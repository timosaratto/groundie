"""A dev server for the view: serves view.html as the MCP host would show it, with a stand-in for the host
that renders a file of blocks, and reloads the page when the file or the package's source changes."""

import json
import re
import subprocess
import sys
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

from groundie.query import SpecError

PACKAGE = Path(__file__).parent
IMPORT_APP = re.compile(r'import \{ App \} from "[^"]+";')

# Stands in for the ext-apps App: hands the view the render of the file, as the host hands it the tool
# result, then reloads the page when a watched file changes.
HOST = """
class App {
  async connect() {
    const version = await (await fetch("/version")).text();
    this.ontoolresult(await (await fetch("/render")).json());
    setInterval(async () => {
      try { if (await (await fetch("/version")).text() !== version) location.reload(); } catch {}
    }, 500);
  }
}
"""


def render(path):
    """The tool result the view gets for the file of blocks, or the error as its text."""
    from groundie import cli
    from groundie.mcp_app import view_fragment

    try:
        fragment, manifest, con = cli.build_report(path)
    except SpecError as e:
        return {"content": [{"type": "text", "text": f"error: {cli.one_line(e)}"}]}
    con.close()
    return {"structuredContent": {"html": view_fragment(fragment, manifest)}}


def version(path):
    """Changes when the file of blocks or the package's source does."""
    watched = [path, *(p for p in PACKAGE.rglob("*") if p.suffix in (".py", ".css", ".html"))]
    return str(max(p.stat().st_mtime_ns for p in watched if p.exists()))


def serve(path: Path, port: int):
    path = path.resolve()

    class Handler(BaseHTTPRequestHandler):
        def do_GET(self):
            if self.path == "/":
                # Read fresh, as from a subprocess, so that edits to the package show without a restart.
                done = subprocess.run([sys.executable, "-c", "from groundie.mcp_app import view_html; print(view_html())"],
                                      capture_output=True, text=True)
                if done.returncode == 0:
                    self.send("text/html", IMPORT_APP.sub(HOST, done.stdout))
                else:
                    self.send("text/plain", done.stderr)
            elif self.path == "/render":
                done = subprocess.run([sys.executable, "-m", "groundie.dev", str(path)], capture_output=True, text=True)
                result = done.stdout if done.returncode == 0 else json.dumps(
                    {"content": [{"type": "text", "text": done.stderr}]})
                self.send("application/json", result)
            elif self.path == "/version":
                self.send("text/plain", version(path))
            else:
                self.send_error(404)

        def send(self, content_type, body):
            data = body.encode()
            self.send_response(200)
            self.send_header("Content-Type", f"{content_type}; charset=utf-8")
            self.send_header("Content-Length", str(len(data)))
            self.end_headers()
            self.wfile.write(data)

        def log_message(self, *args):
            pass

    server = ThreadingHTTPServer(("127.0.0.1", port), Handler)
    print(f"serving {path.name} at http://127.0.0.1:{port}/")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass


if __name__ == "__main__":
    print(json.dumps(render(sys.argv[1])))
