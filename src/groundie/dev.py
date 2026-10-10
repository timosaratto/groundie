"""A dev server for the report: serves the page `groundie build` would write for a file of blocks, and
reloads it when the file or the package's source changes."""

import html
import subprocess
import sys
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

from groundie.query import SpecError

PACKAGE = Path(__file__).parent

# Reloads the page when a watched file changes.
RELOAD = """<script>
(async () => {
  const version = await (await fetch("/version")).text();
  setInterval(async () => {
    try { if (await (await fetch("/version")).text() !== version) location.reload(); } catch {}
  }, 500);
})();
</script>
"""


def version(path):
    """Changes when the file of blocks or the package's source does."""
    watched = [path, *(p for p in PACKAGE.rglob("*") if p.suffix in (".py", ".css", ".html"))]
    return str(max(p.stat().st_mtime_ns for p in watched if p.exists()))


def serve(path: Path, port: int):
    path = path.resolve()

    class Handler(BaseHTTPRequestHandler):
        def do_GET(self):
            if self.path == "/":
                # Build in a subprocess, so that edits to the package show without a restart.
                done = subprocess.run([sys.executable, "-m", "groundie.dev", str(path)], capture_output=True, text=True)
                if done.returncode == 0:
                    self.send("text/html", done.stdout.replace("</body>", RELOAD + "</body>"))
                else:
                    self.send("text/html", f"<pre>{html.escape(done.stderr)}</pre>\n{RELOAD}")
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
    from groundie import cli

    try:
        print(cli.build(Path(sys.argv[1])), end="")
    except SpecError as e:
        sys.exit(f"error: {cli.one_line(e)}")
