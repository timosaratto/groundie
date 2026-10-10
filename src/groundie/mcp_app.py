"""An MCP server that renders blocks inline in hosts that support MCP Apps."""

from importlib.resources import files
from pathlib import Path

from mcp.server.apps import Apps, ResourceCsp, client_supports_apps
from mcp.server.mcpserver import Context, MCPServer
from mcp.server.mcpserver.exceptions import ToolError
from mcp_types import CallToolResult, TextContent

from groundie import cli
from groundie.query import SpecError, ident

VIEW = "ui://groundie/blocks.html"
MAX_ROWS = 200  # rows of a block given to the model as text

apps = Apps()


def view_html():
    """The view: view.html with style.css inlined."""
    pkg = files("groundie")
    return pkg.joinpath("view.html").read_text().replace("/* STYLE */", pkg.joinpath("style.css").read_text())


def view_fragment(fragment, manifest):
    """What the view shows: the blocks' HTML with the provenance folded away below it."""
    provenance = f'<details class="provenance"><summary>Provenance</summary>\n{cli.render_provenance(manifest)}</details>\n'
    return fragment + provenance


def tsv(columns, rows, total):
    """At most MAX_ROWS rows as TSV, then how many of 'total' are left out."""
    cell = lambda v: "" if v is None else str(v)
    lines = ["\t".join(columns), *("\t".join(map(cell, row)) for row in rows[:MAX_ROWS])]
    if total > MAX_ROWS:
        lines.append(f"... {total - MAX_ROWS} more rows")
    return "\n".join(lines)


@apps.tool(resource_uri=VIEW, description=(
    "Run groundie blocks and show the user every datatable, matrix and tree among them, with the prose "
    "around them. Returns one summary line per block, the warnings, and the rows of each datatable and "
    "matrix, and of the last block when it is an input, sql or query block. Pass either path, the absolute "
    "path of a Markdown file of blocks, or blocks, the Markdown itself, with base_dir, the absolute "
    "directory that relative input paths resolve against. Call read_guide first for the format."
))
def render_blocks(
    ctx: Context, path: str | None = None, blocks: str | None = None, base_dir: str | None = None,
) -> CallToolResult:
    if (path is None) == (blocks is None):
        raise ToolError("pass exactly one of path or blocks")
    if path is not None and base_dir is not None:
        raise ToolError("base_dir goes with blocks; with path, input paths resolve against the file's directory")
    if path is not None and not Path(path).is_absolute():
        raise ToolError(f"path {path!r} is relative; pass the file's absolute path")
    try:
        if blocks is None:
            fragment, manifest, con = cli.build_report(path)
        else:
            fragment, manifest, con = cli.build_report("blocks", blocks, base_dir)
    except SpecError as e:
        raise ToolError(cli.one_line(e)) from e
    listed = manifest["blocks"]
    text = [cli.summary_line(manifest["file"]["path"], b) for b in listed]
    text += [f"warning: {w}" for w in manifest["warnings"]]
    for i, b in enumerate(listed):
        if b["kind"] in (*cli.TABLES, *cli.MATRICES) or (i == len(listed) - 1 and b["kind"] != "tree"):
            rows = con.execute(f"select * from {ident(b['name'])} limit {MAX_ROWS}").fetchall()
            text.append(f"\n{b['kind']} {b['name']}:\n{tsv(b['output_columns'], rows, b['rows'])}")
    con.close()
    # Hosts without MCP Apps (e.g. Claude Code) hand structuredContent to the model in place of the text.
    structured = {"html": view_fragment(fragment, manifest)} if client_supports_apps(ctx) else None
    return CallToolResult(content=[TextContent(type="text", text="\n".join(text))], structured_content=structured)


apps.add_html_resource(
    VIEW,
    view_html(),
    title="groundie blocks",
    csp=ResourceCsp(resource_domains=["https://unpkg.com"]),
)

server = MCPServer("groundie", extensions=[apps])


@server.tool(description="The guide to writing groundie blocks. Read it before the first render_blocks.")
def read_guide() -> str:
    return files("groundie").joinpath("GUIDE.md").read_text()


def main():
    server.run()
