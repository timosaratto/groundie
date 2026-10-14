import asyncio
import re
from pathlib import Path

from mcp import Client
from mcp.client import advertise
from mcp.server.apps import APP_MIME_TYPE, EXTENSION_ID

from groundie.mcp_app import MAX_ROWS, VIEW, server

COHORT = Path(__file__).parent.parent / "examples" / "cohort.md"
DATA = (Path(__file__).parent.parent / "data").resolve()
APPS_HOST = [advertise(EXTENSION_ID, {"mimeTypes": [APP_MIME_TYPE]})]
AMRFINDER = "```input amrfinder\nadapter: amrfinder\npaths:\n" + "".join(
    f"  - amrfinder_results/{path.name}\n" for path in sorted((DATA / "amrfinder_results").glob("*.tsv"))) + "```\n"


def call(args, extensions=None):
    async def go():
        async with Client(server, extensions=extensions) as client:
            tools = (await client.list_tools()).tools
            view = (await client.read_resource(VIEW)).contents[0]
            return tools, view, await client.call_tool("render_blocks", args)
    return asyncio.run(go())


def test_render_blocks_for_apps_host():
    tools, view, result = call({"path": str(COHORT)}, APPS_HOST)
    assert {t.name: t for t in tools}["render_blocks"].meta["ui"]["resourceUri"] == VIEW
    assert view.mime_type == APP_MIME_TYPE and "table.matrix" in view.text  # style.css is inlined
    assert not result.is_error
    text = result.content[0].text
    assert text.startswith("cohort.md:3: input mlst: 10 rows; columns: sample, scheme")
    assert "\ncohort.md:52: Tree phylogeny: 10 tips, 10 matched\n" in text
    for shown in ("MlstTable st", "DataTable summary_table", "Matrix resistome"):
        assert f"\n\n{shown}:\n" in text
    assert "\nsql summary:\n" not in text and "\ninput mlst:\n" not in text  # neither shown nor last
    html = result.structured_content["html"]
    kinds = re.findall(r'<section class="block (\w+)">', html)
    assert kinds == ["datatable", "datatable", "tree", "matrix"]  # every display block
    assert "<h1>Staphylococcus aureus cohort</h1>" in html  # and the prose
    assert html.endswith("</footer>\n</details>\n") and "<summary>Provenance</summary>" in html
    assert "<html" not in html


def test_render_blocks_text_only_client():
    _, _, result = call({"path": str(COHORT)})
    assert result.structured_content is None
    assert "Tree phylogeny: 10 tips" in result.content[0].text


def test_render_blocks_inline():
    blocks = AMRFINDER + "\n```DataTable classes\nsource: amrfinder\ngroup_by: [element_class]\n```\n"
    _, _, result = call({"blocks": blocks, "base_dir": str(DATA)}, APPS_HOST)
    assert not result.is_error
    text = result.content[0].text
    assert text.startswith("blocks:1: input amrfinder: 334 rows; columns: sample, contig_id")
    assert "\n\nDataTable classes:\nelement_class\n" in text
    assert '<section class="block datatable">' in result.structured_content["html"]


def test_render_blocks_display_keeps_raw_names():
    blocks = AMRFINDER + ("\n```DataTable classes\nsource: amrfinder\ngroup_by: [element_class]\n"
                          "display: {element_class: Class}\n```\n")
    _, _, result = call({"blocks": blocks, "base_dir": str(DATA)}, APPS_HOST)
    assert "\n\nDataTable classes:\nelement_class\n" in result.content[0].text
    assert '<th title="element_class">Class</th>' in result.structured_content["html"]


def test_render_blocks_inline_paths():
    block = "```input amrfinder\nadapter: amrfinder\npath: {}\n```\n"
    _, _, result = call({"blocks": block.format("amrfinder_results/SAMD00654274.tsv")})
    assert result.is_error and result.content[0].text.endswith(
        "blocks:1: input 'amrfinder': relative path 'amrfinder_results/SAMD00654274.tsv' in inline blocks; "
        "use absolute paths or pass base_dir")
    _, _, result = call({"blocks": block.format(DATA / "amrfinder_results" / "SAMD00654274.tsv")})
    assert not result.is_error


def test_render_blocks_ending_in_an_input():
    _, _, result = call({"blocks": AMRFINDER, "base_dir": str(DATA)}, APPS_HOST)
    text = result.content[0].text
    assert "\n\ninput amrfinder:\nsample\tcontig_id\t" in text
    assert text.endswith(f"\n... {334 - MAX_ROWS} more rows")
    assert text.count("\n") == 1 + 2 + 1 + MAX_ROWS  # summary line, label, header, rows, the rest counted
    assert result.structured_content["html"].startswith('<details class="provenance">')  # nothing to show


def test_render_blocks_ending_in_sql():
    blocks = AMRFINDER + "\n```sql n\nselect count(*) as n from amrfinder\n```\n"
    _, _, result = call({"blocks": blocks, "base_dir": str(DATA)})
    assert result.content[0].text.endswith("\n\nsql n:\nn\n334")
    assert "\ninput amrfinder:\n" not in result.content[0].text  # not the last block


def test_render_blocks_arguments():
    for args, message in [
        ({}, "pass exactly one of path or blocks"),
        ({"path": str(COHORT), "blocks": "# x"}, "pass exactly one of path or blocks"),
        ({"path": str(COHORT), "base_dir": str(DATA)},
         "base_dir goes with blocks; with path, input paths resolve against the file's directory"),
        ({"path": "cohort.md"}, "path 'cohort.md' is relative; pass the file's absolute path"),
        ({"blocks": "```DataTable t\nsource: nope\n```\n"},
         "blocks:1: DataTable 't': source 'nope' names no block; readable here: none"),
    ]:
        _, _, result = call(args)
        assert result.is_error and result.content[0].text.endswith(": " + message)


def test_read_guide():
    async def go():
        async with Client(server) as client:
            return await client.call_tool("read_guide", {})
    text = asyncio.run(go()).content[0].text
    assert text.startswith("# Writing groundie blocks")
    assert "```input " in text and "```DataTable " in text and "render_blocks" in text
