import json
import re
from pathlib import Path

import duckdb
import pytest

from groundie import cli
from groundie import query as q
from groundie.adapters import amrfinder, mash, mlst, utils
from groundie.cli import build
from groundie.mcp_app import tsv

SPEC = Path(__file__).parent.parent / "examples" / "SAMD00654274.md"
COHORT = Path(__file__).parent.parent / "examples" / "cohort.md"
AMRFINDER = (Path(__file__).parent.parent / "data" / "amrfinder_results" / "SAMD00654274.tsv").resolve()
MASH = (Path(__file__).parent.parent / "data" / "tree_results" / "cohort.mash_dist.tsv").resolve()
MLST = (Path(__file__).parent.parent / "data" / "mlst_results" / "mlst.tsv").resolve()
FIXTURES = Path(__file__).parent / "fixtures"


def block(fence, text):
    return f"```{fence}\n{text}```\n"


INPUT = block("input amrfinder", f"adapter: amrfinder\npath: {AMRFINDER}\n")  # 4 lines


def write(tmp_path, text):
    path = tmp_path / "spec.md"
    path.write_text(text)
    return path


def rows(path):
    """The row count of each DataTable and Matrix, an MlstTable included, from the manifest."""
    return {b["name"]: b["rows"] for b in cli.build_report(path)[1]["blocks"]
            if b["kind"] in ("DataTable", "Matrix", "MlstTable")}


def details(out, kind, i=0):
    """What 'View query' expands to under the i-th display block of a kind."""
    return out.split(f'<section class="block {kind}">')[i + 1].split("</summary>")[1].split("</details>")[0]


def spec(tmp_path, body):
    """A file with the amrfinder input, a blank line, then body, which starts on line 6."""
    return write(tmp_path, f"{INPUT}\n{body}")


def test_build_amr_example():
    out = build(SPEC)
    assert out.count("<tr>") == 49  # 1 + 1 + 13 + 5 + 21 datatable rows + 5 header rows + provenance header + 2 inputs
    assert "<title>SAMD00654274</title>" in out
    assert out.count("<h1>") == 1 and "<h1>AMR hits for SAMD00654274</h1>" in out  # the author's, from prose
    assert "<h2>Provenance</h2>" in out
    assert "<h2>AMR hits</h2>" in out
    assert "amrfinder preliminary" in out
    assert rows(SPEC) == {"st": 1, "st_2": 1, "amr": 13, "classes": 5, "virulence": 21}
    assert ('<tr><th title="The genome&#x27;s file, without directories or FASTA extension.">Sample</th>'
            '<th class="num" title="The ST; null when none was assigned.">ST</th>'
            '<th title="calls: arcC">arcC</th>') in out  # display: mlst, then calls spread
    assert '<tr><th title="element_class">Class</th><th class="num" title="hits">Hits</th>' in out
    assert "<h2>Virulence and stress</h2>" in out
    assert "<p>From hits:</p>\n<pre>source: amrfinder\ncolumns: [element_symbol," in out
    assert "<p>From amrfinder:</p>\n<pre>adapter: amrfinder\npath: ../data/amrfinder_results/SAMD00654274.tsv\n" in out
    assert "filters:\n  element_type: [AMR]" in out
    assert "?" not in out.split("<footer")[0]
    assert "<caption>saureus</caption>" in out
    assert '<tr><td>SAMD00654274</td><td class="num">8</td><td>3</td><td>3</td><td>1</td><td>1</td><td>4</td><td>4</td><td>3</td></tr>' in out
    assert "hits (" not in out  # the query block shows nothing
    strip = lambda s: re.sub(r"Built \S+ with", "Built <t> with", s)
    assert strip(build(SPEC)) == strip(out)


def test_build_cohort_example():
    out = build(COHORT)
    assert out.count("<tr>") == 55  # st, summary, resistome: 3 * (1 + 10); provenance: 1 + 21
    assert "<title>cohort</title>" in out
    assert rows(COHORT) == {"st": 10, "summary_table": 10, "resistome": 129}
    assert "<p>From summary (ad hoc SQL):</p>\n<pre>select\n    m.sample," in out
    assert ('<tr><td>SAMD00654274</td><td class="num">8</td><td class="num">13</td><td class="num">19</td>'
            '<td class="num">2</td></tr>') in out  # 13 AMR, 19 virulence, 2 stress in its amrfinder file
    assert '<tr><td>SAMD00654278</td><td class="num null">–</td><td>3</td><td>3</td><td>1</td><td>1</td><td>4</td><td>~4</td>' in out
    assert "<th>arcC</th>" in out and "<th>yqiL</th>" in out
    assert "mlst_results/SAMD00654356.tsv" in out
    assert "<p>From mlst:</p>\n<pre>adapter: mlst\npath: ../data/mlst_results/SAMD*.tsv\n" in out


def test_cohort_resistome_matrix():
    out = build(COHORT)
    # SAMD00654274 is the first row and blaI its first gene column: 98.41% identity, one copy
    assert ('<tr><td>SAMD00654274</td><td class="cell" title="element_class BETA-LACTAM, best_identity 98.4, '
            'lowest_coverage 100.0, copies 1"') in out
    assert out.count('copies 2"') == 2  # the two blaI duplicates
    assert 'best_identity 100.0, lowest_coverage 50.0, copies 2" style="background:#1f4e96;color:#fff"><span class="mark">◆</span>' in out
    query = details(out, "matrix")
    assert "group_by: [sample, element_class, element_symbol]\nmeasures:\n  best_identity: max(pct_identity)\n" in query
    assert "rows: [sample]\ncols: element_symbol\nvalue: best_identity\nscale: [90, 100]\nmark:\n  column: lowest_coverage\n" in query
    assert "<p>From amrfinder:</p>\n<pre>adapter: amrfinder\npath: ../data/amrfinder_results/*.tsv\n" in query
    grid = out.split('<table class="matrix">')[1].split("</table>")[0]
    assert grid.count('<th class="col">') == 32 and grid.count('<td class="cell"') == 129


def test_matrix_layout_names_output_column(tmp_path):
    text = "source: amrfinder\ngroup_by: [sample, element_symbol]\nrows: [sample]\ncols: element_class\n"
    with pytest.raises(ValueError, match="Matrix 'm': 'element_class' not an output column; "
                                         "output columns: sample, element_symbol"):
        build(spec(tmp_path, block("Matrix m", text)))


def test_matrix_value_as_text(tmp_path):
    body = block("input mash", f"adapter: mash\npath: {MASH}\n")
    text = "source: mash\nrows: [reference]\ncols: query\nvalue: distance\n"
    out = build(write(tmp_path, body + block("Matrix m", text + "digits: 4\n")))
    grid = out.split('<table class="matrix">')[1].split("</table>")[0]
    assert re.search(r'<td class="cell num" title="[^"]*">0\.0081</td>', grid)
    assert grid.count('class="cell num"') == 45  # each pair once: one triangle
    assert "cells: distance; blank: no row" in out
    out = build(write(tmp_path, body + block("Matrix m", text)))
    assert ">0.00814823</td>" in out  # without digits, the value as it is


def test_matrix_symmetric(tmp_path):
    body = block("input mash", f"adapter: mash\npath: {MASH}\n")
    text = "source: mash\nrows: [reference]\ncols: query\nvalue: distance\ndigits: 4\nsymmetric: true\n"
    out = build(write(tmp_path, body + block("Matrix m", text)))
    grid = out.split('<table class="matrix">')[1].split("</table>")[0]
    assert grid.count('class="cell num"') == 90  # both triangles; the diagonal is blank
    heads = re.findall(r'<th class="col"><span>([^<]+)</span>', grid)
    assert len(heads) == 10 and heads == sorted(heads) == re.findall(r"<tr><td>([^<]+)</td>", grid)
    assert "symmetric: true" in out  # under View query
    with pytest.raises(ValueError, match="Matrix 'm': 'symmetric' is true or false, and true needs one"):
        build(write(tmp_path, body + block("Matrix m", text.replace("[reference]", "[reference, p_value]"))))


@pytest.mark.parametrize("keys, error", [
    ("scale: [0, 1]\n", "'scale' needs 'value'"),
    ("value: best\nscale: [90, 100]\ndigits: 2\n", "'digits' is a whole number of decimal places"),
    ("value: best\ndigits: -1\n", "'digits' is a whole number of decimal places"),
])
def test_matrix_value_keys(tmp_path, keys, error):
    text = ("source: amrfinder\ngroup_by: [sample, element_symbol]\nmeasures:\n  best: max(pct_identity)\n"
            "rows: [sample]\ncols: element_symbol\n")
    with pytest.raises(ValueError, match=f"Matrix 'm': {re.escape(error)}"):
        build(spec(tmp_path, block("Matrix m", text + keys)))


def test_matrix_rows_sharing_a_cell(tmp_path):
    # SAMD00654319 has two blaI rows, coverage 100 and 50; no grouping needed
    amr = AMRFINDER.with_name("SAMD00654319.tsv")
    text = ("source: amrfinder\ncolumns: [sample, element_symbol, pct_coverage]\nfilters:\n  element_type: AMR\n"
            "sort: [pct_coverage]\nrows: [element_symbol]\ncols: sample\n")
    head = block("input amrfinder", f"adapter: amrfinder\npath: {amr}\n")
    out = build(write(tmp_path, head + block("Matrix m", text)))
    assert '<td class="cell" title="2 rows; shown: pct_coverage 50.0">●</td>' in out  # the first
    out = build(write(tmp_path, head + block("Matrix m", text + "value: pct_coverage\nscale: [0, 100]\n")))
    assert 'title="2 rows; shown: pct_coverage 50.0"' in out  # still the first: a value does not pick
    out = build(write(tmp_path, head + block("Matrix m", text + "pick: max(pct_coverage)\n")))
    assert 'title="2 rows; shown: pct_coverage 100.0"' in out
    with pytest.raises(cli.SpecError, match="Matrix 'm': 'pick' must be max\\(<column>\\) or min"):
        build(write(tmp_path, head + block("Matrix m", text + "pick: pct_coverage\n")))


def test_amrfinder_matrix(tmp_path):
    amr = AMRFINDER.with_name("SAMD00654319.tsv")
    head = block("input amrfinder", f"adapter: amrfinder\npaths: [{amr}, {AMRFINDER}]\n")
    out = build(write(tmp_path, head + block("AmrfinderMatrix r", "source: amrfinder\nfilters: {element_type: AMR}\n")))
    grid = out.split('<table class="matrix">')[1].split("</table>")[0]
    assert grid.startswith("\n<thead>\n<tr><th>element_symbol</th>")  # one type: no element_type column
    # SAMD00654319's two blaI hits: the full-length one is shown, so no ◆
    assert 'title="2 rows; shown: pct_identity 100.0, pct_coverage 100.0">●</td>' in grid
    assert "pick: max(pct_coverage)" in out and "columns: [element_symbol, sample, pct_identity, pct_coverage]" in out
    out = build(write(tmp_path, head + block("AmrfinderMatrix r", "source: amrfinder\n")))  # every type
    assert "<tr><th>element_type</th><th>element_symbol</th>" in out and "<tr><td>STRESS</td>" in out
    marked = lambda out: out.split('<table class="matrix">')[1].split("</table>")[0].count("◆")
    default = marked(build(write(tmp_path, head + block("AmrfinderMatrix r", "source: amrfinder\n"))))
    out = build(write(tmp_path, head + block("AmrfinderMatrix r", "source: amrfinder\npartial_below: 0\n")))
    assert default > 0 and marked(out) == 0 and "◆ pct_coverage below 0" in out


@pytest.mark.parametrize("body, error", [
    ("source: amrfinder\nrows: [sample]\n", "AmrfinderMatrix 'r': unknown key 'rows'; accepted: source, filters, sort, partial_below"),
    ("source: amrfinder\npartial_below: high\n", "AmrfinderMatrix 'r': 'partial_below' must be a number"),
    ("source: t\n", "AmrfinderMatrix 'r': source 't' is loaded by adapter 'tsv'; an AmrfinderMatrix needs an input"),
])
def test_amrfinder_matrix_shape(tmp_path, body, error):
    (tmp_path / "t.tsv").write_text("n\n1\n")
    text = INPUT + block("input t", "adapter: tsv\npath: t.tsv\n") + block("AmrfinderMatrix r", body)
    with pytest.raises(cli.SpecError, match=re.escape(error)):
        build(write(tmp_path, text))


SAUREUS = "arcC(3)\taroE(3)\tglpF(1)\tgmk({})\tpta({})\ttpi({})\tyqiL(3)"
ECOLI = "adk(53)\tfumC(40)\tgyrB(47)\ticd(13)\tmdh(36)\tpurA(28)\trecA(29)"
MLST_ROWS = [  # two schemes and a genome with none, mlst's allele flags, and not in name order
    f"a/E2.fa\tecoli_achtman_4\t131\t{ECOLI}",
    f"S1.fa\tsaureus\t8\t{SAUREUS.format(1, 4, '~4')}",
    f"S2.fa.gz\tsaureus\t-\t{SAUREUS.format('-', '4?', '4,5')}",
    f"E1.fa\tecoli_achtman_4\t10\t{ECOLI}",
    f"S0.fa\tsaureus\t2725\t{SAUREUS.format(1, 4, 4)}",
    "unknown.fa\t-\t-",
]


def mlst_spec(tmp_path, body, lines=MLST_ROWS):
    (tmp_path / "m.tsv").write_text("".join(f"{line}\n" for line in lines))
    return write(tmp_path, block("input mlst", "adapter: mlst\npath: m.tsv\n") + body)


def captions(out):
    return re.findall(r"<caption>([^<]+)</caption>", out)


def first_cells(out):
    return re.findall(r"<tr><td>([^<]+)</td>", out.split("<footer")[0])


def test_mlst_table(tmp_path):
    path = mlst_spec(tmp_path, block("MlstTable st", "source: mlst\n"))
    out = build(path)
    # by size, then by name; no scheme last; rows by sample
    assert captions(out) == ["saureus", "ecoli_achtman_4", "No scheme found"]
    assert first_cells(out) == ["S0", "S1", "S2", "E1", "E2", "unknown"]
    assert ('<tr><th title="sample">Sample</th><th class="num" title="sequence_type">ST</th><th>arcC</th>'
            '<th>aroE</th><th>glpF</th><th>gmk</th>') in out
    assert "<th>adk</th>" in out and "<th>recA</th>" in out
    assert '<tr><td>S1</td><td class="num">8</td><td>3</td><td>3</td><td>1</td><td>1</td><td>4</td><td>~4</td><td>3</td></tr>' in out
    assert '<tr><td>S2</td><td class="num null">–</td><td>3</td><td>3</td><td>1</td><td class="null">–</td><td>4?</td><td>4,5</td><td>3</td></tr>' in out
    no_scheme = out.split("<caption>No scheme found</caption>")[1].split("</table>")[0]
    assert '<th title="sample">Sample</th><th title="sequence_type">ST</th></tr>' in no_scheme  # no loci, no ST
    assert "<tr><td>unknown</td><td class=\"null\">–</td></tr>" in no_scheme
    assert out.count('<section class="block datatable">') == 1 and out.count("View query") == 1
    assert "Tables by how many genomes they hold" in out
    page, manifest, con = cli.build_report(path)
    assert manifest["blocks"][1] == {"line": 5, "kind": "MlstTable", "name": "st", "rows": 6, "tables": 3,
                                     "output_columns": ["sample", "scheme", "sequence_type", "calls"]}
    assert [(t["caption"], t["rows"], len(t["columns"])) for t in manifest["queries"][0]["tables"]] == [
        ("saureus", 3, 9), ("ecoli_achtman_4", 2, 9), ("No scheme found", 1, 2)]
    assert cli.summary_line("spec.md", manifest["blocks"][1]).startswith("spec.md:5: MlstTable st: 6 rows in 3 tables; ")


def test_mlst_table_loci_per_scheme(tmp_path):
    lines = [*MLST_ROWS, "B1.fa\tborrelia\t1\t" + "\t".join(f"g{i}(1)" for i in range(1, 9))]  # eight loci
    out = build(mlst_spec(tmp_path, block("MlstTable st", "source: mlst\n"), lines))
    borrelia = out.split("<caption>borrelia</caption>")[1].split("</table>")[0]
    assert borrelia.count("</th>") == 10 and "<th>g8</th>" in borrelia
    saureus = out.split("<caption>saureus</caption>")[1].split("</table>")[0]
    assert saureus.count("</th>") == 9 and "g8" not in saureus  # its own seven loci only


def test_mlst_table_sort_and_filters(tmp_path):
    out = build(mlst_spec(tmp_path, block("MlstTable st", "source: mlst\nsort: [sequence_type]\n")))
    assert captions(out) == ["saureus", "ecoli_achtman_4", "No scheme found"]  # first seen, no scheme still last
    assert first_cells(out) == ["S1", "S0", "S2", "E1", "E2", "unknown"]  # by number: 8, 10, 131, 2725, then no ST
    assert "Tables by how many" not in out
    out = build(mlst_spec(tmp_path, block("MlstTable st", "source: mlst\nfilters: {scheme: ecoli_achtman_4}\n")))
    assert captions(out) == ["ecoli_achtman_4"] and first_cells(out) == ["E1", "E2"]
    out = build(mlst_spec(tmp_path, block("MlstTable st", "source: mlst\nfilters: {sample: nope}\n")))
    assert captions(out) == [] and '<td class="empty" colspan="4">No rows</td>' in out


@pytest.mark.parametrize("filters, samples", [
    ("{sequence_type: null}", ["S2", "unknown"]),
    ("{sequence_type: ~}", ["S2", "unknown"]),
    ("{sequence_type: [8, null]}", ["S1", "S2", "unknown"]),
    ("{scheme: null}", ["unknown"]),
    ("{sequence_type: None}", None),  # the text 'None', not a null
])
def test_filter_null(tmp_path, filters, samples):
    body = block("DataTable t", f"source: mlst\nfilters: {filters}\nsort: [sample]\n")
    if samples is None:
        with pytest.raises(cli.SpecError):
            build(mlst_spec(tmp_path, body))
        return
    page, manifest, con = cli.build_report(mlst_spec(tmp_path, body))
    assert [r[0] for r in con.execute("select sample from t").fetchall()] == samples


def test_filter_null_text():
    query = q.normalize("t", {"source": "mlst", "filters": {"sequence_type": [8, None]}})
    assert q.to_sql(query) == ('select * from "mlst" where ("sequence_type" in (?) or "sequence_type" is null)', ["8"])
    assert q.filter_text(query) == "sequence_type in 8, null"
    assert "sequence_type: ['8', null]" in q.to_yaml(query)


def test_mlst_table_as_source(tmp_path):
    # a later block reads the query's rows, raw, not the split tables
    body = block("MlstTable st", "source: mlst\n") + block("DataTable t", "source: st\nfilters: {sample: S1}\n")
    page, manifest, con = cli.build_report(mlst_spec(tmp_path, body))
    assert con.execute("select * from t").fetchall() == [
        ("S1", "saureus", 8, [{"locus": locus, "call": call} for locus, call in zip(
            ["arcC", "aroE", "glpF", "gmk", "pta", "tpi", "yqiL"], ["3", "3", "1", "1", "4", "~4", "3"])])]


@pytest.mark.parametrize("body, error", [
    ("source: mlst\ncolumns: [sample]\n", "MlstTable 'st': unknown key 'columns'; accepted: source, filters, sort"),
    ("source: t\n", "MlstTable 'st': source 't' is loaded by adapter 'tsv'; an MlstTable needs an input loaded by the mlst adapter"),
    ("source: q\n", "MlstTable 'st': source 'q' is a query block; an MlstTable needs an input"),
])
def test_mlst_table_shape(tmp_path, body, error):
    (tmp_path / "t.tsv").write_text("n\n1\n")
    head = block("input t", "adapter: tsv\npath: t.tsv\n") + block("query q", "source: mlst\n")
    with pytest.raises(cli.SpecError, match=re.escape(error)):
        build(mlst_spec(tmp_path, head + block("MlstTable st", body)))


def test_datatable_display(tmp_path):
    body = "source: mlst\ncolumns: [sample, sequence_type]\nsort: [sample]\ndisplay:\n  sequence_type: ST\n"
    path = mlst_spec(tmp_path, block("DataTable st", body))
    out = build(path)
    assert '<tr><th>sample</th><th class="num" title="sequence_type">ST</th></tr>' in out
    assert details(out, "datatable").startswith(
        "\n<pre>source: mlst\ncolumns: [sample, sequence_type]\nsort: [sample]\ndisplay:\n  sequence_type: ST\n</pre>")
    page, manifest, con = cli.build_report(path)
    assert manifest["queries"][0]["display"] == {"sequence_type": "ST"}
    assert manifest["blocks"][1]["output_columns"] == ["sample", "sequence_type"]  # raw names
    assert [c[0] for c in con.execute("describe st").fetchall()] == ["sample", "sequence_type"]


def test_datatable_display_descriptions(tmp_path):
    body = ("source: mlst\ncolumns: [sample, scheme, sequence_type]\nsort: [sample]\ndisplay:\n"
            "  sample: {description: The <genome>}\n  scheme: Scheme\n"
            "  sequence_type: {title: ST, description: The sequence type}\n")
    path = mlst_spec(tmp_path, block("DataTable st", body))
    out = build(path)
    # a description replaces the column's name on hover; a title alone keeps it there
    assert ('<tr><th title="The &lt;genome&gt;">sample</th><th title="scheme">Scheme</th>'
            '<th class="num" title="The sequence type">ST</th></tr>') in out
    assert ("  scheme: Scheme\n  sequence_type:\n    title: ST\n    description: The sequence type\n"
            in details(out, "datatable"))  # as written
    page, manifest, con = cli.build_report(path)
    assert manifest["queries"][0]["display"] == {"scheme": "Scheme", "sequence_type": "ST"}
    assert manifest["queries"][0]["descriptions"] == {"sample": "The <genome>", "sequence_type": "The sequence type"}


def test_datatable_display_every_column(tmp_path):
    # without 'columns', display may name any column of the source
    out = build(mlst_spec(tmp_path, block("DataTable st", "source: mlst\ndisplay: {scheme: Scheme}\n")))
    assert '<th title="scheme">Scheme</th>' in out


@pytest.mark.parametrize("body, error", [
    ("columns: [sample, sequence_type]\ndisplay: {st: ST}\n",
     "spec.md:5: DataTable 'st': display names 'st', not an output column; output columns: sample, sequence_type"),
    ("display: {st: ST}\n", "display names 'st', not an output column; output columns: sample, scheme, sequence_type, calls"),
    ("display: {sequence_type: 7}\n", "DataTable 'st': display 'sequence_type' must be text, the column's title, "
                                      "or {title, description}, either or both, each text"),
    ("display: {sequence_type: ''}\n", "display 'sequence_type' must be text, the column's title"),
    ("display: {sequence_type: {}}\n", "display 'sequence_type' must be text, the column's title"),
    ("display: {sequence_type: {title: ST, digits: 1}}\n", "display 'sequence_type' must be text, the column's title"),
    ("display: {sequence_type: {description: 7}}\n", "display 'sequence_type' must be text, the column's title"),
    ("display: {sequence_type: {title: ' '}}\n", "display 'sequence_type' must be text, the column's title"),
    ("display: {sequence_type: [ST]}\n", "display 'sequence_type' must be text, the column's title"),
    ("display: {st: {title: ST}}\n", "display names 'st', not an output column"),
    ("display: [sequence_type]\n", "DataTable 'st': 'display' maps columns to their titles, or to "
                                   "{title, description}, e.g. {sequence_type: ST}, or names an input whose "
                                   "adapter's titles and descriptions it shows"),
])
def test_datatable_display_shape(tmp_path, body, error):
    with pytest.raises(cli.SpecError, match=re.escape(error)):
        build(mlst_spec(tmp_path, block("DataTable st", "source: mlst\n" + body)))


def test_datatable_display_from_input(tmp_path):
    # the adapter's titles and descriptions, by column name, through a sql block; a column it has neither
    # for keeps its name
    body = (block("sql n", "select sample, sequence_type, 1 as one from mlst order by sample\n")
            + block("DataTable st", "source: n\ndisplay: mlst\n"))
    path = mlst_spec(tmp_path, body)
    out = build(path)
    assert ('<tr><th title="The genome&#x27;s file, without directories or FASTA extension.">'
            'Sample</th><th class="num" title="The ST; null when none was assigned.">ST</th>'
            '<th class="num">one</th></tr>') in out
    assert details(out, "datatable").startswith("\n<pre>source: n\ndisplay: mlst\n</pre>")  # as written
    page, manifest, con = cli.build_report(path)
    assert manifest["queries"][0]["display"] == {"sample": "Sample", "sequence_type": "ST"}  # as resolved
    assert manifest["queries"][0]["descriptions"] == {
        "sample": "The genome's file, without directories or FASTA extension.",
        "sequence_type": "The ST; null when none was assigned."}
    assert [c[0] for c in con.execute("describe st").fetchall()] == ["sample", "sequence_type", "one"]


@pytest.mark.parametrize("display, error", [
    ("nope", "spec.md:12: DataTable 'st': display 'nope' names no block; inputs: mlst, t"),
    ("q", "DataTable 'st': display 'q' is a query block; display names an input, whose adapter's titles and "
          "descriptions it shows"),
    ("t", "DataTable 'st': display 't' is loaded by adapter 'tsv', which declares no titles or descriptions; "
          "adapters with them: amrfinder, mlst"),
])
def test_datatable_display_from_input_shape(tmp_path, display, error):
    (tmp_path / "t.tsv").write_text("n\n1\n")
    head = block("input t", "adapter: tsv\npath: t.tsv\n") + block("query q", "source: mlst\n")
    with pytest.raises(cli.SpecError, match=re.escape(error)):
        build(mlst_spec(tmp_path, head + block("DataTable st", f"source: mlst\ndisplay: {display}\n")))


@pytest.mark.parametrize("fence, body", [
    ("query q", "source: mlst\n"),
    ("Matrix q", "source: mlst\nrows: sample\ncols: scheme\n"),
    ("MlstTable q", "source: mlst\n"),
])
def test_display_only_on_datatable(tmp_path, fence, body):
    with pytest.raises(cli.SpecError, match="unknown key 'display'"):
        build(mlst_spec(tmp_path, block(fence, body + "display: {sample: Sample}\n")))


SPREAD = "spread: {column: calls, names: locus, values: call}\n"
ECOLI_LOCI = ["adk", "fumC", "gyrB", "icd", "mdh", "purA", "recA"]
SAUREUS_LOCI = ["arcC", "aroE", "glpF", "gmk", "pta", "tpi", "yqiL"]


def test_datatable_spread(tmp_path):
    body = "source: mlst\ncolumns: [sample, calls, sequence_type]\nsort: [sample]\n" + SPREAD
    path = mlst_spec(tmp_path, block("DataTable st", body))
    out = build(path)
    # in place of calls, a column per locus, by first appearance: E1's loci, then S0's
    heads = re.findall(r"<th[^>]*>([^<]+)</th>", out.split("<footer")[0])
    assert heads == ["sample", *ECOLI_LOCI, *SAUREUS_LOCI, "sequence_type"]
    assert '<th title="calls: adk">adk</th>' in out
    blank = "<td></td>" * 7
    assert ("<tr><td>E1</td><td>53</td><td>40</td><td>47</td><td>13</td><td>36</td><td>28</td><td>29</td>"
            f'{blank}<td class="num">10</td></tr>') in out
    # a null call is a dash; a locus the genome's scheme lacks is blank
    assert (f'<tr><td>S2</td>{blank}<td>3</td><td>3</td><td>1</td><td class="null">–</td><td>4?</td>'
            '<td>4,5</td><td>3</td><td class="num null">–</td></tr>') in out
    assert f'<tr><td>unknown</td>{blank * 2}<td class="num null">–</td></tr>' in out  # no calls: all blank
    assert details(out, "datatable").startswith(
        "\n<pre>source: mlst\ncolumns: [sample, calls, sequence_type]\nsort: [sample]\n"
        "spread:\n  column: calls\n  names: locus\n  values: call\n</pre>")
    page, manifest, con = cli.build_report(path)
    assert manifest["queries"][0]["spread"] == {"column": "calls", "names": "locus", "values": "call",
                                                "columns": [*ECOLI_LOCI, *SAUREUS_LOCI]}
    # the SQL table, the summary and the TSV keep calls as a list
    assert manifest["blocks"][1]["output_columns"] == ["sample", "calls", "sequence_type"]
    assert con.execute("select calls from st where sample = 'S0'").fetchone()[0][0] == {"locus": "arcC", "call": "3"}
    data = con.execute("select * from st").fetchall()
    text = tsv(manifest["blocks"][1]["output_columns"], data, len(data))
    assert text.startswith("sample\tcalls\tsequence_type\n") and "\nS0\t[{'locus': 'arcC', 'call': '3'}, " in text


def test_datatable_spread_descriptions(tmp_path):
    body = "source: mlst\ncolumns: [sample, sequence_type, calls]\nfilters: {sample: S1}\n" + SPREAD
    out = build(mlst_spec(tmp_path, block("DataTable st", body + "display: mlst\n")))
    # the hover names the list column and the name; the list column's description is not shown
    assert '<th title="calls: arcC">arcC</th>' in out and "One call per locus" not in out
    assert details(out, "datatable").startswith(
        "\n<pre>source: mlst\nfilters:\n  sample: [S1]\ncolumns: [sample, sequence_type, calls]\n"
        "spread:\n  column: calls\n  names: locus\n  values: call\ndisplay: mlst\n</pre>")
    # a title or description on the list column is not shown, and not an error
    path = mlst_spec(tmp_path, block("DataTable st", body + "display:\n  calls: {title: Calls, description: The calls}\n"))
    out = build(path)
    assert '<th title="calls: arcC">arcC</th>' in out and "Calls<" not in out and 'title="The calls' not in out
    entry = cli.build_report(path)[1]["queries"][0]
    assert "display" not in entry and "descriptions" not in entry


def test_datatable_spread_null_list_and_sql_source(tmp_path):
    sql = ("select sample, case when sample = 'S1' then null else calls end as calls from mlst "
           "where scheme = 'saureus' order by sample\n")
    out = build(mlst_spec(tmp_path, block("sql n", sql) + block("DataTable st", "source: n\n" + SPREAD)))
    assert '<tr><td>S1</td>' + '<td class="null">–</td>' * 7 + "</tr>" in out  # a null list: a dash in each
    assert "<tr><td>S0</td><td>3</td><td>3</td><td>1</td><td>1</td><td>4</td><td>4</td><td>3</td></tr>" in out


def test_datatable_spread_fixed_size_array(tmp_path):
    sql = ("select 'x' as sample, [{'locus': 'a', 'call': '1'}, {'locus': 'b', 'call': '2'}]"
           "::STRUCT(locus VARCHAR, call VARCHAR)[2] as calls\n")
    out = build(mlst_spec(tmp_path, block("sql n", sql) + block("DataTable st", "source: n\n" + SPREAD)))
    assert "<tr><td>x</td><td>1</td><td>2</td></tr>" in out


def test_datatable_spread_query_source(tmp_path):
    body = (block("query q", "source: mlst\ncolumns: [sample, calls]\nfilters: {scheme: ecoli_achtman_4}\n")
            + block("DataTable st", "source: q\nsort: [sample]\n" + SPREAD))
    out = build(mlst_spec(tmp_path, body))
    assert re.findall(r"<th[^>]*>([^<]+)</th>", out.split("<footer")[0]) == ["sample", *ECOLI_LOCI]


def test_datatable_spread_no_entries(tmp_path):
    # every list empty: no spread columns, and calls is gone
    path = mlst_spec(tmp_path, block("DataTable st", "source: mlst\ncolumns: [sample, calls]\n"
                                                     "filters: {scheme: null}\n" + SPREAD))
    out = build(path)
    assert "<tr><th>sample</th></tr>" in out and "<tr><td>unknown</td></tr>" in out
    assert cli.build_report(path)[1]["queries"][0]["spread"]["columns"] == []


def test_spread_missing_is_not_numeric():
    columns = [("k", "k"), ("l", "l")]
    rows = [("a", [{"n": "x", "v": 1}, {"n": "y", "v": None}]), ("b", [{"n": "y", "v": 2}]), ("c", None)]
    shown, spread = cli.tables.spread("DataTable 't'", columns, rows, "l", "n", "v")
    assert shown == [("k", "k"), ("x", "x"), ("y", "y")]
    assert spread == [("a", 1, None), ("b", cli.tables.MISSING, 2), ("c", None, None)]
    out = cli.render_datatable(shown, spread)
    assert '<th class="num">x</th><th class="num">y</th>' in out
    assert '<tr><td>b</td><td></td><td class="num">2</td></tr>' in out


@pytest.mark.parametrize("body, error", [
    ("spread: calls\n", "spec.md:5: DataTable 'st': 'spread' must be {column, names, values}, each text: the list "
                        "column, and the fields of its structs naming the columns it spreads into and holding their "
                        "values, e.g. {column: calls, names: locus, values: call}"),
    ("spread: {column: calls, names: locus}\n", "DataTable 'st': 'spread' must be {column, names, values}"),
    ("spread: {column: calls, names: locus, values: call, fill: x}\n", "DataTable 'st': 'spread' must be"),
    ("spread: {column: calls, names: locus, values: 1}\n", "DataTable 'st': 'spread' must be"),
    ("columns: [sample, sequence_type]\n" + SPREAD,
     "spec.md:5: DataTable 'st': spread names 'calls', not an output column; output columns: sample, sequence_type"),
    ("sort: [sample]\nspread: {column: sample, names: locus, values: call}\n",
     "DataTable 'st': spread 'sample': row 1 holds 'E1', not a list of {locus, call}"),
    ("sort: [sample]\nspread: {column: calls, names: name, values: call}\n",
     "DataTable 'st': spread 'calls': row 1 holds [{'locus': 'adk', 'call': '53'}"),
    ("sort: [sample]\nspread: {column: calls, names: call, values: locus}\n",
     "DataTable 'st': spread 'calls': row 3 has call '3' twice"),  # S0: arcC(3), aroE(3)
])
def test_datatable_spread_shape(tmp_path, body, error):
    with pytest.raises(cli.SpecError, match=re.escape(error)):
        build(mlst_spec(tmp_path, block("DataTable st", "source: mlst\n" + body)))


@pytest.mark.parametrize("calls, error", [
    ("[{'locus': 'a', 'call': '1'}, {'locus': null, 'call': '2'}]", "spread 'calls': row 1 has a null locus"),
    ("[{'locus': 'sample', 'call': '1'}]",
     "spread 'calls': locus 'sample' is the name of another column; other columns: sample"),
    ("[1, 2]", "spread 'calls': row 1 holds [1, 2], not a list of {locus, call}"),
])
def test_datatable_spread_cells(tmp_path, calls, error):
    body = block("sql n", f"select 'x' as sample, {calls} as calls\n") + block("DataTable st", "source: n\n" + SPREAD)
    with pytest.raises(cli.SpecError, match=re.escape(f"DataTable 'st': {error}")):
        build(mlst_spec(tmp_path, body))


@pytest.mark.parametrize("fence, body", [
    ("query q", "source: mlst\n"),
    ("Matrix q", "source: mlst\nrows: sample\ncols: scheme\n"),
    ("MlstTable q", "source: mlst\n"),
])
def test_spread_only_on_datatable(tmp_path, fence, body):
    with pytest.raises(cli.SpecError, match=f"{fence.split()[0]} 'q': unknown key 'spread'"):
        build(mlst_spec(tmp_path, block(fence, body + SPREAD)))


def test_mlst_table_loci_disagree(tmp_path):
    lines = [MLST_ROWS[1], MLST_ROWS[4].replace("tpi(4)", "tpx(4)")]
    with pytest.raises(cli.SpecError, match=re.escape(
            "MlstTable 'st': scheme 'saureus': genomes call different loci: arcC, aroE, glpF, gmk, pta, tpx, yqiL "
            "in S0; arcC, aroE, glpF, gmk, pta, tpi, yqiL in S1")):
        build(mlst_spec(tmp_path, block("MlstTable st", "source: mlst\nsort: [sample]\n"), lines))


def test_matrix_default_order(tmp_path):
    amr = AMRFINDER.with_name("SAMD00654319.tsv")
    head = (block("input amrfinder", f"adapter: amrfinder\npaths: [{amr}, {AMRFINDER}]\n")
            + block("sql hits", "select * from amrfinder where element_type = 'AMR' order by element_symbol desc\n"))
    out = build(write(tmp_path, head + block("Matrix m", "source: hits\nrows: [element_symbol]\ncols: sample\n")))
    grid = out.split('<table class="matrix">')[1].split("</table>")[0]
    cols = re.findall(r'<th class="col"><span>(\w+)', grid)
    rows = re.findall(r"<tr><td>([^<]+)</td>", grid)
    assert cols == ["SAMD00654274", "SAMD00654319"]  # by name, not as the query ordered them
    both = [r for r in rows if grid.split(f"<tr><td>{r}</td>")[1].split("</tr>")[0].count("●") == 2]
    assert rows[:len(both)] == sorted(both) and len(both) < len(rows)  # in both genomes first, then by name
    assert rows[len(both):] == sorted(rows[len(both):])
    assert "Rows by how many cells they fill, most first" in out
    out = build(write(tmp_path, head + block("Matrix m", "source: hits\nsort: [element_symbol]\n"
                                                          "rows: [element_symbol]\ncols: sample\n")))
    grid = out.split('<table class="matrix">')[1].split("</table>")[0]
    assert re.findall(r"<tr><td>([^<]+)</td>", grid) == sorted(rows)  # the query's order
    assert "Rows by how many cells" not in out


def test_matrix_without_value(tmp_path):
    text = "source: amrfinder\ngroup_by: [element_type, sample]\nsort: [element_type]\nrows: [element_type]\ncols: sample\n"
    out = build(spec(tmp_path, block("Matrix m", text)))
    assert '<td class="cell" title="">●</td>' in out


def test_unknown_block(tmp_path):
    with pytest.raises(ValueError, match="unknown block 'chart'.*known: input, sql, query, clusters, checks, calls, DataTable, Matrix, "
                                         "AmrfinderMatrix, MlstTable, Tree"):
        build(spec(tmp_path, block("chart c", "source: amrfinder\n")))


@pytest.mark.parametrize("old, new", [("datatable", "DataTable"), ("mlst_table", "MlstTable"), ("tree", "Tree")])
def test_lowercase_display_block(tmp_path, old, new):
    with pytest.raises(ValueError, match=f"unknown block '{old}' \\({old} x\\); did you mean {new}\\? known: "):
        build(spec(tmp_path, block(f"{old} x", "source: amrfinder\n")))


@pytest.mark.parametrize("fence", ["table t", "data t"])
def test_old_block_kinds(tmp_path, fence):
    with pytest.raises(ValueError, match=f"^spec.md:6: unknown block '{fence.split()[0]}'"):
        build(spec(tmp_path, block(fence, "source: amrfinder\n")))


def test_tsv_over_amrfinder(tmp_path):
    body = (block("input amr", f"adapter: tsv\npath: {AMRFINDER}\noptions: {{null: NA}}\n")
            + block("DataTable hits", "source: amr\nfilters:\n  Type: AMR\ncolumns: [file, Element symbol, Class]\n"))
    path = write(tmp_path, body)
    out = build(path)
    assert rows(path) == {"hits": 13}
    assert "<tr><td>SAMD00654274.tsv</td>" in out
    assert "<th>Element symbol</th>" in out
    assert "tsv preliminary" in out
    assert "<pre>adapter: tsv\npath: " in out and "options:\n  &#x27;null&#x27;: NA\n</pre>" in out


def test_tsv_paths(tmp_path):
    second = AMRFINDER.with_name("SAMD00654276.tsv")
    body = (block("input amr", f"adapter: tsv\npaths: [{AMRFINDER}, {second}]\n")
            + block("DataTable files", "source: amr\ncolumns: [file, Name]\nsort: [file]\n"))
    out = build(write(tmp_path, body))
    assert "<tr><td>SAMD00654274.tsv</td><td>SAMD00654274</td></tr>" in out
    assert "<tr><td>SAMD00654276.tsv</td><td>SAMD00654276</td></tr>" in out
    assert "<p>From amr:</p>\n<pre>adapter: tsv\npaths: [" in out


def test_csv(tmp_path):
    (tmp_path / "t.csv").write_text("gene,drug class\nmecA,beta-lactam\nblaZ,beta-lactam\nnorA,efflux\n")
    body = block("input t", "adapter: csv\npath: t.csv\n") + block("DataTable x", "source: t\nfilters:\n  gene: norA\n")
    path = write(tmp_path, body)
    assert "<tr><td>t.csv</td><td>norA</td><td>efflux</td></tr>" in build(path)
    assert rows(path) == {"x": 1}


@pytest.mark.parametrize("text", ['{"gene": "mecA", "n": 1}\n{"gene": "norA", "n": 2}\n',
                                  '[{"gene": "mecA", "n": 1}, {"gene": "norA", "n": 2}]'])
def test_json(tmp_path, text):
    (tmp_path / "t.json").write_text(text)
    body = block("input t", "adapter: json\npath: t.json\n") + block("DataTable x", "source: t\nsort: [gene]\n")
    out = build(write(tmp_path, body))
    assert '<tr><td>t.json</td><td>mecA</td><td class="num">1</td></tr>' in out
    assert '<tr><td>t.json</td><td>norA</td><td class="num">2</td></tr>' in out


def test_unknown_option(tmp_path):
    text = f"adapter: tsv\npath: {AMRFINDER}\noptions: {{delim: ','}}\n"
    with pytest.raises(ValueError, match="^spec.md:1: input 'amr': unknown option 'delim'.*accepted: null, header"):
        build(write(tmp_path, block("input amr", text)))


def test_options_on_builtin(tmp_path):
    text = f"adapter: amrfinder\npath: {AMRFINDER}\noptions: {{header: false}}\n"
    with pytest.raises(ValueError, match="input 'amr': unknown option 'header'.*accepted: none"):
        build(write(tmp_path, block("input amr", text)))


def test_existing_file_column(tmp_path):
    (tmp_path / "t.csv").write_text("file,n\na,1\n")
    with pytest.raises(ValueError, match="input 't':.*column \"file\""):
        build(write(tmp_path, block("input t", "adapter: csv\npath: t.csv\n")))


@pytest.mark.parametrize("text, error", [
    ("path: x.tsv\n", "needs 'adapter'; known: amrfinder, mlst, mash, tsv, csv, json, newick"),
    ("adapter: xlsx\npath: x.tsv\n", "unknown adapter 'xlsx'; known: amrfinder"),
    ("adapter: tsv\n", "give exactly one of 'path' or 'paths'"),
    ("adapter: tsv\npath: x.tsv\npaths: [x.tsv]\n", "give exactly one of 'path' or 'paths'"),
    ("adapter: tsv\npaths: []\n", "'paths' must be a non-empty list"),
    ("adapter: tsv\npaths: x.tsv\n", "'paths' must be a non-empty list"),
    ("adapter: tsv\npath: x.tsv\nsource: y\n", "unknown key 'source'; accepted: adapter, path, paths, options"),
    ("- adapter\n", "the body must be YAML keys and values"),
])
def test_input_block_shape(tmp_path, text, error):
    with pytest.raises(cli.SpecError, match=f"^spec.md:3: input 'x': {re.escape(error)}"):
        build(write(tmp_path, "prose\n\n" + block("input x", text)))


@pytest.mark.parametrize("path, problem", [(".", "is not a file"), ("nope.tsv", "is not a file"),
                                           ("*.csv", "matches no file"), ("*.d", "matches no file")])
def test_input_names_files(tmp_path, path, problem):
    (tmp_path / "a.tsv").write_text("n\n1\n")
    (tmp_path / "a.tsv.d").mkdir()  # a directory a glob matches: not a file
    with pytest.raises(cli.SpecError, match=f"^spec.md:1: input 'x': {re.escape(repr(path))} {problem}$"):
        build(write(tmp_path, block("input x", f"adapter: tsv\npaths: ['{path}', a.tsv]\n")))


def test_input_glob(tmp_path):
    (tmp_path / "sub").mkdir()
    for name in ("b.tsv", "a.tsv", "c.csv"):
        (tmp_path / "sub" / name).write_text("n\n1\n")
    text = block("input t", "adapter: tsv\npaths: ['sub/*.tsv', sub/a.tsv]\n") + block("DataTable d", "source: t\n")
    page, manifest, _ = cli.build_report(write(tmp_path, text))
    files = manifest["inputs"]["t"]["files"]
    assert [f["path"] for f in files] == ["sub/a.tsv", "sub/b.tsv"]  # name order, each file once
    assert all(len(f["sha256"]) == 64 for f in files) and manifest["blocks"][1]["rows"] == 2
    assert "paths: [sub/*.tsv, sub/a.tsv]" in page  # 'View query' shows the block as written
    absolute = text.replace("'sub/*.tsv'", f"'{tmp_path}/sub/*.tsv'")
    files = cli.build_report("blocks", absolute, tmp_path)[1]["inputs"]["t"]["files"]
    assert [f["path"] for f in files] == [f"{tmp_path}/sub/a.tsv", f"{tmp_path}/sub/b.tsv"]


def test_input_paths_resolve_against_file(tmp_path):
    (tmp_path / "sub").mkdir()
    (tmp_path / "sub" / "a.tsv").write_text("n\n1\n")
    path = tmp_path / "sub" / "spec.md"
    path.write_text(block("input a", "adapter: tsv\npath: a.tsv\n") + block("DataTable t", "source: a\n"))
    manifest = cli.build_report(path)[1]
    [f] = manifest["inputs"]["a"]["files"]
    assert f["path"] == "a.tsv" and f["resolved_path"] == str((tmp_path / "sub" / "a.tsv").resolve())


def test_inline_blocks_and_base_dir(tmp_path):
    (tmp_path / "a.tsv").write_text("n\n1\n2\n")
    text = block("input a", "adapter: tsv\npath: a.tsv\n") + block("DataTable t", "source: a\n")
    fragment, manifest, _ = cli.build_report("blocks", text, tmp_path)
    assert manifest["file"]["path"] == "blocks" and manifest["blocks"][1]["rows"] == 2
    assert '<section class="block datatable">' in fragment
    with pytest.raises(cli.SpecError, match="^blocks:1: input 'a': relative path 'a.tsv' in inline blocks; "
                                            "use absolute paths or pass base_dir$"):
        cli.build_report("blocks", text)
    absolute = text.replace("path: a.tsv", f"path: {tmp_path / 'a.tsv'}")
    assert cli.build_report("blocks", absolute)[1]["blocks"][1]["rows"] == 2


def test_input_below_its_reader(tmp_path):
    body = (block("DataTable t", "source: late\ncolumns: [element_symbol]\nfilters:\n  element_type: STRESS\n")
            + block("input late", f"adapter: amrfinder\npath: {AMRFINDER}\n"))
    page, manifest, _ = cli.build_report(write(tmp_path, body))
    assert [(b["kind"], b["name"], b["rows"]) for b in manifest["blocks"]] == [("DataTable", "t", 2),
                                                                              ("input", "late", 34)]
    assert "<p>From late:</p>\n<pre>adapter: amrfinder\n" in page


def test_sql_consumed_with_filter(tmp_path):
    sql = "select element_type, count(*) as n from amrfinder group by element_type\n"
    body = block("sql counts", sql) + block("DataTable t", "source: counts\nfilters:\n  element_type: AMR\n")
    out = build(spec(tmp_path, body))
    assert '<tr><td>AMR</td><td class="num">13</td></tr>' in out
    assert "<p>From counts (ad hoc SQL):</p>\n<pre>select element_type, count(*) as n" in out
    assert "counts</summary>" not in out  # the sql block renders nothing


def test_sql_reads_query(tmp_path):
    body = (block("query hits", "source: amrfinder\nfilters:\n  element_type: STRESS\n")
            + block("sql symbols", "select element_symbol from hits order by element_symbol\n")
            + block("DataTable t", "source: symbols\n"))
    path = spec(tmp_path, body)
    out = build(path)
    assert rows(path) == {"t": 2}
    assert "From hits" not in out  # ad hoc SQL ends the lineage


def test_chain_lineage(tmp_path):
    body = (block("query hits", "source: amrfinder\nfilters:\n  element_type: [AMR, STRESS]\n")
            + block("DataTable amr", "source: hits\nfilters:\n  element_type: AMR\n"
                    "columns: [element_type, element_symbol, element_class]\n")
            + block("query classes", "source: amr\ngroup_by: element_class\nmeasures:\n  n: count(*)\n")
            + block("DataTable top", "source: classes\nsort: [n, element_class]\n"))
    path = spec(tmp_path, body)
    query = details(build(path), "datatable", 1)
    assert query.startswith("\n<pre>source: classes\nsort: [n, element_class]\n</pre>\n")
    assert re.findall(r"<p>From (\w+)", query) == ["classes", "amr", "hits", "amrfinder"]
    assert "<p>From hits:</p>\n<pre>source: amrfinder\nfilters:\n  element_type: [AMR, STRESS]\n</pre>" in query
    assert rows(path) == {"amr": 13, "top": 5}  # amr is shown, and read by classes


def test_sql_ancestor_in_lineage(tmp_path):
    body = (block("sql amr", "select * from amrfinder where element_type = 'AMR'\n")
            + block("query genes", "source: amr\ncolumns: [element_symbol]\n")
            + block("DataTable t", "source: genes\nsort: element_symbol\n"))
    out = build(spec(tmp_path, body))
    assert re.findall(r"<p>From (\w+)", details(out, "datatable")) == ["genes", "amr"]
    assert "<p>From amr (ad hoc SQL):</p>\n<pre>select * from amrfinder where element_type = &#x27;AMR&#x27;</pre>" in out
    assert "<p>From amrfinder" not in out


def test_matrix_as_source(tmp_path):
    grid = ("source: amrfinder\ngroup_by: [element_type, sample]\nmeasures:\n  n: count(*)\n"
            "rows: element_type\ncols: sample\n")
    body = block("Matrix m", grid) + block("DataTable t", "source: m\nfilters:\n  element_type: STRESS\n")
    path = spec(tmp_path, body)
    out = build(path)
    assert '<tr><td>STRESS</td><td>SAMD00654274</td><td class="num">2</td></tr>' in out
    assert rows(path)["t"] == 1
    assert re.findall(r"<p>From (\w+)", details(out, "datatable")) == ["m", "amrfinder"]


def test_tree_is_not_a_source(tmp_path):
    body = block("input t", "adapter: newick\npath: t.nwk\n") + block("Tree p", "source: t\n")
    (tmp_path / "t.nwk").write_text("(A,B);")
    with pytest.raises(cli.SpecError, match="^spec.md:8: DataTable 'x': source 'p' is a Tree block, which no block can read$"):
        build(write(tmp_path, body + block("DataTable x", "source: p\n")))
    with pytest.raises(cli.SpecError, match="^spec.md:8: Tree 'q': tips 'p' is a Tree block"):
        build(write(tmp_path, body + block("Tree q", "source: t\ntips: p\nmatch: label\n")))



@pytest.mark.parametrize("second, error", [
    (block("DataTable amrfinder", "source: amrfinder\n"), "DataTable 'amrfinder': the name 'amrfinder' is taken by the "
                                                           "block on line 1"),
    (block("input hits", f"adapter: amrfinder\npath: {AMRFINDER}\n"), "input 'hits': the name 'hits' is taken by "
                                                                       "the block on line 6"),
    (block("Tree hits", "source: amrfinder\n"), "Tree 'hits': the name 'hits' is taken by the block on line 6"),
    (block("sql hits", "select 1\n"), "sql 'hits': the name 'hits' is taken"),
    (block("Matrix hits", "source: amrfinder\nrows: sample\ncols: sample\n"), "Matrix 'hits': the name 'hits' is taken"),
])
def test_unique_names(tmp_path, second, error):
    body = block("query hits", "source: amrfinder\n") + second
    with pytest.raises(cli.SpecError, match=f"^spec.md:9: {re.escape(error)}"):
        build(spec(tmp_path, body))


def test_tree_name_is_taken(tmp_path):
    (tmp_path / "t.nwk").write_text("(A,B);")
    body = (block("input t", "adapter: newick\npath: t.nwk\n") + block("Tree p", "source: t\n")
            + block("DataTable p", "source: t\n"))
    with pytest.raises(cli.SpecError, match="^spec.md:8: DataTable 'p': the name 'p' is taken by the block on line 5$"):
        build(write(tmp_path, body))


def test_source_names_no_block(tmp_path):
    body = (block("DataTable t", "source: later\n") + block("query later", "source: amrfinder\n"))
    with pytest.raises(cli.SpecError, match="^spec.md:6: DataTable 't': source 'later' is not above this block; "
                                            "a block reads inputs and the blocks above it$"):
        build(spec(tmp_path, body))
    with pytest.raises(cli.SpecError, match="^spec.md:6: DataTable 't': source 'nope' names no block; "
                                            "readable here: amrfinder$"):
        build(spec(tmp_path, block("DataTable t", "source: nope\n")))


@pytest.mark.parametrize("sql, found", [
    ("select 1; select 2\n", "SELECT, SELECT"),
    ("create table x as select 1\n", "CREATE"),
    (f"copy amrfinder to '{{path}}'\n", "COPY"),
    ("set enable_external_access = true\n", "SET"),
])
def test_sql_not_one_select(tmp_path, sql, found):
    sql = sql.replace("{path}", str(tmp_path / "out.csv"))
    with pytest.raises(ValueError, match=f"sql 'x': the body must be exactly one SELECT statement; found: {found}"):
        build(spec(tmp_path, block("sql x", sql)))


def test_sql_cannot_read_files(tmp_path):
    with pytest.raises(ValueError, match="sql 'x':.*file system operations are disabled"):
        build(spec(tmp_path, block("sql x", f"select * from read_csv('{AMRFINDER}')\n")))


def test_sql_in_manifest(tmp_path):
    body = block("sql a", "select 1 as n\n") + block("sql b", "select n + 1 as n from a\n")
    manifest = cli.build_report(spec(tmp_path, body))[1]
    assert manifest["sql"] == [{"name": "a", "sql": "select 1 as n", "rows": 1},
                               {"name": "b", "sql": "select n + 1 as n from a", "rows": 1}]
    assert [(b["kind"], b["name"]) for b in manifest["blocks"]] == [("input", "amrfinder"), ("sql", "a"), ("sql", "b")]


def test_cohort_tree():
    out = build(COHORT)
    assert details(out, "tree").startswith(
        "\n<p>Tree from source tree: 10 tips; rooted as written; branch lengths as given<br>"
        "Tips from summary: ad hoc SQL; matched on sample; 10 of 10 tips matched<br>"
        "Coloured by sequence_type; labelled by sample, sequence_type</p>\n<pre>source: tree\ntips: summary\n")
    figure = out.split('<svg class="tree"')[1].split("</svg>")[0]
    tips = re.findall(r"<title>(\w+)</title>", figure)
    assert len(tips) == 10
    # ladderised: the ST8 clade, with the untyped SAMD00654278 inside it, is together at the bottom
    assert tips[-4:] == ["SAMD00654306", "SAMD00654296", "SAMD00654278", "SAMD00654274"]
    assert "#" not in figure  # colours come from style.css
    assert '<p class="legend">sequence_type: ' in out
    assert "other values, sharing one grey: 45 <svg" in out and '<circle class="tree-dot null" cx="5" cy="5" r="4"/></svg>null</p>' in out
    assert "<td>tree</td><td>newick preliminary</td><td><code>../data/tree_results/cohort.nwk</code>" in out


def test_cohort_builds_identical():
    strip = lambda s: re.sub(r"Built \S+ with", "Built <t> with", s)
    assert strip(build(COHORT)) == strip(build(COHORT))


def tree_spec(tmp_path, nwk, body):
    """A file with input 't', the newick nwk, then body, which starts on line 5."""
    (tmp_path / "t.nwk").write_text(nwk)
    return write(tmp_path, block("input t", "adapter: newick\npath: t.nwk\n") + body)


def test_tree_without_tips(tmp_path, monkeypatch):
    manifests = []
    monkeypatch.setattr(cli, "render_provenance", lambda m: manifests.append(m) or "")
    out = build(tree_spec(tmp_path, "((A:1,B:2):1,C:1);\n", block("Tree p", "source: t\n")))
    assert "<p>Tree from source t: 3 tips; rooted as written; branch lengths as given<br>Labelled by tip label</p>" in out
    assert "<title>C</title>" in out and "<title>A</title>" in out
    assert '<p class="legend">' not in out
    assert '<text class="scale"' in out
    entry = manifests[0]["trees"][0]
    assert entry["file"] == "t.nwk" and entry["tree"] == 1 and entry["tips"] == 3 and entry["matched"] is None


def test_tree_sql_reaches_tips(tmp_path):
    body = block("sql tips", "select label from t where is_tip order by label\n") + block("DataTable x", "source: tips\n")
    out = build(tree_spec(tmp_path, "((A,B),C);", body))
    assert "<tr><td>A</td></tr>\n<tr><td>B</td></tr>\n<tr><td>C</td></tr>" in out


def test_tree_mismatches(tmp_path, monkeypatch):
    manifests = []
    monkeypatch.setattr(cli, "render_provenance", lambda m: manifests.append(m) or "")
    sql = "select * from (values ('A', 'x'), ('B', 'y'), ('B', 'z'), ('Q', 'x')) v(sample, st)\n"
    body = block("sql meta", sql) + block("Tree p", "source: t\ntips: meta\nmatch: sample\nlabel: st\ncolor: st\n")
    out = build(tree_spec(tmp_path, "((A:1,B:1):1,C:2);", body))
    assert ("Tips from meta: ad hoc SQL; matched on sample; 2 of 3 tips matched; tips with no row: C; "
            "rows with no tip: Q; match values on more than one row, first row used: B<br>") in out
    assert '<text class="unmatched"' in out and ">C</text>" in out
    assert "grey italic, no dot: the tip label of a tip with no row in meta" in out
    entry = manifests[0]["trees"][0]
    assert entry["matched"] == 2 and entry["unmatched_tips"] == ["C"] and entry["unmatched_rows"] == ["Q"]


def test_tree_tips_from_input(tmp_path):
    (tmp_path / "meta.csv").write_text("sample,st\nA,8\nB,5\n")
    body = (block("Tree p", "source: t\ntips: meta\nmatch: sample\ncolor: st\n")
            + block("input meta", "adapter: csv\npath: meta.csv\n"))
    out = build(tree_spec(tmp_path, "((A:1,B:1):1,C:2);", body))
    assert "Tips from meta: adapter csv, file meta.csv; matched on sample; 2 of 3 tips matched; tips with no row: C<br>" in out
    query = details(out, "tree")
    assert "</p>\n<pre>source: t\ntips: meta\nmatch: sample\ncolor: st\n</pre>\n" in query
    assert re.findall(r"<p>From (\w+)", query) == ["t", "meta"]
    assert "<p>From meta:</p>\n<pre>adapter: csv\npath: meta.csv\n</pre>" in query


def test_tree_tips_lineage(tmp_path):
    (tmp_path / "meta.csv").write_text("sample,st\nA,8\nB,5\n")
    body = (block("input meta", "adapter: csv\npath: meta.csv\n")
            + block("query st8", "source: meta\nfilters:\n  st: 8\n")
            + block("Tree p", "source: t\ntips: st8\nmatch: sample\nlabel: [sample, st]\n"))
    out = build(tree_spec(tmp_path, "(A,B);", body))
    assert "Tips from st8: filters st = 8; matched on sample; 1 of 2 tips matched" in out
    assert re.findall(r"<p>From (\w+)", details(out, "tree")) == ["t", "st8", "meta"]


def test_tree_cladogram_and_notes(tmp_path):
    out = build(tree_spec(tmp_path, "((A,B),(C,A));", block("Tree p", "source: t\n")))
    assert ("Tree from source t: 4 tips; rooted as written; no branch lengths, drawn as a cladogram with tips "
            "aligned; 1 tip label on more than one tip<br>") in out
    assert '<text class="scale"' not in out
    out = build(tree_spec(tmp_path, "((A:-1,B:1):1,C);", block("Tree p", "source: t\n")))
    assert "branch lengths as given; 1 negative length drawn as zero; 1 missing length drawn as zero<br>" in out


def test_tree_several_trees(tmp_path):
    with pytest.raises(ValueError, match=r"Tree 'p': source 't' holds 2 trees \(t.nwk tree 1, t.nwk tree 2\); "
                                         r"pick one with 'tree'"):
        build(tree_spec(tmp_path, "(A,B);\n(C,D);\n", block("Tree p", "source: t\n")))
    out = build(tree_spec(tmp_path, "(A,B);\n(C,D);\n", block("Tree p", "source: t\ntree: 2\n")))
    assert "Tree from source t (t.nwk, tree 2): 2 tips" in out and "<title>D</title>" in out
    with pytest.raises(ValueError, match="holds 2 trees .*; 0 match tree 3; pick one"):
        build(tree_spec(tmp_path, "(A,B);\n(C,D);\n", block("Tree p", "source: t\ntree: 3\n")))


def test_tree_needs_newick_source(tmp_path):
    with pytest.raises(ValueError, match="^spec.md:6: Tree 'p': source 'amrfinder' is loaded by adapter 'amrfinder'; "
                                         "a Tree block needs an input loaded by the newick adapter$"):
        build(spec(tmp_path, block("Tree p", "source: amrfinder\n")))
    body = block("sql s", "select * from t\n") + block("Tree p", "source: s\n")
    with pytest.raises(ValueError, match="^spec.md:8: Tree 'p': source 's' is a sql block; "
                                         "a Tree block needs an input loaded by the newick adapter$"):
        build(tree_spec(tmp_path, "(A,B);", body))


@pytest.mark.parametrize("text, error", [
    ("source: t\nshape: circular\n", "unknown key 'shape'; accepted: source, tips, match, label, color, file, tree"),
    ("tips: meta\nmatch: sample\n", "needs 'source'"),
    ("source: t\nlabel: sample\n", "'label' need 'tips'"),
    ("source: t\ntips: meta\n", "'tips' names a block and needs 'match'"),
    ("source: t\ntree: 0\n", "'tree' must be a whole number from 1"),
])
def test_tree_block_shape(tmp_path, text, error):
    with pytest.raises(ValueError, match=f"Tree 'p': .*{re.escape(error)}"):
        build(tree_spec(tmp_path, "(A,B);", block("Tree p", text)))


def test_tree_column_not_in_tips(tmp_path):
    body = block("sql meta", "select 'A' as sample\n") + block("Tree p", "source: t\ntips: meta\nmatch: sample\ncolor: st\n")
    with pytest.raises(ValueError, match="Tree 'p': 'st' not a column of meta; columns: sample"):
        build(tree_spec(tmp_path, "(A,B);", body))


def test_tree_parse_error_names_source(tmp_path):
    with pytest.raises(ValueError, match="^spec.md:1: input 't': t.nwk, character 4: unbalanced"):
        build(tree_spec(tmp_path, "(A,B", ""))


def test_tree_label_escaped(tmp_path):
    out = build(tree_spec(tmp_path, "('<b>x</b>&',B);", block("Tree p", "source: t\n")))
    assert "<b>x</b>" not in out
    assert "&lt;b&gt;x&lt;/b&gt;&amp;" in out


def test_tree_deep(tmp_path):
    depth = 3000
    out = build(tree_spec(tmp_path, "(" * depth + "A:1" + ",B:1):1" * depth + ";", block("Tree p", "source: t\n")))
    assert f"Tree from source t: {depth + 1} tips" in out


def test_error_names_spec_line(tmp_path):
    body = "prose\n\n" + block("DataTable ok", "source: amrfinder\n") + block("Matrix m", "source: amrfinder\n")
    with pytest.raises(ValueError, match=r"^spec.md:11: Matrix 'm': needs 'rows'"):  # parse time
        build(spec(tmp_path, body))
    body = "prose\n\n" + block("DataTable ok", "source: amrfinder\n") + block("query d", "source: nope\n")
    with pytest.raises(ValueError, match=r"^spec.md:11: query 'd': source 'nope' names no block"):  # parse time
        build(spec(tmp_path, body))
    body = "prose\n\n" + block("DataTable ok", "source: amrfinder\n") + block("query d", "source: ok\nsort: nope\n")
    with pytest.raises(ValueError, match=r"^spec.md:11: query 'd': Binder Error"):  # build time
        build(spec(tmp_path, body))


def test_block_kind_in_messages(tmp_path):
    text = "source: amrfinder\ncolumns: [sample]\ngroup_by: [sample]\nrows: sample\ncols: sample\n"
    with pytest.raises(ValueError, match="Matrix 'm': give 'columns' or 'group_by', not both"):
        build(spec(tmp_path, block("Matrix m", text)))
    with pytest.raises(ValueError, match="query 'hits': needs 'source'"):
        build(spec(tmp_path, block("query hits", "columns: [sample]\n")))


def test_unknown_query_key(tmp_path):
    with pytest.raises(ValueError, match=r"^spec.md:6: DataTable 't': unknown key 'filter'; accepted: source, "
                                         "columns, filters, group_by, measures, sort, display, spread$"):
        build(spec(tmp_path, block("DataTable t", "source: amrfinder\nfilter:\n  sample: x\n")))
    with pytest.raises(ValueError, match=r"^spec.md:6: query 't': unknown key 'result'; accepted: source, "):
        build(spec(tmp_path, block("query t", "result: amrfinder\n")))
    text = "source: amrfinder\ngroup_by: [sample, element_symbol]\nrows: sample\ncols: element_symbol\nlimit: 3\n"
    with pytest.raises(ValueError, match=r"Matrix 'm': unknown key 'limit'; accepted: source, columns, "
                                         "filters, group_by, measures, sort, rows, cols, value, scale, digits, mark, pick, symmetric$"):
        build(spec(tmp_path, block("Matrix m", text)))


def test_scalar_sort_and_columns(tmp_path):
    body = (block("DataTable t", "source: amrfinder\ncolumns: element_symbol\nsort: element_symbol\n")
            + block("DataTable g", "source: amrfinder\ngroup_by: element_type\nsort: element_type\n"))
    page, manifest, _ = cli.build_report(spec(tmp_path, body))
    assert [b["output_columns"] for b in manifest["blocks"][1:]] == [["element_symbol"], ["element_type"]]


@pytest.mark.parametrize("text, error", [
    ("columns: [sample, nope]\n", 'DataTable \'t\': Binder Error: Referenced column "nope" not found'),
    ("filters:\n  nope: x\n", 'DataTable \'t\': Binder Error: Referenced column "nope" not found'),
])
def test_query_duckdb_error(tmp_path, text, error):
    with pytest.raises(cli.SpecError, match=r"^spec.md:6: " + re.escape(error)):
        build(spec(tmp_path, block("DataTable t", "source: amrfinder\n" + text)))


def test_block_yaml_error(tmp_path):
    with pytest.raises(ValueError, match=r"^spec.md:6: DataTable 'x': invalid YAML on line 9: "):
        build(spec(tmp_path, block("DataTable x", "source: amrfinder\nfilters: [a\n")))


def test_cannot_read_file(tmp_path):
    with pytest.raises(cli.SpecError, match="^cannot read '.*nope.md': No such file or directory$"):
        build(tmp_path / "nope.md")


def test_page_title_and_no_heading(tmp_path):
    path = tmp_path / "my report.md"
    path.write_text("Some prose.\n\n" + INPUT + block("DataTable t", "source: amrfinder\ncolumns: sample\n"))
    out = build(path)
    assert "<title>my report</title>" in out
    assert "<h1" not in out and "<main>\n<p>Some prose.</p>\n" in out


def test_block_body_not_a_mapping(tmp_path):
    with pytest.raises(cli.SpecError, match="spec.md:6: Matrix 'm': the body must be YAML keys and values"):
        build(spec(tmp_path, block("Matrix m", "- source\n")))


@pytest.mark.parametrize("fence, problem", [
    ("```DataTable\nsource: amrfinder\n```\n", "'```DataTable' is not read as a block: expected ```DataTable <name>, "
                                               "one space, the name letters, digits and underscores"),
    ("```DataTable  two\nsource: amrfinder\n```\n", "'```DataTable  two' is not read as a block: expected"),
    ("```query my-hits\nsource: amrfinder\n```\n", "'```query my-hits' is not read as a block: expected"),
    ("```sql\nselect 1\n```\n", "'```sql' is not read as a block: expected ```sql <name>"),
    ("```input\nadapter: tsv\n```\n", "'```input' is not read as a block: expected ```input <name>"),
    ("```DataTable open\nsource: amrfinder\n", "'```DataTable open' is never closed"),
])
def test_unparsed_fence_warns(tmp_path, fence, problem):
    page, manifest, _ = cli.build_report(spec(tmp_path, "prose\n\n" + fence))
    [warning] = manifest["warnings"]
    assert warning.startswith(f"spec.md:8: fence {problem}")
    assert warning.endswith("; rendered as a code block")
    assert [b["kind"] for b in manifest["blocks"]] == ["input"] and "rendered as a code block" not in page


def test_no_warnings(tmp_path):
    assert cli.build_report(COHORT)[1]["warnings"] == []
    assert cli.build_report(SPEC)[1]["warnings"] == []
    body = "```python\nx = 1\n```\n" + block("sql x", "select '\n```DataTable\n' as s\n")  # inside a parsed block
    assert cli.build_report(spec(tmp_path, body))[1]["warnings"] == []


def test_empty_result_warns(tmp_path):
    (tmp_path / "e.csv").write_text("a,b\n")
    body = (block("DataTable t", "source: amrfinder\nfilters:\n  element_type: AMRR\n")
            + block("sql s", "select 1 as n where false\n") + block("query u", "source: amrfinder\n")
            + block("input e", "adapter: csv\npath: e.csv\n"))
    page, manifest, _ = cli.build_report(spec(tmp_path, body))
    assert manifest["warnings"] == ["spec.md:6: DataTable 't': 0 rows; check the filter values: element_type = AMRR",
                                    "spec.md:11: sql 's': 0 rows", "spec.md:17: input 'e': 0 rows"]
    assert "check the filter values" not in page


def run_main(monkeypatch, *args):
    monkeypatch.setattr("sys.argv", ["groundie", "build", *map(str, args)])
    try:
        cli.main()
    except SystemExit as e:
        return e.code
    return 0


def test_main_success(tmp_path, monkeypatch, capsys):
    out, manifest = tmp_path / "r.html", tmp_path / "m.json"
    assert run_main(monkeypatch, COHORT, "-o", out, "--manifest", manifest) == 0
    stdout, stderr = capsys.readouterr()
    assert stderr == ""
    lines = stdout.splitlines()
    assert lines[0] == f"built {out} from {COHORT}"
    assert lines[1] == "cohort.md:3: input mlst: 10 rows; columns: sample, scheme, sequence_type, calls"
    assert lines[2].startswith("cohort.md:8: input amrfinder: 334 rows; columns: sample, protein_id, contig_id, ")
    assert lines[3].startswith("cohort.md:13: input tree: 18 rows; columns: file, tree, node, parent, label, branch_length, "
                               "support, is_tip, comment")
    assert lines[4:] == [
        "cohort.md:22: MlstTable st: 10 rows in 1 table; columns: sample, scheme, sequence_type, calls",
        "cohort.md:31: sql summary: 10 rows; columns: sample, sequence_type, amr_hits, virulence_hits, stress_hits",
        "cohort.md:44: DataTable summary_table: 10 rows; columns: sample, sequence_type, amr_hits, virulence_hits, "
        "stress_hits",
        "cohort.md:58: Tree phylogeny: 10 tips, 10 matched",
        "cohort.md:70: Matrix resistome: 129 rows; columns: sample, element_class, element_symbol, best_identity, "
        "lowest_coverage, copies",
    ]
    assert out.read_text().startswith("<!doctype html>")
    data = json.loads(manifest.read_text())
    assert data["warnings"] == []
    assert data["blocks"][3]["tables"] == 1
    assert data["blocks"][6] == {"line": 58, "kind": "Tree", "name": "phylogeny", "tips": 10, "matched": 10}
    assert data["blocks"][7]["output_columns"][-1] == "copies"
    assert data["queries"][-1]["layout"]["scale"] == [90, 100]
    assert [q["kind"] for q in data["queries"]] == ["MlstTable", "DataTable", "Matrix"]
    assert data["queries"][0]["tables"] == [{"caption": "saureus", "rows": 10, "columns": [
        "sample", "sequence_type", "arcC", "aroE", "glpF", "gmk", "pta", "tpi", "yqiL"]}]
    assert list(data["inputs"]) == ["mlst", "amrfinder", "tree"]
    assert data["inputs"]["tree"]["files"][0]["path"] == "../data/tree_results/cohort.nwk"


def test_main_warnings(tmp_path, monkeypatch, capsys):
    path = spec(tmp_path, "```DataTable\nx\n```\n" + block("DataTable t", "source: amrfinder\nfilters:\n  sample: nope\n"))
    assert run_main(monkeypatch, path, "-o", tmp_path / "r.html") == 0
    stdout, stderr = capsys.readouterr()
    assert stderr.splitlines() == [
        "groundie: warning: spec.md:6: fence '```DataTable' is not read as a block: expected ```DataTable <name>, "
        "one space, the name letters, digits and underscores; rendered as a code block",
        "groundie: warning: spec.md:9: DataTable 't': 0 rows; check the filter values: sample = nope",
    ]
    assert "spec.md:9: DataTable t: 0 rows; columns: sample, " in stdout


@pytest.mark.parametrize("body, error", [
    (block("Matrix m", "source: amrfinder\n"), "spec.md:6: Matrix 'm': needs 'rows'"),  # parse time
    (block("sql q", "select * from nope\n"), "spec.md:6: sql 'q': Catalog Error: Table with name nope does not exist!"),
    (block("DataTable t", "source: amrfinder\ncolumns: [sample, nope]\n"),
     "spec.md:6: DataTable 't': Binder Error: Referenced column \"nope\" not found in FROM clause!"),
    (block("query d", "source: amrfinder\nsort: [nope]\n"),
     "spec.md:6: query 'd': Binder Error: Referenced column \"nope\" not found in FROM clause!"),
    (block("input x", "adapter: tsv\npath: nope.tsv\n"), "spec.md:6: input 'x': 'nope.tsv' is not a file"),
])
def test_main_error(tmp_path, monkeypatch, capsys, body, error):
    out = tmp_path / "r.html"
    assert run_main(monkeypatch, spec(tmp_path, body), "-o", out) == 1
    stdout, stderr = capsys.readouterr()
    assert stdout == "" and stderr.count("\n") == 1
    assert stderr.startswith(f"groundie: error: {error}")
    assert not out.exists()


@pytest.mark.parametrize("args", [["nope.md"], ["spec.md"]])
def test_main_error_without_traceback(tmp_path, monkeypatch, capsys, args):
    (tmp_path / "spec.md").write_text("```DataTable t\nsource: [t\n```\n")
    assert run_main(monkeypatch, *(tmp_path / a for a in args), "-o", tmp_path / "r.html") == 1
    stderr = capsys.readouterr().err
    assert stderr.startswith("groundie: error: ") and stderr.count("\n") == 1 and "Traceback" not in stderr


@pytest.mark.parametrize("flag", ["-o", "--manifest"])
def test_main_cannot_write(tmp_path, monkeypatch, capsys, flag):
    args = {"-o": tmp_path / "r.html", flag: tmp_path / "nodir" / "out"}
    assert run_main(monkeypatch, spec(tmp_path, ""), *(a for pair in args.items() for a in pair)) == 1
    stdout, stderr = capsys.readouterr()
    assert stdout == "" and stderr.count("\n") == 1
    assert re.match(r"groundie: error: cannot write '.*nodir/out': No such file or directory", stderr)


def test_main_lets_a_bug_through(tmp_path, monkeypatch):
    def bug(path):
        raise ValueError("not about the spec")
    monkeypatch.setattr(cli, "build_report", bug)
    with pytest.raises(ValueError, match="not about the spec"):
        run_main(monkeypatch, spec(tmp_path, ""), "-o", tmp_path / "r.html")


@pytest.mark.parametrize("path", [SPEC, COHORT])
def test_manifest_is_json(path):
    manifest = cli.build_report(path)[1]
    assert json.loads(json.dumps(manifest)) == manifest


def test_mash_adapter(tmp_path):
    (tmp_path / "d.tsv").write_text("genomes/A.fa.gz\tsketches/B.msh\t0.001\t0\t9000/10000\n")
    text = block("input d", "adapter: mash\npath: d.tsv\n") + block("DataTable t", "source: d\n")
    page, manifest, con = cli.build_report(write(tmp_path, text))
    assert con.execute("select * from d").fetchall() == [("A", "B", 0.001, 0.0, 9000, 10000)]


def test_mlst_adapter(tmp_path):
    # Real mlst 2.23 output: `mlst genomes/SAMD00654274.fa.gz rerun/SAMD00654274.fa.gz genomes/unknown.fa`,
    # where the rerun is another genome under the same name and unknown.fa is random sequence.
    text = block("input m", f"adapter: mlst\npaths: [{MLST}, {FIXTURES / 'mlst_edge.tsv'}]\n") + block("DataTable t", "source: m\n")
    page, manifest, con = cli.build_report(write(tmp_path, text))
    got = con.execute("select * from m").fetchall()
    assert [f["rows"] for f in manifest["inputs"]["m"]["files"]] == [10, 3]
    assert manifest["blocks"][0]["output_columns"] == ["sample", "scheme", "sequence_type", "calls"]
    calls = lambda tpi: [{"locus": locus, "call": call} for locus, call in zip(
        ["arcC", "aroE", "glpF", "gmk", "pta", "tpi", "yqiL"], ["3", "3", "1", "1", "4", tpi, "3"])]
    assert got[0] == got[10] == ("SAMD00654274", "saureus", 8, calls("4"))
    assert got[11] == ("SAMD00654274", "saureus", None, calls("~4"))  # same name, its own call; '-' ST is null
    assert got[12] == ("unknown", None, None, [])  # no scheme found: kept, with no calls


@pytest.mark.parametrize("adapter, paths", [
    (mlst.Mlst, [MLST, FIXTURES / "mlst_edge.tsv"]),
    (amrfinder.Amrfinder, [AMRFINDER]),
])
def test_rows_match_schema(adapter, paths):
    from jsonschema import Draft202012Validator
    Draft202012Validator.check_schema(adapter.SCHEMA)
    names = [n for n, _ in adapter.COLUMNS]
    for row in [row for path in paths for row in adapter.rows(path, {})]:
        Draft202012Validator(adapter.SCHEMA).validate(dict(zip(names, row)))


def test_adapter_titles():
    assert mlst.Mlst.TITLES == {"sample": "Sample", "scheme": "Scheme", "sequence_type": "ST"}
    properties = mlst.Mlst.SCHEMA["properties"]
    assert properties["sequence_type"]["title"] == "ST"
    assert properties["sequence_type"]["description"] == "The ST; null when none was assigned."
    assert "title" not in properties["calls"]
    assert mlst.Mlst.DESCRIPTIONS["calls"] == properties["calls"]["description"]  # the docstrings, as in the schema
    assert set(mlst.Mlst.DESCRIPTIONS) == set(properties)
    titles = amrfinder.Amrfinder.TITLES  # by default, the header's name
    assert titles["sample"] == "Sample" and titles["pct_identity"] == "% Identity to reference"
    assert len(titles) == len(amrfinder.Amrfinder.COLUMNS)


def test_amrfinder_rows(tmp_path):
    hit = ("SAMD00654274", None, "SAMD00654274.contig00001", 72659, 75058, "-", "parC_S80Y",
           "Staphylococcus aureus quinolone resistant ParC", "core", "AMR", "POINT", "QUINOLONE", "QUINOLONE", "POINTX",
           800, 800, 100.0, 99.5, 800, "WP_001289586.1", "DNA topoisomerase IV subunit A ParC", None, None, None)
    assert amrfinder.Amrfinder.rows(AMRFINDER, {})[0] == hit
    # the same file as if run without --name and with --print_node
    lines = [line.split("\t", 1)[1] for line in AMRFINDER.read_text().splitlines()[:2]]
    (tmp_path / "KP01.tsv").write_text(f"{lines[0]}\tHierarchy node\n{lines[1]}\tparC_S80Y\n")
    assert amrfinder.Amrfinder.rows(tmp_path / "KP01.tsv", {}) == [("KP01", *hit[1:-1], "parC_S80Y")]
    # a byte order mark, or a header in another case, still finds Name rather than falling back to the file's
    lines = AMRFINDER.read_text().splitlines()[:2]
    (tmp_path / "KP02.tsv").write_text("﻿" + "\n".join(lines) + "\n", encoding="utf-8")
    (tmp_path / "KP03.tsv").write_text(lines[0].lower() + "\n" + lines[1] + "\n")
    for name in ("KP02.tsv", "KP03.tsv"):
        assert amrfinder.Amrfinder.rows(tmp_path / name, {}) == [hit]


def test_amrfinder_missing_column(tmp_path):
    (tmp_path / "a.tsv").write_text("Name\tProtein id\n")
    with pytest.raises(ValueError, match=r"^a\.tsv: no column 'Contig id'; the header has 'Name', 'Protein id'$"):
        amrfinder.Amrfinder.rows(tmp_path / "a.tsv", {})
    (tmp_path / "a.tsv").write_text("Protein identifier\tContig id\tStart\n")  # version 3, not read
    with pytest.raises(ValueError, match=r"^a\.tsv: no column 'Protein id'; the header has 'Protein identifier', "):
        amrfinder.Amrfinder.rows(tmp_path / "a.tsv", {})
    (tmp_path / "a.tsv").write_text("")
    assert amrfinder.Amrfinder.rows(tmp_path / "a.tsv", {}) == []


def test_amrfinder_sample_from_file_names(tmp_path):
    header, hit = [line.split("\t", 1)[1] for line in AMRFINDER.read_text().splitlines()[:2]]  # without --name
    for sample in ("KP01", "KP02"):
        (tmp_path / sample).mkdir()
        (tmp_path / sample / "amrfinder.tsv").write_text(f"{header}\n{hit}\n")
        (tmp_path / f"{sample}.tsv").write_text(f"{header}\n{hit}\n")
    text = block("input a", "adapter: amrfinder\npath: '*.tsv'\n") + block("DataTable t", "source: a\n")
    assert cli.build_report(write(tmp_path, text))[1]["blocks"][1]["rows"] == 2  # KP01 and KP02
    first, second = (tmp_path / s / "amrfinder.tsv" for s in ("KP01", "KP02"))
    error = (f"input 'a': {first} and {second} both lack column 'Name', and their names both give sample "
             f"'amrfinder'")
    with pytest.raises(cli.SpecError, match=re.escape(error) + "$"):
        build(write(tmp_path, text.replace("'*.tsv'", "'*/amrfinder.tsv'")))
    (tmp_path / "KP02" / "amrfinder.tsv").write_text(AMRFINDER.read_text())  # one with --name: no clash
    build(write(tmp_path, text.replace("'*.tsv'", "'*/amrfinder.tsv'")))


def test_mlst_rows(tmp_path):
    path = tmp_path / "a.tsv"
    # any number of loci, mlst's flags kept in the call, Windows line ends, blank lines skipped
    path.write_text("x/y.fasta.gz\tbig\t1" + "".join(f"\tl{i}({i})" for i in range(1, 9))
                    + "\r\n\r\nz.fna\tsaureus\t-\tarcC(-)\taroE(~4)\tglpF(4?)\tgmk(4,5)\r\n")
    assert mlst.Mlst.rows(path, {}) == [
        ("y", "big", 1, [{"locus": f"l{i}", "call": str(i)} for i in range(1, 9)]),
        ("z", "saureus", None, [{"locus": "arcC", "call": None}, {"locus": "aroE", "call": "~4"},
                                {"locus": "glpF", "call": "4?"}, {"locus": "gmk", "call": "4,5"}])]


@pytest.mark.parametrize("line, error", [
    ("x.fa\tsaureus", "a.tsv, line 2: has 2 columns; needs at least 3"),
    ("-\tsaureus\t8", "a.tsv, line 2: sample: '-' is null, which this column can't be"),
    ("x.fa\tsaureus\tST8\tarcC(3)", "a.tsv, line 2: sequence_type: 'ST8' is not an integer"),
    ("x.fa\tsaureus\t8\tarcC3", "a.tsv, line 2: calls: 'arcC3' is not 'locus(call)'"),
    ("x.fa\tsaureus\t8\t-(3)", "a.tsv, line 2: calls: locus: '-' is null, which this column can't be"),
])
def test_mlst_bad_line(tmp_path, line, error):
    (tmp_path / "a.tsv").write_text(f"y.fa\t-\t-\n{line}\n")
    with pytest.raises(cli.SpecError, match=f"^spec.md:1: input 'm': {re.escape(error)}$"):
        build(write(tmp_path, block("input m", "adapter: mlst\npath: a.tsv\n")))


class Part(utils.Struct):
    PATTERN = "{name}:{n}"

    name: str
    n: int


@pytest.mark.parametrize("hint, source, problem", [
    (int, None, "needs a source: column(...) or columns(...)"),
    (str, utils.column("Name"), "names column 'Name', but the format has no header"),
    (str, utils.column(0), "columns are counted from 1, first to last"),
    (list[str], utils.columns(2, 0), "columns are counted from 1, first to last"),
    (list[str], utils.columns(3, 2), "columns are counted from 1, first to last"),
    (list[str], utils.columns(2, ...), "columns(...) ends at a position: -1 for each line's last"),
    (str, utils.columns(2, -1), "columns(...) reads a list: its type is list[...]"),
    (list[str], utils.column(2), "a list is read from columns(...)"),
    (bool, utils.column(2), "no reader for <class 'bool'>"),
    (str, utils.column(2, title=""), "a title is text, the header a datatable shows"),
    (str, utils.column(2, title=7), "a title is text, the header a datatable shows"),
    (str, utils.column(2, absent=utils.no_value), "absent is for a column named in the header"),
])
def test_adapter_checks(hint, source, problem):
    assert utils.source_problem(utils.Tsv(header=False), utils.Field(hint, source, None)) == problem
    assert (utils.source_problem(utils.Tsv(header=True), utils.Field(str, utils.column("N", absent=utils.no_value), None))
            == "absent=no_value is null, which only a type with '| None' allows")
    with pytest.raises(TypeError, match=r"^Bad.n: needs a source: column\(\.\.\.\) or columns\(\.\.\.\)$"):
        class Bad(utils.Adapter):
            FORMAT = utils.Tsv(header=False)
            n: int


def test_struct_checks():
    with pytest.raises(TypeError, match=r"^NoPattern: needs a PATTERN naming its parts, such as '\{a\}:\{b\}'$"):
        class NoPattern(utils.Struct):
            a: str
    with pytest.raises(TypeError, match=r"^Mismatch: its PATTERN has parts name, count; its fields are name, n$"):
        class Mismatch(utils.Struct):
            PATTERN = "{name}:{count}"
            name: str
            n: int
    with pytest.raises(TypeError, match=r"^Nested.p: a part is str, int or float, not .*Part"):
        class Nested(utils.Struct):
            PATTERN = "{p}"
            p: Part


class Hits(utils.Adapter):
    FORMAT = utils.Tsv(header=True, null="NA")

    sample: str = utils.column("Name")
    identity: float | None = utils.column("% Identity")
    span: Part = utils.column("Span")


def test_adapter_with_header(tmp_path):
    (tmp_path / "h.tsv").write_text("Span\tName\t% Identity\nc1:7\tS1\t99.5\nc2:-3\tS2\tNA\n")
    assert Hits.rows(tmp_path / "h.tsv", {}) == [
        ("S1", 99.5, {"name": "c1", "n": 7}), ("S2", None, {"name": "c2", "n": -3})]
    (tmp_path / "h.tsv").write_text("Span\tName\nc1:7\tS1\n")
    with pytest.raises(ValueError, match=r"^h\.tsv: no column '% Identity'; the header has 'Span', 'Name'$"):
        Hits.rows(tmp_path / "h.tsv", {})


def test_readers():
    assert utils.Pattern("{locus}({call})").split("tpi(~4)") == {"locus": "tpi", "call": "~4"}
    assert utils.Pattern("{locus}({call})").split("gmk()") == {"locus": "gmk", "call": ""}
    with pytest.raises(ValueError, match=r"^'tpi' is not 'locus\(call\)'$"):
        utils.Pattern("{locus}({call})").split("tpi")
    with pytest.raises(ValueError, match=r"^'8\.0' is not an integer$"):
        utils.integer("8.0")


@pytest.mark.parametrize("file", ["genomes/A.fa.gz", "A.fasta", "a/b/A.fna", "A.fas.gz", "A", "A.fa.bak"])
def test_mlst_sample_names_match_mash(file):
    sql = mash.SAMPLE.format(col="?")
    assert utils.extract_sample_name(file) == duckdb.execute(f"select {sql}", [file]).fetchone()[0]


def test_clusters(tmp_path):
    text = (block("input mash", f"adapter: mash\npath: {MASH}\n") + block("clusters c", "source: mash\nthreshold: 0.005\n")
            + block("DataTable t", "source: c\n"))
    page, manifest, con = cli.build_report(write(tmp_path, text))
    got = {}
    for sample, cluster, n in con.execute("select * from c").fetchall():
        got.setdefault(cluster, []).append(sample)
        assert n == len([1 for s, c, _ in con.execute("select * from c").fetchall() if c == cluster])
    assert got == {"C1": ["SAMD00654276", "SAMD00654277", "SAMD00654301", "SAMD00654312", "SAMD00654356"],
                   "C2": ["SAMD00654274", "SAMD00654278", "SAMD00654296", "SAMD00654306"],
                   "C3": ["SAMD00654319"]}  # largest first; the singleton is kept
    assert "From c:</p>\n<pre>source: mash\npairs: [reference, query]\ndistance: distance\nthreshold: 0.005\nlinkage: single\n" in page
    # the threshold is inclusive; below every distance, each sample is its own cluster, named in order
    (tmp_path / "p.tsv").write_text("x\ty\td\nb\ta\t0.1\nc\ta\t0.2\nc\tb\t0.1\n")
    pairs = block("input p", "adapter: tsv\npath: p.tsv\n")
    con = cli.build_report(write(tmp_path, pairs + block("clusters c", "source: p\npairs: [x, y]\ndistance: d\nthreshold: 0.1\n")))[2]
    assert con.execute("select * from c").fetchall() == [("a", "C1", 3), ("b", "C1", 3), ("c", "C1", 3)]
    con = cli.build_report(write(tmp_path, pairs + block("clusters c", "source: p\npairs: [x, y]\ndistance: d\nthreshold: 0\n")))[2]
    assert con.execute("select * from c").fetchall() == [("a", "C1", 1), ("b", "C2", 1), ("c", "C3", 1)]


@pytest.mark.parametrize("body, error", [
    ("source: p\n", "needs 'threshold', a distance of 0 or more"),
    ("source: p\nthreshold: -1\n", "needs 'threshold', a distance of 0 or more"),
    ("source: p\nthreshold: 1\nlinkage: average\n", "'linkage' must be single"),
    ("source: p\nthreshold: 1\npairs: [x]\n", "'pairs' must be the two columns"),
    ("source: p\nthreshold: 1\nsort: [x]\n", "unknown key 'sort'; accepted: source, threshold, pairs, distance, linkage"),
    ("source: p\nthreshold: 1\n", "'reference', 'query', 'distance' not a column of 'p'; its columns: file, x, y, d"),
])
def test_clusters_shape(tmp_path, body, error):
    (tmp_path / "p.tsv").write_text("x\ty\td\na\tb\t0.1\n")
    text = block("input p", "adapter: tsv\npath: p.tsv\n") + block("clusters c", body)
    with pytest.raises(cli.SpecError, match=re.escape(f"clusters 'c': {error}")):
        build(write(tmp_path, text))


def test_checks(tmp_path):
    (tmp_path / "q.tsv").write_text("sample\tlen\tN50\tgc\n"
                                    "a\t100\t50\t33\n"   # every bound met exactly: min and max are inclusive
                                    "b\t99\t10\t33\n"    # len below min; N50 at 'above', which is exclusive
                                    "c\t101\t\t32\n")    # len above max; N50 missing; gc at 'below'
    text = (block("input q", "adapter: tsv\npath: q.tsv\n")
            + block("checks qc", "source: q\nchecks:\n  len: {min: 100, max: 100}\n  N50: {above: 10}\n"
                    "  gc: {below: 33.5, above: 32}\n")
            + block("DataTable t", "source: qc\n"))
    page, manifest, con = cli.build_report(write(tmp_path, text))
    assert [d[0] for d in con.execute("select * from qc").description] == ["sample", "len", "N50", "gc", "qc", "qc_failed"]
    assert con.execute("select sample, qc, qc_failed from qc").fetchall() == [
        ("a", "PASS", None), ("b", "FAIL", "len 99; N50 10"), ("c", "FAIL", "len 101; N50 missing; gc 32")]
    assert "From qc:</p>\n<pre>source: q\nkey: sample\nchecks:\n  len:\n    min: 100\n" in page
    assert manifest["queries"][0]["kind"] == "checks" and "coalesce" in manifest["queries"][0]["sql"]


@pytest.mark.parametrize("body, error", [
    ("checks:\n  x: {min: 1}\n", "needs 'source'"),
    ("source: p\n", "needs 'checks'"),
    ("source: p\nchecks:\n  x: {min: 1}\nsort: [x]\n", "unknown key 'sort'; accepted: source, key, checks"),
    ("source: p\nchecks:\n  x: {least: 1}\n", "check 'x': unknown bound 'least'; accepted: min, max, above, below"),
    ("source: p\nchecks:\n  x: {}\n", "check 'x' must be {min, max, above or below: <number>}"),
    ("source: p\nchecks:\n  x: 5\n", "check 'x' must be {min, max, above or below: <number>}"),
    ("source: p\nchecks:\n  x: {min: one}\n", "check 'x': bounds must be numbers"),
    ("source: p\nchecks:\n  x: {min: 1, above: 0}\n", "check 'x': give min or above, not both"),
    ("source: p\nchecks:\n  x: {max: 1, below: 2}\n", "check 'x': give min or above, not both, and max or below"),
    ("source: p\nchecks:\n  sample: {min: 1}\n", "check 'sample' is on the key column"),
    ("source: p\nkey: id\nchecks:\n  z: {min: 1}\n", "'id', 'z' not a column of 'p'; its columns: file, sample, x"),
])
def test_checks_shape(tmp_path, body, error):
    (tmp_path / "p.tsv").write_text("sample\tx\na\t1\n")
    text = block("input p", "adapter: tsv\npath: p.tsv\n") + block("checks c", body)
    with pytest.raises(cli.SpecError, match=re.escape(f"checks 'c': {error}")):
        build(write(tmp_path, text))


def test_calls(tmp_path):
    (tmp_path / "s.tsv").write_text("sample\nSAMD00654274\nnohits\n")
    calls = "  MRSA: {any: [mecA, mecC]}\n  PVL: {all: [lukF-PV, lukS-PV]}\n  both: {all: [mecA, mecC]}\n"
    text = (block("input amrfinder", f"adapter: amrfinder\npath: {AMRFINDER}\n") + block("input s", "adapter: tsv\npath: s.tsv\n")
            + block("calls k", f"source: amrfinder\nsamples: s\ncalls:\n{calls}") + block("DataTable t", "source: k\n"))
    page, manifest, con = cli.build_report(write(tmp_path, text))
    assert [d[0] for d in con.execute("select * from k").description] == ["sample", "MRSA", "PVL", "both"]
    # a sample with no hits is kept, false throughout; 'all' needs every value
    assert con.execute("select * from k").fetchall() == [("SAMD00654274", True, True, False), ("nohits", False, False, False)]
    assert "From k:</p>\n<pre>source: amrfinder\nsamples: s\ncolumn: element_symbol\nkey: sample\ncalls:\n  MRSA:\n    any: [mecA, mecC]\n" in page
    # without samples, the samples in source; YAML-typed values compare as text
    (tmp_path / "h.tsv").write_text("id\tgene\na\t7\nb\tx\n")
    text = block("input h", "adapter: tsv\npath: h.tsv\n") + block(
        "calls k", "source: h\nkey: id\ncolumn: gene\ncalls:\n  seven: {any: [7]}\n")
    con = cli.build_report(write(tmp_path, text))[2]
    assert con.execute("select * from k").fetchall() == [("a", True), ("b", False)]


@pytest.mark.parametrize("body, error", [
    ("calls:\n  X: {any: [a]}\n", "needs 'source'"),
    ("source: p\n", "needs 'calls'"),
    ("source: p\ncalls:\n  X: {any: [a]}\nfilters: {}\n", "unknown key 'filters'; accepted: source, samples, column, key, calls"),
    ("source: p\ncalls:\n  X: [a]\n", "call 'X' must be {any: [values]} or {all: [values]}"),
    ("source: p\ncalls:\n  X: {any: [a], all: [b]}\n", "call 'X' must be {any: [values]} or {all: [values]}"),
    ("source: p\ncalls:\n  X: {some: [a]}\n", "call 'X' must be {any: [values]} or {all: [values]}"),
    ("source: p\ncalls:\n  X: {any: []}\n", "call 'X': 'any' must be a non-empty list of values"),
    ("source: p\ncalls:\n  X: {all: a}\n", "call 'X': 'all' must be a non-empty list of values"),
    ("source: p\ncalls:\n  Sample: {any: [a]}\n", "call 'Sample': the name is taken by the key or another call"),
    ("source: p\ncalls:\n  X: {any: [a]}\n  x: {any: [b]}\n", "call 'x': the name is taken by the key or another call"),
    ("source: p\ncolumn: [g]\ncalls:\n  X: {any: [a]}\n", "'column' must be one column"),
    ("source: p\ncolumn: g\ncalls:\n  X: {any: [a]}\n", "'g' not a column of 'p'; its columns: file, sample, gene"),
    ("source: p\nsamples: s\ncolumn: gene\ncalls:\n  X: {any: [a]}\n", "'sample' not a column of 's'; its columns: file, id"),
    ("source: p\nsamples: nope\ncalls:\n  X: {any: [a]}\n", "samples 'nope' names no block; readable here: p, s"),
])
def test_calls_shape(tmp_path, body, error):
    (tmp_path / "p.tsv").write_text("sample\tgene\na\tx\n")
    (tmp_path / "s.tsv").write_text("id\na\n")
    text = (block("input p", "adapter: tsv\npath: p.tsv\n") + block("input s", "adapter: tsv\npath: s.tsv\n")
            + block("calls c", body))
    with pytest.raises(cli.SpecError, match=re.escape(f"calls 'c': {error}")):
        build(write(tmp_path, text))


def test_list_and_distinct_measures(tmp_path):
    (tmp_path / "h.tsv").write_text("g\tv\tn\na\tz\t10\na\ty\t9\na\tz\t9\na\t\t\nb\tx\t1\n")
    measures = "  genes: list(v)\n  unique: distinct(v)\n  nums: list(n)\n  unums: distinct(n)\n"
    text = (block("input h", "adapter: tsv\npath: h.tsv\n")
            + block("query m", f"source: h\ngroup_by: g\nmeasures:\n{measures}sort: [g]\n"))
    con = cli.build_report(write(tmp_path, text))[2]
    # sorted, nulls skipped; numbers sort as numbers
    assert con.execute("select * from m").fetchall() == [("a", "y, z, z", "y, z", "9, 9, 10", "9, 10"), ("b", "x", "x", "1", "1")]
    with pytest.raises(cli.SpecError, match=re.escape("expected fn(column) with fn in count, min, max, sum, avg, list, distinct")):
        build(write(tmp_path, block("input h", "adapter: tsv\npath: h.tsv\n")
                    + block("query m", "source: h\ngroup_by: g\nmeasures:\n  x: median(n)\n")))
