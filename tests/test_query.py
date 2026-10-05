import pytest

from groundie import query as q


@pytest.mark.parametrize("block, error", [
    ({"source": "a", "columns": ["x"], "group_by": ["x"]}, "datatable 't': give 'columns' or 'group_by', not both"),
    ({"source": "a", "group_by": ["x"], "measures": {"n": "count(x"}}, "measure 'n': expected fn(column)"),
    ({"source": "a", "group_by": ["x"], "measures": {"n": "median(y)"}}, "measure 'n': expected fn(column)"),
    ({"source": "a", "group_by": ["x"], "measures": {"n": "max(*)"}}, "measure 'n': expected fn(column)"),
    ({"source": "a", "group_by": ["x"], "measures": {"n": "count( y )"}}, "measure 'n': expected fn(column)"),
    ({"source": "a", "measures": {"n": "count(*)"}}, "'measures' needs 'group_by'"),
    ({"source": "a", "group_by": ["x"], "measures": {"x": "count(*)"}}, "measure 'x' has the name of a group_by"),
    ({"source": "a", "group_by": ["x"], "measures": {"n": "count(*)"}, "sort": ["y"]},
     "sort names 'y', not an output column; output columns: x, n"),
    ({"source": "a", "filter": {"x": "1"}, "limit": 3},
     "datatable 't': unknown key 'filter', 'limit'; accepted: source, columns, filters, group_by, measures, sort"),
    ({"source": "a", "sort": 3}, "datatable 't': 'sort' must be one column or a list of columns"),
    ({"source": "a", "columns": {"x": 1}}, "datatable 't': 'columns' must be one column or a list of columns"),
    ({"source": "a", "group_by": ["x", 1]}, "datatable 't': 'group_by' must be one column or a list of columns"),
])
def test_normalize_errors(block, error):
    with pytest.raises(ValueError, match=error.replace("(", r"\(").replace(")", r"\)").replace("*", r"\*")):
        q.normalize("t", block)


def test_grouped_to_sql():
    query = q.normalize("t", {"source": "amr", "filters": {"element_type": "AMR"},
                              "group_by": ["sample", "gene", "sample"],
                              "measures": {"best": "max(pct_identity)", "n": "count(*)"}, "sort": ["n", "sample"]})
    assert q.to_sql(query) == (
        'select "sample", "gene", max("pct_identity") as "best", count(*) as "n" from "amr" '
        'where "element_type" in (?) group by "sample", "gene" order by "n", "sample"', ["AMR"])
    amr = q.Input("amr", "amrfinder", ["/data/a.tsv"])
    assert q.about(query) == ("t: filters element_type = AMR; one row per sample, gene; "
                              "best = max(pct_identity), n = count(*); sorted by n, sample")
    assert q.about(amr) == "amr: adapter amrfinder, file a.tsv"
    assert q.to_yaml(query) == ("source: amr\nfilters:\n  element_type: [AMR]\ngroup_by: [sample, gene]\n"
                                "measures:\n  best: max(pct_identity)\n  n: count(*)\nsort: [n, sample]\n")


def test_scalar_column_keys():
    query = q.normalize("t", {"source": "a", "columns": "x", "sort": "x"})
    assert (query.columns, query.sort) == (["x"], ["x"])
    assert q.to_sql(query) == ('select "x" from "a" order by "x"', [])
    assert q.normalize("t", {"source": "a", "group_by": "x"}).group_by == ["x"]


def test_needs_source():
    with pytest.raises(ValueError, match="query 't': needs 'source', the name of the block it reads"):
        q.normalize("t", {"columns": ["x"]}, "query")
    with pytest.raises(ValueError, match=r"query 't': unknown key 'result'"):
        q.normalize("t", {"result": "a"}, "query")


def test_lineage_and_about():
    blocks = {
        "a": q.Input("a", "tsv", ["x/one.tsv", "y/two.tsv"]),
        "s": q.Sql("s", "select * from a"),
        "b": q.normalize("b", {"source": "a", "filters": {"k": ["1", "2"]}}, "query"),
        "c": q.normalize("c", {"source": "b", "sort": "k"}, "datatable"),
        "d": q.normalize("d", {"source": "s"}, "datatable"),
    }
    assert [b.name for b in q.lineage(blocks, "c")] == ["c", "b", "a"]
    assert [b.name for b in q.lineage(blocks, "d")] == ["d", "s"]  # ad hoc SQL ends the walk
    assert [q.about(b) for b in q.lineage(blocks, "c")] == [
        "c: no filters; sorted by k", "b: filters k in 1, 2", "a: adapter tsv, 2 files: one.tsv, two.tsv"]
    assert q.about(blocks["s"]) == "s: ad hoc SQL"
    assert q.to_yaml(blocks["a"]) == "adapter: tsv\npaths: [x/one.tsv, y/two.tsv]\n"
