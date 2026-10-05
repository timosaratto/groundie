import pytest

from groundie.adapters import newick

# Rows are (tree, node, parent, label, branch_length, support, is_tip, comment).
CASES = {
    "lengths and nesting": ("((A:1,B:2)C:0.5,D:3);", [
        (1, 1, None, None, None, None, False, None),
        (1, 2, 1, "C", 0.5, None, False, None),
        (1, 3, 2, "A", 1.0, None, True, None),
        (1, 4, 2, "B", 2.0, None, True, None),
        (1, 5, 1, "D", 3.0, None, True, None)]),
    "quoted labels, doubled apostrophes, underscores kept": ("('it''s here':1,B_c,'x_(y)');", [
        (1, 1, None, None, None, None, False, None),
        (1, 2, 1, "it's here", 1.0, None, True, None),
        (1, 3, 1, "B_c", None, None, True, None),
        (1, 4, 1, "x_(y)", None, None, True, None)]),
    "NHX and BEAST comments, anywhere": ("[&R] (A[&rate=1]:1[&&NHX:S=x],B)[c];", [
        (1, 1, None, None, None, None, False, "[&R][c]"),
        (1, 2, 1, "A", 1.0, None, True, "[&rate=1][&&NHX:S=x]"),
        (1, 3, 1, "B", None, None, True, None)]),
    "scientific notation and negative lengths": ("(A:1e-5,B:-2.5E+1,C:.5);", [
        (1, 1, None, None, None, None, False, None),
        (1, 2, 1, "A", 1e-5, None, True, None),
        (1, 3, 1, "B", -25.0, None, True, None),
        (1, 4, 1, "C", 0.5, None, True, None)]),
    "polytomy, single child, empty labels": ("((A),,(B,C,D));", [
        (1, 1, None, None, None, None, False, None),
        (1, 2, 1, None, None, None, False, None),
        (1, 3, 2, "A", None, None, True, None),
        (1, 4, 1, None, None, None, True, None),
        (1, 5, 1, None, None, None, False, None),
        (1, 6, 5, "B", None, None, True, None),
        (1, 7, 5, "C", None, None, True, None),
        (1, 8, 5, "D", None, None, True, None)]),
    "support on internal nodes, text kept": ("((A,B)95:1,(C,D)0.87,E)x;", [
        (1, 1, None, "x", None, None, False, None),
        (1, 2, 1, "95", 1.0, 95.0, False, None),
        (1, 3, 2, "A", None, None, True, None),
        (1, 4, 2, "B", None, None, True, None),
        (1, 5, 1, "0.87", None, 0.87, False, None),
        (1, 6, 5, "C", None, None, True, None),
        (1, 7, 5, "D", None, None, True, None),
        (1, 8, 1, "E", None, None, True, None)]),
    "numeric tip label is not support": ("(1,2);", [
        (1, 1, None, None, None, None, False, None),
        (1, 2, 1, "1", None, None, True, None),
        (1, 3, 1, "2", None, None, True, None)]),
    "duplicate tip labels": ("(A,A);", [
        (1, 1, None, None, None, None, False, None),
        (1, 2, 1, "A", None, None, True, None),
        (1, 3, 1, "A", None, None, True, None)]),
    "several trees, whitespace and newlines anywhere, trailing newline": ("( A :\n 1 ,\tB ) ;\n\n(C);\n", [
        (1, 1, None, None, None, None, False, None),
        (1, 2, 1, "A", 1.0, None, True, None),
        (1, 3, 1, "B", None, None, True, None),
        (2, 1, None, None, None, None, False, None),
        (2, 2, 1, "C", None, None, True, None)]),
    "a single tip": ("A;", [(1, 1, None, "A", None, None, True, None)]),
}


@pytest.mark.parametrize("text, expected", CASES.values(), ids=CASES.keys())
def test_parse(text, expected):
    assert newick.parse(text, "t.nwk") == expected


@pytest.mark.parametrize("text, error", [
    ("(A,B", "t.nwk, character 4: unbalanced '(': 1 not closed"),
    ("(A,B));", "t.nwk, character 6: unbalanced ')'"),
    ("('A,B);", "t.nwk, character 2: quoted label is never closed"),
    ("(A,B)[x;", "t.nwk, character 6: comment '[' is never closed"),
    ("(A,B)", "t.nwk, character 5: missing ';' at the end of the last tree"),
    ("(A,B);(C)", "t.nwk, character 9: missing ';'"),
    ("", "t.nwk: no tree found"),
    ("  \n", "t.nwk: no tree found"),
    (";", "t.nwk, character 1: empty tree before ';'"),
    ("(A:,B);", "t.nwk, character 4: expected a branch length after ':'"),
    ("(A:x,B);", "t.nwk, character 4: expected a branch length after ':'"),
    ("(A B);", "t.nwk, character 4: unexpected label 'B'"),
    ("(A,B)(C);", "t.nwk, character 6: unexpected '('"),
    ("A,B;", "t.nwk, character 2: ',' outside parentheses"),
])
def test_parse_errors(text, error):
    with pytest.raises(ValueError, match=error.replace("(", r"\(").replace(")", r"\)").replace("[", r"\[")):
        newick.parse(text, "t.nwk")


def test_deep_tree():
    depth = 5000
    rows = newick.parse("(" * depth + "A" + ",B)" * depth + ";", "t.nwk")
    assert len(rows) == 2 * depth + 1
    assert rows[depth] == (1, depth + 1, depth, "A", None, None, True, None)
    assert rows[-1] == (1, 2 * depth + 1, 1, "B", None, None, True, None)


def test_rows_adds_file(tmp_path):
    path = tmp_path / "t.nwk"
    path.write_text("(A,B);\n")
    assert newick.rows(path, {})[1] == ("t.nwk", 1, 2, 1, "A", None, None, True, None)
