import argparse
import glob
import hashlib
import html
import json
import platform
import re
import sys
from contextlib import contextmanager
from dataclasses import dataclass
from importlib.resources import files
from datetime import datetime, timezone
from importlib.metadata import version as pkg_version
from pathlib import Path

import duckdb
import yaml
from markdown_it import MarkdownIt

from groundie.adapters import ADAPTERS
from groundie import query as q
from groundie import tables, tree
from groundie.query import SpecError, ident

BLOCK = re.compile(r"^```(\w+) (\w+)\n(.*?)^```[ \t]*$", re.MULTILINE | re.DOTALL)
# Data blocks are lowercase; display blocks, which show something, are PascalCase, like JSX components.
KINDS = ("input", "sql", "query", "clusters", "checks", "calls", "DataTable", "Matrix", "AmrfinderMatrix",
         "MlstTable", "Tree")
MATRICES = ("Matrix", "AmrfinderMatrix")
TABLES = ("DataTable", "MlstTable")
SUGAR = {"AmrfinderMatrix": "amrfinder", "MlstTable": "mlst"}  # a block kind for one adapter's input
INPUT_KEYS = ("adapter", "path", "paths", "options")
LAYOUT_KEYS = ("rows", "cols", "value", "scale", "digits", "mark", "pick", "symmetric")
PICK = re.compile(r"(max|min)\((\w+)\)")
CLUSTERS_KEYS = ("source", "threshold", "pairs", "distance", "linkage")
CHECKS_KEYS = ("source", "key", "checks")
CALLS_KEYS = ("source", "samples", "column", "key", "calls")
TREE_KEYS = ("source", "tips", "match", "label", "color", "file", "tree")


@dataclass
class Layout:
    """How a matrix block pivots its query's rows: not part of the query."""
    rows: list[str]
    cols: str
    value: str | None = None  # the column a cell shows: as a shade with scale, as text without
    scale: list | None = None  # [low, high]: the value range the shade spans; values outside are clamped
    digits: int | None = None  # decimal places of a value shown as text
    mark: dict | None = None  # {column: <output column>, below: <number>}
    pick: str | None = None  # 'max(column)' or 'min(column)': which of several rows in one cell it shows
    symmetric: bool = False  # a row for (a, b) also fills (b, a); rows and columns are the same keys


def layout(name, block):
    """Take a matrix block's layout keys out of the block (the rest is the query) and check their shape."""
    keys = {k: block.pop(k) for k in LAYOUT_KEYS if k in block}
    rows = keys.get("rows")
    if isinstance(rows, str):
        rows = [rows]
    if not rows or not isinstance(keys.get("cols"), str):
        raise SpecError(f"Matrix {name!r}: needs 'rows' (one or more columns) and 'cols' (one column)")
    mark = keys.get("mark")
    if mark is not None and (not isinstance(mark, dict) or set(mark) != {"column", "below"}
                             or not is_num(mark["below"])):
        raise SpecError(f"Matrix {name!r}: 'mark' must be {{column: <column>, below: <number>}}")
    value, scale, digits = keys.get("value"), keys.get("scale"), keys.get("digits")
    if scale is not None and value is None:
        raise SpecError(f"Matrix {name!r}: 'scale' needs 'value', the column it shades")
    if digits is not None and (value is None or scale is not None or not isinstance(digits, int)
                               or isinstance(digits, bool) or digits < 0):
        raise SpecError(f"Matrix {name!r}: 'digits' is a whole number of decimal places, for a 'value' "
                        "shown as text (no 'scale')")
    if scale is not None and not (isinstance(scale, list) and len(scale) == 2 and all(map(is_num, scale))
                                  and scale[0] < scale[1]):
        raise SpecError(f"Matrix {name!r}: 'scale' must be [low, high] with low below high")
    if "pick" in keys and not PICK.fullmatch(str(keys["pick"])):
        raise SpecError(f"Matrix {name!r}: 'pick' must be max(<column>) or min(<column>)")
    symmetric = keys.get("symmetric", False)
    if not isinstance(symmetric, bool) or symmetric and len(rows) != 1:
        raise SpecError(f"Matrix {name!r}: 'symmetric' is true or false, and true needs one 'rows' column")
    return Layout(rows, keys["cols"], value, scale, digits, mark, keys.get("pick"), symmetric)


def display(name, block):
    """Take a datatable's display key out of the block (the rest is the query): each column's title, by
    column, or the name of an input whose adapter's titles it shows. Whether the titles name output columns,
    and the input is one, are checked later, when they are known."""
    value = block.pop("display", None) or {}
    if isinstance(value, str):
        return value
    if not isinstance(value, dict):
        raise SpecError(f"DataTable {name!r}: 'display' maps columns to their titles, e.g. {{sequence_type: ST}}, "
                        "or names an input whose adapter's titles it shows")
    for column, title in value.items():
        if not isinstance(title, str) or not title.strip():
            raise SpecError(f"DataTable {name!r}: display {str(column)!r} must be text, the column's title")
    return {str(column): title for column, title in value.items()}


AMRFINDER_MATRIX_KEYS = ("source", "filters", "sort", "partial_below")


def amrfinder_matrix(name, block):
    """An AmrfinderMatrix block as the query and layout of the matrix it stands for: genes by sample, the
    element type beside the gene unless a filter picks one, and in a cell with several hits the most
    complete, marked when even that one is partial: coverage below partial_below, default 100."""
    who = f"AmrfinderMatrix {name!r}"
    unknown = [k for k in block if k not in AMRFINDER_MATRIX_KEYS]
    if unknown:
        raise SpecError(f"{who}: unknown key {', '.join(map(repr, unknown))}; accepted: "
                         f"{', '.join(AMRFINDER_MATRIX_KEYS)}. For another layout, use a Matrix block")
    block = dict(block)
    partial_below = block.pop("partial_below", 100)
    if not is_num(partial_below):
        raise SpecError(f"{who}: 'partial_below' must be a number, the coverage (%) below which a hit is partial")
    query = q.normalize(name, block, "AmrfinderMatrix")
    rows = ["element_symbol"] if len(query.filters.get("element_type", [])) == 1 else ["element_type", "element_symbol"]
    query.columns = [*rows, "sample", "pct_identity", "pct_coverage"]
    return query, Layout(rows, "sample", mark={"column": "pct_coverage", "below": partial_below}, pick="max(pct_coverage)")


MLST_TABLE_KEYS = ("source", "filters", "sort")


def mlst_table(name, block):
    """An MlstTable block as its query: every column of the mlst input. The tables are split at build."""
    unknown = [k for k in block if k not in MLST_TABLE_KEYS]
    if unknown:
        raise SpecError(f"MlstTable {name!r}: unknown key {', '.join(map(repr, unknown))}; accepted: "
                         f"{', '.join(MLST_TABLE_KEYS)}. For other columns, use a DataTable")
    return q.normalize(name, block, "MlstTable")


@dataclass
class TreeBlock:
    """What a tree block draws: one tree of a newick input, its tips annotated from another block. Not a
    query, and it makes no SQL table: no block can read it."""
    name: str
    source: str
    tips: str | None = None  # a block with one row per tip
    match: str | None = None  # the column of tips holding the raw tip label
    label: list[str] | None = None  # columns of tips shown right of the tips; none means the raw tip label
    color: str | None = None  # a column of tips, categorical
    file: str | None = None  # which tree, when the source holds several: the file's basename
    tree: int | None = None  # and the tree's position in that file, from 1


def tree_block(name, block):
    """Check a tree block's shape. Its source and tips are checked by check_reads, its columns at build."""
    unknown = [k for k in block if k not in TREE_KEYS]
    if unknown:
        raise SpecError(f"Tree {name!r}: unknown key {', '.join(map(repr, unknown))}; accepted: {', '.join(TREE_KEYS)}")
    if not isinstance(block.get("source"), str):
        raise SpecError(f"Tree {name!r}: needs 'source', an input loaded by the newick adapter")
    label = [block["label"]] if isinstance(block.get("label"), str) else block.get("label")
    given = [k for k in ("match", "label", "color") if k in block]
    if "tips" not in block and given:
        raise SpecError(f"Tree {name!r}: {', '.join(map(repr, given))} need 'tips', the block to read them from")
    if "tips" in block and not (isinstance(block["tips"], str) and isinstance(block.get("match"), str)):
        raise SpecError(f"Tree {name!r}: 'tips' names a block and needs 'match', the column holding the tip label")
    if label is not None and not (isinstance(label, list) and label and all(isinstance(c, str) for c in label)):
        raise SpecError(f"Tree {name!r}: 'label' must be one column or a list of columns")
    if "color" in block and not isinstance(block["color"], str):
        raise SpecError(f"Tree {name!r}: 'color' must be one column")
    if "file" in block and not isinstance(block["file"], str):
        raise SpecError(f"Tree {name!r}: 'file' must be a file name, as in the input's 'file' column")
    if "tree" in block and not (type(block["tree"]) is int and block["tree"] >= 1):
        raise SpecError(f"Tree {name!r}: 'tree' must be a whole number from 1, the tree's position in its file")
    return TreeBlock(name, block["source"], block.get("tips"), block.get("match"), label, block.get("color"),
                     block.get("file"), block.get("tree"))


def clusters_block(name, block):
    """Check a clusters block's shape. Its columns are checked at build."""
    who = f"clusters {name!r}"
    unknown = [k for k in block if k not in CLUSTERS_KEYS]
    if unknown:
        raise SpecError(f"{who}: unknown key {', '.join(map(repr, unknown))}; accepted: {', '.join(CLUSTERS_KEYS)}")
    if not isinstance(block.get("source"), str):
        raise SpecError(f"{who}: needs 'source', the block holding one row per pair of samples")
    if not is_num(block.get("threshold")) or block["threshold"] < 0:
        raise SpecError(f"{who}: needs 'threshold', a distance of 0 or more; pairs this close or closer are linked")
    pairs = block.get("pairs", ["reference", "query"])
    if not (isinstance(pairs, list) and len(pairs) == 2 and all(isinstance(c, str) for c in pairs)):
        raise SpecError(f"{who}: 'pairs' must be the two columns naming the samples of a pair")
    if not isinstance(block.get("distance", ""), str):
        raise SpecError(f"{who}: 'distance' must be one column")
    if block.get("linkage", "single") != "single":
        raise SpecError(f"{who}: 'linkage' must be single, the only linkage there is so far")
    return q.Clusters(name, block["source"], block["threshold"], pairs, block.get("distance", "distance"))


def checks_block(name, block):
    """Check a checks block's shape. Its columns are checked at build."""
    who = f"checks {name!r}"
    unknown = [k for k in block if k not in CHECKS_KEYS]
    if unknown:
        raise SpecError(f"{who}: unknown key {', '.join(map(repr, unknown))}; accepted: {', '.join(CHECKS_KEYS)}")
    if not isinstance(block.get("source"), str):
        raise SpecError(f"{who}: needs 'source', the block holding the rows to check")
    key = block.get("key", "sample")
    if not isinstance(key, str):
        raise SpecError(f"{who}: 'key' must be one column, the one identifying a row")
    checks = block.get("checks")
    if not (isinstance(checks, dict) and checks):
        raise SpecError(f"{who}: needs 'checks': column: {{min, max, above or below: <number>}}, one or more")
    for column, bounds in checks.items():
        what = f"{who}: check {column!r}"
        if column == key:
            raise SpecError(f"{what} is on the key column; the key is shown as it is")
        if not (isinstance(bounds, dict) and bounds):
            raise SpecError(f"{what} must be {{min, max, above or below: <number>}}, one or more")
        unknown = [b for b in bounds if b not in q.BOUNDS]
        if unknown:
            raise SpecError(f"{what}: unknown bound {', '.join(map(repr, unknown))}; accepted: {', '.join(q.BOUNDS)}")
        if not all(map(is_num, bounds.values())):
            raise SpecError(f"{what}: bounds must be numbers")
        if {"min", "above"} <= set(bounds) or {"max", "below"} <= set(bounds):
            raise SpecError(f"{what}: give min or above, not both, and max or below, not both")
    return q.Checks(name, block["source"], {str(c): b for c, b in checks.items()}, key)


def calls_block(name, block):
    """Check a calls block's shape. Its source and samples are checked by check_reads, its columns at build."""
    who = f"calls {name!r}"
    unknown = [k for k in block if k not in CALLS_KEYS]
    if unknown:
        raise SpecError(f"{who}: unknown key {', '.join(map(repr, unknown))}; accepted: {', '.join(CALLS_KEYS)}")
    if not isinstance(block.get("source"), str):
        raise SpecError(f"{who}: needs 'source', the block holding one row per sample and value, such as amrfinder")
    for k in ("samples", "column", "key"):
        if k in block and not isinstance(block[k], str):
            raise SpecError(f"{who}: {k!r} must be one {'block' if k == 'samples' else 'column'}")
    key = block.get("key", "sample")
    calls = block.get("calls")
    if not (isinstance(calls, dict) and calls):
        raise SpecError(f"{who}: needs 'calls': name: {{any or all: [values]}}, one or more")
    names = [key.lower()]
    for call, rule in calls.items():
        what = f"{who}: call {str(call)!r}"
        if str(call).lower() in names:  # SQL column names ignore case
            raise SpecError(f"{what}: the name is taken by the key or another call")
        names.append(str(call).lower())
        if not (isinstance(rule, dict) and len(rule) == 1 and next(iter(rule)) in ("any", "all")):
            raise SpecError(f"{what} must be {{any: [values]}} or {{all: [values]}}")
        values = next(iter(rule.values()))
        if not (isinstance(values, list) and values) or any(v is None or isinstance(v, (list, dict)) for v in values):
            raise SpecError(f"{what}: {next(iter(rule))!r} must be a non-empty list of values")
    # YAML types bare values (007 -> 7); values are always strings, as in filters.
    calls = {str(c): {how: [str(v) for v in values]} for c, rule in calls.items() for how, values in rule.items()}
    return q.Calls(name, block["source"], calls, block.get("samples"), block.get("column", "element_symbol"), key)


def parse_spec(text, name="spec"):
    """The segments: prose and blocks, each with the line it starts on. 'name' is the file's name, for
    messages. Checks that block names are unique and that each block reads a block it can."""
    segments, lines = [], {}
    pos = 0
    for m in BLOCK.finditer(text):
        segments.append({"type": "prose", "text": text[pos : m.start()], "line": 1 + text.count("\n", 0, pos)})
        line = 1 + text.count("\n", 0, m.start())
        kind, block = m.group(1), m.group(2)
        with at(f"{name}:{line}"):
            if block in lines:
                raise SpecError(f"{kind} {block!r}: the name {block!r} is taken by the block on line {lines[block]}")
            segments.append({**parse_block(kind, block, m.group(3), line), "line": line})
        lines[block] = line
        pos = m.end()
    segments.append({"type": "prose", "text": text[pos:], "line": 1 + text.count("\n", 0, pos)})
    check_reads(name, [seg for seg in segments if seg["type"] != "prose"])
    return segments


def parse_block(kind, name, text, line):
    """One block's segment, without its line. 'line' is the line of its opening fence."""
    if kind not in KINDS:
        meant = [k for k in KINDS if k.lower() == kind.lower().replace("_", "")]
        raise SpecError(f"unknown block {kind!r} ({kind} {name}); "
                         + (f"did you mean {meant[0]}? " if meant else "") + f"known: {', '.join(KINDS)}")
    if kind == "sql":
        return {"type": kind, "block": sql_block(name, text)}
    try:
        block = yaml.safe_load(text) or {}
    except yaml.YAMLError as e:
        raise SpecError(f"{kind} {name!r}: invalid YAML on line {yaml_line(e, line + 1)}: {yaml_problem(e)}") from e
    if not isinstance(block, dict):
        raise SpecError(f"{kind} {name!r}: the body must be YAML keys and values")
    if kind == "input":
        return {"type": kind, "block": input_block(name, block)}
    if kind == "Tree":  # not a query: it reads its input as it is
        return {"type": kind, "block": tree_block(name, block)}
    if kind == "clusters":
        return {"type": kind, "block": clusters_block(name, block)}
    if kind == "checks":
        return {"type": kind, "block": checks_block(name, block)}
    if kind == "calls":
        return {"type": kind, "block": calls_block(name, block)}
    if kind == "AmrfinderMatrix":
        query, grid = amrfinder_matrix(name, block)
        return {"type": kind, "block": query, "layout": grid}
    if kind == "MlstTable":
        return {"type": kind, "block": mlst_table(name, block)}
    if kind == "DataTable":
        titles = display(name, block)
        return {"type": kind, "block": q.normalize(name, block, kind, ("display",)), "display": titles}
    grid =layout(name, block) if kind == "Matrix" else None
    return {"type": kind, "block": q.normalize(name, block, kind, LAYOUT_KEYS if kind == "Matrix" else ()),
            "layout": grid}


def check_reads(name, segments):
    """Each block reads a block it can: any input, or a block above it that is not a tree. A tree's
    source is an input loaded by the newick adapter."""
    blocks = {seg["block"].name: seg for seg in segments}
    above = set()
    for seg in segments:
        block = seg["block"]
        who = f"{seg['type']} {block.name!r}"
        with at(f"{name}:{seg['line']}"):
            for key in ("source", "tips", "samples"):
                target = getattr(block, key, None)
                if target is None:
                    continue
                if target not in blocks:
                    readable = [n for n, s in blocks.items() if s["type"] == "input" or n in above]
                    raise SpecError(f"{who}: {key} {target!r} names no block; readable here: "
                                     f"{', '.join(readable) or 'none'}")
                if blocks[target]["type"] == "Tree":
                    raise SpecError(f"{who}: {key} {target!r} is a Tree block, which no block can read")
                if blocks[target]["type"] != "input" and target not in above:
                    raise SpecError(f"{who}: {key} {target!r} is not above this block; "
                                     "a block reads inputs and the blocks above it")
            if isinstance(seg.get("display"), str):
                check_display(who, seg["display"], blocks)
            if seg["type"] == "Tree":
                source = blocks[block.source]
                if source["type"] != "input" or source["block"].adapter != "newick":
                    what = (f"is loaded by adapter {source['block'].adapter!r}" if source["type"] == "input"
                            else f"is a {source['type']} block")
                    raise SpecError(f"{who}: source {block.source!r} {what}; "
                                     "a Tree block needs an input loaded by the newick adapter")
            if seg["type"] in SUGAR:
                source, adapter = blocks[block.source], SUGAR[seg["type"]]
                if source["type"] != "input" or source["block"].adapter != adapter:
                    what = (f"is loaded by adapter {source['block'].adapter!r}" if source["type"] == "input"
                            else f"is a {source['type']} block")
                    raise SpecError(f"{who}: source {block.source!r} {what}; an {seg['type']} needs an "
                                     f"input loaded by the {adapter} adapter; narrow it with filters")
        if seg["type"] != "Tree":
            above.add(block.name)


def check_display(who, target, blocks):
    """A datatable's 'display: <name>' names an input loaded by an adapter that declares titles."""
    titled = ", ".join(a for a, adapter in ADAPTERS.items() if getattr(adapter, "TITLES", None))
    if target not in blocks:
        inputs = [n for n, s in blocks.items() if s["type"] == "input"]
        raise SpecError(f"{who}: display {target!r} names no block; inputs: {', '.join(inputs) or 'none'}")
    if blocks[target]["type"] != "input":
        raise SpecError(f"{who}: display {target!r} is a {blocks[target]['type']} block; display names an "
                         "input, whose adapter's titles it shows")
    adapter = blocks[target]["block"].adapter
    if not getattr(ADAPTERS[adapter], "TITLES", None):
        raise SpecError(f"{who}: display {target!r} is loaded by adapter {adapter!r}, which declares no titles; "
                         f"adapters with titles: {titled}")


def yaml_line(e, first):
    """The file line of a YAML error, given the file line of the snippet's first line."""
    mark = getattr(e, "problem_mark", None) or getattr(e, "context_mark", None)
    return first + mark.line if mark else first


def yaml_problem(e):
    return getattr(e, "problem", None) or str(e).splitlines()[0]


@contextmanager
def at(where):
    """Prefix a SpecError raised inside with where in the file it arose, e.g. 'cohort.md:42'."""
    try:
        yield
    except SpecError as e:
        raise SpecError(f"{where}: {e}") from e


FENCE = re.compile(r"^```(\w+)(.*)$", re.MULTILINE)


def unparsed_fences(seg):
    """(line, message) for each fence in a prose segment that names a block kind but was not read as a block."""
    for m in FENCE.finditer(seg["text"]):
        if m.group(1) not in KINDS:
            continue
        line = seg["line"] + seg["text"].count("\n", 0, m.start())
        if re.fullmatch(r" \w+", m.group(2)):
            problem = "is never closed; close it with a ``` line"
        else:
            problem = (f"is not read as a block: expected ```{m.group(1)} <name>, one space, "
                       "the name letters, digits and underscores")
        yield line, f"fence {m.group(0)!r} {problem}; rendered as a code block"


def sql_block(name, text):
    """A sql block's body, checked to be exactly one SELECT statement."""
    try:
        statements = duckdb.extract_statements(text)
    except duckdb.Error as e:
        raise SpecError(f"sql {name!r}: {e}") from e
    if len(statements) != 1 or statements[0].type != duckdb.StatementType.SELECT:
        found = ", ".join(s.type.name for s in statements) or "none"
        raise SpecError(f"sql {name!r}: the body must be exactly one SELECT statement; found: {found}")
    return q.Sql(name, text.strip())


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def input_block(name, block):
    """Check an input block's keys. Its files are read at build."""
    who = f"input {name!r}"
    unknown = [k for k in block if k not in INPUT_KEYS]
    if unknown:
        raise SpecError(f"{who}: unknown key {', '.join(map(repr, unknown))}; accepted: {', '.join(INPUT_KEYS)}")
    adapter = block.get("adapter")
    if adapter is None:
        raise SpecError(f"{who}: needs 'adapter'; known: {', '.join(ADAPTERS)}")
    if adapter not in ADAPTERS:
        raise SpecError(f"{who}: unknown adapter {adapter!r}; known: {', '.join(ADAPTERS)}")
    if ("path" in block) == ("paths" in block):
        raise SpecError(f"{who}: give exactly one of 'path' or 'paths'")
    paths = [block["path"]] if "path" in block else block["paths"]
    if not isinstance(paths, list) or not paths:
        raise SpecError(f"{who}: 'paths' must be a non-empty list")
    # YAML reads a bare 'null' key as None; the option is spelled null.
    options = {"null" if k is None else k: v for k, v in (block.get("options") or {}).items()}
    unknown = [k for k in options if k not in ADAPTERS[adapter].OPTIONS]
    if unknown:
        raise SpecError(f"{who}: unknown option {', '.join(map(repr, unknown))} for adapter "
                         f"{adapter!r}; accepted: {', '.join(ADAPTERS[adapter].OPTIONS) or 'none'}")
    return q.Input(name, adapter, paths, options)


def resolve(block, written, base_dir):
    """An input path as (path as written, file): relative to base_dir, which is None for inline blocks. A
    glob gives one pair per matching file, in name order, each written as the file's path."""
    if base_dir is None and not Path(written).is_absolute():
        raise SpecError(f"input {block.name!r}: relative path {written!r} in inline blocks; "
                         "use absolute paths or pass base_dir")
    if not glob.has_magic(written):
        path = Path(base_dir or "") / written
        if not path.is_file():
            raise SpecError(f"input {block.name!r}: {written!r} is not a file")
        return [(written, path.resolve())]
    matches = sorted(m for m in glob.glob(written, root_dir=base_dir, recursive=True)
                     if (Path(base_dir or "") / m).is_file())
    if not matches:
        raise SpecError(f"input {block.name!r}: {written!r} matches no file")
    return [(m, (Path(base_dir or "") / m).resolve()) for m in matches]


def load_input(con, block, base_dir):
    """Load an input's files into an SQL table named after the block. Returns its provenance record."""
    adapter = ADAPTERS[block.adapter]
    paths = {}  # a file once, as the first pattern that matches it writes it
    for pattern in block.paths:
        for written, path in resolve(block, pattern, base_dir):
            paths.setdefault(path, (written, path))
    paths = list(paths.values())
    try:
        if hasattr(adapter, "rows"):  # the adapter parses in Python; DuckDB only stores
            counts = insert_rows(con, block.name, adapter, [path for _, path in paths], block.options)
        else:
            literal = "[" + ", ".join(q.literal(path) for _, path in paths) + "]"
            sql = adapter.canonical_sql(block.options)
            counts = [con.execute(f"select count(*) from ({sql.format(path=q.literal(path))})").fetchone()[0]
                      for _, path in paths]
            con.execute(f"create table {ident(block.name)} as {sql.format(path=literal)}")
    except duckdb.Error as e:
        raise SpecError(f"input {block.name!r}: {e}") from e
    files = [{"path": written, "resolved_path": str(path), "sha256": sha256(path),
              "bytes": path.stat().st_size, "rows": rows} for (written, path), rows in zip(paths, counts)]
    return {"adapter": adapter.NAME, "adapter_version": adapter.VERSION,
            **({"options": block.options} if block.options else {}), "files": files}


def insert_rows(con, label, adapter, paths, options):
    """Create the input's SQL table from the adapter's COLUMNS and insert each file's rows. Returns rows per file."""
    con.execute(f"create table {ident(label)} ({', '.join(f'{ident(n)} {t}' for n, t in adapter.COLUMNS)})")
    counts = []
    for path in paths:
        try:
            rows = adapter.rows(path, options)
        except (OSError, ValueError) as e:
            raise SpecError(f"input {label!r}: {e}") from e
        # The rows go in as one JSON string: executemany, or a list parameter per column, takes seconds
        # for a tree of ten thousand tips; this takes a tenth of one.
        names = [n for n, _ in adapter.COLUMNS]
        shape = json.dumps([{n: t for n, t in adapter.COLUMNS}])
        con.execute(f"insert into {ident(label)} select unnest(json_transform(?, ?), recursive := true)",
                    [json.dumps([dict(zip(names, row)) for row in rows]), shape])
        counts.append(len(rows))
    return counts


def make_clusters(con, block):
    """Create the clusters block's SQL table: sample, cluster, cluster_size. Clusters are numbered C1, C2, ...
    by size, largest first, then by their first sample; a sample linked to none is a cluster of its own."""
    who = f"clusters {block.name!r}"
    missing_columns(con, who, block.source, [*block.pairs, block.distance])
    a, b, d = map(ident, [*block.pairs, block.distance])
    rows = con.execute(f"select {a}, {b}, {d} from {ident(block.source)}").fetchall()
    parent = {}
    def root(x):
        while parent.setdefault(x, x) != x:
            parent[x] = x = parent[parent[x]]
        return x
    for x, y, dist in rows:
        if x is None or y is None:
            continue
        rx, ry = root(x), root(y)  # every sample of a pair is a member, linked or not
        if dist is not None and not is_num(dist):
            raise SpecError(f"{who}: distance {block.distance!r} holds {dist!r}, not a number")
        if dist is not None and dist <= block.threshold:
            parent[rx] = ry
    groups = {}
    for s in sorted(parent, key=str):
        groups.setdefault(root(s), []).append(s)
    ordered = sorted(groups.values(), key=lambda g: (-len(g), str(g[0])))
    members = [(s, f"C{i}", len(g)) for i, g in enumerate(ordered, 1) for s in g]
    create(con, who, block.name, "select unnest(?) as sample, unnest(?) as cluster, unnest(?) as cluster_size",
           [list(col) for col in zip(*members)] if members else [[], [], []])


def missing_columns(con, who, source, wanted):
    """Raise unless every wanted column is a column of the SQL table 'source'."""
    columns = size(con, source)[1]
    missing = [c for c in dict.fromkeys(wanted) if c not in columns]
    if missing:
        raise SpecError(f"{who}: {', '.join(map(repr, missing))} not a column of {source!r}; "
                         f"its columns: {', '.join(columns)}")


def make_checks(con, block):
    """Create the checks block's SQL table. Returns its SQL."""
    who = f"checks {block.name!r}"
    missing_columns(con, who, block.source, [block.key, *block.checks])
    sql = q.checks_sql(block)
    create(con, who, block.name, sql)
    return sql


def make_calls(con, block):
    """Create the calls block's SQL table. Returns its SQL and parameters."""
    who = f"calls {block.name!r}"
    missing_columns(con, who, block.source, [block.key, block.column])
    if block.samples:
        missing_columns(con, who, block.samples, [block.key])
    sql, params = q.calls_sql(block)
    create(con, who, block.name, sql, params)
    return sql, params


def create(con, who, name, sql, params=()):
    """Run a block's SELECT into an SQL table named after the block."""
    try:
        con.execute(f"create table {ident(name)} as\n{sql}", params)
    except duckdb.Error as e:
        raise SpecError(f"{who}: {e}") from e


def size(con, name):
    """An SQL table's row count and column names."""
    rows = con.execute(f"select count(*) from {ident(name)}").fetchone()[0]
    return rows, [d[0] for d in con.execute(f"select * from {ident(name)} limit 0").description]


MARKDOWN = MarkdownIt("commonmark").enable("table")


def render_prose(text):
    return MARKDOWN.render(text)


def is_num(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool)


def cell(value, tag="td", num=False):
    """A table cell. A null data cell shows a muted en dash."""
    null = value is None and tag == "td"
    text = "–" if null else "" if value is None else html.escape(str(value))
    classes = " ".join(c for c, on in (("num", num), ("null", null)) if on)
    cls = f' class="{classes}"' if classes else ""
    return f"<{tag}{cls}>{text}</{tag}>"


def header(name, label, num):
    """A column header. A label other than the column's name keeps the name in its title."""
    if label == name:
        return cell(name, "th", num)
    cls = ' class="num"' if num else ""
    return f'<th{cls} title="{html.escape(name)}">{html.escape(label)}</th>'


def render_datatable(columns, rows, caption=None):
    """The datatable. columns are (name, label): the column, and the header shown."""
    numeric = [all(v is None or is_num(v) for v in col) and any(v is not None for v in col)
               for col in zip(*rows)] if rows else [False] * len(columns)
    out = ["<table>\n" + (f"<caption>{html.escape(caption)}</caption>\n" if caption else "") + "<thead>\n"]
    out.append("<tr>" + "".join(header(c, label, n) for (c, label), n in zip(columns, numeric)) + "</tr>\n")
    out.append("</thead>\n<tbody>\n")
    for row in rows:
        out.append("<tr>" + "".join(cell(v, "td", n) for v, n in zip(row, numeric)) + "</tr>\n")
    if not rows:
        out.append(f'<tr><td class="empty" colspan="{len(columns)}">No rows</td></tr>\n')
    out.append("</tbody>\n</table>\n")
    return "".join(out)


def details(spec_yaml, chain, notes=()):
    """The footer of a datatable, matrix or tree: the grounded note, and beside it 'View query', which
    expands to any notes on the block, its body, then each block behind it."""
    out = [f'<details class="query"><summary>{GROUNDED}<span class="toggle">View query</span></summary>\n']
    if notes:
        out.append(f'<p>{"<br>".join(html.escape(n) for n in notes)}</p>\n')
    out.append(f"<pre>{html.escape(spec_yaml)}</pre>\n")
    for block in chain:
        if isinstance(block, q.Sql):
            out.append(f"<p>From {html.escape(block.name)} (ad hoc SQL):</p>\n<pre>{html.escape(block.sql)}</pre>\n")
        else:
            out.append(f"<p>From {html.escape(block.name)}:</p>\n<pre>{html.escape(q.to_yaml(block))}</pre>\n")
    out.append("</details>\n")
    return "".join(out)


DEFAULT_ORDER = ("Rows by how many cells they fill, most first, then by name; columns by name. "
                 "Give 'sort' to order both by the query instead.")


def by_name(key):
    """A sort key for a row or column key: by value, nulls last."""
    return tuple((v is None, v) for v in (key if isinstance(key, tuple) else (key,)))


def pivot(query, grid, columns, rows):
    """Row keys and column keys, and for each filled cell the query row it shows and how many rows fell in
    it. With the query's sort, keys are in first-seen order; without, in DEFAULT_ORDER. Of several rows, a
    cell shows the first, or with pick the one with the highest (max) or lowest (min) value of its column."""
    fn, by = PICK.fullmatch(grid.pick).groups() if grid.pick else (None, None)
    missing = [c for c in [*grid.rows, grid.cols, grid.value, (grid.mark or {}).get("column"), by]
               if c is not None and c not in columns]
    if missing:
        raise SpecError(f"Matrix {query.name!r}: {', '.join(map(repr, missing))} not an output column; "
                         f"output columns: {', '.join(columns)}")
    at = {c: i for i, c in enumerate(columns)}
    def better(row, than):  # a row with no number in the pick column never wins
        if not by or not is_num(row[by]):
            return False
        return not is_num(than[by]) or (row[by] > than[by] if fn == "max" else row[by] < than[by])

    row_keys, col_keys, cells, counts = {}, {}, {}, {}
    for row in rows:
        r, c = tuple(row[at[c]] for c in grid.rows), row[at[grid.cols]]
        keys = [(r, c), ((c,), r[0])] if grid.symmetric and (c,) != r else [(r, c)]
        row = dict(zip(columns, row))
        for key in keys:
            row_keys[key[0]] = col_keys[key[1]] = None
            counts[key] = counts.get(key, 0) + 1
            if key not in cells or better(row, cells[key]):
                cells[key] = row
    if grid.symmetric:  # one order for both: first seen with sort, else by name
        order = list(col_keys) if query.sort else sorted(col_keys, key=by_name)
        return [(k,) for k in order], order, cells, counts
    if query.sort:
        return list(row_keys), list(col_keys), cells, counts
    filled = {}
    for rk, _ in cells:
        filled[rk] = filled.get(rk, 0) + 1
    return (sorted(row_keys, key=lambda rk: (-filled[rk], by_name(rk))), sorted(col_keys, key=by_name),
            cells, counts)


def fmt(value):
    return f"{value:.1f}" if isinstance(value, float) else str(value)


def shade(t):
    """Light to dark blue for t in [0, 1], with a text colour that stays readable."""
    light, dark = (219, 231, 246), (31, 78, 150)
    r, g, b = (round(a + (z - a) * t) for a, z in zip(light, dark))
    return f"background:#{r:02x}{g:02x}{b:02x};color:{'#fff' if t > 0.55 else '#1a1a1a'}"


def render_matrix(query, grid, columns, rows):
    """The grid and its legend."""
    row_keys, col_keys, cells, counts = pivot(query, grid, columns, rows)
    low, high = grid.scale or (0, 1)
    out = ['<table class="matrix">\n<thead>\n<tr>' + "".join(cell(c, "th") for c in grid.rows)]
    out.append("".join(f'<th class="col"><span>{html.escape(str(k))}</span></th>' for k in col_keys))
    out.append("</tr>\n</thead>\n<tbody>\n")
    for rk in row_keys:
        out.append("<tr>" + "".join(cell(v) for v in rk))
        for ck in col_keys:
            c = cells.get((rk, ck))
            if c is None:
                out.append("<td></td>")
                continue
            title = ", ".join(f"{k} {fmt(v)}" for k, v in c.items() if k not in grid.rows and k != grid.cols)
            if counts[(rk, ck)] > 1:
                title = f"{counts[(rk, ck)]} rows; shown: {title}" if title else f"{counts[(rk, ck)]} rows"
            text, cls, style = "●", "cell", ""
            if grid.value and grid.scale:
                v = c[grid.value]
                text = ""  # colour only: the value is in the title
                if is_num(v):
                    style = f' style="{shade(min(1, max(0, (v - low) / (high - low))))}"'
            elif grid.value:
                v = c[grid.value]
                text = "–" if v is None else html.escape(
                    f"{v:.{grid.digits}f}" if is_num(v) and grid.digits is not None else str(v))
                cls = "cell num null" if v is None else "cell num"
            if grid.mark and is_num(c[grid.mark["column"]]) and c[grid.mark["column"]] < grid.mark["below"]:
                text += '<span class="mark">◆</span>'
            out.append(f'<td class="{cls}" title="{html.escape(title)}"{style}>{text}</td>')
        out.append("</tr>\n")
    out.append("</tbody>\n</table>\n")
    if not grid.value:
        legend = ["● a hit"]
    elif grid.scale:
        legend = [f"shade: {grid.value} from {low} (light) to {high} (dark); blank: no hit"]
    else:
        legend = [f"cells: {grid.value}; blank: no row"]
    if grid.mark:
        legend.append(f"◆ {grid.mark['column']} below {grid.mark['below']}")
    out.append(f'<p class="legend">{html.escape("; ".join(legend))}</p>\n')
    return "".join(out)


def render_provenance(manifest):
    out = ['<footer class="provenance">\n<h2>Provenance</h2>\n<table>\n<thead>\n']
    out.append("<tr><th>Input</th><th>Adapter</th><th>Path</th><th>SHA-256</th>"
               '<th class="num">Rows</th><th class="num">Bytes</th></tr>\n</thead>\n<tbody>\n')
    for label, src in manifest["inputs"].items():
        for f in src["files"]:
            out.append(
                f"<tr><td>{html.escape(label)}</td>"
                f"<td>{html.escape(src['adapter'])} {html.escape(src['adapter_version'])}</td>"
                f"<td><code>{html.escape(f['path'])}</code></td>"
                f"<td><code title=\"{f['sha256']}\">{f['sha256'][:12]}</code></td>"
                f"<td class=\"num\">{f['rows']}</td><td class=\"num\">{f['bytes']}</td></tr>\n"
            )
    out.append("</tbody>\n</table>\n")
    env = manifest["environment"]
    out.append(
        f"<p>File <code>{html.escape(manifest['file']['path'])}</code> "
        f"<code title=\"{manifest['file']['sha256']}\">{manifest['file']['sha256'][:12]}</code>. "
        f"Built {html.escape(manifest['built_at'])} with groundie {env['groundie']}, "
        f"DuckDB {env['duckdb']}, Python {env['python']}.</p>\n</footer>\n"
    )
    return "".join(out)


def build(path: Path) -> str:
    """The report page of a file of blocks."""
    body, manifest, _ = build_report(path)
    return page(Path(path).stem, body + render_provenance(manifest))


def page(title, body):
    css = files("groundie").joinpath("style.css").read_text()
    return (
        "<!doctype html>\n<html lang=\"en\">\n<head>\n<meta charset=\"utf-8\">\n"
        "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n"
        f"<title>{html.escape(title)}</title>\n"
        f"<style>\n{css}</style>\n"
        f"</head>\n<body>\n<main>\n{body}</main>\n</body>\n</html>\n"
    )


def build_report(path, text=None, base_dir=None):
    """The blocks' HTML, prose included and without page or provenance; the manifest, which lists the
    blocks in file order and the build's warnings; and the DuckDB connection, holding an SQL table per
    block but tree. With 'text', the blocks are that text, 'path' only names them in messages, and relative
    input paths resolve against base_dir, or are an error without it."""
    path = Path(path)
    if text is None:
        try:
            text = path.read_text()
        except OSError as e:
            raise SpecError(f"cannot read {str(path)!r}: {e.strerror or e}") from e
        base_dir = path.parent
    segments = parse_spec(text, path.name)
    blocks = {seg["block"].name: seg["block"] for seg in segments if seg["type"] != "prose"}
    manifest = {
        "file": {"path": path.name, "sha256": hashlib.sha256(text.encode()).hexdigest()},
        "inputs": {},
        "queries": [],  # one per query, datatable and matrix block
        "sql": [],
        "trees": [],
        "blocks": [],  # one per block in file order: what the build summary prints
        "warnings": [],
        "environment": {
            "groundie": pkg_version("groundie"),
            "duckdb": duckdb.__version__,
            "python": platform.python_version(),
        },
        "built_at": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
    }
    con = duckdb.connect()
    for seg in segments:  # every input first, so any block can read any input
        if seg["type"] == "input":
            with at(f"{path.name}:{seg['line']}"):
                manifest["inputs"][seg["block"].name] = load_input(con, seg["block"], base_dir)
    # Inputs are in. From here no block may touch files, extensions or the network, or undo that.
    con.execute("set enable_external_access = false")
    con.execute("set lock_configuration = true")

    def warn(line, message):
        manifest["warnings"].append(f"{path.name}:{line}: {message}")

    parts = []
    for seg in segments:
        if seg["type"] == "prose":
            for line, message in unparsed_fences(seg):
                warn(line, message)
            parts.append(render_prose(seg["text"]))
            continue
        with at(f"{path.name}:{seg['line']}"):
            summary = build_block(con, blocks, manifest, parts, seg)
        manifest["blocks"].append({"line": seg["line"], "kind": seg["type"], **summary})
        if summary.get("rows") == 0:
            filters = isinstance(seg["block"], q.Query) and q.filter_text(seg["block"])
            warn(seg["line"], f"{seg['type']} {summary['name']!r}: 0 rows"
                 + (f"; check the filter values: {filters}" if filters else ""))
    return "".join(parts), manifest, con


def build_block(con, blocks, manifest, parts, seg):
    """Run one block, append its HTML to parts and its entries to the manifest. Returns its summary: name,
    rows and output columns, or for a tree, tips and matched tips."""
    block, who = seg["block"], f"{seg['type']} {seg['block'].name!r}"
    if seg["type"] == "Tree":
        figure, notes, entry = tree.render(con, blocks, block)
        manifest["trees"].append(entry)
        chain = [blocks[block.source], *(q.lineage(blocks, block.tips) if block.tips else [])]
        parts.append(framed("tree", figure, q.dump(entry["block"]), chain, notes))
        return {"name": block.name, "tips": entry["tips"], "matched": entry["matched"]}
    summary = lambda rows, columns: {"name": block.name, "rows": rows, "output_columns": columns}
    if seg["type"] == "input":  # loaded before the lock
        return summary(*size(con, block.name))
    if seg["type"] == "clusters":
        make_clusters(con, block)
        rows, columns = size(con, block.name)
        manifest["queries"].append({"kind": "clusters", **q.as_dict(block), "rows": rows})
        return summary(rows, columns)
    if seg["type"] in ("checks", "calls"):
        sql, params = (make_checks(con, block), []) if seg["type"] == "checks" else make_calls(con, block)
        rows, columns = size(con, block.name)
        manifest["queries"].append({"kind": seg["type"], **q.as_dict(block), "sql": sql, "params": params, "rows": rows})
        return summary(rows, columns)
    if seg["type"] == "sql":
        create(con, who, block.name, block.sql)
        rows, columns = size(con, block.name)
        manifest["sql"].append({"name": block.name, "sql": block.sql, "rows": rows})
        return summary(rows, columns)
    sql, params = q.to_sql(block)
    create(con, who, block.name, sql, params)
    rows, columns = size(con, block.name)
    entry = {"kind": seg["type"], **q.as_dict(block), "sql": sql, "params": params, "rows": rows}
    manifest["queries"].append(entry)
    if seg["type"] != "query":
        data = con.execute(f"select * from {ident(block.name)}").fetchall()
    if seg["type"] == "DataTable":
        written = seg["display"]
        if isinstance(written, str):  # an input: its adapter's titles, for the columns that have one
            declared = ADAPTERS[blocks[written].adapter].TITLES
            titles = {c: declared[c] for c in columns if c in declared}
        else:
            titles = written
            outside = [c for c in titles if c not in columns]
            if outside:
                raise SpecError(f"{who}: display names {', '.join(map(repr, outside))}, not an output column; "
                                f"output columns: {', '.join(columns)}")
        if titles:
            entry["display"] = titles
        parts.append(framed("datatable", render_datatable([(c, titles.get(c, c)) for c in columns], data),
                            q.to_yaml(block) + (q.dump({"display": written}) if written else ""),
                            q.lineage(blocks, block.source)))
    if seg["type"] == "MlstTable":
        split = tables.mlst_tables(block.name, columns, data, bool(block.sort))
        entry["tables"] = [{"caption": t.caption, "columns": [c for c, _ in t.columns], "rows": len(t.rows)}
                           for t in split]
        body = "".join(render_datatable(t.columns, t.rows, t.caption) for t in split) or render_datatable(
            [(c, c) for c in columns], [])
        parts.append(framed("datatable", body, q.to_yaml(block), q.lineage(blocks, block.source),
                            () if block.sort else (tables.DEFAULT_ORDER,)))
        return {**summary(rows, columns), "tables": len(split)}
    if seg["type"] in MATRICES:
        entry["layout"] = vars(seg["layout"])
        grid_yaml = q.dump({k: v for k, v in vars(seg["layout"]).items() if v not in (None, False)})
        parts.append(framed("matrix", render_matrix(block, seg["layout"], columns, data),
                            q.to_yaml(block) + grid_yaml, q.lineage(blocks, block.source),
                            () if block.sort else (DEFAULT_ORDER,)))
    return summary(rows, columns)


GROUNDED = ('<span class="grounded"><b>Grounded data</b>: this block is rendered deterministically from the '
            "source files.</span>")


def framed(kind, body, spec_yaml, chain, notes=()):
    """A display block in its bordered frame: the body, then the grounded note and 'View query'."""
    return f'<section class="block {kind}">\n{body}{details(spec_yaml, chain, notes)}</section>\n'


def summary_line(name, block):
    """One block of the build summary, e.g. 'cohort.md:36: DataTable st: 10 rows; columns: sample, scheme'."""
    where = f"{name}:{block['line']}: {block['kind']} {block['name']}: "
    if block["kind"] == "Tree":
        matched = "" if block["matched"] is None else f", {block['matched']} matched"
        return where + q.plural(block["tips"], "tip") + matched
    split = f" in {q.plural(block['tables'], 'table')}" if "tables" in block else ""
    return where + f"{q.plural(block['rows'], 'row')}{split}; columns: {', '.join(block['output_columns'])}"


def one_line(message):
    """A message on one line: DuckDB's run over several, with a caret line under the SQL."""
    return " ".join(l.strip() for l in str(message).splitlines() if l.strip(" ^"))


def write(path, text):
    try:
        path.write_text(text)
    except OSError as e:
        raise SpecError(f"cannot write {str(path)!r}: {e.strerror or e}") from e


def main():
    parser = argparse.ArgumentParser(prog="groundie")
    sub = parser.add_subparsers(dest="command", required=True)
    b = sub.add_parser("build", help="build an HTML report from a Markdown file of blocks")
    b.add_argument("file", type=Path)
    b.add_argument("-o", "--out", type=Path, required=True)
    b.add_argument("--manifest", type=Path, help="also write the build's manifest, as JSON, to this path")
    d = sub.add_parser("dev", help="serve the report of a file of blocks, reloading on changes")
    d.add_argument("file", type=Path)
    d.add_argument("--port", type=int, default=8000)
    args = parser.parse_args()
    if args.command == "dev":
        from groundie import dev
        return dev.serve(args.file, args.port)
    try:
        body, manifest, _ = build_report(args.file)
        write(args.out, page(args.file.stem, body + render_provenance(manifest)))
        if args.manifest:
            write(args.manifest, json.dumps(manifest, indent=2) + "\n")
    except SpecError as e:
        print(f"groundie: error: {one_line(e)}", file=sys.stderr)
        sys.exit(1)
    for warning in manifest["warnings"]:
        print(f"groundie: warning: {warning}", file=sys.stderr)
    print(f"built {args.out} from {args.file}")
    for block in manifest["blocks"]:
        print(summary_line(manifest["file"]["path"], block))
