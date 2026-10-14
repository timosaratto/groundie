"""The blocks that make SQL tables, independent of SQL, and how the page describes them."""

import re
from dataclasses import dataclass, field
from pathlib import PurePath

import yaml

SHOWN = 5  # names listed in a message before "and N more"


class SpecError(ValueError):
    """Something wrong with the blocks or their files: reported to the author in one line, without a traceback."""


@dataclass
class Query:
    """What a query, datatable or matrix block asks for."""
    name: str
    source: str  # the block read from
    columns: list[str] = field(default_factory=list)  # empty means every column
    filters: dict[str, list[str | None]] = field(default_factory=dict)  # column -> accepted values; None: null
    group_by: list[str] = field(default_factory=list)  # set: one row per distinct combination of these
    measures: dict[str, str] = field(default_factory=dict)  # output alias -> fn(column), per group
    sort: list[str] = field(default_factory=list)


@dataclass
class Sql:
    """An ad hoc ```sql block: one SELECT statement, as written."""
    name: str
    sql: str


@dataclass
class Input:
    """An input block: files read by one adapter into one SQL table. Paths are as written."""
    name: str
    adapter: str
    paths: list[str]
    options: dict = field(default_factory=dict)


@dataclass
class Clusters:
    """A clusters block: single-linkage clusters of the samples in a table of pairwise distances."""
    name: str
    source: str  # the block read from: one row per pair
    threshold: float  # pairs at this distance or closer are linked
    pairs: list[str] = field(default_factory=lambda: ["reference", "query"])  # the two sample columns
    distance: str = "distance"
    linkage: str = "single"


@dataclass
class Checks:
    """A checks block: each row of its source passes or fails thresholds on its columns."""
    name: str
    source: str
    checks: dict[str, dict[str, float]]  # column -> {bound: number}, bound in min, max (inclusive), above, below
    key: str = "sample"  # the column identifying a row


@dataclass
class Calls:
    """A calls block: per sample, whether its rows hold any, or all, of some values of one column."""
    name: str
    source: str
    calls: dict[str, dict[str, list[str]]]  # output column -> {'any' or 'all': values}
    samples: str | None = None  # a block listing every sample in its key column; none means the samples in source
    column: str = "element_symbol"
    key: str = "sample"


QUERY_KEYS = ("source", "columns", "filters", "group_by", "measures", "sort")
MEASURE = re.compile(r"(count)\((\*)\)|(count|min|max|sum|avg|list|distinct)\((\w+)\)")
FUNCTIONS = "count, min, max, sum, avg, list, distinct"


def measure(who, alias, text):
    """Parse 'fn(column)' into (fn, column), or raise. 'who' names the block in the message: "datatable 'x'"."""
    m = MEASURE.fullmatch(str(text))
    if not m:
        raise SpecError(f"{who}: measure {alias!r}: expected fn(column) with fn in {FUNCTIONS}")
    return m.group(1, 2) if m.group(1) else m.group(3, 4)


def column_list(who, block, key):
    """A key naming columns: one column or a list of them. Returns the list."""
    value = block.get(key)
    if value is None:
        return []
    if isinstance(value, str):
        value = [value]
    if not (isinstance(value, list) and all(isinstance(c, str) for c in value)):
        raise SpecError(f"{who}: {key!r} must be one column or a list of columns")
    return value


def normalize(name, block, component="datatable", taken=()):
    """Turn a parsed block into a Query: scalars become lists, values become strings, columns dedupe.
    'component' is the block's kind (query, datatable, matrix), for messages. 'taken' are keys the caller
    already took out of the block, such as a matrix's layout: accepted, but not part of the query."""
    who = f"{component} {name!r}"
    unknown = [k for k in block if k not in QUERY_KEYS]
    if unknown:
        raise SpecError(f"{who}: unknown key {', '.join(map(repr, unknown))}; "
                         f"accepted: {', '.join((*QUERY_KEYS, *taken))}")
    if not isinstance(block.get("source"), str):
        raise SpecError(f"{who}: needs 'source', the name of the block it reads")
    if "columns" in block and "group_by" in block:
        raise SpecError(f"{who}: give 'columns' or 'group_by', not both")
    columns = list(dict.fromkeys(column_list(who, block, "columns")))
    group_by = list(dict.fromkeys(column_list(who, block, "group_by")))
    filters = {}
    for column, values in (block.get("filters") or {}).items():
        if not isinstance(values, list):
            values = [values]
        if not values:
            raise SpecError(f"{who}: filter {column!r} needs at least one value")
        # YAML types bare values (007 -> 7, yes -> True); filter values are always strings, but null (or ~)
        # stays None and matches a null.
        filters[column] = [None if v is None else str(v) for v in values]
    measures = {str(alias): str(text) for alias, text in (block.get("measures") or {}).items()}
    if measures and not group_by:
        raise SpecError(f"{who}: 'measures' needs 'group_by'")
    for alias, text in measures.items():
        measure(who, alias, text)
        if alias in group_by:
            raise SpecError(f"{who}: measure {alias!r} has the name of a group_by column")
    query = Query(name, block["source"], columns, filters, group_by, measures, column_list(who, block, "sort"))
    outside = [c for c in query.sort if c not in outputs(query)]
    if group_by and outside:
        raise SpecError(f"{who}: sort names {', '.join(map(repr, outside))}, not an output column; "
                         f"output columns: {', '.join(outputs(query))}")
    return query


def outputs(q):
    """The output columns, when the query names them; empty means every column of its source."""
    return q.group_by + list(q.measures) or q.columns


def as_dict(block):
    """A query or input as its block writes it, with its name."""
    if isinstance(block, Input):
        files = {"path": block.paths[0]} if len(block.paths) == 1 else {"paths": block.paths}
        return {"name": block.name, "adapter": block.adapter, **files, "options": block.options}
    if isinstance(block, Clusters):
        return {"name": block.name, "source": block.source, "pairs": block.pairs, "distance": block.distance,
                "threshold": block.threshold, "linkage": block.linkage}
    if isinstance(block, Checks):
        return {"name": block.name, "source": block.source, "key": block.key, "checks": block.checks}
    if isinstance(block, Calls):
        return {"name": block.name, "source": block.source, "samples": block.samples, "column": block.column,
                "key": block.key, "calls": block.calls}
    return {"name": block.name, "source": block.source, "filters": block.filters, "columns": block.columns,
            "group_by": block.group_by, "measures": block.measures, "sort": block.sort}


def ident(name):
    """Quote a SQL identifier."""
    return '"' + str(name).replace('"', '""') + '"'


def to_sql(q):
    """Return (sql, params) for DuckDB."""
    columns = [ident(c) for c in q.group_by or q.columns]
    for alias, text in q.measures.items():
        fn, column = measure(repr(q.name), alias, text)
        if fn == "list":  # every value, sorted; as text, so numbers join too
            columns.append(f"string_agg(cast({ident(column)} as varchar), ', ' order by {ident(column)}) as {ident(alias)}")
        elif fn == "distinct":  # string_agg(distinct ...) can only sort as text: 10 before 9
            columns.append(f"array_to_string(nullif(list_sort(list_distinct(list({ident(column)}))), []), ', ') "
                           f"as {ident(alias)}")
        else:
            columns.append(f"{fn}({'*' if column == '*' else ident(column)}) as {ident(alias)}")
    columns = ", ".join(columns) or "*"
    sql = f"select {columns} from {ident(q.source)}"
    params = []
    if q.filters:
        clauses = []
        for column, values in q.filters.items():
            known = [v for v in values if v is not None]
            tests = [f"{ident(column)} in ({', '.join('?' * len(known))})"] if known else []
            if None in values:  # 'in' never matches a null
                tests.append(f"{ident(column)} is null")
            clauses.append(tests[0] if len(tests) == 1 else f"({' or '.join(tests)})")
            params.extend(known)
        sql += " where " + " and ".join(clauses)
    if q.group_by:
        sql += " group by " + ", ".join(ident(c) for c in q.group_by)
    if q.sort:
        sql += " order by " + ", ".join(ident(c) for c in q.sort)
    return sql, params


BOUNDS = {"min": ">=", "max": "<=", "above": ">", "below": "<"}


def checks_sql(block):
    """The SQL of a checks block: the key, each checked column, qc (PASS or FAIL) and qc_failed, the failed
    checks with their values in the order written, or null. A null value fails its check."""
    passes, failed = [], []
    for column, bounds in block.checks.items():
        test = " and ".join(f"{ident(column)} {BOUNDS[b]} {n!r}" for b, n in bounds.items())
        passes.append(f"coalesce({test}, false)")
        failed.append(f"case when not coalesce({test}, false) then {literal(column + ' ')} || "
                      f"coalesce(cast({ident(column)} as varchar), 'missing') end")
    columns = ", ".join(ident(c) for c in [block.key, *block.checks])
    return (f"select {columns},\n  case when {' and '.join(passes)} then 'PASS' else 'FAIL' end as qc,\n"
            f"  nullif(concat_ws('; ', {', '.join(failed)}), '') as qc_failed\nfrom {ident(block.source)}")


def calls_sql(block):
    """Return (sql, params) of a calls block: one row per sample, in key order, then a boolean per call. A
    sample with no rows in source gets false."""
    key, column, params, calls = ident(block.key), ident(block.column), [], []
    for name, rule in block.calls.items():
        (how, values), = rule.items()
        values = list(dict.fromkeys(values))
        params.extend(values)
        within = f"{column} in ({', '.join('?' * len(values))})"
        test = "> 0" if how == "any" else f"= {len(values)}"
        calls.append(f"count(distinct {column}) filter (where {within}) {test} as {ident(name)}")
    samples = ident(block.samples or block.source)
    return (f"select s.{key}, {', '.join(calls)}\n"
            f"from (select distinct {key} from {samples} where {key} is not null) s\n"
            f"left join {ident(block.source)} h using ({key})\ngroup by s.{key}\norder by s.{key}"), params


def literal(text):
    """Quote a SQL string."""
    return "'" + str(text).replace("'", "''") + "'"


def lineage(blocks, name):
    """The block 'name' and the blocks behind it, nearest first: each query's source in turn, down to an
    input or a sql block. 'blocks' maps names to blocks."""
    chain = []
    while name is not None:
        chain.append(blocks[name])
        name = getattr(chain[-1], "source", None)
    return chain


def about(block):
    """What a block is, in one phrase: 'hits: filters element_type = AMR', 'summary: ad hoc SQL', or
    'amrfinder: adapter amrfinder, 2 files: a.tsv, b.tsv'."""
    if isinstance(block, Sql):
        return f"{block.name}: ad hoc SQL"
    if isinstance(block, Input):
        names = [PurePath(p).name for p in block.paths]
        files = f"file {names[0]}" if len(names) == 1 else f"{len(names)} files: {listing(names)}"
        return f"{block.name}: adapter {block.adapter}, {files}"
    if isinstance(block, Clusters):
        return f"{block.name}: {block.linkage}-linkage clusters of {block.source} at {block.distance} ≤ {block.threshold}"
    if isinstance(block, Checks):
        return f"{block.name}: checks of {block.source} on {', '.join(block.checks)}"
    if isinstance(block, Calls):
        return f"{block.name}: calls {', '.join(block.calls)} from {block.source} {block.column}"
    return f"{block.name}: " + ("filters " + filters_and_shape(block) if block.filters else "no filters" + shape(block))


def filters_and_shape(q):
    return (filter_text(q) or "none") + shape(q)


def filter_text(q):
    """The filters as one phrase, e.g. 'element_type = AMR, sample in A, B'; empty when there are none."""
    shown = {column: ["null" if v is None else v for v in values] for column, values in q.filters.items()}
    return ", ".join(f"{column} = {values[0]}" if len(values) == 1 else f"{column} in {', '.join(values)}"
                     for column, values in shown.items())


def shape(q):
    """What follows the filters: the grouping, the measures, then the sort."""
    out = "; one row per " + ", ".join(q.group_by) if q.group_by else ""
    if q.measures:
        out += "; " + ", ".join(f"{alias} = {text}" for alias, text in q.measures.items())
    return out + ("; sorted by " + ", ".join(q.sort) if q.sort else "")


class Dumper(yaml.SafeDumper):
    """Lists on one line, maps as blocks: 'element_type: [AMR]', one measure per line."""


Dumper.add_representer(list, lambda d, v: d.represent_sequence("tag:yaml.org,2002:seq", v, flow_style=True))


def dump(d):
    return yaml.dump(d, Dumper=Dumper, sort_keys=False, allow_unicode=True)


def to_yaml(block):
    """A query or input expanded, valid as its block's body."""
    return dump({k: v for k, v in as_dict(block).items() if k != "name" and v not in (None, "", [], {})})


def plural(n, word):
    return f"{n} {word}{'' if n == 1 else 's'}"


def listing(names):
    """Up to SHOWN names, then how many more."""
    more = len(names) - SHOWN
    return ", ".join(names[:SHOWN]) + (f" and {more} more" if more > 0 else "")
