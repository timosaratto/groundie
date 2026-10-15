"""Tables as data, before any markup: what a display block shows when it shows several tables, or a list
column spread out."""
from dataclasses import dataclass

from groundie.query import SpecError

LABELS = {"sample": "Sample", "sequence_type": "ST"}
DEFAULT_ORDER = ("Tables by how many genomes they hold, most first, then by scheme; rows by sample. "
                 "Give 'sort' to order both by the query instead.")


@dataclass
class Table:
    """One table: a caption, its columns as (name, label), where name is the column in the query's rows (or
    what stands for one, such as an mlst locus) and label the header shown, and its rows, one value per
    column."""
    caption: str
    columns: list[tuple[str, str]]
    rows: list[tuple]


def by_name(value):
    return (value is None, value)


def split(columns, rows, by, ordered):
    """The rows grouped by the value of one column, as {value: rows}. With ordered (the query has a sort),
    groups and their rows keep the query's order; without, groups go largest first, then by value."""
    at = columns.index(by)
    groups = {}
    for row in rows:
        groups.setdefault(row[at], []).append(row)
    if ordered:
        return groups
    return {k: groups[k] for k in sorted(groups, key=lambda k: (-len(groups[k]), by_name(k)))}


def mlst_tables(name, columns, rows, ordered):
    """An MlstTable's tables: one per scheme, captioned with it, its loci as headers and the calls in the
    cells. Genomes with no scheme found (a null scheme) come last. Without ordered, rows go by sample."""
    groups = split(columns, rows, "scheme", ordered)
    if None in groups:
        groups[None] = groups.pop(None)
    fixed = ["sample", "sequence_type"]
    at = [columns.index(c) for c in [*fixed, "calls"]]
    tables = []
    for scheme, group in groups.items():
        group = [tuple(row[i] for i in at) for row in group]  # sample, sequence_type, calls
        if not ordered:
            group.sort(key=lambda row: by_name(row[0]))
        loci = scheme_loci(f"MlstTable {name!r}: scheme {scheme!r}", group)
        tables.append(Table(
            "No scheme found" if scheme is None else scheme,
            [(c, LABELS[c]) for c in fixed] + [(locus, locus) for locus in loci],
            [(sample, st, *(c["call"] for c in calls)) for sample, st, calls in group],
        ))
    return tables


class Missing:
    """A spread cell whose row's list has no entry for its column: blank, unlike a null value."""
    def __repr__(self):
        return "MISSING"


MISSING = Missing()


def spread(who, columns, rows, column, names, values):
    """columns and rows with a list column spread out: in its place, one column per distinct value of its
    structs' 'names' field, in first-seen order (row order, then list order), each holding the 'values' field
    of the struct with that name, or MISSING where the row's list has none. A null list gives nulls. Returns
    (columns, rows); a spread column is (name, name)."""
    at = [c for c, _ in columns]
    if column not in at:
        raise SpecError(f"{who}: spread names {column!r}, not an output column; output columns: {', '.join(at)}")
    i = at.index(column)
    shape = f"a list of {{{names}, {values}}}"
    found, spread_names = [], {}
    for n, row in enumerate(rows, 1):
        if row[i] is None:
            found.append(None)
            continue
        if not (isinstance(row[i], (list, tuple)) and all(  # DuckDB gives a fixed-size ARRAY as a tuple
                isinstance(e, dict) and names in e and values in e and not isinstance(e[names], (list, dict))
                for e in row[i])):
            raise SpecError(f"{who}: spread {column!r}: row {n} holds {row[i]!r}, not {shape}")
        entries = {}
        for e in row[i]:
            if e[names] is None:
                raise SpecError(f"{who}: spread {column!r}: row {n} has a null {names}")
            if e[names] in entries:
                raise SpecError(f"{who}: spread {column!r}: row {n} has {names} {e[names]!r} twice")
            entries[e[names]] = e[values]
            spread_names.setdefault(e[names], None)
        found.append(entries)
    others = at[:i] + at[i + 1:]
    taken = [str(k) for k in spread_names if str(k) in others]
    if taken:
        raise SpecError(f"{who}: spread {column!r}: {names} {taken[0]!r} is the name of another column; "
                        f"other columns: {', '.join(others)}")
    blank = (None,) * len(spread_names)
    return ([*columns[:i], *((str(k), str(k)) for k in spread_names), *columns[i + 1:]],
            [(*row[:i], *(blank if entries is None else (entries.get(k, MISSING) for k in spread_names)),
              *row[i + 1:]) for row, entries in zip(rows, found)])


def scheme_loci(who, group):
    """The loci every genome of one scheme calls, in order. Genomes calling different loci are an error."""
    loci = {}
    for sample, _, calls in group:
        loci.setdefault(tuple(c["locus"] for c in calls), sample)
    if len(loci) > 1:
        found = "; ".join(f"{', '.join(k) or 'none'} in {sample}" for k, sample in loci.items())
        raise SpecError(f"{who}: genomes call different loci: {found}")
    return next(iter(loci))
