"""Tables as data, before any markup: what a display block shows when it shows several tables."""
import re
from dataclasses import dataclass

from groundie.adapters.mlst import LOCI
from groundie.query import SpecError

LOCUS = re.compile(r"(.+)\((.*)\)")  # mlst's 'arcC(3)': the locus name, then the allele with its flags
LABELS = {"sample": "Sample", "sequence_type": "ST"}
NO_SCHEME = "-"  # mlst's scheme for a genome it could not type
DEFAULT_ORDER = ("Tables by how many genomes they hold, most first, then by scheme; rows by sample. "
                 "Give 'sort' to order both by the query instead.")


@dataclass
class Table:
    """One table: a caption, its columns as (name, label), where name is the column in the query's rows and
    label the header shown, and its rows, one value per column."""
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


def names_from_cells(who, columns, rows, parse, keys):
    """Move a name shared by every cell of a column into the column's header: the columns as (name, label)
    and the rows with only what is left of each cell. parse splits a cell into (name, rest), or gives None.
    A column whose cells are all null is dropped. Cells naming different names are an error, which names
    the rows by keys, one per row."""
    keep, labels = [], {}
    for i, column in enumerate(columns):
        names = {}
        for row, key in zip(rows, keys):
            if row[i] is None:
                continue
            parsed = parse(row[i])
            if parsed is None:
                raise SpecError(f"{who}: {column} of {key} is {row[i]!r}, not 'name(value)'")
            names.setdefault(parsed[0], f"{row[i]!r} in {key}")
        if len(names) > 1:
            raise SpecError(f"{who}: {column} names different loci: {', '.join(names.values())}")
        if names:
            keep.append(i)
            labels[i] = next(iter(names))
    cells = lambda row: tuple(row[i] if row[i] is None else parse(row[i])[1] for i in keep)
    return [(columns[i], labels[i]) for i in keep], [cells(row) for row in rows]


def mlst_tables(name, columns, rows, ordered):
    """An mlst_table's tables: one per scheme, captioned with it, the locus names as headers and the alleles
    in the cells. Genomes with no scheme found come last. Without ordered, rows go by sample."""
    groups = split(columns, rows, "scheme", ordered)
    if NO_SCHEME in groups:
        groups[NO_SCHEME] = groups.pop(NO_SCHEME)
    fixed = ["sample", "sequence_type"]
    at = [columns.index(c) for c in [*fixed, *LOCI]]
    tables = []
    for scheme, group in groups.items():
        group = [tuple(row[i] for i in at) for row in group]  # sample, sequence_type, then the loci
        if not ordered:
            group.sort(key=lambda row: by_name(row[0]))
        who = f"mlst_table {name!r}: scheme {scheme!r}"
        loci, cells = names_from_cells(who, LOCI, [row[2:] for row in group], parse_locus, [row[0] for row in group])
        tables.append(Table(
            "No scheme found" if scheme == NO_SCHEME else scheme,
            [(c, LABELS[c]) for c in fixed] + loci,
            [row[:2] + rest for row, rest in zip(group, cells)],
        ))
    return tables


def parse_locus(text):
    m = LOCUS.fullmatch(text)
    return m.groups() if m else None
