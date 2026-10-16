"""Tables as data, before any markup: a display block's rows split into several tables, or a list column
spread out."""
from groundie.query import SpecError


def split(who, columns, rows, by):
    """The rows grouped by the value of one column, as {value: rows}: groups, and the rows in each, in
    first-seen order. A null is a value like any other."""
    if by not in columns:
        raise SpecError(f"{who}: split names {by!r}, not an output column; output columns: {', '.join(columns)}")
    at = columns.index(by)
    groups = {}
    for n, row in enumerate(rows, 1):
        if isinstance(row[at], (list, tuple, dict)):  # DuckDB gives a fixed-size ARRAY as a tuple
            raise SpecError(f"{who}: split {by!r}: row {n} holds {row[at]!r}, not a single value")
        groups.setdefault(row[at], []).append(row)
    return groups


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
