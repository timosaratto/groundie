"""Adapters declared as classes: the file's format once, then one field per column of the SQL table, saying
which column of the file it comes from, its type, and what it means.

    class Pair(Adapter):
        \"\"\"One pair of genomes.\"\"\"

        FORMAT = Tsv(header=False, null="NA")

        reference: str = column(1, transform=extract_sample_name)
        "The first genome."

        distance: float | None = column(3)
        "The distance; null where the tool wrote NA."

A field reads one cell, column(...), or a run of cells into a list, columns(...). It sees nothing else: no
other field, no other line. So a declaration is plain data, a position, a type and a transform, that could as
well be compiled to a DuckDB select; joins and logic belong in the report's sql blocks. A cell is read in
three steps:

1. The format's null token is null, which only a type with '| None' allows.
2. The transform, if any, turns the text into other text: extract_sample_name strips directories and FASTA
   extensions.
3. The field's type reads the text: str as is, int and float strictly, and a Struct by its PATTERN, such as
   '{locus}({call})', which splits the text into the Struct's fields, each read as a cell is, from step 1.

A name matches the header ignoring case, as in DuckDB; a UTF-8 byte order mark before the header is dropped.
A column named in the header can say what to read where a file lacks it: column("Name",
absent=sample_from_file_name). That too is one value per file, from its path, as DuckDB's
read_csv(filename = true) gives; so two files of an input that give the same value, other than null, fail it.

A field may also have a title, column(3, title="ST"): the header a datatable shows for it when its display
names the input, 'display: mlst'. A column named in the header is titled by that name unless it says
otherwise.

The class gives the SQL COLUMNS, the JSON SCHEMA, the TITLES and the rows that cli.load_input expects."""

import ast
import inspect
import re
import textwrap
import types
from dataclasses import dataclass
from pathlib import Path
from typing import NamedTuple, get_args, get_origin, get_type_hints

FASTA = re.compile(r"\.(fa|fasta|fna|fas)(\.gz)?$")
INTEGER = re.compile(r"[+-]?\d+")
SCALARS = {str: ("string", "varchar"), int: ("integer", "bigint"), float: ("number", "double")}


# Transforms: text in, text out.

def extract_sample_name(text):
    """A genome's file as its sample name: the directories and a FASTA extension stripped, as in mash.SAMPLE."""
    return FASTA.sub("", text.rsplit("/", 1)[-1])


# Absents: a file's path in, the text to read in place of a column the file lacks out; None for null.

def no_value(path):
    """Null, on every line of the file."""
    return None


def sample_from_file_name(path):
    """The file's name without its directories and last extension: 'results/KP01.tsv' is 'KP01'."""
    return Path(path).stem


# Readers, one per type: text in, value out; a ValueError says what the text is not.

def integer(text):
    if not INTEGER.fullmatch(text):
        raise ValueError(f"{text!r} is not an integer")
    return int(text)


def number(text):
    try:
        return float(text)
    except ValueError:
        raise ValueError(f"{text!r} is not a number") from None


READERS = {str: str, int: integer, float: number}


@dataclass(frozen=True)
class Tsv:
    """A tab-separated file. header: whether its first line names the columns. null: the text that stands for
    no value, in any cell; None for no such text."""
    header: bool
    null: str | None = None
    delimiter = "\t"


class Pattern:
    """A cell's parts by name: '{locus}({call})' splits 'tpi(~4)' into locus 'tpi' and call '~4'. Text
    outside braces must appear as written; each part takes as little text as it can."""

    def __init__(self, text):
        parts = re.split(r"\{(\w+)\}", text)  # literal, name, literal, ..., literal
        self.names = parts[1::2]
        self.shown = "".join(parts)
        self.regex = re.compile("".join(f"(?P<{p}>.*?)" if i % 2 else re.escape(p) for i, p in enumerate(parts)))

    def split(self, text):
        m = self.regex.fullmatch(text)
        if not m:
            raise ValueError(f"{text!r} is not {self.shown!r}")
        return m.groupdict()


@dataclass(frozen=True)
class Source:
    """Where a field's text comes from, and how it is read before its type reads it; and its title, the
    header a datatable with 'display: <input>' shows for the field."""
    first: int | str  # a position, counted from 1, or a name in the header
    last: int | None  # None: the one cell at first; a position, negative from the line's end: a list of cells
    transform: object
    title: str | None = None
    absent: object = None  # where the header lacks a named column: a function of the file's path; None: an error

    @property
    def many(self):
        return self.last is not None


def column(at, transform=None, title=None, absent=None):
    """One cell: the column at a position, counted from 1, or of a name in the header. absent: where the
    header lacks the name, a function of the file's path giving the text to read in its place, such as
    no_value; by default it's an error."""
    return Source(at, None, transform, title, absent)


def columns(first, last, transform=None, title=None):
    """A list of cells: the columns from position first through last, both counted from 1, as in DuckDB's
    list_slice. A negative last counts from each line's end: -1 is its last cell, -2 the one before."""
    return Source(first, last, transform, title)


class Absent(NamedTuple):
    """The value, the same on every line, of a column the file lacks."""
    value: object


class Field(NamedTuple):
    hint: object
    source: Source | None
    doc: str | None

    @property
    def title(self):
        """The source's title; by default, a column named in the header is titled by that name."""
        if not isinstance(self.source, Source):
            return None
        if self.source.title is None and isinstance(self.source.first, str):
            return self.source.first
        return self.source.title


class Struct:
    """A cell of several named, typed, documented parts, which its PATTERN splits it into."""

    def __init_subclass__(cls):
        super().__init_subclass__()
        cls.FIELDS = declared_fields(cls)
        for name, f in cls.FIELDS.items():
            if optional(f.hint)[0] not in READERS:
                raise TypeError(f"{cls.__name__}.{name}: a part is str, int or float, not {f.hint}")
        if not isinstance(vars(cls).get("PATTERN"), str):
            raise TypeError(f"{cls.__name__}: needs a PATTERN naming its parts, such as '{{a}}:{{b}}'")
        cls.parts = Pattern(cls.PATTERN)
        if cls.parts.names != list(cls.FIELDS):
            raise TypeError(f"{cls.__name__}: its PATTERN has parts {', '.join(cls.parts.names)}; "
                            f"its fields are {', '.join(cls.FIELDS)}")

    @classmethod
    def read(cls, text, null):
        out = {}
        for name, part in cls.parts.split(text).items():
            try:
                out[name] = read(part, cls.FIELDS[name].hint, null)
            except ValueError as e:
                raise ValueError(f"{name}: {e}") from None
        return out


class Adapter:
    """A file format and the columns of the SQL table read from it; see the module's docstring."""

    VERSION = "preliminary"
    OPTIONS = ()

    def __init_subclass__(cls):
        super().__init_subclass__()
        cls.NAME = cls.__module__.rsplit(".", 1)[-1]
        cls.FIELDS = declared_fields(cls)
        for name, f in cls.FIELDS.items():
            if problem := source_problem(cls.FORMAT, f):
                raise TypeError(f"{cls.__name__}.{name}: {problem}")
        cls.COLUMNS = tuple((name, sql_type(f.hint)) for name, f in cls.FIELDS.items())
        cls.TITLES = {name: f.title for name, f in cls.FIELDS.items() if f.title}
        cls.SCHEMA = {"$schema": "https://json-schema.org/draft/2020-12/schema", "title": f"{cls.NAME} result",
                      **object_type(cls)}

    @classmethod
    def rows(cls, path, options):
        path = Path(path)
        lines = [(n, line.split(cls.FORMAT.delimiter))
                 for n, line in enumerate(path.read_text(encoding="utf-8-sig").splitlines(), 1) if line.strip()]
        if cls.FORMAT.header and not lines:  # not even a header: no rows
            return []
        header = lines.pop(0)[1] if cls.FORMAT.header else None
        try:
            where = cls.cells(header, path)
        except ValueError as e:
            raise ValueError(f"{path.name}: {e}") from None
        out = []
        for n, cells in lines:
            try:
                out.append(cls.row(cells, where, header))
            except ValueError as e:
                raise ValueError(f"{path.name}, line {n}: {e}") from None
        return out

    @classmethod
    def path_values(cls, path):
        """What the file's path gives in place of the columns it lacks, by field, nulls left out. Like
        read_csv(filename = true), each tells the file apart: two files of an input mustn't give the same."""
        if not cls.FORMAT.header:
            return {}
        with open(path, encoding="utf-8-sig") as f:
            header = next((line.rstrip("\r\n").split(cls.FORMAT.delimiter) for line in f if line.strip()), None)
        if header is None:
            return {}
        return {name: at.value for name, at in zip(cls.FIELDS, cls.cells(header, path))
                if isinstance(at, Absent) and at.value is not None}

    @classmethod
    def cells(cls, header, path):
        """Each field's cells in a line: an index, a slice for a list, or the Absent value of a column the
        file lacks. A name matches the header ignoring case, as DuckDB's identifiers do."""
        folded = [h.casefold() for h in header] if header is not None else None
        out = []
        for name, f in cls.FIELDS.items():
            s = f.source
            if isinstance(s.first, str):
                if s.first.casefold() in folded:
                    out.append(folded.index(s.first.casefold()))
                elif s.absent:
                    try:
                        out.append(Absent(read(s.absent(path), f.hint, None)))
                    except ValueError as e:
                        raise ValueError(f"{name}: {e}") from None
                else:
                    raise ValueError(f"no column {s.first!r}; the header has {', '.join(map(repr, header))}")
            else:
                out.append(slice(s.first - 1, s.last if s.last > 0 else s.last + 1 or None) if s.many else s.first - 1)
        return out

    @classmethod
    def row(cls, cells, where, header):
        if header is not None:
            if len(cells) != len(header):
                raise ValueError(f"has {len(cells)} columns; the header has {len(header)}")
        else:
            sources = [f.source for f in cls.FIELDS.values()]
            open_ended = any(s.many and s.last < 0 for s in sources)
            needed = max(s.first - 2 - s.last if s.many and s.last < 0 else s.last or s.first for s in sources)
            if len(cells) < needed or (len(cells) > needed and not open_ended):
                raise ValueError(f"has {len(cells)} columns; needs {'at least ' if open_ended else ''}{needed}")
        out = []
        for (name, f), at in zip(cls.FIELDS.items(), where):
            try:
                if isinstance(at, Absent):
                    out.append(at.value)
                elif f.source.many:
                    out.append([read(text, get_args(f.hint)[0], cls.FORMAT.null, f.source.transform) for text in cells[at]])
                else:
                    out.append(read(cells[at], f.hint, cls.FORMAT.null, f.source.transform))
            except ValueError as e:
                raise ValueError(f"{name}: {e}") from None
        return tuple(out)


def read(text, hint, null, transform=None):
    """A cell's or a part's value: null, or the text cleaned, then read by the type."""
    t, nullable = optional(hint)
    if text == null:
        if not nullable:
            raise ValueError(f"{text!r} is null, which this column can't be")
        return None
    if transform:
        text = transform(text)
    return t.read(text, null) if is_struct(t) else READERS[t](text)


def source_problem(format, f):
    """What is wrong with a field's declaration, if anything."""
    s = f.source
    if not isinstance(s, Source):
        return "needs a source: column(...) or columns(...)"
    t, nullable = optional(f.hint)
    if s.title is not None and not (isinstance(s.title, str) and s.title.strip()):
        return "a title is text, the header a datatable shows"
    if isinstance(s.first, str) and not format.header:
        return f"names column {s.first!r}, but the format has no header"
    if s.absent and not isinstance(s.first, str):
        return "absent is for a column named in the header"
    if s.absent is no_value and not nullable:
        return "absent=no_value is null, which only a type with '| None' allows"
    if s.many and not isinstance(s.last, int):
        return "columns(...) ends at a position: -1 for each line's last"
    if isinstance(s.first, int) and s.first < 1 or s.last == 0 or isinstance(s.last, int) and 0 < s.last < s.first:
        return "columns are counted from 1, first to last"
    if s.many:
        if get_origin(t) is not list or nullable:
            return "columns(...) reads a list: its type is list[...]"
        if isinstance(s.first, str):
            return "columns(...) takes positions, not names"
        t, nullable = optional(get_args(t)[0])
    elif get_origin(t) is list:
        return "a list is read from columns(...)"
    if t not in READERS and not is_struct(t):
        return f"no reader for {t}"


def declared_fields(cls):
    docs = attribute_docs(cls)
    return {name: Field(hint, vars(cls).get(name), docs.get(name)) for name, hint in get_type_hints(cls).items()}


def is_struct(t):
    return isinstance(t, type) and issubclass(t, Struct)


def optional(hint):
    """(the type without '| None', whether None is allowed)"""
    if isinstance(hint, types.UnionType) and type(None) in get_args(hint):
        others = [a for a in get_args(hint) if a is not type(None)]
        if len(others) == 1:
            return others[0], True
    return hint, False


def sql_type(hint):
    t, _ = optional(hint)
    if t in SCALARS:
        return SCALARS[t][1]
    if get_origin(t) is list:
        return sql_type(get_args(t)[0]) + "[]"
    if is_struct(t):
        return f"struct({', '.join(f'{name} {sql_type(f.hint)}' for name, f in t.FIELDS.items())})"
    raise TypeError(f"no SQL type for {t}")


def json_type(hint, doc=None, title=None):
    t, nullable = optional(hint)
    if t in SCALARS:
        out = {"type": SCALARS[t][0]}
    elif get_origin(t) is list:
        out = {"type": "array", "items": json_type(get_args(t)[0])}
    elif is_struct(t):
        out = object_type(t)
    else:
        raise TypeError(f"no JSON type for {t}")
    if nullable:
        out["type"] = [out["type"], "null"]
    return out | ({"title": title} if title else {}) | ({"description": doc} if doc else {})


def object_type(cls):
    """A Struct or an Adapter's row as a JSON object."""
    return {"type": "object",
            **({"description": inspect.cleandoc(cls.__doc__)} if cls.__doc__ else {}),
            "properties": {name: json_type(f.hint, f.doc, f.title) for name, f in cls.FIELDS.items()},
            "required": list(cls.FIELDS), "additionalProperties": False}


def attribute_docs(cls):
    """The string written under each annotated attribute in the class's source, by attribute name."""
    body = ast.parse(textwrap.dedent(inspect.getsource(cls))).body[0].body
    return {a.target.id: b.value.value for a, b in zip(body, body[1:])
            if isinstance(a, ast.AnnAssign) and isinstance(b, ast.Expr)
            and isinstance(b.value, ast.Constant) and isinstance(b.value.value, str)}
