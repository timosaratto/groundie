"""Records: the rows an adapter emits, declared as a class, and read from text in two explicit steps.

A parse function binds: it splits its text into named pieces and hands them to Record.read. Each field
interprets: its value is the reader read applies to its piece, and the string under it is its description.

    class Pair(Record):
        \"\"\"One pair of genomes.\"\"\"

        reference: str = sample_name
        "The first genome."

        distance: float | None = none_if("NA", number)
        "The distance; null where the tool wrote NA."

    def parse(line):
        reference, distance = line.split("\\t")
        return Pair.read(reference=reference, distance=distance)

    ADAPTER = Lines(Pair, parse)

A field without a reader keeps its text, so only a str field (or list[str]) can go without one. read checks
each value against its field's type. The class gives the SQL COLUMNS and the JSON SCHEMA; Lines makes the
adapter cli.load_input expects."""

import ast
import inspect
import re
import textwrap
import types
from pathlib import Path
from typing import NamedTuple, get_args, get_origin, get_type_hints

FASTA = re.compile(r"\.(fa|fasta|fna|fas)(\.gz)?$")
INTEGER = re.compile(r"[+-]?\d+")
SCALARS = {str: ("string", "varchar"), int: ("integer", "bigint"), float: ("number", "double")}


# Readers: text in, value out; a ValueError says what the text is not.

def sample_name(text):
    """A genome's file as its sample name: the directories and a FASTA extension stripped, as in mash.SAMPLE."""
    return FASTA.sub("", text.rsplit("/", 1)[-1])


def integer(text):
    if not INTEGER.fullmatch(text):
        raise ValueError(f"{text!r} is not an integer")
    return int(text)


def number(text):
    try:
        return float(text)
    except ValueError:
        raise ValueError(f"{text!r} is not a number") from None


def none_if(token, then=None):
    """A reader: None where the text is the token, else the text, or then(text)."""
    return lambda text: None if text == token else then(text) if then else text


def each(reader):
    """A reader of a list of texts: the reader applied to each."""
    return lambda texts: [reader(text) for text in texts]


def fullmatch(pattern, text, shown=None):
    """The regex's groups, matched against the whole text; shown names the expected form in the error."""
    m = re.fullmatch(pattern, text)
    if not m:
        raise ValueError(f"{text!r} is not {shown or pattern!r}")
    return m.groups()


class Field(NamedTuple):
    hint: object
    reader: object  # text -> value; None keeps the text
    doc: str | None


class Record:
    def __init_subclass__(cls):
        super().__init_subclass__()
        docs = attribute_docs(cls)
        cls.FIELDS = {name: Field(hint, vars(cls).get(name), docs.get(name))
                      for name, hint in get_type_hints(cls).items()}
        for name, f in cls.FIELDS.items():
            if f.reader is None and not takes_text(f.hint):
                raise TypeError(f"{cls.__name__}.{name}: a {f.hint} field needs a reader; only text is kept as is")

    @classmethod
    def read(cls, **texts):
        """The record as a dict: each field's reader applied to its text, the value checked against its type."""
        if texts.keys() != cls.FIELDS.keys():
            raise TypeError(f"{cls.__name__}.read: needs {', '.join(cls.FIELDS)}; got {', '.join(texts)}")
        out = {}
        for name, f in cls.FIELDS.items():
            try:
                out[name] = f.reader(texts[name]) if f.reader else texts[name]
            except ValueError as e:
                raise ValueError(f"{name}: {e}") from None
            if not fits(out[name], f.hint):
                raise TypeError(f"{cls.__name__}.{name}: {out[name]!r} is not {f.hint}")
        return out


class Lines:
    """An adapter for text with one record per line: parse is given each non-blank line."""

    VERSION = "preliminary"
    OPTIONS = ()

    def __init__(self, record, parse):
        self.record, self.parse = record, parse
        self.NAME = record.__module__.rsplit(".", 1)[-1]
        self.COLUMNS = tuple((name, sql_type(f.hint)) for name, f in record.FIELDS.items())
        self.SCHEMA = {"$schema": "https://json-schema.org/draft/2020-12/schema", "title": f"{self.NAME} result",
                       **json_type(record)}

    def rows(self, path, options):
        path = Path(path)
        out = []
        for n, line in enumerate(path.read_text().splitlines(), 1):
            if line.strip():
                try:
                    out.append(tuple(self.parse(line).values()))
                except ValueError as e:
                    raise ValueError(f"{path.name}, line {n}: {e}") from None
        return out


def is_record(t):
    return isinstance(t, type) and issubclass(t, Record)


def optional(hint):
    """(the type without '| None', whether None is allowed)"""
    if isinstance(hint, types.UnionType) and type(None) in get_args(hint):
        others = [a for a in get_args(hint) if a is not type(None)]
        if len(others) == 1:
            return others[0], True
    return hint, False


def takes_text(hint):
    t, _ = optional(hint)
    return t is str or (get_origin(t) is list and get_args(t)[0] is str)


def fits(value, hint):
    t, nullable = optional(hint)
    if value is None:
        return nullable
    if get_origin(t) is list:
        return isinstance(value, list) and all(fits(v, get_args(t)[0]) for v in value)
    if is_record(t):
        return isinstance(value, dict) and value.keys() == t.FIELDS.keys()
    if t in (int, float) and isinstance(value, bool):
        return False
    return isinstance(value, (int, float) if t is float else t)


def sql_type(hint):
    t, _ = optional(hint)
    if t in SCALARS:
        return SCALARS[t][1]
    if get_origin(t) is list:
        return sql_type(get_args(t)[0]) + "[]"
    if is_record(t):
        return f"struct({', '.join(f'{name} {sql_type(f.hint)}' for name, f in t.FIELDS.items())})"
    raise TypeError(f"no SQL type for {t}")


def json_type(hint, doc=None):
    t, nullable = optional(hint)
    if t in SCALARS:
        out = {"type": SCALARS[t][0]}
    elif get_origin(t) is list:
        out = {"type": "array", "items": json_type(get_args(t)[0])}
    elif is_record(t):
        out = {"type": "object",
               **({"description": inspect.cleandoc(t.__doc__)} if t.__doc__ else {}),
               "properties": {name: json_type(f.hint, f.doc) for name, f in t.FIELDS.items()},
               "required": list(t.FIELDS), "additionalProperties": False}
    else:
        raise TypeError(f"no JSON type for {t}")
    if nullable:
        out["type"] = [out["type"], "null"]
    return out | ({"description": doc} if doc else {})


def attribute_docs(cls):
    """The string written under each annotated attribute in the class's source, by attribute name."""
    body = ast.parse(textwrap.dedent(inspect.getsource(cls))).body[0].body
    return {a.target.id: b.value.value for a, b in zip(body, body[1:])
            if isinstance(a, ast.AnnAssign) and isinstance(b, ast.Expr)
            and isinstance(b.value, ast.Constant) and isinstance(b.value.value, str)}
