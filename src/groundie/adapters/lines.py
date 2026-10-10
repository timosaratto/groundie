"""Adapters for line-based text, declared as a class. Each annotated attribute is one column, read from the
next field of the line, in order; the string under it is its description, and the class docstring the row's:

    class Mash(Lines, sep="\\t", null="-"):
        \"\"\"One pair of genomes.\"\"\"

        reference: str = sample_name
        "The first genome."

A field's text is converted by its type (str, int, float; '| None' allows the null token), or read by the
step given as its value: a function from the text to the value, such as sample_name. A Match type reads
itself. rest(...) on the last field takes the line's remaining fields as a list. The class becomes the adapter: NAME (the
module's name), VERSION, OPTIONS, COLUMNS, SCHEMA and rows(path, options), as cli.load_input expects."""

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


def sample_name(text):
    """A genome's file as its sample name: the directories and a FASTA extension stripped, as in mash.SAMPLE."""
    return FASTA.sub("", text.rsplit("/", 1)[-1])


class Match:
    """A struct read from text by a regex whose named groups are its annotated fields, each converted by its
    type. A field annotated with it needs no step:

        class Call(Match, pattern=r"(?P<locus>.+)\\((?P<call>.*)\\)"):
            locus: str
            call: str
    """

    def __init_subclass__(cls, pattern):
        super().__init_subclass__()
        cls.REGEX = re.compile(pattern)
        # shown in errors as the groups' names: '(?P<locus>.+)\((?P<call>.*)\)' is 'locus(call)'
        cls.SHOWN = re.sub(r"\(\?P<(\w+)>[^()]*\)", r"\1", pattern).replace("\\", "")
        cls.FIELDS = {k: optional(v)[0] for k, v in get_type_hints(cls).items()}
        if set(cls.REGEX.groupindex) != set(cls.FIELDS):
            raise TypeError(f"{cls.__name__}: the pattern's groups {sorted(cls.REGEX.groupindex)} are not its "
                            f"fields {sorted(cls.FIELDS)}")

    @classmethod
    def parse(cls, text):
        m = cls.REGEX.fullmatch(text)
        if not m:
            raise ValueError(f"{text!r} is not {cls.SHOWN!r}")
        return {k: None if v is None else convert(cls.FIELDS[k], v) for k, v in m.groupdict().items()}


def is_struct(t):
    return isinstance(t, type) and issubclass(t, Match)


class rest:  # noqa: N801 - reads as a function where it is used
    """The last field: the line's remaining fields as a list, each read by the step or converted by type."""

    def __init__(self, step=None):
        self.step = step


class Field(NamedTuple):
    name: str
    type: object  # the Python type, without '| None'
    nullable: bool
    step: object  # text -> value; None converts by type
    rest: bool
    doc: str | None


class Lines:
    VERSION = "preliminary"
    OPTIONS = ()

    def __init_subclass__(cls, sep="\t", null=None):
        super().__init_subclass__()
        cls.NAME = cls.__module__.rsplit(".", 1)[-1]
        cls.SEP, cls.NULL = sep, null
        cls.FIELDS = fields(cls)
        cls.COLUMNS = tuple((f.name, sql_type(f.type)) for f in cls.FIELDS)
        cls.SCHEMA = {
            "$schema": "https://json-schema.org/draft/2020-12/schema",
            "title": f"{cls.NAME} result",
            "description": inspect.cleandoc(cls.__doc__ or ""),
            "type": "object",
            "properties": {f.name: json_type(f.type, f.nullable, f.doc) for f in cls.FIELDS},
            "required": [f.name for f in cls.FIELDS],
            "additionalProperties": False,
        }
        tail = f" and any number of {cls.FIELDS[-1].name}" if cls.FIELDS[-1].rest else ""
        cls.EXPECTED = ", ".join(f.name for f in cls.FIELDS if not f.rest) + tail

    @classmethod
    def rows(cls, path, options):
        path = Path(path)
        out = []
        for n, line in enumerate(path.read_text().splitlines(), 1):
            if line.strip():
                try:
                    out.append(cls.parse(line.split(cls.SEP)))
                except ValueError as e:
                    raise ValueError(f"{path.name}, line {n}: {e}") from None
        return out

    @classmethod
    def parse(cls, texts):
        """One line's fields as a row, in COLUMNS order."""
        fixed = [f for f in cls.FIELDS if not f.rest]
        if len(texts) < len(fixed) or (len(texts) > len(fixed) and not cls.FIELDS[-1].rest):
            raise ValueError(f"{cls.SEP.join(texts)!r} has {len(texts)} fields, not {cls.EXPECTED}")
        row = [cls.read(f, text) for f, text in zip(fixed, texts)]
        if cls.FIELDS[-1].rest:
            f = cls.FIELDS[-1]
            row.append([value(f, get_args(f.type)[0], text) for text in texts[len(fixed):]])
        return tuple(row)

    @classmethod
    def read(cls, f, text):
        if cls.NULL is not None and text == cls.NULL:
            if not f.nullable:
                raise ValueError(f"{f.name} cannot be null ({cls.NULL!r})")
            return None
        return value(f, f.type, text)


def value(f, t, text):
    """The text read by the field's step, or converted to the type t."""
    try:
        return f.step(text) if f.step else convert(t, text)
    except ValueError as e:
        raise ValueError(f"{f.name} {e}") from None


def convert(t, text):
    if t is str:
        return text
    if t is int and INTEGER.fullmatch(text):
        return int(text)
    if t is float:
        try:
            return float(text)
        except ValueError:
            pass
    raise ValueError(f"{text!r} is not {'an integer' if t is int else 'a number'}")


def fields(cls):
    """The class's annotated attributes as Fields, checked: rest() last and on a list, and each read by a step,
    a scalar type, or a Match."""
    docs = attribute_docs(cls)
    out = []
    for name, hint in get_type_hints(cls).items():
        t, nullable = optional(hint)
        default = vars(cls).get(name)
        tail = isinstance(default, rest)
        step = default.step if tail else default
        where = f"{cls.__name__}.{name}"
        if tail and get_origin(t) is not list:
            raise TypeError(f"{where}: rest() needs a list type, not {t}")
        target = get_args(t)[0] if tail else t
        if step is None and is_struct(target):
            step = target.parse
        if step is None and target not in SCALARS:
            raise TypeError(f"{where}: {target} is read from text by a step; give one")
        out.append(Field(name, t, nullable, step, tail, docs.get(name)))
    if any(f.rest for f in out[:-1]):
        raise TypeError(f"{cls.__name__}: only the last field can be rest()")
    return tuple(out)


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
        return f"struct({', '.join(f'{k} {sql_type(v)}' for k, v in get_type_hints(t).items())})"
    raise TypeError(f"no SQL type for {t}")


def json_type(t, nullable=False, doc=None):
    if t in SCALARS:
        out = {"type": [SCALARS[t][0], "null"] if nullable else SCALARS[t][0]}
    elif get_origin(t) is list:
        out = {"type": ["array", "null"] if nullable else "array", "items": json_type(get_args(t)[0])}
    elif is_struct(t):
        hints, docs = get_type_hints(t), attribute_docs(t)
        out = {"type": ["object", "null"] if nullable else "object",
               "properties": {k: json_type(*optional(v), docs.get(k)) for k, v in hints.items()},
               "required": list(hints), "additionalProperties": False}
    else:
        raise TypeError(f"no JSON type for {t}")
    return out | ({"description": doc} if doc else {})


def attribute_docs(cls):
    """The string written under each annotated attribute in the class's source, by attribute name."""
    body = ast.parse(textwrap.dedent(inspect.getsource(cls))).body[0].body
    return {a.target.id: b.value.value for a, b in zip(body, body[1:])
            if isinstance(a, ast.AnnAssign) and isinstance(b, ast.Expr)
            and isinstance(b.value, ast.Constant) and isinstance(b.value.value, str)}
