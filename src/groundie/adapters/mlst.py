"""mlst's default output: no header; the file, the scheme, the ST, then 'locus(call)' for each of the
scheme's loci, as many as the scheme has. One row per genome, as SCHEMA describes it.

The file is named by the path mlst was given, so the directories and a FASTA extension are stripped, as in
the mash adapter. mlst writes '-' for no scheme found and for no ST assigned; both are read as null. The
loci are kept in mlst's order as calls: the locus name, and the call as mlst wrote it, flags included
('3', '~4', '4?', '-', '4,5'). A genome with no scheme found has no calls."""

import re
from pathlib import Path

NAME = "mlst"
VERSION = "preliminary"
OPTIONS = ()
COLUMNS = (("sample", "varchar"), ("scheme", "varchar"), ("sequence_type", "integer"),
           ("calls", "struct(locus varchar, call varchar)[]"))

SCHEMA = {
    "$schema": "https://json-schema.org/draft/2020-12/schema",
    "title": "mlst result",
    "description": "One genome's mlst result.",
    "type": "object",
    "properties": {
        "sample": {"type": "string", "description": "The genome's file, without directories or FASTA extension."},
        "scheme": {"type": ["string", "null"], "description": "The MLST scheme; null when none was found."},
        "sequence_type": {"type": ["integer", "null"], "minimum": 0,
                          "description": "The ST; null when none was assigned."},
        "calls": {
            "type": "array",
            "description": "One call per locus of the scheme, in mlst's order; empty when no scheme was found.",
            "items": {
                "type": "object",
                "properties": {
                    "locus": {"type": "string"},
                    "call": {"type": "string", "description": "As mlst wrote it, flags included: '3', '~4', '4?', '-'."},
                },
                "required": ["locus", "call"],
                "additionalProperties": False,
            },
        },
    },
    "required": ["sample", "scheme", "sequence_type", "calls"],
    "additionalProperties": False,
}

FASTA = re.compile(r"\.(fa|fasta|fna|fas)(\.gz)?$")
LOCUS = re.compile(r"(.+)\((.*)\)")  # 'tpi(~4)': the locus, then the call
NONE = "-"


def rows(path, options):
    path = Path(path)
    out = []
    for n, line in enumerate(path.read_text().splitlines(), 1):
        if line.strip():
            out.append(parse(line, f"{path.name}, line {n}"))
    return out


def parse(line, where):
    """One line of mlst's output as (sample, scheme, sequence_type, calls)."""
    fields = line.split("\t")
    if len(fields) < 3:
        raise ValueError(f"{where}: {line!r} is not the file, the scheme and the ST, tab-separated")
    file, scheme, st, *loci = fields
    if not (st == NONE or st.isdigit()):
        raise ValueError(f"{where}: ST {st!r} is not a number or '-'")
    calls = []
    for text in loci:
        m = LOCUS.fullmatch(text)
        if not m:
            raise ValueError(f"{where}: {text!r} is not 'locus(call)'")
        calls.append({"locus": m[1], "call": m[2]})
    sample = FASTA.sub("", file.rsplit("/", 1)[-1])
    return sample, None if scheme == NONE else scheme, None if st == NONE else int(st), calls
