"""mlst's default output: tab-separated, no header, '-' for nothing."""

from groundie.adapters.utils import Lines, Record, each, fullmatch, integer, none_if, sample_name


class MlstCall(Record):
    locus: str
    call: str
    "As mlst wrote it, flags included: '3', '~4', '4?', '-'."


def parse_call(text):
    locus, call = fullmatch(r"(.+)\((.*)\)", text, shown="locus(call)")
    return MlstCall.read(locus=locus, call=call)


class Mlst(Record):
    """One genome's mlst result."""

    sample: str = sample_name
    "The genome's file, without directories or FASTA extension."

    scheme: str | None = none_if("-")
    "The MLST scheme; null when none was found."

    sequence_type: int | None = none_if("-", integer)
    "The ST; null when none was assigned."

    calls: list[MlstCall] = each(parse_call)
    "One call per locus, in mlst's order; empty when no scheme was found."


def parse(line):
    file, scheme, st, *loci = line.split("\t")
    return Mlst.read(sample=file, scheme=scheme, sequence_type=st, calls=loci)


ADAPTER = Lines(Mlst, parse)
