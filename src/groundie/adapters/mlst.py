"""mlst's default output: tab-separated, no header, '-' for nothing."""

from groundie.adapters.lines import Lines, Match, rest, sample_name


class MlstCall(Match, pattern=r"(?P<locus>.+)\((?P<call>.*)\)"):
    locus: str
    call: str
    "As mlst wrote it, flags included: '3', '~4', '4?', '-'."


class Mlst(Lines, sep="\t", null="-"):
    """One genome's mlst result."""

    sample: str = sample_name
    "The genome's file, without directories or FASTA extension."

    scheme: str | None
    "The MLST scheme; null when none was found."

    sequence_type: int | None
    "The ST; null when none was assigned."

    calls: list[MlstCall] = rest()
    "One call per locus, in mlst's order; empty when no scheme was found."
