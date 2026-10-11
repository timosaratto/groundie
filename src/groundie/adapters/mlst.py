"""mlst's default output: the file, the scheme, the ST, then 'locus(call)' for each of the scheme's loci."""

from groundie.adapters.utils import Adapter, Struct, Tsv, column, columns, extract_sample_name


class MlstCall(Struct):
    PATTERN = "{locus}({call})"

    locus: str
    call: str | None
    "As mlst wrote it, flags included: '3', '~4', '4?'; null where mlst wrote '-', no allele found."


class Mlst(Adapter):
    """One genome's mlst result."""

    FORMAT = Tsv(header=False, null="-")

    sample: str = column(1, transform=extract_sample_name)
    "The genome's file, without directories or FASTA extension."

    scheme: str | None = column(2)
    "The MLST scheme; null when none was found."

    sequence_type: int | None = column(3)
    "The ST; null when none was assigned."

    calls: list[MlstCall] = columns(4, -1)
    "One call per locus, in mlst's order; empty when no scheme was found."
