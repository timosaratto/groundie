"""AMRFinderPlus output: a header, then one hit per line, read by header name. A file without --name's Name
column takes its sample from the file's name. Each column's title is AMRFinderPlus's own header.

Only version 4's format is read. Versions before 4.0 name most columns differently ('Gene symbol' for
'Element symbol', 'Protein identifier' for 'Protein id', 'HMM id' for 'HMM accession', ...), so their files
fail with the column that's missing."""

from groundie.adapters.utils import Adapter, Tsv, column, no_value, sample_from_file_name


class Amrfinder(Adapter):
    """One AMRFinderPlus hit: a gene or point mutation found in one genome."""

    FORMAT = Tsv(header=True, null="NA")
    # Every field is '| None': NA is null in any column, as read_csv(nullstr = 'NA') read it before.

    sample: str | None = column("Name", title="Sample", absent=sample_from_file_name)
    "The genome, as given to --name; without --name, the file's name without its directories and extension."

    protein_id: str | None = column("Protein id")
    "The protein's id in the input; null for nucleotide input."

    contig_id: str | None = column("Contig id")
    "The contig the hit is on; null for protein input."

    start_pos: int | None = column("Start")
    "Null for protein input."

    stop_pos: int | None = column("Stop")
    "Null for protein input."

    strand: str | None = column("Strand")
    "'+' or '-'; null for protein input."

    element_symbol: str | None = column("Element symbol")
    "The gene or mutation: 'tet(38)', 'parC_S80F'."

    element_name: str | None = column("Element name")
    "What it is: 'tetracycline efflux MFS transporter Tet(38)'."

    scope: str | None = column("Scope")
    "'core' or 'plus'."

    element_type: str | None = column("Type")
    "AMR, STRESS or VIRULENCE."

    element_subtype: str | None = column("Subtype")
    "AMR, POINT, BIOCIDE, METAL, VIRULENCE, ..."

    element_class: str | None = column("Class")
    "The drug class, such as QUINOLONE; null for virulence hits."

    element_subclass: str | None = column("Subclass")
    "Null for virulence hits."

    method: str | None = column("Method")
    "How it was found: EXACTX, BLASTX, POINTX, PARTIALX, HMM, ..."

    target_length: int | None = column("Target length")

    ref_length: int | None = column("Reference sequence length")
    "Null where no reference was aligned, as for a hit found by HMM alone; so are the fields to ref_name."

    pct_coverage: float | None = column("% Coverage of reference")

    pct_identity: float | None = column("% Identity to reference")

    alignment_length: int | None = column("Alignment length")

    ref_accession: str | None = column("Closest reference accession")

    ref_name: str | None = column("Closest reference name")

    hmm_accession: str | None = column("HMM accession")
    "Null unless an HMM matched."

    hmm_description: str | None = column("HMM description")
    "Null unless an HMM matched."

    hierarchy_node: str | None = column("Hierarchy node", absent=no_value)
    "The node of NCBI's Reference Gene Hierarchy the hit was placed at; null without --print_node."
