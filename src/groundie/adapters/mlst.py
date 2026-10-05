from groundie.adapters.mash import SAMPLE

NAME = "mlst"
VERSION = "0.3"
OPTIONS = ()

# mlst's default output: no header; the file, the scheme, the ST, then 'locus(allele)' for each of the
# scheme's seven loci. The file is named by the path mlst was given, so the directories and a FASTA
# extension are stripped, as in the mash adapter. ST is text: mlst writes '-' when no type is assigned.
# The loci are named by position, locus_1 to locus_7, in mlst's order; each holds the allele as text with
# mlst's flags ('~4', '4?', '-'). A genome with no scheme found ('-' '-' and no loci) has null loci.
LOCI = [f"locus_{i}" for i in range(1, 8)]
COLUMNS = ", ".join(f"'{c}': 'VARCHAR'" for c in ["file", "scheme", "ST", *LOCI])
CANONICAL_SQL = f"""
select
    {SAMPLE.format(col="file")} as sample,
    scheme,
    "ST" as sequence_type,
    regexp_extract(columns(* exclude (file, scheme, "ST")), '\\((.*)\\)', 1)
from read_csv({{path}}, delim = '\\t', header = false, auto_detect = false, null_padding = true, columns = {{{{{COLUMNS}}}}})
"""


def canonical_sql(options):
    return CANONICAL_SQL
