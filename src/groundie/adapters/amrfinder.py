NAME = "amrfinder"
VERSION = "0.1"
OPTIONS = ()

CANONICAL_SQL = """
select
    "Name"                          as sample,
    "Contig id"                     as contig_id,
    "Start"                         as start_pos,
    "Stop"                          as stop_pos,
    "Strand"                        as strand,
    "Element symbol"                as element_symbol,
    "Element name"                  as element_name,
    "Scope"                         as scope,
    "Type"                          as element_type,
    "Subtype"                       as element_subtype,
    "Class"                         as element_class,
    "Subclass"                      as element_subclass,
    "Method"                        as method,
    "Target length"                 as target_length,
    "Reference sequence length"     as ref_length,
    "% Coverage of reference"       as pct_coverage,
    "% Identity to reference"       as pct_identity,
    "Alignment length"              as alignment_length,
    "Closest reference accession"   as ref_accession,
    "Closest reference name"        as ref_name,
    "HMM accession"::varchar        as hmm_accession,
    "HMM description"::varchar      as hmm_description
from read_csv({path}, delim = '\\t', header = true, nullstr = 'NA', union_by_name = true)
"""


def canonical_sql(options):
    return CANONICAL_SQL
