NAME = "mash"
VERSION = "preliminary"
OPTIONS = ()

# `mash dist` or `mash triangle -E` output: no header; reference, query, distance, p-value, shared hashes
# as 'shared/total'. mash names genomes by the path it was given, so the directories and a FASTA or sketch
# extension are stripped ('genomes/A.fa.gz' -> 'A'), to match the sample names of the other adapters.
SAMPLE = "regexp_replace(regexp_replace({col}, '^.*/', ''), '\\.(fa|fasta|fna|fas|msh)(\\.gz)?$', '')"
CANONICAL_SQL = f"""
select
    {SAMPLE.format(col="reference")} as reference,
    {SAMPLE.format(col="query")} as query,
    distance,
    p_value,
    split_part(shared, '/', 1)::integer as shared_hashes,
    split_part(shared, '/', 2)::integer as sketch_hashes
from read_csv({{path}}, delim = '\\t', header = false, columns = {{{{
    'reference': 'VARCHAR', 'query': 'VARCHAR', 'distance': 'DOUBLE', 'p_value': 'DOUBLE', 'shared': 'VARCHAR'}}}})
"""


def canonical_sql(options):
    return CANONICAL_SQL
