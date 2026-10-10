from groundie.adapters.mash import SAMPLE

NAME = "mlst"
VERSION = "preliminary"
OPTIONS = ()

# mlst's default output: no header; the file, the scheme, the ST, then 'locus(allele)' for each of the
# scheme's loci, as many as the scheme has. The file is named by the path mlst was given, so the
# directories and a FASTA extension are stripped, as in the mash adapter. ST is text: mlst writes '-' when
# no type is assigned. The loci are named by position, locus_1 to locus_<n>, in mlst's order, since the
# names depend on the scheme; n is the most loci of any genome in the input, and a genome with fewer has
# null loci after its own, as does a genome with no scheme found ('-' '-' and no loci). Each holds mlst's
# text as is, the locus name with the allele and its flags ('arcC(3)', 'tpi(~4)', 'pta(4?)', 'gmk(-)').
# Each line is split on tabs here, not read as CSV, because the number of fields varies by scheme.
FIELDS = """
select string_split(rtrim(line, chr(13)), chr(9)) as f
from (select unnest(string_split(content, chr(10))) as line from read_text({path}))
where trim(line) <> ''
"""


def shape(con, paths):
    """What canonical_sql takes from the files themselves: how many loci. paths is a DuckDB list literal."""
    sql = f"select coalesce(max(len(f)), 3) - 3 from ({FIELDS.format(path=paths)})"
    return {"loci": con.execute(sql).fetchone()[0]}


def canonical_sql(options, loci):
    columns = "".join(f",\n    f[{i + 3}] as locus_{i}" for i in range(1, loci + 1))
    return f"""
select
    {SAMPLE.format(col="f[1]")} as sample,
    f[2] as scheme,
    f[3] as sequence_type{columns}
from ({FIELDS})
"""
