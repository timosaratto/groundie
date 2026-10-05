"""Any tab-separated file. Columns are kept as in the file; one column is added first: 'file', the basename
of the file each row came from."""

NAME = "tsv"
VERSION = "0.1"
OPTIONS = ("null", "header")  # null: the string read as null (DuckDB's default if absent); header: default true


def canonical_sql(options):
    return delimited_sql("\\t", options)


def delimited_sql(delim, options):
    """Shared with csv, which differs only in the delimiter."""
    args = f"delim = '{delim}', header = {str(options.get('header', True)).lower()}"
    if "null" in options:
        # Doubled braces survive the {path} format.
        args += ", nullstr = '" + str(options["null"]).replace("'", "''").replace("{", "{{").replace("}", "}}") + "'"
    return f"""
select parse_filename(file) as file, * exclude (file)
from read_csv({{path}}, {args}, filename = 'file', union_by_name = true)
"""
