"""Any JSON file, newline-delimited or a top-level array of objects. Keys become columns as they are;
one column is added first: 'file', the basename of the file each row came from."""

NAME = "json"
VERSION = "0.1"
OPTIONS = ()


def canonical_sql(options):
    return """
select parse_filename(file) as file, * exclude (file)
from read_json({path}, format = 'auto', filename = 'file', union_by_name = true)
"""
