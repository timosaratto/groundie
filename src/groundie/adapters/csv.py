"""Any comma-separated file, read as tsv is."""

from groundie.adapters import tsv

NAME = "csv"
VERSION = "0.1"
OPTIONS = tsv.OPTIONS


def canonical_sql(options):
    return tsv.delimited_sql(",", options)
