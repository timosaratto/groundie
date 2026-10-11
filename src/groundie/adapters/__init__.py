from groundie.adapters import amrfinder, csv, json, mash, mlst, newick, tsv

ADAPTERS = {m.NAME: m for m in (amrfinder, mlst.Mlst, mash, tsv, csv, json, newick)}
