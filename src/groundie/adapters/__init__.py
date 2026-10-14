from groundie.adapters import amrfinder, csv, json, mash, mlst, newick, tsv

ADAPTERS = {m.NAME: m for m in (amrfinder.Amrfinder, mlst.Mlst, mash, tsv, csv, json, newick)}
