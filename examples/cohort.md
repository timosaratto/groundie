# Staphylococcus aureus cohort

```input mlst
adapter: mlst
path: ../data/mlst_results/SAMD*.tsv
```

```input amrfinder
adapter: amrfinder
path: ../data/amrfinder_results/*.tsv
```

```input tree
adapter: newick
path: ../data/tree_results/cohort.nwk
```

## Sequence types

One row per genome, from mlst, one table per scheme, headed by the scheme's loci. A – means no ST was assigned, or, in a call, no allele was found. In a call, `~` marks a novel allele and `?` a partial match.

```MlstTable st
source: mlst
sort: [sample]
```

## Cohort summary

One row per genome: its sequence type and how many AMR, virulence and stress hits AMRFinderPlus reported.

```sql summary
select
    m.sample,
    m.sequence_type,
    count(*) filter (where a.element_type = 'AMR')       as amr_hits,
    count(*) filter (where a.element_type = 'VIRULENCE') as virulence_hits,
    count(*) filter (where a.element_type = 'STRESS')    as stress_hits
from mlst m
left join amrfinder a on a.sample = m.sample
group by m.sample, m.sequence_type
order by m.sample
```

```DataTable summary_table
source: summary
```

## Phylogeny

A neighbour-joining tree of mash distances between the ten genomes: coarse, good for seeing which genomes sit together and not for dating or fine branching order. The genome with no sequence type sits inside the ST8 clade, which suggests it is an ST8 whose typing failed rather than a new lineage.

```Tree phylogeny
source: tree
tips: summary
match: sample
label: [sample, sequence_type]
color: sequence_type
```

## Resistome

AMRFinderPlus hits of type AMR, genes and point mutations: one row per genome, one column per gene. The shade is the best identity to the reference, from 90 % (light) to 100 % (dark); ◆ marks a gene with at least one partial copy (coverage below 100 %). Hover a cell for its values and how many copies were found.

```Matrix resistome
source: amrfinder
filters:
  element_type: AMR
group_by: [sample, element_class, element_symbol]
measures:
  best_identity: max(pct_identity)
  lowest_coverage: min(pct_coverage)
  copies: count(*)
sort: [sample, element_class, element_symbol]
rows: [sample]
cols: element_symbol
value: best_identity
scale: [90, 100]
mark: {column: lowest_coverage, below: 100}
```
