# Writing groundie blocks

A block is a fenced snippet in Markdown that loads data files or asks a question of them. `groundie` runs the
blocks in order and shows the result of each: a datatable, a grid or a phylogeny. Every value shown comes from
the data files, never from the text you write around the blocks.

There are two ways to run blocks:

- A report: a Markdown file of prose and blocks, built into one self-contained HTML page with
  `uv run groundie build`. See [Building a report](#building-a-report).
- In chat: the MCP tool `render_blocks` runs blocks and shows them to the user in the conversation. See
  [Rendering blocks in chat](#rendering-blocks-in-chat).

The blocks are the same in both.

## Blocks

A block opens with ```` ```<kind> <name> ````: one space, then a name of letters, digits and underscores. It
closes with a line of three backticks. Everything outside blocks is Markdown prose. There is no front matter.

````markdown
```datatable st
source: mlst
columns: [sample, sequence_type]
```
````

| Kind | Shows | Can be a source for later blocks |
| --- | --- | --- |
| `input` | Nothing | Yes |
| `sql` | Nothing | Yes |
| `query` | Nothing | Yes |
| `clusters` | Nothing | Yes: one row per sample, its cluster |
| `checks` | Nothing | Yes: one row per source row, PASS or FAIL |
| `calls` | Nothing | Yes: one row per sample, true or false per call |
| `datatable` | A datatable | Yes: the rows it shows |
| `matrix` | A grid | Yes: its query's rows, not the pivoted grid |
| `amrfinder_matrix` | A gene × sample grid of AMRFinderPlus hits | Yes: as a `matrix` |
| `mlst_table` | One datatable of mlst results per scheme, headed by its loci | Yes: its query's rows, not the split tables |
| `tree` | A phylogeny | No |

The source of a block is the block it reads, named by its `source:` key. Every block except a `tree` makes an
SQL table with the block's name, holding the rows described in the last column above. Later blocks read it
with `source: <name>`, or by that name inside a `sql` block.

The rules:

1. Every block has a kind and a name. A name is used once in the file, across blocks of all kinds.
2. `source:` names a block above it, or any `input`, wherever it is. So do a tree's `tips` and a calls
   block's `samples`.
3. A `query` is a datatable that shows nothing. Use it for rows that later blocks read but the reader need not
   see.

Other fences are shown as code. Give them at most one word after the backticks. A fence that starts with a
block kind but is not read as a block, such as ```` ```datatable ```` with no name, is shown as code with a
warning.

### Order

All `input` blocks load first, before any other block runs. After that no block can read files or reach the
network. So a block can read an input from anywhere in the file, but put the inputs first, where a reader
looks for them. The other blocks run top to bottom.

## Input blocks

An `input` block loads data files into an SQL table named after the block. Its body is YAML:

```yaml
adapter: amrfinder      # how to read the files: see Adapters
paths:                  # several files are stacked into one SQL table
  - data/amrfinderplus/KP01.tsv
  - data/amrfinderplus/KP02.tsv
```

```yaml
adapter: amrfinder
path: data/amrfinderplus/*.tsv  # a glob: every matching file, in name order
```

```yaml
adapter: tsv
path: data/metadata.tsv # one file
options: {null: NA}     # optional; only tsv and csv take options
```

- `adapter` is required. Give exactly one of `path` or `paths`.
- A path can be a glob (`*`, `?`, `[...]`, and `**` for any depth of directories). It must match at least
  one file, and directories are skipped. A file matched by several paths is loaded once. The provenance
  lists each matched file with its hash; 'View query' shows the paths as written.
- A path that is not a glob must name a file, not a directory.
- A relative path resolves against the directory of the report file. For blocks sent inline to
  `render_blocks`, it resolves against `base_dir`. An absolute path always works.
- An input shows nothing. To show its rows, read it with a `datatable`.

## Adapters

| Adapter | Files | SQL table |
| --- | --- | --- |
| `amrfinder` | AMRFinderPlus TSV, run with `--name` | The columns listed below. `NA` is read as null |
| `mlst` | mlst's default output: no header; the file, the scheme, the ST, then `locus(allele)` for each of seven loci | `sample` (the file, its directories and a `.fa`, `.fasta`, `.fna` or `.fas` extension, with or without `.gz`, stripped), `scheme`, `sequence_type`, then `locus_1` to `locus_7`, in mlst's order. All are text. `sequence_type` is `-` when none is assigned. A locus holds mlst's text as is: the locus name with the allele and mlst's flags, as in `arcC(3)`, `tpi(~4)`, `pta(4?)` and `gmk(-)`. A genome with no scheme found has null loci. A scheme with more than seven loci is an error |
| `mash` | `mash dist` or `mash triangle -E` output, no header | `reference`, `query` (sample names: the directories and a `.fa`, `.fasta`, `.fna`, `.fas` or `.msh` extension, with or without `.gz`, are stripped), `distance`, `p_value`, `shared_hashes`, `sketch_hashes` |
| `tsv`, `csv` | Any delimited file | `file` (the file's name), then the file's columns as named. Numbers are read as numbers |
| `json` | Newline-delimited objects, or one array of objects | `file`, then one column per key |
| `newick` | Newick trees | One row per node: `file`, `tree`, `node`, `parent`, `label`, `branch_length`, `support`, `is_tip`, `comment`. Drawn by a `tree` block |

The column names come from the adapter, not from the file. The `amrfinder` columns:

```
text:    sample, contig_id, strand, element_symbol, element_name, scope, element_type, element_subtype,
         element_class, element_subclass, method, ref_accession, ref_name, hmm_accession, hmm_description
numbers: start_pos, stop_pos, target_length, ref_length, pct_coverage, pct_identity, alignment_length
```

`sample` is the file's `Name` column. `element_type`, `element_subtype`, `element_class` and
`element_subclass` are its `Type`, `Subtype`, `Class` and `Subclass`.

Only `tsv` and `csv` take `options`:

- `null`: the string read as null
- `header`: default true. With `header: false` the columns are named `column0`, `column1`, and so on

## Query keys: query, datatable, matrix

`query`, `datatable` and `matrix` blocks share these keys. `source` is required; the rest are optional.

Plain rows:

```yaml
source: amrfinder                  # the block read
filters:                           # equality only; columns are ANDed, listed values are ORed
  element_class: [AMINOGLYCOSIDE, QUINOLONE]
columns: [sample, element_symbol]  # default: every column
sort: [sample]                     # ascending only
```

Grouped rows:

```yaml
source: amrfinder
group_by: [sample, element_class]  # one row per combination; used in place of columns
measures:                          # count(*) or fn(column): see below
  hits: count(*)
  best: max(pct_identity)
  genes: list(element_symbol)
sort: [sample]                     # a group_by column or a measure
```

Measures: `count`, `min`, `max`, `sum` and `avg`, and two that give text: `list(column)`, every value of the
group, sorted and joined with `, `, and `distinct(column)`, the same with each value once. Both skip nulls.

- These are all the query keys. Any other key is an error, except a matrix's layout keys.
- Filter values are compared as text: `sequence_type: 258` matches the text `258`.
- Joins, expressions, ranges and descending order need a `sql` block.
- A block reading a `sql` block, with no `sort`, keeps the order of the statement's `order by`.

## sql

One SELECT statement in DuckDB SQL. It reads the SQL tables of inputs and of blocks above it, by name, and has
no access to files or the network. Any DuckDB SELECT works, including `with` and `with recursive`, window
functions and `filter`. It shows nothing, so follow it with a `datatable` that reads it:

````markdown
```sql summary
select m.sample, m.sequence_type, count(a.sample) as amr_hits
from mlst m left join amr a on a.sample = m.sample
group by all
order by amr_hits desc, m.sample
```

```datatable summary_table
source: summary
```
````

## clusters

Single-linkage clusters of samples from a table of pairwise distances, such as a `mash` input, SNP distances
or cgMLST allele differences in long form.

```yaml
source: mash                    # one row per pair of samples
threshold: 0.001                # pairs at this distance or closer are linked (inclusive)
pairs: [reference, query]       # optional: the two sample columns; default the mash adapter's
distance: distance              # optional: the distance column; default 'distance'
linkage: single                 # optional; single is the only linkage so far
```

- It makes an SQL table with `sample`, `cluster` and `cluster_size`, one row per sample that appears in a
  pair. A sample linked to none is a cluster of its own.
- Clusters are named `C1`, `C2`, ... by size, largest first, then by their first sample in name order.
- It shows nothing. To show it, read it with a `datatable`, or join it in a `sql` block.

## checks

Pass or fail each row of a block against thresholds on its columns.

```yaml
source: amr                          # the rows to check
key: element_symbol                  # optional: the column identifying a row; default 'sample'
checks:                              # one or more columns, each with min, max, above or below
  pct_coverage: {min: 90, max: 100}  # min and max are inclusive
  pct_identity: {above: 95}          # above and below are exclusive
```

- It makes an SQL table with the key, each checked column, `qc` (`PASS` or `FAIL`) and `qc_failed`: the
  failed checks with their values, in the order written, such as `pct_coverage 85.2; pct_identity 93.1`. It is
  null on a PASS.
- A null value fails its check: `pct_identity missing`.
- Give min or above, not both, and max or below, not both.
- It shows nothing. To show it, read it with a `datatable`, or join it in a `sql` block.

## calls

Call a marker per sample from one column of its hits: true when the sample has any, or all, of some values.

```yaml
source: amrfinder               # one row per hit
samples: mlst                   # optional: a block whose key column lists every sample; default the samples in source
column: element_symbol          # optional: the column holding the values; default 'element_symbol'
key: sample                     # optional: the sample column, in source and samples; default 'sample'
calls:                          # one or more: name: {any: [values]} or {all: [values]}
  carbapenemase: {any: [blaKPC-2, blaKPC-3, blaNDM-1, blaOXA-48]}
  aerobactin: {all: [iucA, iucB, iucC, iutA]}
```

- It makes an SQL table with the key, then one true/false column per call, one row per sample, in name order.
- A sample with no rows in source is false for every call. AMRFinderPlus writes no rows for a genome with no
  hits, so give `samples` to keep such genomes.
- Values are compared as text, as in filters. A call's name may not be the key's.
- It shows nothing. To show it, read it with a `datatable`, or join it in a `sql` block.

## matrix

The query keys, plus the layout keys below. The query need not be grouped: when several of its rows fall in
one cell, the cell shows the first, or the one `pick` chooses, and its hover text says how many rows there
were. To summarise them some other way, group the query by the `rows` and `cols` columns and use `measures`.

Without `sort`, rows are ordered by how many cells they fill, most first, then by name, and columns by name:
a gene × sample matrix lists the most common genes first. With `sort`, each row and each column appears where
its first cell appears in the sorted result, so `sort` sets both orders together.

```yaml
rows: [sample]                               # one or more columns
cols: element_symbol                         # one column
value: best_identity                         # optional: the column a cell shows; without, a cell shows ●
scale: [95, 100]                             # optional: show value as a shade, low end light, high end dark;
                                             # values outside are clamped. Without, value is shown as text
digits: 2                                    # optional, for value as text: decimal places; default: as is
mark: {column: lowest_coverage, below: 80}   # optional: ◆ in cells where the column is below the number
pick: max(pct_coverage)                      # optional: of several rows in a cell, the one with the highest
                                             # (max) or lowest (min) value of the column; default: the first
symmetric: true                              # optional, with one rows column: a row for (a, b) also fills
                                             # (b, a); rows and columns list the same keys, by name
```

`rows`, `cols`, `value`, the `mark` column and the `pick` column are columns of the query's result. A cell's
values are in its `title` attribute in the HTML, shown on hover. A block that reads a matrix gets the query's
rows, not the grid.

A distance matrix from pairs, one number per cell, the full square from a file that lists each pair once:

```yaml
source: mash
rows: [reference]
cols: query
value: distance
digits: 4
symmetric: true
```

Without `symmetric`, each pair is one cell, so such a file gives one triangle. The diagonal stays blank unless
the query has rows that pair a key with itself.

## amrfinder_matrix

A `matrix` with its layout set for AMRFinderPlus hits: one row per gene, one column per sample.

```yaml
source: amrfinder                         # an input loaded by the amrfinder adapter
filters: {element_class: AMINOGLYCOSIDE}  # optional: any filters, e.g. element_type, element_class, scope, sample
partial_below: 80                         # optional: ◆ marks coverage below this (%); default 100
```

- It takes `source`, `filters`, `sort` and `partial_below`. For any other layout, write a `matrix`.
- Rows are `element_symbol`, with `element_type` before it unless a filter picks one type. Columns are
  `sample`. The order is a matrix's default: the most common genes first.
- A cell with several hits shows the most complete one (`pick: max(pct_coverage)`); hovering gives its
  identity and coverage. ◆ marks a cell whose most complete hit is still partial (coverage below
  `partial_below`).
- 'View query' shows the `matrix` it stands for, every setting written out.

## mlst_table

mlst results as one datatable per scheme, captioned with the scheme. The scheme's locus names are the
headers and the cells hold the alleles: the input's `arcC(3)` shows as `3` under `arcC`.

```yaml
source: mlst                  # an input loaded by the mlst adapter
filters: {scheme: saureus}    # optional: any filters on the input's columns
sort: [sample]                # optional
```

- It takes `source`, `filters` and `sort`. For other columns, write a `datatable`.
- The columns are `Sample`, `ST`, then the scheme's loci. Each header keeps its column's name in its
  `title`, such as `locus_1` for `arcC`.
- Allele numbers compare only within a scheme, so each scheme gets its own table. Genomes with no scheme
  found get the last table, `No scheme found`, with no loci.
- Without `sort`, the tables go by how many genomes they hold, most first, then by scheme, and the rows by
  sample. With `sort`, both follow the sorted result.
- Filters and sort read the input's columns as they are: a locus is the text `tpi(~4)`, not `~4`.
- Cells of one locus naming different loci within a scheme are an error.
- A block that reads it gets its query's rows, as the input has them, not the split tables.

## tree

```yaml
source: tree                    # an input loaded by the newick adapter
tips: mlst                      # optional: a block with one row per tip; any kind but tree, an input included
match: sample                   # needed with tips: the column that equals the tip label, as exact text
label: [sample, sequence_type]  # columns shown at the tips; default: the tip label
color: sequence_type            # the 3 most frequent values get a colour each, the rest share a grey
```

- These are all the tree keys. A tree takes no query keys: to filter or join the tip rows, do it in a `query`
  or `sql` block and name that block in `tips`.
- `label` and `color` need `tips`.
- When the input holds several trees, pick one with `file: core.nwk` and `tree: 1`, the tree's position in
  that file.
- A tree makes no SQL table. No block can read it.

A tip with no matching row is drawn and reported, not an error. Check the summary line:
`tree phylogeny: 12 tips, 12 matched`.

## What the output adds

The output adds these itself, so leave them out of the prose:

- under each datatable (an `mlst_table` included), matrix (an `amrfinder_matrix` included) and tree: its filters, its query, and its
  lineage, a walk back through its sources to the inputs, each with its settings or, for a `sql` block, its
  statement
- in a built report: a provenance footer listing every data file with its hash

The page gets no heading from the tool. Write your own, starting with a `# Heading`.

## A complete report

`cohort.md`, next to a `data` directory:

````markdown
# Klebsiella pneumoniae cohort

```input mlst
adapter: mlst
path: data/mlst/mlst.tsv
```

```input amrfinder
adapter: amrfinder
path: data/amrfinderplus/*.tsv
```

```input tree
adapter: newick
path: data/phylogeny/core.nwk
```

## Sequence types

One row per genome, from mlst. `-` means no sequence type was assigned.

```datatable st
source: mlst
columns: [sample, sequence_type]
sort: [sample]
```

```tree phylogeny
source: tree
tips: mlst
match: sample
label: [sample, sequence_type]
color: sequence_type
```

## AMR genes

AMRFinderPlus hits of type AMR: the count per genome, then the β-lactam genes per genome.

```query amr
source: amrfinder
filters:
  element_type: AMR
```

```sql summary
select m.sample, m.sequence_type, count(a.sample) as amr_hits
from mlst m left join amr a on a.sample = m.sample
group by all
order by amr_hits desc, m.sample
```

```datatable summary_table
source: summary
```

```amrfinder_matrix beta_lactams
source: amrfinder
filters: {element_class: BETA-LACTAM}
```
````

## Building a report

1. Learn each data file's columns and values: build a short file with just its `input` block and a
   `datatable` over it, or send the `input` block to `render_blocks`. The column names come from the adapter,
   not from the file: see [Adapters](#adapters).
2. Write the report.
3. Build it, from the groundie checkout:

   ```
   uv run groundie build cohort.md -o cohort.html
   ```

   - It reads a file only: there is no inline input.
   - Each block prints one line: `cohort.md:58: amrfinder_matrix beta_lactams: 64 rows; columns:
     element_symbol, ...`. A tree prints its tips: `cohort.md:28: tree phylogeny: 12 tips, 12 matched`.
   - A warning goes to stderr and the report is still built. `0 rows` usually means a wrong filter value.
   - An error stops the build with exit code 1 and one line on stderr:
     `groundie: error: cohort.md:47: sql 'summary': ...`
4. Read the HTML to check the values. Each datatable is a `<table>` under its heading; skip `<style>` and
   `<svg>`. Take every number and claim in the prose from a built datatable, never from counting rows in the
   data files.

To look something up in the data while writing, add a `sql` block and a `datatable` over it, build, read the
datatable, then remove both.

## Rendering blocks in chat

The MCP tool `render_blocks` runs blocks and shows every datatable, matrix and tree among them to the user in
the chat, in order, with the prose between them. Pass one of:

- `path`: the absolute path of a Markdown file of blocks. Relative input paths resolve against its directory.
- `blocks`: the Markdown text itself, with optional `base_dir`, the absolute directory that relative input
  paths resolve against. Without `base_dir`, input paths must be absolute.

The text result, which only you see, holds:

- one summary line per block, as the CLI prints them, with `blocks` as the file name for inline text
- the warnings
- the rows, as TSV, of each datatable and matrix (an `mlst_table` and an `amrfinder_matrix` included: their
  query's rows), and of the last block
  unless it is a tree: an `input`, `sql`, `query`, `clusters`, `checks` or `calls`. At most 200 rows per
  block; the rest are counted, not listed.

Each call is self-contained: nothing is remembered between calls. Every call carries the `input` blocks and
any other blocks the displayed blocks read, such as a `query`, `sql` or `calls`. To learn a data file's
columns and values, send just its `input` block: its rows come back as TSV.

In chat, write the prose around the blocks as your own chat message, not inside `blocks`. The same rule
holds: take every number from a rendered result, not from your memory of the files.

A short call, with `base_dir` set to the directory that holds `data`:

````markdown
```input mlst
adapter: mlst
path: data/mlst/mlst.tsv
```

```datatable st258
source: mlst
filters:
  sequence_type: 258
columns: [sample, sequence_type]
sort: [sample]
```
````

The text result:

```
blocks:1: input mlst: 12 rows; columns: sample, scheme, sequence_type, locus_1, locus_2, locus_3, locus_4, locus_5, locus_6, locus_7
blocks:6: datatable st258: 4 rows; columns: sample, sequence_type

datatable st258:
sample	sequence_type
KP02	258
KP05	258
KP09	258
KP11	258
```
