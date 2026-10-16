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
```DataTable st
source: mlst
columns: [sample, sequence_type]
```
````

Blocks that show something on the page are PascalCase, like JSX components; blocks that only make data are
lowercase.

| Kind | Shows | Can be a source for later blocks |
| --- | --- | --- |
| `input` | Nothing | Yes |
| `sql` | Nothing | Yes |
| `query` | Nothing | Yes |
| `clusters` | Nothing | Yes: one row per sample, its cluster |
| `checks` | Nothing | Yes: one row per source row, PASS or FAIL |
| `calls` | Nothing | Yes: one row per sample, true or false per call |
| `DataTable` | A datatable | Yes: the rows it shows |
| `Matrix` | A grid | Yes: its query's rows, not the pivoted grid |
| `AmrfinderMatrix` | A gene × sample grid of AMRFinderPlus hits | Yes: as a `Matrix` |
| `Tree` | A phylogeny | No |

The source of a block is the block it reads, named by its `source:` key. Every block except a `Tree` makes an
SQL table with the block's name, holding the rows described in the last column above. Later blocks read it
with `source: <name>`, or by that name inside a `sql` block.

The rules:

1. Every block has a kind and a name. A name is used once in the file, across blocks of all kinds.
2. `source:` names a block above it, or any `input`, wherever it is. So do a `Tree`'s `tips` and a calls
   block's `samples`.
3. A `query` is a `DataTable` that shows nothing. Use it for rows that later blocks read but the reader need
   not see.

Other fences are shown as code. Give them at most one word after the backticks. A fence that starts with a
block kind but is not read as a block, such as ```` ```DataTable ```` with no name, is shown as code with a
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
- An input shows nothing. To show its rows, read it with a `DataTable`.

## Adapters

| Adapter | Files | SQL table |
| --- | --- | --- |
| `amrfinder` | AMRFinderPlus TSV, version 4 | One row per hit: the columns listed below. `NA` is read as null |
| `mlst` | mlst's default output: no header; the file, the scheme, the ST, then `locus(call)` for each of the scheme's loci | One row per genome: `sample` (the file, its directories and a `.fa`, `.fasta`, `.fna` or `.fas` extension, with or without `.gz`, stripped), `scheme`, `sequence_type` (a number), and `calls`, a list with one `{locus, call}` per locus of the scheme, in mlst's order. `scheme` and `sequence_type` are null where mlst wrote `-`: no scheme found, no ST assigned. A call is text, as mlst wrote it with its flags: mlst's `tpi(~4)` is `{locus: tpi, call: ~4}`; other calls look like `3`, `4?` and `4,5`. A call is null where mlst wrote `-`, no allele found: `gmk(-)` is `{locus: gmk, call: null}`. A genome with no scheme found has no calls. See the SQL section to read `calls` |
| `mash` | `mash dist` or `mash triangle -E` output, no header | `reference`, `query` (sample names: the directories and a `.fa`, `.fasta`, `.fna`, `.fas` or `.msh` extension, with or without `.gz`, are stripped), `distance`, `p_value`, `shared_hashes`, `sketch_hashes` |
| `tsv`, `csv` | Any delimited file | `file` (the file's name), then the file's columns as named. Numbers are read as numbers |
| `json` | Newline-delimited objects, or one array of objects | `file`, then one column per key |
| `newick` | Newick trees | One row per node: `file`, `tree`, `node`, `parent`, `label`, `branch_length`, `support`, `is_tip`, `comment`. Drawn by a `Tree` block |

The column names come from the adapter, not from the file.

Titles and descriptions, which a `DataTable` shows with `display: <input>`: `mlst` titles `sample` Sample,
`scheme` Scheme and `sequence_type` ST. `amrfinder` titles each column by AMRFinderPlus's own header, such as
`Element symbol` and `% Identity to reference`, and `sample` Sample. Both describe each of their columns. No
other adapter declares titles or descriptions yet.

The `amrfinder` columns:

```
text:    sample, protein_id, contig_id, strand, element_symbol, element_name, scope, element_type,
         element_subtype, element_class, element_subclass, method, ref_accession, ref_name, hmm_accession,
         hmm_description, hierarchy_node
numbers: start_pos, stop_pos, target_length, ref_length, pct_coverage, pct_identity, alignment_length
```

`sample` is the file's `Name` column, written by `--name`. A file without it takes `sample` from its file name,
without directories and the last extension: `amrfinder/KP01.tsv` is `KP01`. Two such files that give the
same sample, as `KP01/amrfinder.tsv` and `KP02/amrfinder.tsv` do, fail the input: run with `--name`. `element_type`,
`element_subtype`, `element_class` and `element_subclass` are its `Type`, `Subtype`, `Class` and `Subclass`.
Files from versions before 4.0, which name most columns differently, are not read: the input fails, naming
the column it can't find.

Null where AMRFinderPlus writes `NA`: `protein_id` for nucleotide input; `contig_id`, `start_pos`, `stop_pos`
and `strand` for protein input; `element_class` and `element_subclass` for virulence hits; `ref_length`
through `ref_name` where no reference was aligned; `hmm_accession` and `hmm_description` unless an HMM
matched. `hierarchy_node` is null also where AMRFinderPlus ran without `--print_node`.

Only `tsv` and `csv` take `options`:

- `null`: the string read as null
- `header`: default true. With `header: false` the columns are named `column0`, `column1`, and so on

## Query keys: query, DataTable, Matrix

`query`, `DataTable` and `Matrix` blocks share these keys. `source` is required; the rest are optional.

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

- These are all the query keys. Any other key is an error, except a `Matrix`'s layout keys and a
  `DataTable`'s `display`, `spread` and `split`.
- Filter values are written as text and read as the column's type: `sequence_type: 258` matches the number
  258, and `sample: 007` the text `007`, not `7`.
- `null` (or `~`) matches a missing value: `sequence_type: null` keeps the genomes with no ST, and
  `sequence_type: [8, null]` those with ST 8 or none. `None` is the text `None`.
- Joins, expressions, ranges and descending order need a `sql` block.
- A block reading a `sql` block, with no `sort`, keeps the order of the statement's `order by`.

### display

A `DataTable`, and only a `DataTable`, takes `display`: for some of its columns, a title, shown as the header in
place of the column's name, and a description, shown on hovering the header.

```yaml
source: mlst
columns: [sample, scheme, sequence_type]
display:
  sequence_type: ST        # the header shows ST; hovering it shows sequence_type, the column's name
  scheme: {title: Scheme, description: The PubMLST scheme mlst matched}
  sample: {description: The genome's file name, without its extension}   # the header stays sample
```

- Each key is one of the block's output columns. Its value is the title, as text, or `{title, description}`,
  either or both, each text.
- Hovering a header with a description shows the description. Hovering one with a title and no description
  shows the column's name.
- `display` may instead name an input, whose adapter's titles and descriptions it shows: `display: mlst`
  titles `sample`, `scheme` and `sequence_type` as `Sample`, `Scheme` and `ST`, and describes each column as
  the adapter's schema does. It matches by column name, so it works through a `query` or `sql` block too, as
  long as the name is unchanged; a column the adapter has no title for keeps its name. Make sure a column
  computed under an adapter's name still means what its description says. Only adapters that declare titles
  or descriptions can be named: see [Adapters](#adapters).
- Only the page changes. The SQL table the block makes, later blocks that read it, and the rows
  `render_blocks` returns keep the column's name: `sequence_type`, not `ST`.

### spread

A `DataTable`, and only a `DataTable`, takes `spread`: a list column of structs, such as an `mlst` input's
`calls`, shown as one column per name, each holding that name's value.

```yaml
source: mlst
columns: [sample, sequence_type, calls]
spread: {column: calls, names: locus, values: call}   # a column per locus, holding its call
display: mlst
```

- `column` is one of the block's output columns, a list of structs. `names` is the field of each struct that
  names its column, `values` the field the cell shows. All three are required. One `spread` per table.
- Each row stays one row. In place of `column` come the names, one column each, in the order they first appear:
  row by row, then along each list. Above, a cohort of two schemes shows the first scheme's loci, then the
  other's; to give each scheme its own table, add a [`split`](#split).
- Each header is the name, `arcC`. Hovering it shows `calls: arcC`. A `display` title or description on
  `calls` is not shown. `display` names output columns, so it can't title the spread columns.
- A cell:

  | The row's list | The cell |
  | --- | --- |
  | Has the name, with a value | The value |
  | Has the name, with a null value | – |
  | Lacks the name | Blank |
  | Is empty | Blank, in every spread column |
  | Is null | –, in every spread column |

- When no row's list has an entry, there are no spread columns, and `column` is simply gone.
- Errors: `column` not an output column; a cell that is neither null nor a list of structs with both fields; a
  null name; a name twice in one row's list; a name that is another column's name.
- Only the page changes. The SQL table the block makes, the summary line, later blocks that read it, and the
  rows `render_blocks` returns keep `column` as the list. 'View query' shows `spread` as written.

### split

A `DataTable`, and only a `DataTable`, takes `split`: one of its output columns, whose values split its rows into
one table each, captioned with the value. mlst results, one table per scheme, each headed by its own loci:

```yaml
source: mlst
columns: [sample, scheme, sequence_type, calls]
split: scheme                                         # a table per scheme
spread: {column: calls, names: locus, values: call}   # in each, a column per locus of that scheme
display: mlst
sort: [scheme, sample]
```

- The tables, and the rows in each, come in the order their first row appears: give `sort` to set both. A
  null is a value like any other.
- Each table's caption is the column's title and the value: `Scheme: saureus`. The title is the column's
  `display` title, or its name. A null shows as `Scheme: –`.
- The split column stays in each table, as an ordinary column.
- With `spread`, each table is spread on its own, so it gets only its own names: a cohort of two schemes
  gives two tables, each with its scheme's loci. A row that lacks a name another row of its table has shows
  blank, as `spread` does.
- With no rows, there is one empty table, with no caption.
- Errors: `split` not one output column; a cell that is a list or struct.
- The SQL table the block makes, later blocks that read it, and the rows `render_blocks` returns stay one
  table. Besides the page, only the summary line changes: it counts the tables, `10 rows in 2 tables`. 'View
  query' shows `split` as written.

## sql

One SELECT statement in DuckDB SQL. It reads the SQL tables of inputs and of blocks above it, by name, and has
no access to files or the network. Any DuckDB SELECT works, including `with` and `with recursive`, window
functions and `filter`. It shows nothing, so follow it with a `DataTable` that reads it:

````markdown
```sql summary
select m.sample, m.sequence_type, count(a.sample) as amr_hits
from mlst m left join amr a on a.sample = m.sample
group by all
order by amr_hits desc, m.sample
```

```DataTable summary_table
source: summary
```
````

A list column, such as an `mlst` input's `calls`, is read with `unnest`, one row per element. Every call
that is not a plain allele number:

```sql
select m.sample, c.locus, c.call
from mlst m, unnest(m.calls) as u(c)
where c.call is null or not regexp_full_match(c.call, '\d+')
```

A `DataTable` shows a list or struct cell as raw Python, e.g. `[{'locus': 'arcC', 'call': '3'}, ...]`, so a
`DataTable` over an `mlst` input with no columns chosen shows `calls` that way. To show a list of structs as
one column per name, give the `DataTable` a [`spread`](#spread). Otherwise select its fields, as above, or
format it as text in SQL. Each genome's calls as one cell, `arcC:3, aroE:1, ...`:

```sql
select sample, scheme, sequence_type,
       array_to_string(list_transform(calls, c -> c.locus || ':' || coalesce(c.call, '-')), ', ') as calls
from mlst
```

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
- It shows nothing. To show it, read it with a `DataTable`, or join it in a `sql` block.

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
- It shows nothing. To show it, read it with a `DataTable`, or join it in a `sql` block.

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
- It shows nothing. To show it, read it with a `DataTable`, or join it in a `sql` block.

## Matrix

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
values are in its `title` attribute in the HTML, shown on hover. A block that reads a `Matrix` gets the query's
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

## AmrfinderMatrix

A `Matrix` with its layout set for AMRFinderPlus hits: one row per gene, one column per sample.

```yaml
source: amrfinder                         # an input loaded by the amrfinder adapter
filters: {element_class: AMINOGLYCOSIDE}  # optional: any filters, e.g. element_type, element_class, scope, sample
partial_below: 80                         # optional: ◆ marks coverage below this (%); default 100
```

- It takes `source`, `filters`, `sort` and `partial_below`. For any other layout, write a `Matrix`.
- Rows are `element_symbol`, with `element_type` before it unless a filter picks one type. Columns are
  `sample`. The order is a `Matrix`'s default: the most common genes first.
- A cell with several hits shows the most complete one (`pick: max(pct_coverage)`); hovering gives its
  identity and coverage. ◆ marks a cell whose most complete hit is still partial (coverage below
  `partial_below`).
- 'View query' shows the `Matrix` it stands for, every setting written out.

## Tree

```yaml
source: tree                    # an input loaded by the newick adapter
tips: mlst                      # optional: a block with one row per tip; any kind but Tree, an input included
match: sample                   # needed with tips: the column that equals the tip label, as exact text
label: [sample, sequence_type]  # columns shown at the tips; default: the tip label
color: sequence_type            # the 3 most frequent values get a colour each, the rest share a grey
```

- These are all the `Tree` keys. A `Tree` takes no query keys: to filter or join the tip rows, do it in a `query`
  or `sql` block and name that block in `tips`.
- `label` and `color` need `tips`.
- When the input holds several trees, pick one with `file: core.nwk` and `tree: 1`, the tree's position in
  that file.
- A `Tree` makes no SQL table. No block can read it.

A tip with no matching row is drawn and reported, not an error. Check the summary line:
`Tree phylogeny: 12 tips, 12 matched`.

## What the output adds

The output adds these itself, so leave them out of the prose:

- under each `DataTable`, `Matrix` (an `AmrfinderMatrix` included) and `Tree`: its filters, its query, and its
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

One row per genome, from mlst. A – in the ST column means none was assigned.

```DataTable st
source: mlst
columns: [sample, sequence_type]
sort: [sample]
display: mlst
```

```Tree phylogeny
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

```DataTable summary_table
source: summary
```

```AmrfinderMatrix beta_lactams
source: amrfinder
filters: {element_class: BETA-LACTAM}
```
````

## Building a report

1. Learn each data file's columns and values: build a short file with just its `input` block and a
   `DataTable` over it, or send the `input` block to `render_blocks`. The column names come from the adapter,
   not from the file: see [Adapters](#adapters).
2. Write the report.
3. Build it, from the groundie checkout:

   ```
   uv run groundie build cohort.md -o cohort.html
   ```

   - It reads a file only: there is no inline input.
   - Each block prints one line: `cohort.md:58: AmrfinderMatrix beta_lactams: 64 rows; columns:
     element_symbol, ...`. A `Tree` prints its tips: `cohort.md:28: Tree phylogeny: 12 tips, 12 matched`.
   - A warning goes to stderr and the report is still built. `0 rows` usually means a wrong filter value.
   - An error stops the build with exit code 1 and one line on stderr:
     `groundie: error: cohort.md:47: sql 'summary': ...`
4. Read the HTML to check the values. Each datatable is a `<table>` under its heading; skip `<style>` and
   `<svg>`. Take every number and claim in the prose from a built datatable, never from counting rows in the
   data files.

To look something up in the data while writing, add a `sql` block and a `DataTable` over it, build, read the
datatable, then remove both.

## Rendering blocks in chat

The MCP tool `render_blocks` runs blocks and shows every `DataTable`, `Matrix` and `Tree` among them to the user in
the chat, in order, with the prose between them. Pass one of:

- `path`: the absolute path of a Markdown file of blocks. Relative input paths resolve against its directory.
- `blocks`: the Markdown text itself, with optional `base_dir`, the absolute directory that relative input
  paths resolve against. Without `base_dir`, input paths must be absolute.

The text result, which only you see, holds:

- one summary line per block, as the CLI prints them, with `blocks` as the file name for inline text
- the warnings
- the rows, as TSV, of each `DataTable` and `Matrix` (an `AmrfinderMatrix` included: its query's rows), and
  of the last block
  unless it is a `Tree`: an `input`, `sql`, `query`, `clusters`, `checks` or `calls`. At most 200 rows per
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

```DataTable st258
source: mlst
filters:
  sequence_type: 258
columns: [sample, sequence_type]
sort: [sample]
```
````

The text result:

```
blocks:1: input mlst: 12 rows; columns: sample, scheme, sequence_type, calls
blocks:6: DataTable st258: 4 rows; columns: sample, sequence_type

DataTable st258:
sample	sequence_type
KP02	258
KP05	258
KP09	258
KP11	258
```
