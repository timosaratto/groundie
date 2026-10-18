# groundie: declarative reporting framework for bacterial genomics

Groundie makes turning your pipeline outputs into user-friendly reports a breeze, for both LLMs and humans.

## How it works

Reports are written as declarative templates, with YAML blocks declaring which files to read and how to present the results. Groundie has prebuilt components for common tasks like AMR and phylogeny, and generic query and UI components for complex queries or other tools.

Templates are then compiled into user-friendly HTML reports, with full provenance of the queries that produced them, down to input file hashes. Reports can also be directly embedded into LLM apps like Claude Cowork using groundie's MCP App, enabling LLMs to report their findings in an easily reviewable format.

**Note: groundie is still in its early development stages – functionality is minimalistic, and bugs and breaking changes are to be expected. All feedback is highly appreciated.**

## Examples

For example, this template

````
## Phylogeny

```input tree
adapter: newick
path: ../data/tree_results/cohort.nwk
```

```Tree phylogeny
source: tree
```
````

...will be rendered as this report<img width="1120" height="479" alt="CleanShot 2026-10-18 at 17 21 59" src="https://github.com/user-attachments/assets/fe2aaafb-70b6-4968-b229-d9c3834b11a7" />
<img width="1135" height="444" alt="CleanShot 2026-10-18 at 17 12 33 1" src="https://github.com/user-attachments/assets/9cc25a44-e802-47ae-8575-e8b977c6e7fa" />

While this template

````
## Phylogeny

```input tree
adapter: newick
path: ../data/tree_results/cohort.nwk
```

```input mlst
adapter: mlst
path: ../data/mlst_results/SAMD*.tsv
```

```Tree phylogeny
source: tree
tips: mlst
match: sample
label: [sample, sequence_type]
color: sequence_type
```
````

...will be rendered like this
<img width="1120" height="479" alt="CleanShot 2026-10-18 at 17 21 59" src="https://github.com/user-attachments/assets/1a94c82b-71b0-483f-9947-560bbb0737de" />

## Using groundie

See [GUIDE.md](src/groundie/GUIDE.md) for instructions for creating groundie templates.

## Supported tools

groundie currently supports a minimalistic set of tools. It has built-in query connectors for amrfinder, mlst, mash and newick plus generic connectors for any tsv, csv or json file. It supports rendering tables, heatmaps and phylogenetic trees.

## Quick start with groundie CLI

```
brew install uv   # or: curl -LsSf https://astral.sh/uv/install.sh | sh
git clone https://github.com/timosaratto/groundie.git
cd groundie
uv run groundie build examples/cohort.md -o cohort.html
open cohort.html
```

This builds the example report in [examples/cohort.md](examples/cohort.md) into a self-contained HTML page. The first run builds the environment from `uv.lock`, so it may take a moment.

While writing a template, `uv run groundie dev examples/cohort.md` serves the report at http://localhost:8000 and reloads it whenever the file changes.

## Quick start with Claude Code

```
brew install uv   # or: curl -LsSf https://astral.sh/uv/install.sh | sh
git clone https://github.com/timosaratto/groundie.git
cd groundie
claude mcp add groundie -- uv run --directory "$PWD" groundie-mcp
claude
```

Then give Claude the example prompt below. The first run builds the environment from `uv.lock`, so it may take a moment.

> Run outbreak investigation for the 10 S. aureus genomes in the `data` directory. Format your findings using groundie

Then try the same prompt without groundie, and you'll notice a difference.

## Quick start with Claude Desktop

1. Install uv and clone the repo as in the Claude Code quick start above.
2. In Claude Desktop, open **Settings → Developer → Edit Config** and add:
   ```json
   {
     "mcpServers": {
       "groundie": {
         "command": "</absolute/path/to/uv>",
         "args": [
           "run",
           "--directory",
           "<absolute/path/to/groundie>",
           "groundie-mcp"
         ]
       }
     }
   }
   ```
   Both paths must be absolute, since Desktop doesn't use your shell's `PATH`. Find uv with `which uv`.
3. Quit Claude Desktop completely and reopen it.
4. Run the example prompt from the Claude Code quick start in Claude Cowork. (the regular chat has no file system access)
