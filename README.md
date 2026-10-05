# groundie – declarative reporting framework for microbial genomics

Agents are becoming good at bioinformatics, but do you trust their results aren't hallucinated? I don't. That's why I built groundie.

Groundie lets your LLM generate HTML reports where the tables and figures are guaranteed to match the underlying data files.

Note: groundie is still in its early development stages, so bugs and breaking changes are to be expected. All feedback is highly appreciated.

## How it works

A human or an LLM writes a report in Markdown: prose plus query blocks that define how the underlying files (e.g. amrfinder and mlst TSVs) are queried and displayed. The groundie CLI builds that report into an HTML page by running the blocks against the files. All data are queried deterministically at compile-time, and a provenance chain is included in the report. You can also check the queries to confirm that the LLM hasn't messed them up, although it seldom does.

The nice thing is that the LLM also sees the query results, which allows it to write its remarks around the data tables. These are not fact-checked by groundie and can contain mistakes.

As a result, this allows you to vibe yourself bioinformatics results that you can actually trust.

Or, if you don't like AI, you can also write reports by hand, and add it at the end of your bioinformatics pipeline to get a neat report.

See [GUIDE.md](src/groundie/GUIDE.md) for instructions on using groundie.

## Supported tools

Groundie currently supports a minimalistic set of tools. It has built-in query connectors for amrfinder, mlst, mash and newick plus generic connectors for any tsv, csv or json file. It supports rendering tables, heatmaps and phylogenetic trees.

## Quick start with Claude Code

```
brew install uv   # or: curl -LsSf https://astral.sh/uv/install.sh | sh
git clone https://github.com/timosaratto/groundie.git
cd groundie
claude mcp add groundie -- uv run --directory "$PWD" groundie-mcp
claude
```

Then give Claude the example prompt above. The first run builds the environment from `uv.lock`, so it may take a moment.

> Run outbreak investigation for the 10 S. aureus genomes in the `data` directory. Format your findings using groundie

Then try the same prompt without groundie, and you'll notice a difference.

## Quick start with Claude Desktop

1. Install uv and clone the repo as in the Claude Code quick start below.
2. In Claude Desktop, open **Settings → Developer → Edit Config** and add:
   ```json
   {
     "mcpServers": {
       "groundie": {
         "command": </absolute/path/to/uv>,
         "args": ["run", "--directory", <absolute/path/to/groundie>, "groundie-mcp"]
       }
     }
   }
   ```
   Both paths must be absolute, since Desktop doesn't use your shell's `PATH`. Find uv with `which uv`.
3. Quit Claude Desktop completely and reopen it.
4. Run an example prompt in Claude Cowork. (the regular chat has no file system access)
