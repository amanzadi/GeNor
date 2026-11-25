# Architecture

One package, three layers. The CLI is the core; the MCP server and skills sit on top of it.

```
src/genor/
  cli.py            Typer app; every command can print JSON and writes a run directory
  commands.py       command implementations shared by the CLI and MCP server
  mcp_server.py     fastmcp server; analyses run as background jobs
  runs.py           run directory + manifest.json (inputs, sha256, params, versions, citations)
  jobs.py           detached subprocess jobs: submit / status / result / logs / cancel
  errors.py         GenorError (code + hint), NeedsUserDecision (exit 3), structured warnings
  guards.py         summary-statistics-only refusal, build/ancestry checks
  config.py         platform data/config/state dirs (GENOR_HOME overrides)
  skills.py         locate, serve (get_guidance) and install skills
  sumstats/         columns.py, load.py, qc.py, harmonise.py, liftover.py, sources.toml
  refdata/          registry.toml + checksum-pinned lazy downloads
  backends/         registry.toml + resolver/installer for external tools
  methods/          meta, ldsc (+ vendored _ldsc), clump, mr, pleio, coloc (ABF), magma
  lookup/           http.py (rate limits, retries, test transport), services.py
  report/           plots.py, strobe_mr.py, builder.py (Markdown + HTML)
  data/             build_markers.tsv.gz, example/ (simulated GWAS pair)
skills/             SKILL.md trees (Claude Code plugin reads them here; copied into the wheel)
.claude-plugin/     plugin.json + marketplace.json
.mcp.json, bin/, commands/   plugin MCP config, PATH wrapper, slash commands
```

## Conventions

- `genor <area> <verb> ... --out <run_dir> --json`. Stdout carries one short JSON object
  (`status`, `summary`, `warnings`, `ask_user`, `files`, `next_steps`); tables and figures go to
  the run directory.
- Exit codes: 0 ok, 2 actionable error (`error.code`, `error.hint`), 3 the human must decide
  (`error.ask_user`), 1 unexpected failure.
- A loaded trait is a directory with `sumstats.parquet` (GWAS-SSF column names) and
  `sumstats.json` (build, ancestry, N, cases/controls, prevalences, column mapping, QC report).
  All analyses take trait directories, never raw files.
- Warnings have a stable `code` and a `severity` (`info`, `warning`, `serious`); skills map codes
  to meanings and fixes.
- The MCP server never blocks longer than ~50 s: it submits a job, waits, and returns either the
  result or a job id (Codex's default tool timeout is 60 s).

## Decisions

- **Summary statistics only** and **no authenticated fetching**: GeNor never handles credentials
  or individual-level data; it points users to sources.
- **uv core + opt-in backends**: everything in core works from `uvx genor`; tools that are R-only,
  container-only or not redistributable are optional and detected by `genor doctor`.
- **Native ports only with golden tests** against the reference implementation.
- **GPL-3** because LDSC is vendored.
