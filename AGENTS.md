# AGENTS.md

Guidance for coding agents working **on** GeNor. (Agents *using* GeNor should read
`skills/genor/SKILL.md`.)

## Setup and checks

```bash
uv sync                       # creates .venv with dev deps
uv run pytest                 # fast suite (offline, ~15 s)
uv run pytest -m slow         # golden LDSC test (downloads 25 MB)
uv run pytest -m live         # hits real APIs; nightly only
uv run ruff check src tests scripts
uv run python scripts/gen_cli_docs.py   # regenerate docs/cli.md, docs/mcp.md, docs/data-sources.md
```

Run all four (pytest, ruff, docs regeneration, `git diff docs/`) before finishing a change.

## Rules

- **Summary statistics only.** Never add code that reads individual-level genotypes as study input.
  The refusal in `guards.py` stays.
- **No authenticated fetching.** No tokens, accounts or click-through terms in code. New data
  sources go in `sumstats/sources.toml` as pointers.
- **New native methods need a golden test** against the reference implementation
  (`tests/golden/`, generator in `scripts/golden/`). Deterministic quantities must match to
  numerical precision.
- **Pinned downloads.** New refdata or backend entries carry a checksum. Unpinned entries need
  `--allow-unpinned` and must say why.
- **Don't bundle non-redistributable tools** (MAGMA, METAL, LAVA). Drive user-installed copies.
- **Every command goes through `Run`** (`runs.py`): record inputs, refdata, backends, citations,
  warnings (stable `code`, `severity`), and `ask_user` for decisions that belong to the human.
- **Errors are actionable**: raise `GenorError(message, code=..., hint=...)`, or `NeedsUserDecision`
  when only the user can answer. Never guess build, ancestry, N, prevalence or sample overlap.
- **CLI, MCP and skills stay in sync.** A new CLI option or warning code needs: the MCP mapping
  (`mcp_server.METHODS` or a tool), the regenerated docs, and the relevant skill's
  "Failure signatures" / "Defaults" sections.
- Keep skills under 500 lines; depth goes in `references/`.
- Python ≥ 3.11 (no 3.12-only syntax). Linux and macOS only.

## Layout

See `docs/architecture.md`. The vendored LDSC code in `src/genor/methods/_ldsc/` is unmodified;
don't edit it (excluded from ruff).
