# Contributing

Thanks for helping. GeNor is research software; correctness comes first.

1. Open an issue describing the method or fix, with references.
2. `uv sync`, make the change, add tests. Native statistical methods need a golden test against
   the reference implementation (see `docs/validation.md`).
3. `uv run pytest && uv run ruff check src tests scripts && uv run python scripts/gen_cli_docs.py`
4. Update the relevant skill in `skills/` (defaults, decision points, failure signatures).
5. Open a pull request. By contributing you agree your code is released under GPL-3.0-or-later.

Coding agents: see `AGENTS.md`.
