# Reproducibility with GeNor manifests

Every analysis command writes `<run_dir>/manifest.json` and `<run_dir>/result.json` (source:
`src/genor/runs.py`), including failed runs.

## manifest.json fields

| Field | Content |
|---|---|
| `genor_version` | GeNor version |
| `command`, `argv` | Command name (e.g. `coloc.abf`) and the full argument vector |
| `status` | `ok` or `error` |
| `started`, `finished` | UTC timestamps |
| `platform` | Python version, OS, machine |
| `params` | All parameters including defaults |
| `inputs` | Each input file: `role`, absolute `path`, `size`, `sha256`, and trait name where relevant. Files over 4 GB are not hashed (`sha256` null) |
| `refdata` | Reference datasets used, with build, ancestry, licence |
| `backends` | External tools used, with path and version |
| `meta` | Trait metadata (build, ancestry, N, ...) where recorded |
| `summary` | Same summary as `result.json` |
| `warnings` | Structured warnings (`code`, `message`, `severity`) |
| `ask_user` | Decisions left to the human |
| `outputs` | Files written |
| `citations` | Citation keys and full texts |
| `error` | Error code and message if the run failed |
| `notice` | Research-use notice |

## Re-running

- `argv` is the exact command line; run it again in the same working directory with the same
  inputs. Relative paths in `params` are relative to where GeNor was run.
- Verify inputs before re-running: `shasum -a 256 <file>` (macOS) or `sha256sum <file>` (Linux)
  must match `inputs[].sha256`. Note that downstream runs hash the trait's `sumstats.parquet`;
  the original file's hash is in the `sumstats load` run directory (`source_sha256` in trait
  metadata and `inputs`).
- Stochastic steps have seeds: MR bootstrap (`--seed`, `--nboot`) and pleioFDR random pruning
  (`--seed`, `--randprune-n`). Fine-mapping, coloc, LDSC, MAGMA and meta-analysis are
  deterministic given inputs, reference data, and tool versions.
- Reference data is checksum-verified on download (`genor refdata fetch`).
- Live lookups are not reproducible: Open Targets, Ensembl, STRING and g:Profiler data change
  between releases. Save the JSON output and record the date and any version fields.

## What to hand over with results

1. The run directories (or at least `manifest.json`, `result.json`, and output tables).
2. `report.md` / `report.html`.
3. Saved lookup outputs with dates.
4. The GeNor version and how it was installed.
5. The source of every input GWAS (URL or accession, download date, terms of use). GeNor records
   paths and hashes, not provenance; ask the user.

Do not hand over summary statistics the user is not allowed to redistribute; share hashes and
access instructions instead.
