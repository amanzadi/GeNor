# What `genor report build` produces

Source: `src/genor/report/builder.py`.

```
# <title>

> GeNor is research software. It is not a medical device and must not be used for clinical
> decisions, diagnosis, or commercial services.

## <command>                       one section per run directory, in the order given
Run directory: `<dir>` · status **ok|error** · GeNor <version> · <finished or started time>

**Inputs**
- <role>: `<absolute path>` (sha256 <first 12 hex>…)

**Reference data and tools**
- <refdata name>: <build> <ancestry> (<licence>)
- <backend name> <version>

**Parameters:** `key=value`, ...    (None values omitted)

| Result | Value |                  flattened run summary, nested two levels; lists/dicts
|---|---|                           longer than 200 characters are truncated; first 60 rows

**Warnings**
- [severity] `code`: message

**Decisions for the human**
- question (why)

![key](path.png)                    every PNG output (Manhattan, QQ, MR scatter/forest/funnel, ...)

### STROBE-MR checklist             MR runs only (from strobe_mr.md)

**Error:** `code`: message          failed runs only

## References
1. GeNor: post-GWAS analysis toolkit for AI agents. ...
2. ...                              deduplicated across runs
```

The HTML version is a direct conversion with minimal styling (tables, lists, images, code).

## Practical notes

- **Paths**: figure paths are relative only when the PNG is under the report's directory.
  Otherwise they are absolute, and the HTML breaks when moved. Write the report into the parent
  directory of the runs.
- **Summary truncation**: the table shows the first 60 flattened values. Full numbers are in each
  run's `result.json` and output files (e.g. `finemap.tsv`, `coloc_snp_pp.tsv`, `mr_results.tsv`).
- **coloc sections**: `trait1` / `trait2` hold the trait names; `prior_trait1` / `prior_trait2`
  hold prior details (`type`, `sd_prior`, `sdY_estimated`).
- **Lookups** (Open Targets, VEP, genes in region, g:Profiler, STRING) are not run directories and
  do not appear. Add a section by hand.
- **Order**: sections follow the order of the `RUN_DIRS` arguments. Put data loading first, then
  genome-wide analyses, then locus-level ones.

## Severity levels in warnings

| Severity | Meaning | In the write-up |
|---|---|---|
| `info` | Expected property of the method or a note (e.g. `single_causal_assumption`, `build_inferred`) | Mention in methods where relevant |
| `warning` | May affect results (e.g. `weak_instruments`, `mr_heterogeneity`, `sdy_assumed`) | Discuss in results or limitations |
| `serious` | Likely to invalidate or strongly bias a result (e.g. `build_conflict`, `allele_mismatch`, `mr_directional_pleiotropy`, `ldsc_intercept_high`, `rg_out_of_bounds`) | Must be in limitations; consider not interpreting the result |
