---
name: pleiotropy-condfdr
description: Conditional and conjunctional false discovery rate (condFDR, conjFDR) between two GWAS with GeNor (genor pleio run), using the pleioFDR Python port. Use when the user wants to boost locus discovery in one trait using a genetically related trait, find loci shared by two traits (pleiotropy, shared loci, polygenic overlap), or asks about pleioFDR, condFDR or conjFDR setup, reference data and memory.
---

# condFDR and conjFDR (pleioFDR)

## When to use

- **condFDR** (trait1 conditional on trait2): more discoveries for trait1 by using trait2's
  association signal as prior information (Andreassen et al. 2013).
- **conjFDR**: loci associated with **both** traits (Smeland et al. 2020). Useful when rg is near
  zero but the traits may still share loci with mixed effect directions.

Not for causal direction (use `mendelian-randomisation`) and not for showing a shared causal
variant at a locus (use `colocalisation`).

## Inputs and prerequisites

- Two genome-wide traits loaded with `genor sumstats load`, **EUR ancestry**, with **rsIDs**
  (matching to the template is by rsID; build is not used, but rsIDs must be dbSNP rsIDs).
- The `pleiofdr` backend: the pleioFDR Python port (espenhgn/pleiofdr, branch `python-port`,
  pinned commit 45e03b9), which needs **Python >= 3.13**. GeNor must run in an environment that
  can import it:

```bash
uvx --python 3.13 \
  --with 'pleiofdr @ git+https://github.com/espenhgn/pleiofdr@45e03b9c98c48ce9ef1496c4173a5581b23a853d' \
  genor doctor --json
# or install once: uv tool install genor --python 3.13 --with '<same spec>'
```

  `genor doctor` / `genor backend list` then show `pleiofdr` as available. For MCP, the server must
  be started from that environment.
- Reference data, both large (checksum-pinned in GeNor's registry):
  - `pleiofdr-template`: 9,545,380-SNP template, 262 MB.
  - `pleiofdr-ref-eur`: 1000G EUR LD matrix (r2 > 0.1) and annotations, 2.2 GiB (2.38 GB).
    Loading it needs about **12 GB of RAM**.

```bash
genor refdata fetch pleiofdr-template pleiofdr-ref-eur --yes
```

  MCP: `refdata_fetch(name, confirmed_by_user=true)`, only after the user agrees.

## Steps

1. Run `genor doctor --json`. Confirm the backend, the two datasets and the machine's memory.
2. Ask before fetching (see Decision points).
3. Optional dry run that only writes inputs (still requires both datasets, not the backend):

```bash
genor pleio run runs/scz runs/bip --stat condfdr --prepare-only --out runs/condfdr_scz_bip --json
```

4. Run (long; submit as a job):

```bash
genor job submit -- pleio run runs/scz runs/bip --stat condfdr --out runs/condfdr_scz_bip
genor job submit -- pleio run runs/scz runs/bip --stat conjfdr --out runs/conjfdr_scz_bip
genor job status <job_id>; genor job result <job_id>
```

   MCP: `run_analysis("pleio", {"trait1": ..., "trait2": ..., "stat": "conjfdr", "out": ...})`,
   then `job_status` / `job_result`.
5. For condFDR in both directions, run again with the traits swapped.

Options: `--stat condfdr|conjfdr`, `--fdr` (default 0.01 condFDR, 0.05 conjFDR),
`--randprune-n 500`, `--seed 1`, `--keep-mhc`, `--prepare-only`.

## Defaults and why

| Setting | Default | Why |
|---|---|---|
| condFDR threshold | 0.01 | Convention in the pleioFDR literature for conditional discoveries. |
| conjFDR threshold | 0.05 | Convention for shared loci (Smeland et al. 2020); conjFDR is the maximum of the two condFDRs, so it is already conservative. |
| Random pruning | 500 iterations | FDR estimates average over random LD-pruned SNP subsets to limit LD-driven inflation; published analyses use >= 500 (GeNor warns below 100). |
| MHC | excluded (chr6:25,119,106-33,854,733, GRCh37) | Extreme, long-range LD distorts the enrichment estimates. `--keep-mhc` keeps it. |
| Seed | 1 | Random pruning is stochastic; fixed for reproducibility. |

Background and interpretation detail: [references/method-and-interpretation.md](references/method-and-interpretation.md).

## Decision points

- **Ask the user** before downloading: 262 MB + 2.2 GiB, and about 12 GB RAM at run time. If the machine has less memory, stop.
- **Ask the user** whether the two GWAS share participants or controls. Overlap can create
  apparent shared signal; GeNor does not correct for it. Show the LDSC cross-trait intercept.
- **Ask the user** which trait is primary for condFDR (trait1 | trait2), or whether both
  directions are wanted.
- **Ask the user** before changing `--fdr` or keeping the MHC; report the thresholds used.
- **Ask the user** which shared loci to follow up (coloc, Open Targets).
- Non-EUR GWAS: the only reference is EUR. Stop and say so rather than substituting.

## Reading the output

`summary`: `stattype`, `fdrthresh`, `randprune_n`, `matched` (SNPs matched to the template per
trait), `config` (path to the pleioFDR `config.txt`), `n_loci`, `result_files`, `loci_table`
(path to the `*_loci.csv` table, when written). Everything is under `<out>/pleiofdr/`: the MAT
inputs, `config.txt`, and `results/` with pleioFDR's tables and plots. With `--prepare-only`,
`prepared_only: true` and no results.

- **condFDR** value per SNP: the posterior probability that the SNP is null for trait1 given its
  p-values in both traits. condFDR < 0.01 is a conditional discovery for trait1.
- **conjFDR** < 0.05: the SNP is associated with both traits. Check the sign of the two z-scores
  at each lead SNP: concordant and discordant shared loci are both common.
- Compare `matched` with the trait's `n_variants`. Low overlap (`pleio_low_overlap`) means many
  SNPs were not used.

## Failure signatures

| Code | Meaning | Fix |
|---|---|---|
| `refdata_missing` | template or reference not downloaded | Ask; fetch with `--yes`. |
| `confirm_download` | fetch without `--yes` | Ask; add `--yes`. |
| `backend_missing` | pleioFDR port not importable | Run through `uvx --python 3.13 --with 'pleiofdr @ ...'` as above. |
| `pleio_needs_rsid` | trait has no `variant_id` | Load a file with rsIDs. |
| `needs_user_decision` (exit 3) | ancestry unknown or multi | Ask; reload with `--ancestry`. |
| `ancestry_mismatch` | ancestry is not EUR | No matching reference; stop. |
| `pleio_same_trait_name` | both traits give the same file name (e.g. `MDD-2019` vs `MDD_2019`) | Reload one trait with a different `--trait`. |
| `bad_param` | `--stat` not condfdr/conjfdr | Fix the option. |
| `internal_error` / killed job | often out of memory loading the 2.2 GiB reference | Check `genor job logs`; needs about 12 GB RAM. |

Warnings: `pleio_low_overlap` (serious: < 30% of a trait's variants matched the template; rsIDs
missing or a non-dbSNP ID scheme), `pleio_allele_mismatch` (alleles not matching the template set
to missing), `pleio_few_prunes` (`--randprune-n` < 100).

## What not to conclude

- conjFDR loci are shared at the locus level, not necessarily the same causal variant. Use
  `colocalisation` before claiming a shared variant.
- Shared loci are not evidence of causation in either direction.
- condFDR discoveries for trait1 still need replication like any GWAS hit; they depend on the
  choice of conditioning trait.
- conjFDR loci with discordant effect directions are expected and do not mean an error.
- Overlapping samples or shared controls can produce shared loci that are artefacts.

## Related skills

`sumstats-qc`, `heritability-rg` (rg and overlap), `colocalisation`, `fine-mapping`,
`open-targets-triage`, `interpretation-and-reporting`.
