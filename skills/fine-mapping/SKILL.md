---
name: fine-mapping
description: Statistical fine-mapping of a GWAS locus with Wakefield approximate Bayes factors (genor finemap abf) to get per-variant posterior probabilities and a 95% credible set. Use when the user asks which variant in a locus is likely causal, wants a credible set, PIPs or posterior probabilities, or mentions fine-mapping, SuSiE, FINEMAP or COJO.
---
# Fine-mapping (Wakefield ABF, single causal variant)

GeNor fine-maps one region at a time with Wakefield's approximate Bayes factor (ABF), a port of
`coloc::finemap.abf`. It uses only the marginal summary statistics in the region and assumes
**exactly one causal variant (or none) in the region**. It needs no LD panel, which is both its
convenience and its weakness.

Not in GeNor yet: SuSiE / SuSiE-RSS (multiple causal variants with an LD matrix), GCTA-COJO
conditional and joint analysis, FINEMAP, PolyFun, functional priors. They are planned for a later
milestone. (`genor backend install gcta` downloads GCTA, but no GeNor command drives COJO yet;
the `r` backend lists `r.susier` but nothing is wired to it.) If the user needs multi-signal
fine-mapping now, say so and point them to running SuSiE-RSS or COJO outside GeNor.

## When to use

- After `genor clump run` gives index variants and the user wants to narrow a locus to a set of
  candidate variants.
- Before `colocalisation` or `gene-pathway-celltype`, to know which variants carry the signal.
- To check whether a locus has a sharp signal (small credible set) or a diffuse one.

## Inputs and prerequisites

- A loaded trait directory from `genor sumstats load` (see `sumstats-qc`).
- Best: `beta` and `standard_error` for every variant in the region. Then, for quantitative
  traits, `effect_allele_frequency` and `n` let GeNor estimate sdY.
- Fallback without beta/SE: `p_value`, `effect_allele_frequency` and `n`; for case-control traits
  the case fraction too (load with `--n-cases` / `--n-controls`).
- The region in the trait's own genome build. Check `build` in the trait's `sumstats.json` or the
  load result.
- No backend or reference data needed.

## Steps

1. Get candidate loci (needs PLINK 2 and an LD panel for `--mode plink`):
   ```bash
   genor clump run runs/bmi --preset loci --out runs/clump_bmi --json
   ```
   `index_variants.tsv` lists the index variants; `summary.top` shows the 15 strongest.
2. **Ask the user** which loci to fine-map (do not fine-map every locus silently; see Decision
   points).
3. Check for multiple signals in the window before trusting a single-causal result. Count index
   variants from step 1 that fall inside the planned window. More than one means the
   single-causal assumption is likely broken. See
   [references/multiple-signals.md](references/multiple-signals.md).
4. Run ABF fine-mapping around the lead variant, or on an explicit region:
   ```bash
   genor finemap abf runs/bmi --lead rs796513 --window-kb 500 --out runs/fm_bmi_rs796513 --json
   genor finemap abf runs/bmi --region 12:69565626-70565626 --out runs/fm_bmi_chr12 --json
   ```
   MCP: `run_analysis(method="finemap.abf", params={"trait": "runs/bmi", "lead": "rs796513",
   "window_kb": 500, "p1": 1e-4, "level": 0.95, "out": "runs/fm_bmi_rs796513"})`. If it returns a
   `job_id`, poll `job_status` then `job_result`.
5. Read `pp_null` first, then the credible set (see Reading the output).
6. Optional sensitivity: rerun with a narrower window (e.g. `--window-kb 250`) and with `--p1 1e-5`.
   A stable top variant and credible set is reassuring; a set that changes a lot is fragile.
7. Hand off: candidate variants to `gene-pathway-celltype` (VEP, genes in region) and to
   `open-targets-triage` (`genor lookup variant` for existing credible sets and L2G).

## Defaults and why

| Option | Default | Why |
|---|---|---|
| `--window-kb` | 500 (lead +/- 500 kb) | Covers most LD blocks in European-ancestry data. Wider windows add variants and raise the chance of a second signal. |
| `--p1` | 1e-4 | Prior probability that any given variant is causal (Giambartolomei et al. 2014; coloc default). |
| `--level` | 0.95 | Conventional credible set coverage. |
| Prior SD of effect | 0.15 x sdY (quantitative); 0.2 (case-control, log-OR) | Wakefield 2009; coloc defaults. Encodes the belief that true per-allele effects are small. |
| sdY | Estimated from SE, MAF and N (coloc `sdY.est`); 1 if not estimable | sdY scales the prior to the trait's units. |

Citations recorded in the manifest: Wakefield 2009 (`wakefield-abf`), Giambartolomei 2014 (`coloc`).
Method details and formulas: [references/abf-method.md](references/abf-method.md).

## Decision points

- **Ask the user which loci to follow up.** Offer the clumped index variants with p-values, and
  suggest prioritising genome-wide significant loci (p < 5e-8). Fine-mapping a sub-threshold
  locus mostly returns `finemap_no_signal`.
- **Ask the user whether the trait is standardised** (e.g. inverse-normal transformed) when
  `sdy_assumed` appears or the beta scale is unclear. Wrong sdY shifts the prior.
- **Ask the user what to do with multi-signal loci**: narrow the window around each index variant
  (approximate), or run SuSiE-RSS / COJO outside GeNor. Do not decide for them.
- **Ask the user before relying on an imputation-sparse file** (e.g. HapMap3-only variants): the
  causal variant may simply be absent.

## Reading the output

The JSON `summary` has:

- `region`, `n_variants`: what was fine-mapped. A region with very few variants gives coarse results.
- `pp_null`: posterior that no variant in the region is causal. **Read this first.**
  - Near 0: there is a signal; the credible set is meaningful.
  - Above 0.5: GeNor adds the `finemap_no_signal` warning.
  - The credible set is conditional on one causal variant being present (PPs renormalised
    without the null). When `pp_null` is high, the set exists but describes a signal that is
    probably not there: report "no credible signal", not the set.
- `credible_set`, `credible_set_size`: smallest set of variants whose PPs sum to >= `level`.
  The JSON lists at most 50; the full ranking is in `finemap.tsv`.
- `top_variant`, `top_pp`: the highest-PP variant and its posterior.
- `type` (`quant` or `cc`), `sd_prior`, `sdY_estimated`: the prior actually used.

`finemap.tsv` columns: `variant_id`, `chromosome`, `base_pair_location`, `p_value`, `lABF`
(log ABF), `PP`, sorted by PP.

Rough reading guide (assuming one signal and good coverage):

- Credible set of 1-5 variants with `top_pp` > 0.5: well resolved.
- Tens of variants with flat PPs: the signal sits in a long LD block; fine-mapping did not resolve it.
- Top PP close to 1 often just reflects a very strong signal; it still assumes the causal variant
  was genotyped or imputed.

## Failure signatures

| Code | Meaning | Fix |
|---|---|---|
| `bad_param` | Neither `--region` nor `--lead` given | Pass one of them. |
| `bad_region` | Region string not `chr:start-end` | e.g. `12:69565626-70565626` (a `chr` prefix is accepted). |
| `not_found` | `--lead` ID not in the trait | Use the trait's `variant_id` exactly (often rsID); check clumping output. |
| `region_empty` | Fewer than 2 variants in the region | Wrong chromosome naming or build; widen the window. |
| `coloc_inputs` | No beta/SE, and no MAF + N for the p-value route; or case-control without case fraction | Reload with `--n`, or `--n-cases` and `--n-controls`. |
| `not_a_trait_dir` | Path is not a loaded trait | Point to the `sumstats load` output directory. |
| warning `finemap_no_signal` | `pp_null` > 0.5 | No clear signal; do not interpret the credible set. |
| warning `sdy_assumed` | sdY not estimable (no MAF or N) | Prior assumes a standardised trait; ask the user. |
| info `single_causal_assumption` | Always emitted | Confirm a single signal (step 3) and say so in the report. |

Prior sanity: GeNor's null prior is `1 - n_variants x p1`. Keep `n_variants x p1` well below 1
(for example under 0.5). With very dense imputed data (thousands of variants per Mb) use a smaller
window or smaller `--p1`; GeNor does not check this for you.

## What not to conclude

- The top variant is not proven causal. PP is conditional on the single-causal model, the prior,
  and the causal variant being present in the file.
- A credible set is not a functional annotation. Pair it with VEP and QTL evidence before naming a
  mechanism.
- A large credible set does not mean "many causal variants"; it means the data cannot distinguish
  variants in LD.
- Results from one ancestry do not transfer to another: LD differs, so credible sets differ.
- Do not compare PPs across loci or traits as if they were effect sizes.

## Related skills

- `genor` (router), `sumstats-qc` (loading, build), `colocalisation` (shared causal variants),
  `gene-pathway-celltype` (variants to genes), `open-targets-triage` (existing credible sets, L2G),
  `interpretation-and-reporting`, `gwas-concepts` (LD, priors, Bayes factors explained).
- References: [abf-method.md](references/abf-method.md),
  [multiple-signals.md](references/multiple-signals.md).
