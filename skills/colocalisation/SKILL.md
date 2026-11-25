---
name: colocalisation
description: Bayesian colocalisation of two association signals in one region with coloc ABF (genor coloc abf), giving posterior probabilities PP.H0-PP.H4 for no signal, one trait only, distinct causal variants, or a shared causal variant. Use when the user asks whether two GWAS, or a GWAS and an eQTL/pQTL/sQTL, share a causal variant at a locus, or mentions coloc, PP.H4, shared signal, or QTL colocalisation.
---
# Colocalisation (coloc ABF)

`genor coloc abf` is a port of `coloc::coloc.abf` (coloc 5.x). For one region it compares five
hypotheses about two traits, assuming **at most one causal variant per trait in the region**:

| Hypothesis | Meaning |
|---|---|
| H0 | No causal variant for either trait |
| H1 | Causal variant for trait 1 only |
| H2 | Causal variant for trait 2 only |
| H3 | Both traits have a causal variant, and they are different variants |
| H4 | Both traits share one causal variant |

GeNor does not query eQTL or pQTL APIs. QTL data must be downloaded by the user and loaded as a
trait. Open Targets' precomputed colocalisations are available with `genor lookup coloc`.

## When to use

- A GWAS locus and a QTL for a nearby gene in a relevant tissue: does the GWAS signal act through
  that gene's expression (or splicing, or protein level)?
- Two GWAS traits at a shared locus (e.g. after `pleiotropy-condfdr` conjFDR finds a shared locus).
- Not for genome-wide questions: use `heritability-rg` for genome-wide genetic correlation.

## Inputs and prerequisites

- Two loaded trait directories (`genor sumstats load`, see `sumstats-qc`), **on the same genome
  build**. GeNor stops with `build_mismatch` if both builds are known and differ.
- For each trait: beta + SE (preferred), or p-value + allele frequency + N (plus case fraction for
  case-control). Correct `trait_type` (quantitative vs binary) matters for the prior.
- For a QTL: one molecular trait (one gene, or one gene-tissue pair) loaded as its own trait with
  its sample size. See [references/qtl-data.md](references/qtl-data.md).
- Enough shared variants in the region (at least 10; ideally hundreds). Dense, overlapping
  coverage matters more than anything else.

## Steps

1. Choose the region: a GWAS lead variant (from `genor clump run`) +/- window, or an explicit
   region. The `--lead` variant is looked up in **trait 1**, so put the trait that contains it first.
2. **Ask the user** which QTL dataset and tissue to use, and which loci (Decision points).
3. If the QTL is on another build, lift one trait: `genor sumstats liftover runs/eqtl_GENE --to GRCh37 --out runs/eqtl_GENE_b37`.
4. Run coloc:
   ```bash
   genor coloc abf runs/bmi runs/mdd --lead rs2237324 --window-kb 500 --out runs/coloc_bmi_mdd_rs2237324 --json
   genor coloc abf runs/bmi runs/eqtl_GENE --region 12:69565626-70565626 --out runs/coloc_bmi_GENE --json
   ```
   MCP: `run_analysis(method="coloc.abf", params={"trait1": "runs/bmi", "trait2": "runs/mdd",
   "lead": "rs2237324", "window_kb": 500, "p1": 1e-4, "p2": 1e-4, "p12": 1e-5, "out": "..."})`.
5. Prior sensitivity on `p12` (GeNor has no built-in sweep; run it again):
   ```bash
   genor coloc abf runs/bmi runs/mdd --lead rs2237324 --p12 5e-6 --out runs/coloc_bmi_mdd_rs2237324_p12_5e-6 --json
   genor coloc abf runs/bmi runs/mdd --lead rs2237324 --p12 1e-6 --out runs/coloc_bmi_mdd_rs2237324_p12_1e-6 --json
   ```
6. Check the single-causal-variant assumption for both traits (see `fine-mapping`,
   [../fine-mapping/references/multiple-signals.md](../fine-mapping/references/multiple-signals.md)).
   Run `genor finemap abf` on each trait in the same region and compare top variants.
7. Optional cross-check with public data: `genor lookup variant <GRCh38 id>` to get a
   `studyLocusId`, then `genor lookup coloc <studyLocusId>` (see `open-targets-triage`).
8. Report with `genor report build` (see `interpretation-and-reporting`).

## Defaults and why

| Option | Default | Why |
|---|---|---|
| `--p1` | 1e-4 | Prior that a variant is causal for trait 1 only (Giambartolomei et al. 2014). |
| `--p2` | 1e-4 | Same for trait 2. |
| `--p12` | 1e-5 | Prior that a variant is causal for both. coloc's long-standing default. Wallace (2020) argues 1e-5 is permissive and suggests 5e-6 as a more conservative choice. |
| `--window-kb` | 500 | Lead +/- 500 kb. |
| Prior SD of effect | 0.15 x sdY quantitative; 0.2 case-control (log-OR) | As in coloc and `fine-mapping`. |

Citation recorded: Giambartolomei 2014 (`coloc`). Details, formulas, and prior reasoning:
[references/coloc-method.md](references/coloc-method.md).

## Decision points

- **Ask the user which QTL dataset and tissue** fit the biology (e.g. brain for psychiatric
  traits, adipose or hypothalamus for BMI, blood for immune traits). Show them
  `genor sumstats sources gtex`, `genor sumstats sources eqtl-catalogue`, `genor sumstats sources eqtlgen`.
  GeNor never downloads them.
- **Ask the user which loci and genes** to test. Testing every gene in every window inflates
  the number of chances to find a high PP.H4 by luck; record how many tests were run.
- **Ask the user about the p12 prior** if PP.H4 changes conclusion between 1e-5 and 5e-6.
  Report both.
- **Ask the user what to do with multi-signal loci** (coloc-SuSiE is not in GeNor).
- **Ask the user to confirm the QTL sample ancestry** matches the GWAS.

## Reading the output

`summary` fields:

- `nsnps`: shared variants after harmonisation. Low counts (tens) make results fragile.
- `PP.H0` ... `PP.H4`: sum to 1.
- `interpretation`: GeNor's rule of thumb: PP.H4 > 0.8 "shared causal variant"; PP.H3 > 0.8
  "distinct causal variants"; PP.H0 + PP.H1 + PP.H2 > 0.5 "at least one trait lacks a clear signal";
  else "Inconclusive".
- `top_shared_variant`: variant with the highest `SNP.PP.H4`. Meaningful only when PP.H4 is high.
- `priors`, `trait1` / `trait2` (names), and `prior_trait1` / `prior_trait2` blocks with `type`,
  `sd_prior`, `sdY_estimated`.

`coloc_snp_pp.tsv`: `variant_id`, `SNP.PP.H4` (posterior that each variant is the shared causal
variant, given H4), sorted.

How to read it:

1. Is there power? If PP.H1 or PP.H2 dominates, one trait has no signal in the region (common for
   QTLs of the wrong gene or tissue, and for weak GWAS loci). This is "no evidence", not
   "evidence against".
2. Given both traits have signals, compare PP.H4 with PP.H3. PP.H4 / (PP.H3 + PP.H4) is useful
   when PP.H1 or PP.H2 absorbs much of the mass.
3. PP.H4 > 0.8 is the common convention for "likely shared". It is a convention, not a calibrated
   threshold.
4. Check `top_shared_variant` against the GWAS lead and the fine-mapping credible set.

Example from GeNor's simulated data: at rs796513 (BMI_sim lead), PP.H1 = 0.997, so MDD_sim has
no signal there. At rs2237324 (MDD_sim lead), PP.H4 = 1.0.

## Failure signatures

| Code | Meaning | Fix |
|---|---|---|
| `build_mismatch` | Traits on different builds | `genor sumstats liftover <trait> --to <build>`. |
| `coloc_too_few` | Fewer than 10 shared variants | Different builds or ID types, a QTL file filtered to significant pairs only, or wrong region. |
| `cannot_match_variants` | No shared rsIDs and no chr/pos in both | Reload with position columns, or use consistent IDs. |
| `bad_param` / `bad_region` / `not_found` | Missing region/lead, malformed region, lead not in trait 1 | Pass `--region chr:start-end`, or put the trait with the lead first. |
| `coloc_inputs` | No beta/SE and no MAF+N; or case-control without case fraction | Reload with `--n` / `--n-cases` / `--n-controls`. |
| warning `allele_mismatch` (serious) | > 5% of shared variants have incompatible alleles | Likely build or ID mismatch; stop and check. |
| warning `palindromic_dropped_no_freq` | A/T and C/G variants dropped (no frequencies) | Load allele frequencies if available. |
| warning `sdy_assumed` | sdY not estimable | Ask whether the trait is standardised. |
| info `single_causal_assumption` | Always emitted | Check for multiple signals (step 6). |

Note: if either trait's build is unknown, GeNor cannot check builds. Confirm both builds before
running.

## What not to conclude

- High PP.H4 does not show that the gene causes the disease. It shows the two signals are
  consistent with one shared variant. LD with an unobserved variant, or a variant acting on two
  nearby genes, gives the same result.
- PP.H3 high does not prove the signals are distinct if either trait has several signals.
- Low PP.H4 from an underpowered QTL is not evidence against the gene.
- A PP.H4 for one tissue does not exclude other tissues; QTLs are widely shared across tissues.
- Do not treat PPs as p-values or correct them with Bonferroni; report the number of tests.

## Related skills

- `fine-mapping` (single-causal check, credible sets), `open-targets-triage` (`genor lookup coloc`),
  `gene-pathway-celltype` (genes in the region), `mendelian-randomisation` (causal direction is a
  different question), `pleiotropy-condfdr` (finding shared loci genome-wide),
  `interpretation-and-reporting`, `gwas-concepts`.
- References: [coloc-method.md](references/coloc-method.md), [qtl-data.md](references/qtl-data.md).
