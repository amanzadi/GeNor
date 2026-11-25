# Checking for multiple independent signals

ABF fine-mapping and coloc ABF both assume at most one causal variant per trait in the region.
Many GWAS loci have two or more independent signals, especially in large studies. If the
assumption fails, ABF tends to put most posterior on the strongest signal and ignore the others;
the credible set is then too confident and incomplete.

## What GeNor can do now

### 1. LD clumping at r2 0.1 inside the window

```bash
genor clump run runs/bmi --preset loci --out runs/clump_bmi --json
```

`--preset loci` uses r2 0.1 and a 250 kb window, with PLINK 2 against a reference panel
(`--panel`, e.g. `g1000-eur`; fetch with `genor refdata fetch`, install PLINK with
`genor backend install plink2`). Count the index variants in `index_variants.tsv` that fall within
your fine-mapping window.

- One index variant: consistent with a single signal (not proof).
- Two or more genome-wide significant index variants within the window: suspect multiple signals.

Limits:
- Clumping is greedy and marginal. Two variants with r2 < 0.1 can both be marginally significant
  because of a third variant, and a true secondary signal masked by the primary can be missed.
  Conditional analysis is the proper test.
- `--mode distance` uses no LD and cannot answer this question.
- The reference panel must match the GWAS ancestry. GeNor's panels (`g1000-eur`, `g1000-eas`) are
  GRCh37 and are checked against the trait's build and ancestry; a GRCh38 trait needs
  `genor sumstats liftover <trait> --to GRCh37` first.

### 2. Narrow windows around each index variant (approximation)

If the index variants are far enough apart (hundreds of kb, low LD between them), fine-map each
with a window that contains only one of them:

```bash
genor finemap abf runs/bmi --lead rsA --window-kb 100 --out runs/fm_rsA --json
genor finemap abf runs/bmi --lead rsB --window-kb 100 --out runs/fm_rsB --json
```

This is approximate: marginal statistics near rsA still contain some signal from rsB through LD.
Report it as such. Ask the user before taking this route.

### 3. Look at existing fine-mapping

Open Targets runs SuSiE-based fine-mapping on many public GWAS. For a GRCh38 variant ID,
`genor lookup variant <CHR_POS_REF_ALT>` lists credible sets containing that variant across studies
(see `open-targets-triage`). Several credible sets for the same study near the locus indicate
multiple signals in that study.

## Not in GeNor yet

- **GCTA-COJO** conditional and joint analysis (Yang et al. 2012). The GCTA binary can be
  installed with `genor backend install gcta`, and the citation key `gcta-cojo` exists, but there
  is no `genor cojo` command. Planned for a later milestone.
- **SuSiE / SuSiE-RSS** (Wang et al. 2020; Zou et al. 2022): multiple causal variants with an LD
  matrix. Listed as an optional R backend (`r.susier`) in the registry but not wired. Planned.
- **coloc-SuSiE** (Wallace 2021) for colocalisation with multiple signals. Not in GeNor.

If the user needs these, say they are not in GeNor and that the usual route is to run them
outside GeNor with an LD matrix from in-sample genotypes or a well-matched reference panel.
LD mismatch between the GWAS and the reference panel is the main source of false signals in
SuSiE-RSS and COJO; in-sample LD is preferred.

## How to report

State in the methods: "Fine-mapping used Wakefield ABF assuming a single causal variant per
region. We checked for additional signals by LD clumping (r2 < 0.1, 250 kb, <panel>); loci with
more than one index variant in the window were <excluded | fine-mapped in narrowed windows>."
