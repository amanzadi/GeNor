# Design pitfalls

## Sample overlap

Two GWAS sharing participants (e.g. both include UK Biobank) have correlated sampling errors.

- **Genetic correlation (LDSC rg)**: handled by the cross-trait intercept, which absorbs overlap.
  A non-zero intercept is evidence of overlap or shared confounding (GeNor warning
  `sample_overlap_signal`).
- **Two-sample MR**: overlap combined with weak instruments biases estimates towards the
  observational association. GeNor asks for the overlap (`--overlap none|partial|substantial|unknown`)
  and records it; the user must state it, not the agent.
- **Meta-analysis**: overlapping samples counted twice give falsely small SEs.
- **condFDR / conjFDR**: overlap inflates cross-trait enrichment; random pruning and the
  pleioFDR procedure handle some of it, but large overlap should be reported.
- **Colocalisation**: overlap mostly matters for power and calibration, less for bias.

## Collider and selection bias

Conditioning on a variable affected by both the exposure and the outcome (a collider) creates
spurious associations.

- **Covariate adjustment**: GWAS of trait A adjusted for heritable trait B (e.g. waist adjusted for
  BMI) can produce associations driven by B. Using such GWAS in MR can bias results.
- **Selection**: cohorts that select on health, education, or disease status (e.g. UK Biobank's
  healthy volunteer bias, hospital-based cases) can induce associations between unrelated traits.
- **Index-event bias**: studying disease progression only among people who already have the
  disease.

## Ascertainment and case-control design

- Case-control studies over-sample cases (sample prevalence, e.g. 30%) relative to the population
  prevalence (e.g. 15% for depression). Heritability on the observed 0-1 scale depends on the
  sample prevalence; conversion to the liability scale needs the population prevalence
  (`--pop-prev` in `genor ldsc h2`, or `--population-prevalence` on load). GeNor asks for it.
- **Proxy phenotypes** ("by-proxy" cases from family history, self-report, broad definitions
  like "depression" from a single question) increase N but change what is measured; genetic
  correlation with the strict phenotype may be below 1.
- **Control contamination**: undiagnosed cases among controls attenuate effects.

## Cryptic relatedness

Unreported relatives in a sample violate independence and inflate test statistics. Mixed models
handle it; in summary statistics it shows as an elevated LDSC intercept.

## Population stratification

See population-genetics.md. Also relevant for MR: within-family GWAS estimates can differ from
population estimates for traits influenced by assortative mating, dynastic effects (parental
genotype influencing offspring environment), or stratification (e.g. educational attainment).

## Heterogeneity in meta-analysis

Different cohorts may define the phenotype differently, use different ancestries, or have
different environments. Cochran's Q and I2 quantify heterogeneity. Fixed-effect estimates are
invalid under substantial heterogeneity. See `meta-analysis`.

## Explaining to a user

"Both studies include UK Biobank. That shared sample makes their errors correlated. For MR this
can pull the estimate towards the ordinary observational correlation, so please tell me how much
the samples overlap."
