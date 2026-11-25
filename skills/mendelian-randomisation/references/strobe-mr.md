# STROBE-MR checklist and GeNor's coverage

Skrivankova VW, et al. STROBE-MR Statement. JAMA 2021;326:1614-1621. GeNor writes
`strobe_mr.md` and `strobe_mr.json` in every MR run directory. Status values:

- **filled**: GeNor writes it from the run manifest and results. Check it, do not rewrite it.
- **partial**: GeNor supplies what it knows; the text says "HUMAN:" for the rest.
- **human**: only the authors can write it. Ask the user; never invent rationale, funding or conflicts.

| # | Item | What it asks | GeNor status |
|---|---|---|---|
| 1 | Title and abstract | Indicate MR as the design in title and/or abstract | partial (suggested title; authors write the abstract) |
| 2 | Background | Scientific background and rationale | human |
| 3 | Objectives | Objectives, including pre-specified causal hypotheses | human |
| 4 | Study design and data sources | Design; data sources (setting, participants, exposures, outcomes, variants, covariates) | partial (traits, N, cases/controls, ancestry, build, source paths; HUMAN: settings, recruitment, phenotype definitions, sample overlap) |
| 5 | Assumptions | The three core IV assumptions and any additional ones | partial (states relevance, independence, exclusion restriction and how they were probed; HUMAN: extra assumptions such as homogeneity or monotonicity) |
| 6 | Statistical methods: main analysis | Variable handling, variant selection, MR estimator | filled (threshold, clumping mode, r2, window, panel, harmonisation, IVW random effects, scale) |
| 7 | Assessment of assumptions | Methods used to assess assumptions | filled (F, Cochran's Q, Rucker's Q', Egger intercept, Steiger, LDSC intercept if run) |
| 8 | Sensitivity and additional analyses | Sensitivity analyses performed | filled (Egger, weighted median, weighted/simple mode, IVW fixed, leave-one-out, Wald ratios) |
| 9 | Software and pre-registration | Software, versions, settings; pre-registration | partial (GeNor version, seed, bootstraps, PLINK 2; HUMAN: pre-registration) |
| 10 | Descriptive data | Numbers of individuals and variants; descriptive statistics | filled (candidates, after clumping, after harmonisation, palindromic and mismatch drops, missing in outcome) |
| 11 | Main results | Main estimates with precision; scale | filled (IVW beta, 95% CI, p, number of SNPs; OR note for binary outcomes) |
| 12 | Assessment of assumptions (results) | Heterogeneity, F statistics, I2_GX, etc. | filled (F mean/min, weak count, Q, Egger intercept, I2_GX) |
| 13 | Sensitivity analyses (results) | Robustness to assumption violations, directionality | filled (all sensitivity estimates, Steiger; leave-one-out files) |
| 14 | Key results | Summary with reference to objectives | human |
| 15 | Limitations | Validity of IV assumptions, other biases, imprecision | partial (lists machine-detected warnings; HUMAN: pleiotropy, stratification, selection, overlap, power) |
| 16 | Interpretation | Triangulation, clinical relevance, mechanism | human |
| 17 | Generalizability | Other populations, exposure timing and levels | partial (ancestry of both GWAS; HUMAN: timing and levels) |
| 18 | Funding | Funding sources and roles | human |
| 19 | Data and data sharing | Data access, source references | partial (source files; manifest with SHA-256; HUMAN: URLs, accessions, data-use terms) |
| 20 | Conflicts of interest | Declarations | human |

Totals per run: 7 filled, 7 partial, 6 human (`summary.strobe_mr`).

What to ask the user to complete: items 2, 3, 14, 16, 18, 20 in full; the HUMAN parts of 1, 4,
5, 9, 15, 17, 19. Item 4's sample-overlap description should match the `--overlap` answer.

Note: GeNor runs only the Steiger directionality test, not SNP-level Steiger filtering.
