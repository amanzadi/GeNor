---
description: Two-sample Mendelian randomisation with sensitivity analyses and a STROBE-MR checklist
argument-hint: "<exposure sumstats> <outcome sumstats> [notes]"
---
Use the `mendelian-randomisation` skill for: $ARGUMENTS

Before running, confirm with the user: the exposure and outcome definitions, the ancestry of both
GWAS, and whether the samples overlap. Load both files if they are not loaded yet, run
`genor mr run --exposure <dir> --outcome <dir> --overlap <answer> --json`, then report the IVW
estimate (OR for binary outcomes) with its 95% CI, the sensitivity estimators, heterogeneity,
Egger intercept, Steiger direction and F statistics, and point to `strobe_mr.md` for the items the
authors must complete.
