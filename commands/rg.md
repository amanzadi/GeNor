---
description: SNP heritability and genetic correlation between two traits with LDSC
argument-hint: "<trait A sumstats or run dir> <trait B sumstats or run dir>"
---
Use the `heritability-rg` skill for: $ARGUMENTS

Load the traits if needed, check `genor doctor` for the LDSC reference data (ask before fetching),
run `genor ldsc rg <A> <B> --json`, and report h2 for each trait (liability scale for binary traits
when the population prevalence is known, so ask for it), rg with SE and p, and the cross-trait
intercept as evidence on sample overlap.
