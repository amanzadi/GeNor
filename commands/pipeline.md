---
description: Run the standard post-GWAS pipeline on the user's summary statistics
argument-hint: "<sumstats files and what you want to know>"
---
Use the `genor` skill, then run the standard pipeline for: $ARGUMENTS

1. `genor doctor --json`.
2. Load and QC each file (`sumstats-qc` skill). Ask the user for ancestry, sample size or case/control
   counts, and the build if it is not inferred. Never guess these.
3. LDSC h2 for each trait and rg between traits (`heritability-rg`), if the reference data is
   available or the user agrees to fetch it.
4. Clump into loci. Then do only what the user's question needs (MR, condFDR, fine-mapping, coloc,
   gene mapping).
5. Build a report with `genor report build` and summarise: key numbers with uncertainty, every
   serious warning, and the decisions that were left to the user.

Stop and ask whenever a result contains `ask_user` or the exit code is 3.
