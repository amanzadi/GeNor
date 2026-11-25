"""STROBE-MR checklist (Skrivankova et al., JAMA 2021) filled from an MR run.

Each of the 20 items gets a status:
- ``filled``: GeNor knows this from the run manifest and results.
- ``partial``: GeNor supplies the parts it knows; a human must complete the rest.
- ``human``: only the authors can write this (rationale, funding, interpretation).
"""

from __future__ import annotations

from typing import Any

ITEMS = [
    (1, "Title and abstract", "Indicate Mendelian randomization (MR) as the study design in the title and/or abstract."),
    (2, "Background", "Explain the scientific background and rationale for the reported study."),
    (3, "Objectives", "State specific objectives clearly, including pre-specified causal hypotheses."),
    (4, "Study design and data sources", "Present key elements of the study design; describe the data sources (setting, participants, exposures, outcomes, genetic variants, covariates)."),
    (5, "Assumptions", "Explicitly state the three core IV assumptions for the main analysis (relevance, independence and exclusion restriction) and any additional assumptions."),
    (6, "Statistical methods: main analysis", "Describe statistical methods and statistics used: how quantitative variables were handled, how genetic variants were selected, the MR estimator and related statistics."),
    (7, "Assessment of assumptions", "Describe any methods or prior knowledge used to assess the assumptions or justify their validity."),
    (8, "Sensitivity analyses and additional analyses", "Describe any sensitivity analyses or additional analyses performed."),
    (9, "Software and pre-registration", "Name statistical software and package(s), including version and settings used; state whether the study protocol and details were pre-registered."),
    (10, "Descriptive data", "Report the numbers of individuals and genetic variants, and descriptive statistics of the data sources."),
    (11, "Main results", "Report the results of the main analysis with estimates and precision (e.g. 95% CI); report the scale of the estimates."),
    (12, "Assessment of assumptions", "Assess the validity of the assumptions; report any additional statistics (e.g. assessments of heterogeneity across genetic variants, F statistics, I^2_GX)."),
    (13, "Sensitivity analyses and additional analyses", "Report any sensitivity analyses to assess the robustness of the main results to violations of the assumptions; report any additional analyses (e.g. directionality); report any assessment of direction of effect."),
    (14, "Key results", "Summarize key results with reference to study objectives."),
    (15, "Limitations", "Discuss limitations of the study, taking into account the validity of the IV assumptions, other sources of potential bias, and imprecision."),
    (16, "Interpretation", "Meta-synthesis of results and triangulation; clinical relevance; mechanism."),
    (17, "Generalizability", "Discuss the generalizability of the study results to other populations, different exposure periods/timings, and different levels of exposure."),
    (18, "Funding", "Describe sources of funding and the role of funders in the present study."),
    (19, "Data and data sharing", "Provide the data used to perform all analyses or report where and how the data can be accessed; reference the source data."),
    (20, "Conflicts of interest", "All authors should declare all potential conflicts of interest."),
]


def _fmt(x: Any, d: int = 3) -> str:
    try:
        return f"{float(x):.{d}g}"
    except (TypeError, ValueError):
        return str(x)


def build(manifest: dict[str, Any], results: list[dict[str, Any]]) -> list[dict[str, Any]]:
    s = manifest.get("summary", {})
    p = manifest.get("params", {})
    meta = manifest.get("meta", {})
    ex, out = meta.get("exposure", {}), meta.get("outcome", {})
    main = next((r for r in results if r["method"] == "IVW"), results[0] if results else {})
    sens = [r for r in results if r["method"] != main.get("method")]
    filled: dict[int, tuple[str, str]] = {}

    filled[1] = ("partial", f"Suggested title: 'Effect of {ex.get('trait')} on {out.get('trait')}: "
                 "a two-sample Mendelian randomization study'. Authors write the abstract.")
    filled[4] = ("partial",
                 f"Two-sample MR with summary statistics. Exposure: {ex.get('trait')} "
                 f"(N={ex.get('n')}, ancestry {ex.get('ancestry')}, build {ex.get('build')}, source "
                 f"{ex.get('source')}). Outcome: {out.get('trait')} (N={out.get('n')}, cases="
                 f"{out.get('n_cases')}, controls={out.get('n_controls')}, ancestry {out.get('ancestry')}, "
                 f"source {out.get('source')}). HUMAN: describe settings, recruitment, phenotype "
                 "definitions, and sample overlap between the two GWAS.")
    filled[5] = ("partial", "Relevance: variants associated with the exposure at the selection "
                 "threshold. Independence: no confounding of variant-outcome association (assumed; "
                 "same-ancestry GWAS, PC-adjusted). Exclusion restriction: variants affect the outcome "
                 "only via the exposure (probed by Egger intercept, heterogeneity, robust estimators). "
                 "HUMAN: state any additional assumptions (e.g. homogeneity/monotonicity for "
                 "point estimates).")
    filled[6] = ("filled",
                 f"Instruments: exposure p < {p.get('p_threshold')}, clumped with {p.get('clump_mode')} "
                 f"(r2 < {s.get('clump', {}).get('r2')}, window {s.get('clump', {}).get('kb')} kb, "
                 f"panel {s.get('clump', {}).get('panel', 'none')}). Alleles harmonised to the exposure "
                 "effect allele; palindromic SNPs with MAF > 0.42 removed. Main estimator: "
                 "multiplicative random-effects IVW. Estimates are per unit of the exposure's GWAS "
                 "scale on the outcome's scale (log-odds for binary outcomes).")
    filled[7] = ("filled", "Instrument strength: per-SNP F statistics. Heterogeneity: Cochran's Q "
                 "(IVW) and Rucker's Q' (Egger). Directional pleiotropy: MR-Egger intercept. "
                 "Direction: Steiger directionality test (no SNPs filtered). Overlap: cross-trait LDSC intercept if run.")
    filled[8] = ("filled", "MR-Egger, weighted median, weighted mode, simple mode, IVW fixed effects, "
                 "leave-one-out IVW, single-SNP Wald ratios.")
    filled[9] = ("partial", f"GeNor {manifest.get('genor_version')} (Python), estimators validated "
                 f"against TwoSampleMR; PLINK 2 for clumping where used; seed {p.get('seed')}, "
                 f"{p.get('nboot')} bootstrap replicates. HUMAN: state pre-registration.")
    f = s.get("F", {})
    hz = s.get("harmonisation", {})
    filled[10] = ("filled", f"{s.get('n_candidates')} exposure variants passed the threshold; "
                  f"{s.get('n_instruments_clumped')} after clumping; {s.get('nsnp')} after "
                  f"harmonisation with the outcome (palindromic dropped: {hz.get('palindromic_dropped')}, "
                  f"allele mismatches: {hz.get('allele_mismatch_dropped')}, missing in outcome: "
                  f"{s.get('n_missing_in_outcome')}).")
    if main:
        lo = main["b"] - 1.96 * main["se"]
        hi = main["b"] + 1.96 * main["se"]
        filled[11] = ("filled", f"{main['method']}: beta = {_fmt(main['b'])} (95% CI {_fmt(lo)} to "
                      f"{_fmt(hi)}), p = {_fmt(main['pval'])}, {main.get('nsnp')} SNPs. For a binary "
                      "outcome, OR per unit = exp(beta).")
    het = s.get("heterogeneity_ivw", {})
    ple = s.get("pleiotropy_egger", {})
    st = s.get("steiger", {})
    filled[12] = ("filled", f"F mean {_fmt(f.get('F_mean'))}, min {_fmt(f.get('F_min'))} "
                  f"({f.get('n_weak_F_lt_10')} SNPs with F < 10). Cochran's Q = {_fmt(het.get('Q'))} "
                  f"(df {het.get('Q_df')}, p = {_fmt(het.get('Q_pval'))}). Egger intercept "
                  f"{_fmt(ple.get('intercept'))} (p = {_fmt(ple.get('pval'))}). I2_GX = "
                  f"{_fmt(s.get('I2_GX'))}.")
    filled[13] = ("filled", "; ".join(f"{r['method']}: {_fmt(r['b'])} (SE {_fmt(r['se'])}, p {_fmt(r['pval'])})"
                                      for r in sens)
                  + (f". Steiger: exposure explains more variance than outcome = "
                     f"{st.get('correct_causal_direction')} (p = {_fmt(st.get('steiger_pval'))})." if st else "")
                  + " Leave-one-out and single-SNP results are in the run directory.")
    filled[15] = ("partial", "Machine-detected issues: " + ("; ".join(
        w["message"] for w in manifest.get("warnings", []) if w.get("severity") in ("serious", "warning"))
        or "none") + ". HUMAN: discuss residual pleiotropy, population stratification, "
        "selection bias, sample overlap, and power.")
    filled[17] = ("partial", f"Both GWAS are {ex.get('ancestry')}/{out.get('ancestry')} ancestry; "
                  "results may not transfer to other ancestries. HUMAN: discuss timing and levels of exposure.")
    filled[19] = ("partial", "Summary statistics sources: " + f"exposure {ex.get('source')}; outcome "
                  f"{out.get('source')}. GeNor run manifest (inputs with SHA-256) accompanies results. "
                  "HUMAN: give access URLs / accession numbers and data-use terms.")
    items = []
    for num, name, desc in ITEMS:
        status, text = filled.get(num, ("human", "To be completed by the authors."))
        items.append({"item": num, "name": name, "recommendation": desc, "status": status, "text": text})
    return items


def to_markdown(items: list[dict[str, Any]]) -> str:
    lines = ["# STROBE-MR checklist", "",
             "Generated by GeNor from the run manifest. Items marked **HUMAN** need the authors.", "",
             "| # | Item | Status | Content |", "|---|---|---|---|"]
    for it in items:
        status = {"filled": "filled", "partial": "partial (HUMAN)", "human": "**HUMAN**"}[it["status"]]
        lines.append(f"| {it['item']} | {it['name']} | {status} | {it['text'].replace('|', '/')} |")
    lines += ["", "Reference: Skrivankova VW et al. JAMA 2021;326(16):1614-1621."]
    return "\n".join(lines) + "\n"
