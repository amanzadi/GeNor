"""Column-name detection for summary statistics files.

Internal names follow GWAS-SSF. Header names are normalised (lower case, only
letters and digits) before lookup, so ``P-value``, ``p_value`` and ``PVALUE`` all match.
"""

from __future__ import annotations

import re

# Order matters inside each list only for documentation; lookup is by set membership.
ALIASES: dict[str, list[str]] = {
    "variant_id": ["variantid", "snp", "rsid", "rs", "rsids", "markername", "snpid", "rsnumber",
                   "marker", "variant", "snpname", "id", "rsnumbers", "hm_rsid"],
    "chromosome": ["chromosome", "chr", "chrom", "chrname", "hg19chrc", "chromosomenumber",
                   "hm_chrom", "chrid"],
    "base_pair_location": ["basepairlocation", "bp", "pos", "position", "bphg19", "bpos", "genpos",
                           "hm_pos", "posb37", "posb38", "basepair"],
    "effect_allele": ["effectallele", "a1", "ea", "allele1", "testedallele", "codedallele",
                      "incallele", "riskallele", "alt", "altallele", "hm_effect_allele",
                      "alleleb", "effectalleles"],
    "other_allele": ["otherallele", "a2", "oa", "allele2", "allele0", "noneffectallele",
                     "noncodedallele", "decallele", "nea", "ref", "refallele",
                     "hm_other_allele", "allelea", "otheralleles"],
    "beta": ["beta", "b", "effect", "effects", "logor", "logodds", "betasnp", "est", "hm_beta",
             "effectsize", "allelicbeta"],
    "odds_ratio": ["oddsratio", "or", "hm_odds_ratio"],
    "standard_error": ["standarderror", "se", "stderr", "sebeta", "sesnp", "stderror", "logorse"],
    "p_value": ["pvalue", "p", "pval", "pbolt_lmm", "pboltlmm", "pboltlmminf", "pvalnominal",
                "frequentistaddpvalue", "pvalgc", "gcpvalue", "pdgc", "pvalues", "pv"],
    "neg_log10_p": ["neglog10pvalue", "log10p", "mlog10p", "neglog10pval", "neglog10p",
                    "minuslog10p", "lp", "mlogp"],
    "effect_allele_frequency": ["effectallelefrequency", "eaf", "freq", "frq", "af", "a1freq",
                                "freq1", "freqa1", "a1frq", "af1", "afalt", "freqeffect",
                                "hm_effect_allele_frequency", "frequency", "effectallelefreq",
                                "allelefreq", "eafukb"],
    "n": ["n", "ntotal", "totaln", "samplesize", "nsample", "nsamples", "obsct", "neff", "nobs",
          "nmiss", "weight", "totalsamplesize"],
    "n_cases": ["ncases", "ncas", "nca", "ncase", "cases"],
    "n_controls": ["ncontrols", "ncon", "nco", "ncontrol", "controls"],
    "z": ["z", "zscore", "zstat", "zstatistic", "tstat"],
    "info": ["info", "imputationinfo", "infoscore", "impinfo", "rsq", "r2"],
}

AMBIGUOUS_NOTES = {
    "a1": "In PLINK and LDSC output A1 is the effect allele; GeNor assumes the same.",
    "rsq": "Treated as imputation quality.",
    "neff": "Effective sample size used as N.",
    "weight": "METAL 'Weight' column used as N.",
}



def normalise(name: str) -> str:
    return re.sub(r"[^a-z0-9]", "", name.lower())


_LOOKUP = {normalise(alias): canon for canon, aliases in ALIASES.items() for alias in aliases}


def detect(header: list[str], overrides: dict[str, str] | None = None) -> tuple[dict[str, str], list[str]]:
    """Map file columns to canonical names.

    Returns (mapping file_column -> canonical, notes). Explicit ``overrides``
    (canonical -> file column) win.
    """
    overrides = dict(overrides or {})
    norm = {normalise(h): h for h in header}
    # SAIGE: BETA refers to Allele2 (AF_Allele2 / AC_Allele2 columns identify the format)
    if "afallele2" in norm or "acallele2" in norm:
        if "allele2" in norm and "allele1" in norm:
            overrides.setdefault("effect_allele", norm["allele2"])
            overrides.setdefault("other_allele", norm["allele1"])
        if "afallele2" in norm:
            overrides.setdefault("effect_allele_frequency", norm["afallele2"])
    # PLINK 2 .glm: A1 is the tested (effect) allele; REF/ALT are genome alleles. The other
    # allele is derived later (see load.standardise), so don't map REF/ALT here.
    plink2_glm = "a1" in norm and "ref" in norm and "alt" in norm
    if plink2_glm:
        overrides.setdefault("effect_allele", norm["a1"])
    mapping: dict[str, str] = {}
    notes: list[str] = []
    taken: set[str] = set()
    for canon, col in overrides.items():
        if col not in header:
            raise KeyError(f"Column '{col}' (given for {canon}) is not in the file header")
        mapping[col] = canon
        taken.add(canon)
    # regenie: ALLELE1 is the effect allele, ALLELE0 the other.
    for col in header:
        if col in mapping:
            continue
        key = normalise(col)
        canon = _LOOKUP.get(key)
        if plink2_glm and key in ("ref", "alt", "a1freq"):
            continue
        if canon and canon not in taken:
            mapping[col] = canon
            taken.add(canon)
            if key in AMBIGUOUS_NOTES:
                notes.append(f"{col}: {AMBIGUOUS_NOTES[key]}")
    if plink2_glm:
        notes.append("PLINK 2 .glm format: A1 is the effect allele; other allele = whichever of "
                     "REF/ALT is not A1.")
    if "afallele2" in norm or "acallele2" in norm:
        notes.append("SAIGE format: Allele2 is the effect allele (BETA refers to Allele2).")
    return mapping, notes
