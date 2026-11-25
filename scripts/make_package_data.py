"""Regenerate GeNor's bundled data (maintainers only).

Inputs are the public reference files listed in src/genor/refdata/registry.toml:
  LDscore/ (ldsc-ldscores-eur), 1000G_Phase3_weights_hm3_no_MHC/ (ldsc-weights-eur),
  1000G_Phase3_frq/ (Zenodo 10515792 1000G_Phase3_frq.tgz), hg19ToHg38.over.chain.gz

Usage: python scripts/make_package_data.py <dir-with-those-files>

Writes:
  src/genor/data/build_markers.tsv.gz   rsID with GRCh37 and GRCh38 positions (build inference)
  src/genor/data/example/*.tsv.gz        simulated exposure/outcome GWAS (see example/README.md)
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from genor.sumstats.liftover import lift, load_chain  # noqa: E402

SEED = 20261125


def read_chr(d: Path) -> pd.DataFrame:
    ld = pd.concat(pd.read_csv(d / f"LDscore/LDscore.{c}.l2.ldscore.gz", sep="\t") for c in range(1, 23))
    wl = pd.concat(pd.read_csv(d / f"1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.{c}.l2.ldscore.gz",
                               sep="\t") for c in range(1, 23))
    frq = pd.concat(pd.read_csv(d / f"1000G_Phase3_frq/1000G.EUR.QC.{c}.frq", sep=r"\s+")
                    for c in range(1, 23))
    m = ld.merge(wl[["SNP"]], on="SNP").merge(frq[["SNP", "A1", "A2", "MAF"]], on="SNP")
    m = m[m["A1"].isin(list("ACGT")) & m["A2"].isin(list("ACGT")) & (m["MAF"] >= 0.01)]
    return m.reset_index(drop=True)


def build_markers(m: pd.DataFrame, d: Path, rng: np.random.Generator) -> None:
    sub = m.sample(4000, random_state=1).copy()
    chain = load_chain(d / "hg19ToHg38.over.chain.gz")
    c38, p38 = lift(chain, sub["CHR"].astype(str), sub["BP"])
    sub["pos38"] = p38
    sub = sub[p38.notna().to_numpy() & (sub["BP"].to_numpy() != sub["pos38"].to_numpy())]
    out = pd.DataFrame({"rsid": sub["SNP"], "chromosome": sub["CHR"].astype(str),
                        "pos37": sub["BP"].astype(int), "pos38": sub["pos38"].astype(int)})
    out.sort_values(["chromosome", "pos37"]).to_csv(
        ROOT / "src/genor/data/build_markers.tsv.gz", sep="\t", index=False, compression="gzip")
    print("build markers:", len(out))


def simulate(m: pd.DataFrame, d: Path, rng: np.random.Generator) -> None:
    M = sum(float(open(d / f"LDscore/LDscore.{c}.l2.M_5_50").read().split()[0]) for c in range(1, 23))
    snps = m.iloc[::12].reset_index(drop=True)  # ~1 in 12 HapMap3 SNPs
    n = len(snps)
    l = snps["L2"].clip(lower=1).to_numpy()  # noqa: E741
    f = snps["MAF"].to_numpy()  # frequency of A1 (minor allele) in 1000G EUR

    truth = {"exposure": "BMI_sim", "outcome": "MDD_sim", "h2_exposure": 0.20,
             "h2_outcome_observed": 0.06, "rg": 0.15, "causal_effect": 0.15,
             "n_exposure": 300_000, "n_cases": 60_000, "n_controls": 140_000,
             "n_instruments_planted": 40, "n_pleiotropic_instruments": 3, "M_5_50": M,
             "seed": SEED, "n_variants": n, "build": "GRCh37", "ancestry": "EUR"}
    n1 = truth["n_exposure"]
    n2 = truth["n_cases"] + truth["n_controls"]
    h1, h2, rg = truth["h2_exposure"], truth["h2_outcome_observed"], truth["rg"]
    v1 = 1 + n1 * h1 * l / M
    v2 = 1 + n2 * h2 * l / M
    cov = np.sqrt(n1 * n2) * rg * np.sqrt(h1 * h2) * l / M
    u1, u2 = rng.standard_normal(n), rng.standard_normal(n)
    z1 = np.sqrt(v1) * u1
    r = cov / np.sqrt(v1 * v2)
    z2 = np.sqrt(v2) * (r * u1 + np.sqrt(1 - r**2) * u2)

    # plant instruments: >2 Mb apart, common
    cand = snps[(snps["MAF"] > 0.1)].sample(frac=1, random_state=SEED)
    chosen, last = [], {}
    for i, row in cand.iterrows():
        key = row["CHR"]
        if all(abs(row["BP"] - bp) > 2_000_000 for bp in last.get(key, [])):
            chosen.append(i)
            last.setdefault(key, []).append(row["BP"])
        if len(chosen) == truth["n_instruments_planted"]:
            break
    chosen = np.array(chosen)
    neff2 = 4 / (1 / truth["n_cases"] + 1 / truth["n_controls"])
    se1 = 1 / np.sqrt(2 * f * (1 - f) * n1)
    se2 = 1 / np.sqrt(2 * f * (1 - f) * neff2)
    bx = rng.choice([-1, 1], len(chosen)) * rng.uniform(0.025, 0.06, len(chosen))
    by = truth["causal_effect"] * bx
    pleio = chosen[:3]
    by[:3] += rng.choice([-1, 1], 3) * 0.03
    z1[chosen] = bx / se1[chosen] + rng.standard_normal(len(chosen))
    z2[chosen] = by / se2[chosen] + rng.standard_normal(len(chosen))
    truth["pleiotropic_snps"] = snps.loc[pleio, "SNP"].tolist()
    truth["instrument_snps"] = snps.loc[chosen, "SNP"].tolist()

    from scipy import stats

    b1, b2 = z1 * se1, z2 * se2
    p1, p2 = 2 * stats.norm.sf(np.abs(z1)), 2 * stats.norm.sf(np.abs(z2))
    # exposure: PLINK/GIANT-like, effect allele = A1 (minor)
    expo = pd.DataFrame({"SNP": snps["SNP"], "CHR": snps["CHR"], "BP": snps["BP"],
                         "A1": snps["A1"], "A2": snps["A2"], "FRQ": np.round(f, 4),
                         "BETA": np.round(b1, 6), "SE": np.round(se1, 6),
                         "P": [f"{x:.4g}" for x in p1], "N": n1})
    # outcome: METAL-like names, lower-case alleles, effect allele swapped for ~40% of SNPs
    swap = rng.random(n) < 0.4
    ea = np.where(swap, snps["A2"], snps["A1"])
    oa = np.where(swap, snps["A1"], snps["A2"])
    outc = pd.DataFrame({"MarkerName": snps["SNP"], "chrom": snps["CHR"], "pos": snps["BP"],
                         "Allele1": pd.Series(ea).str.lower(), "Allele2": pd.Series(oa).str.lower(),
                         "Freq1": np.round(np.where(swap, 1 - f, f), 4),
                         "Effect": np.round(np.where(swap, -b2, b2), 6), "StdErr": np.round(se2, 6),
                         "P-value": [f"{x:.4g}" for x in p2],
                         "Ncases": truth["n_cases"], "Ncontrols": truth["n_controls"]})
    # a few malformed rows to exercise QC
    bad = outc.sample(25, random_state=2).index
    outc.loc[bad[:10], "P-value"] = "NA"
    outc.loc[bad[10:], "StdErr"] = 0
    out = ROOT / "src/genor/data/example"
    out.mkdir(parents=True, exist_ok=True)
    expo.to_csv(out / "bmi_sim.tsv.gz", sep="\t", index=False, compression="gzip")
    outc.to_csv(out / "mdd_sim.tsv.gz", sep="\t", index=False, compression="gzip")
    (out / "truth.json").write_text(json.dumps(truth, indent=2))
    print("example variants:", n, "instruments:", len(chosen))


def main() -> None:
    d = Path(sys.argv[1])
    rng = np.random.default_rng(SEED)
    m = read_chr(d)
    build_markers(m, d, rng)
    simulate(m, d, rng)


if __name__ == "__main__":
    main()
