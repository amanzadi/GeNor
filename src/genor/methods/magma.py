"""MAGMA driver (gene, gene-set, gene-property analysis).

MAGMA's licence forbids redistribution, so the user downloads it once and registers
it (``genor backend add magma /path/to/magma``). The auxiliary files (gene locations,
reference genotypes ``g1000_<anc>``) come from the MAGMA website as well and are passed
by path.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from genor import backends
from genor.errors import GenorError, Warnings
from genor.sumstats.load import TraitMeta


def gene_analysis(df: pd.DataFrame, meta: TraitMeta, work: Path, w: Warnings, *, gene_loc: Path,
                  ref_bfile: str, window_kb: tuple[int, int] = (35, 10)) -> dict[str, Any]:
    backends.require("magma")
    if meta.n is None:
        raise GenorError("MAGMA needs N.", code="n_missing")
    for p in (gene_loc, Path(f"{ref_bfile}.bed")):
        if not Path(p).exists():
            raise GenorError(f"{p} not found", code="file_not_found",
                             hint="Download MAGMA auxiliary files from https://cncr.nl/research/magma/")
    work.mkdir(parents=True, exist_ok=True)
    snploc = work / "snp.loc"
    df[["variant_id", "chromosome", "base_pair_location"]].dropna().to_csv(
        snploc, sep="\t", index=False, header=False)
    pv = work / "pvals.txt"
    df[["variant_id", "p_value"]].assign(N=df["n"] if "n" in df else meta.n).rename(
        columns={"variant_id": "SNP", "p_value": "P"}).to_csv(pv, sep="\t", index=False)
    up, down = window_kb
    backends.run("magma", ["--annotate", f"window={up},{down}", "--snp-loc", str(snploc),
                           "--gene-loc", str(gene_loc), "--out", str(work / "annot")],
                 log=work / "magma_annotate.log")
    backends.run("magma", ["--bfile", ref_bfile, "--pval", str(pv), "use=SNP,P", "ncol=N",
                           "--gene-annot", str(work / "annot.genes.annot"), "--out", str(work / "genes")],
                 log=work / "magma_gene.log")
    res = pd.read_csv(work / "genes.genes.out", sep=r"\s+")
    n = len(res)
    bonf = 0.05 / n if n else float("nan")
    return {"n_genes": n, "bonferroni_p": bonf, "n_sig": int((res["P"] < bonf).sum()),
            "top_genes": res.nsmallest(10, "P")[["GENE", "P"]].to_dict("records"),
            "genes_raw": str(work / "genes.genes.raw"), "genes_out": str(work / "genes.genes.out")}


def gene_set(genes_raw: Path, set_annot: Path, work: Path, covar: Path | None = None) -> dict[str, Any]:
    backends.require("magma")
    args = ["--gene-results", str(genes_raw)]
    if covar:
        args += ["--gene-covar", str(covar)]
    else:
        args += ["--set-annot", str(set_annot)]
    args += ["--out", str(work / ("geneprop" if covar else "geneset"))]
    work.mkdir(parents=True, exist_ok=True)
    backends.run("magma", args, log=work / "magma_sets.log")
    out = work / (("geneprop" if covar else "geneset") + ".gsa.out")
    res = pd.read_csv(out, sep=r"\s+", comment="#")
    return {"n_sets": len(res), "top": res.nsmallest(10, "P").to_dict("records"), "gsa_out": str(out)}
