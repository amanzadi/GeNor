"""condFDR / conjFDR through the pleioFDR Python port (no MATLAB).

GeNor's part: align each trait to the pleioFDR 9,545,380-SNP template (``logpvec`` /
``zvec`` MAT files, z signed to the template's A1), write a ``config.txt`` and call
``pleiofdr.cli.run_from_config``. The statistics are pleioFDR's own code.

Requirements: the ``pleiofdr`` backend (Python >= 3.13), the template
(``pleiofdr-template``, 262 MB) and the LD reference (``pleiofdr-ref-eur``, 2.3 GB,
~12 GB RAM when loaded). The reference is EUR / GRCh37.
"""

from __future__ import annotations

import glob
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from genor import backends, refdata
from genor.errors import GenorError, Warnings
from genor.guards import check_ancestry
from genor.sumstats.load import TraitMeta

COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}


def read_template(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", usecols=["CHR", "SNP", "BP", "A1", "A2"],
                       dtype={"SNP": str, "A1": str, "A2": str})


def align_to_template(df: pd.DataFrame, tmpl: pd.DataFrame, w: Warnings, label: str) -> dict[str, np.ndarray]:
    if "variant_id" not in df:
        raise GenorError(f"{label}: pleioFDR matches the template by rsID; none found.",
                         code="pleio_needs_rsid")
    d = df.drop_duplicates("variant_id").set_index("variant_id")
    idx = d.index.get_indexer(tmpl["SNP"])
    hit = idx >= 0
    ea = np.full(len(tmpl), "", dtype=object)
    oa = np.full(len(tmpl), "", dtype=object)
    z = np.full(len(tmpl), np.nan)
    p = np.full(len(tmpl), np.nan)
    ea[hit] = d["effect_allele"].to_numpy()[idx[hit]]
    oa[hit] = d["other_allele"].to_numpy()[idx[hit]]
    z[hit] = d["z"].to_numpy(float)[idx[hit]]
    p[hit] = d["p_value"].to_numpy(float)[idx[hit]]
    a1, a2 = tmpl["A1"].to_numpy(), tmpl["A2"].to_numpy()
    cea = pd.Series(ea).map(COMP).to_numpy()
    coa = pd.Series(oa).map(COMP).to_numpy()
    # direct matches win over strand complements (palindromes match both ways)
    same_d = (ea == a1) & (oa == a2)
    swap_d = (ea == a2) & (oa == a1)
    same = same_d | ((cea == a1) & (coa == a2) & ~swap_d)
    swap = (swap_d | ((cea == a2) & (coa == a1) & ~same_d)) & ~same
    z = np.where(swap, -z, z)
    bad = hit & ~(same | swap)
    z[bad] = np.nan
    p[bad] = np.nan
    n_match = int((hit & (same | swap)).sum())
    if bad.sum():
        w.add("pleio_allele_mismatch", f"{label}: {int(bad.sum())} variants with alleles not "
              "matching the template set to missing.", "info")
    if n_match < 0.3 * len(df):
        w.add("pleio_low_overlap", f"{label}: only {n_match} of {len(df)} variants matched the "
              "pleioFDR template.", "serious")
    return {"logpvec": (-np.log10(p)).reshape(-1, 1), "zvec": z.reshape(-1, 1), "n_matched": n_match}


def write_mat(path: Path, arrays: dict[str, Any]) -> None:
    import scipy.io

    scipy.io.savemat(path, {"logpvec": arrays["logpvec"], "zvec": arrays["zvec"]}, do_compression=False)


def write_config(path: Path, *, reffile: Path, refinfo: Path, trait_dir: Path, t1: str, n1: str,
                 t2: str, n2: str, outputdir: Path, stattype: str, fdrthresh: float,
                 randprune_n: int, exclude_mhc: bool) -> None:
    reffile, refinfo, trait_dir, outputdir = (Path(x).resolve() for x in (reffile, refinfo, trait_dir, outputdir))
    lines = [
        f"reffile={reffile}", f"refinfo={refinfo}", f"traitfolder={trait_dir}",
        f"traitfile1={t1}", f"traitname1={n1}", f"traitfiles={{'{t2}'}}", f"traitnames={{'{n2}'}}",
        f"outputdir={outputdir}", f"stattype={stattype}", f"fdrthresh={fdrthresh}",
        f"randprune_n={randprune_n}", "onscreen=false",
        "exclude_chr_pos=[6 25119106 33854733]" if exclude_mhc else "exclude_chr_pos=[]",
    ]
    path.write_text("\n".join(lines) + "\n")


def run(d1: pd.DataFrame, m1: TraitMeta, d2: pd.DataFrame, m2: TraitMeta, work: Path, w: Warnings, *,
        stattype: str = "condfdr", fdrthresh: float | None = None, randprune_n: int = 500,
        exclude_mhc: bool = True, seed: int | None = 1, prepare_only: bool = False,
        ref_mat: Path | None = None, template: Path | None = None) -> dict[str, Any]:
    if stattype not in ("condfdr", "conjfdr"):
        raise GenorError("stattype must be condfdr or conjfdr", code="bad_param")
    fdrthresh = fdrthresh if fdrthresh is not None else (0.01 if stattype == "condfdr" else 0.05)
    for m in (m1, m2):
        check_ancestry(m.ancestry, "EUR", "the pleioFDR reference (1000G EUR)")
    tmpl_path = Path(template) if template else refdata.file_in("pleiofdr-template")
    if ref_mat:
        ref_path = Path(ref_mat)
    elif prepare_only and not refdata.is_installed("pleiofdr-ref-eur"):
        ref_path = refdata.dataset_dir("pleiofdr-ref-eur") / refdata.entry("pleiofdr-ref-eur")["filename"]
        w.add("pleio_ref_not_downloaded", "MAT files and config written; the 2.4 GB reference is "
              "not downloaded yet, so config.txt points at where it will be.", "info")
    else:
        ref_path = refdata.file_in("pleiofdr-ref-eur")
    for p in (tmpl_path, ref_path) if not prepare_only else (tmpl_path,):
        if not p.exists():
            raise GenorError(f"{p} not found", code="file_not_found")
    if not prepare_only:
        backends.require("pleiofdr")
    tmpl = read_template(tmpl_path)
    work.mkdir(parents=True, exist_ok=True)
    names = []
    for d, m in ((d1, m1), (d2, m2)):
        arr = align_to_template(d, tmpl, w, m.trait)
        safe = "".join(ch if ch.isalnum() else "_" for ch in m.trait)
        write_mat(work / f"{safe}.mat", arr)
        names.append((safe, arr["n_matched"]))
    if exclude_mhc:
        in_mhc = (tmpl["CHR"] == 6) & tmpl["BP"].between(25_119_106, 33_854_733)
        if not in_mhc.any():
            exclude_mhc = False
            w.add("pleio_mhc_not_in_template", "Template has no MHC variants; no exclusion applied.", "info")
    cfg = work / "config.txt"
    outdir = work / "results"
    write_config(cfg, reffile=ref_path, refinfo=tmpl_path, trait_dir=work, t1=f"{names[0][0]}.mat",
                 n1=names[0][0], t2=f"{names[1][0]}.mat", n2=names[1][0], outputdir=outdir,
                 stattype=stattype, fdrthresh=fdrthresh, randprune_n=randprune_n,
                 exclude_mhc=exclude_mhc)
    summary: dict[str, Any] = {"stattype": stattype, "fdrthresh": fdrthresh, "randprune_n": randprune_n,
                               "matched": {n: k for n, k in names}, "config": str(cfg)}
    if prepare_only:
        summary["prepared_only"] = True
        return summary
    from pleiofdr.cli import run_from_config  # type: ignore[import-not-found]

    out = run_from_config(cfg, seed=seed, make_plots=True)
    summary["n_loci"] = float(out.locusnum) if out.locusnum is not None else None
    summary["result_files"] = sorted(Path(p).name for p in glob.glob(str(outdir / "*")))
    loci = sorted(outdir.glob("*_loci.csv"))
    if loci:
        summary["loci_table"] = str(loci[0])
    if randprune_n < 100:
        w.add("pleio_few_prunes", f"randprune_n={randprune_n}; publications use >= 500.", "warning")
    return summary
