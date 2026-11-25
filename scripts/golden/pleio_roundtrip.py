"""Round-trip check: GeNor's template alignment + pleioFDR port == pleioFDR on the original MAT files.

Needs the pleioFDR demo data (pleioFDR_demo_data.tar.gz, 21 MB) and the 9545380.ref template
(262 MB) in the working directory, and the `pleiofdr` package installed (Python >= 3.13).
Verified 2026-11-25 on the chr21 demo: identical *_all.csv (287 rows, 8 loci, max diff 0).

    python scripts/golden/pleio_roundtrip.py
"""

import glob
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.io
from pleiofdr import io
from pleiofdr.cli import run_from_config

from genor.commands import pleio_run
from genor.sumstats.load import TraitMeta, save_trait

t = pd.read_csv("9545380.ref", sep="\t", dtype={"SNP": str, "A1": str, "A2": str})
t21 = t[t.CHR == 21].reset_index(drop=True)
t21.to_csv("template_chr21.ref", sep="\t", index=False)
rng = np.random.default_rng(3)
for f, name in (("CTG_COG_2018_DEMO.mat", "COG"), ("SSGAC_EDU_2016_DEMO.mat", "EDU")):
    m = io.load_mat(f)
    lp, z = np.ravel(m["logpvec"]), np.ravel(m["zvec"])
    ok = np.isfinite(lp) & np.isfinite(z)
    flip = rng.random(len(z)) < 0.5  # report half the SNPs on the other allele
    d = pd.DataFrame({"variant_id": t21.SNP, "chromosome": "21", "base_pair_location": t21.BP,
                      "effect_allele": np.where(flip, t21.A2, t21.A1),
                      "other_allele": np.where(flip, t21.A1, t21.A2),
                      "z": np.where(flip, -z, z), "p_value": 10 ** -lp})[ok]
    save_trait(d, TraitMeta(trait=name, build="GRCh37", ancestry="EUR"), Path(f"tr_{name}"))
    lp = lp.copy()
    lp[~np.isfinite(z)] = np.nan  # same SNP set for the direct run
    scipy.io.savemat(f"masked_{f}", {"logpvec": lp.reshape(-1, 1), "zvec": z.reshape(-1, 1)})

pleio_run(Path("tr_COG"), Path("tr_EDU"), Path("g_run"), stattype="condfdr", randprune_n=20, seed=7,
          ref_mat=Path("ref_1kgPhase3eur_LDr2p1_DEMO.mat"), template=Path("template_chr21.ref"))
gd = str(Path("g_run/pleiofdr").resolve())
cfg = (Path("g_run/pleiofdr/config.txt").read_text()
       .replace(f"traitfolder={gd}", f"traitfolder={Path('.').resolve()}")
       .replace("traitfile1=COG.mat", "traitfile1=masked_CTG_COG_2018_DEMO.mat")
       .replace("{'EDU.mat'}", "{'masked_SSGAC_EDU_2016_DEMO.mat'}")
       .replace(f"outputdir={gd}/results", f"outputdir={Path('direct_results').resolve()}"))
Path("direct.txt").write_text(cfg)
run_from_config("direct.txt", seed=7, make_plots=False)
A = pd.read_csv(glob.glob("g_run/pleiofdr/results/*_all.csv")[0])
B = pd.read_csv(glob.glob("direct_results/*_all.csv")[0])
assert A.shape == B.shape, (A.shape, B.shape)
for c in [c for c in A.columns if A[c].dtype.kind == "f" and c != "geneid"]:
    assert np.allclose(A[c], B[c], equal_nan=True), c
print("identical:", A.shape)
