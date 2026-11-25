import json
import sys

sys.path.insert(0, sys.argv[1])
import ldsc
from ldscore import sumstats as S

parser = ldsc.parser
R = sys.argv[2]
out = {}
for mode, val in (("h2", "bmi.sumstats"), ("rg", "bmi.sumstats,mdd.sumstats")):
    args = parser.parse_args([f"--{mode}", val, "--ref-ld-chr", f"{R}/ldsc-ldscores-eur/LDscore/LDscore.",
        "--w-ld-chr", f"{R}/ldsc-weights-eur/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.", "--out", f"g_{mode}"])
    log = ldsc.Logger(f"g_{mode}.log")
    if mode == "h2":
        h = S.estimate_h2(args, log)
        out["h2"] = {"h2": float(h.tot), "se": float(h.tot_se), "intercept": float(h.intercept), "intercept_se": float(h.intercept_se)}
    else:
        r = S.estimate_rg(args, log)[0]
        out["rg"] = {"rg": float(r.rg_ratio), "se": float(r.rg_se), "h2_1": float(r.hsq1.tot), "h2_1_se": float(r.hsq1.tot_se),
                     "h2_2": float(r.hsq2.tot), "h2_2_se": float(r.hsq2.tot_se), "gencov": float(r.gencov.tot),
                     "gencov_se": float(r.gencov.tot_se), "gencov_intercept": float(r.gencov.intercept)}
print(json.dumps(out, indent=1))
