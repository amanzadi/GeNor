"""Command implementations shared by the CLI and the MCP server.

Each analysis function opens a :class:`~genor.runs.Run`, does the work, records
inputs/parameters/warnings/citations, and returns the short result dict that the
CLI prints and MCP tools return.
"""

from __future__ import annotations

import json
import platform
import shutil
import sys
from dataclasses import asdict
from importlib import resources
from pathlib import Path
from typing import Any

import pandas as pd

from genor import RESEARCH_USE_NOTICE, __version__, backends, refdata
from genor.errors import GenorError, NeedsUserDecision
from genor.runs import Run
from genor.sumstats.load import TraitMeta, load, open_trait, save_trait

# ============================================================== environment


def doctor() -> dict[str, Any]:
    bk = backends.status()
    rd = refdata.status()
    blocked: dict[str, list[str]] = {}
    for b in bk:
        if not b["available"]:
            for m in b["provides"]:
                blocked.setdefault(m, []).append(f"backend {b['name']}: {b['how_to_enable']}")
    for r in rd:
        if not r["installed"]:
            for m in r["provides"]:
                blocked.setdefault(m, []).append(f"refdata {r['name']} ({r['size']}): genor refdata fetch {r['name']}")
    ready_now = ["sumstats.load", "sumstats.qc", "sumstats.liftover*", "meta", "mr (with --clump-mode distance)",
                 "finemap.abf", "coloc.abf", "report", "lookup.*", "example"]
    from genor.config import data_dir, jobs_dir
    return {
        "status": "ok",
        "genor_version": __version__,
        "python": sys.version.split()[0],
        "platform": backends.platform_key(),
        "supported_platform": platform.system() in ("Linux", "Darwin"),
        "data_dir": str(data_dir()),
        "jobs_dir": str(jobs_dir()),
        "backends": [{k: v for k, v in b.items() if k in ("name", "available", "version", "how_to_enable", "mode")}
                     for b in bk],
        "refdata": [{k: v for k, v in r.items() if k in ("name", "installed", "size", "ancestry", "build")}
                    for r in rd],
        "works_without_downloads": ready_now,
        "blocked_until_setup": blocked,
        "policy": "Summary statistics only. Users supply their own GWAS files; GeNor never "
                  "downloads GWAS or handles credentials.",
        "notice": RESEARCH_USE_NOTICE,
    }


def sources(name: str | None = None) -> dict[str, Any]:
    import tomllib

    src = tomllib.loads(resources.files("genor.sumstats").joinpath("sources.toml").read_text())
    if name:
        if name not in src:
            raise GenorError(f"Unknown source '{name}'", code="unknown_source", hint=", ".join(src))
        return {"status": "ok", "source": {"key": name, **src[name]}}
    return {"status": "ok", "note": "GeNor does not download GWAS. The user downloads the file "
            "(accepting the provider's terms) and passes a local path.",
            "sources": [{"key": k, "name": v["name"], "url": v["url"], "build": v.get("build"),
                         "access": v["access"]} for k, v in src.items()]}


def example(out: Path) -> dict[str, Any]:
    out.mkdir(parents=True, exist_ok=True)
    pkg = resources.files("genor.data").joinpath("example")
    files = {}
    for name in ("bmi_sim.tsv.gz", "mdd_sim.tsv.gz", "truth.json"):
        with resources.as_file(pkg.joinpath(name)) as p:
            shutil.copy(p, out / name)
            files[name] = str(out / name)
    truth = json.loads((out / "truth.json").read_text())
    return {
        "status": "ok", "files": files,
        "description": ("Simulated GWAS (not real data) on ~98k HapMap3 SNPs, GRCh37, EUR. "
                        "BMI_sim (quantitative, N=300k) and MDD_sim (case-control 60k/140k) with "
                        "a planted causal effect for MR and a genetic correlation for LDSC."),
        "truth": {k: truth[k] for k in ("h2_exposure", "h2_outcome_observed", "rg", "causal_effect",
                                        "n_instruments_planted")},
        "next_steps": [
            f"genor sumstats load {out}/bmi_sim.tsv.gz --trait BMI_sim --ancestry EUR --out {out}/bmi",
            f"genor sumstats load {out}/mdd_sim.tsv.gz --trait MDD_sim --ancestry EUR --population-prevalence 0.15 --out {out}/mdd",
            f"genor ldsc rg {out}/bmi {out}/mdd --out {out}/rg   (needs: genor refdata fetch ldsc-ldscores-eur ldsc-weights-eur)",
            f"genor mr run --exposure {out}/bmi --outcome {out}/mdd --clump-mode distance --overlap none --out {out}/mr",
        ],
    }


# ============================================================== sumstats


def sumstats_load(file: Path, trait: str, out: Path, *, build: str | None = None,
                  ancestry: str | None = None, n: float | None = None, n_cases: float | None = None,
                  n_controls: float | None = None, population_prevalence: float | None = None,
                  effect_type: str | None = None, columns: dict[str, str] | None = None,
                  maf: float = 0.01, info: float = 0.8, drop_mhc: bool = False,
                  plots: bool = True) -> dict[str, Any]:
    from genor.report import plots as P
    from genor.sumstats.qc import run_qc

    params = dict(file=file, trait=trait, build=build, ancestry=ancestry, n=n, n_cases=n_cases,
                  n_controls=n_controls, population_prevalence=population_prevalence,
                  effect_type=effect_type, columns=columns, maf=maf, info=info, drop_mhc=drop_mhc)
    with Run("sumstats.load", out, params) as run:
        run.add_input(file, "sumstats")
        run.cite("gwas-ssf")
        df, meta, _ = load(Path(file), trait, run.dir, build=build, ancestry=ancestry, n=n,
                           n_cases=n_cases, n_controls=n_controls,
                           population_prevalence=population_prevalence, effect_type=effect_type,
                           columns=columns, w=run.warnings)
        meta.source_sha256 = run.inputs[0].get("sha256")
        clean, report, _ = run_qc(df, meta.build, maf=maf, info=info, drop_mhc=drop_mhc, w=run.warnings)
        meta.qc = report
        save_trait(clean, meta, run.dir)
        run.output("sumstats", run.dir / "sumstats.parquet")
        run.output("trait_meta", run.dir / "sumstats.json")
        if plots and len(clean):
            if "chromosome" in clean:
                run.output("manhattan", P.manhattan(clean, run.path("manhattan.png"), trait))
            run.output("qq", P.qq(clean["p_value"], run.path("qq.png"), trait, report["lambda_gc"]))
        run.meta = {"trait": asdict(meta)}
        run.summary = {"trait": trait, "n_variants": report["n_output"], "build": meta.build,
                       "build_inferred": meta.build_inferred, "ancestry": meta.ancestry,
                       "n": meta.n, "n_cases": meta.n_cases, "n_controls": meta.n_controls,
                       "trait_type": meta.trait_type, "lambda_gc": report["lambda_gc"],
                       "mean_chi2": report["mean_chi2"], "n_genome_wide_sig": report["n_genome_wide_sig"],
                       "columns_detected": meta.columns, "qc": report}
        if meta.build is None:
            run.ask("Which genome build is this file on (GRCh37 or GRCh38)?", ["GRCh37", "GRCh38"],
                    "Build could not be inferred; LD-based methods need it.")
        if meta.ancestry is None:
            run.ask("What is the genetic ancestry of this GWAS sample?", ["EUR", "EAS", "AFR", "AMR", "SAS"],
                    "Needed to choose an LD reference panel.")
        if meta.n is None:
            run.ask("What is the sample size (or case/control counts)?", [], "LDSC and MR diagnostics need N.")
        run.next_steps = [f"genor ldsc h2 {run.dir}", f"genor clump run {run.dir}"]
    return run.result()


def sumstats_liftover(trait_dir: Path, to: str, out: Path) -> dict[str, Any]:
    from genor.guards import norm_build
    from genor.sumstats.liftover import lift, load_chain

    target = norm_build(to)
    with Run("sumstats.liftover", out, {"trait_dir": trait_dir, "to": target}) as run:
        df, meta = open_trait(trait_dir)
        run.add_input(Path(trait_dir) / "sumstats.parquet", "sumstats")
        if meta.build is None:
            raise NeedsUserDecision("Source build unknown.", "Which build is the source on?",
                                    ["GRCh37", "GRCh38"])
        if meta.build == target:
            raise GenorError(f"Already on {target}.", code="nothing_to_do")
        name = "chain-hg19-to-hg38" if target == "GRCh38" else "chain-hg38-to-hg19"
        chain = load_chain(refdata.file_in(name))
        run.add_refdata(name, **refdata.info_for_manifest(name))
        run.cite("liftover")
        c, p = lift(chain, df["chromosome"].astype(str), df["base_pair_location"])
        ok = p.notna()
        df = df[ok.to_numpy()].copy()
        df["chromosome"] = c[ok].to_numpy()
        df["base_pair_location"] = p[ok].astype(int).to_numpy()
        meta.build = target
        meta.build_inferred = False
        save_trait(df, meta, run.dir)
        run.output("sumstats", run.dir / "sumstats.parquet")
        run.summary = {"from": "GRCh37" if target == "GRCh38" else "GRCh38", "to": target,
                       "n_mapped": int(ok.sum()), "n_unmapped": int((~ok).sum())}
        if (~ok).mean() > 0.05:
            run.warn("liftover_many_unmapped", f"{(~ok).mean():.1%} of variants did not map.", "warning")
    return run.result()


# ============================================================== meta-analysis


def meta_run(traits: list[Path], out: Path, scheme: str = "ivw", trait_name: str = "meta") -> dict[str, Any]:
    from genor.methods.meta import meta_analyse

    with Run("meta", out, {"traits": traits, "scheme": scheme}) as run:
        tabs, metas = [], []
        for t in traits:
            d, m = open_trait(t)
            run.add_input(Path(t) / "sumstats.parquet", "sumstats", trait=m.trait)
            tabs.append(d)
            metas.append(m)
        res, summary = meta_analyse(tabs, metas, scheme, run.warnings)
        run.cite("metal", "higgins-i2", "dersimonian-laird")
        nm = TraitMeta(trait=trait_name, build=metas[0].build,
                       ancestry=metas[0].ancestry if len({m.ancestry for m in metas}) == 1 else "multi",
                       n=float(res["n"].max()) if res["n"].notna().any() else None,
                       source="genor meta: " + ", ".join(m.trait for m in metas))
        keep = [c for c in ("variant_id", "chromosome", "base_pair_location", "effect_allele", "other_allele",
                            "beta", "standard_error", "p_value", "z", "n", "effect_allele_frequency") if c in res]
        save_trait(res[keep], nm, run.dir)
        res.to_csv(run.path("meta_results.tsv.gz"), sep="\t", index=False)
        run.output("sumstats", run.dir / "sumstats.parquet")
        run.output("table", run.path("meta_results.tsv.gz"))
        if scheme == "samplesize":
            run.warn("meta_no_beta", "The sample-size scheme gives Z and P only (no beta/SE): the "
                     "result can be clumped or used for LDSC, but not as MR or coloc input.", "warning")
        if any(m.trait_type == "binary" for m in metas):
            run.warn("meta_binary_counts", "Case/control counts are not carried into the meta-analysed "
                     "trait; give --n-cases/--n-controls when reusing it for liability-scale h2.", "info")
        run.summary = summary
        run.next_steps = [f"genor clump run {run.dir}"]
    return run.result()


# ============================================================== LDSC


def _ld_inputs(run: Run, ref: str, w: str):
    from genor.methods.ldsc import load_reference

    ld = load_reference(ref, w)
    for name in (ref, w):
        run.add_refdata(name, **refdata.info_for_manifest(name))
    run.cite("ldsc", "ldsc-ref-zenodo")
    return ld


def ldsc_h2(trait: Path, out: Path, *, ref: str = "ldsc-ldscores-eur", w: str = "ldsc-weights-eur",
            pop_prev: float | None = None, intercept: float | None = None) -> dict[str, Any]:
    from genor.methods import ldsc

    with Run("ldsc.h2", out, {"trait": trait, "ref": ref, "w": w, "pop_prev": pop_prev,
                              "intercept": intercept}) as run:
        df, meta = open_trait(trait)
        run.add_input(Path(trait) / "sumstats.parquet", "sumstats", trait=meta.trait)
        ld = _ld_inputs(run, ref, w)
        res = ldsc.h2(df, meta, ld, run.warnings, pop_prev=pop_prev, intercept=intercept)
        run.summary = {"trait": meta.trait, **res}
        if meta.trait_type == "binary" and res.get("h2_liability") is None:
            if meta.sample_prevalence is None:
                run.ask(f"How many cases and controls does the {meta.trait} GWAS have?", [],
                        "Sample prevalence is needed for liability-scale h2; reload the trait "
                        "with --n-cases/--n-controls.")
            else:
                run.ask(f"What is the population prevalence of {meta.trait}?", [],
                        "Needed to convert h2 to the liability scale (--pop-prev).")
        (run.path("h2.json")).write_text(json.dumps(res, indent=2))
        run.output("h2", run.path("h2.json"))
    return run.result()


def ldsc_rg(trait1: Path, trait2: Path, out: Path, *, ref: str = "ldsc-ldscores-eur",
            w: str = "ldsc-weights-eur", pop_prev1: float | None = None,
            pop_prev2: float | None = None) -> dict[str, Any]:
    from genor.methods import ldsc

    with Run("ldsc.rg", out, {"trait1": trait1, "trait2": trait2, "ref": ref, "w": w}) as run:
        d1, m1 = open_trait(trait1)
        d2, m2 = open_trait(trait2)
        run.add_input(Path(trait1) / "sumstats.parquet", "sumstats", trait=m1.trait)
        run.add_input(Path(trait2) / "sumstats.parquet", "sumstats", trait=m2.trait)
        ld = _ld_inputs(run, ref, w)
        run.cite("ldsc-rg")
        res = ldsc.rg(d1, m1, d2, m2, ld, run.warnings, pop_prev=(pop_prev1, pop_prev2))
        (run.path("rg.json")).write_text(json.dumps(res, indent=2))
        run.output("rg", run.path("rg.json"))
        run.summary = res
    return run.result()


# ============================================================== clumping


def _record_plink(run: Run) -> None:
    p = backends.find("plink2")
    if p:
        run.add_backend("plink2", path=p, version=backends._version("plink2", p))
    run.cite("plink2", "1000g")


def clump_run(trait: Path, out: Path, *, mode: str = "plink", preset: str = "loci",
              p1: float | None = None, r2: float | None = None, kb: int | None = None,
              panel: str | None = None) -> dict[str, Any]:
    from genor.methods.clump import clump

    with Run("clump", out, {"trait": trait, "mode": mode, "preset": preset, "p1": p1, "r2": r2,
                            "kb": kb, "panel": panel}) as run:
        df, meta = open_trait(trait)
        run.add_input(Path(trait) / "sumstats.parquet", "sumstats", trait=meta.trait)
        res, info = clump(df, meta, run.path("work"), mode=mode, preset=preset, p1=p1, r2=r2, kb=kb,
                          panel=panel, w=run.warnings)
        if mode == "plink":
            _record_plink(run)
            run.add_refdata(info["panel"], **refdata.info_for_manifest(info["panel"]))
        res.to_csv(run.path("index_variants.tsv"), sep="\t", index=False)
        run.output("index_variants", run.path("index_variants.tsv"))
        top = res.nsmallest(15, "p_value")
        run.summary = {**info, "top": [
            {"variant_id": r.get("variant_id"), "chr": r.get("chromosome"), "pos": r.get("base_pair_location"),
             "p": r["p_value"]} for _, r in top.iterrows()]}
        run.next_steps = ["genor finemap abf <trait> --lead <variant>",
                          f"genor lookup genes --region chr:start-end --build {meta.build or 'GRCh37'}"]
    return run.result()


# ============================================================== MR


def mr_run(exposure: Path, outcome: Path, out: Path, *, p_threshold: float = 5e-8,
           clump_mode: str = "plink", panel: str | None = None, clump_r2: float = 0.001,
           clump_kb: int = 10_000, instruments: Path | None = None, exclude_mhc: bool = False,
           overlap: str | None = None, nboot: int = 1000, seed: int = 1) -> dict[str, Any]:
    from genor.methods import mr
    from genor.methods.clump import clump
    from genor.report import plots as P
    from genor.report import strobe_mr
    from genor.sumstats.harmonise import harmonise
    from genor.sumstats.qc import mhc_mask

    params = dict(exposure=exposure, outcome=outcome, p_threshold=p_threshold, clump_mode=clump_mode,
                  panel=panel, clump_r2=clump_r2, clump_kb=clump_kb, instruments=instruments,
                  exclude_mhc=exclude_mhc, overlap=overlap, nboot=nboot, seed=seed)
    with Run("mr", out, params) as run:
        ex, mx = open_trait(exposure)
        oc, mo = open_trait(outcome)
        run.add_input(Path(exposure) / "sumstats.parquet", "exposure", trait=mx.trait)
        run.add_input(Path(outcome) / "sumstats.parquet", "outcome", trait=mo.trait)
        run.cite("mr-ivw", "mr-egger", "mr-median", "mr-mode", "mr-steiger", "twosamplemr", "strobe-mr")
        if mx.ancestry and mo.ancestry and mx.ancestry != mo.ancestry:
            run.warn("mr_ancestry_mismatch", f"Exposure is {mx.ancestry}, outcome {mo.ancestry}: "
                     "variant effects may not transport across ancestries.", "serious")
        summary: dict[str, Any] = {"exposure": mx.trait, "outcome": mo.trait}
        if instruments:
            ids = set(pd.read_csv(instruments, sep=None, engine="python").iloc[:, 0].astype(str))
            cand = ex[ex["variant_id"].astype(str).isin(ids)]
            summary["n_candidates"] = len(cand)
            inst = cand
            summary["clump"] = {"mode": "user-supplied instruments"}
        else:
            cand = ex[ex["p_value"] < p_threshold]
            summary["n_candidates"] = len(cand)
            if len(cand) == 0:
                raise GenorError(f"No exposure variants with p < {p_threshold}.", code="mr_no_instruments",
                                 hint="A relaxed threshold (e.g. 5e-6) weakens instruments; ask the user.")
            inst, cinfo = clump(ex, mx, run.path("clump"), mode=clump_mode, preset="mr", p1=p_threshold,
                                r2=clump_r2, kb=clump_kb, panel=panel, w=run.warnings)
            summary["clump"] = cinfo
            if clump_mode == "plink":
                _record_plink(run)
        if exclude_mhc:
            mm = mhc_mask(inst, mx.build)
            inst = inst[~mm]
            summary["n_mhc_excluded"] = int(mm.sum())
        summary["n_instruments_clumped"] = len(inst)
        h, hrep = harmonise(inst, oc, suffixes=("_exp", "_out"), w=run.warnings)
        summary["harmonisation"] = hrep
        summary["n_missing_in_outcome"] = int(len(inst) - hrep["n_shared"])
        if summary["n_missing_in_outcome"]:
            run.warn("instruments_missing_in_outcome", f"{summary['n_missing_in_outcome']} instruments "
                     "are absent from the outcome GWAS (GeNor does not substitute LD proxies).", "warning")
        tab = pd.DataFrame({
            "SNP": h["variant_id_exp"] if "variant_id_exp" in h else h["variant_id"],
            "effect_allele": h["effect_allele"], "other_allele": h["other_allele"],
            "beta_exp": h["beta_exp"], "se_exp": h["standard_error_exp"], "p_exp": h["p_value_exp"],
            "beta_out": h["beta_out"], "se_out": h["standard_error_out"], "p_out": h["p_value_out"],
            "n_exp": h["n_exp"] if "n_exp" in h else mx.n,
            "n_out": h["n_out"] if "n_out" in h else mo.n,
        }).reset_index(drop=True)
        tab.to_csv(run.path("harmonised.tsv"), sep="\t", index=False)
        run.output("harmonised", run.path("harmonised.tsv"))
        if len(tab) < 3:
            run.warn("mr_few_instruments", f"Only {len(tab)} instruments: sensitivity estimators "
                     "(Egger, median, mode) need >= 3.", "serious")
        res = mr.run_all(tab, nboot=nboot, seed=seed)
        results = res.pop("results")
        single, loo = res.pop("single"), res.pop("loo")
        results.to_csv(run.path("mr_results.tsv"), sep="\t", index=False)
        single.to_csv(run.path("single_snp.tsv"), sep="\t", index=False)
        run.output("results", run.path("mr_results.tsv"))
        run.output("single_snp", run.path("single_snp.tsv"))
        if len(loo):
            loo.to_csv(run.path("leave_one_out.tsv"), sep="\t", index=False)
            run.output("leave_one_out", run.path("leave_one_out.tsv"))
        summary.update(res)
        recs = results.to_dict("records")
        summary["estimates"] = [{k: r.get(k) for k in ("method", "nsnp", "b", "se", "pval")} for r in recs]
        main = next((r for r in recs if r["method"] in ("IVW", "Wald ratio")), None)
        if main and mo.trait_type == "binary":
            import math
            summary["main_or"] = math.exp(main["b"])
            summary["main_or_95ci"] = [math.exp(main["b"] - 1.96 * main["se"]),
                                       math.exp(main["b"] + 1.96 * main["se"])]
        # diagnostics -> warnings
        F = res.get("F", {})
        if F and F.get("n_weak_F_lt_10"):
            run.warn("weak_instruments", f"{F['n_weak_F_lt_10']} instruments have F < 10.", "warning")
        if res.get("heterogeneity_ivw", {}).get("Q_pval", 1) < 0.05:
            run.warn("mr_heterogeneity", "Cochran's Q p < 0.05: heterogeneous ratio estimates "
                     "(possible pleiotropy). Rely on random-effects IVW and robust methods.", "warning")
        if res.get("pleiotropy_egger", {}).get("pval", 1) < 0.05:
            run.warn("mr_directional_pleiotropy", "MR-Egger intercept p < 0.05: directional "
                     "pleiotropy; IVW may be biased.", "serious")
        if res.get("I2_GX") is not None and res["I2_GX"] < 0.9:
            run.warn("egger_low_i2gx", f"I2_GX = {res['I2_GX']:.2f} < 0.9: MR-Egger is diluted "
                     "towards the null (NOME violation).", "info")
        st = res.get("steiger")
        if st and (mx.trait_type == "binary" or mo.trait_type == "binary"):
            run.warn("steiger_binary_approx", "Steiger r2 uses b/SE/N for binary traits too; "
                     "TwoSampleMR's liability-scale r2 (needs prevalence) can differ. Treat the "
                     "direction test as approximate.", "info")
        if st and not st["correct_causal_direction"]:
            run.warn("steiger_reverse", "Instruments explain more variance in the outcome than the "
                     "exposure: possible reverse causation or weak instruments.", "serious")
        if recs and len(recs) > 2:
            signs = {(r["b"] > 0) for r in recs if r["method"] in ("IVW", "MR-Egger", "Weighted median", "Weighted mode")}
            summary["direction_consistent_across_methods"] = len(signs) == 1
        if mx.trait_type == "binary":
            run.warn("binary_exposure", "Binary exposure: the estimate is per unit log-odds of "
                     "liability to the exposure, not the effect of having it.", "info")
        # decisions for the human
        if overlap is None:
            run.ask("Do the exposure and outcome GWAS share participants (sample overlap)?",
                    ["none", "partial", "substantial", "unknown"],
                    "Overlap biases two-sample MR towards the observational estimate when "
                    "instruments are weak. The cross-trait LDSC intercept is evidence.")
        summary["sample_overlap"] = overlap or "not stated"
        # figures
        if len(tab) >= 2:
            run.output("scatter", P.mr_scatter(tab, results, run.path("scatter.png"), mx.trait, mo.trait))
            run.output("forest", P.forest(single["SNP"].tolist(), single["b"].to_numpy(),
                                          single["se"].to_numpy(), run.path("forest.png"),
                                          f"Wald ratio ({mo.trait} per {mx.trait})"))
            iv = next(r for r in recs if r["method"] == "IVW")
            run.output("funnel", P.funnel(single["b"].to_numpy(), single["se"].to_numpy(), iv["b"],
                                          run.path("funnel.png")))
        if len(loo):
            run.output("leave_one_out_plot", P.forest(loo["SNP_left_out"].tolist(), loo["b"].to_numpy(),
                                                      loo["se"].to_numpy(), run.path("leave_one_out.png"),
                                                      "IVW leaving out SNP"))
        run.meta = {"exposure": asdict(mx), "outcome": asdict(mo)}
        for side in ("exposure", "outcome"):
            run.meta[side].pop("qc", None)
            run.meta[side].pop("columns", None)
        run.summary = summary
        # STROBE-MR after the summary exists
        items = strobe_mr.build(run.manifest(), recs)
        run.path("strobe_mr.md").write_text(strobe_mr.to_markdown(items))
        run.path("strobe_mr.json").write_text(json.dumps(items, indent=2))
        run.output("strobe_mr", run.path("strobe_mr.md"))
        summary["strobe_mr"] = {s: sum(1 for i in items if i["status"] == s) for s in ("filled", "partial", "human")}
        run.next_steps = [f"genor report build {run.dir} --out {run.dir}/report.md"]
    return run.result()


# ============================================================== pleioFDR


def pleio_run(trait1: Path, trait2: Path, out: Path, *, stattype: str = "condfdr",
              fdr: float | None = None, randprune_n: int = 500, exclude_mhc: bool = True,
              seed: int = 1, prepare_only: bool = False, ref_mat: Path | None = None,
              template: Path | None = None) -> dict[str, Any]:
    from genor.methods import pleio

    with Run(f"pleio.{stattype}", out, {"trait1": trait1, "trait2": trait2, "stattype": stattype,
                                        "fdr": fdr, "randprune_n": randprune_n, "seed": seed}) as run:
        d1, m1 = open_trait(trait1)
        d2, m2 = open_trait(trait2)
        run.add_input(Path(trait1) / "sumstats.parquet", "primary", trait=m1.trait)
        run.add_input(Path(trait2) / "sumstats.parquet", "conditional", trait=m2.trait)
        for name, override in (("pleiofdr-template", template), ("pleiofdr-ref-eur", ref_mat)):
            if override:
                run.add_input(override, name)
            else:
                run.add_refdata(name, **refdata.info_for_manifest(name))
        p = backends.find("pleiofdr")
        if p:
            run.add_backend("pleiofdr", path=p, version=backends._version("pleiofdr", p))
        run.cite("pleiofdr", "conjfdr")
        run.summary = pleio.run(d1, m1, d2, m2, run.path("pleiofdr"), run.warnings, stattype=stattype,
                                fdrthresh=fdr, randprune_n=randprune_n, exclude_mhc=exclude_mhc,
                                seed=seed, prepare_only=prepare_only, ref_mat=ref_mat, template=template)
        run.output("pleiofdr_dir", run.path("pleiofdr"))
    return run.result()


# ============================================================== locus level


def _region_from(df: pd.DataFrame, region: str | None, lead: str | None, window_kb: int) -> str:
    if region:
        return region
    if lead:
        r = df[df["variant_id"].astype(str) == lead]
        if r.empty:
            raise GenorError(f"Lead variant {lead} not in the trait", code="not_found")
        c, pos = r["chromosome"].iloc[0], int(r["base_pair_location"].iloc[0])
        return f"{c}:{max(1, pos - window_kb * 1000)}-{pos + window_kb * 1000}"
    raise GenorError("Give --region chr:start-end or --lead <variant>.", code="bad_param")


def finemap_abf(trait: Path, out: Path, *, region: str | None = None, lead: str | None = None,
                window_kb: int = 500, p1: float = 1e-4, level: float = 0.95) -> dict[str, Any]:
    from genor.methods.coloc import finemap_abf as fm
    from genor.methods.coloc import subset_region

    with Run("finemap.abf", out, {"trait": trait, "region": region, "lead": lead,
                                  "window_kb": window_kb, "p1": p1, "level": level}) as run:
        df, meta = open_trait(trait)
        run.add_input(Path(trait) / "sumstats.parquet", "sumstats", trait=meta.trait)
        reg = _region_from(df, region, lead, window_kb)
        sub = subset_region(df, reg)
        res, summary = fm(sub, meta, run.warnings, p1=p1, level=level)
        res.to_csv(run.path("finemap.tsv"), sep="\t", index=False)
        run.output("finemap", run.path("finemap.tsv"))
        run.cite("wakefield-abf", "coloc")
        run.warn("single_causal_assumption", "ABF fine-mapping assumes one causal variant in the "
                 "region; multiple independent signals make the credible set misleading.", "info")
        run.summary = {"region": reg, **summary}
    return run.result()


def coloc_run(trait1: Path, trait2: Path, out: Path, *, region: str | None = None,
              lead: str | None = None, window_kb: int = 500, p1: float = 1e-4, p2: float = 1e-4,
              p12: float = 1e-5) -> dict[str, Any]:
    from genor.methods.coloc import coloc_abf, subset_region

    with Run("coloc.abf", out, {"trait1": trait1, "trait2": trait2, "region": region, "lead": lead,
                                "window_kb": window_kb, "p1": p1, "p2": p2, "p12": p12}) as run:
        d1, m1 = open_trait(trait1)
        d2, m2 = open_trait(trait2)
        if m1.build and m2.build and m1.build != m2.build:
            raise GenorError(f"Build mismatch ({m1.build} vs {m2.build}).", code="build_mismatch",
                             hint="genor sumstats liftover one trait first.")
        if not (m1.build and m2.build):
            run.warn("coloc_build_unchecked", "Build unknown for at least one trait; regions are "
                     "matched by position, so make sure both are on the same build.", "serious")
        run.add_input(Path(trait1) / "sumstats.parquet", "trait1", trait=m1.trait)
        run.add_input(Path(trait2) / "sumstats.parquet", "trait2", trait=m2.trait)
        reg = _region_from(d1, region, lead, window_kb)
        from genor.sumstats.harmonise import harmonise
        s1, s2 = subset_region(d1, reg), subset_region(d2, reg)
        h, _ = harmonise(s1, s2, suffixes=("", "_b"), w=run.warnings)
        s2a = pd.DataFrame({c: h[f"{c}_b"] for c in ("beta", "standard_error", "p_value", "n",
                                                    "effect_allele_frequency") if f"{c}_b" in h})
        s2a["variant_id"] = h["variant_id"].to_numpy()
        res, summary = coloc_abf(s1, m1, s2a, m2, run.warnings, p1, p2, p12)
        res.to_csv(run.path("coloc_snp_pp.tsv"), sep="\t", index=False)
        run.output("snp_pp", run.path("coloc_snp_pp.tsv"))
        run.cite("coloc")
        run.warn("single_causal_assumption", "coloc ABF assumes at most one causal variant per "
                 "trait; with multiple signals H3 can be favoured wrongly.", "info")
        run.summary = {"region": reg, "trait1": m1.trait, "trait2": m2.trait, **summary}
    return run.result()


# ============================================================== MAGMA


def magma_gene(trait: Path, out: Path, *, gene_loc: Path, ref_bfile: str,
               window_up_kb: int = 35, window_down_kb: int = 10) -> dict[str, Any]:
    from genor.methods import magma

    with Run("magma.gene", out, {"trait": trait, "gene_loc": gene_loc, "ref_bfile": ref_bfile,
                                 "window": [window_up_kb, window_down_kb]}) as run:
        df, meta = open_trait(trait)
        run.add_input(Path(trait) / "sumstats.parquet", "sumstats", trait=meta.trait)
        run.add_input(gene_loc, "gene_loc")
        p = backends.require("magma")
        run.add_backend("magma", path=p, version=backends._version("magma", p))
        run.cite("magma")
        run.summary = magma.gene_analysis(df, meta, run.path("magma"), run.warnings, gene_loc=gene_loc,
                                          ref_bfile=ref_bfile, window_kb=(window_up_kb, window_down_kb))
    return run.result()


def magma_sets(genes_raw: Path, out: Path, *, set_annot: Path | None = None,
               covar: Path | None = None) -> dict[str, Any]:
    from genor.methods import magma

    if not set_annot and not covar:
        raise GenorError("Give --set-annot (gene sets) or --covar (gene properties).", code="bad_param")
    with Run("magma.sets", out, {"genes_raw": genes_raw, "set_annot": set_annot, "covar": covar}) as run:
        run.add_input(genes_raw, "genes_raw")
        p = backends.require("magma")
        run.add_backend("magma", path=p, version=backends._version("magma", p))
        run.cite("magma")
        run.summary = magma.gene_set(genes_raw, set_annot or Path(), run.path("magma"), covar)
    return run.result()


# ============================================================== reports & misc


def report_build(run_dirs: list[Path], out: Path, title: str = "GeNor analysis report") -> dict[str, Any]:
    from genor.report.builder import build

    return {"status": "ok", **build(run_dirs, out, title)}


def list_runs(root: Path = Path("genor_runs")) -> dict[str, Any]:
    runs = []
    for m in sorted(Path(".").rglob("manifest.json")):
        if any(part.startswith(".") for part in m.parts):
            continue
        try:
            d = json.loads(m.read_text())
        except json.JSONDecodeError:
            continue
        runs.append({"run_dir": str(m.parent), "command": d.get("command"), "status": d.get("status"),
                     "finished": d.get("finished")})
    return {"status": "ok", "runs": runs[-100:]}


def guidance(topic: str | None = None) -> dict[str, Any]:
    from genor import skills

    return skills.guidance(topic)
