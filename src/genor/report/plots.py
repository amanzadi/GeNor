"""Figures: Manhattan, QQ, MR scatter / forest / funnel / leave-one-out."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

COLORS = ("#3b6ea8", "#8fb3d9")


def _thin(df: pd.DataFrame, keep_p: float = 1e-3, frac: float = 0.1, seed: int = 1) -> pd.DataFrame:
    strong = df[df["p_value"] < keep_p]
    weak = df[df["p_value"] >= keep_p]
    if len(weak) > 200_000:
        weak = weak.sample(frac=frac, random_state=seed)
    return pd.concat([strong, weak])


def manhattan(df: pd.DataFrame, path: Path, title: str = "", sig: float = 5e-8) -> Path:
    d = _thin(df[["chromosome", "base_pair_location", "p_value"]].dropna())
    d = d.assign(chr=pd.to_numeric(d["chromosome"], errors="coerce")).dropna(subset=["chr"])
    d = d.sort_values(["chr", "base_pair_location"])
    offsets, ticks, x = {}, [], np.zeros(len(d))
    off = 0
    for c, g in d.groupby("chr", sort=True):
        offsets[c] = off
        ticks.append((off + g["base_pair_location"].max() / 2, int(c)))
        off += g["base_pair_location"].max() + 5e6
    x = d["base_pair_location"].to_numpy() + d["chr"].map(offsets).to_numpy()
    y = -np.log10(d["p_value"].clip(lower=1e-300).to_numpy())
    col = np.where(d["chr"].to_numpy() % 2 == 0, COLORS[0], COLORS[1])
    fig, ax = plt.subplots(figsize=(12, 4))
    ax.scatter(x, y, c=col, s=2, rasterized=True)
    ax.axhline(-np.log10(sig), color="#b03a2e", lw=0.8, ls="--")
    ax.set_xticks([t[0] for t in ticks], [str(t[1]) for t in ticks], fontsize=7)
    ax.set_xlabel("Chromosome")
    ax.set_ylabel(r"$-\log_{10} P$")
    ax.set_title(title)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def qq(p: pd.Series, path: Path, title: str = "", lam: float | None = None) -> Path:
    p = np.sort(p.dropna().clip(lower=1e-300).to_numpy())
    n = len(p)
    exp = -np.log10((np.arange(1, n + 1) - 0.5) / n)
    obs = -np.log10(p)
    if n > 100_000:  # thin the bulk
        idx = np.unique(np.concatenate([np.arange(0, 10_000), np.linspace(0, n - 1, 50_000).astype(int)]))
        exp, obs = exp[idx], obs[idx]
    fig, ax = plt.subplots(figsize=(4.5, 4.5))
    ax.scatter(exp, obs, s=3, color=COLORS[0], rasterized=True)
    mx = max(exp.max(), 1)
    ax.plot([0, mx], [0, mx], color="#888", lw=0.8)
    ax.set_xlabel(r"Expected $-\log_{10} P$")
    ax.set_ylabel(r"Observed $-\log_{10} P$")
    ax.set_title(title + (f"  (λ = {lam:.3f})" if lam else ""))
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def mr_scatter(h: pd.DataFrame, results: pd.DataFrame, path: Path, exposure: str, outcome: str) -> Path:
    bx, by = h["beta_exp"].to_numpy(), h["beta_out"].to_numpy()
    sx, sy = h["se_exp"].to_numpy(), h["se_out"].to_numpy()
    sgn = np.sign(bx)
    bx, by = bx * sgn, by * sgn
    fig, ax = plt.subplots(figsize=(5.5, 4.5))
    ax.errorbar(bx, by, xerr=sx, yerr=sy, fmt="o", ms=3, color=COLORS[0], ecolor="#bbb", lw=0.6)
    xs = np.array([0, bx.max() * 1.05])
    for _, r in results.iterrows():
        if r["method"] in ("IVW", "MR-Egger", "Weighted median", "Weighted mode"):
            icpt = r.get("intercept", 0.0) if r["method"] == "MR-Egger" else 0.0
            icpt = 0.0 if pd.isna(icpt) else icpt
            ax.plot(xs, icpt + r["b"] * xs, lw=1, label=r["method"])
    ax.axhline(0, color="#888", lw=0.5)
    ax.set_xlabel(f"SNP effect on {exposure}")
    ax.set_ylabel(f"SNP effect on {outcome}")
    ax.legend(fontsize=7, frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def forest(labels: list[str], b: np.ndarray, se: np.ndarray, path: Path, xlabel: str,
           highlight: set[str] | None = None) -> Path:
    n = len(labels)
    fig, ax = plt.subplots(figsize=(5.5, max(2.5, 0.18 * n + 1)))
    y = np.arange(n)[::-1]
    for i, lab in enumerate(labels):
        c = "#b03a2e" if highlight and lab in highlight else COLORS[0]
        ax.errorbar(b[i], y[i], xerr=1.96 * se[i], fmt="o", ms=3, color=c, ecolor=c, lw=0.8)
    ax.axvline(0, color="#888", lw=0.6)
    ax.set_yticks(y, labels, fontsize=6)
    ax.set_xlabel(xlabel)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def funnel(ratio: np.ndarray, ratio_se: np.ndarray, ivw: float, path: Path) -> Path:
    fig, ax = plt.subplots(figsize=(4.5, 4.5))
    ax.scatter(ratio, 1 / ratio_se, s=8, color=COLORS[0])
    ax.axvline(ivw, color="#b03a2e", lw=0.8)
    ax.set_xlabel("Wald ratio estimate")
    ax.set_ylabel("1 / SE")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path
