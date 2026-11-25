"""Two-sample Mendelian randomisation from summary statistics.

Estimators follow TwoSampleMR 0.7.x exactly where the method is deterministic
(IVW, MR-Egger, Wald ratio, weighted median point estimate, mode point estimates
via a port of R's ``density()``); bootstrap standard errors match in distribution.
See tests/golden/ for the reference values.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
from scipy import stats

# --------------------------------------------------------------------------- helpers


def _wls(y: np.ndarray, X: np.ndarray, w: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    """Weighted least squares like R's lm(weights=w): coef, se, residual sigma."""
    sw = np.sqrt(w)
    Xw, yw = X * sw[:, None], y * sw
    coef, *_ = np.linalg.lstsq(Xw, yw, rcond=None)
    resid = yw - Xw @ coef
    dof = len(y) - X.shape[1]
    sigma = float(np.sqrt(resid @ resid / dof)) if dof > 0 else float("nan")
    cov = np.linalg.inv(Xw.T @ Xw) * sigma**2
    return coef, np.sqrt(np.diag(cov)), sigma


def _p_norm(b: float, se: float) -> float:
    return float(2 * stats.norm.sf(abs(b / se)))


def r_mad(x: np.ndarray) -> float:
    return float(1.4826 * np.median(np.abs(x - np.median(x))))


def r_density_mode(x: np.ndarray, weights: np.ndarray, bw: float, n: int = 512) -> float:
    """x at the maximum of R's stats::density(x, weights, bw) (gaussian kernel, R >= 4.4)."""
    cut, ext = 3, 4
    frm, to = x.min() - cut * bw, x.max() + cut * bw
    lo, up = frm - ext * bw, to + ext * bw
    # C_BinDist: linear binning into 2n cells
    y = np.zeros(2 * n)
    xdelta = (up - lo) / (n - 1)
    xpos = (x - lo) / xdelta
    ix = np.floor(xpos).astype(int)
    fx = xpos - ix
    for i in range(len(x)):
        k, f, wi = ix[i], fx[i], weights[i]
        if 0 <= k <= n - 2:
            y[k] += wi * (1 - f)
            y[k + 1] += wi * f
        elif k == -1:
            y[0] += wi * f
        elif k == n - 1:
            y[k] += wi * (1 - f)
    kords = np.linspace(0, (2 * n - 1) / (n - 1) * (up - lo), 2 * n)
    kords[n + 1:] = -kords[n - 1:0:-1]
    kords = stats.norm.pdf(kords, scale=bw)
    dens = np.fft.ifft(np.fft.fft(y) * np.conj(np.fft.fft(kords))).real[:n]
    dens = np.maximum(0, dens)  # R: Re()/length(y) then pmax; ifft already divides by 2n
    xords = np.linspace(lo, up, n)
    xs = np.linspace(frm, to, n)
    ys = np.interp(xs, xords, dens)
    return float(xs[int(np.argmax(ys))])


def _mode_beta(beta_iv: np.ndarray, se_iv: np.ndarray, phi: float = 1.0) -> float:
    s = 0.9 * min(np.std(beta_iv, ddof=1), r_mad(beta_iv)) / len(beta_iv) ** 0.2
    w = se_iv**-2 / np.sum(se_iv**-2)
    h = max(1e-8, s * phi)
    return r_density_mode(beta_iv, w, h)


def weighted_median(b_iv: np.ndarray, weights: np.ndarray) -> float:
    o = np.argsort(b_iv, kind="stable")
    b, w = b_iv[o], weights[o]
    cs = (np.cumsum(w) - 0.5 * w) / np.sum(w)
    below = np.max(np.where(cs < 0.5)[0])
    return float(b[below] + (b[below + 1] - b[below]) * (0.5 - cs[below]) / (cs[below + 1] - cs[below]))


# --------------------------------------------------------------------------- estimators


def wald_ratio(bx: float, by: float, sx: float, sy: float) -> dict[str, float]:
    b = by / bx
    se = sy / abs(bx)
    return {"b": b, "se": se, "pval": _p_norm(b, se), "nsnp": 1}


def ivw(bx: np.ndarray, by: np.ndarray, sx: np.ndarray, sy: np.ndarray,
        fixed: bool = False) -> dict[str, float]:
    coef, se, sigma = _wls(by, bx[:, None], 1 / sy**2)
    b = float(coef[0])
    se_b = float(se[0]) / (sigma if fixed else min(1.0, sigma))
    q_df = len(bx) - 1
    Q = sigma**2 * q_df
    return {"b": b, "se": se_b, "pval": _p_norm(b, se_b), "nsnp": len(bx),
            "Q": Q, "Q_df": q_df, "Q_pval": float(stats.chi2.sf(Q, q_df))}


def egger(bx: np.ndarray, by: np.ndarray, sx: np.ndarray, sy: np.ndarray) -> dict[str, float]:
    sign = np.where(bx == 0, 1.0, np.sign(bx))
    y, x = by * sign, np.abs(bx)
    X = np.column_stack([np.ones_like(x), x])
    coef, se, sigma = _wls(y, X, 1 / sy**2)
    k = len(bx)
    f = min(1.0, sigma)
    b, se_b = float(coef[1]), float(se[1]) / f
    bi, se_i = float(coef[0]), float(se[0]) / f
    Q = sigma**2 * (k - 2)
    return {"b": b, "se": se_b, "pval": float(2 * stats.t.sf(abs(b / se_b), k - 2)),
            "intercept": bi, "intercept_se": se_i,
            "intercept_pval": float(2 * stats.t.sf(abs(bi / se_i), k - 2)),
            "nsnp": k, "Q": Q, "Q_df": k - 2, "Q_pval": float(stats.chi2.sf(Q, k - 2)),
            "I2_GX": i2_gx(np.abs(bx), sx)}


def i2_gx(bx_abs: np.ndarray, sx: np.ndarray) -> float:
    """Bowden's I^2_GX: NOME violation for MR-Egger (want > 0.9)."""
    w = 1 / sx**2
    mu = np.sum(bx_abs * w) / np.sum(w)
    Q = np.sum(w * (bx_abs - mu) ** 2)
    k = len(bx_abs)
    return float(max(0.0, (Q - (k - 1)) / Q)) if Q > 0 else 0.0


def weighted_median_est(bx, by, sx, sy, nboot: int = 1000, seed: int = 1) -> dict[str, float]:
    b_iv = by / bx
    vbj = sy**2 / bx**2 + by**2 * sx**2 / bx**4
    b = weighted_median(b_iv, 1 / vbj)
    rng = np.random.default_rng(seed)
    med = np.empty(nboot)
    for i in range(nboot):
        bxb = rng.normal(bx, sx)
        byb = rng.normal(by, sy)
        med[i] = weighted_median(byb / bxb, 1 / vbj)
    se = float(np.std(med, ddof=1))
    return {"b": b, "se": se, "pval": _p_norm(b, se), "nsnp": len(bx)}


def mode_est(bx, by, sx, sy, weighted: bool = True, nboot: int = 1000, seed: int = 1,
             phi: float = 1.0) -> dict[str, float]:
    b_iv = by / bx
    se_iv = np.sqrt(sy**2 / bx**2 + by**2 * sx**2 / bx**4)
    use_se = se_iv if weighted else np.ones_like(se_iv)
    b = _mode_beta(b_iv, use_se, phi)
    rng = np.random.default_rng(seed)
    boots = np.empty(nboot)
    for i in range(nboot):
        bb = rng.normal(b_iv, se_iv)
        boots[i] = _mode_beta(bb, use_se, phi)
    se = r_mad(boots)
    k = len(bx)
    return {"b": b, "se": se, "pval": float(2 * stats.t.sf(abs(b / se), k - 1)), "nsnp": k}


def steiger(bx, sx, nx, by, sy, ny) -> dict[str, Any]:
    """Steiger directionality using r from b, se, N for both traits (TwoSampleMR get_r_from_bsen)."""
    def r2(b, se, n):
        F = (b / se) ** 2
        return F / (n - 2 + F)
    rx = float(np.sqrt(np.sum(r2(bx, sx, nx))))
    ry = float(np.sqrt(np.sum(r2(by, sy, ny))))
    nx_, ny_ = float(np.median(nx)), float(np.median(ny))
    z = (np.arctanh(rx) - np.arctanh(ry)) / np.sqrt(1 / (nx_ - 3) + 1 / (ny_ - 3))
    return {"r2_exposure": rx**2, "r2_outcome": ry**2, "correct_causal_direction": bool(rx > ry),
            "steiger_pval": float(2 * stats.norm.sf(abs(z)))}


def f_stats(bx: np.ndarray, sx: np.ndarray) -> dict[str, float]:
    F = (bx / sx) ** 2
    return {"F_mean": float(F.mean()), "F_min": float(F.min()), "n_weak_F_lt_10": int((F < 10).sum())}


def run_all(h: pd.DataFrame, nboot: int = 1000, seed: int = 1) -> dict[str, Any]:
    """h: columns SNP, beta_exp, se_exp, beta_out, se_out, n_exp, n_out."""
    bx, by = h["beta_exp"].to_numpy(float), h["beta_out"].to_numpy(float)
    sx, sy = h["se_exp"].to_numpy(float), h["se_out"].to_numpy(float)
    k = len(h)
    rows: list[dict[str, Any]] = []
    out: dict[str, Any] = {"nsnp": k}
    if k == 0:
        return {"nsnp": 0, "results": pd.DataFrame(), "single": pd.DataFrame(), "loo": pd.DataFrame()}
    if k == 1:
        rows.append({"method": "Wald ratio", **wald_ratio(bx[0], by[0], sx[0], sy[0])})
    if k >= 2:
        iv = ivw(bx, by, sx, sy)
        rows.append({"method": "IVW", **iv})
        rows.append({"method": "IVW (fixed effects)", **ivw(bx, by, sx, sy, fixed=True)})
        out["heterogeneity_ivw"] = {"Q": iv["Q"], "Q_df": iv["Q_df"], "Q_pval": iv["Q_pval"]}
    if k >= 3:
        eg = egger(bx, by, sx, sy)
        rows.append({"method": "MR-Egger", **{kk: v for kk, v in eg.items()}})
        rows.append({"method": "Weighted median", **weighted_median_est(bx, by, sx, sy, nboot, seed)})
        rows.append({"method": "Weighted mode", **mode_est(bx, by, sx, sy, True, nboot, seed)})
        rows.append({"method": "Simple mode", **mode_est(bx, by, sx, sy, False, nboot, seed)})
        out["pleiotropy_egger"] = {"intercept": eg["intercept"], "se": eg["intercept_se"],
                                   "pval": eg["intercept_pval"]}
        out["heterogeneity_egger"] = {"Q": eg["Q"], "Q_df": eg["Q_df"], "Q_pval": eg["Q_pval"]}
        out["I2_GX"] = eg["I2_GX"]
    res = pd.DataFrame(rows)
    single = pd.DataFrame([{"SNP": s, **wald_ratio(a, b, c, d)}
                           for s, a, b, c, d in zip(h["SNP"], bx, by, sx, sy, strict=True)])
    loo = pd.DataFrame()
    if k >= 3:
        loo = pd.DataFrame([{"SNP_left_out": h["SNP"].iloc[i],
                             **ivw(np.delete(bx, i), np.delete(by, i), np.delete(sx, i), np.delete(sy, i))}
                            for i in range(k)])
    out.update({"results": res, "single": single, "loo": loo, "F": f_stats(bx, sx)})
    if "n_exp" in h and "n_out" in h and h["n_exp"].notna().all() and h["n_out"].notna().all():
        out["steiger"] = steiger(bx, sx, h["n_exp"].to_numpy(float), by, sy, h["n_out"].to_numpy(float))
    return out
