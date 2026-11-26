"""Empirical null calibration of DE p-values as a function of n.

At small cell counts, DE p-values can be miscalibrated (typically
anti-conservative), so a fixed padj < 0.05 cutoff does not carry the same
false-discovery rate at every sample size. This calibrates *raw* p-values
against the label-permuted arm of run_downsampling (k_permute > 0):
permuted labels give DE runs with no true group difference, so their raw
p-values form an empirical null at each n.

apply_adaptive_alpha() maps each observed raw p-value to a calibrated
p-value (the fraction of null p-values at or below it, with a linear tail
below the target_fpr null quantile), applies BH to the calibrated p-values,
and counts genes with calibrated padj < alpha. The FDR meaning of alpha is
kept at every n; only the p-values feeding BH are corrected.

Calibration uses raw p-values, not BH-adjusted ones: under the null nearly
all adjusted p-values sit near 1, so quantiles of null padj are not a
usable threshold (versions <= 0.2.0 did this and produced thresholds around
0.75-0.95). See docs/DESIGN.md (repo root) for details.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .de_wrappers import _bh_adjust


def calibrate_alpha(
    downsampling_result: dict, target_fpr: float = 0.05, cap_at_target: bool = False
) -> pd.DataFrame:
    """Collect the empirical null of raw p-values at each grid n.

    Returns DataFrame(n, null_quantile, n_null_tests, null_p). null_quantile
    is the target_fpr quantile of the null raw p-values (below target_fpr
    means the DE method is anti-conservative at that n); null_p holds the
    sorted null p-values. target_fpr and cap_at_target are stored in
    ``.attrs``.
    """
    perm = downsampling_result["permuted"]
    if perm is None:
        raise ValueError(
            "downsampling_result has no permuted arm; rerun run_downsampling() "
            "with k_permute > 0."
        )

    rows = []
    for n, group in perm.groupby("n"):
        null_p = np.concatenate([r["pval"].to_numpy(dtype=float) for r in group["de_res"]])
        null_p = np.sort(null_p[~np.isnan(null_p)])
        q = float(np.quantile(null_p, target_fpr, method="inverted_cdf")) if len(null_p) else np.nan
        rows.append({"n": int(n), "null_quantile": q, "n_null_tests": len(null_p), "null_p": null_p})
    out = pd.DataFrame(rows).sort_values("n").reset_index(drop=True)
    out.attrs["target_fpr"] = target_fpr
    out.attrs["cap_at_target"] = cap_at_target
    return out


def _calibrate_p(p, null_p, target_fpr: float = 0.05, cap_at_target: bool = False) -> np.ndarray:
    """Map raw p-values to calibrated p-values against an empirical null."""
    p = np.asarray(p, dtype=float)
    null_p = np.sort(np.asarray(null_p, dtype=float))
    null_p = null_p[~np.isnan(null_p)]
    if len(null_p) == 0:
        return p.copy()
    q = float(np.quantile(null_p, target_fpr, method="inverted_cdf"))
    f0 = np.searchsorted(null_p, p, side="right") / len(null_p)
    if q > 0:
        p_cal = np.where(p < q, p * target_fpr / q, f0)
    else:
        p_cal = f0
    p_cal = np.minimum(p_cal, 1.0)
    p_cal = np.where(np.isnan(p), np.nan, p_cal)
    if cap_at_target:
        p_cal = np.maximum(p_cal, p)
    return p_cal


def predict_alpha(alpha_curve: pd.DataFrame, n) -> np.ndarray:
    """Linearly interpolate null_quantile to arbitrary cell counts.

    Endpoints are held constant outside the observed range. Diagnostic only
    (how anti-conservative the DE method is at n); not a significance
    threshold.
    """
    n = np.atleast_1d(np.asarray(n, dtype=float))
    return np.interp(n, alpha_curve["n"].to_numpy(), alpha_curve["null_quantile"].to_numpy())


def apply_adaptive_alpha(
    downsampling_result: dict, alpha_curve: pd.DataFrame, alpha: float = 0.05
) -> pd.DataFrame:
    """Recompute DEG counts from null-calibrated p-values.

    Reuses the DE result tables already stored in downsampling_result, so no
    DE refitting is needed. For each replicate, raw p-values are calibrated
    against the null from the nearest grid n, BH-adjusted, and genes with
    calibrated padj < alpha are counted. Returns the `downsampled`
    DataFrame with n_degs replaced and new columns `alpha_used` (= alpha)
    and `null_quantile`.
    """
    ds = downsampling_result["downsampled"].copy()
    target_fpr = alpha_curve.attrs.get("target_fpr", 0.05)
    cap = alpha_curve.attrs.get("cap_at_target", False)
    grid_n = alpha_curve["n"].to_numpy(dtype=float)
    idx = [int(np.argmin(np.abs(grid_n - n))) for n in ds["n"].to_numpy(dtype=float)]

    counts = []
    for res, i in zip(ds["de_res"], idx):
        p_cal = _calibrate_p(res["pval"].to_numpy(dtype=float), alpha_curve["null_p"].iloc[i], target_fpr, cap)
        padj = _bh_adjust(p_cal)
        counts.append(int(np.sum(padj[~np.isnan(padj)] < alpha)))
    ds["n_degs"] = counts
    ds["alpha_used"] = alpha
    ds["null_quantile"] = alpha_curve["null_quantile"].to_numpy()[idx]
    return ds
