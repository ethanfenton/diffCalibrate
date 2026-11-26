import numpy as np
import pandas as pd
import pytest

from diffCalibrate.adaptive_alpha import (
    _calibrate_p,
    apply_adaptive_alpha,
    calibrate_alpha,
    predict_alpha,
)
from diffCalibrate.downsample import run_downsampling


def _ds(sim, seed, k_permute=3, **kw):
    d = sim(seed=seed, **kw)
    return d, run_downsampling(
        d["counts"], d["group"], "wilcoxon", grid_args={"n_points": 3}, k=2, k_permute=k_permute, seed=1
    )


def test_calibrate_alpha_requires_permuted_arm(sim):
    _, res = _ds(sim, 1, k_permute=0, n_genes=40, n_per_group=30)
    with pytest.raises(ValueError):
        calibrate_alpha(res)


def test_calibrate_alpha_shape_and_bounds(sim):
    _, res = _ds(sim, 2, n_genes=60, n_per_group=40)
    curve = calibrate_alpha(res, target_fpr=0.05)
    assert set(curve.columns) == {"n", "null_quantile", "n_null_tests", "null_p"}
    assert curve["null_quantile"].between(0, 1).all()
    assert (curve["n_null_tests"] == 60 * 3).all()
    assert list(curve["n"]) == sorted(curve["n"])


def test_calibration_uses_raw_pvalues_not_padj(sim):
    # Regression: <= 0.2.0 took quantiles of null padj (all near 1), giving
    # thresholds ~0.75-0.95 that called nearly every gene a DEG.
    d, res = _ds(sim, 3, k_permute=5, n_genes=100, n_per_group=40, n_de=5)
    curve = calibrate_alpha(res)
    assert (curve["null_quantile"] < 0.3).all()
    out = apply_adaptive_alpha(res, curve, alpha=0.05)
    assert (out["n_degs"] < 50).all()


def test_calibrate_p_uniform_null_is_identity_like():
    rng = np.random.default_rng(1)
    p = np.array([1e-6, 1e-3, 0.01, 0.05, 0.2, 0.5, 0.9])
    np.testing.assert_allclose(_calibrate_p(p, rng.uniform(size=100_000)), p, rtol=0.1, atol=0.01)


def test_calibrate_p_corrects_inflated_null():
    rng = np.random.default_rng(2)
    p = np.array([1e-6, 1e-3, 0.01, 0.05, 0.2, 0.5, 0.9])
    p_cal = _calibrate_p(p, rng.uniform(size=100_000) ** 2)
    assert (p_cal >= p).all()
    assert (np.diff(p_cal) >= 0).all()
    assert ((p_cal > 0) & (p_cal <= 1)).all()


def test_cap_at_target_never_shrinks_pvalues():
    rng = np.random.default_rng(3)
    conservative = np.sqrt(rng.uniform(size=10_000))
    p = np.array([0.001, 0.01, 0.1, 0.5])
    assert (_calibrate_p(p, conservative) < p).any()
    assert (_calibrate_p(p, conservative, cap_at_target=True) >= p).all()


def test_predict_alpha_interpolates_and_clamps():
    curve = pd.DataFrame({"n": [20, 50, 100], "null_quantile": [0.01, 0.03, 0.05]})
    out = predict_alpha(curve, [20, 35, 100, 500, 1])
    assert out[0] == pytest.approx(0.01)
    assert out[1] == pytest.approx(0.02)
    assert out[3] == pytest.approx(0.05)
    assert out[4] == pytest.approx(0.01)


def test_apply_adaptive_alpha_adds_columns(sim):
    _, res = _ds(sim, 4, n_genes=60, n_per_group=40)
    curve = calibrate_alpha(res)
    out = apply_adaptive_alpha(res, curve, alpha=0.05)
    assert {"alpha_used", "null_quantile"} <= set(out.columns)
    assert len(out) == len(res["downsampled"])
