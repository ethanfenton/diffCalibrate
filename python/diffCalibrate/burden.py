"""End-to-end sample-size-adjusted DE burden pipeline."""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from typing import Optional, Union

import numpy as np
import pandas as pd

from .adaptive_alpha import apply_adaptive_alpha, calibrate_alpha
from .concordance import concordance_curve
from .downsample import n_degs, run_downsampling_dataset
from .model_fit import fit_deg_model, predict_degs, test_pooling


@dataclass
class BurdenResult:
    """Result of :func:`calculate_burden`.

    Attributes
    ----------
    burden : pd.DataFrame(cell_type, n_full, reference_n, n_degs_observed,
        burden, burden_se, extrapolated) -- the main output.
    downsampling : dict, per-cell-type run_downsampling() results.
    concordance : dict, per-cell-type concordance_curve() results.
    alpha_curves : dict or None, per-cell-type calibrate_alpha() results.
    pooling_test : test_pooling() output, always computed.
    pooling_used : the pooling mode actually used for prediction.
    model : the fitted model object(s) actually used for prediction.
    """

    burden: pd.DataFrame
    downsampling: dict
    concordance: dict
    alpha_curves: Optional[dict]
    pooling_test: dict
    pooling_used: str
    model: object = field(repr=False)

    def __repr__(self) -> str:
        lines = [f"<BurdenResult> {len(self.burden)} cell types, pooling: {self.pooling_used}"]
        lines.append(self.burden.to_string(index=False))
        extrapolated = self.burden.loc[self.burden["extrapolated"], "cell_type"]
        if len(extrapolated):
            lines.append(
                "Note: burden is extrapolated beyond observed cell counts for: "
                + ", ".join(extrapolated)
            )
        return "\n".join(lines)


def calculate_burden(
    data: dict,
    de_fun: Union[str, callable],
    model_type: str = "auto",
    pooling: str = "shared_shape",
    adaptive_alpha: bool = False,
    alpha: float = 0.05,
    reference_n: Optional[float] = None,
    min_concordance: Optional[float] = None,
    concordance_metric: str = "cat_overlap",
    weight_floor: float = 0.1,
    sig_level: float = 0.05,
    grid_args: Optional[dict] = None,
    k: int = 5,
    k_permute: Optional[int] = None,
    seed: Optional[int] = 1,
    **kwargs,
) -> BurdenResult:
    """Calculate sample-size-adjusted DE burden across cell types.

    End-to-end pipeline: downsamples each cell type across a grid of cell
    counts, runs DE at each, computes concordance against the full-data DE
    result, optionally calibrates a per-n adaptive alpha from label-permuted
    null runs, tests whether a spline is needed and whether cell types share
    a common DEGs~n curve, fits the resulting model(s), and predicts the
    expected DEG count ("burden") for every cell type at a common reference
    cell count. See docs/DESIGN.md (repo root) for the full methodology.

    Parameters
    ----------
    data : dict, cell_type -> {"counts": genes x cells matrix, "group": ...}.
    de_fun : DE backend: "wilcoxon", "edger", "deseq2", or a user function
        following the de_wrappers contract.
    model_type : "auto", "spline", or "linear" -- see fit_deg_model().
    pooling : "shared_shape" (default: a single DEGs~n curve *shape* fit
        across all cell types plus a per-cell-type offset, i.e.
        ``degs ~ cr(n) + cell_type``), "auto" (run test_pooling() and use
        its recommendation), "stratified" (fully separate curve per cell
        type), or "pooled" (``degs ~ cr(n)`` with no cell-type term).
        "shared_shape" keeps the stability of pooling every cell type's
        replicates into one curve while still letting cell types differ in
        level. "pooled" has no cell-type term, so it returns the *same*
        burden for every cell type by construction; it is kept only for
        inspecting the shared curve and emits a warning. Under "auto", a
        "pooled" recommendation is predicted with the shared_shape model
        instead, so burden stays cell-type-specific; `pooling_used` records
        the model actually used. test_pooling()'s result is always computed
        and returned as `pooling_test` regardless of this argument.
    adaptive_alpha : if True, calibrate each cell type's raw DE p-values
        against its label-permuted null runs at each n (see
        calibrate_alpha()) before BH adjustment, and count calibrated
        padj < alpha as DEGs for the model. Requires k_permute > 0.
        Default False (fixed alpha on the DE method's own padj).
    alpha : FDR threshold for counting DEGs, applied to the DE method's padj
        or, with adaptive_alpha=True, to the calibrated padj.
    reference_n : cell count at which burden is predicted. Default None
        uses the median, across cell types, of each cell type's full
        per-group cell count.
    min_concordance : downsampled replicates with a concordance_metric value
        below this are excluded from model fitting entirely (None disables
        hard exclusion).
    concordance_metric : which concordance() measure to use for
        min_concordance and for weighting. Default "cat_overlap".
    weight_floor : minimum weight assigned to any retained replicate.
    sig_level : significance threshold for the nonlinearity and pooling tests.
    grid_args, k, k_permute, seed : passed to run_downsampling(). k_permute
        defaults to 5 when adaptive_alpha=True, else 0 (the permuted runs are
        only used for calibration).
    **kwargs : passed to de_fun.

    Returns
    -------
    BurdenResult
    """
    if pooling not in ("pooled", "auto", "shared_shape", "stratified"):
        raise ValueError('pooling must be one of "shared_shape", "auto", "stratified", "pooled"')
    if model_type not in ("auto", "spline", "linear"):
        raise ValueError('model_type must be one of "auto", "spline", "linear"')
    if k_permute is None:
        k_permute = 5 if adaptive_alpha else 0
    if adaptive_alpha and k_permute <= 0:
        raise ValueError("adaptive_alpha=True requires k_permute > 0.")

    cell_types = list(data.keys())

    downsampling = run_downsampling_dataset(
        data, de_fun, grid_args=grid_args or {}, k=k, k_permute=k_permute, alpha=alpha, seed=seed, **kwargs
    )

    concordance = {ct: concordance_curve(downsampling[ct], alpha=alpha) for ct in cell_types}

    alpha_curves = None
    if adaptive_alpha:
        alpha_curves = {ct: calibrate_alpha(downsampling[ct]) for ct in cell_types}
        for ct in cell_types:
            downsampling[ct]["downsampled"] = apply_adaptive_alpha(downsampling[ct], alpha_curves[ct], alpha=alpha)

    combined_rows = []
    for ct in cell_types:
        ds = downsampling[ct]["downsampled"]
        conc = concordance[ct]
        w = np.maximum(conc[concordance_metric].to_numpy(), weight_floor)
        w = np.where(np.isnan(w), weight_floor, w)
        combined_rows.append(
            pd.DataFrame({"cell_type": ct, "n": ds["n"], "degs": ds["n_degs"], "weight": w})
        )
    combined = pd.concat(combined_rows, ignore_index=True)

    if min_concordance is not None:
        conc_all = pd.concat(
            [concordance[ct].assign(cell_type=ct) for ct in cell_types], ignore_index=True
        )
        keep = conc_all[concordance_metric].to_numpy() >= min_concordance
        combined = combined.loc[keep].reset_index(drop=True)

    pooling_test = test_pooling(
        combined["n"], combined["degs"], combined["cell_type"], weights=combined["weight"], sig_level=sig_level
    )
    pooling_used = pooling_test["recommendation"] if pooling == "auto" else pooling
    if pooling == "auto" and pooling_used == "pooled":
        # A model with no cell-type term predicts the same burden for every
        # cell type; use the shared-shape model so burden stays cell-type-specific.
        pooling_used = "shared_shape"
    if pooling_used == "pooled" and len(cell_types) > 1:
        warnings.warn(
            'pooling="pooled" fits no cell-type term, so burden is identical for '
            'every cell type by construction. Use pooling="shared_shape" (default) '
            "to compare cell types."
        )

    if reference_n is None:
        n_full_per_ct = {ct: downsampling[ct]["downsampled"]["n"].max() for ct in cell_types}
        reference_n = float(np.median(list(n_full_per_ct.values())))

    if pooling_used == "pooled":
        model = fit_deg_model(
            combined["n"], combined["degs"], weights=combined["weight"], model_type=model_type, sig_level=sig_level
        )
        pred = predict_degs(model, np.repeat(reference_n, len(cell_types)), se=True)
        burden_rows = pd.DataFrame({"cell_type": cell_types, "burden": pred["fit"], "burden_se": pred["se"]})
    elif pooling_used == "shared_shape":
        model = fit_deg_model(
            combined["n"], combined["degs"], weights=combined["weight"], cell_type=combined["cell_type"],
            model_type=model_type, sig_level=sig_level,
        )
        pred = predict_degs(model, np.repeat(reference_n, len(cell_types)), cell_type_new=cell_types, se=True)
        burden_rows = pd.DataFrame({"cell_type": cell_types, "burden": pred["fit"], "burden_se": pred["se"]})
    else:
        model = {}
        rows = []
        for ct in cell_types:
            sub = combined.loc[combined["cell_type"] == ct]
            model[ct] = fit_deg_model(
                sub["n"], sub["degs"], weights=sub["weight"], model_type=model_type, sig_level=sig_level
            )
            pred = predict_degs(model[ct], reference_n, se=True)
            rows.append({"cell_type": ct, "burden": pred["fit"].iloc[0], "burden_se": pred["se"].iloc[0]})
        burden_rows = pd.DataFrame(rows)

    n_full = {ct: downsampling[ct]["downsampled"]["n"].max() for ct in cell_types}
    n_degs_observed = {ct: n_degs(downsampling[ct]["full"], alpha) for ct in cell_types}

    burden = pd.DataFrame(
        {
            "cell_type": cell_types,
            "n_full": [n_full[ct] for ct in cell_types],
            "reference_n": reference_n,
            "n_degs_observed": [n_degs_observed[ct] for ct in cell_types],
        }
    )
    burden = burden.merge(burden_rows, on="cell_type", how="left")
    burden["extrapolated"] = reference_n > burden["n_full"]

    if burden["extrapolated"].any():
        flagged = ", ".join(burden.loc[burden["extrapolated"], "cell_type"])
        warnings.warn(
            f"reference_n ({reference_n}) exceeds the full cell count for: {flagged}. "
            "Burden for these cell types is extrapolated beyond observed data "
            "(see the `extrapolated` column and treat burden_se accordingly) "
            "-- this is most severe for stratified/per-cell-type models with a "
            'spline, which can diverge sharply outside the training range; '
            'consider model_type="linear" or excluding these cell types.'
        )

    return BurdenResult(
        burden=burden,
        downsampling=downsampling,
        concordance=concordance,
        alpha_curves=alpha_curves,
        pooling_test=pooling_test,
        pooling_used=pooling_used,
        model=model,
    )
