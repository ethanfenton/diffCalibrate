#' Empirical null calibration of DE p-values as a function of n
#'
#' At small cell counts, DE p-values can be miscalibrated (typically
#' anti-conservative), so a fixed `padj < 0.05` cutoff does not carry the
#' same false-discovery rate at every sample size. This calibrates *raw*
#' p-values against the label-permuted arm of [run_downsampling()]
#' (`k_permute > 0`): permuted labels give DE runs with no true group
#' difference, so their raw p-values form an empirical null at each `n`.
#'
#' [apply_adaptive_alpha()] then maps each observed raw p-value to a
#' calibrated p-value (the fraction of null p-values at or below it; see
#' Details), applies BH to the calibrated p-values, and counts genes with
#' calibrated `padj < alpha`. The FDR meaning of `alpha` is therefore kept
#' at every `n`; only the p-values feeding BH are corrected.
#'
#' Calibration is done on raw p-values, not BH-adjusted ones: under the
#' null nearly all adjusted p-values sit near 1, so quantiles of null
#' `padj` are not a usable significance threshold (versions <= 0.2.0 did
#' this and produced thresholds around 0.75-0.95).
#'
#' @details For an observed raw p-value `p` and null p-values with
#'   empirical CDF `F0`, the calibrated p-value is `F0(p)` when `p` is at or
#'   above `q`, the `target_fpr` quantile of the null, and `p * target_fpr /
#'   q` below it. The linear tail avoids the resolution floor of an
#'   empirical CDF (otherwise every very small p would tie at `1 / N`) and
#'   is continuous and monotone with `F0` at `q`. For well-calibrated
#'   p-values `q ~ target_fpr`, so calibrated p ~ raw p.
#'
#' @param downsampling_result A single cell type's result from
#'   [run_downsampling()], run with `k_permute > 0`.
#' @param target_fpr Null quantile at which the empirical CDF hands over to
#'   the linear tail (see Details). Default 0.05.
#' @param cap_at_target If `TRUE`, calibration may only make p-values
#'   larger (correct anti-conservative bias, never spend extra power on
#'   conservative bias). Default `FALSE`.
#' @return `data.frame(n, null_quantile, n_null_tests, null_p)`, one row per
#'   grid point. `null_quantile` is the `target_fpr` quantile of the null
#'   raw p-values (below `target_fpr` means anti-conservative at that `n`);
#'   `n_null_tests` is the number of null (gene x permutation replicate)
#'   p-values; `null_p` is a list column holding them. `target_fpr` and
#'   `cap_at_target` are stored as attributes.
#' @export
calibrate_alpha <- function(downsampling_result, target_fpr = 0.05,
                             cap_at_target = FALSE) {
  perm <- downsampling_result$permuted
  if (is.null(perm)) {
    stop("downsampling_result has no permuted arm; rerun run_downsampling() ",
         "with k_permute > 0.", call. = FALSE)
  }

  ns <- sort(unique(perm$n))
  null_p <- lapply(ns, function(n) {
    rows <- perm[perm$n == n, , drop = FALSE]
    p <- unlist(lapply(rows$de_res, function(r) r$pval))
    sort(p[!is.na(p)])
  })
  null_quantile <- vapply(null_p, function(p) {
    if (length(p) == 0) NA_real_ else
      as.numeric(stats::quantile(p, probs = target_fpr, type = 1))
  }, numeric(1))

  out <- data.frame(n = ns, null_quantile = null_quantile,
                    n_null_tests = lengths(null_p))
  out$null_p <- null_p
  attr(out, "target_fpr") <- target_fpr
  attr(out, "cap_at_target") <- cap_at_target
  out
}

#' Calibrate raw p-values against an empirical null
#'
#' @param p Raw p-values.
#' @param null_p Null raw p-values (sorted or not).
#' @param target_fpr,cap_at_target See [calibrate_alpha()].
#' @return Calibrated p-values, same length as `p`.
#' @keywords internal
.calibrate_p <- function(p, null_p, target_fpr = 0.05, cap_at_target = FALSE) {
  null_p <- sort(null_p[!is.na(null_p)])
  if (length(null_p) == 0) return(p)
  q <- as.numeric(stats::quantile(null_p, probs = target_fpr, type = 1))
  f0 <- findInterval(p, null_p, rightmost.closed = FALSE) / length(null_p)
  p_cal <- if (q > 0) ifelse(p < q, p * target_fpr / q, f0) else f0
  p_cal <- pmin(p_cal, 1)
  if (cap_at_target) p_cal <- pmax(p_cal, p)
  p_cal
}

#' Interpolate the calibration quantile to arbitrary cell counts
#'
#' Linear interpolation of `null_quantile` from [calibrate_alpha()], with
#' endpoints held constant outside the observed range. Diagnostic only:
#' shows how anti-conservative (`< target_fpr`) or conservative the DE
#' method is at a given `n`. It is not a significance threshold.
#'
#' @param alpha_curve Output of [calibrate_alpha()].
#' @param n Vector of cell counts to evaluate at.
#' @return Numeric vector, same length as `n`.
#' @export
predict_alpha <- function(alpha_curve, n) {
  stats::approx(alpha_curve$n, alpha_curve$null_quantile, xout = n, rule = 2)$y
}

#' Recompute DEG counts from null-calibrated p-values
#'
#' Reuses the DE result tables already stored in a [run_downsampling()]
#' result (`downsampled$de_res`), so no DE refitting is needed. For each
#' replicate, raw p-values are calibrated against the null p-values from
#' the nearest grid `n` in `alpha_curve` (see [calibrate_alpha()]),
#' BH-adjusted, and genes with calibrated `padj < alpha` are counted.
#'
#' @param downsampling_result A single cell type's result from
#'   [run_downsampling()].
#' @param alpha_curve Output of [calibrate_alpha()] for the same cell type.
#' @param alpha FDR threshold applied to the calibrated, BH-adjusted
#'   p-values. Default 0.05.
#' @return The `downsampled` `data.frame` from `downsampling_result`, with
#'   `n_degs` replaced by calibrated counts, plus columns `alpha_used`
#'   (= `alpha`) and `null_quantile` (the calibration quantile used).
#' @export
apply_adaptive_alpha <- function(downsampling_result, alpha_curve, alpha = 0.05) {
  ds <- downsampling_result$downsampled
  target_fpr <- attr(alpha_curve, "target_fpr") %||% 0.05
  cap <- attr(alpha_curve, "cap_at_target") %||% FALSE
  idx <- vapply(ds$n, function(n) which.min(abs(alpha_curve$n - n)), integer(1))

  ds$n_degs <- vapply(seq_len(nrow(ds)), function(i) {
    p_cal <- .calibrate_p(ds$de_res[[i]]$pval, alpha_curve$null_p[[idx[i]]],
                          target_fpr = target_fpr, cap_at_target = cap)
    sum(stats::p.adjust(p_cal, method = "BH") < alpha, na.rm = TRUE)
  }, numeric(1))
  ds$alpha_used <- alpha
  ds$null_quantile <- alpha_curve$null_quantile[idx]
  ds
}
