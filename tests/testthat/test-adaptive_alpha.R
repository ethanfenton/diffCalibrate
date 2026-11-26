test_that("calibrate_alpha returns one row per grid point and errors without permuted arm", {
  d <- sim_cell_type(40, n_genes = 40, seed = 8)
  out <- run_downsampling(d$counts, d$group, "wilcoxon",
                           grid_args = list(n_min = 10, n_points = 3),
                           k = 2, k_permute = 3, seed = 12)
  ac <- calibrate_alpha(out, target_fpr = 0.05)
  expect_equal(nrow(ac), length(unique(out$downsampled$n)))
  expect_true(all(ac$null_quantile >= 0 & ac$null_quantile <= 1))
  expect_equal(ac$n_null_tests, lengths(ac$null_p))
  expect_true(all(ac$n_null_tests == 40 * 3))

  out_noperm <- run_downsampling(d$counts, d$group, "wilcoxon",
                                  grid_args = list(n_min = 10, n_points = 3),
                                  k = 2, k_permute = 0, seed = 13)
  expect_error(calibrate_alpha(out_noperm), "no permuted arm")
})

test_that("calibration uses raw p-values, not BH-adjusted ones", {
  # Regression: <= 0.2.0 took quantiles of null padj (all near 1), giving
  # thresholds ~0.75-0.95 that called nearly every gene a DEG.
  d <- sim_cell_type(40, n_genes = 80, n_de = 5, seed = 20)
  out <- run_downsampling(d$counts, d$group, "wilcoxon",
                           grid_args = list(n_min = 10, n_points = 3),
                           k = 2, k_permute = 5, seed = 21)
  ac <- calibrate_alpha(out)
  expect_true(all(ac$null_quantile < 0.3))
  ds2 <- apply_adaptive_alpha(out, ac, alpha = 0.05)
  expect_true(all(ds2$n_degs < 0.5 * nrow(d$counts)))
})

test_that(".calibrate_p is identity-like for uniform nulls and corrects inflated ones", {
  set.seed(1)
  p <- c(1e-6, 1e-3, 0.01, 0.05, 0.2, 0.5, 0.9)
  unif <- stats::runif(1e5)
  expect_equal(.calibrate_p(p, unif), p, tolerance = 0.05)

  inflated <- unif^2  # anti-conservative null: too many small p-values
  p_cal <- .calibrate_p(p, inflated)
  expect_true(all(p_cal >= p))
  expect_false(is.unsorted(p_cal))
  expect_true(all(p_cal > 0 & p_cal <= 1))
})

test_that("cap_at_target never makes p-values smaller", {
  set.seed(2)
  conservative <- sqrt(stats::runif(1e4))
  p <- c(0.001, 0.01, 0.1, 0.5)
  expect_true(any(.calibrate_p(p, conservative) < p))
  expect_true(all(.calibrate_p(p, conservative, cap_at_target = TRUE) >= p))
})

test_that("predict_alpha interpolates null_quantile and holds endpoints constant", {
  ac <- data.frame(n = c(10, 50, 100), null_quantile = c(0.01, 0.05, 0.1))
  expect_equal(predict_alpha(ac, 50), 0.05)
  expect_equal(predict_alpha(ac, 5), 0.01)
  expect_equal(predict_alpha(ac, 200), 0.1)
})

test_that("apply_adaptive_alpha recomputes n_degs without rerunning DE", {
  d <- sim_cell_type(40, n_genes = 40, seed = 10)
  out <- run_downsampling(d$counts, d$group, "wilcoxon",
                           grid_args = list(n_min = 10, n_points = 3),
                           k = 2, k_permute = 3, seed = 15)
  ac <- calibrate_alpha(out)
  ds2 <- apply_adaptive_alpha(out, ac, alpha = 0.05)
  expect_true(all(c("alpha_used", "null_quantile") %in% names(ds2)))
  expect_true(all(ds2$alpha_used == 0.05))
  expect_equal(nrow(ds2), nrow(out$downsampled))
})
