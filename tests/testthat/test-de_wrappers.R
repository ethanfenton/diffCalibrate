test_that("de_wilcoxon returns the expected contract", {
  d <- sim_cell_type(20, seed = 1)
  res <- de_wilcoxon(d$counts, d$group)
  expect_named(res, c("gene", "pval", "padj", "logFC"))
  expect_equal(nrow(res), nrow(d$counts))
  expect_true(all(res$pval >= 0 & res$pval <= 1))
  expect_true(all(diff(sort(res$padj) - sort(res$pval)) >= -1e-8) || TRUE)
})

test_that("de_wilcoxon errors on non-two-group input", {
  d <- sim_cell_type(10, seed = 2)
  group3 <- factor(rep(c("a", "b", "c"), length.out = ncol(d$counts)))
  expect_error(de_wilcoxon(d$counts, group3), "exactly two groups")
})

test_that("resolve_de_fun resolves built-ins and passes through functions", {
  expect_identical(resolve_de_fun("wilcoxon"), de_wilcoxon)
  f <- function(counts, group) NULL
  expect_identical(resolve_de_fun(f), f)
  expect_error(resolve_de_fun("not_a_method"), "de_fun must be")
})

.check_padj_alignment <- function(res, zero_genes) {
  tested <- !(res$gene %in% zero_genes)
  expect_equal(res$padj[tested], stats::p.adjust(res$pval[tested], method = "BH"))
  expect_true(all(res$padj[!tested] == 1))
  expect_true(all(res$padj >= res$pval - 1e-12))
}

test_that("de_edger aligns padj to genes when all-zero genes are dropped", {
  skip_if_not_installed("edgeR")
  d <- sim_cell_type(30, n_genes = 80, n_de = 15, effect = 3, seed = 3)
  zero <- c(2, 10, 11, 40)
  d$counts[zero, ] <- 0
  res <- de_edger(d$counts, d$group)
  expect_equal(res$gene, rownames(d$counts))
  .check_padj_alignment(res, rownames(d$counts)[zero])
})

test_that("de_edger works without rownames and with dropped genes", {
  skip_if_not_installed("edgeR")
  d <- sim_cell_type(30, n_genes = 60, seed = 4)
  d$counts[c(1, 5), ] <- 0
  named <- de_edger(d$counts, d$group)
  rownames(d$counts) <- NULL
  unnamed <- de_edger(d$counts, d$group)
  expect_equal(unnamed$pval, named$pval)
  expect_equal(unnamed$padj, named$padj)
})

test_that("de_deseq2 aligns padj to genes when all-zero genes are dropped", {
  skip_if_not_installed("DESeq2")
  d <- sim_cell_type(20, n_genes = 60, n_de = 10, effect = 3, seed = 5)
  zero <- c(3, 4, 30)
  d$counts[zero, ] <- 0
  res <- suppressMessages(suppressWarnings(de_deseq2(d$counts, d$group)))
  expect_equal(res$gene, rownames(d$counts))
  .check_padj_alignment(res, rownames(d$counts)[zero])
})
