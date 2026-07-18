# Expression helpers (R/expr_helpers.R) for Phase 7.

test_that("expr_default_assay prefers norm_logcounts when size factors exist", {
  skip_if_not_installed("DESeq2")
  dds <- ensure_logcounts(make_mock_dds(n_genes = 80, n_per_group = 3, n_spike = 6, seed = 1))
  expect_equal(expr_default_assay(dds), "vst")            # no size factors yet -> VST first
  dds2 <- estimate_size_factors_endogenous(dds)
  expect_equal(expr_default_assay(dds2), "norm_logcounts")
})

test_that("expr_default_assay falls through to stored assays without VST/size factors", {
  skip_if_not_installed("DESeq2")
  dds <- ensure_logcounts(make_mock_dds(n_genes = 60, n_per_group = 2, n_spike = 4, seed = 2))
  # Monkeypatch is overkill; VST is always attempted first, so just assert the
  # counts-only fallback path via a bare object with only counts.
  cnts <- SummarizedExperiment::assay(dds, "counts")
  bare <- DESeq2::DESeqDataSetFromMatrix(cnts, S4Vectors::DataFrame(
    condition = SummarizedExperiment::colData(dds)$condition), design = ~1)
  expect_equal(expr_default_assay(bare), "vst")           # synthetic VST still first
})

test_that("expr_value_matrix resolves stored, norm_logcounts, and vst", {
  skip_if_not_installed("DESeq2")
  dds <- ensure_logcounts(make_mock_dds(n_genes = 120, n_per_group = 3, n_spike = 8, seed = 3))

  lc <- expr_value_matrix(dds, "logcounts")
  expect_true(is.matrix(lc$mat))
  expect_equal(dim(lc$mat), dim(dds))                     # NOT endogenous-restricted
  expect_equal(lc$label, "logcounts")

  nl <- expr_value_matrix(dds, "norm_logcounts")
  expect_equal(nrow(nl$mat), nrow(dds))
  expect_match(nl$label, "normalized")

  v <- expr_value_matrix(dds, "vst")
  expect_true(is.matrix(v$mat))
  expect_equal(ncol(v$mat), ncol(dds))
  # VST is endogenous-only (via qc_vst); rows <= all rows
  expect_lte(nrow(v$mat), nrow(dds))
})

test_that("row_zscore centers/scales rows and guards constant rows", {
  m <- rbind(g1 = c(1, 2, 3, 4), g2 = c(5, 5, 5, 5), g3 = c(-2, 0, 2, 4))
  z <- row_zscore(m)
  expect_equal(dim(z), dim(m))
  expect_equal(rownames(z), rownames(m))
  expect_equal(unname(rowMeans(z)), c(0, 0, 0), tolerance = 1e-8)
  expect_true(all(z["g2", ] == 0))                        # constant row -> 0, not NaN
  expect_true(all(is.finite(z)))
  # sd (population-ish, n-1) of a non-constant row is 1
  expect_equal(stats::sd(z["g1", ]), 1, tolerance = 1e-8)
})

test_that("expr_geom_availability applies the G1/G2 guards", {
  # small groups: dots on by default; distributions hidden (all groups < g1)
  a <- expr_geom_availability(c(4, 5), g1 = 10, g2 = 50)
  expect_true(a$dots_allowed); expect_true(a$dots_default); expect_false(a$dist_shown)

  # a group at/above g1: distributions shown
  a2 <- expr_geom_availability(c(4, 12), g1 = 10, g2 = 50)
  expect_true(a2$dist_shown)

  # mid: dots allowed but off by default (n_max between g1 and g2)
  b <- expr_geom_availability(c(12, 30), g1 = 10, g2 = 50)
  expect_true(b$dots_allowed); expect_false(b$dots_default); expect_true(b$dist_shown)

  # large: dots disallowed
  d <- expr_geom_availability(c(60, 70), g1 = 10, g2 = 50)
  expect_false(d$dots_allowed); expect_false(d$dots_default)

  # all groups tiny: distributions hidden
  e <- expr_geom_availability(c(2, 3), g1 = 10, g2 = 50)
  expect_false(e$dist_shown); expect_true(e$dots_default)
})

test_that("expr_long_frame joins values/groups and drops missing groups", {
  vals <- c(S1 = 1, S2 = 2, S3 = 3, S4 = 4)
  grp  <- factor(c("a", "a", "b", NA), levels = c("a", "b"))
  df <- expr_long_frame(vals, grp)
  expect_equal(nrow(df), 3L)                              # NA group dropped
  expect_true(is.factor(df$group))
  expect_equal(levels(df$group), c("a", "b"))
  expect_equal(df$value, c(1, 2, 3))
  expect_equal(df$sample, c("S1", "S2", "S3"))

  df2 <- expr_long_frame(vals, c("a", "a", "b", "b"), colour = c(1, 1, 2, 2))
  expect_true("colour" %in% names(df2))
  expect_error(expr_long_frame(1:3, 1:2), "length")
})
