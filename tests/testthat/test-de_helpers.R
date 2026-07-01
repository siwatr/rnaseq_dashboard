# DE engine (R/de_helpers.R) for Phase 5a.

# --- classification --------------------------------------------------------

test_that("classify applies the padj + |LFC| rule and fixed DEG levels", {
  lfc  <- c(3, -3, 0.5, 3,   NA)
  padj <- c(0.01, 0.01, 0.001, 0.2, 0.001)
  cl <- classify(lfc, padj, padj_cut = 0.05, lfc_cut = log2(2))
  expect_equal(cl$sig, c(TRUE, TRUE, FALSE, FALSE, FALSE))     # #3 small LFC, #4 padj, #5 NA lfc
  expect_equal(as.character(cl$DEG), c("up", "down", "no_change", "no_change", "no_change"))
  expect_equal(levels(cl$DEG), c("up", "down", "no_change"))
})

test_that("classify treats NA padj as not significant", {
  cl <- classify(c(5, 5), c(NA, 0.001), 0.05, log2(2))
  expect_equal(cl$sig, c(FALSE, TRUE))
})

test_that("de_classify_table adds both LFC variants when the shrunk column exists", {
  df <- data.frame(
    log2FoldChange = c(3, -3, 0.2),
    log2FoldChange_shrunk = c(1.5, -0.4, 0.2),
    padj = c(0.001, 0.001, 0.001))
  out <- de_classify_table(df, 0.05, log2(2))
  expect_true(all(c("sig", "DEG", "sig_shrunk", "DEG_shrunk") %in% names(out)))
  expect_equal(as.character(out$DEG), c("up", "down", "no_change"))
  # shrinkage pulls #2 below the LFC cut -> no_change under the shrunk classification
  expect_equal(as.character(out$DEG_shrunk), c("up", "no_change", "no_change"))
})

test_that("de_summary counts up/down/total", {
  df <- data.frame(DEG = factor(c("up", "up", "down", "no_change"),
                                levels = c("up", "down", "no_change")))
  s <- de_summary(df)
  expect_equal(unname(s[c("up", "down", "total")]), c(2, 1, 3))
})

# --- design + contrast -----------------------------------------------------

test_that("de_full_rank flags a confounded design", {
  cd <- data.frame(
    condition = factor(rep(c("ctrl", "trt"), each = 2)),
    batch     = factor(rep(c("b1", "b2"), each = 2)))   # batch == condition -> confounded
  ok  <- de_full_rank(~ condition, cd)
  bad <- de_full_rank(~ condition + batch, cd)
  expect_true(ok$ok)
  expect_false(bad$ok)
  expect_lt(bad$rank, bad$ncoef)
  expect_match(bad$msg, "full rank")
})

test_that("de_coef_name matches the DESeq2 convention", {
  expect_equal(de_coef_name(c("condition", "treated", "control")),
               "condition_treated_vs_control")
})

test_that("de_design_factors / de_contrast_levels / de_relevel work on a dds", {
  skip_if_not_installed("DESeq2")
  dds <- make_mock_dds(n_genes = 40, n_per_group = 3, n_spike = 4, seed = 1)
  expect_true(all(c("condition", "bio_rep", "group") %in% de_design_factors(dds)))
  expect_equal(de_contrast_levels(dds, "condition"), c("control", "treated"))
  re <- de_relevel(dds, "condition", "treated")
  expect_equal(levels(SummarizedExperiment::colData(re)$condition)[1], "treated")
  expect_error(de_relevel(dds, "condition", "nope"), "not in")
})

# --- fitting + results -----------------------------------------------------

test_that("de_run + de_results carry both LFC variants and classify", {
  skip_if_not_installed("DESeq2")
  dds <- ensure_logcounts(make_mock_dds(n_genes = 120, n_per_group = 4, n_spike = 6, seed = 2))
  fit <- de_run(dds)
  df  <- de_results(fit, c("condition", "treated", "control"), shrink_type = "none")
  expect_true(all(c("baseMean", "log2FoldChange", "padj", "log2FoldChange_shrunk") %in% names(df)))
  expect_equal(nrow(df), nrow(dds))
  # shrink_type = "none" -> the shrunk column is all NA
  expect_true(all(is.na(df$log2FoldChange_shrunk)))
  cdf <- de_classify_table(df)
  expect_s3_class(cdf$DEG, "factor")
  expect_true(is.logical(cdf$sig))
})

test_that("de_shrink returns an aligned vector + method attr; apeglm/ashr populate it", {
  skip_if_not_installed("DESeq2")
  dds <- make_mock_dds(n_genes = 120, n_per_group = 4, n_spike = 6, seed = 3)
  fit <- de_run(dds)
  sh_none <- de_shrink(fit, c("condition", "treated", "control"), type = "none")
  expect_length(sh_none, nrow(fit))
  expect_true(all(is.na(sh_none)))
  expect_equal(attr(sh_none, "method"), "none")
  if (requireNamespace("apeglm", quietly = TRUE) || requireNamespace("ashr", quietly = TRUE)) {
    sh <- de_shrink(fit, c("condition", "treated", "control"), type = "apeglm")
    expect_length(sh, nrow(fit))
    expect_true(any(is.finite(sh)))
    expect_true(attr(sh, "method") %in% c("apeglm", "ashr", "normal"))
    df <- de_results(fit, c("condition", "treated", "control"), shrink_type = "apeglm")
    expect_false(is.null(attr(df, "shrink_method")))
  }
})

test_that("de_group_means guards empty / unknown groups", {
  skip_if_not_installed("DESeq2")
  dds <- ensure_logcounts(make_mock_dds(n_genes = 30, n_per_group = 3, n_spike = 2, seed = 5))
  ctrl <- colnames(dds)[SummarizedExperiment::colData(dds)$condition == "control"]
  trt  <- colnames(dds)[SummarizedExperiment::colData(dds)$condition == "treated"]
  gm <- de_group_means(dds, "logcounts", ctrl, trt)
  expect_equal(nrow(gm), nrow(dds))
  expect_error(de_group_means(dds, "logcounts", character(0), trt), "at least one")
  expect_error(de_group_means(dds, "logcounts", "NOPE", trt), "Unknown sample")
})

# --- size-factor edge case (poscounts fallback) ----------------------------

test_that("estimate_size_factors_endogenous falls back to poscounts on sparse data", {
  skip_if_not_installed("DESeq2")
  set.seed(7)
  m <- matrix(rpois(6 * 4, 5) + 1L, nrow = 6,
              dimnames = list(paste0("g", 1:6), paste0("S", 1:4)))
  for (i in seq_len(6)) m[i, (i %% 4) + 1] <- 0L        # every gene gets >=1 zero
  cd <- S4Vectors::DataFrame(condition = factor(rep(c("a", "b"), each = 2)),
                             row.names = colnames(m))
  rd <- S4Vectors::DataFrame(
    feature_class = factor(rep("endogenous", 6),
                           levels = c("endogenous", "spike_in", "exogenous")),
    row.names = rownames(m))
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = m, colData = cd, rowData = rd,
                                        design = ~ condition)
  # the standard estimator errors here; our helper must recover via poscounts
  expect_error(DESeq2::estimateSizeFactors(dds), "zero")
  suppressMessages(sf <- DESeq2::sizeFactors(estimate_size_factors_endogenous(dds)))
  expect_length(sf, 4)
  expect_true(all(is.finite(sf) & sf > 0))
})

# --- plot helpers ----------------------------------------------------------

test_that("de_clamp pulls out-of-range values to the limit and flags them", {
  cl <- de_clamp(c(-5, 0, 5, NA), lo = -2, hi = 2)
  expect_equal(cl$value, c(-2, 0, 2, NA))
  expect_equal(cl$clamped, c(TRUE, FALSE, TRUE, FALSE))
})

test_that("de_colour_resolve gives a discrete DEG scale and a continuous scale", {
  df <- data.frame(DEG = factor(c("up", "down"), levels = c("up", "down", "no_change")),
                   padj = c(0.01, 0.2))
  expect_null(de_colour_resolve(df, "__none__"))
  deg <- de_colour_resolve(df, "DEG")
  expect_equal(deg$kind, "discrete")
  expect_s3_class(deg$scale, "ScaleDiscrete")
  num <- de_colour_resolve(df, "padj")
  expect_equal(num$kind, "continuous")
})

test_that("de_ma_gg / de_volcano_gg / de_direct_gg build ggplots and validate inputs", {
  df <- data.frame(
    baseMean = c(10, 100, 1000),
    log2FoldChange = c(-3, 0.1, 4),
    padj = c(0.001, 0.5, 0.001),
    row.names = c("g1", "g2", "g3"))
  df <- de_classify_table(df)
  col <- de_colour_resolve(df, "DEG")
  expect_s3_class(de_ma_gg(df, colour = col, y_range = c(-2, 2)), "ggplot")
  expect_s3_class(de_volcano_gg(df, colour = col), "ggplot")
  means <- data.frame(id = c("g1", "g2"), control = c(1, 2), test = c(3, 4))
  expect_s3_class(de_direct_gg(means), "ggplot")
  expect_error(de_ma_gg(data.frame(x = 1)), "Missing column")
})
