# Expression helpers (R/expr_helpers.R) for Phase 7.

test_that("expr_default_assay prefers norm_logcounts when size factors exist", {
  skip_if_not_installed("DESeq2")
  dds <- ensure_logcounts(make_mock_dds(n_genes = 80, n_per_group = 3, n_spike = 6, seed = 1))
  expect_equal(expr_default_assay(dds), "vst")            # no size factors yet -> VST first
  dds2 <- estimate_size_factors_endogenous(dds)
  expect_equal(expr_default_assay(dds2), "norm_logcounts")
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

test_that("expr_geom_availability splits the distribution + dots thresholds", {
  # small groups: dots on by default; distributions hidden (all groups < dist_min)
  a <- expr_geom_availability(c(4, 5), dist_min = 10, dots_max = 100, dots_hard = 500)
  expect_true(a$dots_allowed); expect_true(a$dots_default); expect_false(a$dist_shown)

  # a group at/above dist_min: distributions shown
  a2 <- expr_geom_availability(c(4, 12), dist_min = 10)
  expect_true(a2$dist_shown)

  # groups < dots_max (100): dots still on by default (and distributions shown)
  b <- expr_geom_availability(c(12, 90), dist_min = 10, dots_max = 100)
  expect_true(b$dots_allowed); expect_true(b$dots_default); expect_true(b$dist_shown)

  # >= dots_max but < dots_hard: allowed, off by default
  c1 <- expr_geom_availability(c(120), dots_max = 100, dots_hard = 500)
  expect_true(c1$dots_allowed); expect_false(c1$dots_default)

  # >= dots_hard: disallowed (overplotted)
  d <- expr_geom_availability(c(600), dots_hard = 500)
  expect_false(d$dots_allowed); expect_false(d$dots_default)

  # all groups tiny: distributions hidden, dots default on
  e <- expr_geom_availability(c(2, 3))
  expect_false(e$dist_shown); expect_true(e$dots_default)
})

test_that("expr_set_aggregate reduces a set to one value per sample + accounting", {
  set.seed(11)
  mat <- matrix(rpois(6 * 5, 20), nrow = 6,
                dimnames = list(paste0("g", 1:6), paste0("S", 1:5)))
  counts <- mat
  counts["g6", ] <- 0                                   # g6 not expressed anywhere
  ids <- c("g1", "g2", "g3", "g6", "absent1", "absent2")

  a <- expr_set_aggregate(mat, ids, counts = counts, method = "mean",
                          zscore = FALSE, only_expressed = TRUE)
  expect_equal(a$n_total, 6L)                           # includes absent + g6
  expect_equal(a$n_present, 4L)                          # g1,g2,g3,g6
  expect_equal(a$n_absent, 2L)
  expect_equal(a$n_used, 3L)                             # g6 dropped by only_expressed
  expect_false("g6" %in% a$ids_used)
  expect_length(a$values, 5L)
  expect_equal(names(a$values), colnames(mat))
  # mean of the 3 used rows, per column
  expect_equal(unname(a$values),
               unname(colMeans(mat[c("g1", "g2", "g3"), ])), tolerance = 1e-8)

  # only_expressed = FALSE keeps g6 (all-zero counts) -> 4 genes used
  b <- expr_set_aggregate(mat, ids, counts = counts, zscore = FALSE,
                          only_expressed = FALSE)
  expect_equal(b$n_used, 4L)
})

test_that("expr_set_aggregate z-scores per gene and counts non-varying genes", {
  mat <- rbind(varies = c(1, 2, 3, 4), flat = c(5, 5, 5, 5))
  colnames(mat) <- paste0("S", 1:4)
  counts <- matrix(1, nrow = 2, ncol = 4, dimnames = dimnames(mat))
  a <- expr_set_aggregate(mat, rownames(mat), counts = counts, method = "mean",
                          zscore = TRUE, only_expressed = TRUE)
  expect_equal(a$n_used, 2L)
  expect_equal(a$n_nonvar, 1L)                           # the flat row
  # z-scored mean: flat -> all 0, varies -> its own z-score; averaged
  z <- row_zscore(mat)
  expect_equal(unname(a$values), unname(colMeans(z)), tolerance = 1e-8)
  expect_equal(unname(mean(a$values)), 0, tolerance = 1e-8)
})

test_that("expr_set_aggregate returns NULL values when nothing survives", {
  mat <- matrix(1:4, nrow = 2, dimnames = list(c("g1", "g2"), c("S1", "S2")))
  a <- expr_set_aggregate(mat, c("nope1", "nope2"), zscore = FALSE)
  expect_null(a$values)
  expect_equal(a$n_used, 0L)
  expect_equal(a$n_total, 2L)
  expect_equal(a$n_present, 0L)
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

# --- P7c: gene-set heatmap helpers -----------------------------------------

test_that("expr_heatmap_matrix selects present genes, drops all-zero, keeps shape", {
  mat <- matrix(c(1, 2, 3, 4, 5, 6, 0, 0), nrow = 4, byrow = TRUE,
                dimnames = list(c("g1", "g2", "g3", "g0"), c("S1", "S2")))
  counts <- matrix(c(5, 6, 7, 8, 9, 10, 0, 0), nrow = 4, byrow = TRUE,
                   dimnames = list(c("g1", "g2", "g3", "g0"), c("S1", "S2")))
  h <- expr_heatmap_matrix(mat, c("g1", "g2", "g0", "absent"), counts = counts,
                           zscore = FALSE, only_expressed = TRUE)
  expect_equal(h$n_total, 4L)                 # 4 authored (deduped)
  expect_equal(h$n_present, 3L)               # g1,g2,g0 present; "absent" not
  expect_equal(h$n_absent, 1L)
  expect_equal(h$n_used, 2L)                  # g0 dropped (all-zero counts)
  expect_equal(rownames(h$mat), c("g1", "g2"))
  expect_equal(ncol(h$mat), 2L)
})

test_that("expr_heatmap_matrix z-scores rows and counts non-varying ones", {
  mat <- matrix(c(1, 3, 5,      # varies
                  2, 2, 2),     # flat -> z-score 0, counted as non-varying
                nrow = 2, byrow = TRUE, dimnames = list(c("g1", "g2"), c("A", "B", "C")))
  h <- expr_heatmap_matrix(mat, c("g1", "g2"), counts = NULL,
                           zscore = TRUE, only_expressed = FALSE)
  expect_equal(h$n_used, 2L)
  expect_equal(h$n_nonvar, 1L)                # the flat row
  expect_true(all(h$mat["g2", ] == 0))        # constant -> 0, not NaN
  expect_equal(unname(rowMeans(h$mat)), c(0, 0), tolerance = 1e-8)
})

test_that("expr_heatmap_matrix returns NULL mat when no gene survives", {
  mat <- matrix(1:4, 2, dimnames = list(c("g1", "g2"), c("S1", "S2")))
  h <- expr_heatmap_matrix(mat, c("x", "y"), zscore = TRUE)
  expect_null(h$mat)
  expect_equal(h$n_used, 0L)
})

test_that("expr_heatmap_labels maps to a column with duplicate + NA-safe fallback", {
  meta <- data.frame(gene_name = c("Actb", "Actb", NA, ""),
                     row.names = c("id1", "id2", "id3", "id4"),
                     stringsAsFactors = FALSE)
  lab <- expr_heatmap_labels(c("id1", "id2", "id3", "id4"), "gene_name",
                             meta = meta, meta_keys = rownames(meta))
  expect_equal(lab, c("Actb", "Actb", "id3", "id4"))   # duplicates kept; NA/"" -> id
  # "__id__" and an unknown column both return the keys verbatim
  expect_equal(expr_heatmap_labels(c("id1", "id2"), "__id__", meta = meta), c("id1", "id2"))
  expect_equal(expr_heatmap_labels(c("id1"), "nope", meta = meta), "id1")
})

test_that("heatmap_label_default follows the size threshold", {
  expect_equal(heatmap_label_default(50, 50), "all")
  expect_equal(heatmap_label_default(51, 50), "none")
})

test_that("expr_label_coverage counts searched labels that cannot be shown", {
  cov <- expr_label_coverage(c("a", "b", "c", "a"), present_keys = c("a", "b", "z"))
  expect_equal(cov$n_selected, 3L)            # deduped
  expect_equal(cov$n_shown, 2L)               # a, b
  expect_equal(cov$n_hidden, 1L)              # c
})

test_that("expr_symmetric_limits centres on zero (and guards empty/zero)", {
  expect_equal(expr_symmetric_limits(c(-1, 0.5, 2)), c(-2, 2))
  expect_equal(expr_symmetric_limits(c(0, 0, 0)), c(-1, 1))
  expect_equal(expr_symmetric_limits(c(NA, NaN, Inf)), c(-1, 1))
})

# --- P7d: k-means split helpers --------------------------------------------

test_that("expr_kmeans is reproducible, RNG-safe, and separates clear groups", {
  m <- rbind(matrix(rnorm(20, 0), 5, 4), matrix(rnorm(20, 8), 5, 4))
  rownames(m) <- paste0("g", 1:10)
  set.seed(42); before <- runif(1)
  cl <- expr_kmeans(m, 2, seed = 1)
  expect_named(cl, paste0("g", 1:10))
  expect_equal(sort(as.integer(table(cl))), c(5L, 5L))   # two clean groups of 5
  expect_equal(cl, expr_kmeans(m, 2, seed = 1))          # reproducible
  set.seed(42); expect_equal(runif(1), before)           # global RNG untouched
})

test_that("expr_kmeans with a blank/NA seed clusters and advances the global RNG", {
  m <- rbind(matrix(rnorm(20, 0), 5, 4), matrix(rnorm(20, 8), 5, 4))
  rownames(m) <- paste0("g", 1:10)
  set.seed(7)
  before <- get(".Random.seed", envir = .GlobalEnv)
  cl <- expr_kmeans(m, 2, seed = NA)                    # no seed
  expect_named(cl, paste0("g", 1:10))
  expect_equal(sort(as.integer(table(cl))), c(5L, 5L))  # still a valid partition
  after <- get(".Random.seed", envir = .GlobalEnv)
  expect_false(identical(before, after))                # RNG advanced (not restored)
})

test_that("expr_kmeans handles the degenerate k / n edge cases", {
  m <- matrix(rnorm(8), 4, 2, dimnames = list(paste0("g", 1:4), c("A", "B")))
  expect_null(expr_kmeans(m, 1))                          # k < 2 -> no split
  expect_null(expr_kmeans(m[1, , drop = FALSE], 2))       # < 2 rows -> no split
  expect_equal(length(unique(expr_kmeans(m, 10))), 4L)    # k >= n -> each its own
})

test_that("expr_kmeans relabels clusters by decreasing size (1 = largest)", {
  m <- rbind(matrix(0, 2, 3), matrix(9, 5, 3))            # groups of 2 and 5
  rownames(m) <- paste0("g", 1:7)
  cl <- expr_kmeans(m, 2, seed = 1)
  expect_equal(sum(cl == 1), 5L)                          # cluster 1 = the larger
  expect_equal(sum(cl == 2), 2L)
})

test_that("split_with_counts labels slices 'C<id>\\n(count)', ordered", {
  cl <- c(g1 = 2, g2 = 1, g3 = 1, g4 = 2, g5 = 1)
  f <- split_with_counts(cl)
  expect_s3_class(f, "factor")
  expect_equal(levels(f), c("C1\n(3)", "C2\n(2)"))        # numeric order, two-line labels
  expect_equal(as.character(f[1]), "C2\n(2)")
  expect_equal(levels(split_with_counts(cl, prefix = "K")),
               c("K1\n(3)", "K2\n(2)"))
})

test_that("cluster_membership groups ids by cluster in order", {
  cl <- c(g1 = 2, g2 = 1, g3 = 1, g4 = 2, g5 = 1)
  mem <- cluster_membership(cl)
  expect_equal(names(mem), c("1", "2"))
  expect_equal(mem[["1"]], c("g2", "g3", "g5"))
  expect_equal(mem[["2"]], c("g1", "g4"))
})
