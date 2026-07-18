# Phase 7: Expression page. Pure, exported helpers the mod_expression page calls --
# choosing/resolving the value matrix (incl. VST and normalized log-counts), the
# per-gene z-score used by the gene-set heatmap, the single-gene long frame, and
# the sample-count guards that decide which geoms (violin/box/dots) are offered.
# Unlike the PCA helpers these do NOT restrict to endogenous rows (a gene set or a
# searched single gene may be a spike-in/exogenous feature) -- except VST, which is
# inherently computed on the endogenous diagnostic rows via qc_vst().

#' Default expression value key for a dataset
#'
#' Picks the most interpretable value assay for the Expression page. Prefers
#' size-factor-normalized log-counts when size factors exist, then the
#' variance-stabilized transform, then the linear abundance assays (which the
#' plot's transform control can log), then log-counts / counts.
#'
#' @param dds A `DESeqDataSet`.
#' @return A single value key: `"norm_logcounts"`, `"vst"`, or a stored assay name.
#' @export
expr_default_assay <- function(dds) {
  an <- SummarizedExperiment::assayNames(dds)
  if (!is.null(DESeq2::sizeFactors(dds))) return("norm_logcounts")
  # No size factors: VST first (always attemptable, with a logcounts fallback at
  # resolve time), then the stored linear/log assays in preference order.
  for (a in c("vst", "TPM", "FPKM", "CPM", "logcounts")) {
    if (identical(a, "vst") || a %in% an) return(a)
  }
  if ("counts" %in% an) return("counts")
  if (length(an)) an[1] else "counts"
}

#' Resolve an Expression value matrix (genes x samples) + an honest label
#'
#' Handles the two computed inputs (`"vst"`, `"norm_logcounts"`) and any stored
#' assay. VST is computed on the endogenous diagnostic rows (via [qc_vst()]) and
#' falls back to log-counts on any failure; `"norm_logcounts"` estimates
#' endogenous size factors when absent. Stored assays are returned as-is (the
#' caller's transform/pseudocount then applies on top).
#'
#' @param dds A `DESeqDataSet`.
#' @param assay A value key (`"vst"`, `"norm_logcounts"`, or a stored assay name).
#' @return A list with `mat` (numeric genes x samples) and `label` (for the axis /
#'   subtitle, reflecting any fallback).
#' @export
expr_value_matrix <- function(dds, assay = "norm_logcounts") {
  fallback_logcounts <- function(why)
    list(mat = .qc_assay_matrix(dds, "logcounts"),
         label = paste0("logcounts (", why, ")"))

  if (identical(assay, "vst")) {
    res <- tryCatch(
      list(mat = as.matrix(SummarizedExperiment::assay(qc_vst(dds, blind = TRUE))),
           label = "VST (blind)"),
      error = function(e) NULL)
    return(res %||% fallback_logcounts("VST unavailable"))
  }

  if (identical(assay, "norm_logcounts")) {
    res <- tryCatch({
      d <- if (is.null(DESeq2::sizeFactors(dds))) estimate_size_factors_endogenous(dds) else dds
      nc <- as.matrix(DESeq2::counts(d, normalized = TRUE))
      list(mat = log2(nc + 1), label = "normalized log-counts (log2)")
    }, error = function(e) NULL)
    return(res %||% fallback_logcounts("normalized counts unavailable"))
  }

  list(mat = .qc_assay_matrix(dds, assay), label = assay)
}

#' Per-gene (row-wise) z-score of a value matrix
#'
#' Centers and scales each row to mean 0 / sd 1 -- the heatmap default so genes on
#' very different expression scales are comparable. Zero-variance rows (constant
#' expression) would divide by zero; they are returned as all-zero rows rather than
#' `NaN`. Non-finite inputs are treated as the row mean before scaling.
#'
#' @param mat A numeric genes x samples matrix.
#' @return A matrix the same shape, each row z-scored (constant rows -> 0).
#' @export
row_zscore <- function(mat) {
  mat <- as.matrix(mat)
  if (!nrow(mat) || !ncol(mat)) return(mat)
  mu <- rowMeans(mat, na.rm = TRUE)
  centered <- mat - mu
  sdv <- sqrt(rowSums(centered^2, na.rm = TRUE) / max(ncol(mat) - 1L, 1L))
  sdv[!is.finite(sdv) | sdv == 0] <- 1  # constant row -> centered values are 0
  z <- centered / sdv
  z[!is.finite(z)] <- 0
  dimnames(z) <- dimnames(mat)
  z
}

#' Which distribution geoms the single-gene plot should offer
#'
#' Given the per-group sample counts and the two thresholds, decides whether dots
#' are allowed / on by default and whether the violin/box distributions are shown.
#' Dots default on for small groups (individual points are informative) and are
#' disallowed once groups get large (overplotting, and the distribution geoms carry
#' the signal); the distribution geoms hide only when every group is too small to
#' summarize.
#'
#' @param group_sizes Integer vector of per-group sample counts (>= 1 each).
#' @param g1 Small-group threshold (default `getOption("ddsdashboard.expr_dots_max", 10)`).
#' @param g2 Hard cap above which dots are never drawn
#'   (default `getOption("ddsdashboard.expr_dots_hard", 50)`).
#' @return A list: `n_max` (largest group), `dots_allowed`, `dots_default`,
#'   `dist_shown` (whether violin/box are offered).
#' @export
expr_geom_availability <- function(group_sizes,
                                   g1 = getOption("ddsdashboard.expr_dots_max", 10L),
                                   g2 = getOption("ddsdashboard.expr_dots_hard", 50L)) {
  sizes <- as.integer(group_sizes[is.finite(group_sizes)])
  n_max <- if (length(sizes)) max(sizes) else 0L
  dots_allowed <- n_max < g2
  list(
    n_max = n_max,
    dots_allowed = dots_allowed,
    dots_default = dots_allowed && n_max <= g1,
    dist_shown = length(sizes) > 0 && any(sizes >= g1)
  )
}

#' Assemble the single-gene long data frame
#'
#' Joins one gene's per-sample values to a grouping vector (and optional colour
#' vector), dropping samples with a missing group. The grouping keeps its factor
#' level order when supplied as a factor.
#'
#' @param values Named (or same-order) numeric vector of the gene's values.
#' @param groups Per-sample grouping (character/factor), same length/order as
#'   `values`.
#' @param samples Sample ids (same length/order); defaults to `names(values)`.
#' @param colour Optional per-sample colour attribute (same length/order).
#' @return A data.frame with `sample`, `group` (factor), `value`, and (if given)
#'   `colour`.
#' @export
expr_long_frame <- function(values, groups, samples = names(values), colour = NULL) {
  n <- length(values)
  if (length(groups) != n) stop("groups must match values in length", call. = FALSE)
  if (is.null(samples)) samples <- as.character(seq_len(n))
  df <- data.frame(
    sample = as.character(samples),
    group  = if (is.factor(groups)) groups else factor(groups),
    value  = as.numeric(values),
    stringsAsFactors = FALSE)
  if (!is.null(colour)) df$colour <- colour
  df[!is.na(df$group), , drop = FALSE]
}
