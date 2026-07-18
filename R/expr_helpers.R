# Phase 7: Expression page. Pure, exported helpers the mod_expression page calls --
# choosing/resolving the value matrix (incl. VST and normalized log-counts), the
# per-gene z-score used by the gene-set heatmap, the single-gene long frame, and
# the sample-count guards that decide which geoms (violin/box/dots) are offered.
# Unlike the PCA helpers these do NOT restrict to endogenous rows (a gene set or a
# searched single gene may be a spike-in/exogenous feature) -- except VST, which is
# inherently computed on the endogenous diagnostic rows via qc_vst().

#' Default expression value key for a dataset
#'
#' Picks the most interpretable value for the Expression page: size-factor
#' normalized log-counts when size factors exist, otherwise the variance-
#' stabilized transform. Both are computed value keys (resolved by
#' [expr_value_matrix()]); VST is always attemptable with a log-counts fallback
#' at resolve time, so no stored-assay fallthrough is needed here.
#'
#' @param dds A `DESeqDataSet`.
#' @return A single value key: `"norm_logcounts"` or `"vst"`.
#' @export
expr_default_assay <- function(dds) {
  if (!is.null(DESeq2::sizeFactors(dds))) "norm_logcounts" else "vst"
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
  n_row <- rowSums(!is.na(mat))                          # per-row non-NA count
  sdv <- sqrt(rowSums(centered^2, na.rm = TRUE) / pmax(n_row - 1L, 1L))
  sdv[!is.finite(sdv) | sdv == 0] <- 1  # constant row -> centered values are 0
  z <- centered / sdv
  z[!is.finite(z)] <- 0
  dimnames(z) <- dimnames(mat)
  z
}

#' Which distribution geoms the single-gene plot should offer
#'
#' Decides whether data points are allowed / on by default and whether the
#' violin/box distributions are shown, from the per-group sample counts. The two
#' concerns use separate thresholds: the distribution geoms appear once any group
#' is large enough to summarize (`dist_min`), while data points default on up to a
#' larger group size (`dots_max`) and are disallowed only for genuinely overplotted
#' groups (`dots_hard`).
#'
#' @param group_sizes Integer vector of per-group sample counts (>= 1 each).
#' @param dist_min Minimum group size for the violin/box geoms to be offered
#'   (default `getOption("ddsdashboard.expr_dist_min", 10)`).
#' @param dots_max Data points default ON when the largest group is below this
#'   (default `getOption("ddsdashboard.expr_dots_max", 100)`).
#' @param dots_hard Hard cap above which data points are never drawn
#'   (default `getOption("ddsdashboard.expr_dots_hard", 500)`).
#' @return A list: `n_max` (largest group), `dots_allowed`, `dots_default`,
#'   `dist_shown` (whether violin/box are offered).
#' @export
expr_geom_availability <- function(group_sizes,
                                   dist_min = getOption("ddsdashboard.expr_dist_min", 10L),
                                   dots_max = getOption("ddsdashboard.expr_dots_max", 100L),
                                   dots_hard = getOption("ddsdashboard.expr_dots_hard", 500L)) {
  sizes <- as.integer(group_sizes[is.finite(group_sizes)])
  n_max <- if (length(sizes)) max(sizes) else 0L
  dots_allowed <- n_max < dots_hard
  list(
    n_max = n_max,
    dots_allowed = dots_allowed,
    dots_default = dots_allowed && n_max < dots_max,
    dist_shown = length(sizes) > 0 && any(sizes >= dist_min)
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
  df <- df[!is.na(df$group), , drop = FALSE]
  df$group <- droplevels(df$group)                       # drop groups with no shown samples
  df
}
