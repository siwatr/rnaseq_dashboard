# Pure helpers behind the QC "Filtering" sub-tab (P3c): compute *advisory*
# suggestions of poor-quality samples/features (with per-reason detail) and the
# actual removal that subsets the dds. The app suggests; the user decides. All
# expression-based feature rules run on ENDOGENOUS features only - spike-in /
# exogenous are flagged, never dropped (rnaseq-bioc: they stay plottable and are
# exempt from filtering, mirroring controlGenes / variable-gene selection).
# Removal goes through state_mutate (mod_qc.R); these helpers stay Shiny-free.

# Robust fences for advisory sample flagging. median +/- k*MAD; a zero MAD
# (degenerate spread) yields an infinite fence so nothing is flagged. Advisory
# only - we never auto-drop samples (rnaseq-bioc: small-n MAD is unreliable).
.lower_fence <- function(x, k = 3) {
  m <- stats::median(x, na.rm = TRUE); s <- stats::mad(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) -Inf else m - k * s
}
.upper_fence <- function(x, k = 3) {
  m <- stats::median(x, na.rm = TRUE); s <- stats::mad(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) Inf else m + k * s
}

# Endogenous mask from the always-present feature_class column (all TRUE when
# the column is somehow absent, so the rules still apply).
.endogenous_mask <- function(dds) {
  rd <- SummarizedExperiment::rowData(dds)
  if (!"feature_class" %in% colnames(rd)) return(rep(TRUE, nrow(dds)))
  as.character(rd$feature_class) == "endogenous"
}

# Paste the fired reasons (named logical -> "label; label") per row.
.reasons <- function(flag_list) {
  labs <- names(flag_list)
  apply(do.call(cbind, flag_list), 1L, function(row) {
    hit <- labs[as.logical(row)]
    if (length(hit)) paste(hit, collapse = "; ") else ""
  })
}

#' Flag low-quality features for removal
#'
#' Computes, per feature, expression-based drop suggestions. Always flags
#' all-zero endogenous features; then either [edgeR::filterByExpr()] (the
#' design-aware default, when `edgeR` is installed) or a manual
#' total-count / detected-samples rule. Spike-in and exogenous features are
#' **exempt** (`suggested_drop = FALSE` always) - they are flagged on the
#' Feature tab elsewhere but never dropped by expression filtering.
#'
#' @param dds A `DESeqDataSet` with a `"counts"` assay.
#' @param use_filter_by_expr Use `edgeR::filterByExpr()` when available
#'   (default `TRUE`); falls back to the manual rule otherwise.
#' @param group Optional `colData` column giving the grouping for
#'   `filterByExpr`; defaults to the design / first discrete column.
#' @param min_count Manual rule only: minimum total count across samples.
#'   (Not forwarded to `filterByExpr`, whose `min.count` is per-sample and would
#'   mean something different; `filterByExpr` uses its own design-aware default.)
#' @param min_samples Manual rule only: minimum number of samples with count > 0
#'   (default `NULL` = not applied).
#' @return A `data.frame`, one row per feature, with `feature_id`,
#'   `feature_class`, `total_count`, `n_detected`, `mean_logcounts`,
#'   `drop_all_zero`, `fails_min_count`, `fails_min_samples`,
#'   `fails_filterByExpr`, `suggested_drop`, and `reason`.
#' @export
flag_features <- function(dds, use_filter_by_expr = TRUE, group = NULL,
                          min_count = 10, min_samples = NULL) {
  counts <- as.matrix(SummarizedExperiment::assay(dds, "counts"))
  endo <- .endogenous_mask(dds)
  total <- rowSums(counts)
  detected <- rowSums(counts > 0)

  mean_logcounts <- tryCatch({
    lc <- if ("logcounts" %in% SummarizedExperiment::assayNames(dds)) {
      as.matrix(SummarizedExperiment::assay(dds, "logcounts"))
    } else {
      logcounts_from_counts(counts)
    }
    rowMeans(lc)
  }, error = function(e) rep(NA_real_, nrow(counts)))

  # Manual-rule columns are gated to endogenous so the frame never reports a
  # spike-in/exogenous as "failing" a rule it is exempt from.
  drop_all_zero   <- endo & (total == 0)
  fails_min_count <- endo & (total < min_count)
  fails_min_samples <- if (is.null(min_samples)) rep(FALSE, nrow(counts)) else endo & (detected < min_samples)

  # filterByExpr on the endogenous submatrix only; map keep -> fail back. We do
  # NOT pass the manual min_count: filterByExpr's min.count is a per-sample count
  # (converted to a CPM cutoff), not a total, so reusing it would mis-filter.
  use_fbe <- isTRUE(use_filter_by_expr) && requireNamespace("edgeR", quietly = TRUE)
  fails_fbe <- rep(FALSE, nrow(counts))
  if (use_fbe && any(endo)) {
    grp <- .group_vector(dds, group)
    keep_endo <- tryCatch(
      edgeR::filterByExpr(counts[endo, , drop = FALSE], group = grp),
      error = function(e) rep(TRUE, sum(endo)))
    fails_fbe[endo] <- !keep_endo
  }

  active_fail <- if (use_fbe) fails_fbe else (fails_min_count | fails_min_samples)
  suggested_drop <- endo & (drop_all_zero | active_fail)

  reason <- .reasons(list(
    "all-zero"                 = endo & drop_all_zero,
    "below min total count"    = endo & !drop_all_zero & !use_fbe & fails_min_count,
    "detected in too few samples" = endo & !drop_all_zero & !use_fbe & fails_min_samples,
    "filterByExpr"             = endo & !drop_all_zero & use_fbe & fails_fbe
  ))

  rd <- SummarizedExperiment::rowData(dds)
  data.frame(
    feature_id     = rownames(dds),
    feature_class  = as.character(rd$feature_class %||% "endogenous"),
    total_count    = as.numeric(total),
    n_detected     = as.integer(detected),
    mean_logcounts = as.numeric(mean_logcounts),
    drop_all_zero  = drop_all_zero,
    fails_min_count = fails_min_count,
    fails_min_samples = fails_min_samples,
    fails_filterByExpr = fails_fbe,
    suggested_drop = suggested_drop,
    reason         = reason,
    row.names      = NULL,
    stringsAsFactors = FALSE,
    check.names    = FALSE
  )
}

# Per-sample group vector (named by sample) for a colData column, or NULL.
.group_vector <- function(dds, group = NULL) {
  cd <- as.data.frame(SummarizedExperiment::colData(dds))
  col <- group %||% .qc_default_group(dds)
  if (is.null(col) || !col %in% colnames(cd)) return(NULL)
  factor(cd[[col]])
}

#' Flag low-quality samples (advisory)
#'
#' Per-sample flags from [qc_per_sample_metrics()] plus the within-group
#' correlation outlier check ([qc_within_group_correlation()]). Thresholds
#' default to robust fences (median +/- 3*MAD) so the flags are populated before
#' any tuning; pass explicit values to override. Advisory only - the UI never
#' auto-drops samples.
#'
#' @param dds A `DESeqDataSet`.
#' @param group `colData` column defining replicate groups for the within-group
#'   check; defaults to a design variable.
#' @param lib_size_min,detected_min Lower thresholds; `NULL` uses a robust fence.
#' @param pct_mito_max Upper threshold for % mitochondrial; `NULL` uses a fence.
#' @param within_group_z Flag a sample when its within-group correlation z-score
#'   (within its group) is below `-within_group_z`.
#' @param method Correlation method for the within-group check.
#' @return A `data.frame`, one row per sample, with the metric columns plus
#'   `within_group_corr`, the per-reason logicals `low_lib_size`,
#'   `low_detected`, `high_mito`, `within_group_outlier`, `flagged`, and `reason`.
#' @export
flag_samples <- function(dds, group = NULL, lib_size_min = NULL,
                         detected_min = NULL, pct_mito_max = NULL,
                         within_group_z = 2, method = c("spearman", "pearson")) {
  method <- match.arg(method)
  m <- qc_per_sample_metrics(dds)
  wg <- qc_within_group_correlation(dds, method = method, group = group)
  wg_corr <- wg$mean_corr[match(m$sample, wg$sample)]
  wg_grp  <- wg$group[match(m$sample, wg$sample)]

  lib_thr <- lib_size_min %||% .lower_fence(m$library_size)
  det_thr <- detected_min %||% .lower_fence(m$detected)
  mito_thr <- pct_mito_max %||% .upper_fence(m$pct_mito)

  low_lib   <- m$library_size < lib_thr
  low_det   <- m$detected < det_thr
  high_mito <- m$pct_mito > mito_thr

  # Within-group outlier: per-group z-score of mean_corr below -within_group_z.
  wg_out <- rep(FALSE, nrow(m))
  for (g in unique(stats::na.omit(as.character(wg_grp)))) {
    idx <- which(as.character(wg_grp) == g & !is.na(wg_corr))
    if (length(idx) < 3L) next                      # z-score needs >= 3 points
    z <- (wg_corr[idx] - mean(wg_corr[idx])) / stats::sd(wg_corr[idx])
    wg_out[idx] <- is.finite(z) & z < -abs(within_group_z)
  }

  reason <- .reasons(list(
    "low library size"     = low_lib,
    "few detected features" = low_det,
    "high % mito"          = high_mito,
    "within-group outlier" = wg_out
  ))
  data.frame(
    sample = m$sample, library_size = m$library_size, detected = m$detected,
    pct_mito = m$pct_mito, pct_spike = m$pct_spike, within_group_corr = wg_corr,
    low_lib_size = low_lib, low_detected = low_det, high_mito = high_mito,
    within_group_outlier = wg_out,
    flagged = low_lib | low_det | high_mito | wg_out, reason = reason,
    row.names = NULL, stringsAsFactors = FALSE, check.names = FALSE
  )
}

#' Three-level removal status for plotting
#'
#' Maps a flagged vector (and, optionally, the hits for one specific reason) to a
#' factor for the reason-aware colour scheme: `pass` (not suggested),
#' `suggested_other` (suggested for a different reason), `suggested_this`
#' (suggested for the reason of the current plot).
#'
#' @param flagged Logical vector: is each item suggested for removal?
#' @param this_reason Optional logical vector (same length): does each item fire
#'   the reason of interest? `NULL` collapses to pass / suggested_other.
#' @return A factor with levels `pass`, `suggested_other`, `suggested_this`.
#' @export
removal_status <- function(flagged, this_reason = NULL) {
  flagged <- as.logical(flagged)
  flagged[is.na(flagged)] <- FALSE
  out <- ifelse(flagged, "suggested_other", "pass")
  if (!is.null(this_reason)) {
    tr <- as.logical(this_reason); tr[is.na(tr)] <- FALSE
    out[flagged & tr] <- "suggested_this"
  }
  factor(out, levels = c("pass", "suggested_other", "suggested_this"))
}

# Subset + keep downstream consistent: recompute library-size-dependent assays
# and (only if they were already set) re-estimate endogenous size factors.
.refit_after_subset <- function(dds) {
  dds <- refresh_assays(dds)
  if (!is.null(DESeq2::sizeFactors(dds))) dds <- estimate_size_factors_endogenous(dds)
  dds
}

#' Remove features from the dataset
#'
#' Subsets out `drop_ids` (ignored when absent), then refreshes normalized assays
#' and re-estimates size factors if they were set. A no-op when `drop_ids` is
#' empty. `feature_class` / `rowData` / `colData` are preserved by subsetting.
#'
#' @param dds A `DESeqDataSet`.
#' @param drop_ids Feature ids (rownames) to remove.
#' @return The subset `DESeqDataSet`.
#' @export
drop_features <- function(dds, drop_ids) {
  drop_ids <- intersect(as.character(drop_ids), rownames(dds))
  if (!length(drop_ids)) return(dds)
  .refit_after_subset(dds[setdiff(rownames(dds), drop_ids), , drop = FALSE])
}

#' Remove samples from the dataset
#'
#' Subsets out `drop_ids` (ignored when absent), then refreshes normalized assays
#' and re-estimates size factors if they were set. A no-op when `drop_ids` is empty.
#'
#' @param dds A `DESeqDataSet`.
#' @param drop_ids Sample ids (colnames) to remove.
#' @return The subset `DESeqDataSet`.
#' @export
drop_samples <- function(dds, drop_ids) {
  drop_ids <- intersect(as.character(drop_ids), colnames(dds))
  if (!length(drop_ids)) return(dds)
  dds <- dds[, setdiff(colnames(dds), drop_ids), drop = FALSE]
  # Drop now-empty factor levels so a later DESeq()/full-rank check on the DE
  # page does not trip over a level that no remaining sample uses.
  cd <- SummarizedExperiment::colData(dds)
  for (j in seq_along(cd)) if (is.factor(cd[[j]])) cd[[j]] <- droplevels(cd[[j]])
  SummarizedExperiment::colData(dds) <- cd
  .refit_after_subset(dds)
}

#' Before/after expression density data for the feature filter
#'
#' Long-format log-expression values over endogenous features, labelled `before`
#' (all endogenous) and `after` (only `keep_ids`), for the standard
#' filtering-density QC plot (limma/edgeR RNAseq123 workflow). Both curves are
#' sliced from the *same* value matrix, so they share one normalization basis
#' and the plot shows purely which low-expression mass the filter removes (the
#' "after" curve is not renormalized on the kept subset).
#'
#' @param dds A `DESeqDataSet`.
#' @param keep_ids Feature ids retained by the proposed filter.
#' @param assay Value assay; falls back to `log2(CPM + 1)` if absent.
#' @return A `data.frame` with `sample` (factor), `value`, and `status`
#'   (factor `before`/`after`).
#' @export
qc_filter_density <- function(dds, keep_ids, assay = "logcounts") {
  m <- .qc_diagnostic_matrix(dds, assay)                 # endogenous, full-library basis
  keep <- intersect(rownames(m), as.character(keep_ids))
  long <- function(mat, status) data.frame(
    sample = factor(rep(colnames(mat), each = nrow(mat)), levels = colnames(mat)),
    value  = as.numeric(mat), status = status, stringsAsFactors = FALSE)
  out <- rbind(long(m, "before"), long(m[keep, , drop = FALSE], "after"))
  out$status <- factor(out$status, levels = c("before", "after"))
  out
}
