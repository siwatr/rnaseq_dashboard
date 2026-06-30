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
#' correlation outlier check ([qc_within_group_correlation()]). Each threshold is
#' opt-in: a `NULL` threshold means that rule is **not applied** (no flag), so
#' blank inputs in the UI simply disable a check. [suggest_sample_thresholds()]
#' supplies sensible robust-fence values for the UI "Auto" button. Advisory
#' only - the UI never auto-drops samples.
#'
#' @param dds A `DESeqDataSet`.
#' @param group `colData` column defining replicate groups for the within-group
#'   check; defaults to a design variable.
#' @param lib_size_min,detected_min Lower thresholds; `NULL` disables the rule.
#' @param pct_mito_max Upper threshold for % mitochondrial; `NULL` disables it.
#' @param within_group_z Flag a sample when its within-group correlation z-score
#'   (within its group) is below `-within_group_z`; `NULL` disables the check.
#' @param method Correlation method for the within-group check.
#' @param spike Optional precomputed per-sample spike-in summary - the
#'   `per_sample` element of [spike_dose_response()] (columns `sample`,
#'   `n_spike_detected`, `n_points`, `slope`, `r_squared`, `lod`). When `NULL`,
#'   0-row, or sharing no sample with `dds`, every spike rule is silently
#'   disabled (back-compatible). Passed in (not recomputed) so the module reuses
#'   the cached dose-response artifact and this helper stays free of the ERCC
#'   concentration-resolution parameters.
#' @param pct_spike_min,pct_spike_max Two-sided fence on % spike-in (canonical
#'   `pct_spike` from [qc_per_sample_metrics()], not the `spike` frame's copy):
#'   under-spiked (`< pct_spike_min`) or over-spiked (`> pct_spike_max`). Each
#'   `NULL` disables that side.
#' @param min_spike_detected Flag when fewer than this many spike features are
#'   detected; `NULL` disables.
#' @param dose_r2_min Flag when the dose-response R^2 is below this; `NULL`
#'   disables.
#' @param dose_slope_min,dose_slope_max Flag when the dose-response slope falls
#'   outside `[dose_slope_min, dose_slope_max]` (ideal titration slope ~ 1); each
#'   `NULL` disables that side. (The `lod` column - lowest detected concentration -
#'   is carried for display context only and is intentionally not a filter rule:
#'   its units depend on the concentration source, so there is no portable cutoff.)
#' @return A `data.frame`, one row per sample, with the metric columns
#'   (`library_size`, `detected`, `pct_mito`, `pct_spike`, `within_group_corr`,
#'   `n_spike_detected`, `dose_r2`, `dose_slope`, `lod`), the per-reason logicals
#'   (`low_lib_size`, `low_detected`, `high_mito`, `within_group_outlier`,
#'   `low_spike`, `high_spike`, `low_spike_detected`, `low_dose_r2`, `bad_slope`,
#'   `few_spike_points`), `flagged`, and `reason`.
#' @export
flag_samples <- function(dds, group = NULL, lib_size_min = NULL,
                         detected_min = NULL, pct_mito_max = NULL,
                         within_group_z = 2, method = c("spearman", "pearson"),
                         spike = NULL, pct_spike_min = NULL, pct_spike_max = NULL,
                         min_spike_detected = NULL, dose_r2_min = NULL,
                         dose_slope_min = NULL, dose_slope_max = NULL) {
  method <- match.arg(method)
  m <- qc_per_sample_metrics(dds)
  wg <- qc_within_group_correlation(dds, method = method, group = group)
  wg_corr <- wg$mean_corr[match(m$sample, wg$sample)]
  wg_grp  <- wg$group[match(m$sample, wg$sample)]
  n_samp  <- nrow(m)
  no_flag <- rep(FALSE, n_samp)

  # Each rule is opt-in: a NULL threshold disables it (no flag).
  low_lib   <- if (is.null(lib_size_min)) no_flag else m$library_size < lib_size_min
  low_det   <- if (is.null(detected_min)) no_flag else m$detected < detected_min
  high_mito <- if (is.null(pct_mito_max)) no_flag else m$pct_mito > pct_mito_max

  # Within-group outlier: per-group z-score of mean_corr below -within_group_z.
  wg_out <- no_flag
  if (!is.null(within_group_z)) {
    for (g in unique(stats::na.omit(as.character(wg_grp)))) {
      idx <- which(as.character(wg_grp) == g & !is.na(wg_corr))
      if (length(idx) < 3L) next                    # z-score needs >= 3 points
      z <- (wg_corr[idx] - mean(wg_corr[idx])) / stats::sd(wg_corr[idx])
      wg_out[idx] <- is.finite(z) & z < -abs(within_group_z)
    }
  }

  # --- Spike-in (ERCC) criteria (P3e) --------------------------------------
  # All opt-in. The fit metrics come from the precomputed `spike` summary;
  # when it is absent or shares no sample, every spike rule is FALSE and the
  # spike columns are NA (so the frame never reports false data).
  spike_ok <- !is.null(spike) && is.data.frame(spike) && nrow(spike) > 0 &&
    any(m$sample %in% spike$sample)
  si <- if (spike_ok) match(m$sample, spike$sample) else rep(NA_integer_, n_samp)
  sp_detected <- if (spike_ok) as.numeric(spike$n_spike_detected[si]) else rep(NA_real_, n_samp)
  sp_npoints  <- if (spike_ok) as.numeric(spike$n_points[si])         else rep(NA_real_, n_samp)
  sp_r2       <- if (spike_ok) as.numeric(spike$r_squared[si])        else rep(NA_real_, n_samp)
  sp_slope    <- if (spike_ok) as.numeric(spike$slope[si])            else rep(NA_real_, n_samp)
  sp_lod      <- if (spike_ok) as.numeric(spike$lod[si])              else rep(NA_real_, n_samp)

  # `is.finite(x) & cond` is non-NA for NA x; the helper just gates on the
  # threshold being set. pct_spike is canonical from qc_per_sample_metrics.
  rule <- function(thr, cond) if (is.null(thr)) no_flag else cond
  low_spike  <- rule(pct_spike_min, is.finite(m$pct_spike) & m$pct_spike < pct_spike_min)
  high_spike <- rule(pct_spike_max, is.finite(m$pct_spike) & m$pct_spike > pct_spike_max)
  low_spike_detected <- rule(min_spike_detected,
                             is.finite(sp_detected) & sp_detected < min_spike_detected)
  low_dose_r2 <- rule(dose_r2_min, is.finite(sp_r2) & sp_r2 < dose_r2_min)
  bad_slope <- rule(dose_slope_min, is.finite(sp_slope) & sp_slope < dose_slope_min) |
               rule(dose_slope_max, is.finite(sp_slope) & sp_slope > dose_slope_max)
  # A dose rule (R^2/slope) is enabled but a sample has < 3 usable points, so its
  # slope/R^2 are NA and cannot fire - flag it explicitly so failed samples are
  # not silently passed (the detected-count rule independently catches zeros).
  dose_rule_on <- !is.null(dose_r2_min) || !is.null(dose_slope_min) || !is.null(dose_slope_max)
  few_points <- if (dose_rule_on) (is.finite(sp_npoints) & sp_npoints < 3L) else no_flag

  reason <- .reasons(list(
    "low library size"     = low_lib,
    "few detected features" = low_det,
    "high % mito"          = high_mito,
    "within-group outlier" = wg_out,
    "low % spike-in"       = low_spike,
    "high % spike-in"      = high_spike,
    "few detected spikes"  = low_spike_detected,
    "low dose-response R2" = low_dose_r2,
    "dose-response slope out of range" = bad_slope,
    "insufficient spike points" = few_points
  ))
  data.frame(
    sample = m$sample, library_size = m$library_size, detected = m$detected,
    pct_mito = m$pct_mito, pct_spike = m$pct_spike, within_group_corr = wg_corr,
    n_spike_detected = sp_detected, dose_r2 = sp_r2, dose_slope = sp_slope, lod = sp_lod,
    low_lib_size = low_lib, low_detected = low_det, high_mito = high_mito,
    within_group_outlier = wg_out, low_spike = low_spike, high_spike = high_spike,
    low_spike_detected = low_spike_detected, low_dose_r2 = low_dose_r2,
    bad_slope = bad_slope, few_spike_points = few_points,
    flagged = low_lib | low_det | high_mito | wg_out | low_spike | high_spike |
      low_spike_detected | low_dose_r2 | bad_slope | few_points,
    reason = reason,
    row.names = NULL, stringsAsFactors = FALSE, check.names = FALSE
  )
}

#' Suggested auto thresholds for sample flagging
#'
#' Robust fences (median +/- 3*MAD) for the sample-filter inputs, used to
#' populate the UI on load and behind its "Auto" button. A degenerate (zero-MAD)
#' fence yields `NA` so that input stays blank (rule disabled).
#'
#' @param dds A `DESeqDataSet`.
#' @param spike Optional per-sample spike summary ([spike_dose_response()]'s
#'   `per_sample`). When supplied, the list also carries spike suggestions:
#'   `pct_spike_min`/`pct_spike_max` (two-sided MAD fence on % spike-in),
#'   `dose_r2_min` (lower MAD fence on dose-response R^2), and a fixed advisory
#'   slope range `dose_slope_min = 0.5` / `dose_slope_max = 1.5` (ideal slope
#'   ~ 1). The lowest-detected-concentration has no auto default (its units
#'   depend on the concentration source) and is left to the user.
#' @return A list with `lib_size_min`, `detected_min`, `pct_mito_max` (plus the
#'   spike fields above when `spike` is supplied).
#' @export
suggest_sample_thresholds <- function(dds, spike = NULL) {
  m <- qc_per_sample_metrics(dds)
  lo <- function(x) { f <- .lower_fence(x); if (is.finite(f)) max(0, f) else NA_real_ }
  hi <- function(x) { f <- .upper_fence(x); if (is.finite(f)) f else NA_real_ }
  out <- list(lib_size_min = lo(m$library_size), detected_min = lo(m$detected),
              pct_mito_max = hi(m$pct_mito))
  if (!is.null(spike) && is.data.frame(spike) && nrow(spike) > 0) {
    out$pct_spike_min  <- lo(m$pct_spike)
    out$pct_spike_max  <- hi(m$pct_spike)
    out$dose_r2_min    <- lo(spike$r_squared)
    out$dose_slope_min <- 0.5
    out$dose_slope_max <- 1.5
  }
  out
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

# Semantic 3-colour scheme for the "Removal status" colour-by (reason-aware):
# green = QC pass, amber = suggested for some other reason, red = suggested for
# the reason of the current plot. A fixed scale (not the qualitative palette);
# the Palette page can override it via state$palette$other$removal_status. Shared
# by the QC plots and the PCA "Suggested removal" aesthetic so they agree.
.removal_palette <- c(pass = "#2CA02C", suggested_other = "#E6B800",
                      suggested_this = "#D62728")
.removal_labels <- c(pass = "QC pass", suggested_other = "Suggested drop (other)",
                     suggested_this = "Suggested drop (this reason)")
# Labels for the metric-free 2-level view (no "this reason" context): used by PCA
# + the RLE/density/spike colour-by, where only pass vs suggested is meaningful.
.removal_labels_2 <- c(pass = "QC pass", suggested_other = "Suggested removal",
                       suggested_this = "Suggested removal")

#' Resolve the removal-status colour vector
#'
#' The named colour vector for [removal_status()]'s levels: the project Palette
#' "Other" config when set, else the built-in green/amber/red. Pure so the QC
#' page and PCA resolve identical colours.
#'
#' @param config Optional `state$palette$other$removal_status` config
#'   (`list(name, colors, custom)`), or `NULL` for the built-in scheme.
#' @return A named character vector (`pass`/`suggested_other`/`suggested_this` -> hex).
#' @export
removal_status_colors <- function(config = NULL) {
  if (is.null(config)) return(.removal_palette)
  palette_discrete(names(.removal_palette), config$colors,
                   config$name %||% "Okabe-Ito", config$custom)
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

# Merge metadata for a restore: result is indexed by `ids`, on the CURRENT
# column schema (`current`). Rows present in `current` keep their (edited) values;
# rows only in `original` (being re-added) take original's values for the columns
# that exist there (columns added after load stay NA). Both args are DataFrames.
.merge_meta <- function(current, original, ids) {
  cols <- colnames(current)
  out  <- current[match(ids, rownames(current)), cols, drop = FALSE]
  rownames(out) <- ids
  miss <- which(!ids %in% rownames(current))
  if (length(miss)) {
    opos <- match(ids[miss], rownames(original))
    for (j in intersect(cols, colnames(original))) {
      out[[j]] <- .fill_rows(out[[j]], miss, original[[j]][opos])
    }
  }
  out
}

# Assign `vals` into rows `idx` of `target`, widening factor levels first so a
# value absent from the (possibly narrowed) target levels is not turned into NA.
.fill_rows <- function(target, idx, vals) {
  if (is.factor(target)) {
    target <- factor(target, levels = union(levels(target), as.character(vals)))
    target[idx] <- as.character(vals)
  } else {
    target[idx] <- vals
  }
  target
}

# Rebuild a dds from original's counts over a target feature/sample set, merging
# metadata so kept items keep their edits. Used by restore_samples/features.
# Raw counts are never edited, so kept-item counts in `original` equal `working`'s;
# pulling all counts from `original` is therefore authoritative. Note: rebuilding
# via DESeqDataSetFromMatrix flattens rowRanges to a plain rowData DataFrame, so
# GTF-derived GRanges/seqinfo (if any) are not carried through a restore.
.reconstruct <- function(working, original, feats, samples) {
  counts <- as.matrix(SummarizedExperiment::assay(original, "counts"))[feats, samples, drop = FALSE]
  rd <- .merge_meta(SummarizedExperiment::rowData(working), SummarizedExperiment::rowData(original), feats)
  cd <- .merge_meta(SummarizedExperiment::colData(working), SummarizedExperiment::colData(original), samples)
  # Build with a placeholder design first: re-added rows can be NA for a
  # user-added design column (or collapse a factor to one level), which
  # DESeqDataSetFromMatrix rejects at construction. Then restore the intended
  # design by slot assignment (not re-validated) so the DE page can full-rank it.
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts, colData = cd, rowData = rd, design = ~ 1)
  dds <- tryCatch({ DESeq2::design(dds) <- DESeq2::design(working); dds },
                  error = function(e) dds)
  # DESeqDataSetFromMatrix keeps only `counts`; restore the normalized assays
  # `working` carried, recomputed from the reconstructed counts.
  dds <- ensure_logcounts(dds)
  want <- intersect(c("CPM", "TPM", "FPKM"), SummarizedExperiment::assayNames(working))
  if (length(want)) dds <- add_normalized_assays(dds, which = want)
  if (!is.null(DESeq2::sizeFactors(working))) dds <- estimate_size_factors_endogenous(dds)
  dds
}

#' Restore removed samples (reset sample filtering)
#'
#' Re-adds every sample that is in `original` but not in `working`, returning the
#' full original sample set while **keeping edits on the samples that stayed**
#' (re-added samples get their original colData; columns added after load are
#' `NA`). Feature filtering and rowData edits are untouched. No-op when nothing
#' was removed.
#'
#' @param working The current `DESeqDataSet`.
#' @param original The originally loaded `DESeqDataSet`.
#' @return `working` with removed samples restored (assays refreshed, size
#'   factors re-estimated when set).
#' @export
restore_samples <- function(working, original) {
  if (!length(setdiff(colnames(original), colnames(working)))) return(working)
  .reconstruct(working, original, feats = rownames(working), samples = colnames(original))
}

#' Restore removed features (reset feature filtering)
#'
#' Re-adds every feature that is in `original` but not in `working`, returning the
#' full original feature set while **keeping edits on the features that stayed**
#' (re-added features get their original rowData; columns added after load are
#' `NA`). Sample filtering and colData edits are untouched. No-op when nothing
#' was removed.
#'
#' @param working The current `DESeqDataSet`.
#' @param original The originally loaded `DESeqDataSet`.
#' @return `working` with removed features restored (assays refreshed, size
#'   factors re-estimated when set).
#' @export
restore_features <- function(working, original) {
  if (!length(setdiff(rownames(original), rownames(working)))) return(working)
  .reconstruct(working, original, feats = rownames(original), samples = colnames(working))
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
