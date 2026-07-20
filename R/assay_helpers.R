# Pure helpers behind the Assay tab: report/compute normalized assays and
# estimate size factors on endogenous genes. Library-size-dependent assays must
# be recomputed when the feature set changes (refresh_assays), which the
# filtering PR will call.

.feature_length <- function(dds) as.numeric(SummarizedExperiment::rowData(dds)$feature_length)

#' Is a complete, positive feature_length available?
#'
#' TPM/FPKM normalize across all features, so a partial length is unusable.
#' @param dds A `DESeqDataSet`.
#' @return Logical.
#' @export
has_feature_length <- function(dds) {
  rd <- SummarizedExperiment::rowData(dds)
  if (!"feature_length" %in% colnames(rd)) return(FALSE)
  len <- rd$feature_length
  length(len) == nrow(dds) && !anyNA(len) && all(len > 0)
}

#' Add normalized assays
#'
#' CPM is always computed; TPM/FPKM only when [has_feature_length()] is TRUE
#' (otherwise skipped with a warning).
#' @param dds A `DESeqDataSet`.
#' @param which Subset of `c("CPM","TPM","FPKM")` to add.
#' @return The `DESeqDataSet` with the requested assays attached.
#' @export
add_normalized_assays <- function(dds, which = c("CPM", "TPM", "FPKM")) {
  which <- intersect(which, c("CPM", "TPM", "FPKM"))
  counts <- as.matrix(SummarizedExperiment::assay(dds, "counts"))
  have_len <- has_feature_length(dds)
  len <- if (have_len) .feature_length(dds) else NULL
  if ("CPM" %in% which) SummarizedExperiment::assay(dds, "CPM") <- cpm(counts)
  if ("TPM" %in% which) {
    if (!have_len) warning("TPM skipped: no complete feature_length.", call. = FALSE)
    else SummarizedExperiment::assay(dds, "TPM") <- tpm(counts, len)
  }
  if ("FPKM" %in% which) {
    if (!have_len) warning("FPKM skipped: no complete feature_length.", call. = FALSE)
    else SummarizedExperiment::assay(dds, "FPKM") <- fpkm(counts, len)
  }
  dds
}

# --- Size-factor normalization config (carried on the dds via metadata) -------
# Size factors are decoupled from assay assignment (they are a DESeq2 median-of-
# ratios normalization, unrelated to the library-size/length CPM/TPM/FPKM assays).
# The config records HOW they were estimated so the choice survives edits (undo /
# subset carry metadata) and structural re-estimations reuse it, and the Size-
# factors tab can display/restore it. `control`: "endogenous" (feature_class ==
# endogenous, the default) / "spike_in" / "custom" (an explicit control-gene id
# set). `type`: the DESeq2::estimateSizeFactors estimator. `provenance`: "auto" (a
# default we materialized), "loaded" (the object arrived with size factors), or
# "user" (chosen on the tab).
.SF_CONFIG_KEY <- "sizefactor_config"

#' Default size-factor normalization config
#' @return `list(control, custom_ids, type, provenance)`.
#' @export
default_sizefactor_config <- function() {
  list(control = "endogenous", custom_ids = character(0),
       type = "ratio", provenance = "auto")
}

#' Read / write the size-factor config carried on a dds
#'
#' Stored under `metadata(dds)$sizefactor_config`; [sizefactor_config()] fills any
#' missing field from [default_sizefactor_config()]. (`metadata()` is an S4Vectors
#' accessor -- SummarizedExperiment re-exports the setter but not the getter, so we
#' import both from S4Vectors and call them bare.)
#' @param dds A `DESeqDataSet`.
#' @return `sizefactor_config()`: the config list. `set_sizefactor_config()`: `dds`.
#' @importFrom S4Vectors metadata metadata<-
#' @export
sizefactor_config <- function(dds) {
  cfg <- tryCatch(metadata(dds)[[.SF_CONFIG_KEY]], error = function(e) NULL)
  if (is.null(cfg)) return(default_sizefactor_config())
  utils::modifyList(default_sizefactor_config(), cfg)
}

#' @rdname sizefactor_config
#' @param config A config list (see [default_sizefactor_config()]).
#' @export
set_sizefactor_config <- function(dds, config) {
  md <- tryCatch(metadata(dds), error = function(e) list())
  if (!is.list(md)) md <- list()
  md[[.SF_CONFIG_KEY]] <- config
  metadata(dds) <- md
  dds
}

# Control-gene row indices for a config. Returns integer(0) when the requested set
# resolves to nothing (e.g. "spike_in" with no spike-ins, "custom" whose ids match
# no feature) -- estimate_size_factors() decides what to do rather than silently
# widening to all genes and mislabelling the result.
.sf_control_index <- function(dds, config) {
  rd <- SummarizedExperiment::rowData(dds)
  fc <- if ("feature_class" %in% colnames(rd)) as.character(rd$feature_class) else NULL
  switch(config$control %||% "endogenous",
    endogenous = if (!is.null(fc)) which(fc == "endogenous") else seq_len(nrow(dds)),
    spike_in   = if (!is.null(fc)) which(fc == "spike_in") else integer(0),
    custom     = which(rownames(dds) %in% (config$custom_ids %||% character(0))),
    all_genes  = seq_len(nrow(dds)),
    seq_len(nrow(dds)))
}

# Number of coefficients the design implies (for the iterate identifiability guard).
.sf_design_rank <- function(dds) {
  tryCatch(
    ncol(stats::model.matrix(DESeq2::design(dds),
                             data = as.data.frame(SummarizedExperiment::colData(dds)))),
    error = function(e) 1L)
}

#' Estimate size factors under a control-gene / type config
#'
#' Generalizes [estimate_size_factors_endogenous()]: the control-gene set is
#' chosen by `config$control` (`"endogenous"`/`"spike_in"`/`"custom"`) and the
#' estimator by `config$type` (`"ratio"`/`"poscounts"`/`"iterate"`). Size factors
#' are per-sample scalars, so they are estimated on the control-gene **row-subset**
#' and the full dds **inherits** them -- this honors the control set for *every*
#' estimator, including `"iterate"` (which ignores DESeq2's `controlGenes`), and is
#' mathematically identical to `controlGenes=` for ratio/poscounts (a gene's
#' geometric mean depends only on its own row). An empty control set is an
#' **error** for `"spike_in"`/`"custom"`; only the degenerate `"endogenous"` case
#' with no endogenous rows falls back to all genes. `"iterate"` (design-based) is
#' **refused** when the control set is too small to identify the model. A ratio run
#' whose control subset is all-zero-containing retries once with `"poscounts"`;
#' non-finite factors (degenerate control set) warn. Records the config on the
#' returned dds so the choice travels with the object.
#' @param dds A `DESeqDataSet`.
#' @param config A config list (default [default_sizefactor_config()]).
#' @return `dds` with `sizeFactors` set and the config recorded.
#' @export
estimate_size_factors <- function(dds, config = default_sizefactor_config()) {
  config <- utils::modifyList(default_sizefactor_config(), config)
  idx <- .sf_control_index(dds, config)
  if (!length(idx)) {
    if (identical(config$control, "endogenous")) {
      idx <- seq_len(nrow(dds))                          # degenerate: no endogenous rows
      message("No endogenous features; size factors estimated on all genes.")
    } else {
      stop(sprintf("No control genes resolved for control = '%s' (%s). Adjust the selection.",
                   config$control,
                   if (identical(config$control, "custom")) "custom ids matched no feature"
                   else "no spike-in features in this dataset"), call. = FALSE)
    }
  }
  type <- config$type %||% "ratio"
  if (identical(type, "iterate")) {                      # design-based -> needs rank
    p <- .sf_design_rank(dds)
    if (length(idx) < p)
      stop(sprintf("The 'iterate' estimator needs at least %d control genes for this design (got %d). Use more control genes, a simpler design, or the ratio/poscounts estimator.",
                   p, length(idx)), call. = FALSE)
  }

  # Estimate on the control-gene ROW-subset (all columns kept, in order), then have
  # the full dds inherit the per-sample factors -- see the note above.
  ctrl <- dds[idx, , drop = FALSE]
  fit <- tryCatch(
    DESeq2::estimateSizeFactors(ctrl, type = type),
    error = function(e) {
      if (identical(type, "ratio") &&
          grepl("every gene contains at least one zero", conditionMessage(e))) {
        message("Size-factor estimation fell back to type='poscounts' (sparse control set).")
        DESeq2::estimateSizeFactors(ctrl, type = "poscounts")
      } else stop(e)
    })
  sf <- tryCatch(DESeq2::sizeFactors(fit), error = function(e) NULL)
  if (is.null(sf))
    stop("Size factors unavailable - the dataset uses normalizationFactors (e.g. a tximport avgTxLength object), which control-gene size factors do not support.",
         call. = FALSE)
  if (any(!is.finite(sf)))
    warning("Some size factors are non-finite - the control set may be degenerate ",
            "(all-zero control genes).", call. = FALSE)
  DESeq2::sizeFactors(dds) <- sf[colnames(dds)]          # name-aligned per-sample inherit
  set_sizefactor_config(dds, config)
}

#' @rdname estimate_size_factors
#' @return `estimate_size_factors_endogenous()`: `dds` with endogenous-control
#'   size factors (a thin wrapper kept for the many internal callers).
#' @export
estimate_size_factors_endogenous <- function(dds) {
  estimate_size_factors(dds, config = list(control = "endogenous"))
}

# Re-estimate under the dds's OWN stored config, after a structural edit. No-op
# when size factors were never set, AND when they were "loaded" (externally
# provided): we cannot reproduce an upstream method, so we keep the (auto-
# subsetted) provided vector rather than silently swapping in a fresh estimate.
reestimate_size_factors <- function(dds) {
  if (is.null(tryCatch(DESeq2::sizeFactors(dds), error = function(e) NULL))) return(dds)
  cfg <- sizefactor_config(dds)
  if (identical(cfg$provenance, "loaded")) return(dds)
  estimate_size_factors(dds, config = cfg)
}

#' Ensure a dds carries size factors (materialize a default when absent)
#'
#' Respects size factors the object already carries (records provenance
#' `"loaded"`); otherwise estimates them under `config` (provenance `"auto"`).
#' Called at load so downstream consumers (DE/PCA/Expression) read a consistent,
#' visible normalization instead of each re-deriving it transiently.
#' @param dds A `DESeqDataSet`.
#' @param config Config to use when estimating (default endogenous).
#' @return `dds` with size factors + config recorded.
#' @export
ensure_size_factors <- function(dds, config = default_sizefactor_config()) {
  sf <- tryCatch(DESeq2::sizeFactors(dds), error = function(e) NULL)
  if (!is.null(sf)) {
    cur <- sizefactor_config(dds)
    if (!identical(cur$provenance, "user")) cur$provenance <- "loaded"
    return(set_sizefactor_config(dds, cur))
  }
  cfg <- utils::modifyList(default_sizefactor_config(), config)
  cfg$provenance <- "auto"
  estimate_size_factors(dds, cfg)
}

#' Recompute every present normalized/log assay from current counts
#'
#' Call after the feature set changes so library-size-dependent assays
#' (logcounts/CPM/TPM/FPKM) stay correct.
#' @param dds A `DESeqDataSet`.
#' @return The `DESeqDataSet` with its normalized assays refreshed.
#' @export
refresh_assays <- function(dds) {
  present <- intersect(c("logcounts", "CPM", "TPM", "FPKM"),
                       SummarizedExperiment::assayNames(dds))
  if (!length(present)) return(dds)
  counts <- as.matrix(SummarizedExperiment::assay(dds, "counts"))
  len <- if (has_feature_length(dds)) .feature_length(dds) else NULL
  for (a in present) {
    val <- switch(a,
      logcounts = logcounts_from_counts(counts),
      CPM = cpm(counts),
      TPM = if (!is.null(len)) tpm(counts, len) else NULL,
      FPKM = if (!is.null(len)) fpkm(counts, len) else NULL)
    if (!is.null(val)) SummarizedExperiment::assay(dds, a) <- val
  }
  dds
}
