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

# Control-gene row indices for a config (empty resolution -> all-genes fallback).
.sf_control_index <- function(dds, config) {
  rd <- SummarizedExperiment::rowData(dds)
  fc <- if ("feature_class" %in% colnames(rd)) as.character(rd$feature_class) else NULL
  idx <- switch(config$control %||% "endogenous",
    endogenous = if (!is.null(fc)) which(fc == "endogenous") else seq_len(nrow(dds)),
    spike_in   = if (!is.null(fc)) which(fc == "spike_in") else integer(0),
    custom     = which(rownames(dds) %in% (config$custom_ids %||% character(0))),
    seq_len(nrow(dds)))
  if (!length(idx)) idx <- seq_len(nrow(dds))            # never an empty control set
  idx
}

#' Estimate size factors under a control-gene / type config
#'
#' Generalizes [estimate_size_factors_endogenous()]: the control-gene set is
#' chosen by `config$control` (`"endogenous"`/`"spike_in"`/`"custom"`) and the
#' estimator by `config$type`. Keeps the poscounts fallback -- a `type = "ratio"`
#' run that **fails when every control gene has at least one zero** (all geometric
#' means zero, common on sparse data) retries with `type = "poscounts"`. Records
#' the config on the returned dds so the choice travels with the object.
#' @param dds A `DESeqDataSet`.
#' @param config A config list (default [default_sizefactor_config()]).
#' @return `dds` with `sizeFactors` set and the config recorded.
#' @export
estimate_size_factors <- function(dds, config = default_sizefactor_config()) {
  config <- utils::modifyList(default_sizefactor_config(), config)
  idx <- .sf_control_index(dds, config)
  type <- config$type %||% "ratio"
  out <- tryCatch(
    DESeq2::estimateSizeFactors(dds, controlGenes = idx, type = type),
    error = function(e) {
      if (!identical(type, "ratio")) stop(e)
      message("Size-factor estimation fell back to type='poscounts' (sparse data: ",
              conditionMessage(e), ")")
      DESeq2::estimateSizeFactors(dds, controlGenes = idx, type = "poscounts")
    })
  set_sizefactor_config(out, config)
}

#' @rdname estimate_size_factors
#' @return `estimate_size_factors_endogenous()`: `dds` with endogenous-control
#'   size factors (a thin wrapper kept for the many internal callers).
#' @export
estimate_size_factors_endogenous <- function(dds) {
  estimate_size_factors(dds, config = list(control = "endogenous"))
}

# Re-estimate under the dds's OWN stored config, after a structural edit. A no-op
# when size factors were never set (respects the "not estimated yet" state).
reestimate_size_factors <- function(dds) {
  if (is.null(tryCatch(DESeq2::sizeFactors(dds), error = function(e) NULL))) return(dds)
  estimate_size_factors(dds, config = sizefactor_config(dds))
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
