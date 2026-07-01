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

#' Estimate size factors on endogenous genes
#'
#' Uses DESeq2's median-of-ratios on the endogenous `controlGenes`. That estimator
#' **fails when every control gene has at least one zero** (all geometric means are
#' zero -- common on sparse/shallow data): we catch that and retry with the
#' positive-counts estimator (`type = "poscounts"`), which computes geometric means
#' over the non-zero entries. Protects the assay/QC/PCA paths and `DESeq()`.
#' @param dds A `DESeqDataSet`.
#' @return The `DESeqDataSet` with `sizeFactors` set (endogenous `controlGenes`).
#' @export
estimate_size_factors_endogenous <- function(dds) {
  rd <- SummarizedExperiment::rowData(dds)
  endo <- if ("feature_class" %in% colnames(rd)) which(rd$feature_class == "endogenous") else seq_len(nrow(dds))
  if (!length(endo)) endo <- seq_len(nrow(dds))
  tryCatch(
    DESeq2::estimateSizeFactors(dds, controlGenes = endo),
    error = function(e) {
      message("Size-factor estimation fell back to type='poscounts' (sparse data: ",
              conditionMessage(e), ")")
      DESeq2::estimateSizeFactors(dds, controlGenes = endo, type = "poscounts")
    })
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
