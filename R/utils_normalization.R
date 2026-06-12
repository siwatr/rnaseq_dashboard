# Normalization helpers. Pure functions on a counts matrix (features x
# samples) so they are reusable outside the app and unit-testable. Formulas
# follow the rnaseq-bioc skill: length-normalize before per-sample scaling.

.as_count_matrix <- function(counts) {
  if (is.data.frame(counts)) counts <- as.matrix(counts)
  if (!is.matrix(counts) || !is.numeric(counts)) {
    stop("`counts` must be a numeric matrix (features x samples).", call. = FALSE)
  }
  if (anyNA(counts) || any(counts < 0)) {
    stop("`counts` must be non-negative and contain no NA.", call. = FALSE)
  }
  counts
}

.check_feature_length <- function(feature_length, counts) {
  if (length(feature_length) != nrow(counts)) {
    stop("`feature_length` must have one value per feature (nrow(counts)).",
         call. = FALSE)
  }
  if (anyNA(feature_length) || any(feature_length <= 0)) {
    stop("`feature_length` must be positive and contain no NA.", call. = FALSE)
  }
  feature_length
}

.check_lib_size <- function(lib) {
  if (any(lib == 0)) {
    bad <- if (is.null(names(lib))) which(lib == 0) else names(lib)[lib == 0]
    stop("Sample(s) with zero total counts: ", paste(bad, collapse = ", "),
         call. = FALSE)
  }
  lib
}

#' Counts per million (CPM)
#'
#' @param counts Numeric matrix of raw counts (features x samples).
#' @return Matrix of CPM values, same dimensions as `counts`.
#' @export
cpm <- function(counts) {
  counts <- .as_count_matrix(counts)
  lib <- .check_lib_size(colSums(counts))
  sweep(counts, 2L, lib, "/") * 1e6
}

#' Transcripts per million (TPM)
#'
#' Length-normalizes first, then scales each sample to sum to 1e6.
#'
#' @param counts Numeric matrix of raw counts (features x samples).
#' @param feature_length Effective length per feature (bases), length nrow(counts).
#' @return Matrix of TPM values.
#' @export
tpm <- function(counts, feature_length) {
  counts <- .as_count_matrix(counts)
  feature_length <- .check_feature_length(feature_length, counts)
  rate <- counts / feature_length
  denom <- .check_lib_size(colSums(rate))
  sweep(rate, 2L, denom, "/") * 1e6
}

#' Fragments per kilobase per million (FPKM)
#'
#' @param counts Numeric matrix of raw counts (features x samples).
#' @param feature_length Effective length per feature (bases), length nrow(counts).
#' @return Matrix of FPKM values.
#' @export
fpkm <- function(counts, feature_length) {
  counts <- .as_count_matrix(counts)
  feature_length <- .check_feature_length(feature_length, counts)
  lib <- .check_lib_size(colSums(counts))
  per_kb <- counts / (feature_length / 1e3)
  sweep(per_kb, 2L, lib / 1e6, "/")
}

#' Log-counts from raw counts
#'
#' Default `logcounts` assay: `log2(CPM + 1)`.
#'
#' @param counts Numeric matrix of raw counts (features x samples).
#' @return Matrix of log2(CPM + 1) values.
#' @export
logcounts_from_counts <- function(counts) {
  log2(cpm(counts) + 1)
}
