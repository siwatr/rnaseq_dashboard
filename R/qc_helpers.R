# Pure helpers behind the QC page (Sample QC tab): per-sample QC metrics and the
# feature-subset detection they need. Metrics are computed on raw counts. scater
# (a Suggests package) is used when available; a base-R path returns the same
# schema so the table still renders without it. See the rnaseq-bioc skill.

# First rowData "<...>_name" column (e.g. gene_name), else the rownames.
.feature_names <- function(dds) {
  rd <- SummarizedExperiment::rowData(dds)
  nm_cols <- grep("_name$", colnames(rd), value = TRUE)
  if (length(nm_cols)) as.character(rd[[nm_cols[1]]]) else rownames(dds)
}

# A chromosome-like vector from rowData (chromosome/seqnames/chr) when present.
.feature_chrom <- function(dds) {
  rd <- SummarizedExperiment::rowData(dds)
  for (col in c("chromosome", "seqnames", "chr")) {
    if (col %in% colnames(rd)) return(as.character(rd[[col]]))
  }
  NULL
}

#' Flag mitochondrial features
#'
#' Identifies mito rows from a chromosome-like `rowData` column
#' (`MT`/`chrM`/`chrMT`/`M`, case-insensitive) or a `<...>_name` column matching
#' `^MT-` (covers human `MT-` and mouse `mt-`).
#' @param dds A `DESeqDataSet` / `SummarizedExperiment`.
#' @return A logical vector of length `nrow(dds)`.
#' @export
detect_mito_features <- function(dds) {
  n <- nrow(dds)
  is_mito <- logical(n)
  chrom <- .feature_chrom(dds)
  if (!is.null(chrom) && length(chrom) == n) {
    is_mito <- is_mito | toupper(chrom) %in% c("MT", "CHRM", "CHRMT", "M")
  }
  nm <- .feature_names(dds)
  if (length(nm) == n) is_mito <- is_mito | grepl("^MT-", nm, ignore.case = TRUE)
  is_mito
}

# Rows tagged spike-in via the always-present feature_class column.
.detect_spike_features <- function(dds) {
  rd <- SummarizedExperiment::rowData(dds)
  if (!"feature_class" %in% colnames(rd)) return(logical(nrow(dds)))
  as.character(rd$feature_class) == "spike_in"
}

# scater path: percentages come straight from perCellQCMetrics subsets.
.qc_metrics_scater <- function(dds, is_mito, is_spike) {
  n <- ncol(dds)
  subsets <- list()
  if (any(is_mito))  subsets$mito  <- which(is_mito)
  if (any(is_spike)) subsets$spike <- which(is_spike)
  qc <- scater::perCellQCMetrics(dds, subsets = subsets)
  list(library_size = qc$sum,
       detected     = qc$detected,
       pct_mito     = qc$subsets_mito_percent  %||% rep(0, n),
       pct_spike    = qc$subsets_spike_percent %||% rep(0, n))
}

# Base-R path: identical schema, no scater dependency.
.qc_metrics_base <- function(dds, is_mito, is_spike) {
  counts <- as.matrix(SummarizedExperiment::assay(dds, "counts"))
  lib <- colSums(counts)
  n <- ncol(dds)
  pct_of <- function(mask) {
    if (!any(mask)) return(rep(0, n))
    ifelse(lib > 0, 100 * colSums(counts[mask, , drop = FALSE]) / lib, 0)
  }
  list(library_size = lib,
       detected     = colSums(counts > 0),
       pct_mito     = pct_of(is_mito),
       pct_spike    = pct_of(is_spike))
}

#' Per-sample QC metrics
#'
#' Computes, for each sample (column): library size (total counts), number of
#' detected features (count > 0), and the percentage of counts in mitochondrial
#' and in spike-in features. Uses [scater::perCellQCMetrics()] when `scater` is
#' installed; otherwise a base-R fallback returning the identical schema.
#'
#' @param dds A `DESeqDataSet` with a `"counts"` assay.
#' @return A `data.frame`, one row per sample (row names = `colnames(dds)`), with
#'   columns `sample`, `library_size`, `detected`, `pct_mito`, `pct_spike`.
#' @export
qc_per_sample_metrics <- function(dds) {
  samples  <- colnames(dds)
  is_mito  <- detect_mito_features(dds)
  is_spike <- .detect_spike_features(dds)
  vals <- if (requireNamespace("scater", quietly = TRUE)) {
    .qc_metrics_scater(dds, is_mito, is_spike)
  } else {
    .qc_metrics_base(dds, is_mito, is_spike)
  }
  data.frame(
    sample       = samples,
    library_size = as.numeric(vals$library_size),
    detected     = as.integer(vals$detected),
    pct_mito     = as.numeric(vals$pct_mito),
    pct_spike    = as.numeric(vals$pct_spike),
    row.names    = samples,
    stringsAsFactors = FALSE,
    check.names  = FALSE
  )
}

# ---- Dataset-level diagnostics (P3b) ---------------------------------------

# A samples-by-? value matrix: the named assay when present, else log2(counts+1).
# Returns features x samples with colnames = sample ids.
.qc_assay_matrix <- function(dds, assay = "logcounts") {
  m <- if (assay %in% SummarizedExperiment::assayNames(dds)) {
    SummarizedExperiment::assay(dds, assay)
  } else {
    log2(as.matrix(SummarizedExperiment::assay(dds, "counts")) + 1)
  }
  m <- as.matrix(m)
  colnames(m) <- colnames(dds)
  m
}

#' Variance-stabilizing transform for QC
#'
#' Returns a [DESeq2::DESeqTransform] for the mean-SD diagnostic. Uses the fast
#' [DESeq2::vst()] approximation when the dataset is large enough, falling back
#' to [DESeq2::varianceStabilizingTransformation()] (which works on small or
#' low-count data, where `vst()` errors on its `nsub` requirement).
#'
#' @param dds A `DESeqDataSet`.
#' @param blind Estimate dispersions blind to the design (default `TRUE` for QC).
#' @return A `DESeqTransform`.
#' @export
qc_vst <- function(dds, blind = TRUE) {
  tryCatch(
    DESeq2::vst(dds, blind = blind),
    error = function(e) DESeq2::varianceStabilizingTransformation(dds, blind = blind)
  )
}

#' Sample-to-sample correlation matrix
#'
#' Correlation between samples on a value assay (default `logcounts`, falling
#' back to `log2(counts + 1)`). Drives the sample-correlation heatmap.
#'
#' @param dds A `DESeqDataSet`.
#' @param method Correlation method, `"spearman"` (default) or `"pearson"`.
#' @param assay Assay to correlate on; falls back to `log2(counts + 1)` if absent.
#' @return A symmetric samples-by-samples correlation matrix.
#' @export
qc_sample_correlation <- function(dds, method = c("spearman", "pearson"),
                                  assay = "logcounts") {
  method <- match.arg(method)
  stats::cor(.qc_assay_matrix(dds, assay), method = method)
}

#' Relative log expression (RLE) matrix
#'
#' Each value is its assay value minus that feature's median across samples; a
#' well-normalized sample has RLE values centered on 0 with small spread. All-NA
#' / zero-variance rows contribute 0.
#'
#' @param dds A `DESeqDataSet`.
#' @param assay Assay to use; falls back to `log2(counts + 1)` if absent.
#' @return A features-by-samples matrix of median-centered values.
#' @export
qc_rle_matrix <- function(dds, assay = "logcounts") {
  m <- .qc_assay_matrix(dds, assay)
  med <- if (requireNamespace("matrixStats", quietly = TRUE)) {
    matrixStats::rowMedians(m)
  } else {
    apply(m, 1, stats::median)
  }
  m - med
}

#' Per-sample expression values in long form
#'
#' Long `data.frame` of the value assay (default `logcounts`) for the per-sample
#' expression-density diagnostic.
#'
#' @param dds A `DESeqDataSet`.
#' @param assay Assay to use; falls back to `log2(counts + 1)` if absent.
#' @return A `data.frame` with columns `sample` (factor) and `value`, one row per
#'   feature x sample.
#' @export
qc_expression_long <- function(dds, assay = "logcounts") {
  m <- .qc_assay_matrix(dds, assay)
  data.frame(
    sample = factor(rep(colnames(m), each = nrow(m)), levels = colnames(m)),
    value  = as.numeric(m),
    stringsAsFactors = FALSE
  )
}
