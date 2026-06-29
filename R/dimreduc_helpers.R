# Phase 4: dimensionality reduction (PCA). Pure, exported helpers the
# mod_dimreduc page calls -- choosing/transforming the input matrix, ranking the
# top-variable genes, and computing the PCA embedding (cached at the call site via
# state_derive). Endogenous-only throughout (spike-in/exogenous excluded, mirroring
# the size-factor + variable-gene conventions). Matches DESeq2::plotPCA semantics:
# top-ntop by variance on the (stabilized) matrix -> prcomp(t(mat), center=TRUE,
# scale.=FALSE) -> %variance from sdev^2.

#' PCA-input advice for an assay
#'
#' Guides which assays are sensible PCA inputs. Linear abundance assays
#' (CPM/TPM/FPKM) let a few high-abundance genes dominate PC1 (mean-variance
#' dependence) and should be log-transformed first; raw `counts` is unsuitable
#' (not depth-normalized - a log won't fix it). VST / logcounts / normalized
#' log-counts are ready to use.
#'
#' @param assay Assay/transform name (`"vst"`, `"logcounts"`, `"norm_logcounts"`,
#'   `"CPM"`, `"TPM"`, `"FPKM"`, `"counts"`, or any stored assay).
#' @return A list: `tier` (`"recommended"`/`"log_first"`/`"unsuitable"`),
#'   `recommend_log` (logical, the default log-transform state), and `msg` (a
#'   user-facing warning, or `NULL`).
#' @export
pca_assay_advice <- function(assay) {
  rec <- list(tier = "recommended", recommend_log = FALSE, msg = NULL)
  switch(assay,
    vst            = rec,
    logcounts      = rec,
    norm_logcounts = rec,
    CPM = , TPM = , FPKM = list(
      tier = "log_first", recommend_log = TRUE,
      msg = sprintf(paste0("%s is a linear abundance measure. Running PCA on it lets a ",
        "few highly expressed genes dominate (mean-variance dependence) rather than ",
        "reflecting biological structure. VST or log-counts are recommended; applying ",
        "log2(x+1) before PCA helps."), assay)),
    counts = list(
      tier = "unsuitable", recommend_log = FALSE,
      msg = paste0("counts is not normalized for sequencing depth and is unsuitable for ",
        "PCA - the first components will mostly reflect library-size differences, and a ",
        "log transform does not fix this. Use VST or log-counts instead.")),
    list(tier = "log_first", recommend_log = TRUE,
      msg = sprintf(paste0("'%s' is not a standard PCA input; VST or log-counts are ",
        "recommended. Applying log2(x+1) may help if it is a linear measure."), assay)))
}

#' Build the PCA input matrix for an assay (endogenous-only) + an honest label
#'
#' @param dds A `DESeqDataSet`.
#' @param assay One of `"vst"` (default), `"norm_logcounts"`, or a stored assay
#'   name (e.g. `"logcounts"`, `"CPM"`).
#' @param log_transform Apply `log2(x + 1)` to a *stored* assay (ignored for the
#'   computed `"vst"` / `"norm_logcounts"` inputs, which are already on a log/
#'   stabilized scale).
#' @return A list with `mat` (genesxsamples numeric, endogenous rows only) and
#'   `label` (for the plot subtitle, reflecting any fallback/transform).
#' @export
pca_input <- function(dds, assay = "vst", log_transform = FALSE) {
  endo <- .qc_diagnostic_rows(dds)
  fallback_logcounts <- function(why)
    list(mat = .qc_assay_matrix(dds, "logcounts")[endo, , drop = FALSE],
         label = paste0("logcounts (", why, ")"))

  if (identical(assay, "vst")) {
    # qc_vst() already subsets to endogenous -> do NOT re-subset. Broaden the
    # fallback to *any* VST failure (small/low-count/zero-library data).
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
      list(mat = log2(nc + 1)[.qc_diagnostic_rows(d), , drop = FALSE],
           label = "normalized log-counts (log2)")
    }, error = function(e) NULL)
    return(res %||% fallback_logcounts("normalized counts unavailable"))
  }

  m <- .qc_assay_matrix(dds, assay)[endo, , drop = FALSE]
  if (isTRUE(log_transform)) return(list(mat = log2(m + 1), label = paste0(assay, " (log2)")))
  list(mat = m, label = assay)
}

#' Top-variable feature ids from a value matrix
#'
#' Per-row variance (drops non-finite / zero-variance rows), then the most
#' variable `n_top` (clamped to what's available).
#' @param mat A genesxsamples numeric matrix (already endogenous-only).
#' @param n_top Number of features to keep (default 500).
#' @return Character vector of row ids (possibly fewer than `n_top`).
#' @export
top_variable_features <- function(mat, n_top = 500) {
  v <- if (requireNamespace("matrixStats", quietly = TRUE)) {
    matrixStats::rowVars(as.matrix(mat))
  } else {
    apply(mat, 1, stats::var)
  }
  names(v) <- rownames(mat)
  v <- v[is.finite(v) & v > 0]
  if (!length(v)) return(character(0))
  n <- min(as.integer(n_top), length(v))
  names(sort(v, decreasing = TRUE))[seq_len(n)]
}

#' Compute a PCA embedding from a value matrix
#'
#' Selects the top-variable features, runs `prcomp(t(mat), center = TRUE,
#' scale. = FALSE)` (samples as observations), and applies a deterministic sign
#' convention (each PC's largest-magnitude loading made positive) so re-renders
#' don't mirror-flip.
#' @param mat A genesxsamples numeric matrix (endogenous-only).
#' @param n_top Number of top-variable features (default 500).
#' @return A list: `scores` (data.frame, rows = samples, cols `PC1..PCk`),
#'   `var_pct` (numeric, % variance per PC), `n_genes`, `n_pc`.
#' @export
compute_pca <- function(mat, n_top = 500) {
  top <- top_variable_features(mat, n_top)
  if (length(top) < 2L)
    stop("Not enough variable features for PCA.", call. = FALSE)
  sub <- as.matrix(mat[top, , drop = FALSE])
  pca <- stats::prcomp(t(sub), center = TRUE, scale. = FALSE)
  flip <- apply(pca$rotation, 2, function(col) {
    s <- sign(col[which.max(abs(col))]); if (s == 0) 1 else s
  })
  scores <- sweep(pca$x, 2, flip, `*`)
  var_pct <- pca$sdev^2 / sum(pca$sdev^2) * 100
  scores_df <- as.data.frame(scores)
  rownames(scores_df) <- colnames(sub)
  list(scores = scores_df, var_pct = var_pct, n_genes = length(top), n_pc = ncol(scores))
}
