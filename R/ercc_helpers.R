# Spike-in (ERCC) titration QC. The dose-response checks how well observed
# spike-in expression tracks the *known* input concentration per sample (slope
# ~ 1 and high R^2 on a log-log fit = a healthy titration). NOTE: spike-ins are
# QC-only here -- the app normalizes on endogenous controlGenes
# (estimate_size_factors_endogenous), so these are NOT used for normalization.
# Concentrations come from a rowData `spike_concentration` column or the bundled
# ERCC Mix 1 / Mix 2 reference (inst/extdata/ercc_concentrations.csv). Pure +
# tested; the QC module (mod_qc.R) wraps these.

.ercc_cache <- new.env(parent = emptyenv())

#' Bundled ERCC spike-in reference concentrations
#'
#' The 92 ERCC control transcripts with their nominal Mix 1 / Mix 2 input
#' concentrations (attomoles/uL). See `inst/extdata/ERCC_SOURCE.md` for provenance.
#'
#' @return A `data.frame` with `ercc_id`, `subgroup`, `conc_mix1`, `conc_mix2`.
#' @export
ercc_concentrations <- function() {
  if (is.null(.ercc_cache$tbl)) {
    path <- system.file("extdata", "ercc_concentrations.csv", package = "ddsdashboard")
    if (!nzchar(path)) stop("Bundled ERCC concentration table not found.", call. = FALSE)
    .ercc_cache$tbl <- utils::read.csv(path, stringsAsFactors = FALSE)
  }
  .ercc_cache$tbl
}

#' Resolve known concentrations for the spike-in features
#'
#' @param dds A `DESeqDataSet`.
#' @param source `"column"` (use a `rowData` concentration column), `"mix1"` or
#'   `"mix2"` (join the bundled ERCC reference by id).
#' @param conc_col `rowData` column name for the `"column"` source.
#' @param ref The ERCC reference table (defaults to [ercc_concentrations()]).
#' @return A named numeric vector (spike-in id -> concentration), `NA` where
#'   unknown. Empty when the dataset has no spike-in features.
#' @export
resolve_spike_concentration <- function(dds, source = c("column", "mix1", "mix2"),
                                        conc_col = "spike_concentration",
                                        ref = ercc_concentrations()) {
  source <- match.arg(source)
  spike <- rownames(dds)[.detect_spike_features(dds)]
  if (!length(spike)) return(stats::setNames(numeric(0), character(0)))
  if (source == "column") {
    rd <- SummarizedExperiment::rowData(dds)
    v <- if (conc_col %in% colnames(rd)) as.numeric(rd[[conc_col]][match(spike, rownames(dds))])
         else rep(NA_real_, length(spike))
  } else {
    mixcol <- if (source == "mix1") "conc_mix1" else "conc_mix2"
    v <- as.numeric(ref[[mixcol]][match(spike, ref$ercc_id)])
  }
  stats::setNames(v, spike)
}

# Observed expression matrix (linear, depth-normalized) for the dose-response:
# the named assay when present, else CPM computed from counts. Never counts /
# logcounts (the caller restricts choices; this guards the default).
.spike_expr_matrix <- function(dds, assay = "CPM") {
  an <- SummarizedExperiment::assayNames(dds)
  m <- if (assay %in% an && !assay %in% c("counts", "logcounts")) {
    as.matrix(SummarizedExperiment::assay(dds, assay))
  } else {
    cpm(as.matrix(SummarizedExperiment::assay(dds, "counts")))
  }
  colnames(m) <- colnames(dds)
  m
}

# Per-sample log-log fit of observed vs known concentration. Drops zeros before
# logging (no pseudo-count); needs >= 3 usable points or slope/R^2 are NA.
.spike_fit <- function(conc, expr) {
  ok <- is.finite(conc) & is.finite(expr) & conc > 0 & expr > 0
  n <- sum(ok)
  lod <- if (n > 0L) min(conc[ok]) else NA_real_       # lowest *detected* spiked conc
  if (n < 3L) return(list(n_points = n, slope = NA_real_, r_squared = NA_real_, lod = lod))
  fit <- stats::lm(log10(expr[ok]) ~ log10(conc[ok]))
  list(n_points = n, slope = unname(stats::coef(fit)[2]),
       r_squared = summary(fit)$r.squared, lod = lod)
}

#' Spike-in dose-response (per-sample titration QC)
#'
#' For each sample, relates observed spike-in expression to the known input
#' concentration. Returns the long-form points (for the scatter) and a per-sample
#' summary (content + log-log fit). Advisory QC only.
#'
#' @param dds A `DESeqDataSet` with spike-in features.
#' @param assay Observed abundance assay (linear, depth-normalized): `"CPM"`
#'   (default; computed if absent), `"TPM"` or `"FPKM"` (when present).
#' @param source,conc_col,ref Passed to [resolve_spike_concentration()].
#' @return `list(long, per_sample, assay)`. `long`: `sample` (factor), `feature`,
#'   `concentration`, `expression` (one row per spike feature x sample; zeros
#'   kept -- drop them when plotting/fitting). `per_sample`: `sample`, `pct_spike`,
#'   `n_spike_detected`, `n_points`, `slope`, `r_squared`, `lod`.
#' @export
spike_dose_response <- function(dds, assay = "CPM", source = c("column", "mix1", "mix2"),
                                conc_col = "spike_concentration", ref = ercc_concentrations()) {
  source <- match.arg(source)
  spike <- rownames(dds)[.detect_spike_features(dds)]
  samples <- colnames(dds)
  empty <- list(
    long = data.frame(sample = factor(character(0), levels = samples), feature = character(0),
                      concentration = numeric(0), expression = numeric(0),
                      stringsAsFactors = FALSE),
    per_sample = data.frame(sample = samples, pct_spike = NA_real_, n_spike_detected = 0L,
                            n_points = 0L, slope = NA_real_, r_squared = NA_real_, lod = NA_real_,
                            stringsAsFactors = FALSE),
    assay = assay)
  if (!length(spike)) return(empty)

  conc   <- resolve_spike_concentration(dds, source, conc_col, ref)[spike]
  expr   <- .spike_expr_matrix(dds, assay)[spike, , drop = FALSE]
  counts <- as.matrix(SummarizedExperiment::assay(dds, "counts"))
  lib    <- colSums(counts)
  spike_counts <- counts[spike, , drop = FALSE]

  long <- data.frame(
    sample        = factor(rep(samples, each = length(spike)), levels = samples),
    feature       = rep(spike, times = length(samples)),
    concentration = rep(unname(conc), times = length(samples)),
    expression    = as.numeric(expr),
    stringsAsFactors = FALSE)

  per_sample <- do.call(rbind, lapply(seq_along(samples), function(j) {
    f <- .spike_fit(conc, expr[, j])
    data.frame(sample = samples[j],
               pct_spike = if (lib[j] > 0) 100 * sum(spike_counts[, j]) / lib[j] else NA_real_,
               n_spike_detected = sum(spike_counts[, j] > 0),
               n_points = f$n_points, slope = f$slope, r_squared = f$r_squared, lod = f$lod,
               stringsAsFactors = FALSE)
  }))
  list(long = long, per_sample = per_sample, assay = assay)
}
