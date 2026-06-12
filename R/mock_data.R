# Synthetic test/demo fixtures. Deterministic (seeded) so tests are stable and
# the app has data to load without a real dataset. DESeq2/SummarizedExperiment
# are Suggests, so these functions check for them and error cleanly if absent.

#' Build a mock DESeqDataSet for tests and demos
#'
#' Produces a small bulk RNA-seq `DESeqDataSet` exercising every convention the
#' app relies on: a two-level `condition` design with replicates, a pseudobulk
#' grouping column, endogenous genes (some differentially expressed, a few
#' mitochondrial), ERCC spike-ins with known concentrations, and one exogenous
#' gene — all reflected in `rowData$feature_class`, `feature_length`, and a
#' `gene_name` column.
#'
#' @param n_genes Number of endogenous genes.
#' @param n_per_group Replicates per condition (two conditions are created).
#' @param n_spike Number of ERCC spike-in features.
#' @param frac_de Fraction of endogenous genes that are differentially expressed.
#' @param seed RNG seed for reproducibility.
#' @return A `DESeqDataSet` with a raw `counts` assay and `design = ~ condition`.
#' @export
make_mock_dds <- function(n_genes = 200, n_per_group = 4, n_spike = 10,
                          frac_de = 0.1, seed = 1) {
  for (pkg in c("DESeq2", "SummarizedExperiment", "S4Vectors")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("make_mock_dds() needs the '", pkg, "' package.", call. = FALSE)
    }
  }
  withr_seed(seed)

  ## --- samples (colData) ---
  condition <- factor(rep(c("control", "treated"), each = n_per_group),
                      levels = c("control", "treated"))
  n_samples <- length(condition)
  sample_ids <- sprintf("S%02d", seq_len(n_samples))
  col_data <- S4Vectors::DataFrame(
    condition = condition,
    bio_rep   = factor(rep(seq_len(n_per_group), times = 2)),
    group     = condition,                 # auto-suggested pseudobulk grouping
    row.names = sample_ids
  )

  ## --- endogenous genes ---
  base_mu <- 2^stats::runif(n_genes, 1, 12)            # wide dynamic range
  lfc <- numeric(n_genes)
  de_idx <- sample.int(n_genes, max(1L, round(frac_de * n_genes)))
  lfc[de_idx] <- stats::rnorm(length(de_idx), 0, 2)
  mu <- matrix(base_mu, nrow = n_genes, ncol = n_samples)
  treated <- condition == "treated"
  mu[, treated] <- mu[, treated] * 2^lfc
  counts_endo <- matrix(
    stats::rnbinom(n_genes * n_samples, mu = as.vector(mu), size = 1 / 0.2),
    nrow = n_genes
  )
  n_mito <- min(5L, n_genes)
  endo_chr <- rep("1", n_genes); endo_chr[seq_len(n_mito)] <- "MT"
  endo_name <- c(sprintf("mt-Gene%d", seq_len(n_mito)),
                 sprintf("Gene%d", seq_len(n_genes - n_mito)))
  endo <- list(
    counts    = counts_endo,
    id        = sprintf("ENSMUSG%011d", seq_len(n_genes)),
    gene_name = endo_name,
    class     = rep("endogenous", n_genes),
    length    = round(2^stats::runif(n_genes, 8, 12)),     # ~256–4096 bp
    chr       = endo_chr,
    conc      = rep(NA_real_, n_genes)
  )

  ## --- ERCC spike-ins (counts scale with known concentration) ---
  spike <- NULL
  if (n_spike > 0) {
    conc <- 2^stats::runif(n_spike, -2, 14)              # attomoles/uL, broad range
    spike_mu <- matrix(pmax(conc * 0.5, 1), nrow = n_spike, ncol = n_samples)
    spike <- list(
      counts    = matrix(stats::rnbinom(n_spike * n_samples,
                                        mu = as.vector(spike_mu), size = 1 / 0.1),
                         nrow = n_spike),
      id        = sprintf("ERCC-%05d", seq_len(n_spike)),
      gene_name = sprintf("ERCC-%05d", seq_len(n_spike)),
      class     = rep("spike_in", n_spike),
      length    = round(stats::runif(n_spike, 250, 2000)),
      chr       = rep("ERCC", n_spike),
      conc      = conc
    )
  }

  ## --- one exogenous gene (e.g. an overexpression construct) ---
  exo <- list(
    counts    = matrix(stats::rnbinom(n_samples, mu = 5000, size = 1 / 0.1), nrow = 1),
    id        = "EXO-GFP",
    gene_name = "GFP",
    class     = "exogenous",
    length    = 720L,
    chr       = "exogenous",
    conc      = NA_real_
  )

  parts <- Filter(Negate(is.null), list(endo, spike, exo))
  counts <- do.call(rbind, lapply(parts, `[[`, "counts"))
  storage.mode(counts) <- "integer"
  ids <- unlist(lapply(parts, `[[`, "id"))
  colnames(counts) <- sample_ids
  rownames(counts) <- ids

  row_data <- S4Vectors::DataFrame(
    gene_name      = unlist(lapply(parts, `[[`, "gene_name")),
    feature_class  = factor(unlist(lapply(parts, `[[`, "class")),
                            levels = c("endogenous", "spike_in", "exogenous")),
    feature_length = unlist(lapply(parts, `[[`, "length")),
    chromosome     = unlist(lapply(parts, `[[`, "chr")),
    spike_concentration = unlist(lapply(parts, `[[`, "conc")),
    row.names      = ids
  )

  DESeq2::DESeqDataSetFromMatrix(
    countData = counts,
    colData   = col_data,
    rowData   = row_data,
    design    = ~ condition
  )
}

# Seed without taking a hard dependency on withr.
withr_seed <- function(seed) set.seed(seed)
