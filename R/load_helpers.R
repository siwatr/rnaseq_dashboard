# Pure helpers for turning user input into a normalized DESeqDataSet, and for
# the on-load conventions every dataset must satisfy (feature_class, logcounts).
# Kept free of Shiny so they're unit-testable and reusable. Bulk-first: large
# single-cell objects are rejected with a clear "pseudobulk required" message.

# Read a user-uploaded table by file extension (csv/tsv/txt/xlsx/xls). Shared by
# every page that accepts a tabular upload (Input counts/sample sheet, the
# Sample "Additional metadata" sheet, the Gene Sets table import). `col_names`
# follows readr/readxl (TRUE = first row is a header); when FALSE the columns are
# named Column_1.. (friendlier than readr's X1.. / readxl's ...1).
.read_user_table <- function(path, name, col_names = TRUE) {
  ext <- tolower(tools::file_ext(name))
  df <- switch(ext,
    csv  = as.data.frame(readr::read_csv(path, col_names = col_names, show_col_types = FALSE)),
    tsv  = as.data.frame(readr::read_tsv(path, col_names = col_names, show_col_types = FALSE)),
    txt  = as.data.frame(readr::read_tsv(path, col_names = col_names, show_col_types = FALSE)),
    xlsx = as.data.frame(readxl::read_excel(path, col_names = col_names)),
    xls  = as.data.frame(readxl::read_excel(path, col_names = col_names)),
    stop("Unsupported file type '.", ext, "'. Use CSV, TSV, or XLSX.", call. = FALSE)
  )
  if (isFALSE(col_names) && ncol(df)) names(df) <- paste0("Column_", seq_len(ncol(df)))
  df
}

.counts_assay <- function(obj) {
  an <- SummarizedExperiment::assayNames(obj)
  nm <- if ("counts" %in% an) "counts" else if (length(an)) an[[1L]] else NA_character_
  if (is.na(nm)) stop("Input has no assay to use as counts.", call. = FALSE)
  as.matrix(SummarizedExperiment::assay(obj, nm))
}

.build_dds <- function(counts, col_data, row_data = NULL, design = ~ 1) {
  counts <- as.matrix(counts)
  if (!is.numeric(counts)) stop("Count data must be numeric.", call. = FALSE)
  counts <- round(counts)
  storage.mode(counts) <- "integer"
  if (anyNA(counts) || any(counts < 0)) {
    stop("Count data must be non-negative with no missing values.", call. = FALSE)
  }
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts, colData = col_data,
                                        design = design)
  if (!is.null(row_data) && ncol(row_data) > 0L) {
    SummarizedExperiment::rowData(dds) <- row_data
  }
  dds
}

#' Coerce supported input to a DESeqDataSet
#'
#' A `DESeqDataSet` passes through. A `SingleCellExperiment`/`SummarizedExperiment`
#' is converted by per-cell coercion (`design = ~ 1`) — but a `SingleCellExperiment`
#' with more than `max_cells` columns is rejected (pseudobulk is required, and
#' lands in a later version). Anything else errors.
#'
#' @param obj Input object.
#' @param max_cells Cell-count ceiling for per-cell coercion of single-cell input.
#' @return A `DESeqDataSet`.
#' @export
as_input_dds <- function(obj, max_cells = 1000L) {
  if (methods::is(obj, "DESeqDataSet")) return(obj)
  if (methods::is(obj, "SingleCellExperiment") || methods::is(obj, "SummarizedExperiment")) {
    if (methods::is(obj, "SingleCellExperiment") && ncol(obj) > max_cells) {
      stop(sprintf(paste0("This looks like a large single-cell object (%d cells). ",
                          "Per-cell DESeq2 is not supported at this size; pseudobulk ",
                          "aggregation is required and is coming in a later version."),
                  ncol(obj)), call. = FALSE)
    }
    dds <- .build_dds(.counts_assay(obj), SummarizedExperiment::colData(obj), design = ~ 1)
    rd <- SummarizedExperiment::rowData(obj)
    if (ncol(rd) > 0L) SummarizedExperiment::rowData(dds) <- rd
    return(dds)
  }
  stop("Unsupported input: expected a DESeqDataSet or SingleCellExperiment.", call. = FALSE)
}

#' App-level metadata describing the input (before coercion)
#' @param obj The original loaded object.
#' @param max_cells Same ceiling as [as_input_dds()].
#' @return list(`data_type`, `sce_per_cell`).
#' @export
input_meta <- function(obj, max_cells = 1000L) {
  is_sc <- methods::is(obj, "SingleCellExperiment")
  list(
    data_type    = if (is_sc) "single-cell" else "bulk",
    sce_per_cell = is_sc && ncol(obj) <= max_cells
  )
}

#' Ensure the rowData feature_class column exists
#'
#' Adds a `feature_class` factor (`endogenous` default, `spike_in` for `^ERCC-`
#' ids) when absent. Never overwrites an existing column.
#' @param dds A `DESeqDataSet`.
#' @return The `DESeqDataSet` with `rowData$feature_class`.
#' @export
ensure_feature_class <- function(dds) {
  rd <- SummarizedExperiment::rowData(dds)
  if ("feature_class" %in% colnames(rd)) return(dds)
  cls <- rep("endogenous", nrow(dds))
  cls[grepl("^ERCC-", rownames(dds), ignore.case = TRUE)] <- "spike_in"
  rd$feature_class <- factor(cls, levels = c("endogenous", "spike_in", "exogenous"))
  SummarizedExperiment::rowData(dds) <- rd
  dds
}

#' Ensure the logcounts assay exists
#'
#' Adds `logcounts = log2(CPM + 1)` (via [logcounts_from_counts()]) when absent.
#' @param dds A `DESeqDataSet` with a `counts` assay.
#' @return The `DESeqDataSet` with a `logcounts` assay.
#' @export
ensure_logcounts <- function(dds) {
  if ("logcounts" %in% SummarizedExperiment::assayNames(dds)) return(dds)
  counts <- as.matrix(SummarizedExperiment::assay(dds, "counts"))
  SummarizedExperiment::assay(dds, "logcounts") <- logcounts_from_counts(counts)
  dds
}

# Resolve sample ids from a sample sheet: explicit column, then meaningful
# rownames, then the first column.
.sample_ids <- function(sample_df, sample_id_col = NULL) {
  if (!is.null(sample_id_col)) {
    if (!sample_id_col %in% colnames(sample_df)) {
      stop("sample_id_col '", sample_id_col, "' not found in the sample sheet.", call. = FALSE)
    }
    return(as.character(sample_df[[sample_id_col]]))
  }
  rn <- rownames(sample_df)
  default_rn <- is.null(rn) || identical(rn, as.character(seq_len(nrow(sample_df))))
  if (!default_rn) return(as.character(rn))
  as.character(sample_df[[1L]])
}

#' Build a DESeqDataSet from a counts table and a sample sheet
#'
#' @param counts_df Feature-by-sample counts (data.frame/matrix). Feature ids are
#'   taken from `id_col`, else the first non-numeric column, else rownames.
#' @param sample_df Sample sheet (one row per sample).
#' @param id_col Optional feature-id column name in `counts_df`.
#' @param sample_id_col Optional sample-id column name in `sample_df`.
#' @param design Model formula (default `~ 1`).
#' @return A `DESeqDataSet`. Errors (listing mismatches) if sample ids and counts
#'   columns don't match.
#' @export
tabular_to_dds <- function(counts_df, sample_df, id_col = NULL,
                           sample_id_col = NULL, design = ~ 1) {
  counts_df <- as.data.frame(counts_df, check.names = FALSE)
  is_num <- vapply(counts_df, is.numeric, logical(1))
  if (!is.null(id_col)) {
    ids <- as.character(counts_df[[id_col]]); counts_df[[id_col]] <- NULL
  } else if (length(is_num) && !is_num[[1L]]) {
    ids <- as.character(counts_df[[1L]]); counts_df <- counts_df[-1L]
  } else {
    ids <- rownames(counts_df)
  }
  mat <- as.matrix(counts_df)
  if (!is.numeric(mat)) stop("Counts must be numeric after removing the id column.", call. = FALSE)
  rownames(mat) <- ids

  sample_df <- as.data.frame(sample_df, check.names = FALSE)
  sids <- .sample_ids(sample_df, sample_id_col)
  rownames(sample_df) <- sids

  miss_sheet  <- setdiff(colnames(mat), sids)
  miss_counts <- setdiff(sids, colnames(mat))
  if (length(miss_sheet) || length(miss_counts)) {
    msg <- "Sample IDs in the counts matrix and the sample sheet do not match."
    if (length(miss_sheet))  msg <- paste0(msg, "\n  In counts but not the sheet: ",
                                           paste(utils::head(miss_sheet, 10), collapse = ", "))
    if (length(miss_counts)) msg <- paste0(msg, "\n  In the sheet but not counts: ",
                                           paste(utils::head(miss_counts, 10), collapse = ", "))
    stop(msg, call. = FALSE)
  }
  .build_dds(mat, sample_df[colnames(mat), , drop = FALSE], design = design)
}

#' Best-effort guess of the feature unit (gene / transcript / …)
#'
#' Looks for a `<type>_name` column in `rowData`, then Ensembl-style ids.
#' @param dds A `DESeqDataSet`.
#' @return list(`feature_type`, `confident`).
#' @export
detect_feature_type <- function(dds) {
  cols <- colnames(SummarizedExperiment::rowData(dds))
  m <- regmatches(cols, regexec("^([A-Za-z]+)_name$", cols))
  hits <- vapply(m, function(x) if (length(x) == 2L) x[[2L]] else NA_character_, character(1))
  hits <- hits[!is.na(hits)]
  if (length(hits)) return(list(feature_type = hits[[1L]], confident = TRUE))
  ids <- rownames(dds)
  if (length(ids) && any(grepl("^ENS[A-Z]*G[0-9]+", ids))) return(list(feature_type = "gene", confident = TRUE))
  if (length(ids) && any(grepl("^ENS[A-Z]*T[0-9]+", ids))) return(list(feature_type = "transcript", confident = TRUE))
  list(feature_type = "feature", confident = FALSE)
}
