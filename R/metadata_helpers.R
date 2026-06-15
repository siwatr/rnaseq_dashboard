# Pure helpers for editing sample/feature metadata on a DESeqDataSet. Kept free
# of Shiny so they're unit-testable; the Sample-info module wraps them in
# state_mutate(). Raw counts are never touched.

# Coerce a (string) cell value to the type of an existing column. Errors when a
# numeric/integer/logical target can't accept the value.
.coerce_cell <- function(value, template) {
  value <- as.character(value)
  na_like <- value %in% c("NA", "")
  if (is.factor(template) || is.character(template)) return(value)
  if (is.logical(template)) {
    v <- as.logical(value)
    if (is.na(v) && !na_like) stop("Cannot coerce '", value, "' to logical.", call. = FALSE)
    return(v)
  }
  if (is.integer(template)) {
    v <- suppressWarnings(as.integer(value))
    if (is.na(v) && !na_like) stop("Cannot coerce '", value, "' to integer.", call. = FALSE)
    return(v)
  }
  if (is.numeric(template)) {
    v <- suppressWarnings(as.numeric(value))
    if (is.na(v) && !na_like) stop("Cannot coerce '", value, "' to a number.", call. = FALSE)
    return(v)
  }
  value
}

#' Edit one colData cell
#'
#' Sets `colData(dds)[row, col]` to `value`, coercing to the column's existing
#' type. For factor columns an unseen value is added as a new level.
#'
#' @param dds A `DESeqDataSet`.
#' @param row Sample id (rowname) or 1-based row index.
#' @param col Column name in `colData`.
#' @param value New value (character; coerced to the column type).
#' @return The updated `DESeqDataSet`.
#' @export
edit_coldata_cell <- function(dds, row, col, value) {
  cd <- SummarizedExperiment::colData(dds)
  if (!col %in% colnames(cd)) stop("Unknown colData column: ", col, call. = FALSE)
  ridx <- if (is.numeric(row)) as.integer(row) else match(as.character(row), rownames(cd))
  if (length(ridx) != 1L || is.na(ridx) || ridx < 1L || ridx > nrow(cd)) {
    stop("Unknown sample row: ", row, call. = FALSE)
  }
  column <- cd[[col]]
  new <- .coerce_cell(value, column)
  if (is.factor(column)) {
    column <- factor(as.character(column),
                     levels = union(levels(column), as.character(new)))
  }
  column[ridx] <- new
  cd[[col]] <- column
  SummarizedExperiment::colData(dds) <- cd
  dds
}

# Resolve the sample-id vector of an uploaded table.
.table_ids <- function(table, id_col = NULL) {
  if (!is.null(id_col)) {
    if (!id_col %in% colnames(table)) {
      stop("id_col '", id_col, "' not found in the uploaded table.", call. = FALSE)
    }
    return(as.character(table[[id_col]]))
  }
  rn <- rownames(table)
  if (is.null(rn) || identical(rn, as.character(seq_len(nrow(table))))) {
    return(as.character(table[[1L]]))
  }
  as.character(rn)
}

#' Merge an uploaded sample sheet into colData
#'
#' Joins `table` into `colData(dds)` by sample id (an explicit `id_col`, else
#' meaningful rownames, else the first column). Samples absent from the upload
#' get `NA`; existing columns of the same name are overwritten.
#'
#' @param dds A `DESeqDataSet`.
#' @param table A data.frame-like sample sheet (one row per sample).
#' @param id_col Optional sample-id column name in `table`.
#' @return A list with `dds` (updated) and `report`
#'   (`matched`, `unmatched_in_data`, `unmatched_in_table`).
#' @export
merge_sample_metadata <- function(dds, table, id_col = NULL) {
  table <- as.data.frame(table, check.names = FALSE, stringsAsFactors = FALSE)
  ids <- .table_ids(table, id_col)
  samples <- colnames(dds)
  matched <- intersect(samples, ids)
  if (length(matched) == 0L) {
    stop("No sample IDs in the uploaded table match the dataset.", call. = FALSE)
  }
  # Columns to merge: everything except the id-bearing column.
  newcols <- table
  drop_col <- if (!is.null(id_col)) id_col
              else if (is.null(rownames(table)) ||
                       identical(rownames(table), as.character(seq_len(nrow(table))))) names(table)[1L]
              else NULL
  if (!is.null(drop_col)) newcols[[drop_col]] <- NULL

  aligned <- newcols[match(samples, ids), , drop = FALSE]
  cd <- SummarizedExperiment::colData(dds)
  for (cn in colnames(aligned)) cd[[cn]] <- aligned[[cn]]
  SummarizedExperiment::colData(dds) <- cd

  list(
    dds = dds,
    report = list(
      matched            = length(matched),
      unmatched_in_data  = setdiff(samples, ids),
      unmatched_in_table = setdiff(ids, samples)
    )
  )
}

#' Set the feature_class of selected features
#'
#' @param dds A `DESeqDataSet`.
#' @param ids Feature ids (rownames) to tag.
#' @param class One of `endogenous` / `spike_in` / `exogenous`.
#' @return The updated `DESeqDataSet` (with `feature_class` ensured).
#' @export
set_feature_class <- function(dds, ids, class = c("endogenous", "spike_in", "exogenous")) {
  class <- match.arg(class)
  dds <- ensure_feature_class(dds)
  idx <- match(as.character(ids), rownames(dds))
  if (all(is.na(idx))) stop("None of the given feature ids were found.", call. = FALSE)
  rd <- SummarizedExperiment::rowData(dds)
  fc <- as.character(rd$feature_class)
  fc[idx[!is.na(idx)]] <- class
  rd$feature_class <- factor(fc, levels = c("endogenous", "spike_in", "exogenous"))
  SummarizedExperiment::rowData(dds) <- rd
  dds
}
