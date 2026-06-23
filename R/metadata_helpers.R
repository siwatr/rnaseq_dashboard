# Pure helpers for editing sample (colData) or feature (rowData) metadata on a
# DESeqDataSet. Slot-agnostic so one draft-editor module can serve both tabs.
# Kept free of Shiny; modules wrap them in state_mutate(). Raw counts untouched.

.meta_get <- function(dds, slot) {
  if (slot == "colData") SummarizedExperiment::colData(dds) else SummarizedExperiment::rowData(dds)
}
.meta_set <- function(dds, slot, df) {
  if (slot == "colData") SummarizedExperiment::colData(dds) <- df
  else SummarizedExperiment::rowData(dds) <- df
  dds
}

# Columns constrained to a fixed value set (validated on edit). feature_class is
# the only one today.
.constrained_levels <- function(slot, col) {
  if (slot == "rowData" && col == "feature_class") c("endogenous", "spike_in", "exogenous") else NULL
}

# Coerce a (string) cell value to the type of an existing column.
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

#' Protected metadata columns
#'
#' Columns that cannot be removed: `colData` design variables (renamable — the
#' formula follows) and the `rowData` `feature_class` (also not renamable).
#' @param dds A `DESeqDataSet`.
#' @param slot `"colData"` or `"rowData"`.
#' @return Character vector of protected column names.
#' @export
protected_columns <- function(dds, slot = "colData") {
  if (slot == "rowData") return(intersect("feature_class", colnames(.meta_get(dds, slot))))
  d <- tryCatch(DESeq2::design(dds), error = function(e) NULL)
  if (is.null(d)) character(0) else all.vars(d)
}

#' Edit one metadata cell
#'
#' Coerces `value` to the column's type. A constrained column (`feature_class`)
#' rejects values outside its allowed set; a free factor column gains a new
#' level for an unseen value.
#' @param dds A `DESeqDataSet`.
#' @param slot `"colData"` (rows = samples) or `"rowData"` (rows = features).
#' @param row Row id (rowname) or 1-based index within the slot.
#' @param col Column name.
#' @param value New value (character; coerced/validated).
#' @return The updated `DESeqDataSet`.
#' @export
edit_meta_cell <- function(dds, slot, row, col, value) {
  df <- .meta_get(dds, slot)
  if (!col %in% colnames(df)) stop("Unknown ", slot, " column: ", col, call. = FALSE)
  ridx <- if (is.numeric(row)) as.integer(row) else match(as.character(row), rownames(df))
  if (length(ridx) != 1L || is.na(ridx) || ridx < 1L || ridx > nrow(df)) {
    stop("Unknown row: ", row, call. = FALSE)
  }
  allowed <- .constrained_levels(slot, col)
  column <- df[[col]]
  if (!is.null(allowed)) {
    if (!as.character(value) %in% allowed) {
      stop(col, " must be one of: ", paste(allowed, collapse = ", "), ".", call. = FALSE)
    }
    if (!is.factor(column)) column <- factor(as.character(column), levels = allowed)
    column[ridx] <- as.character(value)
  } else {
    new <- .coerce_cell(value, column)
    if (is.factor(column)) {
      column <- factor(as.character(column), levels = union(levels(column), as.character(new)))
    }
    column[ridx] <- new
  }
  df[[col]] <- column
  .meta_set(dds, slot, df)
}

#' Add a typed metadata column
#' @inheritParams edit_meta_cell
#' @param name New column name (must not already exist).
#' @param type One of character/numeric/integer/logical/factor.
#' @param default Fill value for every row (default `NA`).
#' @return The updated `DESeqDataSet`.
#' @export
add_meta_column <- function(dds, slot, name,
                            type = c("character", "numeric", "integer", "logical", "factor"),
                            default = NA) {
  type <- match.arg(type)
  df <- .meta_get(dds, slot)
  if (!nzchar(name)) stop("Column name must be non-empty.", call. = FALSE)
  if (name %in% colnames(df)) stop("Column '", name, "' already exists.", call. = FALSE)
  n <- nrow(df)
  df[[name]] <- switch(type,
    character = rep(as.character(default), n),
    numeric   = rep(as.numeric(default), n),
    integer   = rep(as.integer(default), n),
    logical   = rep(as.logical(default), n),
    factor    = factor(rep(as.character(default), n)))
  .meta_set(dds, slot, df)
}

#' Remove metadata columns (multi; protected columns are skipped)
#' @inheritParams edit_meta_cell
#' @param names Character vector of columns to remove.
#' @return list(`dds`, `removed`, `skipped` (protected), `unknown`).
#' @export
remove_meta_columns <- function(dds, slot, names) {
  df <- .meta_get(dds, slot)
  names <- as.character(names)
  unknown   <- setdiff(names, colnames(df))
  skipped   <- intersect(names, protected_columns(dds, slot))
  removable <- setdiff(intersect(names, colnames(df)), skipped)
  for (cn in removable) df[[cn]] <- NULL
  list(dds = .meta_set(dds, slot, df), removed = removable, skipped = skipped, unknown = unknown)
}

#' Rename a metadata column
#'
#' For `colData`, renaming a design variable rewrites the design formula. A
#' protected `rowData` column (`feature_class`) cannot be renamed.
#' @inheritParams edit_meta_cell
#' @param old Existing column name.
#' @param new New column name (non-empty, not already used).
#' @return The updated `DESeqDataSet`.
#' @export
rename_meta_column <- function(dds, slot, old, new) {
  df <- .meta_get(dds, slot)
  if (!old %in% colnames(df)) stop("Unknown ", slot, " column: ", old, call. = FALSE)
  if (!nzchar(new)) stop("New column name must be non-empty.", call. = FALSE)
  if (new %in% colnames(df)) stop("Column '", new, "' already exists.", call. = FALSE)
  if (slot == "rowData" && old %in% protected_columns(dds, slot)) {
    stop("Column '", old, "' is required and cannot be renamed.", call. = FALSE)
  }
  was_design <- slot == "colData" && old %in% protected_columns(dds, "colData")
  colnames(df)[match(old, colnames(df))] <- new
  dds <- .meta_set(dds, slot, df)
  if (was_design) {
    f <- DESeq2::design(dds)
    rhs <- paste(deparse(f[[length(f)]]), collapse = " ")
    rhs <- gsub(paste0("\\b", old, "\\b"), new, rhs)
    DESeq2::design(dds) <- stats::as.formula(paste("~", rhs))
  }
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
#' @param dds A `DESeqDataSet`.
#' @param table A data.frame-like sample sheet (one row per sample).
#' @param id_col Optional sample-id column name in `table`.
#' @return A list with `dds` (updated) and `report`
#'   (`matched`, `unmatched_in_data`, `unmatched_in_table`, `overwritten`).
#' @export
merge_sample_metadata <- function(dds, table, id_col = NULL) {
  table <- as.data.frame(table, check.names = FALSE, stringsAsFactors = FALSE)
  ids <- .table_ids(table, id_col)
  samples <- colnames(dds)
  matched <- intersect(samples, ids)
  if (length(matched) == 0L) {
    stop("No sample IDs in the uploaded table match the dataset.", call. = FALSE)
  }
  newcols <- table
  drop_col <- if (!is.null(id_col)) id_col
              else if (is.null(rownames(table)) ||
                       identical(rownames(table), as.character(seq_len(nrow(table))))) names(table)[1L]
              else NULL
  if (!is.null(drop_col)) newcols[[drop_col]] <- NULL

  aligned <- newcols[match(samples, ids), , drop = FALSE]
  cd <- SummarizedExperiment::colData(dds)
  overwritten <- intersect(colnames(aligned), colnames(cd))
  for (cn in colnames(aligned)) cd[[cn]] <- aligned[[cn]]
  SummarizedExperiment::colData(dds) <- cd

  list(
    dds = dds,
    report = list(
      matched            = length(matched),
      unmatched_in_data  = setdiff(samples, ids),
      unmatched_in_table = setdiff(ids, samples),
      overwritten        = overwritten
    )
  )
}

#' Rename samples (colnames), enforcing uniqueness
#' @param dds A `DESeqDataSet`.
#' @param old Existing sample name(s).
#' @param new Replacement name(s), same length as `old`.
#' @return The updated `DESeqDataSet`.
#' @export
rename_samples <- function(dds, old, new) {
  if (length(old) != length(new)) stop("`old` and `new` must be the same length.", call. = FALSE)
  cn <- colnames(dds)
  idx <- match(as.character(old), cn)
  if (anyNA(idx)) stop("Unknown sample(s): ", paste(old[is.na(idx)], collapse = ", "), call. = FALSE)
  if (any(!nzchar(new))) stop("Sample names must be non-empty.", call. = FALSE)
  cn[idx] <- as.character(new)
  if (anyDuplicated(cn)) stop("Sample names must be unique.", call. = FALSE)
  colnames(dds) <- cn
  dds
}

#' Set spike-in concentrations from an existing numeric column
#'
#' Copies `from_col` into the standard `spike_concentration` rowData column,
#' keeping a value only for spike-in features (`NA` elsewhere). Used by the
#' Feature tab to designate which column holds the known ERCC dose.
#' @param dds A `DESeqDataSet`.
#' @param from_col A numeric `rowData` column name.
#' @return The updated `DESeqDataSet`.
#' @export
set_spike_concentration <- function(dds, from_col) {
  rd <- SummarizedExperiment::rowData(dds)
  if (!from_col %in% colnames(rd)) stop("Unknown rowData column: ", from_col, call. = FALSE)
  v <- suppressWarnings(as.numeric(rd[[from_col]]))
  v[!.detect_spike_features(dds)] <- NA_real_     # only spike-ins carry a concentration
  rd$spike_concentration <- v
  SummarizedExperiment::rowData(dds) <- rd
  dds
}

#' Spike-in features lacking a usable concentration
#'
#' @param dds A `DESeqDataSet`.
#' @param conc_col Concentration `rowData` column (default `spike_concentration`).
#' @return Character vector of spike-in feature ids with `NA`/non-positive
#'   concentration (empty when all are usable or there are no spike-ins).
#' @export
spike_features_missing_conc <- function(dds, conc_col = "spike_concentration") {
  spike <- .detect_spike_features(dds)
  if (!any(spike)) return(character(0))
  rd <- SummarizedExperiment::rowData(dds)
  v <- if (conc_col %in% colnames(rd)) suppressWarnings(as.numeric(rd[[conc_col]])) else rep(NA_real_, nrow(dds))
  rownames(dds)[spike & (is.na(v) | v <= 0)]
}

#' Reset one metadata slot to the originally loaded values
#'
#' Slot-scoped "reset to original": restores `colData` or `rowData` of `working`
#' to the values in `original`, **by id, for the rows currently present** in
#' `working` — so it does not undo sample/feature *filtering*. The original
#' column schema wins (user-added columns are dropped). For `colData` the design
#' formula is restored from `original` (so it cannot reference a column that no
#' longer exists); for `rowData`, normalized assays are refreshed and size
#' factors re-estimated when set (feature_length / feature_class feed TPM/FPKM
#' and `controlGenes`). Rows whose id is absent from `original` (e.g. a renamed
#' sample) keep their current values for the original-schema columns.
#'
#' @param working The current `DESeqDataSet`.
#' @param original The originally loaded `DESeqDataSet` (the reset target).
#' @param slot `"colData"` (samples) or `"rowData"` (features).
#' @return `working` with that slot reverted to `original`.
#' @export
reset_metadata_slot <- function(working, original, slot = c("colData", "rowData")) {
  slot <- match.arg(slot)
  ids  <- if (slot == "colData") colnames(working) else rownames(working)
  orig <- .meta_get(original, slot)
  cur  <- .meta_get(working, slot)
  ocols <- colnames(orig)
  pos  <- match(ids, rownames(orig))                  # NA for ids not in original
  out  <- orig[pos, ocols, drop = FALSE]              # original values, original schema
  rownames(out) <- ids
  miss <- which(is.na(pos))                            # renamed / not-in-original rows
  if (length(miss)) {
    for (j in intersect(ocols, colnames(cur))) out[[j]] <- .fill_rows(out[[j]], miss, cur[[j]][miss])
  }
  working <- .meta_set(working, slot, out)
  if (slot == "colData") {
    working <- tryCatch({ DESeq2::design(working) <- DESeq2::design(original); working },
                        error = function(e) working)
  } else {
    working <- refresh_assays(working)
    if (!is.null(DESeq2::sizeFactors(working))) working <- estimate_size_factors_endogenous(working)
  }
  working
}

#' Set the feature_class of selected features (bulk)
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
