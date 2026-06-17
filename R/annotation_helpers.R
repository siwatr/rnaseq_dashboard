# Pure helpers for OrgDb annotation: detect the feature id type, map ids to
# gene symbol / chromosome / description via a Bioconductor OrgDb, and report
# coverage. Org packages are Suggests (loaded on demand). GTF annotation and
# ERCC concentrations come in later PRs.

.orgdb_pkg <- c(mouse = "org.Mm.eg.db", human = "org.Hs.eg.db")

# OrgDb keytype columns offered for import, and the rowData column each becomes.
# SYMBOL/GENENAME keep their established homes; the rest fall back to a lowercase
# name. Used by both the annotate and preview paths so they stay in lockstep.
.orgdb_targets <- function(feature_type = "gene") {
  c(SYMBOL   = paste0(feature_type, "_name"),
    GENENAME = "description",
    ENSEMBL  = "ensembl_id",
    ENTREZID = "entrez_id",
    GENETYPE = "gene_biotype")
}

.orgdb_target <- function(col, feature_type = "gene") {
  m <- .orgdb_targets(feature_type)
  if (col %in% names(m)) unname(m[[col]]) else tolower(col)
}

#' Guess the feature id type
#'
#' @param ids Character vector of feature ids (rownames).
#' @return One of `"ensembl"`, `"entrez"`, `"symbol"` (majority vote).
#' @export
detect_id_type <- function(ids) {
  ids <- ids[!is.na(ids)]
  if (!length(ids)) return("symbol")
  if (mean(grepl("^ENS", ids)) > 0.5) return("ensembl")
  if (mean(grepl("^[0-9]+$", ids)) > 0.5) return("entrez")
  "symbol"
}

# Overlay mapped values onto existing ones: use the mapping where it found a
# hit, otherwise keep whatever was already there (never wipe).
.fill <- function(existing, mapped) {
  out <- if (is.null(existing)) rep(NA_character_, length(mapped)) else as.character(existing)
  ok <- !is.na(mapped)
  out[ok] <- as.character(mapped[ok])
  out
}

# Resolve an OrgDb, returning the package object; errors if the org pkg is absent.
.orgdb_org <- function(organism, who) {
  pkg <- .orgdb_pkg[[organism]]
  for (p in c("AnnotationDbi", pkg)) {
    if (!requireNamespace(p, quietly = TRUE)) {
      stop(who, "() needs the '", p, "' package.", call. = FALSE)
    }
  }
  getExportedValue(pkg, pkg)
}

# mapIds keytype for an id-type string.
.orgdb_keytype <- function(id_type) {
  switch(id_type, ensembl = "ENSEMBL", entrez = "ENTREZID", symbol = "SYMBOL",
         stop("Unknown id_type: ", id_type, call. = FALSE))
}

#' Annotate features from a Bioconductor OrgDb
#'
#' Maps each requested OrgDb column onto `rowData`, filling where a mapping is
#' found and keeping existing values otherwise; `feature_class` is left
#' untouched. Each OrgDb column lands in a fixed `rowData` column: `SYMBOL` ->
#' `<feature_type>_name`, `GENENAME` -> `description`, `ENSEMBL`/`ENTREZID` ->
#' `ensembl_id`/`entrez_id`, `GENETYPE` -> `gene_biotype`. Columns not present in
#' the OrgDb are skipped. Chromosome is intentionally not taken from OrgDb (its
#' `CHR` accessor is deprecated) — it comes from the GTF instead.
#'
#' @param dds A `DESeqDataSet`.
#' @param organism `"mouse"` (org.Mm.eg.db) or `"human"` (org.Hs.eg.db).
#' @param id_type `"ensembl"`/`"entrez"`/`"symbol"`; auto-detected when `NULL`.
#' @param feature_type Feature unit; sets the name column `<feature_type>_name`.
#' @param columns OrgDb keytype columns to import (default `c("SYMBOL",
#'   "GENENAME")`, reproducing the historical name + description behaviour).
#' @param matched_col Optional name of a logical `rowData` column to write,
#'   flagging which features were mapped in the OrgDb (`NULL` = none). The flag
#'   tracks `SYMBOL` when imported, else the first requested column.
#' @return The annotated `DESeqDataSet`.
#' @export
annotate_with_orgdb <- function(dds, organism = c("mouse", "human"),
                                id_type = NULL, feature_type = "gene",
                                columns = c("SYMBOL", "GENENAME"),
                                matched_col = NULL) {
  organism <- match.arg(organism)
  org <- .orgdb_org(organism, "annotate_with_orgdb")
  columns <- intersect(columns, AnnotationDbi::columns(org))
  ids <- rownames(dds)
  id_type <- id_type %||% detect_id_type(ids)
  keytype <- .orgdb_keytype(id_type)
  keys <- if (keytype == "ENSEMBL") sub("\\..*$", "", ids) else ids  # drop version suffix
  pull <- function(col) {
    suppressMessages(tryCatch(
      AnnotationDbi::mapIds(org, keys = keys, column = col, keytype = keytype,
                            multiVals = "first"),
      error = function(e) stats::setNames(rep(NA_character_, length(keys)), keys)
    ))
  }
  rd <- SummarizedExperiment::rowData(dds)
  primary <- if ("SYMBOL" %in% columns) "SYMBOL" else if (length(columns)) columns[[1]] else NULL
  matched <- rep(FALSE, length(ids))
  for (col in columns) {
    vals <- unname(pull(col))
    rd[[.orgdb_target(col, feature_type)]] <- .fill(rd[[.orgdb_target(col, feature_type)]], vals)
    if (identical(col, primary)) matched <- !is.na(vals)
  }
  if (!is.null(matched_col)) rd[[matched_col]] <- matched  # mapped in OrgDb?
  SummarizedExperiment::rowData(dds) <- rd
  dds
}

#' Preview what OrgDb annotation makes available to join
#'
#' A small, non-committing table of the values [annotate_with_orgdb()] would
#' import for the first `n` features: the join key (`id`, the feature ids matched
#' against the OrgDb) as the first column, followed by one column per requested
#' OrgDb column (named by its `rowData` target). OrgDb is keyed by id, so the
#' "available" rows are exactly your dataset's first `n` ids resolved against it.
#' @inheritParams annotate_with_orgdb
#' @param n Number of features to preview (default 20).
#' @return A data.frame with up to `n` rows: `id` plus one column per imported value.
#' @export
orgdb_annotation_preview <- function(dds, organism = c("mouse", "human"),
                                     id_type = NULL, feature_type = "gene",
                                     columns = c("SYMBOL", "GENENAME"), n = 20L) {
  organism <- match.arg(organism)
  org <- .orgdb_org(organism, "orgdb_annotation_preview")
  columns <- intersect(columns, AnnotationDbi::columns(org))
  ids <- utils::head(rownames(dds), as.integer(n))
  id_type <- id_type %||% detect_id_type(ids)
  keytype <- .orgdb_keytype(id_type)
  keys <- if (keytype == "ENSEMBL") sub("\\..*$", "", ids) else ids
  pull <- function(col) suppressMessages(tryCatch(
    AnnotationDbi::mapIds(org, keys = keys, column = col, keytype = keytype, multiVals = "first"),
    error = function(e) stats::setNames(rep(NA_character_, length(keys)), keys)))
  out <- data.frame(id = ids, stringsAsFactors = FALSE, check.names = FALSE)
  for (col in columns) out[[.orgdb_target(col, feature_type)]] <- unname(pull(col))
  out
}

#' Annotation coverage over endogenous features
#'
#' @param dds A `DESeqDataSet`.
#' @param name_col The name column to check (e.g. `"gene_name"`).
#' @return list(`matched`, `total`) counting non-empty names among endogenous features.
#' @export
annotation_coverage <- function(dds, name_col) {
  rd <- SummarizedExperiment::rowData(dds)
  endo <- if ("feature_class" %in% colnames(rd)) rd$feature_class == "endogenous" else rep(TRUE, nrow(dds))
  if (!name_col %in% colnames(rd)) return(list(matched = 0L, total = sum(endo)))
  vals <- as.character(rd[[name_col]])[endo]
  list(matched = sum(!is.na(vals) & nzchar(vals)), total = sum(endo))
}

#' Which target columns already hold values (would be overwritten)
#'
#' Given the columns an annotation step would write, returns those that already
#' exist in `rowData` with at least one non-empty value -- i.e. matched rows
#' would overwrite existing data. Used to warn before applying.
#' @param dds A `DESeqDataSet`.
#' @param target_cols Character vector of column names the step will write.
#' @return Character vector of existing, non-empty target columns.
#' @export
annotation_overwrites <- function(dds, target_cols) {
  rd <- SummarizedExperiment::rowData(dds)
  existing <- intersect(as.character(target_cols), colnames(rd))
  existing[vapply(existing, function(c) {
    v <- rd[[c]]
    any(!is.na(v) & nzchar(as.character(v)))
  }, logical(1))]
}
