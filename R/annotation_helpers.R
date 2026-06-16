# Pure helpers for OrgDb annotation: detect the feature id type, map ids to
# gene symbol / chromosome / description via a Bioconductor OrgDb, and report
# coverage. Org packages are Suggests (loaded on demand). GTF annotation and
# ERCC concentrations come in later PRs.

.orgdb_pkg <- c(mouse = "org.Mm.eg.db", human = "org.Hs.eg.db")

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

#' Annotate features from a Bioconductor OrgDb
#'
#' Writes `rowData`: `<feature_type>_name` (SYMBOL) and `description` (GENENAME).
#' Only fills where a mapping is found (keeps existing values otherwise);
#' `feature_class` is left untouched. Chromosome is intentionally not taken from
#' OrgDb (its `CHR` accessor is deprecated) — it comes from the GTF instead.
#'
#' @param dds A `DESeqDataSet`.
#' @param organism `"mouse"` (org.Mm.eg.db) or `"human"` (org.Hs.eg.db).
#' @param id_type `"ensembl"`/`"entrez"`/`"symbol"`; auto-detected when `NULL`.
#' @param feature_type Feature unit; sets the name column `<feature_type>_name`.
#' @param matched_col Optional name of a logical `rowData` column to write,
#'   flagging which features were mapped in the OrgDb (`NULL` = none).
#' @return The annotated `DESeqDataSet`.
#' @export
annotate_with_orgdb <- function(dds, organism = c("mouse", "human"),
                                id_type = NULL, feature_type = "gene",
                                matched_col = NULL) {
  organism <- match.arg(organism)
  pkg <- .orgdb_pkg[[organism]]
  for (p in c("AnnotationDbi", pkg)) {
    if (!requireNamespace(p, quietly = TRUE)) {
      stop("annotate_with_orgdb() needs the '", p, "' package.", call. = FALSE)
    }
  }
  org <- getExportedValue(pkg, pkg)
  ids <- rownames(dds)
  id_type <- id_type %||% detect_id_type(ids)
  keytype <- switch(id_type,
                    ensembl = "ENSEMBL", entrez = "ENTREZID", symbol = "SYMBOL",
                    stop("Unknown id_type: ", id_type, call. = FALSE))
  keys <- if (keytype == "ENSEMBL") sub("\\..*$", "", ids) else ids  # drop version suffix
  pull <- function(col) {
    suppressMessages(tryCatch(
      AnnotationDbi::mapIds(org, keys = keys, column = col, keytype = keytype,
                            multiVals = "first"),
      error = function(e) stats::setNames(rep(NA_character_, length(keys)), keys)
    ))
  }
  symbol <- unname(pull("SYMBOL"))
  rd <- SummarizedExperiment::rowData(dds)
  name_col <- paste0(feature_type, "_name")
  rd[[name_col]]      <- .fill(rd[[name_col]],      symbol)
  rd[["description"]] <- .fill(rd[["description"]], unname(pull("GENENAME")))
  if (!is.null(matched_col)) rd[[matched_col]] <- !is.na(symbol)  # mapped in OrgDb?
  SummarizedExperiment::rowData(dds) <- rd
  dds
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
