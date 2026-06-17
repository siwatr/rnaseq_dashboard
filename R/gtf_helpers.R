# Pure helpers for GTF annotation: read a GTF, compute feature length over a
# chosen feature `type` (union of ranges per group, so it matches how reads were
# counted -- exonic for mature mRNA, gene/transcript body for nascent reads),
# and overlay GTF attributes onto rowData (authoritative over OrgDb). Heavy
# Bioconductor packages are Suggests, loaded on demand (like the OrgDb path).
# Matching is by dds rownames: extra GTF entries are ignored and dds features
# absent from the GTF keep their existing values (never wiped).

.need <- function(pkgs) {
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      stop("GTF annotation needs the '", p, "' package.", call. = FALSE)
    }
  }
}

#' Import a GTF/GFF file
#'
#' Thin wrapper over `rtracklayer::import()` (accepts `.gtf` / `.gtf.gz`).
#' @param path Path to the annotation file.
#' @return A `GRanges`.
#' @export
import_gtf <- function(path) {
  .need("rtracklayer")
  rtracklayer::import(path)
}

#' Feature `type` values present in a GTF
#' @param gtf A `GRanges` from [import_gtf()].
#' @return Sorted unique `type` values (e.g. `gene`/`transcript`/`exon`/`CDS`).
#' @export
gtf_feature_types <- function(gtf) {
  sort(unique(as.character(gtf$type)))
}

#' Columns a GTF exposes for matching / import
#'
#' `seqnames` plus the attribute columns (`gene_id`, `gene_name`,
#' `gene_biotype`, `transcript_id`, ...); the standard range/bookkeeping columns
#' are dropped.
#' @param gtf A `GRanges` from [import_gtf()].
#' @return Character vector of column names.
#' @export
available_gtf_columns <- function(gtf) {
  nm <- colnames(as.data.frame(gtf[0L]))
  setdiff(nm, c("start", "end", "width", "strand", "source", "score", "phase", "type"))
}

#' Union feature length over a chosen GTF `type`
#'
#' Keeps rows of the given `type`, groups by `group_col`, reduces overlapping /
#' adjacent ranges, and sums their widths -- the effective length for TPM/FPKM.
#' Pick `type` to match how reads were quantified: `exon` for mature mRNA, or
#' `gene`/`transcript` when reads fall outside exons (e.g. nascent transcripts).
#' @param gtf A `GRanges` from [import_gtf()].
#' @param type Feature type to measure (default `"exon"`).
#' @param group_col Attribute to group by (default `"gene_id"`).
#' @return Named numeric vector of lengths (names = `group_col` values).
#' @export
gtf_feature_lengths <- function(gtf, type = "exon", group_col = "gene_id") {
  .need(c("GenomicRanges", "IRanges", "S4Vectors"))
  types <- as.character(gtf$type)
  sub <- gtf[!is.na(types) & types == type]
  if (!length(sub)) stop("No '", type, "' rows in the GTF.", call. = FALSE)
  md <- S4Vectors::mcols(sub)
  if (!group_col %in% colnames(md)) {
    stop("Column '", group_col, "' not found in the GTF.", call. = FALSE)
  }
  grp <- as.character(md[[group_col]])
  keep <- !is.na(grp) & nzchar(grp)
  sub <- sub[keep]; grp <- grp[keep]
  by <- S4Vectors::split(sub, grp)
  len <- sum(IRanges::width(GenomicRanges::reduce(by)))   # named (by group) numeric
  stats::setNames(as.numeric(len), names(len))
}

#' First-rows preview of a GTF as a data.frame
#'
#' A small head() of the parsed `GRanges` for display (choosing import columns),
#' converting only the first `n` rows -- never the whole object.
#' @param gtf A `GRanges` from [import_gtf()].
#' @param n Maximum rows to show (default 20).
#' @return A data.frame with at most `n` rows.
#' @export
gtf_preview <- function(gtf, n = 20L) {
  if (is.null(gtf) || !length(gtf)) return(data.frame())
  as.data.frame(gtf[seq_len(min(as.integer(n), length(gtf)))], stringsAsFactors = FALSE)
}

#' One-row-per-feature attribute table from a GTF
#'
#' The first row seen for each `group_col` value (the gene/transcript row in a
#' well-formed GTF), carrying `seqnames` and every attribute column.
#' @param gtf A `GRanges` from [import_gtf()].
#' @param group_col Attribute keyed on (default `"gene_id"`).
#' @return A data.frame, row names = `group_col` values.
#' @export
gtf_attribute_table <- function(gtf, group_col = "gene_id") {
  .need(c("S4Vectors", "GenomicRanges"))
  md <- S4Vectors::mcols(gtf)
  # Grouping key vector, without materializing the whole GRanges as a data.frame
  # (that would copy every row + expand the Rle seqnames -- the memory spike).
  keys <- if (group_col == "seqnames") as.character(GenomicRanges::seqnames(gtf))
          else if (group_col %in% colnames(md)) as.character(md[[group_col]])
          else stop("Column '", group_col, "' not found in the GTF.", call. = FALSE)
  keep <- !is.na(keys) & nzchar(keys) & !duplicated(keys)
  # Convert only the deduplicated subset (n_groups rows, not n_features).
  out <- as.data.frame(md[keep, , drop = FALSE], stringsAsFactors = FALSE)
  out$seqnames <- as.character(GenomicRanges::seqnames(gtf))[keep]
  rownames(out) <- keys[keep]
  out
}

.strip_version <- function(x) sub("\\..*$", "", as.character(x))

# Pick the GTF column to match dds rownames against.
.resolve_match_col <- function(match_col, ids, gtf) {
  cols <- available_gtf_columns(gtf)
  if (!identical(match_col, "auto")) {
    if (!match_col %in% cols) stop("match_col '", match_col, "' not in the GTF.", call. = FALSE)
    return(match_col)
  }
  if (detect_id_type(ids) == "ensembl" && "gene_id" %in% cols) return("gene_id")
  if ("gene_name" %in% cols) return("gene_name")
  if ("gene_id"   %in% cols) return("gene_id")
  stop("Could not auto-resolve a match column; please choose one.", call. = FALSE)
}

#' Annotate features from a GTF (authoritative over OrgDb)
#'
#' Matches dds rownames to a GTF column, overlays the selected attribute columns
#' onto `rowData` (filling where the GTF has a value, keeping existing values
#' elsewhere), and -- when `compute_length` -- writes union-`length_type`
#' lengths to `rowData$feature_length`. `gene_name` is stored as
#' `<feature_type>_name` and `seqnames` as `chromosome`; other columns keep their
#' names. `feature_class` and counts are untouched.
#'
#' @param dds A `DESeqDataSet`.
#' @param gtf A `GRanges` from [import_gtf()].
#' @param match_col GTF column matched to rownames; `"auto"` (id then name).
#' @param import_cols Attribute columns to copy into `rowData` (character).
#' @param compute_length Whether to compute and store `feature_length`.
#' @param length_type Feature `type` used for length (default `"exon"`).
#' @param feature_type Feature unit; sets the name column `<feature_type>_name`.
#' @param matched_col Optional name of a logical `rowData` column to write,
#'   flagging which features matched the GTF (`NULL` = none).
#' @return list(`dds`, `report` = list(`matched`, `total`, `length_set`, `length_complete`)).
#' @export
annotate_with_gtf <- function(dds, gtf, match_col = "auto", import_cols = NULL,
                              compute_length = FALSE, length_type = "exon",
                              feature_type = "gene", matched_col = NULL) {
  ids <- rownames(dds)
  match_col <- .resolve_match_col(match_col, ids, gtf)
  tab <- gtf_attribute_table(gtf, group_col = match_col)

  is_ens   <- detect_id_type(ids) == "ensembl"
  dds_keys <- if (is_ens) .strip_version(ids) else as.character(ids)
  tab_keys <- if (is_ens) .strip_version(rownames(tab)) else rownames(tab)
  idx      <- match(dds_keys, tab_keys)
  matched  <- !is.na(idx)

  rd <- SummarizedExperiment::rowData(dds)
  name_col <- paste0(feature_type, "_name")
  for (col in import_cols) {
    if (!col %in% colnames(tab)) next
    vals <- rep(NA_character_, length(ids))
    vals[matched] <- as.character(tab[[col]][idx[matched]])
    target <- switch(col, gene_name = name_col, seqnames = "chromosome", col)
    rd[[target]] <- .fill(rd[[target]], vals)         # GTF wins where matched
  }

  length_set <- 0L
  if (compute_length) {
    lens <- gtf_feature_lengths(gtf, type = length_type, group_col = match_col)
    len_keys <- if (is_ens) .strip_version(names(lens)) else names(lens)
    lidx <- match(dds_keys, len_keys)
    new  <- rep(NA_real_, length(ids))
    ok   <- !is.na(lidx)
    new[ok] <- as.numeric(lens[lidx[ok]])
    existing <- if (is.null(rd[["feature_length"]])) rep(NA_real_, length(ids))
                else as.numeric(rd[["feature_length"]])
    existing[ok] <- new[ok]                            # GTF authoritative
    rd[["feature_length"]] <- existing
    length_set <- sum(!is.na(existing))
  }

  if (!is.null(matched_col)) rd[[matched_col]] <- matched   # matched in GTF?
  SummarizedExperiment::rowData(dds) <- rd
  list(dds = dds, report = list(
    matched         = sum(matched),
    total           = length(ids),
    length_set      = length_set,
    length_complete = has_feature_length(dds)
  ))
}

#' Count how many dds features match a GTF
#'
#' Non-committing tally for a coverage banner: how many dataset rownames resolve
#' against the GTF on the chosen match column (Ensembl ids matched
#' version-insensitively, as in [annotate_with_gtf()]).
#' @param dds A `DESeqDataSet`.
#' @param gtf A `GRanges` from [import_gtf()].
#' @param match_col GTF column matched to rownames; `"auto"` (id then name).
#' @return list(`matched`, `total`).
#' @export
gtf_match_count <- function(dds, gtf, match_col = "auto") {
  ids <- rownames(dds)
  match_col <- .resolve_match_col(match_col, ids, gtf)
  tab <- gtf_attribute_table(gtf, group_col = match_col)
  is_ens   <- detect_id_type(ids) == "ensembl"
  dds_keys <- if (is_ens) .strip_version(ids) else as.character(ids)
  tab_keys <- if (is_ens) .strip_version(rownames(tab)) else rownames(tab)
  list(matched = sum(!is.na(match(dds_keys, tab_keys))), total = length(ids))
}

#' Adopt an existing numeric rowData column as feature_length
#'
#' For inputs that already carry length (some quantifiers / imported `dds`).
#' @param dds A `DESeqDataSet`.
#' @param col A numeric `rowData` column name.
#' @return The updated `DESeqDataSet`.
#' @export
set_feature_length_from_column <- function(dds, col) {
  rd <- SummarizedExperiment::rowData(dds)
  if (!col %in% colnames(rd)) stop("Unknown rowData column: ", col, call. = FALSE)
  v <- suppressWarnings(as.numeric(as.character(rd[[col]])))
  if (all(is.na(v))) stop("Column '", col, "' is not numeric.", call. = FALSE)
  if (any(!is.na(v) & v <= 0)) stop("feature_length values must be positive.", call. = FALSE)
  rd[["feature_length"]] <- v
  SummarizedExperiment::rowData(dds) <- rd
  dds
}
