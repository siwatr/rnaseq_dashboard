# Pure helpers for named gene sets of interest (the Gene Sets page). A gene set
# is a structured record so it stays forward-compatible with the P7 "annotated"
# layer (id -> label) and so membership is stored NON-DESTRUCTIVELY: `ids` is the
# full authored set, and "present / absent in the current dataset" is derived
# live (gene_set_present / gene_set_absent) rather than trimmed. These functions
# are Shiny-free and unit-tested; the module (mod_geneset.R) wires them to state.

#' Construct a gene-set record
#'
#' @param ids Character vector of feature ids (deduplicated; NA / blank dropped).
#' @param kind `"simple"` (a plain id list) or `"annotated"` (P7: ids carry a
#'   per-gene label in `annotation`).
#' @param annotation Optional named character vector (id -> label) for an
#'   annotated set; `NULL` for a simple set.
#' @param source Provenance label (e.g. `"paste"`, `"DE: <contrast>"`,
#'   `"top-variable"`), recorded for the reproducibility export.
#' @return A list `list(ids, kind, annotation, source)`.
#' @export
new_gene_set <- function(ids, kind = c("simple", "annotated"),
                         annotation = NULL, source = "manual") {
  kind <- match.arg(kind)
  ids <- as.character(ids)
  ids <- unique(ids[!is.na(ids) & nzchar(ids)])
  list(ids = ids, kind = kind, annotation = annotation, source = source)
}

#' Commit ids into a named gene set (New or Append)
#'
#' `mode = "new"` creates a set under `name`, auto-suffixing (`name_2`, `name_3`,
#' ...) if the name is taken so an existing set is never silently overwritten.
#' `mode = "append"` unions the ids into the existing `name` (deduplicated), or
#' creates it when absent.
#'
#' @param sets The current named list of gene-set records (`state$gene_sets`).
#' @param name Target set name (trimmed; must be non-empty).
#' @param ids Feature ids to add.
#' @param mode `"new"` or `"append"`.
#' @param source Provenance label for a newly created set.
#' @return A list `list(sets, name)` — the updated store and the resolved name
#'   (which may differ from `name` when a New collision was auto-suffixed).
#' @export
gene_set_commit <- function(sets, name, ids, mode = c("new", "append"),
                            source = "manual") {
  mode <- match.arg(mode)
  if (is.null(sets)) sets <- list()
  name <- trimws(as.character(name)[1])
  if (!length(name) || is.na(name) || !nzchar(name))
    stop("A non-empty set name is required.", call. = FALSE)
  ids <- as.character(ids); ids <- unique(ids[!is.na(ids) & nzchar(ids)])

  if (identical(mode, "append") && !is.null(sets[[name]])) {
    cur <- sets[[name]]
    cur$ids <- unique(c(cur$ids, ids))
    sets[[name]] <- cur
    return(list(sets = sets, name = name))
  }
  # New (or Append to a non-existent name): pick a free name.
  final <- name
  if (!is.null(sets[[final]])) {
    i <- 2L
    while (!is.null(sets[[paste0(name, "_", i)]])) i <- i + 1L
    final <- paste0(name, "_", i)
  }
  sets[[final]] <- new_gene_set(ids, source = source)
  list(sets = sets, name = final)
}

#' Ids of a set present in / absent from the current dataset
#'
#' Non-destructive: the set keeps its full authored membership; these derive the
#' live view against the current feature ids (rownames). `present` preserves the
#' set's id order.
#' @param set A gene-set record (or a bare id vector).
#' @param feature_ids The current dataset's feature ids (`rownames(dds)`).
#' @return Character vector of ids.
#' @export
gene_set_present <- function(set, feature_ids) {
  ids <- if (is.list(set)) set$ids else as.character(set)
  ids[ids %in% feature_ids]
}

#' @rdname gene_set_present
#' @export
gene_set_absent <- function(set, feature_ids) {
  ids <- if (is.list(set)) set$ids else as.character(set)
  ids[!(ids %in% feature_ids)]
}

#' Split feature ids into groups by annotation column(s)
#'
#' Powers the Gene Sets table import: an imported table (e.g. a DESeq2 result
#' sheet from another analysis) is sliced into **one gene set per annotation
#' group** -- filter out `no_change`, split by direction, get `up` / `down` sets.
#' The group key is the annotation columns' values joined by `sep`.
#'
#' @param df A data.frame carrying the id column and the annotation column(s).
#'   Subset it to the rows you want *before* calling.
#' @param id_col Name of the column holding feature ids.
#' @param anno_cols Character vector of column names whose combined values key
#'   the groups. Empty (the default) yields a single group named `"All"`.
#' @param sep Separator joining multiple annotation values (default `"."`).
#' @return A named list of unique id vectors, one per group. `NA` / blank ids are
#'   dropped, as are groups left with no ids; `NA` annotation values key as
#'   `"NA"`. Returns an empty list when nothing survives.
#' @export
split_ids_by_group <- function(df, id_col, anno_cols = character(0), sep = ".") {
  stopifnot(is.data.frame(df))
  if (!id_col %in% names(df))
    stop("ID column '", id_col, "' is not in the table.", call. = FALSE)
  ids  <- as.character(df[[id_col]])
  keep <- !is.na(ids) & nzchar(ids)
  anno_cols <- intersect(anno_cols, names(df))
  if (!length(anno_cols)) {
    out <- unique(ids[keep])
    return(if (length(out)) list(All = out) else list())
  }
  key <- do.call(paste, c(lapply(anno_cols, function(cc) {
    v <- as.character(df[[cc]]); v[is.na(v)] <- "NA"; v
  }), list(sep = sep)))
  sp <- split(ids[keep], key[keep])
  sp <- lapply(sp, unique)
  sp[vapply(sp, length, integer(1)) > 0L]
}

# ===========================================================================
# File round-trip (P6d): JSON / GMT / TSV serializers.
#
# JSON is the FAITHFUL format -- it mirrors the structured record (ids + kind +
# source + annotation) and round-trips exactly. GMT and TSV are interchange
# formats that carry what their conventions support: GMT keeps name + a
# description (used for `source`) + genes; the long TSV keeps set / id (+ an
# annotation column, forward-compat with the P7 annotated layer). All are
# simple-set oriented (kind defaults to "simple" on read, matching P6 scope).
# The `to_*` functions return a single string (lines joined by "\n") so a
# download handler writes them uniformly with writeLines(); the `from_*`
# functions accept either a single string (possibly multi-line) or a character
# vector of lines.
# ===========================================================================

# Normalize file content to a character vector of lines.
.gs_lines <- function(x) {
  if (length(x) > 1L) return(as.character(x))
  x <- as.character(x)
  if (!length(x) || is.na(x)) return(character(0))
  strsplit(x, "\r?\n")[[1]]
}

#' Serialize gene sets to JSON
#'
#' Writes a versioned, faithful mirror of the gene-set store (`state$gene_sets`):
#' each set becomes `{ids, kind, source[, annotation]}` with `ids` always an
#' array (even for a single gene). Round-trips with [gene_sets_from_json()].
#'
#' @param sets A named list of gene-set records (`state$gene_sets`).
#' @param pretty Pretty-print the JSON (default `TRUE`).
#' @return A length-1 JSON string.
#' @export
gene_sets_to_json <- function(sets, pretty = TRUE) {
  if (is.null(sets)) sets <- list()
  sets <- lapply(sets, function(s) {
    # I() forces an array under auto_unbox so a 1-gene set stays ["x"], not "x".
    rec <- list(ids = I(as.character(s$ids %||% character(0))),
                kind = s$kind %||% "simple",
                source = s$source %||% "manual")
    if (length(s$annotation)) rec$annotation <- as.list(s$annotation)
    rec
  })
  if (!length(sets)) sets <- stats::setNames(list(), character(0))
  jsonlite::toJSON(list(ddsdashboard_gene_sets_version = 1L, gene_sets = sets),
                   auto_unbox = TRUE, pretty = pretty, null = "null")
}

#' Parse gene sets from JSON
#'
#' Inverse of [gene_sets_to_json()]. Accepts the wrapped object
#' (`{..., gene_sets}`) or a bare `{name: {...}}` / `{name: [ids]}` object (a
#' minimal hand-authored form is tolerated). Each entry is normalized through
#' [new_gene_set()]. Invalid JSON `stop()`s.
#'
#' @param txt A JSON string (or a character vector of lines).
#' @return A named list of gene-set records (possibly empty).
#' @export
gene_sets_from_json <- function(txt) {
  txt <- paste(.gs_lines(txt), collapse = "\n")
  obj <- jsonlite::fromJSON(txt, simplifyVector = TRUE,
                            simplifyDataFrame = FALSE, simplifyMatrix = FALSE)
  gs <- if (is.list(obj) && !is.null(obj$gene_sets)) obj$gene_sets else obj
  if (!is.list(gs) || !length(gs) || is.null(names(gs)))
    return(stats::setNames(list(), character(0)))
  out <- list()
  for (nm in names(gs)) {
    if (!nzchar(nm)) next
    it <- gs[[nm]]
    if (is.list(it)) {
      ids  <- as.character(it$ids)
      kind <- it$kind %||% "simple"
      if (!kind %in% c("simple", "annotated")) kind <- "simple"
      anno <- if (length(it$annotation)) unlist(it$annotation) else NULL
      src  <- it$source %||% "import: json"
    } else {                                   # bare {name: [ids]}
      ids <- as.character(it); kind <- "simple"; anno <- NULL; src <- "import: json"
    }
    out[[nm]] <- new_gene_set(ids, kind = kind, annotation = anno, source = src)
  }
  out
}

#' Serialize gene sets to GMT
#'
#' MSigDB GMT: one tab-delimited line per set --
#' `name <TAB> description <TAB> gene1 <TAB> gene2 ...`. The set `source` is
#' written to the description field. Round-trips ids with [gene_sets_from_gmt()]
#' (`kind`/`annotation` are not represented and default to simple on read).
#'
#' @param sets A named list of gene-set records.
#' @return A length-1 string (lines joined by `"\n"`).
#' @export
gene_sets_to_gmt <- function(sets) {
  if (is.null(sets) || !length(sets)) return("")
  lines <- vapply(names(sets), function(nm) {
    s <- sets[[nm]]
    paste(c(nm, s$source %||% "", as.character(s$ids %||% character(0))), collapse = "\t")
  }, character(1))
  paste(lines, collapse = "\n")
}

#' Parse gene sets from GMT
#'
#' Inverse of [gene_sets_to_gmt()]. Each non-empty line is
#' `name <TAB> description <TAB> genes...`; the description populates `source`.
#'
#' @param txt A GMT string (or a character vector of lines).
#' @return A named list of gene-set records (possibly empty).
#' @export
gene_sets_from_gmt <- function(txt) {
  lines <- .gs_lines(txt); lines <- lines[nzchar(trimws(lines))]
  out <- list()
  for (ln in lines) {
    parts <- strsplit(ln, "\t", fixed = TRUE)[[1]]
    nm <- trimws(parts[1] %||% "")
    if (!nzchar(nm)) next
    desc <- if (length(parts) >= 2L) trimws(parts[2]) else ""
    ids  <- if (length(parts) >= 3L) parts[-(1:2)] else character(0)
    out[[nm]] <- new_gene_set(ids, source = if (nzchar(desc)) desc else "import: gmt")
  }
  out
}

#' Serialize gene sets to a long TSV
#'
#' A tidy long table with columns `set`, `id`, `annotation` -- one row per
#' (set, gene). `annotation` is blank for simple sets (forward-compat with the
#' P7 annotated layer). Round-trips ids with [gene_sets_from_tsv()].
#'
#' @param sets A named list of gene-set records.
#' @return A length-1 string (a TSV, header included).
#' @export
gene_sets_to_tsv <- function(sets) {
  if (is.null(sets)) sets <- list()
  rows <- lapply(names(sets), function(nm) {
    ids <- as.character(sets[[nm]]$ids %||% character(0))
    if (!length(ids)) return(NULL)
    anno <- sets[[nm]]$annotation
    a <- if (length(anno)) unname(anno[ids]) else rep(NA_character_, length(ids))
    data.frame(set = nm, id = ids, annotation = a,
               check.names = FALSE, stringsAsFactors = FALSE)
  })
  df <- do.call(rbind, rows)
  if (is.null(df))
    df <- data.frame(set = character(0), id = character(0), annotation = character(0))
  sub("\n$", "", readr::format_tsv(df, na = ""))    # drop only the trailing newline
}

#' Parse gene sets from a long TSV
#'
#' Inverse of [gene_sets_to_tsv()]. Needs `set` and `id` columns (matched
#' case-insensitively); groups ids by set, preserving first-appearance order.
#' The `annotation` column is ignored (simple sets, P6 scope).
#'
#' @param txt A TSV string (or a character vector of lines).
#' @return A named list of gene-set records (possibly empty).
#' @export
gene_sets_from_tsv <- function(txt) {
  txt <- paste(.gs_lines(txt), collapse = "\n")
  if (!nzchar(trimws(txt))) return(stats::setNames(list(), character(0)))
  df <- as.data.frame(readr::read_tsv(
    I(txt), show_col_types = FALSE,
    col_types = readr::cols(.default = readr::col_character())))
  if (!nrow(df)) return(stats::setNames(list(), character(0)))
  nms <- tolower(names(df))
  sc <- match("set", nms); ic <- match("id", nms)
  if (is.na(sc) || is.na(ic))
    stop("A long gene-set TSV needs 'set' and 'id' columns.", call. = FALSE)
  set_v <- as.character(df[[sc]]); id_v <- as.character(df[[ic]])
  by <- split(id_v, set_v)
  out <- list()
  for (nm in unique(set_v)) {
    if (is.na(nm) || !nzchar(nm)) next
    out[[nm]] <- new_gene_set(by[[nm]], source = "import: tsv")
  }
  out
}

#' Read gene sets from a file (format auto-detected)
#'
#' Dispatches to [gene_sets_from_json()] / [gene_sets_from_gmt()] /
#' [gene_sets_from_tsv()] by file extension (`.json` / `.gmt` / `.tsv` / `.txt`
#' / `.csv`), falling back to a content sniff (leading `{`/`[` -> JSON; a
#' `set`+`id` header -> TSV; otherwise GMT).
#'
#' @param path Path to the file to read.
#' @param format One of `"auto"` (default), `"json"`, `"gmt"`, `"tsv"`.
#' @param name Original file name used for extension detection (defaults to
#'   `path`; supply the upload's display name when `path` is a temp file).
#' @return A named list of gene-set records (possibly empty).
#' @export
gene_sets_from_file <- function(path, format = c("auto", "json", "gmt", "tsv"),
                                name = path) {
  format <- match.arg(format)
  if (identical(format, "auto")) {
    format <- switch(tolower(tools::file_ext(name)),
      json = "json", gmt = "gmt", tsv = "tsv", txt = "tsv", csv = "tsv", "auto")
  }
  txt <- readLines(path, warn = FALSE)
  if (identical(format, "auto")) {
    first <- trimws(paste(txt, collapse = "\n"))
    format <- if (startsWith(first, "{") || startsWith(first, "[")) "json"
      else if (all(c("set", "id") %in% tolower(strsplit(txt[1] %||% "", "\t")[[1]]))) "tsv"
      else "gmt"
  }
  switch(format,
    json = gene_sets_from_json(txt),
    gmt  = gene_sets_from_gmt(txt),
    tsv  = gene_sets_from_tsv(txt))
}
