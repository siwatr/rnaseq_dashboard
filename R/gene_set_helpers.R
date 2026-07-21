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

# A set's authored ids (record or bare vector).
.gs_ids <- function(set) if (is.list(set)) as.character(set$ids) else as.character(set)

#' Ids a gene set contributes to a comparison
#'
#' The Compare tab's "Within this dataset" toggle chooses whether a set
#' contributes its full authored membership or only the ids present in the
#' current dataset. `within = TRUE` -> [gene_set_present()]; else the full ids.
#' @param within When TRUE, keep only ids in `feature_ids`; else the full
#'   authored membership.
#' @rdname gene_set_present
#' @export
gene_set_ids_for <- function(set, feature_ids = NULL, within = FALSE) {
  unique(if (isTRUE(within)) gene_set_present(set, feature_ids) else .gs_ids(set))
}

#' Per-set size breakdown for the Compare Stats bar
#'
#' Powers the Gene Sets Compare tab's set-size bar chart. Non-destructive: sizes
#' derive live from each set's authored membership against the current dataset.
#' @param sets A named list of gene-set records (or bare id vectors).
#' @param feature_ids Current dataset feature ids (`rownames(dds)`).
#' @param within When TRUE, one "present" row per set (present ids only); else a
#'   present + absent breakdown of the full authored membership (stacked bar).
#' @param order Bar order (sets the `set` factor levels): `"none"` (input order),
#'   `"inc"`/`"dec"` (by displayed total size), `"az"`/`"za"` (by name).
#' @return A long data.frame: `set` (factor, ordered per `order`), `status`
#'   (`"present"`/`"absent"`, factor), `n` (integer). Empty when `sets` is empty.
#' @export
gene_set_size_frame <- function(sets, feature_ids = NULL, within = FALSE,
                                order = c("none", "inc", "dec", "az", "za")) {
  order <- match.arg(order)
  if (!length(sets)) {
    return(data.frame(set = factor(character(0)),
                      status = factor(character(0), levels = c("present", "absent")),
                      n = integer(0), stringsAsFactors = FALSE))
  }
  nms <- names(sets)
  present <- vapply(sets, function(s) length(gene_set_present(s, feature_ids)), integer(1))
  rows <- if (isTRUE(within)) {
    data.frame(set = nms, status = "present", n = present, stringsAsFactors = FALSE)
  } else {
    absent <- vapply(sets, function(s) length(gene_set_absent(s, feature_ids)), integer(1))
    rbind(data.frame(set = nms, status = "present", n = present, stringsAsFactors = FALSE),
          data.frame(set = nms, status = "absent",  n = absent,  stringsAsFactors = FALSE))
  }
  # Order the set factor: by displayed total (inc/dec), name (az/za), or input.
  lvls <- switch(order,
    az   = sort(nms),
    za   = rev(sort(nms)),
    none = nms,
    { tot <- vapply(nms, function(s) sum(rows$n[rows$set == s]), integer(1))
      nms[order(tot, decreasing = identical(order, "dec"))] })
  rows$set    <- factor(rows$set, levels = lvls)
  rows$status <- factor(rows$status, levels = c("present", "absent"))
  rows
}

#' Named id-lists for a set-overlap comparison (Euler / Venn / UpSet)
#'
#' @inheritParams gene_set_size_frame
#' @return A named list of unique id vectors, one per set (a set that
#'   contributes nothing is kept as `character(0)` so callers can report it).
#' @export
gene_set_overlap_list <- function(sets, feature_ids = NULL, within = FALSE) {
  if (!length(sets)) return(list())
  lapply(sets, gene_set_ids_for, feature_ids = feature_ids, within = within)
}

# Present / absent stacked-bar colours (present solid, absent a faded tint).
.gene_set_presence_palette <- c(present = "#8b58db", absent = "#d9cef0")

#' Present / absent colours for the Gene Sets Compare stats bar
#'
#' Mirrors [removal_status_colors()]: the default 2-level scheme, overridable by
#' a Palette `other/gene_set_presence` config (edited on the Palette page).
#' @param config A palette config `list(name, colors, custom)`, or `NULL` for the
#'   default scheme.
#' @return A named character vector `c(present = , absent = )`.
#' @export
gene_set_presence_colors <- function(config = NULL) {
  if (is.null(config)) return(.gene_set_presence_palette)
  palette_discrete(names(.gene_set_presence_palette), config$colors,
                   config$name %||% "Custom palette", config$custom)
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

#' Combine several gene sets into one id -> label annotation
#'
#' The P7 "annotated" layer: a second-level object grouping several member gene
#' sets, where each gene's label = the member set it belongs to. Non-destructive
#' -- it reads each set's authored ids and produces one label per gene over the
#' union. Feeds `new_gene_set(kind = "annotated", annotation = <result$annotation>)`.
#'
#' A gene in more than one member set is resolved by `shared`:
#' - `"concat"` (default): the member-set names joined by `sep` (in input order),
#'   e.g. `"up;hypoxia"` -- stays one label per gene, so `row_split` / saving keep
#'   working, and overlaps read as their own level.
#' - `"label"`: a single distinct level named `label` (default `"multiple"`).
#'   **Rejected when `label` collides with a member-set name** (it would silently
#'   merge a real set into the overlap bucket).
#' - `"first"`: the first member set (in input order) the gene appears in.
#'
#' @param sets A **named** list of gene-set records (or bare id vectors); each
#'   name is the label for genes in that set. `NULL` / empty yields empty output.
#' @param shared How to label a gene in >1 set (see above).
#' @param sep Separator for `shared = "concat"` (default `";"`).
#' @param label Level name for `shared = "label"` (default `"multiple"`).
#' @return A list `list(annotation, levels, shared_ids, note)` -- `annotation` a
#'   named character vector (id -> label), `levels` the labels in display order
#'   (member-set names in input order, then any overlap labels in first-appearance
#'   order), `shared_ids` the genes in >1 set, `note` a human summary string.
#' @export
combine_gene_set_annotation <- function(sets, shared = c("concat", "label", "first"),
                                        sep = ";", label = "multiple") {
  shared <- match.arg(shared)
  empty <- list(annotation = stats::setNames(character(0), character(0)),
                levels = character(0), shared_ids = character(0), note = "")
  if (is.null(sets) || !length(sets)) return(empty)
  nms <- names(sets)
  if (is.null(nms)) stop("`sets` must be a named list.", call. = FALSE)
  keep <- nzchar(nms) & !is.na(nms)
  sets <- sets[keep]; nms <- nms[keep]
  if (!length(sets)) return(empty)
  if (identical(shared, "label")) {
    label <- trimws(as.character(label)[1])
    if (!nzchar(label)) stop("`label` must be non-empty.", call. = FALSE)
    if (label %in% nms)
      stop("`label` ('", label, "') collides with a member-set name; pick another.",
           call. = FALSE)
  }
  # (id, set) pairs in input set order; drop NA / blank ids.
  pairs <- do.call(rbind, lapply(nms, function(nm) {
    ids <- .gs_ids(sets[[nm]]); ids <- ids[!is.na(ids) & nzchar(ids)]
    if (!length(ids)) return(NULL)
    data.frame(id = ids, set = nm, stringsAsFactors = FALSE)
  }))
  if (is.null(pairs)) return(empty)
  # Membership per id (set names in input order), keyed by first-appearance id order.
  by <- split(pairs$set, factor(pairs$id, levels = unique(pairs$id)))
  label_for <- function(s) {
    if (length(s) == 1L) return(s)
    switch(shared, concat = paste(s, collapse = sep), label = label, first = s[1])
  }
  labels <- vapply(by, label_for, character(1))
  annotation <- stats::setNames(unname(labels), names(by))
  shared_ids <- names(by)[lengths(by) > 1L]
  singles <- intersect(nms, labels)             # member-set names used as labels
  levels  <- c(singles, setdiff(unique(labels), singles))
  note <- if (length(shared_ids))
    sprintf("%d gene(s) belong to more than one set.", length(shared_ids)) else ""
  # concat can (rarely) produce an overlap label equal to a member-set name (a set
  # literally named like a concatenation, e.g. "a;b"): the overlap genes then merge
  # into that set's level. Same hazard the "label" mode guards against -- surface it.
  if (identical(shared, "concat") && length(shared_ids)) {
    collide <- intersect(nms, unique(annotation[shared_ids]))
    if (length(collide))
      note <- paste0(note, sprintf(" Overlap label(s) %s match a member-set name.",
                                   paste(collide, collapse = ", ")))
  }
  list(annotation = annotation, levels = levels, shared_ids = shared_ids, note = note)
}

#' Per-level gene counts for an annotated set (the composition view)
#'
#' Powers the Annotation tab's composition bar. Non-destructive: counts derive
#' live from the record's `annotation` map (id -> label).
#' @param set An annotated gene-set record (with an `annotation` map).
#' @param feature_ids Current dataset feature ids (`rownames(dds)`); used when
#'   `within = TRUE`.
#' @param within When TRUE, count only genes present in the dataset.
#' @return A data.frame `level` (character, first-appearance order) / `n`
#'   (integer), or an empty frame when there is no annotation.
#' @export
gene_set_annotation_composition <- function(set, feature_ids = NULL, within = FALSE) {
  anno <- if (is.list(set)) set$annotation else NULL
  empty <- data.frame(level = character(0), n = integer(0), stringsAsFactors = FALSE)
  if (is.null(anno) || !length(anno)) return(empty)
  ids <- names(anno)
  if (isTRUE(within) && !is.null(feature_ids)) {
    keep <- ids %in% feature_ids; anno <- anno[keep]
  }
  if (!length(anno)) return(empty)
  lv <- unique(unname(anno))                        # first-appearance level order
  n  <- vapply(lv, function(l) sum(anno == l), integer(1))
  data.frame(level = lv, n = unname(n), stringsAsFactors = FALSE)
}

#' Colours for a gene-set annotation's levels
#'
#' Mirrors [gene_set_presence_colors()] / [removal_status_colors()]: a default
#' qualitative scheme, overridable by a Palette `Gene Set` config (edited on the
#' Palette page, keyed by the annotation name).
#' @param levels The annotation's levels (character).
#' @param config A palette config `list(name, colors, custom)`, or `NULL` for the
#'   default qualitative scheme.
#' @return A named character vector (level -> hex).
#' @export
gene_set_anno_colors <- function(levels, config = NULL) {
  levels <- as.character(levels)
  if (is.null(config))
    return(palette_discrete(levels, NULL, "Okabe-Ito"))
  palette_discrete(levels, config$colors, config$name %||% "Custom palette", config$custom)
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
  # Detect the wrapped form by the version KEY, not the `gene_sets` name -- else a
  # bare {name: [ids]} object holding a set literally named "gene_sets" would be
  # misread as the wrapper. The app always writes the version key.
  gs <- if (is.list(obj) && !is.null(obj[["ddsdashboard_gene_sets_version"]]))
    obj$gene_sets else obj
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
    I(txt), show_col_types = FALSE, na = "",           # symmetric with the writer;
    col_types = readr::cols(.default = readr::col_character())))   # keeps a literal "NA" id
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
    # `.txt` is genuinely ambiguous (a GMT or a long TSV both get saved as .txt),
    # so leave it on "auto" for the content sniff below rather than forcing tsv.
    format <- switch(tolower(tools::file_ext(name)),
      json = "json", gmt = "gmt", tsv = "tsv", csv = "tsv", "auto")
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
