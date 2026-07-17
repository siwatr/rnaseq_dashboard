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
