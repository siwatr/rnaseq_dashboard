# The shared app-state object and its helper API. This is the canonical store
# the whole dashboard reads from and writes to (see the `shiny-module` skill and
# CLAUDE.md "State model"). One mechanism powers invalidation, undo/reset, and
# the reproducibility action log.

#' Internal null-coalescing operator
#' @noRd
`%||%` <- function(x, y) if (is.null(x)) y else x

# Cap on the undo snapshot stack.
.undo_depth <- 5L

#' Create a fresh app-state object
#'
#' @return A [shiny::reactiveValues] with: `original` (immutable loaded object),
#'   `working` (current `dds`), `data_version` (bumped on every data edit),
#'   `history` (action log), `undo_stack`, and `meta` (app-level flags such as
#'   data type / feature type). `derived` is a reference-semantics environment
#'   (not a reactive field) used as the version-stamped cache for heavy results.
#' @export
new_app_state <- function() {
  state <- shiny::reactiveValues(
    original     = NULL,
    working      = NULL,
    data_version = 0L,
    history      = list(),
    undo_stack   = list(),
    n_edits      = 0L,      # net edits applied vs original (load/reset -> 0, undo -1)
    meta         = list()
  )
  # A plain environment so writing cache entries does not trigger reactivity
  # (avoids reactive-write-in-reactive churn); staleness is keyed on data_version.
  state$derived <- new.env(parent = emptyenv())
  state
}

.clear_derived <- function(state) {
  rm(list = ls(state$derived, all.names = TRUE), envir = state$derived)
}

.log <- function(state, entry) {
  state$history <- c(state$history, list(utils::modifyList(list(time = Sys.time()), entry)))
}

# Count features of a given feature_class (no class column -> all endogenous).
.class_count <- function(dds, class) {
  fc <- tryCatch(as.character(SummarizedExperiment::rowData(dds)$feature_class),
                 error = function(e) NULL)
  if (is.null(fc)) return(if (class == "endogenous") nrow(dds) else 0L)
  sum(fc == class)
}

#' Reactive read of the current working `dds`
#' @param state An app-state object from [new_app_state()].
#' @return A reactive expression yielding the working `dds` (or `NULL`).
#' @export
state_dds <- function(state) {
  shiny::reactive(state$working)
}

#' Load an object as the canonical dataset
#'
#' Sets `original` + `working` to `obj` (assumed already normalized by the
#' caller — see [as_input_dds()]/[ensure_feature_class()]/[ensure_logcounts()]),
#' clears derived/undo, resets the history to a single `load` entry, and bumps
#' `data_version`.
#'
#' @param state App-state object.
#' @param obj The `DESeqDataSet` to load.
#' @param source Short label of where it came from ("rds", "tabular", "demo").
#' @param meta Optional named list of app-level flags (data_type, feature_type, …).
#' @return `state`, invisibly.
#' @export
state_load <- function(state, obj, source = "rds", meta = list()) {
  state$original     <- obj
  state$working      <- obj
  state$meta         <- meta
  state$undo_stack   <- list()
  .clear_derived(state)
  state$history      <- list()
  state$n_edits      <- 0L
  state$data_version <- state$data_version + 1L
  .log(state, list(action = "load", source = source))
  invisible(state)
}

#' Apply an edit to the working object
#'
#' The single entry point for any change to samples/features/metadata. Pushes
#' the current `working` onto the undo stack, applies `fn`, **bumps
#' `data_version`** (invalidating all derived results), and logs `action`.
#'
#' @param state App-state object.
#' @param fn Function `(dds) -> dds` producing the new working object.
#' @param action Named list describing the edit (e.g. `list(action = "filter_features", min_count = 10)`).
#' @return `state`, invisibly.
#' @export
state_mutate <- function(state, fn, action = list()) {
  cur <- state$working
  if (is.null(cur)) stop("No dataset loaded; cannot mutate state.", call. = FALSE)
  state$undo_stack   <- utils::head(c(list(cur), state$undo_stack), .undo_depth)
  state$working      <- fn(cur)
  .clear_derived(state)
  state$n_edits      <- state$n_edits + 1L
  state$data_version <- state$data_version + 1L
  .log(state, action)
  invisible(state)
}

#' Get-or-compute a version-stamped derived artifact
#'
#' Returns the cached value for `key` iff it was computed under the current
#' `data_version` with identical `params`; otherwise computes `expr()`, caches
#' it, and returns it. Use for heavy results (VST, PCA, DESeq fit, DE tables).
#'
#' @param state App-state object.
#' @param key Cache key (character).
#' @param params List of parameters the result depends on.
#' @param expr A zero-arg function computing the value.
#' @return The cached or freshly computed value.
#' @export
state_derive <- function(state, key, params = list(), expr) {
  ver <- state$data_version
  if (exists(key, envir = state$derived, inherits = FALSE)) {
    hit <- get(key, envir = state$derived)
    if (identical(hit$version, ver) && identical(hit$params, params)) {
      return(hit$value)
    }
  }
  value <- expr()
  assign(key, list(value = value, version = ver, params = params), envir = state$derived)
  value
}

#' Restore the originally loaded object
#' @param state App-state object.
#' @return `state`, invisibly.
#' @export
state_reset <- function(state) {
  if (is.null(state$original)) return(invisible(state))
  state$working      <- state$original
  state$undo_stack   <- list()
  .clear_derived(state)
  state$n_edits      <- 0L
  state$data_version <- state$data_version + 1L
  .log(state, list(action = "reset"))
  invisible(state)
}

#' Undo the most recent edit
#' @param state App-state object.
#' @return `state`, invisibly.
#' @export
state_undo <- function(state) {
  if (length(state$undo_stack) == 0L) return(invisible(state))
  state$working      <- state$undo_stack[[1L]]
  state$undo_stack   <- state$undo_stack[-1L]
  .clear_derived(state)
  state$n_edits      <- max(0L, state$n_edits - 1L)
  state$data_version <- state$data_version + 1L
  .log(state, list(action = "undo"))
  invisible(state)
}

#' Summary metadata for the status bar
#'
#' @param state App-state object.
#' @return A list: `loaded`, and when loaded `data_type`, `feature_type`,
#'   `n_features`, `n_samples`, `assays`, `design`, `n_edits`, `n_undo` (undo
#'   steps currently available, capped at the snapshot depth), `n_endogenous` /
#'   `n_spike_in` / `n_exogenous` (feature-class counts), `data_version`, and
#'   `sce_per_cell` (TRUE when a single-cell object was coerced per-cell).
#' @export
state_meta <- function(state) {
  dds <- state$working
  if (is.null(dds)) return(list(loaded = FALSE))
  m <- state$meta %||% list()
  list(
    loaded       = TRUE,
    data_type    = m$data_type %||% "bulk",
    feature_type = m$feature_type %||% "feature",
    sce_per_cell = isTRUE(m$sce_per_cell),
    n_features   = nrow(dds),
    n_samples    = ncol(dds),
    assays       = SummarizedExperiment::assayNames(dds),
    design       = tryCatch(paste(deparse(DESeq2::design(dds)), collapse = " "),
                            error = function(e) NA_character_),
    n_edits      = state$n_edits %||% 0L,
    n_undo       = length(state$undo_stack),
    n_endogenous = .class_count(dds, "endogenous"),
    n_spike_in   = .class_count(dds, "spike_in"),
    n_exogenous  = .class_count(dds, "exogenous"),
    data_version = state$data_version
  )
}
