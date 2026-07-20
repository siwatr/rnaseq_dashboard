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
    meta         = list(),
    # Global, session-scoped UI preference (not data-dependent, so NOT in meta):
    # the plot engine toggle (FALSE = static ggplot, TRUE = interactive plotly).
    # Read by plot modules; untouched by load/reset.
    plot_interactive = FALSE,
    # Project colour palette (the Palette page). A UI preference like
    # plot_interactive: per-attribute colour mappings, untouched by load/reset,
    # no data_version impact. Shape: list(colData = list(<col> = list(type,
    # palette, pins)), ...). Empty list = fall back to thematic / default colours.
    palette = list(),
    # Sample-removal selections owned by the QC Filtering page but promoted here
    # so other plot pages (PCA, ...) can colour/shape by them: `samp_pool` is the
    # staged removal pool (sample ids); `samp_flags` is the latest flag_samples()
    # data.frame (or NULL). Session UI state (like palette) - the QC page clears
    # the pool on data load; no data_version impact.
    samp_pool  = character(0),
    samp_flags = NULL,
    # Design-scoped version stamp: bumped by state_set_design() when the model
    # formula / reference levels change, WITHOUT bumping data_version -- so the DE
    # fit (which keys on it) recomputes while design-independent caches (PCA/VST/QC,
    # keyed on data_version) survive. See state_set_design().
    design_version = 0L,
    # Differential-expression results owned by the DE page but held here so the
    # status bar, the Gene Sets page (add-from-DEG), and the Expression page can
    # read them. `contrasts` = list of stored specs; `results` = named list of
    # result data.frames; `active` = the shown contrast; `stamp` = the
    # (data_version, design_version) the fit ran under (staleness); `shrink` = the
    # shrinkage method. Cleared on data load; marked stale (not cleared) on edits.
    de = list(contrasts = list(), results = list(), active = NULL,
              stamp = NULL, shrink = "apeglm"),
    # Named gene sets of interest (the Gene Sets page). A session/UI field like
    # `palette`: NO data_version impact. Each entry is a structured record
    # list(ids, kind = "simple"|"annotated", annotation, source); `ids` is the
    # FULL authored membership (never trimmed by feature edits) -- "present in the
    # dataset" is derived live via gene_set_present()/gene_set_absent(). Unlike
    # `palette` (keyed by column name -> persists), gene sets are keyed by feature
    # id, so state_load() CLEARS them (meaningless on a new object).
    gene_sets = list()
  )
  # A plain environment so writing cache entries does not trigger reactivity
  # (avoids reactive-write-in-reactive churn); staleness is keyed on data_version.
  state$derived <- new.env(parent = emptyenv())
  state
}

# Derived caches that are *content-addressed* (validated by dds_content_fingerprint
# on read: the DESeq fit and the blind VST), so a structure-preserving mutation
# such as adding a normalized assay does not throw them away. state_mutate keeps
# these; a full reset/undo clears everything.
.content_derived_keys <- c("de_fit", "vst")

.clear_derived <- function(state, keep = character(0)) {
  env <- state$derived
  drop <- setdiff(ls(envir = env, all.names = TRUE), keep)
  if (length(drop)) rm(list = drop, envir = env)
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
  state$design_version <- 0L
  state$de           <- list(contrasts = list(), results = list(), active = NULL,
                             stamp = NULL, shrink = "apeglm")
  state$gene_sets    <- list()   # feature-id-keyed -> meaningless on a new object
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
  # Keep the content-addressed fit/VST: they self-invalidate by fingerprint on
  # read, so an assay-add (which leaves counts/samples/features/size factors
  # unchanged) preserves them, while a real structural change stales them.
  .clear_derived(state, keep = .content_derived_keys)
  state$n_edits      <- state$n_edits + 1L
  state$data_version <- state$data_version + 1L
  .log(state, action)
  invisible(state)
}

#' Set the model design (and optional reference-level relevels) -- design-scoped
#'
#' Updates `design(working)` and any factor reference levels **without** bumping
#' `data_version`, so design-independent caches (PCA/VST/QC) survive; bumps a
#' separate `design_version` that the DE fit keys on. Logs to history. Not pushed
#' onto the undo stack (the design is directly re-editable). The design formula
#' does not change sample/feature values, and a relevel only reorders factor
#' levels, so nothing downstream except DE actually depends on it. **Invariant:** no
#' `derived` cache key may depend on factor level order except `de_fit` (which keys on
#' `design_version`); if a future artifact orders/colours by a reference level, revisit this.
#'
#' @param state App-state object.
#' @param design A one-sided model formula (e.g. `~ condition + batch`).
#' @param relevel Optional named list `col = ref` of reference levels to set
#'   (via [de_relevel()]) before applying the design.
#' @param action Extra named fields to merge into the history entry.
#' @return `state`, invisibly.
#' @export
state_set_design <- function(state, design, relevel = NULL, action = list()) {
  dds <- state$working
  if (is.null(dds)) stop("No dataset loaded; cannot set a design.", call. = FALSE)
  for (col in names(relevel)) dds <- de_relevel(dds, col, relevel[[col]])
  DESeq2::design(dds) <- design
  state$working        <- dds
  state$design_version <- (state$design_version %||% 0L) + 1L
  .log(state, utils::modifyList(
    list(action = "set_design",
         design = paste(deparse(design), collapse = " "),
         reference = if (length(relevel)) relevel else NULL), action))
  invisible(state)
}

#' Content fingerprint of a dds's fit inputs
#'
#' The DESeq2 fit and the (blind) VST depend on the count matrix, the sample and
#' feature sets, and the size factors -- not on which *derived* assays (CPM/TPM/
#' FPKM) happen to be present. Raw counts are immutable in this app (edits only
#' ever subset), so `rownames` + `colnames` fully capture "the counts changed",
#' and `sizeFactors` captures a size-factor change directly (e.g. a future
#' control-gene re-estimate that moves them with the counts unchanged). Keying
#' fit/VST staleness on this fingerprint -- rather than the coarse `data_version`
#' that *any* mutation bumps -- means an assay-add preserves them while a real
#' structural or size-factor change invalidates them.
#'
#' @param dds A `DESeqDataSet` (or `NULL`).
#' @return A list `list(rn, cn, sf)`, or `NULL` when `dds` is `NULL`. `sf` is the
#'   **effective** size factors: the stored ones when present, else the endogenous
#'   estimate the fit/normalization would compute -- so materializing that same
#'   estimate later (e.g. adding an assay) is not mistaken for a change, while a
#'   genuinely different size-factor vector is.
#' @export
dds_content_fingerprint <- function(dds) {
  if (is.null(dds)) return(NULL)
  sf <- tryCatch(DESeq2::sizeFactors(dds), error = function(e) NULL)
  if (is.null(sf))                                    # NULL == "will estimate endogenous"
    sf <- tryCatch(DESeq2::sizeFactors(estimate_size_factors_endogenous(dds)),
                   error = function(e) NULL)
  list(rn = rownames(dds), cn = colnames(dds), sf = unname(sf))
}

# The (content-fingerprint, design_version) stamp the DESeq fit + its extracted
# results key on. Design enters via design_version (a relevel/design change must
# invalidate the fit even though the fingerprint is unchanged).
.de_stamp <- function(state) {
  list(fp = dds_content_fingerprint(state$working),
       desv = state$design_version %||% 0L)
}

#' DE result staleness relative to the current data + design
#' @param state App-state object.
#' @return `"none"` (no results), `"current"`, or `"stale"`.
#' @export
de_status <- function(state) {
  de <- state$de %||% list()
  if (is.null(de$stamp) || !length(de$results)) return("none")
  if (identical(de$stamp, .de_stamp(state))) "current" else "stale"
}

#' DESeq2 *fit* status (independent of result extraction)
#'
#' Unlike [de_status()] (which reports on extracted result tables), this reports
#' whether a `DESeq()` fit exists and is up to date: `"none"` (never run on this
#' dataset), `"stale"` (run, but data/design changed since), or `"current"`.
#' Drives the note above the Run DESeq2 button.
#' @param state App-state object.
#' @return `"none"`, `"stale"`, or `"current"`.
#' @export
de_fit_status <- function(state) {
  de <- state$de %||% list()
  if (is.null(de$stamp)) return("none")
  if (identical(de$stamp, .de_stamp(state))) "current" else "stale"
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
#' @param version The staleness token the entry is validated against (default
#'   `state$data_version`). Pass a content fingerprint (e.g.
#'   [dds_content_fingerprint()]) for a *content-addressed* artifact that should
#'   survive an assay-add -- pair it with a key in `.content_derived_keys` so
#'   [state_mutate()] does not clear it.
#' @return The cached or freshly computed value.
#' @export
state_derive <- function(state, key, params = list(), expr,
                         version = state$data_version) {
  if (exists(key, envir = state$derived, inherits = FALSE)) {
    hit <- get(key, envir = state$derived)
    if (identical(hit$version, version) && identical(hit$params, params)) {
      return(hit$value)
    }
  }
  value <- expr()
  assign(key, list(value = value, version = version, params = params), envir = state$derived)
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
    data_version = state$data_version,
    de_status    = de_status(state),
    de_n_results = length((state$de %||% list())$results)   # result tables, not stored specs
  )
}
