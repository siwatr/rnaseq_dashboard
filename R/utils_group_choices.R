# Shared "colour / group by" selectInput choices for the plot pages (PCA, QC).
# A single grouped (<optgroup>) layout so the pages agree, in a fixed order:
#   General      - an optional "(none)" entry
#   This session - session-derived items (QC removal pool / flags, gene
#                  expression, QC metrics) passed in by the caller
#   Data metadata - the colData columns
# Pure (no Shiny), so it is unit-testable and reusable.

# The shared "This session" removal options (QC removal pool / flag status
# promoted to shared state, see new_app_state()), used as `session_items` by the
# PCA and QC per-sample colour/shape selectors.
.session_removal_items <- c("Suggested removal" = "__removal__",
                            "In removal pool"   = "__pool__")

#' Grouped colour/group-by choices for a plot-page selector
#'
#' @param coldata_cols Character vector of `colData` column names (the
#'   "Data metadata" optgroup).
#' @param session_items Named character vector (label -> value) for the
#'   "This session" optgroup (e.g. `c("Suggested removal" = "__removal__")`), or
#'   `NULL`/empty to omit it.
#' @param none Include a "General" optgroup holding a single "(none)" =
#'   `"__none__"` entry (default `TRUE`). Set `FALSE` for selectors that must
#'   always pick a real grouping.
#' @param spike_items Named character vector for the "Spike-in" optgroup (the
#'   per-sample dose-response metrics), or `NULL`/empty to omit it. Placed
#'   between "This session" and "Data metadata".
#' @return A named list suitable for `selectInput(choices = )`, each element an
#'   `<optgroup>`. Empty optgroups are dropped. When the only group would be
#'   "Data metadata" (no `none`, no `session_items`, no `spike_items`), the bare
#'   named vector is returned instead, so a plain colData selector shows no
#'   redundant optgroup.
#' @export
group_field_choices <- function(coldata_cols, session_items = NULL, none = TRUE,
                                spike_items = NULL) {
  groups <- list()
  if (isTRUE(none)) groups[["General"]] <- c("(none)" = "__none__")
  if (length(session_items)) groups[["This session"]] <- session_items
  if (length(spike_items)) groups[["Spike-in"]] <- spike_items
  if (length(coldata_cols))
    groups[["Data metadata"]] <- stats::setNames(as.character(coldata_cols),
                                                 as.character(coldata_cols))
  if (length(groups) == 1L && identical(names(groups), "Data metadata"))
    return(groups[["Data metadata"]])
  groups
}
