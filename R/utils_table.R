# The standard read-only display table for the dashboard. Every non-editable
# table (QC metrics, upload previews, annotation previews) should render through
# dt_table() so they share one look + behaviour: per-column filters, a
# rows-per-page selector, a search box, and horizontal scroll. The editable
# metadata editor (R/mod_meta_editor.R) is the deliberate exception - it carries
# its own editable/selection config.

#' Standard read-only table for the dashboard
#'
#' A thin wrapper over [DT::datatable()] applying the repo conventions for
#' read-only display tables: per-column filters (`filter = "top"`), a
#' rows-per-page selector (`lengthMenu`), a search box, horizontal scrolling,
#' and hidden row names.
#'
#' @param df A data frame to display.
#' @param page_length Initial number of rows per page.
#' @param length_menu Row-count choices offered in the rows-per-page selector.
#' @param filter DT column-filter placement; `"top"` (default) or `"none"`.
#' @param scrollX Enable horizontal scrolling for wide tables.
#' @param selection Row-selection mode passed to [DT::datatable()]. Defaults to
#'   `"none"` so display tables are read-only (DT's own default is `"multiple"`,
#'   which makes rows selectable); pass e.g. `list(mode = "multiple")` for the
#'   actionable filtering tables.
#' @param options Extra DataTables options, merged over (and overriding) the
#'   defaults (`dom`, `pageLength`, `lengthMenu`, `scrollX`).
#' @param ... Further arguments passed to [DT::datatable()].
#' @return A `datatables` htmlwidget.
#' @export
dt_table <- function(df, page_length = 10L, length_menu = c(10, 25, 50, 100),
                     filter = "top", scrollX = TRUE, selection = "none",
                     options = list(), ...) {
  defaults <- list(
    dom        = "lftip",   # length, filter (search), table, info, pagination
    pageLength = page_length,
    lengthMenu = length_menu,
    scrollX    = scrollX
  )
  DT::datatable(df, rownames = FALSE, filter = filter, selection = selection,
                options = utils::modifyList(defaults, options), ...)
}
