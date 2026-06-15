# Page 1, "Sample info" tab: the shared metadata draft editor configured for
# colData (samples) — adds the sheet-merge and sample-rename controls. All
# editing logic lives in R/mod_meta_editor.R.

.sample_editor_opts <- list(
  slot = "colData", title = "Sample information", row_noun = "sample",
  allow_merge = TRUE, allow_row_rename = TRUE, bulk_class = FALSE
)

mod_metadata_ui <- function(id) {
  ns <- NS(id)
  meta_editor_ui(ns("editor"), .sample_editor_opts)
}

#' @param state the shared app-state object (see [new_app_state()]).
#' @return Invisible NULL.
mod_metadata_server <- function(id, state) {
  moduleServer(id, function(input, output, session) {
    meta_editor_server("editor", state, .sample_editor_opts)
    invisible(NULL)
  })
}
