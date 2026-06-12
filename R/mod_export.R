# Page 7: Data export
# Export the processed dds object, DESeq2 result tables, and plots.

mod_export_ui <- function(id) {
  ns <- NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      title = "Export",
      downloadButton(ns("dds"), "Download dds (.rds)")
      # TODO: result tables, plot exports
    ),
    bslib::card(
      bslib::card_header("Export"),
      "Download the processed object, DE tables, and plots."
    )
  )
}

#' @param state the shared app-state object (see [new_app_state()]).
#' @return Invisible NULL.
mod_export_server <- function(id, state) {
  moduleServer(id, function(input, output, session) {
    output$dds <- downloadHandler(
      filename = function() paste0("dds-", Sys.Date(), ".rds"),
      content = function(file) {
        req(state$working)
        saveRDS(state$working, file)
      }
    )
    invisible(NULL)
  })
}
