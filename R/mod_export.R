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

#' @param dds reactive() yielding the current DESeqDataSet.
#' @return Invisible NULL.
mod_export_server <- function(id, dds) {
  moduleServer(id, function(input, output, session) {
    output$dds <- downloadHandler(
      filename = function() paste0("dds-", Sys.Date(), ".rds"),
      content = function(file) {
        req(dds())
        saveRDS(dds(), file)
      }
    )
    invisible(NULL)
  })
}
