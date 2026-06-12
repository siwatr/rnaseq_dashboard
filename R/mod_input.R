# Page 1: Input data
# Load a dds/sce (.rds), convert sce -> dds, show summary, edit colData,
# map rowData/rowRanges from a GTF, tag exogenous / spike-in features.
# This module *produces* the canonical dds for the rest of the app.

mod_input_ui <- function(id) {
  ns <- NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      title = "Load data",
      fileInput(ns("file"), "DESeqDataSet or SingleCellExperiment (.rds)",
                accept = ".rds")
      # TODO: metadata upload, GTF upload, feature-type / name-column pickers
    ),
    bslib::card(
      bslib::card_header("Dataset summary"),
      verbatimTextOutput(ns("summary"))
    )
  )
}

#' @return reactive() yielding the loaded DESeqDataSet (NULL until loaded).
mod_input_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    dds <- reactiveVal(NULL)

    observeEvent(input$file, {
      obj <- readRDS(input$file$datapath)
      # TODO: if SingleCellExperiment, convert to DESeqDataSet (design = ~ 1)
      # TODO: auto-add logcounts assay (see rnaseq-bioc skill)
      dds(obj)
    })

    output$summary <- renderPrint({
      if (is.null(dds())) cat("No dataset loaded.") else methods::show(dds())
    })

    reactive(dds())
  })
}
