# Page 1: Input data
# Load a dataset from one of three sources (demo / .rds / counts + sample sheet),
# normalize it (sce->dds coercion, feature_class, logcounts), and make it the
# canonical working object via state_load(). Shows a summary, a colData preview,
# and a user-correctable feature-type guess.
#
# Later rounds add: colData/rowData editing, annotation (OrgDb/GTF), pseudobulk.

mod_input_ui <- function(id) {
  ns <- NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      title = tags$h4("Load data", class="fs-6 mb-0"),
      radioButtons(ns("source"), "Source", inline = FALSE,
        choices = c("Demo dataset" = "demo", "R object (.rds)" = "rds",
                    "Counts + sample sheet" = "tabular")),
      conditionalPanel(
        sprintf("input['%s'] == 'rds'", ns("source")),
        fileInput(ns("rds"), "DESeqDataSet or SingleCellExperiment (.rds)", accept = ".rds")
      ),
      conditionalPanel(
        sprintf("input['%s'] == 'tabular'", ns("source")),
        fileInput(ns("counts"), "Counts table (features x samples)",
                  accept = c(".csv", ".tsv", ".txt", ".xlsx", ".xls")),
        fileInput(ns("samples"), "Sample sheet (one row per sample)",
                  accept = c(".csv", ".tsv", ".txt", ".xlsx", ".xls"))
      ),
      actionButton(ns("load"), "Load dataset", class = "btn-primary"),
      helpText("Set the feature unit (gene/transcript/...) on the Feature info tab.")
    ),
    bslib::card(
      bslib::card_header(tags$h4("Dataset summary", class = "fs-6")),
      verbatimTextOutput(ns("summary"))
    )
  )
}

#' @param state the shared app-state object (see [new_app_state()])
#' @return Invisible NULL.
mod_input_server <- function(id, state) {
  moduleServer(id, function(input, output, session) {

    observeEvent(input$load, {
      withProgress(message = "Loading dataset...", value = 0.1, {
      obj <- tryCatch(
        switch(input$source,
          demo    = make_mock_dds(),
          rds     = {
            validate(need(!is.null(input$rds), "Upload an .rds file first."))
            readRDS(input$rds$datapath)
          },
          tabular = {
            validate(need(!is.null(input$counts) && !is.null(input$samples),
                          "Upload both a counts table and a sample sheet."))
            tabular_to_dds(
              .read_user_table(input$counts$datapath, input$counts$name),
              .read_user_table(input$samples$datapath, input$samples$name)
            )
          }
        ),
        error = function(e) { showNotification(conditionMessage(e), type = "error", duration = NULL); NULL }
      )
      req(obj)

      setProgress(0.6, message = "Preparing object (feature class, logcounts)...")
      res <- tryCatch({
        meta <- input_meta(obj)
        dds  <- ensure_logcounts(ensure_feature_class(as_input_dds(obj)))
        ft   <- detect_feature_type(dds)
        meta$feature_type <- ft$feature_type
        list(dds = dds, meta = meta, ft = ft)
      }, error = function(e) { showNotification(conditionMessage(e), type = "error", duration = NULL); NULL })
      req(res)

      state_load(state, res$dds, source = input$source, meta = res$meta)
      msg <- sprintf("Loaded %d features x %d samples (%s).",
                     nrow(res$dds), ncol(res$dds), res$meta$data_type)
      if (isTRUE(res$meta$sce_per_cell)) {
        msg <- paste(msg, "Single-cell coerced per-cell: DESeq2 results will be",
                     "statistically unreliable; pseudobulk support is coming.")
      }
      showNotification(msg, type = if (isTRUE(res$meta$sce_per_cell)) "warning" else "message")
      })
    })

    output$summary <- renderPrint({
      if (is.null(state$working)) cat("No dataset loaded.") else methods::show(state$working)
    })

    invisible(NULL)
  })
}
