# Page 1: Input data
# Load a dataset from one of three sources (demo / .rds / counts + sample sheet),
# normalize it (sce->dds coercion, feature_class, logcounts), and make it the
# canonical working object via state_load(). Shows a summary, a colData preview,
# and a user-correctable feature-type guess.
#
# Later rounds add: colData/rowData editing, annotation (OrgDb/GTF), pseudobulk.

# Read a user-uploaded table by file extension (csv/tsv/txt/xlsx/xls).
.read_user_table <- function(path, name) {
  ext <- tolower(tools::file_ext(name))
  switch(ext,
    csv         = as.data.frame(readr::read_csv(path, show_col_types = FALSE)),
    tsv         = as.data.frame(readr::read_tsv(path, show_col_types = FALSE)),
    txt         = as.data.frame(readr::read_tsv(path, show_col_types = FALSE)),
    xlsx        = as.data.frame(readxl::read_excel(path)),
    xls         = as.data.frame(readxl::read_excel(path)),
    stop("Unsupported file type '.", ext, "'. Use CSV, TSV, or XLSX.", call. = FALSE)
  )
}

mod_input_ui <- function(id) {
  ns <- NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      title = "Load data",
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
      hr(),
      conditionalPanel(
        sprintf("output['%s']", ns("loaded")),
        selectInput(ns("feature_type"), "Feature unit",
                    choices = c("gene", "transcript", "exon", "feature")),
        helpText("What each row represents; used to label the app and find name columns.")
      )
    ),
    bslib::card(
      bslib::card_header("Dataset summary"),
      verbatimTextOutput(ns("summary"))
    )
  )
}

#' @param state the shared app-state object (see [new_app_state()])
#' @return Invisible NULL.
mod_input_server <- function(id, state) {
  moduleServer(id, function(input, output, session) {

    output$loaded <- reactive(!is.null(state$working))
    outputOptions(output, "loaded", suspendWhenHidden = FALSE)

    observeEvent(input$load, {
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

      res <- tryCatch({
        meta <- input_meta(obj)
        dds  <- ensure_logcounts(ensure_feature_class(as_input_dds(obj)))
        ft   <- detect_feature_type(dds)
        meta$feature_type <- ft$feature_type
        list(dds = dds, meta = meta, ft = ft)
      }, error = function(e) { showNotification(conditionMessage(e), type = "error", duration = NULL); NULL })
      req(res)

      state_load(state, res$dds, source = input$source, meta = res$meta)
      updateSelectInput(session, "feature_type", selected = res$meta$feature_type)
      msg <- sprintf("Loaded %d features x %d samples (%s).",
                     nrow(res$dds), ncol(res$dds), res$meta$data_type)
      if (isTRUE(res$meta$sce_per_cell)) {
        msg <- paste(msg, "Single-cell coerced per-cell: DESeq2 results will be",
                     "statistically unreliable; pseudobulk support is coming.")
      }
      showNotification(msg, type = if (isTRUE(res$meta$sce_per_cell)) "warning" else "message")
    })

    # Feature-type correction is a labeling change, not a data edit — don't bump
    # data_version; just update app-level meta.
    observeEvent(input$feature_type, {
      req(state$working)
      state$meta <- utils::modifyList(state$meta, list(feature_type = input$feature_type))
    }, ignoreInit = TRUE)

    output$summary <- renderPrint({
      if (is.null(state$working)) cat("No dataset loaded.") else methods::show(state$working)
    })

    invisible(NULL)
  })
}
