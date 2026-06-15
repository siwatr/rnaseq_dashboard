# Page 1, "Sample info" tab: edit colData in place, merge an uploaded sample
# sheet, and tag feature_class. All edits go through state_mutate() (bumping
# data_version, logging history). Backed by R/metadata_helpers.R.

mod_metadata_ui <- function(id) {
  ns <- NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      title = "Edit metadata", width = 340,
      tags$strong("Merge a sample sheet"),
      fileInput(ns("sheet"), NULL,
                accept = c(".csv", ".tsv", ".txt", ".xlsx", ".xls")),
      textInput(ns("id_col"), "Sample-id column (blank = row names / first col)", ""),
      actionButton(ns("merge"), "Merge into sample info"),
      hr(),
      tags$strong("Tag feature class"),
      selectizeInput(ns("feat"), "Find features", choices = NULL, multiple = TRUE,
                     options = list(placeholder = "gene name or id")),
      radioButtons(ns("fclass"), "Set class to",
                   c("endogenous", "spike_in", "exogenous")),
      actionButton(ns("tag"), "Apply tag")
    ),
    bslib::card(
      bslib::card_header("Sample information (double-click a cell to edit)"),
      DT::DTOutput(ns("coldata"))
    ),
    bslib::card(
      bslib::card_header("Feature classes"),
      verbatimTextOutput(ns("fclass_summary"))
    )
  )
}

#' @param state the shared app-state object (see [new_app_state()]).
#' @return Invisible NULL.
mod_metadata_server <- function(id, state) {
  moduleServer(id, function(input, output, session) {

    coldata_df <- reactive({
      req(state$working)
      as.data.frame(SummarizedExperiment::colData(state$working))
    })

    output$coldata <- DT::renderDT(
      DT::datatable(coldata_df(), editable = "cell", selection = "none",
                    options = list(pageLength = 10, scrollX = TRUE)),
      server = TRUE
    )
    proxy <- DT::dataTableProxy("coldata")

    # In-cell edit -> coerce + apply via state_mutate; revert the cell on error.
    observeEvent(input$coldata_cell_edit, {
      info <- input$coldata_cell_edit
      df <- coldata_df()
      col <- colnames(df)[info$col]          # rownames occupy col 0; data cols are 1..p
      sample <- rownames(df)[info$row]
      tryCatch(
        state_mutate(
          state,
          function(d) edit_coldata_cell(d, sample, col, info$value),
          action = list(action = "edit_colData", column = col, row = sample)
        ),
        error = function(e) {
          showNotification(conditionMessage(e), type = "error")
          DT::replaceData(proxy, coldata_df(), resetPaging = FALSE, rownames = TRUE)
        }
      )
    })

    # Server-side feature search, labelled by <feature_type>_name when present.
    observeEvent(state$working, {
      req(state$working)
      ids <- rownames(state$working)
      name_col <- paste0(state_meta(state)$feature_type, "_name")
      rd <- SummarizedExperiment::rowData(state$working)
      labels <- if (name_col %in% colnames(rd)) paste0(rd[[name_col]], " (", ids, ")") else ids
      updateSelectizeInput(session, "feat",
                           choices = stats::setNames(ids, labels), server = TRUE)
    })

    observeEvent(input$merge, {
      req(state$working, input$sheet)
      tbl <- tryCatch(.read_user_table(input$sheet$datapath, input$sheet$name),
                      error = function(e) { showNotification(conditionMessage(e), type = "error"); NULL })
      req(tbl)
      id_col <- if (nzchar(input$id_col)) input$id_col else NULL
      res <- tryCatch(merge_sample_metadata(state$working, tbl, id_col),
                      error = function(e) { showNotification(conditionMessage(e), type = "error"); NULL })
      req(res)
      state_mutate(state, function(d) res$dds, action = list(action = "merge_colData"))
      showNotification(
        sprintf("Merged metadata for %d sample(s); %d row(s) in the upload matched no sample.",
                res$report$matched, length(res$report$unmatched_in_table)),
        type = "message"
      )
    })

    observeEvent(input$tag, {
      req(state$working, input$feat)
      cls <- input$fclass
      ids <- input$feat
      state_mutate(state, function(d) set_feature_class(d, ids, cls),
                   action = list(action = "set_feature_class", class = cls, n = length(ids)))
      showNotification(sprintf("Tagged %d feature(s) as %s.", length(ids), cls), type = "message")
    })

    output$fclass_summary <- renderPrint({
      req(state$working)
      print(table(SummarizedExperiment::rowData(state$working)$feature_class))
    })

    invisible(NULL)
  })
}
