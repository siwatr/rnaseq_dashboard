# Page 1, "Sample info" tab: a navset of two tabs sharing one draft -
#   "Sample Metadata"     - the editable colData table (R/mod_meta_editor.R), and
#   "Additional Metadata" - upload a sample sheet and bind it onto the same draft.
# The merge composes onto the editor draft (via editor$set) and commits on Save in
# the Sample Metadata tab, mirroring how Feature info composes annotation.

.sample_editor_opts <- list(
  slot = "colData", title = tags$h4("Sample metadata", class = "fs-6"), row_noun = "sample",
  allow_row_rename = TRUE, bulk_class = FALSE
)

mod_metadata_ui <- function(id) {
  ns <- NS(id)
  bslib::navset_card_pill(
    title = tags$h4("Sample info", class = "fs-6 mb-0 pe-3"),
    bslib::nav_panel("Sample Metadata", meta_editor_ui(ns("editor"), .sample_editor_opts)),
    bslib::nav_panel(
      "Additional Metadata",
      bslib::layout_sidebar(
        sidebar = bslib::sidebar(
          title = "Bind a sample sheet", width = 320,
          fileInput(ns("sheet"), "Sample sheet",
                    accept = c(".csv", ".tsv", ".txt", ".xlsx", ".xls")),
          textInput(ns("id_col"), "Sample-id column (blank = row names / first col)", ""),
          actionButton(ns("merge"), "Bind into draft", class = "btn-primary"),
          helpText("Joins by sample id; switch to Sample Metadata and Save to keep.")
        ),
        uiOutput(ns("sheet_cov")),
        tags$small(class = "text-muted", "Uploaded sheet preview:"),
        DT::DTOutput(ns("sheet_preview"))
      )
    )
  )
}

#' @param state the shared app-state object (see [new_app_state()]).
#' @return Invisible NULL.
mod_metadata_server <- function(id, state) {
  moduleServer(id, function(input, output, session) {
    editor <- meta_editor_server("editor", state, .sample_editor_opts)
    sheet  <- reactiveVal(NULL)

    observeEvent(input$sheet, {
      tbl <- tryCatch(.read_user_table(input$sheet$datapath, input$sheet$name),
                      error = function(e) { showNotification(conditionMessage(e), type = "error"); NULL })
      sheet(tbl)
    })

    output$sheet_preview <- DT::renderDT({
      tbl <- sheet(); req(tbl)
      DT::datatable(tbl, rownames = FALSE,
                    options = list(dom = "ltp", pageLength = 10, scrollX = TRUE,
                                   lengthMenu = list(c(10, 25, 50, 100),
                                                     c("10", "25", "50", "100"))))
    })

    # How many of the dataset's samples are present in the uploaded sheet.
    output$sheet_cov <- renderUI({
      d <- editor$draft(); tbl <- sheet(); req(d, tbl)
      id_col <- if (nzchar(input$id_col)) input$id_col else NULL
      keys <- tryCatch(.table_ids(tbl, id_col), error = function(e) NULL)
      req(keys)
      ids <- colnames(d)
      .coverage_banner(sum(ids %in% keys), length(ids), "samples", "the uploaded sheet")
    })

    observeEvent(input$merge, {
      req(editor$draft(), sheet())
      id_col <- if (nzchar(input$id_col)) input$id_col else NULL
      res <- tryCatch(merge_sample_metadata(editor$draft(), sheet(), id_col),
                      error = function(e) { showNotification(conditionMessage(e), type = "error", duration = NULL); NULL })
      req(res)
      editor$set(res$dds)
      ow <- res$report$overwritten
      msg <- sprintf("Bound %d sample(s); %d upload row(s) matched nothing.",
                     res$report$matched, length(res$report$unmatched_in_table))
      if (length(ow)) msg <- paste0(msg, " Overwrote: ", paste(ow, collapse = ", "), ".")
      showNotification(paste(msg, "Switch to Sample Metadata and Save to keep."), type = "warning")
    })

    invisible(NULL)
  })
}
