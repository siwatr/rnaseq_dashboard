# Page 1, "Sample info" tab: a draft editor over colData. Edits (cell changes,
# add/remove/rename column, rename sample, merge an uploaded sheet) accumulate in
# a local draft DESeqDataSet and are committed in one state_mutate() on Save.
# Reset discards the draft (to last save, or to the originally loaded object).
# Backed by R/metadata_helpers.R. Feature metadata lives in the Feature-info tab.

mod_metadata_ui <- function(id) {
  ns <- NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      title = "Edit sample info", width = 360,
      div(
        actionButton(ns("save"), "Save changes", class = "btn-primary"),
        actionButton(ns("reset_save"), "Reset to last save", class = "btn-outline-secondary"),
        actionButton(ns("reset_orig"), "Reset to original", class = "btn-outline-danger")
      ),
      uiOutput(ns("protected_note")),
      bslib::accordion(
        open = FALSE,
        bslib::accordion_panel(
          "Columns",
          textInput(ns("add_name"), "Add column: name"),
          selectInput(ns("add_type"), "Type",
                      c("character", "numeric", "integer", "logical", "factor")),
          textInput(ns("add_default"), "Default value", ""),
          actionButton(ns("add_col"), "Add column"),
          hr(),
          selectInput(ns("rm_col"), "Remove column", choices = character(0)),
          actionButton(ns("remove_col"), "Remove"),
          hr(),
          selectInput(ns("rn_col_old"), "Rename column", choices = character(0)),
          textInput(ns("rn_col_new"), "New name"),
          actionButton(ns("rename_col"), "Rename column")
        ),
        bslib::accordion_panel(
          "Rename sample",
          selectInput(ns("rn_samp_old"), "Sample", choices = character(0)),
          textInput(ns("rn_samp_new"), "New name"),
          actionButton(ns("rename_samp"), "Rename sample")
        ),
        bslib::accordion_panel(
          "Merge a sample sheet",
          fileInput(ns("sheet"), NULL,
                    accept = c(".csv", ".tsv", ".txt", ".xlsx", ".xls")),
          textInput(ns("id_col"), "Sample-id column (blank = row names / first col)", ""),
          actionButton(ns("merge"), "Merge into draft")
        )
      )
    ),
    bslib::card(
      bslib::card_header("Sample information (double-click a cell to edit; filters on top)"),
      DT::DTOutput(ns("coldata"))
    )
  )
}

#' @param state the shared app-state object (see [new_app_state()]).
#' @return Invisible NULL.
mod_metadata_server <- function(id, state) {
  moduleServer(id, function(input, output, session) {

    draft  <- reactiveVal(NULL)
    redraw <- reactiveVal(0L)
    bump   <- function() redraw(redraw() + 1L)

    # (Re)initialize the draft from the canonical object whenever it changes
    # (load, Save, or an edit from another tab). Keeps the two consistent.
    observeEvent(state$working, { draft(state$working); bump() }, ignoreNULL = FALSE)

    # Keep column/sample selectors in sync with the draft.
    observeEvent(draft(), {
      d <- draft(); req(d)
      cols <- colnames(SummarizedExperiment::colData(d))
      prot <- protected_columns(d)
      updateSelectInput(session, "rm_col",      choices = setdiff(cols, prot))
      updateSelectInput(session, "rn_col_old",  choices = cols)
      updateSelectInput(session, "rn_samp_old", choices = colnames(d))
    })

    output$protected_note <- renderUI({
      d <- draft(); req(d)
      prot <- protected_columns(d)
      if (length(prot)) {
        helpText(sprintf("Design column(s) cannot be removed: %s.", paste(prot, collapse = ", ")))
      }
    })

    output$coldata <- DT::renderDT({
      redraw()
      d <- isolate(draft())
      req(d)
      DT::datatable(as.data.frame(SummarizedExperiment::colData(d)),
                    filter = "top", editable = "cell", selection = "none",
                    options = list(pageLength = 10, scrollX = TRUE))
    }, server = TRUE)

    # --- structural edits on the draft (each bumps the table redraw) ---
    .apply <- function(expr) {
      tryCatch({ draft(expr); bump() },
               error = function(e) showNotification(conditionMessage(e), type = "error"))
    }

    observeEvent(input$coldata_cell_edit, {
      info <- input$coldata_cell_edit
      d <- draft(); req(d)
      df <- as.data.frame(SummarizedExperiment::colData(d))
      col <- colnames(df)[info$col]
      sample <- rownames(df)[info$row]
      tryCatch(
        draft(edit_coldata_cell(d, sample, col, info$value)),  # silent; client shows it
        error = function(e) { showNotification(conditionMessage(e), type = "error"); bump() }
      )
    })

    observeEvent(input$add_col, {
      req(draft(), nzchar(input$add_name))
      .apply(add_coldata_column(draft(), input$add_name, input$add_type, input$add_default))
      updateTextInput(session, "add_name", value = "")
    })
    observeEvent(input$remove_col, {
      req(draft(), nzchar(input$rm_col))
      .apply(remove_coldata_column(draft(), input$rm_col))
    })
    observeEvent(input$rename_col, {
      req(draft(), nzchar(input$rn_col_old), nzchar(input$rn_col_new))
      .apply(rename_coldata_column(draft(), input$rn_col_old, input$rn_col_new))
      updateTextInput(session, "rn_col_new", value = "")
    })
    observeEvent(input$rename_samp, {
      req(draft(), nzchar(input$rn_samp_old), nzchar(input$rn_samp_new))
      .apply(rename_samples(draft(), input$rn_samp_old, input$rn_samp_new))
      updateTextInput(session, "rn_samp_new", value = "")
    })

    observeEvent(input$merge, {
      req(draft(), input$sheet)
      tbl <- tryCatch(.read_user_table(input$sheet$datapath, input$sheet$name),
                      error = function(e) { showNotification(conditionMessage(e), type = "error"); NULL })
      req(tbl)
      id_col <- if (nzchar(input$id_col)) input$id_col else NULL
      res <- tryCatch(merge_sample_metadata(draft(), tbl, id_col),
                      error = function(e) { showNotification(conditionMessage(e), type = "error"); NULL })
      req(res)
      draft(res$dds); bump()
      ow <- res$report$overwritten
      msg <- sprintf("Merged %d sample(s); %d upload row(s) matched nothing.",
                     res$report$matched, length(res$report$unmatched_in_table))
      if (length(ow)) msg <- paste0(msg, " Overwrote column(s): ", paste(ow, collapse = ", "), ".")
      showNotification(paste(msg, "Click Save to keep."), type = "warning")
    })

    # --- commit / reset ---
    observeEvent(input$save, {
      req(draft(), state$working)
      d <- draft()
      unchanged <-
        identical(as.data.frame(SummarizedExperiment::colData(d)),
                  as.data.frame(SummarizedExperiment::colData(state$working))) &&
        identical(colnames(d), colnames(state$working)) &&
        identical(deparse(DESeq2::design(d)), deparse(DESeq2::design(state$working)))
      if (unchanged) { showNotification("No changes to save.", type = "message"); return() }
      state_mutate(state, function(.) d, action = list(action = "edit_metadata"))
      showNotification("Sample info saved.", type = "message")
    })
    observeEvent(input$reset_save, { draft(state$working);  bump(); showNotification("Reverted to last save.") })
    observeEvent(input$reset_orig, { draft(state$original); bump(); showNotification("Reverted to original sample info.") })

    invisible(NULL)
  })
}
