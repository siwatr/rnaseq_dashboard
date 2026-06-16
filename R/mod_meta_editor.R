# Reusable draft editor for sample (colData) or feature (rowData) metadata.
# Edits accumulate in a local draft DESeqDataSet (cell edits, add/remove/rename
# columns, optional row rename / sheet merge / bulk feature_class) and commit in
# one state_mutate() on Save. Backed by the slot-agnostic helpers in
# R/metadata_helpers.R. Parameterized by `opts`:
#   slot, title, row_noun, allow_merge, allow_row_rename, bulk_class.

# Build the display data.frame for the table: the row id (sample/feature name)
# prepended as a real first column so `filter = "top"` gives it a search box
# (DT skips the rownames column). The id name is made unique against existing
# columns, so an input column literally named like `row_noun` (e.g. "sample")
# keeps its name and the id column becomes "sample.1". Display-only: the id is
# never written back to colData/rowData, and data-column indices are unchanged
# (the id sits at display column 0, where the rowname did).
.meta_display_df <- function(df, row_noun) {
  id_name <- make.unique(c(colnames(df), row_noun))[length(colnames(df)) + 1L]
  out <- data.frame(rownames(df), df, check.names = FALSE, stringsAsFactors = FALSE)
  names(out)[1] <- id_name
  out
}

meta_editor_ui <- function(id, opts, extra_sidebar = NULL, extra_main = NULL) {
  ns <- NS(id)
  panels <- list(
    bslib::accordion_panel(
      "Columns",
      textInput(ns("add_name"), "Add column: name"),
      selectInput(ns("add_type"), "Type",
                  c("character", "numeric", "integer", "logical", "factor")),
      textInput(ns("add_default"), "Default value", ""),
      actionButton(ns("add_col"), "Add column"),
      hr(),
      selectizeInput(ns("rm_cols"), "Remove columns (multi-select)",
                     choices = character(0), multiple = TRUE),
      actionButton(ns("remove_cols"), "Delete selected"),
      hr(),
      selectInput(ns("rn_col_old"), "Rename column", choices = character(0)),
      textInput(ns("rn_col_new"), "New name"),
      actionButton(ns("rename_col"), "Rename column")
    )
  )
  if (isTRUE(opts$allow_row_rename)) {
    panels <- c(panels, list(bslib::accordion_panel(
      paste("Rename", opts$row_noun),
      selectInput(ns("rn_row_old"), opts$row_noun, choices = character(0)),
      textInput(ns("rn_row_new"), "New name"),
      actionButton(ns("rename_row"), paste("Rename", opts$row_noun))
    )))
  }
  if (isTRUE(opts$allow_merge)) {
    panels <- c(panels, list(bslib::accordion_panel(
      "Merge a sample sheet",
      fileInput(ns("sheet"), NULL, accept = c(".csv", ".tsv", ".txt", ".xlsx", ".xls")),
      textInput(ns("id_col"), "Sample-id column (blank = row names / first col)", ""),
      actionButton(ns("merge"), "Merge into draft")
    )))
  }
  if (isTRUE(opts$bulk_class)) {
    panels <- c(panels, list(bslib::accordion_panel(
      "Set feature class (bulk)",
      selectInput(ns("bulk_value"), "Class", c("endogenous", "spike_in", "exogenous")),
      actionButton(ns("bulk_apply"), "Set class on filtered rows"),
      helpText("Filter the table (e.g. type ERCC), then apply to every matching row.")
    )))
  }

  sidebar <- bslib::sidebar(
    title = paste("Edit", opts$row_noun, "info"), width = 360,
    extra_sidebar,
    uiOutput(ns("protected_note")),
    do.call(bslib::accordion, c(list(open = FALSE), panels))
  )
  toolbar <- div(
    class = "d-flex gap-2 mb-2 align-items-center",
    actionButton(ns("save"), tagList(icon("floppy-disk"), "Save"), class = "btn-primary"),
    bslib::tooltip(actionButton(ns("reset_save"), icon("arrow-rotate-left"),
                                class = "btn-outline-secondary"),
                   "Reset to last save"),
    bslib::tooltip(actionButton(ns("reset_orig"), icon("arrows-rotate"),
                                class = "btn-outline-danger"),
                   "Reset to original")
  )
  card <- bslib::card(
    bslib::card_header(paste(opts$title, "(double-click a cell to edit; filters on top)")),
    toolbar,
    DT::DTOutput(ns("table"))
  )
  do.call(bslib::layout_sidebar,
          c(Filter(Negate(is.null), list(extra_main, card)), list(sidebar = sidebar)))
}

#' @param opts editor configuration (slot/title/row_noun/allow_*/bulk_class).
#' @noRd
meta_editor_server <- function(id, state, opts) {
  moduleServer(id, function(input, output, session) {
    slot <- opts$slot
    draft  <- reactiveVal(NULL)
    redraw <- reactiveVal(0L)
    bump   <- function() redraw(redraw() + 1L)

    observeEvent(state$working, { draft(state$working); bump() }, ignoreNULL = FALSE)

    observeEvent(draft(), {
      d <- draft(); req(d)
      cols <- colnames(.meta_get(d, slot))
      prot <- protected_columns(d, slot)
      updateSelectizeInput(session, "rm_cols", choices = setdiff(cols, prot))
      updateSelectInput(session, "rn_col_old", choices = cols)
      if (isTRUE(opts$allow_row_rename)) {
        updateSelectInput(session, "rn_row_old", choices = rownames(.meta_get(d, slot)))
      }
    })

    output$protected_note <- renderUI({
      d <- draft(); req(d)
      prot <- protected_columns(d, slot)
      if (length(prot)) helpText(sprintf("Protected column(s): %s.", paste(prot, collapse = ", ")))
    })

    output$table <- DT::renderDT({
      redraw()
      d <- isolate(draft()); req(d)
      df <- as.data.frame(.meta_get(d, slot))
      display <- .meta_display_df(df, opts$row_noun)
      DT::datatable(display, rownames = FALSE,
                    filter = "top",
                    editable = list(target = "cell", disable = list(columns = 0)),
                    selection = "none",
                    options = list(pageLength = 10, scrollX = TRUE))
    }, server = TRUE)

    .apply <- function(expr) {
      tryCatch({ draft(expr); bump() },
               error = function(e) showNotification(conditionMessage(e), type = "error"))
    }

    observeEvent(input$table_cell_edit, {
      info <- input$table_cell_edit
      if (identical(info$col, 0L)) return()              # id column is read-only
      d <- draft(); req(d)
      df <- as.data.frame(.meta_get(d, slot))
      col <- colnames(df)[info$col]
      rid <- rownames(df)[info$row]
      tryCatch(
        draft(edit_meta_cell(d, slot, rid, col, info$value)),  # silent; client shows it
        error = function(e) { showNotification(conditionMessage(e), type = "error"); bump() }
      )
    })

    observeEvent(input$add_col, {
      req(draft(), nzchar(input$add_name))
      .apply(add_meta_column(draft(), slot, input$add_name, input$add_type, input$add_default))
      updateTextInput(session, "add_name", value = "")
    })

    observeEvent(input$remove_cols, {
      req(draft(), length(input$rm_cols))
      res <- remove_meta_columns(draft(), slot, input$rm_cols)
      draft(res$dds); bump()
      msg <- sprintf("Removed %d column(s).", length(res$removed))
      if (length(res$skipped)) msg <- paste0(msg, " Skipped protected: ",
                                              paste(res$skipped, collapse = ", "), ".")
      showNotification(msg, type = "message")
    })

    observeEvent(input$rename_col, {
      req(draft(), nzchar(input$rn_col_old), nzchar(input$rn_col_new))
      .apply(rename_meta_column(draft(), slot, input$rn_col_old, input$rn_col_new))
      updateTextInput(session, "rn_col_new", value = "")
    })

    if (isTRUE(opts$allow_row_rename)) observeEvent(input$rename_row, {
      req(draft(), nzchar(input$rn_row_old), nzchar(input$rn_row_new))
      .apply(rename_samples(draft(), input$rn_row_old, input$rn_row_new))
      updateTextInput(session, "rn_row_new", value = "")
    })

    if (isTRUE(opts$allow_merge)) observeEvent(input$merge, {
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
      if (length(ow)) msg <- paste0(msg, " Overwrote: ", paste(ow, collapse = ", "), ".")
      showNotification(paste(msg, "Click Save to keep."), type = "warning")
    })

    if (isTRUE(opts$bulk_class)) observeEvent(input$bulk_apply, {
      req(draft())
      d <- draft()
      rows <- input$table_rows_all                       # row indices after filtering
      ids <- rownames(.meta_get(d, slot))[rows]
      if (!length(ids)) { showNotification("No rows in the current filter.", type = "warning"); return() }
      .apply(set_feature_class(d, ids, input$bulk_value))
      showNotification(sprintf("Set %d feature(s) to %s.", length(ids), input$bulk_value),
                       type = "message")
    })

    observeEvent(input$save, {
      req(draft(), state$working)
      d <- draft()
      unchanged <-
        identical(as.data.frame(.meta_get(d, slot)),
                  as.data.frame(.meta_get(state$working, slot))) &&
        identical(colnames(d), colnames(state$working)) &&
        identical(rownames(d), rownames(state$working)) &&
        identical(deparse(DESeq2::design(d)), deparse(DESeq2::design(state$working)))
      if (unchanged) { showNotification("No changes to save.", type = "message"); return() }
      state_mutate(state, function(.) d, action = list(action = paste0("edit_", slot)))
      showNotification("Changes saved.", type = "message")
    })
    observeEvent(input$reset_save, { draft(state$working);  bump(); showNotification("Reverted to last save.") })
    observeEvent(input$reset_orig, { draft(state$original); bump(); showNotification("Reverted to original.") })

    # Expose the draft so a host page (e.g. Feature info) can compose its own
    # edits -- annotation, feature_length -- onto the same buffer instead of
    # committing straight to state$working (which would reset the draft and drop
    # the user's unsaved edits). `set()` replaces the draft and redraws; changes
    # still commit only when the user clicks Save.
    list(draft = draft, set = function(d) { draft(d); bump() })
  })
}
