# Gene Sets page -- define, record, and manage named gene sets of interest (the
# resource DE seeds and the P7 Expression heatmap consumes). Sets live in
# state$gene_sets (a session/UI field, NO data_version impact) as structured,
# NON-DESTRUCTIVE records: `ids` is the full authored membership, and "present /
# absent in the dataset" is derived live (gene_set_present/gene_set_absent).
# Sources are SNAPSHOTS -- their controls preview live, but Add freezes the ids.
#
# P6b builds the Manage tab (Paste / DE DEGs / top-variable sources + set
# management). Tabular import (P6c), file round-trip (P6d), and the Compare tab
# (P6e) land in later sub-PRs; the page is a navset_card_tab so Compare is a
# one-line nav_panel addition.

mod_geneset_ui <- function(id) {
  ns <- NS(id)
  bslib::navset_card_tab(
    id = ns("tabs"),
    bslib::nav_panel(
      tags$h3("Manage sets", class = "fs-6 mb-0"),
      bslib::layout_sidebar(
        sidebar = bslib::sidebar(
          width = 360, title = "Add to set",
          textInput(ns("set_name"), "Set name", placeholder = "e.g. my_genes"),
          radioButtons(ns("add_mode"), NULL,
                       c("New set" = "new", "Append to existing" = "append"),
                       selected = "new", inline = TRUE),
          bslib::accordion(
            id = ns("sources"), open = "Paste names/IDs", multiple = FALSE,
            bslib::accordion_panel(
              "Paste names/IDs", icon = shiny::icon("keyboard"),
              gene_search_ui(ns, "paste", multiple = TRUE,
                             search_modes = c("exact", "contains", "regex")),
              bslib::input_switch(ns("paste_literal"),
                                  "Also add unmatched as literal IDs", value = FALSE),
              actionButton(ns("add_paste"), "Add to set",
                           class = "btn-primary btn-sm", icon = shiny::icon("plus"))),
            bslib::accordion_panel(
              "From DE contrast (DEGs)", icon = shiny::icon("chart-column"),
              uiOutput(ns("deg_controls"))),
            bslib::accordion_panel(
              "From top-variable genes", icon = shiny::icon("arrow-down-wide-short"),
              numericInput(ns("topvar_n"), "Number of genes", value = 100, min = 1, step = 50),
              uiOutput(ns("topvar_note")),
              actionButton(ns("add_topvar"), "Add to set",
                           class = "btn-primary btn-sm", icon = shiny::icon("plus")))
          )
        ),
        tags$h5("Defined sets", class = "fs-6"),
        uiOutput(ns("empty_note")),
        DT::DTOutput(ns("sets_table")),
        tags$div(
          class = "d-flex gap-2 my-2",
          actionButton(ns("rename"), "Rename", class = "btn-sm btn-outline-secondary",
                       icon = shiny::icon("pen")),
          actionButton(ns("delete"), "Delete", class = "btn-sm btn-outline-danger",
                       icon = shiny::icon("trash")),
          actionButton(ns("clear_all"), "Clear all", class = "btn-sm btn-outline-danger",
                       icon = shiny::icon("xmark"))),
        tags$hr(),
        uiOutput(ns("members_header")),
        DT::DTOutput(ns("members_table"))
      )
    )
  )
}

mod_geneset_server <- function(id, state) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    feature_type <- function() (state$meta %||% list())$feature_type %||% "feature"
    working_rn <- function() if (is.null(state$working)) character(0) else rownames(state$working)

    paste_search <- gene_search_server(input, output, session, state, "paste",
                                       multiple = TRUE,
                                       search_modes = c("exact", "contains", "regex"))
    selected_set <- reactiveVal(NULL)   # the set name shown in the members table

    # --- Non-destructive absence tracking (baseline + notify-on-new) --------
    prev_absent <- reactiveVal(list())
    snapshot_absent <- function() {
      rn <- working_rn()
      prev_absent(lapply(state$gene_sets, function(s) gene_set_absent(s, rn)))
    }

    # --- The shared commit path (New / Append) ------------------------------
    do_add <- function(ids, source) {
      ids <- unique(as.character(ids)); ids <- ids[!is.na(ids) & nzchar(ids)]
      nm  <- trimws(input$set_name %||% "")
      if (!nzchar(nm)) { showNotification("Enter a set name first.", type = "warning"); return() }
      if (!length(ids)) { showNotification("No genes to add.", type = "warning"); return() }
      res <- gene_set_commit(state$gene_sets, nm, ids, input$add_mode %||% "new", source = source)
      state$gene_sets <- res$sets
      selected_set(res$name)
      snapshot_absent()   # baseline so later feature edits notify only on NEW absences
      .log(state, list(action = "gene_set_add", set = res$name, source = source,
                       mode = input$add_mode %||% "new", n = length(ids)))
      showNotification(sprintf("Added %d gene(s) to '%s'.", length(ids), res$name),
                       type = "message", duration = 4)
    }

    # --- Source: paste names/IDs -------------------------------------------
    observeEvent(input$add_paste, {
      r <- paste_search(); ids <- r$ids
      if (isTRUE(input$paste_literal)) ids <- c(ids, r$unmatched)   # authored, may be absent
      do_add(ids, "paste")
    })

    # --- Source: DE contrast DEGs ------------------------------------------
    output$deg_controls <- renderUI({
      res <- (state$de %||% list())$results
      if (!length(res))
        return(helpText(class = "small text-muted",
          "Run DESeq2 (on the DE page) to seed gene sets from a contrast's DEGs."))
      tagList(
        selectInput(ns("deg_contrast"), "Contrast", choices = names(res)),
        radioButtons(ns("deg_dir"), "Direction",
                     c("Up" = "up", "Down" = "down", "Both" = "both"),
                     selected = "both", inline = TRUE),
        bslib::layout_columns(
          col_widths = c(6, 6),
          numericInput(ns("deg_padj"), "padj <", value = 0.05, min = 0, max = 1, step = 0.01),
          numericInput(ns("deg_lfc"), "|log2FC| >=", value = 1, min = 0, step = 0.5)),
        bslib::input_switch(ns("deg_shrunk"), "Use shrunk LFC", value = FALSE),
        uiOutput(ns("deg_preview")),
        actionButton(ns("add_deg"), "Add to set",
                     class = "btn-primary btn-sm", icon = shiny::icon("plus")))
    })
    # Classified DEG table for the chosen contrast + thresholds (snapshotted on Add).
    deg_classified <- reactive({
      res <- (state$de %||% list())$results; lab <- input$deg_contrast
      req(length(res), !is.null(lab), lab %in% names(res))
      de_classify_table(res[[lab]], input$deg_padj %||% 0.05, input$deg_lfc %||% log2(2))
    })
    deg_col <- reactive({
      df <- deg_classified()
      if (isTRUE(input$deg_shrunk) && "DEG_shrunk" %in% names(df)) "DEG_shrunk" else "DEG"
    })
    deg_ids <- reactive({
      df <- deg_classified()
      want <- switch(input$deg_dir %||% "both", up = "up", down = "down",
                     both = c("up", "down"))
      rownames(df)[as.character(df[[deg_col()]]) %in% want]
    })
    output$deg_preview <- renderUI({
      df <- tryCatch(deg_classified(), error = function(e) NULL); req(df)
      tab <- table(factor(as.character(df[[deg_col()]]),
                          levels = c("up", "down", "no_change")))
      helpText(class = "small text-muted",
        sprintf("%d up / %d down at padj < %s, |log2FC| >= %s.",
                tab[["up"]], tab[["down"]], format(input$deg_padj %||% 0.05),
                format(input$deg_lfc %||% 1)))
    })
    observeEvent(input$add_deg, {
      ids <- tryCatch(deg_ids(), error = function(e) character(0))
      if (!length(ids)) {
        showNotification("No DEGs in that direction at these thresholds.", type = "warning"); return() }
      do_add(ids, paste0("DE: ", input$deg_contrast %||% ""))
    })

    # --- Source: top-variable genes ----------------------------------------
    output$topvar_note <- renderUI({
      helpText(class = "small text-muted",
        "Adds the most variable endogenous features (VST) -- snapshotted on Add.")
    })
    observeEvent(input$add_topvar, {
      req(state$working)
      n <- suppressWarnings(as.integer(input$topvar_n %||% 100L)); if (is.na(n)) n <- 100L
      n <- max(1L, n)
      ids <- tryCatch(top_variable_features(pca_input(state$working, "vst")$mat, n_top = n),
                      error = function(e) character(0))
      if (!length(ids)) { showNotification("No variable features available.", type = "warning"); return() }
      do_add(ids, "top-variable")
    })

    # --- Manage: defined-sets table + members ------------------------------
    sets_df <- reactive({
      sets <- state$gene_sets
      if (!length(sets)) return(NULL)
      rn <- working_rn()
      data.frame(
        Name  = names(sets),
        Genes = vapply(sets, function(s) length(s$ids), integer(1)),
        `Not in dataset` = vapply(sets, function(s) length(gene_set_absent(s, rn)), integer(1)),
        Source = vapply(sets, function(s) s$source %||% "", character(1)),
        check.names = FALSE, stringsAsFactors = FALSE, row.names = NULL)
    })
    output$empty_note <- renderUI({
      if (length(state$gene_sets)) return(NULL)
      helpText(class = "small text-muted", "No gene sets yet. Add one from the sidebar.")
    })
    output$sets_table <- DT::renderDT({
      df <- sets_df(); req(df)
      dt_table(df, selection = "single", page_length = 10L)
    })
    observeEvent(input$sets_table_rows_selected, {
      df <- sets_df(); i <- input$sets_table_rows_selected
      if (!is.null(df) && length(i)) selected_set(df$Name[i])
    }, ignoreNULL = FALSE)

    display_names <- function(ids) {
      if (is.null(state$working)) return(ids)
      rd <- SummarizedExperiment::rowData(state$working)
      fn <- paste0(feature_type(), "_name")
      if (!(fn %in% colnames(rd))) return(ids)
      nm  <- rep(NA_character_, length(ids))
      hit <- ids %in% rownames(state$working)
      if (any(hit)) nm[hit] <- as.character(rd[ids[hit], fn])
      ifelse(is.na(nm) | !nzchar(nm), ids, nm)
    }
    output$members_header <- renderUI({
      nm <- selected_set()
      if (is.null(nm) || is.null(state$gene_sets[[nm]]))
        return(helpText(class = "small text-muted", "Select a set above to view its members."))
      tags$h5(sprintf("Members of '%s'", nm), class = "fs-6")
    })
    output$members_table <- DT::renderDT({
      nm <- selected_set(); req(!is.null(nm), !is.null(state$gene_sets[[nm]]))
      ids <- state$gene_sets[[nm]]$ids
      df <- data.frame(ID = ids, Name = display_names(ids),
                       `In dataset` = ifelse(ids %in% working_rn(), "yes", "no"),
                       check.names = FALSE, stringsAsFactors = FALSE)
      dt_table(df, page_length = 10L)
    })

    # --- Rename / Delete / Clear -------------------------------------------
    observeEvent(input$rename, {
      nm <- selected_set(); req(!is.null(nm), !is.null(state$gene_sets[[nm]]))
      showModal(modalDialog(
        title = "Rename gene set", size = "s",
        textInput(ns("rename_to"), "New name", value = nm),
        footer = tagList(modalButton("Cancel"),
                         actionButton(ns("rename_ok"), "Rename", class = "btn-primary"))))
    })
    observeEvent(input$rename_ok, {
      nm <- selected_set(); to <- trimws(input$rename_to %||% "")
      req(!is.null(nm), !is.null(state$gene_sets[[nm]]))
      if (!nzchar(to)) { showNotification("Enter a name.", type = "warning"); return() }
      if (!identical(to, nm) && !is.null(state$gene_sets[[to]])) {
        showNotification("A set with that name already exists.", type = "warning"); return() }
      if (!identical(to, nm)) {
        sets <- state$gene_sets
        names(sets)[names(sets) == nm] <- to
        state$gene_sets <- sets
        selected_set(to)
        .log(state, list(action = "gene_set_rename", from = nm, to = to))
      }
      removeModal()
    })
    observeEvent(input$delete, {
      nm <- selected_set(); req(!is.null(nm), !is.null(state$gene_sets[[nm]]))
      sets <- state$gene_sets; sets[[nm]] <- NULL
      state$gene_sets <- sets; selected_set(NULL)
      .log(state, list(action = "gene_set_delete", set = nm))
      showNotification(sprintf("Deleted '%s'.", nm), type = "message", duration = 3)
    })
    observeEvent(input$clear_all, {
      req(length(state$gene_sets))
      showModal(modalDialog(
        title = "Clear all gene sets", size = "s",
        sprintf("Delete all %d gene set(s)? This cannot be undone.", length(state$gene_sets)),
        footer = tagList(modalButton("Cancel"),
                         actionButton(ns("clear_ok"), "Clear all", class = "btn-danger"))))
    })
    observeEvent(input$clear_ok, {
      state$gene_sets <- list(); selected_set(NULL); prev_absent(list())
      .log(state, list(action = "gene_set_clear_all"))
      removeModal()
    })

    # --- Non-destructive reconcile on a data edit: notify on NEW absences ---
    observeEvent(state$data_version, {
      sets <- state$gene_sets
      if (is.null(state$working) || !length(sets)) { prev_absent(list()); return() }
      rn <- rownames(state$working); prev <- prev_absent()
      cur <- list(); newly <- character(0)
      for (nm in names(sets)) {
        ab <- gene_set_absent(sets[[nm]], rn); cur[[nm]] <- ab
        gained <- setdiff(ab, prev[[nm]] %||% character(0))
        if (length(gained)) newly <- c(newly, sprintf("%s (+%d)", nm, length(gained)))
      }
      prev_absent(cur)
      if (length(newly))
        showNotification(
          sprintf("Gene set(s) now include features not in the dataset: %s. They are kept (see 'Not in dataset').",
                  paste(newly, collapse = ", ")),
          type = "warning", duration = 7)
    }, ignoreInit = TRUE)
  })
}
