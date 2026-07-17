# Gene Sets page -- define, record, and manage named gene sets of interest (the
# resource DE seeds and the P7 Expression heatmap consumes). Sets live in
# state$gene_sets (a session/UI field, NO data_version impact) as structured,
# NON-DESTRUCTIVE records: `ids` is the full authored membership, and "present /
# absent in the dataset" is derived live (gene_set_present/gene_set_absent).
#
# The Manage tab is a STAGING workflow (git-metaphor on the backend, friendly
# labels in the UI): a "Build a gene set" card on the left (source pills -> a live
# Preview of the *staged* genes -> a New/Add Save) beside a "Your gene sets" card
# on the right (the session store). `staged()` = a named list of incoming
# id-vectors that follows the ACTIVE source pill; single-set sources yield a
# 1-element list, so the P6c file import (N staged sets) reuses the same layout.
#
# Tabular import (P6c), file round-trip (P6d), and the Compare tab (P6e) land in
# later sub-PRs; the page is a navset_card_tab so Compare is a one-line addition.

# A segmented "pill" control (shinyWidgets when available, radioButtons fallback).
.gs_pills <- function(id, choices, selected = NULL) {
  if (requireNamespace("shinyWidgets", quietly = TRUE))
    shinyWidgets::radioGroupButtons(id, label = NULL, choices = choices,
                                    selected = selected, size = "sm", status = "primary")
  else
    radioButtons(id, label = NULL, choices = choices, selected = selected, inline = TRUE)
}

# The app "primary action" purple, matching the Run DESeq2 button (mod_de.R).
.gs_action_style <- "background-color:#8b58db;border-color:#8b58db;color:#fff;"

mod_geneset_ui <- function(id) {
  ns <- NS(id)
  cond <- function(inp, val) sprintf("input['%s'] == '%s'", ns(inp), val)
  bslib::navset_card_tab(
    id = ns("tabs"),
    bslib::nav_panel(
      tags$h3("Manage sets", class = "fs-6 mb-0"),
      bslib::layout_columns(
        fillable = FALSE, col_widths = c(6, 6),

        # ---- LEFT: Build a gene set (title text + one accordion per step) --
        tags$div(
          tags$h4("Build a gene set", class = "fs-6 mb-2"),
          bslib::accordion(
            open = TRUE,
            bslib::accordion_panel(
              "1 · Select genes",
              .gs_pills(ns("source"),
                        c("Paste" = "paste", "From DE" = "deg", "Top-variable" = "topvar"),
                        selected = "paste"),
              conditionalPanel(
                cond("source", "paste"),
                gene_search_ui(ns, "paste", multiple = TRUE,
                               search_modes = c("exact", "contains", "regex")),
                uiOutput(ns("paste_literal_ui"))),
              conditionalPanel(cond("source", "deg"), uiOutput(ns("deg_controls"))),
              conditionalPanel(
                cond("source", "topvar"),
                numericInput(ns("topvar_n"), "Number of genes", value = 100, min = 1, step = 50),
                helpText(class = "small text-muted",
                         "The most variable endogenous features (VST)."))),
            bslib::accordion_panel(
              "2 · Preview",
              uiOutput(ns("preview_head")),
              DT::DTOutput(ns("preview_table")),
              tags$div(class = "mt-2",
                actionButton(ns("clear_staged"), "Clear", class = "btn-sm btn-outline-danger",
                             icon = shiny::icon("trash")))),
            bslib::accordion_panel(
              "3 · Save",
              uiOutput(ns("save_note")),
              .gs_pills(ns("save_mode"), c("New set" = "new", "Add to existing" = "add"),
                        selected = "new"),
              conditionalPanel(
                cond("save_mode", "new"),
                textInput(ns("new_name"), "Set name",
                          placeholder = "What should this set be called?"),
                uiOutput(ns("new_warn")),
                actionButton(ns("create"), "Create set", class = "btn btn-sm fw-semibold",
                             style = .gs_action_style, icon = shiny::icon("plus"))),
              conditionalPanel(
                cond("save_mode", "add"),
                uiOutput(ns("add_targets_ui")),
                actionButton(ns("add_existing"), "Add to selected", class = "btn btn-sm fw-semibold",
                             style = .gs_action_style, icon = shiny::icon("plus"))))
          )
        ),

        # ---- RIGHT: preview-control (small) + Your gene sets ---------------
        tags$div(
          class = "d-flex flex-column gap-3",
          bslib::accordion(
            open = TRUE,
            bslib::accordion_panel(
              "Gene set preview control",
              helpText(class = "small text-muted mb-2",
                       "Changes here affect all gene set preview tables on this page."),
              uiOutput(ns("show_cols_ui")))),
          tags$div(
            tags$div(
              class = "mb-2",
              tags$h4("Your gene sets", class = "fs-6 mb-0 d-inline"),
              tags$span(class = "small text-muted ms-2",
                        "Available to the other tabs this session")),
            bslib::accordion(
              open = TRUE,
              bslib::accordion_panel(
                "Defined sets",
                uiOutput(ns("empty_note")),
                DT::DTOutput(ns("sets_table")),
                tags$div(
                  class = "d-flex gap-2 my-2",
                  actionButton(ns("rename"), "Rename", class = "btn-sm btn-outline-secondary",
                               icon = shiny::icon("pen")),
                  actionButton(ns("delete"), "Delete", class = "btn-sm btn-outline-danger",
                               icon = shiny::icon("trash")),
                  actionButton(ns("clear_all"), "Clear all", class = "btn-sm btn-outline-danger",
                               icon = shiny::icon("xmark")))),
              bslib::accordion_panel(
                "Set members",
                uiOutput(ns("members_header")),
                DT::DTOutput(ns("members_table"))))
          )
        )
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

    # =====================================================================
    # 1 - Sources: each yields ids; staged() picks the ACTIVE source.
    # =====================================================================
    # Paste: literal-add is offered ONLY for an EXACT search by feature id
    # (rownames). Under another column `unmatched` holds names / class values
    # ("spike-in"), and in contains/regex mode it holds *patterns* ("^ENSMUSG.*")
    # -- neither are ids, and non-destructive membership would never trim them.
    .paste_literal_ok <- function(mode) {
      identical(input$paste_searchby %||% "__rownames__", "__rownames__") &&
        identical(mode %||% "exact", "exact")
    }
    output$paste_literal_ui <- renderUI({
      if (!.paste_literal_ok(input$paste_mode)) return(NULL)
      bslib::input_switch(ns("paste_literal"), "Also add unmatched as literal IDs", value = FALSE)
    })
    paste_ids <- reactive({
      r <- paste_search(); ids <- r$ids
      # Guard on the RESOLVED mode (authoritative), not just the input.
      if (isTRUE(input$paste_literal) && .paste_literal_ok(r$mode))
        ids <- unique(c(ids, r$unmatched))
      ids
    })

    # DE DEGs. Results are marked stale (not cleared) on a data edit, so surface
    # that here rather than silently seeding a set from a fit computed on a
    # different feature set (the project's never-reuse-a-stale-artifact rule).
    output$deg_controls <- renderUI({
      res <- (state$de %||% list())$results
      if (!length(res))
        return(helpText(class = "small text-muted",
          "Run DESeq2 (on the DE page) to seed a set from a contrast's DEGs."))
      tagList(
        if (identical(de_status(state), "stale"))
          tags$div(class = "small text-warning mb-2",
            "These DE results are out of date (the data changed since the fit). Re-run DESeq2 on the DE page before seeding a set."),
        selectInput(ns("deg_contrast"), "Contrast", choices = names(res)),
        radioButtons(ns("deg_dir"), "Direction",
                     c("Up" = "up", "Down" = "down", "Both" = "both"),
                     selected = "both", inline = TRUE),
        bslib::layout_columns(
          col_widths = c(6, 6),
          numericInput(ns("deg_padj"), "padj <", value = 0.05, min = 0, max = 1, step = 0.01),
          numericInput(ns("deg_lfc"), "|log2FC| >=", value = 1, min = 0, step = 0.5)),
        bslib::input_switch(ns("deg_shrunk"), "Use shrunk LFC", value = FALSE),
        uiOutput(ns("deg_preview")))
    })
    deg_classified <- reactive({
      res <- (state$de %||% list())$results; lab <- input$deg_contrast
      req(length(res), !is.null(lab), lab %in% names(res))
      de_classify_table(res[[lab]], input$deg_padj %||% 0.05, input$deg_lfc %||% log2(2))
    })
    # de_results() ALWAYS creates log2FoldChange_shrunk -- all-NA when shrinkage
    # was "none" or every fallback failed. A column-presence check would pass and
    # DEG_shrunk would be uniformly no_change (an empty set with no explanation),
    # so test for real shrunk values instead.
    deg_shrunk_ok <- reactive({
      df <- tryCatch(deg_classified(), error = function(e) NULL)
      !is.null(df) && "log2FoldChange_shrunk" %in% names(df) &&
        any(!is.na(df$log2FoldChange_shrunk))
    })
    deg_col <- reactive({
      if (isTRUE(input$deg_shrunk) && isTRUE(deg_shrunk_ok())) "DEG_shrunk" else "DEG"
    })
    deg_ids <- reactive({
      df <- deg_classified()
      want <- switch(input$deg_dir %||% "both", up = "up", down = "down", both = c("up", "down"))
      rownames(df)[as.character(df[[deg_col()]]) %in% want]
    })
    output$deg_preview <- renderUI({
      df <- tryCatch(deg_classified(), error = function(e) NULL); req(df)
      tab <- table(factor(as.character(df[[deg_col()]]), levels = c("up", "down", "no_change")))
      tagList(
        if (isTRUE(input$deg_shrunk) && !isTRUE(deg_shrunk_ok()))
          tags$div(class = "small text-warning",
            "Shrunken LFCs are not available for this fit -- using standard LFCs."),
        helpText(class = "small text-muted",
          sprintf("%d up / %d down at padj < %s, |log2FC| >= %s.",
                  tab[["up"]], tab[["down"]], format(input$deg_padj %||% 0.05),
                  format(input$deg_lfc %||% 1))))
    })

    # Top-variable: rank once per dataset (cached), slice top-n for the preview.
    topvar_ranked <- reactive({
      req(state$working)
      state_derive(state, "geneset_topvar_rank", params = list(), expr = function() {
        mat <- tryCatch(pca_input(state$working, "vst")$mat, error = function(e) NULL)
        if (is.null(mat)) return(character(0))
        top_variable_features(mat, n_top = nrow(mat))
      })
    })
    topvar_ids <- reactive({
      n <- suppressWarnings(as.integer(input$topvar_n %||% 100L)); if (is.na(n)) n <- 100L
      utils::head(topvar_ranked(), max(1L, n))
    })

    # The active source's provenance label + its staged sets (named list).
    staged_source <- function() switch(input$source %||% "paste",
      paste = "paste", deg = paste0("DE: ", input$deg_contrast %||% ""),
      topvar = "top-variable", "manual")
    # The parameters that DEFINE the staged set. Because a set is an authored
    # snapshot (its thresholds are local to this page and the user may move them
    # afterwards), these are the only record of how it was built -- log them so
    # the Phase 9 reproducibility export can regenerate the set.
    staged_params <- function() switch(input$source %||% "paste",
      paste = list(search_by = input$paste_searchby %||% "__rownames__",
                   match_mode = input$paste_mode %||% "exact",
                   case_insensitive = isTRUE(input$paste_ci),
                   literal_unmatched = isTRUE(input$paste_literal) &&
                     .paste_literal_ok(input$paste_mode)),
      deg = list(contrast = input$deg_contrast %||% "",
                 direction = input$deg_dir %||% "both",
                 padj = input$deg_padj %||% 0.05,
                 lfc = input$deg_lfc %||% 1,
                 shrunk = isTRUE(input$deg_shrunk) && isTRUE(deg_shrunk_ok()),
                 de_status = de_status(state)),
      topvar = list(n_top = input$topvar_n %||% 100L, input = "vst (endogenous)"),
      list())
    staged <- reactive({
      src <- input$source %||% "paste"
      ids <- switch(src,
        paste  = paste_ids(),
        deg    = tryCatch(deg_ids(), error = function(e) character(0)),
        topvar = topvar_ids(),
        character(0))
      ids <- unique(as.character(ids)); ids <- ids[!is.na(ids) & nzchar(ids)]
      if (!length(ids)) return(list())
      nm <- switch(src, paste = "Pasted genes", deg = staged_source(),
                   topvar = "Top-variable", "Staged")
      stats::setNames(list(ids), nm)
    })

    # =====================================================================
    # 2 - Preview of the staged genes.
    # =====================================================================
    output$preview_head <- renderUI({
      st <- staged(); n <- length(st)
      if (!n) return(helpText(class = "small text-muted", "Nothing staged yet."))
      total <- length(unique(unlist(st, use.names = FALSE)))
      txt <- if (n == 1L) sprintf("%d gene(s) staged.", total)
             else sprintf("%d sets staged (%d unique genes).", n, total)
      tagList(tags$div(class = "small text-muted mb-1", txt),
              if (n > 1L) selectizeInput(ns("preview_pick"), NULL, choices = names(st)))
    })
    staged_active_ids <- reactive({
      st <- staged(); if (!length(st)) return(character(0))
      if (length(st) == 1L) return(st[[1]])
      pick <- input$preview_pick %||% names(st)[1]
      st[[if (pick %in% names(st)) pick else names(st)[1]]]
    })
    output$preview_table <- DT::renderDT({
      ids <- staged_active_ids(); req(length(ids))
      dt_table(members_frame(ids), page_length = 5L)
    })
    observeEvent(input$clear_staged, {
      switch(input$source %||% "paste",
        paste  = updateTextAreaInput(session, "paste_q", value = ""),
        deg    = { updateRadioButtons(session, "deg_dir", selected = "both")
                   updateNumericInput(session, "deg_padj", value = 0.05)
                   updateNumericInput(session, "deg_lfc", value = 1) },
        topvar = updateNumericInput(session, "topvar_n", value = 100))
    })

    # =====================================================================
    # 3 - Save: New (reject-on-clash) / Add to existing.
    # =====================================================================
    output$save_note <- renderUI({
      if (length(staged()) > 1L)
        tags$div(class = "small text-warning mb-1",
          "Multiple sets staged -- bulk save (one per staged set) arrives in a later update.")
    })
    output$new_warn <- renderUI({
      nm <- trimws(input$new_name %||% "")
      if (nzchar(nm) && nm %in% names(state$gene_sets))
        tags$div(class = "small text-danger mb-1", "That name is taken -- pick another.")
    })
    output$add_targets_ui <- renderUI({
      if (!length(state$gene_sets))
        return(helpText(class = "small text-muted",
                        "Create a set first, then you can add to it."))
      selectizeInput(ns("add_targets"), "Add to set(s)", choices = names(state$gene_sets),
                     multiple = TRUE)
    })
    observeEvent(input$create, {
      st <- staged()
      if (length(st) != 1L) { showNotification("Stage exactly one set to create.", type = "warning"); return() }
      ids <- st[[1]]; nm <- trimws(input$new_name %||% "")
      if (!nzchar(nm)) { showNotification("Enter a set name.", type = "warning"); return() }
      if (nm %in% names(state$gene_sets)) { showNotification("That name is taken -- pick another.", type = "warning"); return() }
      if (!length(ids)) { showNotification("No genes staged.", type = "warning"); return() }
      res <- gene_set_commit(state$gene_sets, nm, ids, "new", source = staged_source())
      state$gene_sets <- res$sets; snapshot_absent()
      updateTextInput(session, "new_name", value = "")     # prevent an accidental re-add
      .log(state, list(action = "gene_set_add", set = res$name, source = staged_source(),
                       mode = "new", n = length(ids), params = staged_params()))
      showNotification(sprintf("Created '%s' (%d genes).", res$name, length(ids)),
                       type = "message", duration = 4)
    })
    observeEvent(input$add_existing, {
      st <- staged()
      if (length(st) != 1L) { showNotification("Stage exactly one set to add.", type = "warning"); return() }
      ids <- st[[1]]; targets <- input$add_targets
      if (!length(targets)) { showNotification("Select at least one set to add to.", type = "warning"); return() }
      if (!length(ids)) { showNotification("No genes staged.", type = "warning"); return() }
      sets <- state$gene_sets
      for (t in targets) sets <- gene_set_commit(sets, t, ids, "append", source = staged_source())$sets
      state$gene_sets <- sets; snapshot_absent()
      .log(state, list(action = "gene_set_addto", targets = targets, source = staged_source(),
                       n = length(ids), params = staged_params()))
      showNotification(sprintf("Added %d gene(s) to %d set(s).", length(ids), length(targets)),
                       type = "message", duration = 4)
    })

    # =====================================================================
    # Your gene sets: the shared member frame + tables.
    # =====================================================================
    # Extra rowData columns to show alongside id (default <feature_type>_name).
    output$show_cols_ui <- renderUI({
      req(state$working)
      cols <- names(as.data.frame(SummarizedExperiment::rowData(state$working), optional = TRUE))
      ft_name <- paste0(feature_type(), "_name")
      cur <- shiny::isolate(input$show_cols)
      sel <- if (!is.null(cur)) intersect(cur, cols)
             else if (ft_name %in% cols) ft_name else character(0)
      selectizeInput(ns("show_cols"), "Also show columns", choices = cols,
                     selected = sel, multiple = TRUE)
    })
    # id + chosen columns + in-dataset flag. Governs BOTH preview + members tables.
    members_frame <- function(ids) {
      df <- data.frame(ID = ids, check.names = FALSE, stringsAsFactors = FALSE)
      cols <- input$show_cols
      if (length(cols) && !is.null(state$working)) {
        rd <- SummarizedExperiment::rowData(state$working); hit <- ids %in% rownames(state$working)
        for (co in cols) {
          v <- rep(NA_character_, length(ids))
          if (co %in% colnames(rd)) v[hit] <- as.character(rd[ids[hit], co])
          df[[co]] <- v
        }
      }
      df[["In dataset"]] <- ifelse(ids %in% working_rn(), "yes", "no")
      df
    }

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
      helpText(class = "small text-muted", "No gene sets yet. Build one on the left.")
    })
    output$sets_table <- DT::renderDT({
      df <- sets_df(); req(df)
      dt_table(df, selection = "single", page_length = 10L)
    })
    # The selected set = the current table selection (single source of truth; a
    # deselect / delete reverts the members view to a placeholder).
    selected_name <- reactive({
      df <- sets_df(); i <- input$sets_table_rows_selected
      if (!is.null(df) && length(i) && i <= nrow(df)) df$Name[i] else NULL
    })
    output$members_header <- renderUI({
      nm <- selected_name()
      if (is.null(nm) || is.null(state$gene_sets[[nm]]))
        return(helpText(class = "small text-muted", "Select a set above to view its members."))
      tags$h5(sprintf("Members of '%s'", nm), class = "fs-6")
    })
    output$members_table <- DT::renderDT({
      nm <- selected_name(); req(!is.null(nm), !is.null(state$gene_sets[[nm]]))
      dt_table(members_frame(state$gene_sets[[nm]]$ids), page_length = 10L)
    })

    # --- Rename / Delete / Clear -------------------------------------------
    observeEvent(input$rename, {
      nm <- selected_name(); req(!is.null(nm))
      showModal(modalDialog(
        title = "Rename gene set", size = "s",
        textInput(ns("rename_to"), "New name", value = nm),
        footer = tagList(modalButton("Cancel"),
                         actionButton(ns("rename_ok"), "Rename", class = "btn-primary"))))
    })
    observeEvent(input$rename_ok, {
      nm <- selected_name(); to <- trimws(input$rename_to %||% "")
      req(!is.null(nm), !is.null(state$gene_sets[[nm]]))
      if (!nzchar(to)) { showNotification("Enter a name.", type = "warning"); return() }
      if (!identical(to, nm) && !is.null(state$gene_sets[[to]])) {
        showNotification("A set with that name already exists.", type = "warning"); return() }
      if (!identical(to, nm)) {
        sets <- state$gene_sets; names(sets)[names(sets) == nm] <- to; state$gene_sets <- sets
        snapshot_absent()   # re-key the absence baseline, else the rename reads as "newly absent"
        .log(state, list(action = "gene_set_rename", from = nm, to = to))
      }
      removeModal()
    })
    observeEvent(input$delete, {
      nm <- selected_name(); req(!is.null(nm), !is.null(state$gene_sets[[nm]]))
      sets <- state$gene_sets; sets[[nm]] <- NULL; state$gene_sets <- sets
      snapshot_absent()
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
      state$gene_sets <- list(); prev_absent(list())
      .log(state, list(action = "gene_set_clear_all"))
      removeModal()
    })

    # --- Non-destructive reconcile on a data edit: notify on NEW absences ---
    prev_absent <- reactiveVal(list())
    snapshot_absent <- function() {
      rn <- working_rn()
      prev_absent(lapply(state$gene_sets, function(s) gene_set_absent(s, rn)))
    }
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
