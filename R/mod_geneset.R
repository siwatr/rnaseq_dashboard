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

# A namespaced conditionalPanel condition. Module-level (NOT a UI-local closure)
# so the server's renderUI blocks -- e.g. the multi-set Save UI -- can build
# conditions too; `ns` is the UI's NS() in the UI, session$ns in the server.
.gs_cond <- function(ns, inp, val) sprintf("input['%s'] == '%s'", ns(inp), val)

# Which field to match an imported ID column against: the rowData column (or
# rownames) with the most hits for `query`. Ties: rownames win, then first
# rowData column. Returns a feature_search_choices() value ("__rownames__" or a
# column name). Cheap (hashed match over a handful of columns).
.gs_best_match_field <- function(query, dds) {
  q <- unique(trimws(as.character(query))); q <- q[!is.na(q) & nzchar(q)]
  if (!length(q)) return("__rownames__")
  best <- "__rownames__"; best_n <- sum(q %in% rownames(dds))
  rd <- as.data.frame(SummarizedExperiment::rowData(dds), optional = TRUE)
  ok <- vapply(rd, function(x) is.character(x) || is.factor(x) || is.numeric(x), logical(1))
  for (co in names(rd)[ok]) {
    n <- sum(q %in% as.character(rd[[co]]))
    if (n > best_n) { best_n <- n; best <- co }   # strict > keeps the earlier winner on a tie
  }
  best
}

mod_geneset_ui <- function(id) {
  ns <- NS(id)
  cond <- function(inp, val) .gs_cond(ns, inp, val)
  bslib::navset_card_tab(
    id = ns("tabs"),
    bslib::nav_panel(
      tags$h3("Manage sets", class = "fs-6 mb-0"),
      bslib::layout_columns(
        fillable = FALSE, col_widths = c(6, 6),

        # ---- LEFT: Build a gene set (title text + one accordion per step) --
        tags$div(
          tags$h4("Build a gene set", class = "fs-5 mb-2"),
          bslib::accordion(
            open = TRUE,
            bslib::accordion_panel(
              "1 · Select genes",
              .gs_pills(ns("source"),
                        c("Paste" = "paste", "From DE" = "deg", "Top-variable" = "topvar",
                          "Import table" = "file", "Gene set file" = "gsfile"),
                        selected = "paste"),
              conditionalPanel(
                cond("source", "paste"),
                gene_search_ui(ns, "paste", multiple = TRUE,
                               search_modes = c("exact", "contains", "regex"))),
              conditionalPanel(cond("source", "deg"), uiOutput(ns("deg_controls"))),
              conditionalPanel(
                cond("source", "topvar"),
                numericInput(ns("topvar_n"), "Number of genes", value = 100, min = 1, step = 50),
                helpText(class = "small text-muted",
                         "The most variable endogenous features (VST).")),
              conditionalPanel(
                cond("source", "file"),
                uiOutput(ns("tbl_file_ui")),   # rebuilt on Clear so re-uploading the same file works
                bslib::input_switch(ns("tbl_header"), "File has a header row", value = TRUE),
                actionButton(ns("tbl_load"), "Load file", class = "btn-sm",
                             icon = shiny::icon("file-import"), style = .gs_action_style),
                uiOutput(ns("tbl_controls")),
                uiOutput(ns("tbl_loaded_head")),
                DT::DTOutput(ns("tbl_table")),
                uiOutput(ns("tbl_stats"))),
              conditionalPanel(
                cond("source", "gsfile"),
                uiOutput(ns("gsfile_file_ui")),   # rebuilt on Clear so re-uploading works
                actionButton(ns("gsfile_load"), "Load file", class = "btn-sm",
                             icon = shiny::icon("file-import"), style = .gs_action_style),
                helpText(class = "small text-muted",
                         "A previously exported JSON, GMT, or long TSV (set / id) file. ",
                         "Format is detected from the extension. The file's named sets are ",
                         "staged below -- review, then save."),
                uiOutput(ns("gsfile_controls")),
                uiOutput(ns("gsfile_summary")))),
            bslib::accordion_panel(
              "2 · Preview",
              uiOutput(ns("preview_head")),
              # The keep-unmatched toggle lives here (not with the source) because
              # it directly changes what enters the set, next to the drop warning.
              uiOutput(ns("preview_literal_ui")),
              DT::DTOutput(ns("preview_table")),
              tags$div(class = "mt-2",
                actionButton(ns("clear_staged"), "Clear", class = "btn-sm btn-outline-danger",
                             icon = shiny::icon("trash")))),
            bslib::accordion_panel(
              "3 · Save",
              uiOutput(ns("save_note")),
              # One staged set -> the New/Add controls (kept in conditionalPanels,
              # never renderUI, so typing in the name box is never torn down).
              # Many staged sets (an annotation-split import) -> the multi-set UI.
              # `multi_staged` is a server flag exposed to the client condition.
              conditionalPanel(
                sprintf("output['%s'] != '1'", ns("multi_staged")),
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
                  actionButton(ns("add_existing"), "Add to selected",
                               class = "btn btn-sm fw-semibold",
                               style = .gs_action_style, icon = shiny::icon("plus")))),
              conditionalPanel(sprintf("output['%s'] == '1'", ns("multi_staged")),
                               uiOutput(ns("multi_save"))))
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
              tags$h4("Your gene sets", class = "fs-5 mb-0 d-inline"),
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
                DT::DTOutput(ns("members_table"))))),
          # ---- Export: its own flex child (own section, gets the gap-3), with a
          # top divider so it reads as a separate zone from "Your gene sets" -----
          tags$div(
            class = "pt-3 border-top",
            tags$div(
              class = "mb-2",
              tags$h4("Export gene sets", class = "fs-5 mb-0 d-inline"),
              tags$span(class = "small text-muted ms-2",
                        "Download sets as JSON / GMT / TSV")),
            bslib::accordion(
              open = TRUE,
              bslib::accordion_panel("Export format", uiOutput(ns("export_ui"))),
              bslib::accordion_panel("Preview", uiOutput(ns("export_preview_ui")))))
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
    # The active source's "keep unmatched" toggle, rendered in the Preview
    # accordion (next to the drop warning). It keeps each source's own input id
    # (paste_literal / tbl_literal / gsfile_literal), so all resolution logic is
    # unchanged -- only WHERE the switch is drawn moved. isolate() preserves the
    # value across re-renders (a match-field / mode change), without self-looping.
    output$preview_literal_ui <- renderUI({
      mk <- function(id, label)
        bslib::input_switch(ns(id), label, value = isTRUE(shiny::isolate(input[[id]])))
      sw <- switch(input$source %||% "paste",
        paste = if (.paste_literal_ok(input$paste_mode))
          mk("paste_literal", "Also add unmatched as literal IDs"),
        file = if (identical(input$tbl_match_by %||% "__rownames__", "__rownames__"))
          mk("tbl_literal", "Also keep unmatched as literal IDs"),
        gsfile = mk("gsfile_literal",
          if (identical(input$gsfile_match_by %||% "__rownames__", "__rownames__"))
            "Keep IDs not in the dataset" else "Keep unmatched IDs (not in the dataset)"),
        NULL)
      if (is.null(sw)) NULL else tags$div(class = "mb-2", sw)
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

    # ---- Table import: an arbitrary sheet (e.g. a DESeq2 result table from
    # another analysis) -> view/filter/select -> one or MANY staged sets. -----
    # The fileInput is rebuilt under a nonce so Clear truly resets it -- otherwise
    # re-selecting the SAME file emits no browser change event and the importer
    # stays stuck (the classic fileInput footgun; there is no updateFileInput).
    tbl_nonce <- reactiveVal(0L)
    output$tbl_file_ui <- renderUI({
      tbl_nonce()
      fileInput(ns("tbl_file"), "Table file (CSV / TSV / XLSX)",
                accept = c(".csv", ".tsv", ".txt", ".xlsx", ".xls"))
    })
    # Reading is intentional: the file is read only on the Load button (like the
    # app's other uploads). `tbl_src` holds the loaded file REFERENCE, not the
    # parsed frame -- so the header toggle can re-parse reactively after the first
    # load without a re-click. Shiny keeps the upload on disk for the session, so
    # this is a cheap re-read with no in-memory duplicate. NULL src -> NULL raw, so
    # staged() always yields a list (no req cascade escapes the Save observers).
    tbl_src <- reactiveVal(NULL)
    observeEvent(input$tbl_load, {
      f <- input$tbl_file
      if (is.null(f)) { showNotification("Choose a file first.", type = "warning"); return() }
      tbl_src(list(datapath = f$datapath, name = f$name))
    })
    tbl_raw <- reactive({
      src <- tbl_src(); if (is.null(src)) return(NULL)
      tryCatch(.read_user_table(src$datapath, src$name, col_names = isTRUE(input$tbl_header)),
               error = function(e) { showNotification(conditionMessage(e), type = "error"); NULL })
    })

    output$tbl_controls <- renderUI({
      df <- tbl_raw(); req(df, ncol(df) > 0, state$working)
      cols <- names(df)
      tip <- function(x, msg) bslib::tooltip(x, msg, placement = "right")
      tagList(
        bslib::layout_columns(
          col_widths = c(6, 6),
          tip(selectInput(ns("tbl_id_col"), "ID column (in your file)",
                          choices = cols, selected = cols[1]),
              "The column of the uploaded file that holds the gene / feature identifiers."),
          tip(selectInput(ns("tbl_match_by"), "Match against (dataset)",
                          choices = feature_search_choices(
                            SummarizedExperiment::rowData(state$working)),
                          selected = "__rownames__"),
              "Which part of the loaded dataset to match those identifiers to -- the feature IDs (row names) or a rowData column such as gene_name. Auto-detected from your ID column; override if needed.")),
        # A single-column file has nothing to split by.
        if (length(cols) > 1L)
          tip(selectizeInput(ns("tbl_anno_cols"), "Split into sets by column(s) (optional)",
                             choices = cols, multiple = TRUE),
              "Create one gene set per unique value (or value combination) in these columns -- e.g. split a results table by a 'direction' column into up / down sets."),
        conditionalPanel(
          sprintf("input['%s'] && input['%s'].length > 1", ns("tbl_anno_cols"), ns("tbl_anno_cols")),
          tip(textInput(ns("tbl_sep"), "Group name separator", value = "."),
              "Joins the values of multiple split columns into each set's name (e.g. 'up.brain').")),
        radioButtons(ns("tbl_rows"), "Use rows",
                     c("All rows in the current view (filtered)" = "view",
                       "Selected rows only" = "selected"),
                     selected = "view"),
        uiOutput(ns("tbl_match_opts")))
    })
    # Auto-detect the best "Match against" field for the chosen ID column.
    # (Re-detects on a header flip indirectly: toggling the header re-renders
    # tbl_controls and changes the ID-column name, which re-fires input$tbl_id_col.)
    observeEvent(list(input$tbl_id_col, tbl_src()), {
      df <- tbl_raw(); idc <- input$tbl_id_col
      req(df, !is.null(idc), idc %in% names(df), state$working)
      updateSelectInput(session, "tbl_match_by",
                        selected = .gs_best_match_field(df[[idc]], state$working))
    }, ignoreInit = TRUE)
    # Keep-all vs first-only for a name column that can match several features.
    # (The literal-add toggle moved to the Preview accordion.)
    output$tbl_match_opts <- renderUI({
      if (identical(input$tbl_match_by %||% "__rownames__", "__rownames__"))
        return(NULL)
      bslib::tooltip(
        radioButtons(ns("tbl_multi"), "When a name matches several features",
                     c("Keep all matches" = "all", "First match only" = "first"),
                     selected = "all", inline = TRUE),
        "A gene name can map to more than one feature (duplicated / paralogous IDs). Keep all adds every matching feature; First only keeps one.",
        placement = "top")
    })
    output$tbl_loaded_head <- renderUI({
      if (is.null(tbl_raw())) return(NULL)
      tags$div(class = "fw-semibold small mt-2 mb-1", "Loaded data")
    })
    output$tbl_table <- DT::renderDT({
      df <- tbl_raw(); req(df)
      dt_table(df, selection = list(mode = "multiple"))
    })
    # The rows the user is acting on: the filtered view (_rows_all, across pages)
    # or the explicit selection (_rows_selected).
    tbl_rows <- reactive({
      df <- tbl_raw(); if (is.null(df)) return(integer(0))
      if (identical(input$tbl_rows %||% "view", "selected")) input$tbl_table_rows_selected
      else input$tbl_table_rows_all %||% seq_len(nrow(df))
    })
    # Resolve the in-scope rows' ID column to feature id(s). 1-to-MANY when
    # matching a name column with Keep-all (a name can hit several features).
    # Returns per-input-row id vectors + the rows, so the expanded staging frame
    # and the stats both build from it.
    # value -> all dataset row indices for the chosen match column. Its own
    # reactive (keyed on the column + dataset) so it isn't rebuilt over all
    # features on every import-table filter keystroke. NULL for rownames matching.
    tbl_match_index <- reactive({
      field <- input$tbl_match_by %||% "__rownames__"
      if (identical(field, "__rownames__") || is.null(state$working)) return(NULL)
      vals <- as.character(as.data.frame(
        SummarizedExperiment::rowData(state$working), optional = TRUE)[[field]])
      split(seq_along(vals), vals)                        # NA-valued features drop out
    })
    tbl_match <- reactive({
      df <- tbl_raw(); if (is.null(df) || is.null(state$working)) return(NULL)
      idc <- input$tbl_id_col
      if (is.null(idc) || !(idc %in% names(df))) return(NULL)
      rows <- tbl_rows(); if (!length(rows)) return(NULL)
      sub <- df[rows, , drop = FALSE]
      raw <- trimws(as.character(sub[[idc]])); raw[!nzchar(raw)] <- NA   # blanks are not ids
      field <- input$tbl_match_by %||% "__rownames__"
      rn <- rownames(state$working)
      if (identical(field, "__rownames__")) {
        lit <- isTRUE(input$tbl_literal)                 # rownames are unique -> 1:1
        ids <- lapply(raw, function(x)
          if (is.na(x)) character(0) else if (x %in% rn) x else if (lit) x else character(0))
      } else {
        by_val <- tbl_match_index()                       # NULL[[x]] -> NULL, so still safe
        keep_all <- !identical(input$tbl_multi %||% "all", "first")
        ids <- lapply(raw, function(x) {
          if (is.na(x)) return(character(0))
          hit <- rn[by_val[[x]] %||% integer(0)]
          if (!length(hit)) character(0) else if (keep_all) hit else hit[1]
        })
      }
      list(sub = sub, ids = ids)
    })
    # The expanded frame (one row per input-row x matched-id) for the split.
    tbl_resolved <- reactive({
      m <- tbl_match(); if (is.null(m)) return(NULL)
      reps <- pmax(lengths(m$ids), 1L)                    # 0-match rows kept as one NA row
      exp <- m$sub[rep(seq_len(nrow(m$sub)), reps), , drop = FALSE]
      exp$.gs_id <- unlist(lapply(m$ids, function(z) if (length(z)) z else NA_character_))
      exp
    })
    tbl_staged <- reactive({
      sub <- tbl_resolved(); if (is.null(sub)) return(list())
      anno <- intersect(input$tbl_anno_cols %||% character(0), names(sub))
      sep  <- input$tbl_sep %||% "."; if (!nzchar(sep)) sep <- "."
      out <- split_ids_by_group(sub, ".gs_id", anno, sep = sep)
      if (!length(anno) && length(out)) names(out) <- "Imported genes"
      out
    })
    output$tbl_stats <- renderUI({
      df <- tbl_raw(); req(df)
      m <- tbl_match()
      n_view <- length(input$tbl_table_rows_all %||% seq_len(nrow(df)))
      n_sel  <- length(input$tbl_table_rows_selected)
      st <- tbl_staged()
      n_ids   <- if (is.null(m)) 0L else length(unique(unlist(m$ids)))
      n_miss  <- if (is.null(m)) 0L else sum(lengths(m$ids) == 0L)
      n_multi <- if (is.null(m)) 0L else sum(lengths(m$ids) > 1L)
      grp <- if (length(st) > 1L) {
        shown <- utils::head(st, 5L)
        txt <- paste(sprintf("%s (%d)", names(shown), lengths(shown)), collapse = ", ")
        if (length(st) > 5L) txt <- paste0(txt, sprintf(", +%d more", length(st) - 5L))
        tags$div(class = "small text-muted", sprintf("%d sets: %s", length(st), txt))
      }
      multi_note <- if (n_multi > 0L)
        tags$div(class = "small text-warning",
          sprintf("%d name(s) matched multiple features -- %s.", n_multi,
                  if (identical(input$tbl_multi %||% "all", "first")) "kept the first"
                  else "all their IDs added"))
      tagList(
        # "resolved" (not "matched") because with the literal switch on some ids
        # are kept verbatim without being in the dataset.
        tags$div(class = "small text-muted",
          sprintf("%d rows imported | %d in view | %d selected | %d IDs resolved%s.",
                  nrow(df), n_view, n_sel, n_ids,
                  if (n_miss) sprintf(" (%d unmatched)", n_miss) else "")),
        multi_note, grp)
    })

    # ---- Gene set file import (P6d): a serialized JSON / GMT / TSV file ------
    # A file already carries complete NAMED sets, so it stages them directly (the
    # Preview shows them) and flows into the multi-set Save UI (New / Add) -- the
    # same non-destructive path the annotation-split import uses. Like tbl, the
    # fileInput is rebuilt under a nonce so re-selecting the same file re-fires.
    gsfile_nonce <- reactiveVal(0L)
    output$gsfile_file_ui <- renderUI({
      gsfile_nonce()
      fileInput(ns("gsfile_file"), "Gene set file (JSON / GMT / TSV)",
                accept = c(".json", ".gmt", ".tsv", ".txt", ".csv"))
    })
    # Parsed store: list(recs = named gene-set records, name = source file name).
    gsfile_parsed <- reactiveVal(NULL)
    observeEvent(input$gsfile_load, {
      f <- input$gsfile_file
      if (is.null(f)) { showNotification("Choose a file first.", type = "warning"); return() }
      recs <- tryCatch(gene_sets_from_file(f$datapath, name = f$name),
                       error = function(e) { showNotification(conditionMessage(e), type = "error"); NULL })
      if (is.null(recs)) return()
      if (!length(recs)) {
        showNotification("No gene sets found in that file.", type = "warning")
        gsfile_parsed(NULL); return()
      }
      gsfile_parsed(list(recs = recs, name = f$name))
    })
    # A file may use a DIFFERENT id scheme than this dataset's rownames (e.g. a
    # gene_name column from another pipeline), so the same match-field scheme as
    # the table import applies: pick what the file's IDs match against, keep-all
    # vs first when a name hits several features, and a "keep IDs not in the
    # dataset" escape hatch. The match field is auto-detected on load.
    output$gsfile_controls <- renderUI({
      p <- gsfile_parsed(); req(p, state$working)
      q <- unique(unlist(lapply(p$recs, `[[`, "ids"), use.names = FALSE))
      det <- .gs_best_match_field(q, state$working)
      tagList(
        bslib::tooltip(
          selectInput(ns("gsfile_match_by"), "Match against (dataset)",
                      choices = feature_search_choices(
                        SummarizedExperiment::rowData(state$working)),
                      selected = det),
          "Which part of the dataset the file's IDs match -- the feature IDs (row names) or a rowData column such as gene_name. Auto-detected from the file; override if needed.",
          placement = "right"),
        uiOutput(ns("gsfile_match_opts")))
    })
    # Keep-all vs first-only for a name column. (The keep-unmatched toggle moved
    # to the Preview accordion.)
    output$gsfile_match_opts <- renderUI({
      if (identical(input$gsfile_match_by %||% "__rownames__", "__rownames__"))
        return(NULL)
      bslib::tooltip(
        radioButtons(ns("gsfile_multi"), "When a name matches several features",
                     c("Keep all matches" = "all", "First match only" = "first"),
                     selected = "all", inline = TRUE),
        "A gene name can map to several features. Keep all adds every match; First only keeps one.",
        placement = "top")
    })
    # dataset-value -> row indices for the chosen match column (NULL for rownames).
    # trimws both sides (file ids are trimmed at resolve time) so a rowData value
    # with stray whitespace still matches a clean file id.
    gsfile_match_index <- reactive({
      field <- input$gsfile_match_by %||% "__rownames__"
      if (identical(field, "__rownames__") || is.null(state$working)) return(NULL)
      vals <- trimws(as.character(as.data.frame(
        SummarizedExperiment::rowData(state$working), optional = TRUE)[[field]]))
      split(seq_along(vals), vals)                      # NA-valued features drop out
    })
    # Resolve each set's file IDs -> dataset ids (1:many keep-all, or first), with
    # the literal escape hatch. Returns per-set resolved ids + the pooled misses.
    gsfile_resolved <- reactive({
      p <- gsfile_parsed()
      if (is.null(p) || !length(p$recs) || is.null(state$working)) return(NULL)
      field <- input$gsfile_match_by %||% "__rownames__"
      rn <- rownames(state$working)
      keep_all <- !identical(input$gsfile_multi %||% "all", "first")
      literal  <- isTRUE(input$gsfile_literal)
      idx <- gsfile_match_index()
      resolve_one <- function(raw) {                    # one file id -> dataset id(s)
        if (identical(field, "__rownames__"))
          return(if (raw %in% rn) raw else character(0))
        h <- rn[idx[[raw]] %||% integer(0)]
        if (!length(h) || keep_all) h else h[1]
      }
      per <- lapply(p$recs, function(r) {
        raw  <- unique(trimws(as.character(r$ids))); raw <- raw[nzchar(raw)]
        hits <- lapply(raw, resolve_one)
        miss <- raw[lengths(hits) == 0L]
        ids  <- unique(unlist(hits, use.names = FALSE))
        if (literal) ids <- unique(c(ids, miss))
        list(ids = ids, miss = miss)
      })
      list(sets = lapply(per, `[[`, "ids"),
           unmatched = unique(unlist(lapply(per, `[[`, "miss"), use.names = FALSE)))
    })
    # Resolved id-vectors only (the staged() contract). A file's kind/annotation
    # aren't carried and its source is relabelled "import: <file>" on commit --
    # lossless for P6's simple sets; TODO(P7): an annotated import would downgrade
    # to simple here, so the annotated layer must resolve records, not just ids.
    gsfile_staged <- reactive({
      r <- gsfile_resolved(); if (is.null(r)) return(list())
      r$sets[lengths(r$sets) > 0L]                      # drop sets that resolved empty
    })
    output$gsfile_summary <- renderUI({
      p <- gsfile_parsed(); if (is.null(p)) return(NULL)
      r <- gsfile_resolved()
      tot   <- length(unique(unlist(lapply(p$recs, `[[`, "ids"), use.names = FALSE)))
      n_res <- if (is.null(r)) 0L else length(unique(unlist(r$sets, use.names = FALSE)))
      n_mis <- if (is.null(r)) 0L else length(r$unmatched)
      tags$div(class = "small text-muted mt-2",
               sprintf("Loaded %d set(s), %d unique gene(s) from '%s'. %d ID(s) resolved%s.",
                       length(p$recs), tot, p$name, n_res,
                       if (n_mis) sprintf(" (%d unmatched)", n_mis) else ""))
    })

    # The active source's provenance label + its staged sets (named list).
    staged_source <- function() switch(input$source %||% "paste",
      paste = "paste", deg = paste0("DE: ", input$deg_contrast %||% ""),
      topvar = "top-variable",
      file = paste0("import: ", (tbl_src() %||% list(name = "table"))$name),
      gsfile = paste0("import: ", (gsfile_parsed() %||% list(name = "file"))$name),
      "manual")
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
      file = list(file = (tbl_src() %||% list(name = ""))$name,
                  multi_match = if (identical(input$tbl_match_by %||% "__rownames__", "__rownames__"))
                    NA_character_ else (input$tbl_multi %||% "all"),
                  id_col = input$tbl_id_col %||% "",
                  match_by = input$tbl_match_by %||% "__rownames__",
                  rows = input$tbl_rows %||% "view",
                  split_by = input$tbl_anno_cols %||% character(0),
                  sep = input$tbl_sep %||% ".",
                  literal_unmatched = isTRUE(input$tbl_literal) &&
                    identical(input$tbl_match_by %||% "__rownames__", "__rownames__")),
      gsfile = list(file = (gsfile_parsed() %||% list(name = ""))$name,
                    match_by = input$gsfile_match_by %||% "__rownames__",
                    multi = if (identical(input$gsfile_match_by %||% "__rownames__", "__rownames__"))
                      NA_character_ else (input$gsfile_multi %||% "all"),
                    literal_unmatched = isTRUE(input$gsfile_literal)),
      list())
    # The staged sets: a NAMED LIST, so a source may stage one set (paste / DE /
    # top-variable) or many (an annotation-split import).
    staged <- reactive({
      src <- input$source %||% "paste"
      if (identical(src, "file")) return(tbl_staged())
      if (identical(src, "gsfile")) return(gsfile_staged())
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
    # How many queried ids the active source DROPPED as not-in-dataset (only when
    # the literal escape hatch is off), + whether that hatch is available here.
    dropped_info <- function() {
      src <- input$source %||% "paste"
      if (identical(src, "paste")) {
        can <- .paste_literal_ok(input$paste_mode)
        if (isTRUE(input$paste_literal) && can) return(list(n = 0L, can = can))
        return(list(n = length(paste_search()$unmatched), can = can))
      }
      if (identical(src, "file")) {
        can <- identical(input$tbl_match_by %||% "__rownames__", "__rownames__")
        sub <- tbl_resolved(); if (is.null(sub)) return(list(n = 0L, can = can))
        if (isTRUE(input$tbl_literal) && can) return(list(n = 0L, can = can))
        return(list(n = sum(is.na(sub$.gs_id)), can = can))
      }
      if (identical(src, "gsfile")) {
        if (isTRUE(input$gsfile_literal)) return(list(n = 0L, can = TRUE))
        r <- gsfile_resolved()
        return(list(n = if (is.null(r)) 0L else length(r$unmatched), can = TRUE))
      }
      list(n = 0L, can = FALSE)
    }
    output$preview_head <- renderUI({
      st <- staged(); n <- length(st)
      if (!n && !dropped_info()$n)
        return(helpText(class = "small text-muted", "Nothing staged yet."))
      total <- length(unique(unlist(st, use.names = FALSE)))
      txt <- if (n <= 1L) sprintf("%d gene(s) staged.", total)
             else sprintf("%d sets staged (%d unique genes).", n, total)
      d <- dropped_info()
      warn <- if (d$n > 0L)
        tags$div(class = "small text-warning mb-1",
          sprintf("%d ID(s) not in the dataset were left out.%s", d$n,
                  if (d$can) " Turn on the 'keep unmatched IDs' option to include them." else ""))
      tagList(warn,
              tags$div(class = "small text-muted mb-1", txt),
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
      dt_table(members_frame(ids))
    })
    observeEvent(input$clear_staged, {
      switch(input$source %||% "paste",
        paste  = updateTextAreaInput(session, "paste_q", value = ""),
        deg    = { updateRadioButtons(session, "deg_dir", selected = "both")
                   updateNumericInput(session, "deg_padj", value = 0.05)
                   updateNumericInput(session, "deg_lfc", value = 1) },
        topvar = updateNumericInput(session, "topvar_n", value = 100),
        file   = { tbl_src(NULL)                                 # drop the loaded table
                   tbl_nonce(shiny::isolate(tbl_nonce()) + 1L)   # rebuild -> resets the file input
                   updateSelectizeInput(session, "tbl_anno_cols", selected = character(0)) },
        gsfile = { gsfile_parsed(NULL)                           # drop the loaded sets
                   gsfile_nonce(shiny::isolate(gsfile_nonce()) + 1L) })   # reset the file input
    })

    # =====================================================================
    # 3 - Save: New (reject-on-clash) / Add to existing.
    # =====================================================================
    # A client-visible flag so the Save section can swap between the single-set
    # and multi-set UIs without a renderUI (which would reset the name box).
    # An explicit renderText (not a bare reactive) so the value is a plain "1"/"0"
    # for the conditionalPanel JS; suspendWhenHidden = FALSE keeps it live while
    # the panel it drives is hidden.
    # A gene-set file already carries named sets, so it always uses the multi-set
    # Save UI (which commits under those names) even when it holds a single set.
    output$multi_staged <- renderText(
      if (length(staged()) > 1L || identical(input$source %||% "paste", "gsfile")) "1" else "0")
    outputOptions(output, "multi_staged", suspendWhenHidden = FALSE)

    output$save_note <- renderUI({
      n <- length(staged())
      if (n > 1L)
        tags$div(class = "small text-muted mb-1",
                 sprintf("%d sets staged -- save them in one go below.", n))
    })

    # ---- Multi-set save (staged > 1, i.e. an annotation-split import) -------
    # Value-stable staged names, so re-filtering the imported table doesn't reset
    # the picker unless the group names actually change.
    staged_names_rv <- reactiveVal(character(0))
    observe({
      nms <- names(staged())
      if (!identical(nms, shiny::isolate(staged_names_rv()))) staged_names_rv(nms)
    })
    output$multi_save <- renderUI({
      # >1 staged set, or a gene-set-file import of any size (names come from it).
      nms <- staged_names_rv()
      if (!length(nms) ||
          (length(nms) == 1L && !identical(input$source %||% "paste", "gsfile")))
        return(NULL)
      pick_id <- ns("multi_pick")
      tagList(
        # Cap the selected-items box height so many staged sets don't grow the UI
        # unbounded (the dropdown is absolutely positioned, so it isn't clipped).
        tags$style(HTML(sprintf(
          "#%s + .selectize-control .selectize-input { max-height: 110px; overflow-y: auto; }",
          pick_id))),
        selectizeInput(pick_id, "Staged sets to save", choices = nms,
                       selected = nms, multiple = TRUE),
        tags$div(class = "d-flex gap-2 mb-2",
          actionButton(ns("multi_all"), "Select all", class = "btn-sm btn-outline-secondary"),
          actionButton(ns("multi_none"), "Deselect all", class = "btn-sm btn-outline-secondary")),
        .gs_pills(ns("multi_mode"), c("Create new sets" = "new", "Add to existing" = "add"),
                  selected = "new"),
        conditionalPanel(
          .gs_cond(ns, "multi_mode", "new"),
          textInput(ns("multi_prefix"), "Name prefix (optional)", placeholder = "e.g. res_"),
          uiOutput(ns("multi_conflicts")),
          bslib::input_switch(ns("multi_autoname"), "Auto-rename conflicts", value = FALSE),
          actionButton(ns("multi_create"), "Create sets", class = "btn btn-sm fw-semibold",
                       style = .gs_action_style, icon = shiny::icon("plus"))),
        conditionalPanel(
          .gs_cond(ns, "multi_mode", "add"),
          uiOutput(ns("multi_targets_ui")),
          actionButton(ns("multi_add"), "Add to selected", class = "btn btn-sm fw-semibold",
                       style = .gs_action_style, icon = shiny::icon("plus"))))
    })
    multi_pick <- reactive(intersect(input$multi_pick %||% character(0), names(staged())))
    observeEvent(input$multi_all,
      updateSelectizeInput(session, "multi_pick", selected = names(staged())))
    observeEvent(input$multi_none,
      updateSelectizeInput(session, "multi_pick", selected = character(0)))
    multi_final_names <- reactive(paste0(trimws(input$multi_prefix %||% ""), multi_pick()))
    output$multi_conflicts <- renderUI({
      clash <- intersect(multi_final_names(), names(state$gene_sets))
      if (!length(clash)) return(NULL)
      if (isTRUE(input$multi_autoname))
        return(tags$div(class = "small text-muted mb-1",
          sprintf("%d name(s) taken -- they will be auto-renamed (e.g. %s_2).",
                  length(clash), clash[1])))
      tags$div(class = "small text-danger mb-1",
        sprintf("Already taken: %s. Add a prefix or turn on auto-rename.",
                paste(utils::head(clash, 5), collapse = ", ")))
    })
    output$multi_targets_ui <- renderUI({
      if (!length(state$gene_sets))
        return(helpText(class = "small text-muted",
                        "Create a set first, then you can add to it."))
      selectizeInput(ns("multi_targets"), "Add to set(s)", choices = names(state$gene_sets),
                     multiple = TRUE)
    })
    observeEvent(input$multi_create, {
      st <- staged(); pick <- multi_pick()
      if (!length(pick)) { showNotification("Select at least one staged set.", type = "warning"); return() }
      fin <- multi_final_names()
      clash <- intersect(fin, names(state$gene_sets))
      if (length(clash) && !isTRUE(input$multi_autoname)) {
        showNotification(sprintf("These names are taken: %s. Add a prefix or turn on auto-rename.",
                                 paste(utils::head(clash, 5), collapse = ", ")),
                         type = "warning"); return()
      }
      sets <- state$gene_sets; made <- character(0)
      for (i in seq_along(pick)) {                 # New auto-suffixes any clash
        r <- gene_set_commit(sets, fin[i], st[[pick[i]]], "new", source = staged_source())
        sets <- r$sets; made <- c(made, r$name)
      }
      state$gene_sets <- sets; snapshot_absent()
      # Reset the prefix + deselect the just-created sets so the conflict banner
      # doesn't flip to "already taken" right after a successful create.
      updateTextInput(session, "multi_prefix", value = "")
      updateSelectizeInput(session, "multi_pick", selected = character(0))
      .log(state, list(action = "gene_set_add_multi", sets = made, source = staged_source(),
                       n = unname(vapply(st[pick], length, integer(1))),
                       params = staged_params()))
      showNotification(sprintf("Created %d set(s): %s.", length(made),
                               paste(utils::head(made, 5), collapse = ", ")),
                       type = "message", duration = 5)
    })
    observeEvent(input$multi_add, {
      st <- staged(); pick <- multi_pick(); targets <- input$multi_targets
      if (!length(pick)) { showNotification("Select at least one staged set.", type = "warning"); return() }
      if (!length(targets)) { showNotification("Select at least one set to add to.", type = "warning"); return() }
      ids <- unique(unlist(st[pick], use.names = FALSE))
      sets <- state$gene_sets
      for (t in targets) sets <- gene_set_commit(sets, t, ids, "append", source = staged_source())$sets
      state$gene_sets <- sets; snapshot_absent()
      .log(state, list(action = "gene_set_addto", targets = targets, source = staged_source(),
                       n = length(ids), params = staged_params()))
      showNotification(sprintf("Added %d gene(s) to %d set(s).", length(ids), length(targets)),
                       type = "message", duration = 4)
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

    # --- Export: write a SELECTED subset of the store to JSON / GMT / TSV ----
    # JSON is the faithful mirror (round-trips source/kind); GMT and TSV are
    # interchange formats. The Import source pill reads any of them back. A live,
    # capped preview shows exactly what will download (mirrors the Palette tab).
    output$export_ui <- renderUI({
      nms <- names(state$gene_sets)                     # re-renders on add/delete/rename
      if (!length(nms))
        return(helpText(class = "small text-muted", "Create a gene set to enable export."))
      cur <- shiny::isolate(input$export_which)
      sel <- if (!is.null(cur)) intersect(cur, nms) else nms
      tagList(
        # Cap the selected-items box height so many sets don't grow the UI.
        tags$style(HTML(sprintf(
          "#%s + .selectize-control .selectize-input { max-height: 110px; overflow-y: auto; }",
          ns("export_which")))),
        selectInput(ns("export_fmt"), "Format",
                    c("JSON (faithful, recommended)" = "json",
                      "GMT (MSigDB)" = "gmt", "TSV (long: set / id)" = "tsv"),
                    selected = shiny::isolate(input$export_fmt) %||% "json"),
        selectizeInput(ns("export_which"), "Sets to export", choices = nms,
                       selected = sel, multiple = TRUE),
        tags$div(class = "d-flex gap-2 mb-2",
          actionButton(ns("export_all"), "Select all", class = "btn-sm btn-outline-secondary"),
          actionButton(ns("export_none"), "Deselect all", class = "btn-sm btn-outline-secondary")),
        downloadButton(ns("export_dl"), "Download", class = "btn btn-sm fw-semibold",
                       style = .gs_action_style))
    })
    observeEvent(input$export_all,
      updateSelectizeInput(session, "export_which", selected = names(state$gene_sets)))
    observeEvent(input$export_none,
      updateSelectizeInput(session, "export_which", selected = character(0)))
    # The selected sets. An emptied `selectize` reports NULL (not character(0)),
    # so treat absent selection as "nothing selected" -- the initial all-selected
    # state comes from the render default (`selected = sel`), which the browser
    # reflects into input on first paint. (Mirrors the Palette export selector;
    # using NULL as an "all" sentinel would make Deselect-all silently export all.)
    export_sets <- reactive({
      nms <- names(state$gene_sets); if (!length(nms)) return(list())
      state$gene_sets[intersect(input$export_which %||% character(0), nms)]
    })
    export_txt <- reactive({
      switch(input$export_fmt %||% "json",
        json = gene_sets_to_json(export_sets()),
        gmt  = gene_sets_to_gmt(export_sets()),
        tsv  = gene_sets_to_tsv(export_sets()))
    })
    output$export_preview_ui <- renderUI({
      if (!length(state$gene_sets))
        return(helpText(class = "small text-muted", "Create a gene set to preview the export."))
      # Cap the previewed text (not just its box height) so a large store doesn't
      # push thousands of lines into the DOM; the full content still downloads.
      lines <- strsplit(export_txt(), "\n", fixed = TRUE)[[1]]
      cap <- 300L
      shown <- if (length(lines) > cap)
        paste0(paste(utils::head(lines, cap), collapse = "\n"),
               sprintf("\n... (%d more lines; full content downloads)", length(lines) - cap))
      else export_txt()
      tagList(
        helpText(class = "small text-muted mb-1",
                 sprintf("%d of %d set(s) selected.",
                         length(export_sets()), length(state$gene_sets))),
        tags$pre(class = "border rounded p-2 mb-0",
                 style = "max-height: 240px; overflow: auto; font-size: 0.8em; white-space: pre;",
                 shown))
    })
    output$export_dl <- downloadHandler(
      filename = function()
        sprintf("ddsdashboard-gene-sets-%s.%s", Sys.Date(), input$export_fmt %||% "json"),
      content = function(file) writeLines(export_txt(), file))

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
