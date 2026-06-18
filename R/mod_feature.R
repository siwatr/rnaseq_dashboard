# Page 1, "Feature info" tab: annotation (OrgDb, then an authoritative GTF) plus
# the shared metadata draft editor configured for rowData (features). Annotation
# and feature_length operate on the editor's DRAFT (composing with the user's
# unsaved edits) and commit only on Save -- they do not write straight to
# state$working, which would reset the draft and drop pending edits. The editable
# table (R/mod_meta_editor.R) handles feature_class and free edits. feature_length
# can be computed from the GTF (union over a chosen feature type) or adopted from
# an existing numeric column; either unlocks TPM/FPKM on the Assay tab.

.feature_editor_opts <- list(
  slot = "rowData", title = tags$h4("Feature information", class = "fs-6"), row_noun = "feature",
  allow_merge = FALSE, allow_row_rename = FALSE, bulk_class = TRUE
)

# OrgDb columns offered in the "Information to add" selector (values are OrgDb
# keytypes; each maps to a fixed rowData column, see .orgdb_target()).
.orgdb_col_choices <- c(
  "Gene symbol (-> <unit>_name)" = "SYMBOL",
  "Description (-> description)" = "GENENAME",
  "Ensembl id (-> ensembl_id)"  = "ENSEMBL",
  "Entrez id (-> entrez_id)"    = "ENTREZID",
  "Gene type (-> gene_biotype)" = "GENETYPE"
)

# Notification text after an annotation step (how much of the table it touched).
.append_msg <- function(source, matched, total) {
  pct <- if (total > 0L) round(100 * matched / total) else 0L
  sprintf("Appending %s annotation to %d of %d features (%d%%). Click Save to keep.",
          source, matched, total, pct)
}

mod_feature_ui <- function(id) {
  ns <- NS(id)
  # General feature settings (belong to neither annotation source) go in the
  # metadata card's own sidebar.
  feature_settings <- tagList(
    tags$strong("Feature settings"),
    selectInput(ns("feature_type"), "Feature unit",
                choices = c("gene", "transcript", "exon", "feature")),
    helpText("What each row represents; sets the <unit>_name column and labels."),
    tags$strong("Set feature_length from a column"),
    helpText("Unlocks TPM/FPKM from an existing numeric rowData column:"),
    uiOutput(ns("len_col_ui")),
    actionButton(ns("set_len"), "Set length from column"),
    hr()
  )
  # One navset: the metadata table and each annotation source are separate tabs,
  # so only one table + its sidebar is shown at a time.
  bslib::navset_card_pill(
    title = tags$h4("Feature info", class = "fs-6 mb-0 pe-3"),
    bslib::nav_panel(
      "Feature Metadata",
      meta_editor_ui(ns("editor"), .feature_editor_opts,
                     extra_sidebar = feature_settings,
                     extra_main = div(class = "mb-2", textOutput(ns("coverage"))))
    ),
    bslib::nav_panel(
      "OrgDb Annotation",
      bslib::layout_sidebar(
        sidebar = bslib::sidebar(
          title = "OrgDb", width = 320,
          selectInput(ns("organism"), "Organism",
                      c("Mouse (org.Mm.eg.db)" = "mouse", "Human (org.Hs.eg.db)" = "human")),
          selectInput(ns("id_type"), "Feature id type",
                      c("Auto-detect" = "auto", "Ensembl" = "ensembl",
                        "Entrez" = "entrez", "Symbol" = "symbol")),
          helpText(textOutput(ns("detected_id"), inline = TRUE)),
          selectizeInput(ns("orgdb_cols"), "Information to add", multiple = TRUE,
                         choices = .orgdb_col_choices, selected = c("SYMBOL", "GENENAME")),
          checkboxInput(ns("orgdb_flag"), "Flag mapped features (in_orgdb column)", value = TRUE),
          actionButton(ns("annotate"), "Annotate from OrgDb", class = "btn-primary")
        ),
        uiOutput(ns("orgdb_cov")),
        tags$small(class = "text-muted",
                   "Available to join (your feature ids resolved against the OrgDb):"),
        DT::DTOutput(ns("orgdb_preview"))
      )
    ),
    bslib::nav_panel(
      "GTF Annotation",
      bslib::layout_sidebar(
        sidebar = bslib::sidebar(
          title = "GTF", width = 320,
          helpText("GTF values override OrgDb where they match."),
          mod_gtf_reader_ui(ns("gtf"), preview = FALSE),
          uiOutput(ns("gtf_opts")),
          uiOutput(ns("gtf_len_ui"))
        ),
        uiOutput(ns("gtf_cov")),
        mod_gtf_reader_preview_ui(ns("gtf"))
      )
    )
  )
}

#' @param state the shared app-state object (see [new_app_state()]).
#' @return Invisible NULL.
mod_feature_server <- function(id, state) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    editor <- meta_editor_server("editor", state, .feature_editor_opts)
    gtf_obj <- mod_gtf_reader_server("gtf")   # confirmed (trimmed) GRanges, or NULL
    suppress_overwrite <- reactiveVal(FALSE)  # "don't warn again this session"
    pending_apply      <- reactiveVal(NULL)   # deferred op awaiting modal confirm

    # Apply fn(draft, ...) to the editor draft, reporting errors; returns the new
    # dds (invisibly via the editor) or NULL on failure.
    edit_draft <- function(fn) {
      d <- editor$draft(); if (is.null(d)) return(NULL)
      res <- tryCatch(fn(d), error = function(e) {
        showNotification(conditionMessage(e), type = "error", duration = NULL); NULL })
      if (is.null(res)) return(NULL)
      editor$set(res)
      newft <- detect_feature_type(res)
      state$meta <- utils::modifyList(state$meta, list(feature_type = newft$feature_type))
      res
    }

    # Run `do()` immediately, unless it would overwrite columns that already hold
    # values -- then ask first (Proceed/Cancel + a session-wide "don't warn").
    guarded <- function(targets, do) {
      d <- editor$draft(); if (is.null(d)) return(invisible())
      ow <- annotation_overwrites(d, targets)
      if (length(ow) && !isTRUE(suppress_overwrite())) {
        pending_apply(do)
        showModal(modalDialog(
          title = "Overwrite existing annotation?",
          tags$p(sprintf("These column(s) already hold values and will be updated where the source matches: %s.",
                         paste(ow, collapse = ", "))),
          checkboxInput(ns("ow_suppress"), "Don't warn again this session", FALSE),
          footer = tagList(modalButton("Cancel"),
                           actionButton(ns("ow_proceed"), "Proceed", class = "btn-warning"))
        ))
      } else do()
    }

    observeEvent(input$ow_proceed, {
      if (isTRUE(input$ow_suppress)) suppress_overwrite(TRUE)
      op <- pending_apply(); pending_apply(NULL); removeModal()
      if (is.function(op)) op()
    })

    output$detected_id <- renderText({
      d <- editor$draft(); req(d)
      sprintf("Detected: %s", detect_id_type(rownames(d)))
    })

    # Feature-unit selector (moved here from Load): keep it in sync with meta and
    # write back as a labeling change only (no data_version bump).
    observeEvent(state$working, {
      updateSelectInput(session, "feature_type", selected = state_meta(state)$feature_type)
    })
    observeEvent(input$feature_type, {
      req(state$working)
      state$meta <- utils::modifyList(state$meta, list(feature_type = input$feature_type))
    }, ignoreInit = TRUE)

    # Columns the user picked, falling back to the historical default so the page
    # works before the selectize input reports (and in headless tests).
    orgdb_cols <- reactive(input$orgdb_cols %||% c("SYMBOL", "GENENAME"))

    # Match-coverage banners above each preview (how many of the user's feature
    # ids the source resolves), coloured red/amber/green by completeness.
    output$orgdb_cov <- renderUI({
      req(editor$draft())
      id_type <- if (identical(input$id_type, "auto")) NULL else input$id_type
      cnt <- tryCatch(orgdb_match_count(editor$draft(), organism = input$organism,
                                        id_type = id_type, columns = orgdb_cols()),
                      error = function(e) NULL)
      req(cnt)
      .coverage_banner(cnt$matched, cnt$total, "feature IDs", "OrgDb")
    })

    output$gtf_cov <- renderUI({
      g <- gtf_obj(); req(g, editor$draft())
      cnt <- tryCatch(gtf_match_count(editor$draft(), g, match_col = input$gtf_match %||% "auto"),
                      error = function(e) NULL)
      req(cnt)
      .coverage_banner(cnt$matched, cnt$total, "feature IDs", "the GTF")
    })

    # Non-committing preview of what OrgDb makes available to join: the join key
    # (id) plus one column per selected piece of information, for the first rows.
    output$orgdb_preview <- DT::renderDT({
      req(editor$draft())
      id_type <- if (identical(input$id_type, "auto")) NULL else input$id_type
      df <- tryCatch(
        orgdb_annotation_preview(editor$draft(), organism = input$organism,
                                 id_type = id_type, feature_type = state_meta(state)$feature_type,
                                 columns = orgdb_cols(), n = 100L),
        error = function(e) NULL)
      req(df)
      dt_table(df)
    })

    observeEvent(input$annotate, {
      req(editor$draft())
      organism <- input$organism
      id_type <- if (identical(input$id_type, "auto")) NULL else input$id_type
      ft <- state_meta(state)$feature_type
      cols <- orgdb_cols()
      flag <- if (isTRUE(input$orgdb_flag)) "in_orgdb" else NULL
      do <- function() {
        res <- edit_draft(function(d)
          annotate_with_orgdb(d, organism = organism, id_type = id_type,
                              feature_type = ft, columns = cols, matched_col = flag))
        req(res)
        cov <- annotation_coverage(res, paste0(detect_feature_type(res)$feature_type, "_name"))
        showNotification(.append_msg("OrgDb", cov$matched, cov$total), type = "message")
      }
      targets <- vapply(cols, .orgdb_target, character(1), feature_type = ft)
      guarded(unname(targets), do)
    })

    # --- GTF attribute import (draft) ----------------------------------------
    output$gtf_opts <- renderUI({
      g <- gtf_obj(); req(g)
      cols <- available_gtf_columns(g)
      tagList(
        tags$hr(), tags$strong("Import GTF attributes"),
        selectInput(ns("gtf_match"), "Match dds rows by",
                    c("Auto (id, then name)" = "auto", stats::setNames(cols, cols))),
        selectizeInput(ns("gtf_import"), "Import columns", choices = cols, multiple = TRUE,
                       selected = intersect(c("gene_name", "seqnames", "gene_biotype"), cols)),
        checkboxInput(ns("gtf_flag"), "Flag matched features (in_gtf column)", value = TRUE),
        actionButton(ns("apply_gtf"), "Import into features", class = "btn-primary")
      )
    })

    observeEvent(input$apply_gtf, {
      req(editor$draft(), gtf_obj())
      g <- gtf_obj()
      ft <- state_meta(state)$feature_type
      import_cols <- input$gtf_import
      flag <- if (isTRUE(input$gtf_flag)) "in_gtf" else NULL
      targets <- vapply(import_cols, function(c)
        switch(c, gene_name = paste0(ft, "_name"), seqnames = "chromosome", c),
        character(1))
      do <- function() {
        report <- NULL
        res <- edit_draft(function(d) {
          out <- annotate_with_gtf(d, g, match_col = input$gtf_match, import_cols = import_cols,
                                   compute_length = FALSE, feature_type = ft, matched_col = flag)
          report <<- out$report
          out$dds
        })
        req(res)
        showNotification(.append_msg("GTF", report$matched, report$total), type = "message")
      }
      guarded(targets, do)
    })

    # --- feature_length: from an existing numeric column ---------------------
    output$len_col_ui <- renderUI({
      d <- editor$draft(); req(d)
      rd <- SummarizedExperiment::rowData(d)
      num <- colnames(rd)[vapply(as.list(rd), is.numeric, logical(1))]
      selectInput(ns("len_col"), NULL, choices = num)
    })

    observeEvent(input$set_len, {
      req(editor$draft(), input$len_col)
      do <- function() {
        res <- edit_draft(function(d) set_feature_length_from_column(d, input$len_col))
        req(res)
        showNotification(sprintf("feature_length set from '%s'. Click Save to keep.", input$len_col),
                         type = "message")
      }
      guarded("feature_length", do)
    })

    # --- feature_length: computed from the loaded GTF ------------------------
    output$gtf_len_ui <- renderUI({
      g <- gtf_obj(); req(g)
      types <- gtf_feature_types(g)
      tagList(
        tags$hr(), tags$strong("Compute feature length"),
        helpText("Union length over a feature type (exon for mature mRNA):"),
        selectInput(ns("gtf_len_type"), NULL, choices = types,
                    selected = if ("exon" %in% types) "exon" else types[[1]]),
        actionButton(ns("compute_len"), "Compute length from GTF")
      )
    })

    observeEvent(input$compute_len, {
      req(editor$draft(), gtf_obj())
      g <- gtf_obj()
      ft <- state_meta(state)$feature_type
      do <- function() {
        report <- NULL
        res <- edit_draft(function(d) {
          out <- annotate_with_gtf(d, g, match_col = input$gtf_match %||% "auto",
                                   import_cols = NULL, compute_length = TRUE,
                                   length_type = input$gtf_len_type, feature_type = ft)
          report <<- out$report
          out$dds
        })
        req(res)
        showNotification(sprintf("feature_length set for %d/%d features%s. Click Save to keep.",
                                 report$length_set, report$total,
                                 if (report$length_complete) " - TPM/FPKM available"
                                 else " - incomplete, TPM/FPKM need all features"),
                         type = "message", duration = NULL)
      }
      guarded("feature_length", do)
    })

    output$coverage <- renderText({
      d <- editor$draft(); req(d)
      name_col <- paste0(state_meta(state)$feature_type, "_name")
      cov <- annotation_coverage(d, name_col)
      name_msg <- if (cov$matched == 0L) "Names: not annotated yet."
                  else sprintf("Names: %d of %d endogenous features.", cov$matched, cov$total)
      len_msg <- if (has_feature_length(d)) "feature_length: complete (TPM/FPKM available)."
                 else "feature_length: not set (TPM/FPKM unavailable)."
      paste(name_msg, len_msg, "(draft - Save to commit.)")
    })

    invisible(NULL)
  })
}
