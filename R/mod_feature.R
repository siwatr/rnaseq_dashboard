# Page 1, "Feature info" tab: annotation (OrgDb, then an authoritative GTF) plus
# the shared metadata draft editor configured for rowData (features). Annotation
# and feature_length operate on the editor's DRAFT (composing with the user's
# unsaved edits) and commit only on Save -- they do not write straight to
# state$working, which would reset the draft and drop pending edits. The editable
# table (R/mod_meta_editor.R) handles feature_class and free edits. feature_length
# can be computed from the GTF (union over a chosen feature type) or adopted from
# an existing numeric column; either unlocks TPM/FPKM on the Assay tab.

.feature_editor_opts <- list(
  slot = "rowData", title = "Feature information", row_noun = "feature",
  allow_merge = FALSE, allow_row_rename = FALSE, bulk_class = TRUE
)

mod_feature_ui <- function(id) {
  ns <- NS(id)
  annotation_ui <- tagList(
    tags$strong("Annotate from OrgDb"),
    selectInput(ns("organism"), "Organism",
                c("Mouse (org.Mm.eg.db)" = "mouse", "Human (org.Hs.eg.db)" = "human")),
    selectInput(ns("id_type"), "Feature id type",
                c("Auto-detect" = "auto", "Ensembl" = "ensembl",
                  "Entrez" = "entrez", "Symbol" = "symbol")),
    helpText(textOutput(ns("detected_id"), inline = TRUE)),
    actionButton(ns("annotate"), "Annotate from OrgDb", class = "btn-primary"),
    hr(),
    tags$strong("Annotate from GTF"),
    helpText("A GTF overrides OrgDb for matching features. Edits land in the draft - click Save to keep."),
    fileInput(ns("gtf_file"), "Upload GTF/GFF (gzip OK)",
              accept = c(".gtf", ".gff", ".gff3", ".gtf.gz", ".gff.gz", ".gff3.gz",
                         ".gz", "application/gzip", "application/x-gzip", "text/plain")),
    helpText("Large uncompressed GTFs may exceed the upload limit - upload gzipped, ",
             "or read a local file already on this machine:"),
    textInput(ns("gtf_path"), NULL, placeholder = "/path/to/annotation.gtf(.gz)"),
    actionButton(ns("read_gtf"), "Read GTF"),
    uiOutput(ns("gtf_opts")),
    hr(),
    tags$strong("Set feature_length"),
    helpText("Unlocks TPM/FPKM. From an existing numeric rowData column:"),
    uiOutput(ns("len_col_ui")),
    actionButton(ns("set_len"), "Set length from column"),
    uiOutput(ns("gtf_len_ui")),
    hr()
  )
  meta_editor_ui(ns("editor"), .feature_editor_opts,
                 extra_sidebar = annotation_ui,
                 extra_main = div(class = "mb-2", textOutput(ns("coverage"))))
}

#' @param state the shared app-state object (see [new_app_state()]).
#' @return Invisible NULL.
mod_feature_server <- function(id, state) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    editor <- meta_editor_server("editor", state, .feature_editor_opts)
    gtf_obj <- reactiveVal(NULL)

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

    output$detected_id <- renderText({
      d <- editor$draft(); req(d)
      sprintf("Detected: %s", detect_id_type(rownames(d)))
    })

    observeEvent(input$annotate, {
      req(editor$draft())
      organism <- input$organism
      id_type <- if (identical(input$id_type, "auto")) NULL else input$id_type
      ft <- state_meta(state)$feature_type
      res <- edit_draft(function(d)
        annotate_with_orgdb(d, organism = organism, id_type = id_type, feature_type = ft))
      req(res)
      cov <- annotation_coverage(res, paste0(detect_feature_type(res)$feature_type, "_name"))
      showNotification(sprintf("Annotated %d of %d endogenous features. Click Save to keep.",
                               cov$matched, cov$total), type = "message")
    })

    # --- GTF read ------------------------------------------------------------
    observeEvent(input$read_gtf, {
      # Prefer a local path (no browser upload, no size cap); else the upload.
      # Either way a real file extension drives format + gzip handling.
      path <- NULL
      if (nzchar(input$gtf_path %||% "")) {
        if (!file.exists(input$gtf_path)) {
          showNotification("No file at that path.", type = "error"); return()
        }
        path <- input$gtf_path
      } else if (!is.null(input$gtf_file)) {
        path <- file.path(tempdir(), basename(input$gtf_file$name))
        file.copy(input$gtf_file$datapath, path, overwrite = TRUE)
      } else {
        showNotification("Upload a file or enter a local path first.", type = "warning"); return()
      }
      g <- tryCatch(
        withProgress(message = "Reading GTF...", value = 0.5, import_gtf(path)),
        error = function(e) { showNotification(conditionMessage(e), type = "error", duration = NULL); NULL })
      req(g)
      gtf_obj(g)
      showNotification(sprintf("Read GTF: %d records; types: %s.", length(g),
                               paste(gtf_feature_types(g), collapse = ", ")), type = "message")
    })

    # --- GTF attribute import (draft) ----------------------------------------
    output$gtf_opts <- renderUI({
      g <- gtf_obj(); req(g)
      cols <- available_gtf_columns(g)
      tagList(
        selectInput(ns("gtf_match"), "Match dds rows by",
                    c("Auto (id, then name)" = "auto", stats::setNames(cols, cols))),
        selectizeInput(ns("gtf_import"), "Import columns", choices = cols, multiple = TRUE,
                       selected = intersect(c("gene_name", "seqnames", "gene_biotype"), cols)),
        actionButton(ns("apply_gtf"), "Apply GTF annotation", class = "btn-primary")
      )
    })

    observeEvent(input$apply_gtf, {
      req(editor$draft(), gtf_obj())
      g <- gtf_obj()
      ft <- state_meta(state)$feature_type
      report <- NULL
      res <- edit_draft(function(d) {
        out <- annotate_with_gtf(d, g, match_col = input$gtf_match,
                                 import_cols = input$gtf_import,
                                 compute_length = FALSE, feature_type = ft)
        report <<- out$report
        out$dds
      })
      req(res)
      showNotification(sprintf("GTF: matched %d of %d features. Click Save to keep.",
                               report$matched, report$total), type = "message")
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
      res <- edit_draft(function(d) set_feature_length_from_column(d, input$len_col))
      req(res)
      showNotification(sprintf("feature_length set from '%s'. Click Save to keep.", input$len_col),
                       type = "message")
    })

    # --- feature_length: computed from the loaded GTF ------------------------
    output$gtf_len_ui <- renderUI({
      g <- gtf_obj(); req(g)
      types <- gtf_feature_types(g)
      tagList(
        helpText("Or compute union length from the loaded GTF over a feature type:"),
        selectInput(ns("gtf_len_type"), NULL, choices = types,
                    selected = if ("exon" %in% types) "exon" else types[[1]]),
        actionButton(ns("compute_len"), "Compute length from GTF")
      )
    })

    observeEvent(input$compute_len, {
      req(editor$draft(), gtf_obj())
      g <- gtf_obj()
      ft <- state_meta(state)$feature_type
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
