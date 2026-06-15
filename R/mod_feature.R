# Page 1, "Feature info" tab: annotation (OrgDb, then an authoritative GTF) plus
# the shared metadata draft editor configured for rowData (features). Annotation
# applies immediately via state_mutate(); the editable table (R/mod_meta_editor.R)
# handles feature_class and free edits. GTF annotation also sources feature_length
# (union over a chosen feature type, or adopted from an existing numeric column),
# which unlocks TPM/FPKM on the Assay tab.

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
    helpText("A GTF overrides OrgDb for matching features and can supply feature_length."),
    fileInput(ns("gtf_file"), "Upload GTF/GFF (gzip OK)",
              accept = c(".gtf", ".gff", ".gff3", ".gtf.gz", ".gff.gz", ".gff3.gz",
                         ".gz", "application/gzip", "application/x-gzip", "text/plain")),
    helpText("Large uncompressed GTFs may exceed the upload limit - upload gzipped, ",
             "or read a local file already on this machine:"),
    textInput(ns("gtf_path"), NULL, placeholder = "/path/to/annotation.gtf(.gz)"),
    actionButton(ns("read_gtf"), "Read GTF"),
    uiOutput(ns("gtf_opts")),
    hr(),
    tags$strong("Set feature_length from a column"),
    helpText("Use an existing numeric rowData column (e.g. length from your quantifier)."),
    uiOutput(ns("len_col_ui")),
    actionButton(ns("set_len"), "Set length from column"),
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
    meta_editor_server("editor", state, .feature_editor_opts)
    gtf_obj <- reactiveVal(NULL)

    output$detected_id <- renderText({
      req(state$working)
      sprintf("Detected: %s", detect_id_type(rownames(state$working)))
    })

    observeEvent(input$annotate, {
      req(state$working)
      organism <- input$organism
      id_type <- if (identical(input$id_type, "auto")) NULL else input$id_type
      ft <- state_meta(state)$feature_type
      ok <- tryCatch({
        state_mutate(state,
          function(d) annotate_with_orgdb(d, organism = organism, id_type = id_type, feature_type = ft),
          action = list(action = "annotate_orgdb", organism = organism))
        TRUE
      }, error = function(e) { showNotification(conditionMessage(e), type = "error", duration = NULL); FALSE })
      req(ok)
      newft <- detect_feature_type(state$working)
      state$meta <- utils::modifyList(state$meta, list(feature_type = newft$feature_type))
      cov <- annotation_coverage(state$working, paste0(newft$feature_type, "_name"))
      showNotification(sprintf("Annotated %d of %d endogenous features.", cov$matched, cov$total),
                       type = "message")
    })

    # --- GTF -----------------------------------------------------------------
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

    output$gtf_opts <- renderUI({
      g <- gtf_obj(); req(g)
      cols <- available_gtf_columns(g)
      types <- gtf_feature_types(g)
      tagList(
        selectInput(ns("gtf_match"), "Match dds rows by",
                    c("Auto (id, then name)" = "auto", stats::setNames(cols, cols))),
        selectizeInput(ns("gtf_import"), "Import columns", choices = cols, multiple = TRUE,
                       selected = intersect(c("gene_name", "seqnames", "gene_biotype"), cols)),
        checkboxInput(ns("gtf_len"), "Compute feature_length from this GTF", value = TRUE),
        selectInput(ns("gtf_len_type"), "Length feature type", choices = types,
                    selected = if ("exon" %in% types) "exon" else types[[1]]),
        actionButton(ns("apply_gtf"), "Apply GTF annotation", class = "btn-primary")
      )
    })

    observeEvent(input$apply_gtf, {
      req(state$working, gtf_obj())
      g <- gtf_obj()
      ft <- state_meta(state)$feature_type
      report <- NULL
      ok <- tryCatch({
        res <- annotate_with_gtf(state$working, g, match_col = input$gtf_match,
                                 import_cols = input$gtf_import,
                                 compute_length = isTRUE(input$gtf_len),
                                 length_type = input$gtf_len_type, feature_type = ft)
        report <- res$report
        state_mutate(state, function(.) res$dds,
                     action = list(action = "annotate_gtf", match_col = input$gtf_match,
                                   import_cols = input$gtf_import,
                                   length_type = if (isTRUE(input$gtf_len)) input$gtf_len_type else NA))
        TRUE
      }, error = function(e) { showNotification(conditionMessage(e), type = "error", duration = NULL); FALSE })
      req(ok)
      newft <- detect_feature_type(state$working)
      state$meta <- utils::modifyList(state$meta, list(feature_type = newft$feature_type))
      msg <- sprintf("GTF: matched %d of %d features.", report$matched, report$total)
      if (isTRUE(input$gtf_len)) {
        msg <- paste0(msg, sprintf(" feature_length set for %d/%d%s.", report$length_set, report$total,
                      if (report$length_complete) " - TPM/FPKM available"
                      else " - incomplete, TPM/FPKM need all features"))
      }
      showNotification(msg, type = "message", duration = NULL)
    })

    # --- feature_length from an existing numeric column ----------------------
    output$len_col_ui <- renderUI({
      req(state$working)
      rd <- SummarizedExperiment::rowData(state$working)
      num <- colnames(rd)[vapply(as.list(rd), is.numeric, logical(1))]
      selectInput(ns("len_col"), NULL, choices = num)
    })

    observeEvent(input$set_len, {
      req(state$working, input$len_col)
      ok <- tryCatch({
        state_mutate(state, function(d) set_feature_length_from_column(d, input$len_col),
                     action = list(action = "set_feature_length", column = input$len_col))
        TRUE
      }, error = function(e) { showNotification(conditionMessage(e), type = "error", duration = NULL); FALSE })
      if (ok) showNotification(sprintf("feature_length set from '%s'.", input$len_col), type = "message")
    })

    output$coverage <- renderText({
      req(state$working)
      name_col <- paste0(state_meta(state)$feature_type, "_name")
      cov <- annotation_coverage(state$working, name_col)
      name_msg <- if (cov$matched == 0L) "Names: not annotated yet."
                  else sprintf("Names: %d of %d endogenous features.", cov$matched, cov$total)
      len_msg <- if (has_feature_length(state$working)) "feature_length: complete (TPM/FPKM available)."
                 else "feature_length: not set (TPM/FPKM unavailable)."
      paste(name_msg, len_msg)
    })

    invisible(NULL)
  })
}
