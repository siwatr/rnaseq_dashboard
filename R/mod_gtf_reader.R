# GTF reader/filter submodule. Read a GTF/GFF (upload or a local path; gzip-aware),
# let the user keep only the feature types and mcols columns they need, then
# Confirm to trim the GRanges and free the full parse from the session (the parse
# itself still spikes -- rtracklayer is not streaming -- but only the trimmed
# object persists). Returns the confirmed (trimmed) GRanges reactiveVal so the
# host page composes annotation on top of it; tests inject with gtf(<GRanges>).

#' Reader controls UI
#'
#' @param id Module id.
#' @param preview If TRUE, also embed the preview table inline; set FALSE to place
#'   the preview elsewhere via [mod_gtf_reader_preview_ui()] (e.g. a tab's main area).
#' @noRd
mod_gtf_reader_ui <- function(id, preview = TRUE) {
  ns <- NS(id)
  tagList(
    tags$strong("GTF upload"),
    fileInput(ns("file"), "Upload GTF/GFF (gzip OK)",
              accept = c(".gtf", ".gff", ".gff3", ".gtf.gz", ".gff.gz", ".gff3.gz",
                         ".gz", "application/gzip", "application/x-gzip", "text/plain")),
    helpText("Large uncompressed GTFs may exceed the upload limit - upload gzipped, ",
             "or read a local file already on this machine:"),
    textInput(ns("path"), NULL, placeholder = "/path/to/annotation.gtf(.gz)"),
    div(class = "d-flex gap-2 align-items-start",
        actionButton(ns("read"), "Read GTF"),
        uiOutput(ns("remove_ui"), inline = TRUE)),
    uiOutput(ns("status")),
    uiOutput(ns("select")),
    if (isTRUE(preview)) uiOutput(ns("preview_ui"))
  )
}

#' The reader's preview table alone (for placing in a separate region)
#' @param id Module id (same as passed to [mod_gtf_reader_ui()]).
#' @noRd
mod_gtf_reader_preview_ui <- function(id) {
  uiOutput(NS(id)("preview_ui"))
}

#' @param id Module id.
#' @return The `confirmed` reactiveVal yielding the trimmed `GRanges` (or NULL).
#' @noRd
mod_gtf_reader_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    raw       <- reactiveVal(NULL)   # full parse, transient
    confirmed <- reactiveVal(NULL)   # trimmed object, returned

    observeEvent(input$read, {
      # Prefer a local path (no browser upload, no size cap); else the upload.
      path <- NULL
      if (nzchar(input$path %||% "")) {
        if (!file.exists(input$path)) { showNotification("No file at that path.", type = "error"); return() }
        path <- input$path
      } else if (!is.null(input$file)) {
        path <- file.path(tempdir(), basename(input$file$name))
        file.copy(input$file$datapath, path, overwrite = TRUE)
      } else {
        showNotification("Upload a file or enter a local path first.", type = "warning"); return()
      }
      g <- tryCatch(
        withProgress(message = "Reading GTF...", value = 0.5, import_gtf(path)),
        error = function(e) { showNotification(conditionMessage(e), type = "error", duration = NULL); NULL })
      req(g)
      raw(g); confirmed(NULL)
      showNotification(sprintf("Read GTF: %d records; types: %s. Select what to keep, then Confirm.",
                               length(g), paste(gtf_feature_types(g), collapse = ", ")), type = "message")
    })

    # Filtering stays available for the life of the loaded GTF (driven by the
    # current object, raw or already-trimmed), so the user can keep narrowing it.
    # Re-applying can only trim further; to recover a dropped column they Remove
    # the GTF and re-read it.
    output$select <- renderUI({
      g <- raw() %||% confirmed(); req(g)
      types <- gtf_feature_types(g)
      cols  <- available_gtf_columns(g)
      def_types <- intersect(c("gene", "transcript", "exon"), types)
      sel_types <- if (length(def_types)) def_types else types
      suggested <- intersect(c("gene_id", "gene_name", "transcript_id", "gene_biotype", "seqnames"), cols)
      tagList(
        tags$hr(), tags$strong("GTF filtering"),
        selectizeInput(ns("keep_types"), "Feature types to keep", choices = types,
                       selected = sel_types, multiple = TRUE),
        selectizeInput(ns("keep_cols"), "Columns to keep", choices = cols,
                       selected = if (length(suggested)) suggested else cols, multiple = TRUE),
        div(class = "form-text", "Leave a field empty to keep all. type and id/name keys are always kept."),
        actionButton(ns("confirm"), "Apply filter", class = "btn-primary")
      )
    })

    observeEvent(input$confirm, {
      g <- raw() %||% confirmed(); req(g)            # filter whichever is current
      types <- if (length(input$keep_types)) input$keep_types else gtf_feature_types(g)
      g2 <- g[as.character(g$type) %in% types]
      md_cols <- colnames(S4Vectors::mcols(g2))
      want <- if (length(input$keep_cols)) input$keep_cols else md_cols
      # Always retain `type` (needed to subset by feature type for length) and any
      # present binder columns, so the match key survives whatever the user picks
      # at apply time. seqnames is intrinsic to the GRanges (not in mcols).
      protect <- intersect(c("type", "gene_id", "gene_name", "transcript_id"), md_cols)
      keep <- intersect(union(setdiff(want, "seqnames"), protect), md_cols)
      S4Vectors::mcols(g2) <- S4Vectors::mcols(g2)[, keep, drop = FALSE]
      confirmed(g2)
      raw(NULL); gc()                                  # release the full parse
      showNotification(sprintf("Filtered: %d records, %d column(s) kept.",
                               length(g2), length(keep)), type = "message")
    })

    # Drop the loaded GTF entirely, freeing the session (re-read to use it again).
    observeEvent(input$remove, {
      confirmed(NULL); raw(NULL); gc()
      showNotification("GTF removed from session.", type = "message")
    })

    # Remove button lives with the upload group, available whenever a GTF is loaded.
    output$remove_ui <- renderUI({
      req(raw() %||% confirmed())
      bslib::tooltip(
        actionButton(ns("remove"), "Remove GTF", class = "btn-outline-danger"),
        "Free the session; re-read the GTF to use it again.")
    })

    output$status <- renderUI({
      g <- confirmed() %||% raw(); if (is.null(g)) return(NULL)
      ready <- !is.null(confirmed())
      div(class = if (ready) "form-text text-success" else "form-text",
          sprintf("GTF loaded: %d records, %d column(s).",
                  length(g), length(available_gtf_columns(g))))
    })

    # ~20-row preview of the loaded GTF (raw while selecting, else the trimmed
    # object) so users can see columns/values when choosing what to keep/import.
    output$preview_ui <- renderUI({
      g <- raw() %||% confirmed(); req(g)
      tagList(tags$small(class = "text-muted", "GTF preview (first rows):"),
              DT::DTOutput(ns("preview")))
    })
    output$preview <- DT::renderDT({
      g <- raw() %||% confirmed(); req(g)
      DT::datatable(gtf_preview(g, n = 100L), rownames = FALSE,
                    options = list(dom = "ltp", pageLength = 10, scrollX = TRUE,
                                   lengthMenu = list(c(10, 25, 50, 100),
                                                     c("10", "25", "50", "100"))))
    })

    confirmed
  })
}
