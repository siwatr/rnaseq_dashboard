# Page 1, "Assay" tab: inspect assays and add normalized ones. CPM is always
# available; TPM/FPKM require a complete feature_length (populated by annotation
# / GTF). Size factors are NOT set here -- they are a separate DESeq2 median-of-
# ratios normalization (unrelated to these library-size/length assays) managed on
# the "Size factors" tab. Backed by R/assay_helpers.R.

mod_assay_ui <- function(id) {
  ns <- NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      title = tags$h4("Add assays", class = "fs-6 mb-0"), width = 250,
      uiOutput(ns("controls")),
      actionButton(ns("apply"), "Add assays", class = "btn-primary")
    ),
    bslib::card(
      bslib::card_header(tags$h4("Assay information", class = "fs-6 mb-0")),
      verbatimTextOutput(ns("info"))
    )
  )
}

#' @param state the shared app-state object (see [new_app_state()]).
#' @return Invisible NULL.
mod_assay_server <- function(id, state) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    output$controls <- renderUI({
      req(state$working)
      if (has_feature_length(state$working)) {
        checkboxGroupInput(ns("assays"), "Assays to compute",
                           choices = c("CPM", "TPM", "FPKM"), selected = "CPM")
      } else {
        tagList(
          checkboxGroupInput(ns("assays"), "Assays to compute",
                             choices = c("CPM"), selected = "CPM"),
          helpText("TPM/FPKM need a complete feature_length: add it via annotation or a GTF.")
        )
      }
    })

    suppress_overwrite <- reactiveVal(FALSE)   # "don't warn again this session"
    pending_assays     <- reactiveVal(NULL)    # selection awaiting modal confirm

    do_apply <- function(sel) {
      ok <- tryCatch({
        state_mutate(state, function(d) add_normalized_assays(d, sel),
                     action = list(action = "add_assays", assays = sel))
        TRUE
      }, error = function(e) { showNotification(conditionMessage(e), type = "error"); FALSE })
      if (!isTRUE(ok)) return(invisible())
      showNotification(sprintf("Updated assays: %s.",
                               paste(SummarizedExperiment::assayNames(state$working), collapse = ", ")),
                       type = "message")
    }

    # Recomputing an assay that already exists overwrites it -- warn first
    # (Proceed/Cancel + a session-wide "don't warn"), mirroring the annotation tab.
    observeEvent(input$apply, {
      req(state$working, input$assays)
      sel <- input$assays
      existing <- intersect(sel, SummarizedExperiment::assayNames(state$working))
      if (length(existing) && !isTRUE(suppress_overwrite())) {
        pending_assays(sel)
        showModal(modalDialog(
          title = "Overwrite existing assay(s)?",
          tags$p(sprintf("These assay(s) already exist and will be recomputed and overwritten: %s.",
                         paste(existing, collapse = ", "))),
          checkboxInput(ns("assay_ow_suppress"), "Don't warn again this session", FALSE),
          footer = tagList(modalButton("Cancel"),
                           actionButton(ns("assay_ow_proceed"), "Proceed", class = "btn-warning"))
        ))
      } else do_apply(sel)
    })
    observeEvent(input$assay_ow_proceed, {
      if (isTRUE(input$assay_ow_suppress)) suppress_overwrite(TRUE)
      sel <- pending_assays(); pending_assays(NULL); removeModal()
      if (!is.null(sel)) do_apply(sel)
    })

    output$info <- renderPrint({
      if (is.null(state$working)) { cat("No dataset loaded."); return() }
      dds <- state$working
      cat("Assays:", paste(SummarizedExperiment::assayNames(dds), collapse = ", "), "\n")
      cat("Dimensions:", nrow(dds), "features x", ncol(dds), "samples\n")
    })

    invisible(NULL)
  })
}
