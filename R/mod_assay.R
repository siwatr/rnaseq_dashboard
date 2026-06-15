# Page 1, "Assay" tab: inspect assays and add normalized ones. CPM is always
# available; TPM/FPKM require a complete feature_length (populated by annotation
# / GTF). Adding assays also (re)estimates size factors on endogenous genes.
# Backed by R/assay_helpers.R.

mod_assay_ui <- function(id) {
  ns <- NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      title = "Add assays", width = 320,
      uiOutput(ns("controls")),
      actionButton(ns("apply"), "Add assays / update size factors", class = "btn-primary")
    ),
    bslib::card(
      bslib::card_header("Assay information"),
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

    observeEvent(input$apply, {
      req(state$working, input$assays)
      sel <- input$assays
      ok <- tryCatch({
        state_mutate(state, function(d) {
          d <- add_normalized_assays(d, sel)
          estimate_size_factors_endogenous(d)
        }, action = list(action = "add_assays", assays = sel))
        TRUE
      }, error = function(e) { showNotification(conditionMessage(e), type = "error"); FALSE })
      req(ok)
      showNotification(sprintf("Updated assays: %s.",
                               paste(SummarizedExperiment::assayNames(state$working), collapse = ", ")),
                       type = "message")
    })

    output$info <- renderPrint({
      if (is.null(state$working)) { cat("No dataset loaded."); return() }
      dds <- state$working
      cat("Assays:", paste(SummarizedExperiment::assayNames(dds), collapse = ", "), "\n")
      cat("Dimensions:", nrow(dds), "features x", ncol(dds), "samples\n")
      sf <- tryCatch(DESeq2::sizeFactors(dds), error = function(e) NULL)
      cat("Size factors:", if (is.null(sf)) "not set" else "set (endogenous)", "\n")
    })

    invisible(NULL)
  })
}
