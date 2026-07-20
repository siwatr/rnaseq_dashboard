# Page 1, "Size factors" tab (Input group, after Assay). Size factors are a DESeq2
# median-of-ratios normalization -- deliberately separate from the CPM/TPM/FPKM
# assays (which are library-size/length normalizations that do NOT use size
# factors). This tab is where the user chooses the control-gene set (endogenous /
# spike-in / a custom set via the shared gene search) and the estimator `type`,
# then estimates. The config is carried on the dds (metadata) so it survives edits
# and structural re-estimations reuse it. DE / PCA / Expression are consumers.
# Backed by R/assay_helpers.R (estimate_size_factors + sizefactor_config).

mod_sizefactors_ui <- function(id) {
  ns <- NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      title = tags$h4("Estimate size factors", class = "fs-6 mb-0"), width = 300,
      radioButtons(ns("sf_control"), "Control genes",
                   c("Endogenous (default)" = "endogenous",
                     "Spike-in" = "spike_in", "Custom set" = "custom"),
                   selected = "endogenous"),
      conditionalPanel(
        sprintf("input['%s'] == 'custom'", ns("sf_control")),
        gene_search_ui(ns, "sf", multiple = TRUE,
                       search_modes = c("exact", "contains", "regex"),
                       label = "Control-gene set",
                       placeholder = "e.g. Actb, Gapdh, ENSG...")),
      selectInput(ns("sf_type"), "Estimator (type)",
                  c("ratio (median-of-ratios)" = "ratio",
                    "poscounts" = "poscounts", "iterate" = "iterate"),
                  selected = "ratio"),
      actionButton(ns("sf_estimate"), "Estimate size factors", class = "btn-primary"),
      helpText(class = "small text-muted mt-2",
               "Size factors normalize for sequencing depth (used by DE and normalized log-counts). Re-estimating changes only samples' scaling, not the counts.")
    ),
    bslib::card(
      bslib::card_header(tags$h4("Current size factors", class = "fs-6 mb-0")),
      uiOutput(ns("status")),
      DT::DTOutput(ns("table"))
    )
  )
}

#' @param state the shared app-state object (see [new_app_state()]).
#' @return Invisible NULL.
mod_sizefactors_server <- function(id, state) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    sf_search <- gene_search_server(input, output, session, state, "sf",
                                    multiple = TRUE,
                                    search_modes = c("exact", "contains", "regex"))

    # Reflect the dataset's stored config in the inputs (on load / any data change;
    # a no-op update when the config is unchanged, so unrelated edits don't disrupt).
    observeEvent(state$working, {
      req(state$working)
      cfg <- sizefactor_config(state$working)
      updateRadioButtons(session, "sf_control", selected = cfg$control)
      updateSelectInput(session, "sf_type", selected = cfg$type)
    })

    # Build the config the sidebar currently describes.
    pending_config <- reactive({
      ctrl <- input$sf_control %||% "endogenous"
      custom <- if (identical(ctrl, "custom")) sf_search()$ids else character(0)
      list(control = ctrl, custom_ids = custom,
           type = input$sf_type %||% "ratio", provenance = "user")
    })

    observeEvent(input$sf_estimate, {
      req(state$working)
      cfg <- pending_config()
      if (identical(cfg$control, "custom") && !length(cfg$custom_ids)) {
        showNotification("Enter at least one control gene for a custom set.", type = "error")
        return()
      }
      cur <- sizefactor_config(state$working)
      same <- identical(cfg$control, cur$control) &&
              setequal(cfg$custom_ids, cur$custom_ids) &&
              identical(cfg$type, cur$type)
      if (same && !is.null(tryCatch(DESeq2::sizeFactors(state$working), error = function(e) NULL))) {
        showNotification("Size factors already reflect this configuration.", type = "message")
        return()
      }
      ok <- tryCatch({
        state_mutate(state, function(d) estimate_size_factors(d, cfg),
                     action = list(action = "size_factors", control = cfg$control,
                                   type = cfg$type, n_control = length(cfg$custom_ids)))
        TRUE
      }, error = function(e) { showNotification(conditionMessage(e), type = "error"); FALSE })
      if (isTRUE(ok)) showNotification("Size factors updated.", type = "message")
    })

    output$status <- renderUI({
      req(state$working)
      dds <- state$working
      sf <- tryCatch(DESeq2::sizeFactors(dds), error = function(e) NULL)
      if (is.null(sf))
        return(tags$p(class = "text-warning", "Size factors are not set for this dataset."))
      cfg <- sizefactor_config(dds)
      n_ctrl <- length(.sf_control_index(dds, cfg))
      prov <- switch(cfg$provenance %||% "auto",
                     loaded = "provided by the loaded object",
                     user = "set on this tab", "endogenous default (computed on load)")
      ctrl_lab <- switch(cfg$control,
                         endogenous = "endogenous genes", spike_in = "spike-in genes",
                         custom = "a custom gene set", cfg$control)
      tags$p(class = "text-muted mb-2", sprintf(
        "Set from %s using %s (%d control genes), estimator '%s'. Range %.3f-%.3f.",
        ctrl_lab, prov, n_ctrl, cfg$type %||% "ratio", min(sf), max(sf)))
    })

    output$table <- DT::renderDT({
      req(state$working)
      sf <- tryCatch(DESeq2::sizeFactors(state$working), error = function(e) NULL)
      validate(need(!is.null(sf), "No size factors to show."))
      df <- data.frame(sample = names(sf) %||% as.character(seq_along(sf)),
                       size_factor = round(as.numeric(sf), 4),
                       stringsAsFactors = FALSE)
      dt_table(df, rownames = FALSE)
    })

    invisible(NULL)
  })
}
