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
      # All three honor the control set: estimate_size_factors() estimates on the
      # control-gene row-subset (so DESeq2's "iterate", which ignores controlGenes,
      # still respects it) and the full dds inherits the per-sample factors.
      selectInput(ns("sf_type"), "Estimator (type)",
                  c("ratio (median-of-ratios)" = "ratio",
                    "poscounts" = "poscounts", "iterate" = "iterate"),
                  selected = "ratio"),
      conditionalPanel(
        sprintf("input['%s'] != 'endogenous'", ns("sf_control")),
        tags$div(class = "small text-warning mt-1",
                 "This normalization is used by DE, PCA, and Expression (not just QC). Spike-in factors can be noisy with few or low-count spikes.")),
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
      # Always re-estimate (it's cheap) so the user can re-run at will; commit only
      # when the values or the recorded config actually change, so a repeat of the
      # same known config doesn't needlessly invalidate downstream caches. A loaded
      # (unknown) config always commits -- the fresh estimate differs from the
      # inherited vector, and this is how the user replaces it with a known one.
      cand <- tryCatch(estimate_size_factors(state$working, cfg),
                       error = function(e) { showNotification(conditionMessage(e), type = "error"); NULL })
      if (is.null(cand)) return()
      cur <- sizefactor_config(state$working)
      cur_sf <- tryCatch(DESeq2::sizeFactors(state$working), error = function(e) NULL)
      unchanged <- identical(cur$provenance, "user") &&
                   identical(cur$control, cfg$control) && setequal(cur$custom_ids, cfg$custom_ids) &&
                   identical(cur$type, cfg$type) && !is.null(cur_sf) &&
                   isTRUE(all.equal(unname(cur_sf), unname(DESeq2::sizeFactors(cand))))
      if (unchanged) {
        showNotification("Size factors re-estimated - values unchanged.", type = "message")
        return()
      }
      state_mutate(state, function(d) cand,
                   action = list(action = "size_factors", control = cfg$control, type = cfg$type,
                                 n_control = length(.sf_control_index(state$working, cfg))))
      showNotification("Size factors updated.", type = "message")
    })

    output$status <- renderUI({
      req(state$working)
      dds <- state$working
      sf <- tryCatch(DESeq2::sizeFactors(dds), error = function(e) NULL)
      if (is.null(sf))
        return(tags$p(class = "text-warning",
                      "Size factors are not set for this dataset. Estimate them below."))
      cfg <- sizefactor_config(dds)
      rng <- sprintf("Range %.3f-%.3f.", min(sf), max(sf))
      # For a loaded object we don't know the upstream control set / estimator, so
      # don't invent them -- only report what WE computed (auto default / user set).
      txt <- switch(cfg$provenance %||% "auto",
        loaded = sprintf(
          "Inherited from the loaded object; the control genes and estimator used upstream are unknown. %s Re-estimate below to set a known configuration.", rng),
        user = {
          ctrl_lab <- switch(cfg$control, endogenous = "endogenous genes",
                             spike_in = "spike-in genes", custom = "a custom gene set", cfg$control)
          sprintf("Set on this tab from %s (%d control genes), estimator '%s'. %s",
                  ctrl_lab, length(.sf_control_index(dds, cfg)), cfg$type %||% "ratio", rng)
        },
        sprintf("Endogenous default (median-of-ratios), computed on load from %d endogenous genes. %s",
                length(.sf_control_index(dds, cfg)), rng))
      tags$p(class = "text-muted mb-2", txt)
    })

    output$table <- DT::renderDT({
      req(state$working)
      sf <- tryCatch(DESeq2::sizeFactors(state$working), error = function(e) NULL)
      validate(need(!is.null(sf), "No size factors to show."))
      df <- data.frame(sample = names(sf) %||% as.character(seq_along(sf)),
                       size_factor = round(as.numeric(sf), 4),
                       stringsAsFactors = FALSE)
      dt_table(df)
    })

    invisible(NULL)
  })
}
