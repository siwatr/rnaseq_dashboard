# Page 2: Quality control & filtering. A sub-tabbed page (like Input):
#   - "Dataset diagnostics" - VST mean-SD + sample correlation heatmap (P3b; stub).
#   - "Sample QC"           - per-sample QC metric table + plot (this slice, P3a).
#   - "Filtering"           - feature/sample filtering (P3c; later).
# Metrics are pure (R/qc_helpers.R) and computed on raw counts. Removing samples
# or features later must invalidate downstream assays/derived (rnaseq-bioc); the
# metric table is cached via state_derive() keyed on data_version.

# Plot-able metrics: input value -> axis label (library size is shown in millions).
.qc_metric_labels <- c(
  library_size = "Library size (millions)",
  detected     = "Detected features",
  pct_mito     = "% mitochondrial",
  pct_spike    = "% spike-in"
)

# Build the per-metric ggplot from a metrics frame carrying a `group` column.
.qc_metric_plot <- function(tbl, metric, plot_type, group_lab = NULL) {
  y <- if (identical(metric, "library_size")) tbl[[metric]] / 1e6 else tbl[[metric]]
  df <- data.frame(sample = factor(tbl$sample, levels = tbl$sample),
                   group = tbl$group, y = y)
  ylab <- .qc_metric_labels[[metric]] %||% metric
  p <- switch(plot_type,
    box = ggplot2::ggplot(df, ggplot2::aes(x = .data$group, y = .data$y, fill = .data$group)) +
      ggplot2::geom_boxplot(outlier.size = 0.8) +
      ggplot2::labs(x = group_lab %||% "group", y = ylab),
    scatter = ggplot2::ggplot(df, ggplot2::aes(x = .data$sample, y = .data$y, colour = .data$group)) +
      ggplot2::geom_point(size = 2) +
      ggplot2::labs(x = "sample", y = ylab, colour = group_lab %||% "group"),
    # default: bar
    ggplot2::ggplot(df, ggplot2::aes(x = .data$sample, y = .data$y, fill = .data$group)) +
      ggplot2::geom_col() +
      ggplot2::labs(x = "sample", y = ylab, fill = group_lab %||% "group")
  )
  p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

mod_qc_ui <- function(id) {
  ns <- NS(id)
  bslib::navset_card_tab(
    title = tags$h3("QC & filtering", class = "fs-6 mb-0 pe-3"),
    bslib::nav_panel(
      tags$h4("Dataset diagnostics", class = "fs-6"),
      bslib::card_body(
        tags$p(class = "text-muted",
               "Variance-stabilization (mean-SD) and sample-correlation diagnostics arrive in the next slice.")
      )
    ),
    bslib::nav_panel(
      tags$h4("Sample QC", class = "fs-6"),
      bslib::layout_sidebar(
        sidebar = bslib::sidebar(
          title = tags$h4("Per-sample metrics", class = "fs-6 mb-0"), width = 280,
          uiOutput(ns("dtype_badge")),
          selectInput(ns("metric"), "Metric",
                      choices = c("Library size" = "library_size",
                                  "Detected features" = "detected",
                                  "% mitochondrial" = "pct_mito",
                                  "% spike-in" = "pct_spike")),
          selectInput(ns("plot_type"), "Plot type",
                      choices = c("Bar" = "bar", "Box" = "box", "Scatter" = "scatter")),
          uiOutput(ns("group_ui")),
          uiOutput(ns("auto_ui")),
          actionButton(ns("render"), "Render", class = "btn-primary")
        ),
        DT::DTOutput(ns("tbl")),
        plotOutput(ns("plot"))
      )
    )
  )
}

#' @param state the shared app-state object (see [new_app_state()]).
#' @return Invisible NULL.
mod_qc_server <- function(id, state) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Version-stamped per-sample metric table; recomputed when data_version bumps.
    qc_tbl <- reactive({
      req(state$working)
      state_derive(state, "qc_metrics", params = list(),
                   expr = function() qc_per_sample_metrics(state$working))
    })

    output$dtype_badge <- renderUI({
      m <- state_meta(state)
      if (!isTRUE(m$loaded)) return(.badge("no dataset loaded", "text-bg-light"))
      unit <- if (identical(m$data_type, "single-cell")) "per cell" else "per sample"
      tags$div(class = "d-flex gap-1 align-items-center mb-2",
               .badge(m$data_type, "text-bg-info"), .badge(unit))
    })

    output$group_ui <- renderUI({
      req(state$working)
      cols <- colnames(SummarizedExperiment::colData(state$working))
      sel <- if ("condition" %in% cols) "condition" else cols[1]
      selectInput(ns("group"), "Group by", choices = cols, selected = sel)
    })

    # Auto-render default depends on dataset size: live for small, button-gated big.
    output$auto_ui <- renderUI({
      req(state$working)
      checkboxInput(ns("auto"), "Auto-render", value = ncol(state$working) <= 30L)
    })

    # The ready-to-plot frame for the current selectors (metrics + chosen grouping).
    current_spec <- reactive({
      req(state$working, input$metric, input$plot_type)
      tbl <- qc_tbl()
      grp <- input$group
      cd <- as.data.frame(SummarizedExperiment::colData(state$working))
      tbl$group <- if (!is.null(grp) && grp %in% colnames(cd)) {
        as.factor(cd[tbl$sample, grp])
      } else {
        factor(rep("all", nrow(tbl)))
      }
      list(tbl = tbl, metric = input$metric, plot_type = input$plot_type, group_lab = grp)
    })

    # What the plot actually draws: updated live when auto-render is on, else only
    # on the Render button. One reactiveVal so the plot code is not duplicated.
    shown <- reactiveVal(NULL)
    observe({ if (isTRUE(input$auto)) shown(current_spec()) })
    observeEvent(input$render, shown(current_spec()))

    output$tbl <- DT::renderDT({
      validate(need(!is.null(state$working), "No dataset loaded."))
      df <- qc_tbl()
      disp <- data.frame(
        Sample            = df$sample,
        `Library size (M)` = round(df$library_size / 1e6, 3),
        Detected          = df$detected,
        `% mito`          = round(df$pct_mito, 2),
        `% spike`         = round(df$pct_spike, 2),
        check.names = FALSE
      )
      DT::datatable(disp, rownames = FALSE, options = list(pageLength = 10, dom = "tp"))
    })

    output$plot <- renderPlot({
      validate(need(!is.null(shown()),
                    "Click Render (or enable auto-render) to draw the plot."))
      spec <- shown()
      .qc_metric_plot(spec$tbl, spec$metric, spec$plot_type, spec$group_lab)
    })

    invisible(NULL)
  })
}
