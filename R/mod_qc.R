# Page 2: Quality control & filtering. A sub-tabbed page (like Input):
#   - "Dataset diagnostics" - VST mean-SD + sample correlation heatmap (P3b; stub).
#   - "Sample QC"           - per-sample QC metric Plot / Table sub-tabs (P3a).
#   - "Filtering"           - feature/sample filtering (P3c; later).
# Metrics are pure (R/qc_helpers.R) and computed on raw counts. Removing samples
# or features later must invalidate downstream assays/derived (rnaseq-bioc); the
# metric table is cached via state_derive() keyed on data_version.

# Metric -> axis label. Library size is reported in millions on its axis.
.qc_metric_labels <- c(
  library_size = "Library size (millions)",
  detected     = "Detected features",
  pct_mito     = "% mitochondrial",
  pct_spike    = "% spike-in"
)

# Pull a metric column as an axis (value + label); library size -> millions.
.qc_axis <- function(tbl, name) {
  v <- tbl[[name]]
  if (identical(name, "library_size")) {
    list(v = v / 1e6, lab = .qc_metric_labels[["library_size"]])
  } else {
    list(v = v, lab = .qc_metric_labels[[name]] %||% name)
  }
}

# Plot theme. thematic recolors fg/bg/accent to follow the live bslib theme
# (incl. dark mode) at draw time; `dark_theme` is the explicit lever for element
# choices thematic does not manage (here, gridline contrast on a dark panel).
.qc_theme <- function(dark_theme = FALSE) {
  base <- ggplot2::theme_minimal(base_size = 13)
  grid <- if (isTRUE(dark_theme)) "grey35" else "grey85"
  base + ggplot2::theme(
    panel.grid.major = ggplot2::element_line(colour = grid),
    panel.grid.minor = ggplot2::element_blank(),
    legend.position  = "bottom"
  )
}

# Build the QC plot. x_var = "sample" (discrete bar of the metric per sample,
# optionally sorted by value) or another metric name (numeric scatter, metric vs
# metric). Colour/fill by `group`.
.qc_metric_plot <- function(tbl, x_var, metric, group_lab = NULL,
                            sort = "none", dark_theme = FALSE) {
  yy <- .qc_axis(tbl, metric)
  if (identical(x_var, "sample")) {
    lvls <- if (identical(sort, "none")) {
      tbl$sample
    } else {
      tbl$sample[order(yy$v, decreasing = identical(sort, "decreasing"))]
    }
    df <- data.frame(x = factor(tbl$sample, levels = lvls), y = yy$v, group = tbl$group)
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = .data$y, fill = .data$group)) +
      ggplot2::geom_col() +
      ggplot2::labs(x = "sample", y = yy$lab, fill = group_lab %||% "group") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  } else {
    xx <- .qc_axis(tbl, x_var)
    df <- data.frame(x = xx$v, y = yy$v, group = tbl$group)
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = .data$y, colour = .data$group)) +
      ggplot2::geom_point(size = 2) +
      ggplot2::labs(x = xx$lab, y = yy$lab, colour = group_lab %||% "group")
  }
  p + .qc_theme(dark_theme)
}

# Choices shared by the X-axis and metric selectors.
.qc_metric_choices <- c("Library size" = "library_size",
                        "Detected features" = "detected",
                        "% mitochondrial" = "pct_mito",
                        "% spike-in" = "pct_spike")

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
      bslib::navset_pill(
        bslib::nav_panel(
          "Plot",
          bslib::layout_sidebar(
            sidebar = bslib::sidebar(
              title = tags$h4("Per-sample metrics", class = "fs-6 mb-0"), width = 280,
              uiOutput(ns("dtype_badge")),
              selectInput(ns("x_axis"), "X-axis",
                          choices = c("Sample" = "sample", .qc_metric_choices),
                          selected = "sample"),
              selectInput(ns("metric"), "QC metric (Y-axis)",
                          choices = .qc_metric_choices, selected = "library_size"),
              uiOutput(ns("group_ui")),
              conditionalPanel(
                "input.x_axis == 'sample'", ns = ns,
                selectInput(ns("sort"), "Sort by (discrete X only)",
                            choices = c("None" = "none", "Decreasing" = "decreasing",
                                        "Increasing" = "increasing"),
                            selected = "none")
              ),
              uiOutput(ns("auto_ui")),
              actionButton(ns("render"), "Render", class = "btn-primary")
            ),
            plotOutput(ns("plot"))
          )
        ),
        bslib::nav_panel(
          "Table",
          tags$small(class = "text-muted", "Per-sample QC metrics (raw counts):"),
          DT::DTOutput(ns("tbl"))
        )
      )
    )
  )
}

#' @param state the shared app-state object (see [new_app_state()]).
#' @param dark_mode A reactive returning `TRUE` when the app is in dark mode
#'   (wired from the navbar `input_dark_mode`); drives the plot theme.
#' @return Invisible NULL.
mod_qc_server <- function(id, state, dark_mode = reactive(FALSE)) {
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
      selectInput(ns("group"), "Group / colour by", choices = cols, selected = sel)
    })

    # Auto-render default depends on dataset size: live for small, button-gated big.
    output$auto_ui <- renderUI({
      req(state$working)
      checkboxInput(ns("auto"), "Auto-render", value = ncol(state$working) <= 30L)
    })

    # The ready-to-plot frame for the current selectors (metrics + chosen grouping).
    current_spec <- reactive({
      req(state$working, input$x_axis, input$metric)
      tbl <- qc_tbl()
      grp <- input$group
      cd <- as.data.frame(SummarizedExperiment::colData(state$working))
      tbl$group <- if (!is.null(grp) && grp %in% colnames(cd)) {
        as.factor(cd[tbl$sample, grp])
      } else {
        factor(rep("all", nrow(tbl)))
      }
      list(tbl = tbl, x_axis = input$x_axis, metric = input$metric,
           group_lab = grp, sort = input$sort %||% "none")
    })

    # What the plot draws: updated live when auto-render is on, else only on the
    # Render button. One reactiveVal so the plot code is not duplicated.
    shown <- reactiveVal(NULL)
    observe({ if (isTRUE(input$auto)) shown(current_spec()) })
    observeEvent(input$render, shown(current_spec()))

    output$tbl <- DT::renderDT({
      validate(need(!is.null(state$working), "No dataset loaded."))
      df <- qc_tbl()
      disp <- data.frame(
        Sample             = df$sample,
        `Library size (M)` = round(df$library_size / 1e6, 3),
        Detected           = df$detected,
        `% mito`           = round(df$pct_mito, 2),
        `% spike`          = round(df$pct_spike, 2),
        check.names = FALSE
      )
      DT::datatable(disp, rownames = FALSE, options = list(pageLength = 10, dom = "tp"))
    })

    # Reading dark_mode() here (not in `shown`) makes the plot re-theme immediately
    # on a light/dark toggle, while data/encoding still follow the auto/button rule.
    output$plot <- renderPlot({
      validate(need(!is.null(shown()),
                    "Click Render (or enable auto-render) to draw the plot."))
      spec <- shown()
      .qc_metric_plot(spec$tbl, spec$x_axis, spec$metric, spec$group_lab,
                      sort = spec$sort, dark_theme = isTRUE(dark_mode()))
    })

    invisible(NULL)
  })
}
