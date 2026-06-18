# Page 2: Quality control & filtering. A sub-tabbed page (navset_card_tab):
#   - "Dataset diagnostics" - VST mean-SD, RLE, expression density (P3b).
#   - "Sample QC"           - per-sample QC metric Plot / Table sub-tabs (P3a).
#   - "Sample Correlation"  - sample-sample correlation heatmap (P3b).
#   - "Filtering"           - feature/sample filtering (P3c; later).
# Metrics/diagnostics are pure (R/qc_helpers.R). Removing samples or features
# later must invalidate downstream assays/derived (rnaseq-bioc); heavy artifacts
# (the metric table, VST, correlation) are cached via state_derive() keyed on
# data_version. Diagnostic plots are less self-explanatory than the Sample QC
# metrics, so each carries a "How to read this" note below it (.qc_diag_help).

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
# choices thematic does not manage (here, gridline + text contrast).
.qc_theme <- function(dark_theme = FALSE) {
  base <- ggplot2::theme() # Placeholder for now
  grid <- if (isTRUE(dark_theme)) "grey35" else "grey85"
  text <- if (isTRUE(dark_theme)) "grey90" else "grey25"
  base + ggplot2::theme(
    panel.grid.major = ggplot2::element_line(colour = grid),
    panel.grid.minor = ggplot2::element_blank(),
    text = ggplot2::element_text(size = 14, color = text),
    legend.position  = "bottom"
  )
}

# A blank plot carrying a centered message (graceful degradation when an
# optional package is missing).
.qc_msg_plot <- function(msg) {
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0, y = 0, label = msg) +
    ggplot2::theme_void()
}

# Build the per-sample QC plot. x_var = "sample" (discrete bar of the metric per
# sample, optionally sorted by value) or another metric name (numeric scatter,
# metric vs metric). Colour/fill by `group`.
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

# ---- Dataset-level diagnostic plot builders (P3b) --------------------------

# VST mean-SD plot (vsn). vst_mat = features x samples VST matrix. meanSdPlot
# draws with geom_hex, so it needs hexbin; degrade to a message if either is
# missing. (vsn emits a benign upstream aes_string deprecation on build.)
.qc_meansd_plot <- function(vst_mat, dark_theme = FALSE) {
  if (!requireNamespace("vsn", quietly = TRUE) ||
      !requireNamespace("hexbin", quietly = TRUE)) {
    return(.qc_msg_plot("Install 'vsn' and 'hexbin' to show the mean-SD plot."))
  }
  vsn::meanSdPlot(vst_mat, plot = FALSE)$gg + .qc_theme(dark_theme)
}

# Relative log expression boxplots. df: sample (factor), group, value.
.qc_rle_plot <- function(df, dark_theme = FALSE, n_samples = NULL) {
  n <- n_samples %||% length(unique(df$sample))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$sample, y = .data$value, fill = .data$group)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    ggplot2::geom_boxplot(outlier.size = 0.4) +
    ggplot2::labs(x = "sample", y = "relative log expression", fill = "group") +
    .qc_theme(dark_theme)
  if (n > 30) {
    p + ggplot2::theme(axis.text.x = ggplot2::element_blank())
  } else {
    p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  }
}

# Per-sample log-expression density. df: sample (factor), group, value.
.qc_density_plot <- function(df, dark_theme = FALSE, n_samples = NULL) {
  n <- n_samples %||% length(unique(df$sample))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$value, group = .data$sample,
                                        colour = .data$group)) +
    ggplot2::geom_density() +
    ggplot2::labs(x = "log2 expression", y = "density", colour = "group") +
    .qc_theme(dark_theme)
  if (n > 30) p + ggplot2::theme(legend.position = "none") else p
}

# Sample-sample correlation heatmap (ComplexHeatmap). Returns a Heatmap, or NULL
# when ComplexHeatmap is unavailable (the module shows a message instead).
# Placeholder palette - tune candidates via the Themer "Heatmap" sub-tab.
.qc_correlation_heatmap <- function(cor_mat, anno = NULL, anno_lab = NULL,
                                    dark_theme = FALSE, n_samples = ncol(cor_mat)) {
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) return(NULL)
  fg <- if (isTRUE(dark_theme)) "grey90" else "grey10"
  rng <- range(cor_mat, na.rm = TRUE)
  col_fun <- if (requireNamespace("circlize", quietly = TRUE)) {
    circlize::colorRamp2(c(rng[1], mean(rng), rng[2]),
                         c("#fff7fb", "#74a9cf", "#023858"))
  } else {
    NULL
  }
  top <- NULL
  if (!is.null(anno)) {
    adf <- data.frame(x = factor(anno[colnames(cor_mat)]))
    names(adf) <- anno_lab %||% "group"
    top <- ComplexHeatmap::HeatmapAnnotation(df = adf)
  }
  show_names <- n_samples <= 30
  gp <- grid::gpar(col = fg)
  ComplexHeatmap::Heatmap(
    cor_mat, name = "correlation", col = col_fun, top_annotation = top,
    show_row_names = show_names, show_column_names = show_names,
    row_names_gp = gp, column_names_gp = gp,
    column_title = "Sample-to-sample correlation",
    column_title_gp = grid::gpar(col = fg, fontsize = 13)
  )
}

# ---- Diagnostic guideline text ---------------------------------------------

.qc_diag_help <- list(
  meansd = paste(
    "The variance-stabilizing transform (VST) should remove the dependence of",
    "variance on mean. Genes are ranked by mean expression (x) and plotted",
    "against their standard deviation (y); the red line tracks a running sd.",
    "A roughly flat line means variance is well stabilized. A line that rises",
    "at the left (low-expression genes) shows remaining mean-variance trend."),
  rle = paste(
    "Relative Log Expression centers each gene on its median across samples,",
    "then box-plots the residuals per sample. Well-normalized samples have",
    "boxes centered on zero with similar, small spread. A box shifted off zero",
    "or much wider than the rest flags technical bias or a normalization issue."),
  density = paste(
    "Each curve is the distribution of log2 expression for one sample. Samples",
    "from the same protocol should overlap closely; a curve shifted or oddly",
    "shaped flags a problem sample. The tall peak at low values is",
    "undetected/near-zero genes - the motivation for low-count filtering."),
  correlation = paste(
    "Pairwise correlation of samples on log expression, clustered so similar",
    "samples sit together. Replicates of a group should correlate highly and",
    "form a block on the diagonal. A sample that correlates poorly with",
    "everything (even its own group) is a candidate outlier.")
)

# A subtle "How to read this" note placed below a diagnostic plot.
.qc_help_note <- function(key) {
  tags$div(
    class = "mt-2 p-2 rounded bg-body-tertiary border small",
    tags$strong("How to read this: "),
    tags$span(class = "text-muted", .qc_diag_help[[key]])
  )
}

# The data-type badge (bulk / single-cell + per-sample/per-cell unit). Bound to
# more than one output, so it is a plain function not an output expression.
.dtype_badge_ui <- function(state) {
  m <- state_meta(state)
  if (!isTRUE(m$loaded)) return(.badge("no dataset loaded", "text-bg-light"))
  unit <- if (identical(m$data_type, "single-cell")) "per cell" else "per sample"
  tags$div(class = "d-flex gap-1 align-items-center mb-2",
           .badge(m$data_type, "text-bg-info"), .badge(unit))
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

    # ---- Dataset diagnostics -------------------------------------------------
    bslib::nav_panel(
      tags$h4("Dataset diagnostics", class = "fs-6"),
      bslib::layout_sidebar(
        sidebar = bslib::sidebar(
          title = tags$h4("Dataset diagnostics", class = "fs-6 mb-0"), width = 280,
          uiOutput(ns("diag_badge")),
          helpText("Variance and distribution checks across the whole dataset."),
          uiOutput(ns("diag_auto_ui")),
          actionButton(ns("diag_render"), "Render diagnostics", class = "btn-primary")
        ),
        bslib::navset_pill(
          bslib::nav_panel("Mean-SD",
                           plotOutput(ns("diag_meansd")), .qc_help_note("meansd")),
          bslib::nav_panel("RLE",
                           plotOutput(ns("diag_rle")), .qc_help_note("rle")),
          bslib::nav_panel("Expression density",
                           plotOutput(ns("diag_density")), .qc_help_note("density"))
        )
      )
    ),

    # ---- Sample QC (P3a) -----------------------------------------------------
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
    ),

    # ---- Sample Correlation --------------------------------------------------
    bslib::nav_panel(
      tags$h4("Sample Correlation", class = "fs-6"),
      bslib::layout_sidebar(
        sidebar = bslib::sidebar(
          title = tags$h4("Sample correlation", class = "fs-6 mb-0"), width = 280,
          selectInput(ns("cor_method"), "Correlation method",
                      choices = c("Spearman" = "spearman", "Pearson" = "pearson"),
                      selected = "spearman"),
          uiOutput(ns("cor_anno_ui")),
          uiOutput(ns("cor_auto_ui")),
          actionButton(ns("cor_render"), "Render", class = "btn-primary")
        ),
        plotOutput(ns("cor_plot")),
        .qc_help_note("correlation")
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

    # Default grouping column: the design variable / condition, else first col.
    default_group_col <- function() {
      cols <- colnames(SummarizedExperiment::colData(state$working))
      if ("condition" %in% cols) "condition" else cols[1]
    }

    # ---- Sample QC (P3a) ---------------------------------------------------
    qc_tbl <- reactive({
      req(state$working)
      state_derive(state, "qc_metrics", params = list(),
                   expr = function() qc_per_sample_metrics(state$working))
    })

    output$dtype_badge <- renderUI(.dtype_badge_ui(state))
    output$diag_badge  <- renderUI(.dtype_badge_ui(state))

    output$group_ui <- renderUI({
      req(state$working)
      cols <- colnames(SummarizedExperiment::colData(state$working))
      selectInput(ns("group"), "Group / colour by", choices = cols,
                  selected = default_group_col())
    })

    # Auto-render default depends on dataset size: live for small, button-gated big.
    output$auto_ui <- renderUI({
      req(state$working)
      checkboxInput(ns("auto"), "Auto-render", value = ncol(state$working) <= 30L)
    })

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
      dt_table(disp)
    })

    output$plot <- renderPlot({
      validate(need(!is.null(shown()),
                    "Click Render (or enable auto-render) to draw the plot."))
      spec <- shown()
      .qc_metric_plot(spec$tbl, spec$x_axis, spec$metric, spec$group_lab,
                      sort = spec$sort, dark_theme = isTRUE(dark_mode()))
    })

    # ---- Dataset diagnostics (P3b) -----------------------------------------
    output$diag_auto_ui <- renderUI({
      req(state$working)
      checkboxInput(ns("diag_auto"), "Auto-render", value = ncol(state$working) <= 30L)
    })

    # Snapshot all three diagnostics' data at trigger time (VST is cached). The
    # group vector is attached for RLE/density colouring.
    diag_data <- reactive({
      req(state$working)
      dds <- state$working
      vst_mat <- SummarizedExperiment::assay(
        state_derive(state, "vst", params = list(), expr = function() qc_vst(dds)))
      cd <- as.data.frame(SummarizedExperiment::colData(dds))
      grp <- factor(cd[[default_group_col()]])
      n_feat <- nrow(dds)
      rle <- qc_rle_matrix(dds)
      rle_long <- data.frame(
        sample = factor(rep(colnames(dds), each = n_feat), levels = colnames(dds)),
        group  = rep(grp, each = n_feat),
        value  = as.numeric(rle))
      expr_long <- qc_expression_long(dds)
      expr_long$group <- rep(grp, each = n_feat)
      list(vst = vst_mat, rle = rle_long, expr = expr_long, n = ncol(dds))
    })

    diag_shown <- reactiveVal(NULL)
    observe({ if (isTRUE(input$diag_auto)) diag_shown(diag_data()) })
    observeEvent(input$diag_render, diag_shown(diag_data()))

    diag_guard <- function() {
      validate(need(!is.null(diag_shown()),
                    "Click 'Render diagnostics' (or enable auto-render)."))
    }
    output$diag_meansd <- renderPlot({
      diag_guard()
      .qc_meansd_plot(diag_shown()$vst, dark_theme = isTRUE(dark_mode()))
    })
    output$diag_rle <- renderPlot({
      diag_guard()
      .qc_rle_plot(diag_shown()$rle, isTRUE(dark_mode()), diag_shown()$n)
    })
    output$diag_density <- renderPlot({
      diag_guard()
      .qc_density_plot(diag_shown()$expr, isTRUE(dark_mode()), diag_shown()$n)
    })

    # ---- Sample Correlation (P3b) ------------------------------------------
    output$cor_anno_ui <- renderUI({
      req(state$working)
      cols <- colnames(SummarizedExperiment::colData(state$working))
      selectInput(ns("cor_anno"), "Annotate by", choices = cols,
                  selected = default_group_col())
    })
    output$cor_auto_ui <- renderUI({
      req(state$working)
      checkboxInput(ns("cor_auto"), "Auto-render", value = ncol(state$working) <= 30L)
    })

    cor_data <- reactive({
      req(state$working, input$cor_method)
      method <- input$cor_method
      cm <- state_derive(state, "sample_cor", params = list(method = method),
                         expr = function() qc_sample_correlation(state$working, method = method))
      cd <- as.data.frame(SummarizedExperiment::colData(state$working))
      anno_col <- input$cor_anno
      anno <- if (!is.null(anno_col) && anno_col %in% colnames(cd)) {
        stats::setNames(as.character(cd[[anno_col]]), rownames(cd))
      } else {
        NULL
      }
      list(cm = cm, anno = anno, anno_lab = anno_col, n = ncol(state$working))
    })

    cor_shown <- reactiveVal(NULL)
    observe({ if (isTRUE(input$cor_auto)) cor_shown(cor_data()) })
    observeEvent(input$cor_render, cor_shown(cor_data()))

    output$cor_plot <- renderPlot({
      validate(need(!is.null(cor_shown()), "Click Render to draw the correlation heatmap."))
      validate(need(requireNamespace("ComplexHeatmap", quietly = TRUE),
                    "Install 'ComplexHeatmap' to show the correlation heatmap."))
      s <- cor_shown()
      ht <- .qc_correlation_heatmap(s$cm, s$anno, s$anno_lab,
                                    dark_theme = isTRUE(dark_mode()), n_samples = s$n)
      ComplexHeatmap::draw(ht)
    })

    invisible(NULL)
  })
}
