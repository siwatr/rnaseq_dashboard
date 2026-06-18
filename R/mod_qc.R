# Page 2: Quality control & filtering. A sub-tabbed page (navset_card_tab):
#   - "Dataset diagnostics" - VST mean-SD (the one dataset-level variance check).
#   - "Sample QC"           - per-sample pills: General QC / RLE / Expression
#                             density / QC Matrix.
#   - "Sample Correlation"  - pills: Heatmap / Within-group correlation.
#   - "Filtering"           - feature/sample filtering (P3c; later).
# Diagnostics are pure (R/qc_helpers.R). Removing samples/features later must
# invalidate downstream derived (rnaseq-bioc); heavy artifacts (metric table,
# VST, correlation) are cached via state_derive() keyed on data_version. Each
# diagnostic carries a "How to read this" note below it (.qc_diag_help).

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
  grid <- if (isTRUE(dark_theme)) "grey35" else "grey85"
  text <- if (isTRUE(dark_theme)) "grey90" else "grey25"
  ggplot2::theme() +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_line(colour = grid),
      panel.grid.minor = ggplot2::element_blank(),
      text = ggplot2::element_text(size = 14, color = text),
      legend.position  = "bottom"
    )
}

# A blank plot carrying a centered message (graceful degradation / empty state).
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

# ---- Dataset-level diagnostic plot builders --------------------------------

# Per-point 2D kernel density (the colouring geom_pointdensity provides). Uses
# MASS::kde2d (a base "recommended" package) so no extra dependency is needed;
# returns NA for every point when it cannot be computed.
.point_density <- function(x, y, n = 100L) {
  d <- rep(NA_real_, length(x))
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3L || !requireNamespace("MASS", quietly = TRUE)) return(d)
  k <- tryCatch(MASS::kde2d(x[ok], y[ok], n = n), error = function(e) NULL)
  if (is.null(k)) return(d)
  ix <- pmin(pmax(findInterval(x[ok], k$x), 1L), n)
  iy <- pmin(pmax(findInterval(y[ok], k$y), 1L), n)
  d[ok] <- k$z[cbind(ix, iy)]
  d
}

# meanSdPlot-style running median of sd across rank bins (fast, dependency-free).
.qc_meansd_trend <- function(rank, sd, nbins = 50L) {
  if (length(rank) < 2L) return(data.frame(rank = rank, sd = sd))
  br <- unique(stats::quantile(rank, seq(0, 1, length.out = nbins + 1L), na.rm = TRUE))
  if (length(br) < 2L) return(data.frame(rank = rank, sd = sd))
  bin <- cut(rank, br, include.lowest = TRUE)
  data.frame(rank = tapply(rank, bin, mean), sd = tapply(sd, bin, stats::median))
}

# Mean-SD diagnostic on a VST matrix: rank(mean) vs sd, points coloured by local
# density (point-density style) with a running-median trend line. A flat trend
# means variance is well stabilized.
.qc_meansd_plot <- function(vst_mat, dark_theme = FALSE) {
  rmean <- rowMeans(vst_mat)
  rsd <- if (requireNamespace("matrixStats", quietly = TRUE)) {
    matrixStats::rowSds(vst_mat)
  } else {
    apply(vst_mat, 1L, stats::sd)
  }
  df <- data.frame(rank = rank(rmean, ties.method = "first"), sd = rsd)
  df$density <- .point_density(df$rank, df$sd)
  trend <- .qc_meansd_trend(df$rank, df$sd)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$rank, y = .data$sd))
  p <- if (all(is.na(df$density))) {
    p + ggplot2::geom_point(size = 0.6, alpha = 0.5)
  } else {
    p + ggplot2::geom_point(ggplot2::aes(colour = .data$density), size = 0.6) +
      ggplot2::scale_color_viridis_c(name = "point density")
  }
  p +
    ggplot2::geom_line(data = trend, ggplot2::aes(x = .data$rank, y = .data$sd),
                       colour = "red", linewidth = 0.7) +
    ggplot2::labs(x = "rank (mean expression)", y = "standard deviation") +
    .qc_theme(dark_theme)
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

# Within-group correlation: boxplot per group + sample points. df from
# qc_within_group_correlation() (sample, group, mean_corr). Singleton groups
# (NA mean_corr) are dropped.
.qc_within_group_plot <- function(df, dark_theme = FALSE) {
  d <- df[!is.na(df$mean_corr), , drop = FALSE]
  if (!nrow(d)) return(.qc_msg_plot("No multi-sample groups to summarize."))
  ggplot2::ggplot(d, ggplot2::aes(x = .data$group, y = .data$mean_corr)) +
    ggplot2::geom_boxplot(ggplot2::aes(fill = .data$group), alpha = 0.35,
                          outlier.shape = NA) +
    ggplot2::geom_jitter(ggplot2::aes(colour = .data$group),
                         width = 0.15, height = 0, size = 2) +
    ggplot2::labs(x = "group", y = "mean within-group correlation",
                  fill = "group", colour = "group") +
    .qc_theme(dark_theme)
}

# Sample-sample correlation heatmap (ComplexHeatmap). Returns a Heatmap, or NULL
# when ComplexHeatmap is unavailable (the module shows a message instead).
# `anno_df` is a samples-by-columns data.frame of metadata to annotate (or NULL);
# colours come from qc_annotation_colors() so they stay stable across re-renders.
# `method` is shown in the legend title + plot title. Placeholder body palette -
# tune via the Themer "Heatmap" sub-tab.
.qc_correlation_heatmap <- function(cor_mat, anno_df = NULL, dark_theme = FALSE,
                                    n_samples = ncol(cor_mat), method = "spearman") {
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) return(NULL)
  fg <- if (isTRUE(dark_theme)) "grey90" else "grey10"
  legend_lab <- if (identical(method, "pearson")) "Pearson r" else "Spearman rho"
  title_lab  <- if (identical(method, "pearson")) "Pearson" else "Spearman"
  rng <- range(cor_mat, na.rm = TRUE)
  col_fun <- if (requireNamespace("circlize", quietly = TRUE)) {
    circlize::colorRamp2(c(rng[1], mean(rng), rng[2]),
                         c("#fff7fb", "#74a9cf", "#023858"))
  } else {
    NULL
  }
  top <- NULL
  if (!is.null(anno_df) && ncol(as.data.frame(anno_df)) > 0) {
    adf <- as.data.frame(anno_df)[colnames(cor_mat), , drop = FALSE]
    top <- ComplexHeatmap::HeatmapAnnotation(df = adf, col = qc_annotation_colors(adf))
  }
  show_names <- n_samples <= 30
  gp <- grid::gpar(col = fg)
  ComplexHeatmap::Heatmap(
    cor_mat, name = legend_lab, col = col_fun, top_annotation = top,
    show_row_names = show_names, show_column_names = show_names,
    row_names_gp = gp, column_names_gp = gp,
    column_title = paste0("Sample-to-sample correlation (", title_lab, ")"),
    column_title_gp = grid::gpar(col = fg, fontsize = 13)
  )
}

# ---- Diagnostic guideline text ---------------------------------------------

.qc_diag_help <- list(
  meansd = paste(
    "The variance-stabilizing transform should make a gene's variance",
    "independent of its mean. Each point is a gene (x = rank of mean expression,",
    "y = standard deviation), coloured by local point density; the red line is",
    "the running median. A roughly flat line means variance is well stabilized;",
    "a line that rises at the left (low-expression genes) shows a remaining trend."),
  rle = paste(
    "Relative Log Expression centers each gene on its median across samples,",
    "then box-plots the residuals per sample. Well-normalized samples have",
    "boxes centered on zero with similar, small spread. A box shifted off zero",
    "or much wider than the rest flags technical bias or a normalization issue."),
  density = paste(
    "Each curve is the distribution of log2 expression for one sample. Samples",
    "from the same protocol should overlap closely; a curve shifted or oddly",
    "shaped flags a problem sample. The tall peak at low values is",
    "near-zero genes - the motivation for low-count filtering."),
  correlation = paste(
    "Pairwise correlation of samples on log expression, clustered so similar",
    "samples sit together. Replicates of a group should correlate highly and",
    "form a block on the diagonal. A sample that correlates poorly with",
    "everything (even its own group) is a candidate outlier."),
  within_group = paste(
    "For each sample, the mean correlation to the other samples in the same",
    "group (by the chosen column). Points are samples; one sitting clearly",
    "below its group's box correlates poorly with its replicates and is a",
    "candidate outlier. Single-sample groups are omitted.")
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

# Spinner-wrapped plot output (busy indicator for slow renders).
.qc_plot <- function(id) shinycssloaders::withSpinner(plotOutput(id), proxy.height = "300px")

# Choices shared by the X-axis and metric selectors.
.qc_metric_choices <- c("Library size" = "library_size",
                        "Detected features" = "detected",
                        "% mitochondrial" = "pct_mito",
                        "% spike-in" = "pct_spike")

mod_qc_ui <- function(id) {
  ns <- NS(id)
  bslib::navset_card_tab(
    title = tags$h3("QC & filtering", class = "fs-6 mb-0 pe-3"),

    # ---- Dataset diagnostics (Mean-SD only) ----------------------------------
    bslib::nav_panel(
      tags$h4("Dataset diagnostics", class = "fs-6"),
      bslib::layout_sidebar(
        sidebar = bslib::sidebar(
          title = tags$h4("Dataset diagnostics", class = "fs-6 mb-0"), width = 280,
          uiOutput(ns("diag_badge")),
          helpText("Variance-stabilization check across the whole dataset."),
          uiOutput(ns("diag_auto_ui")),
          actionButton(ns("diag_render"), "Render", class = "btn-primary")
        ),
        .qc_plot(ns("diag_meansd")), .qc_help_note("meansd")
      )
    ),

    # ---- Sample QC -----------------------------------------------------------
    bslib::nav_panel(
      tags$h4("Sample QC", class = "fs-6"),
      bslib::navset_pill(
        bslib::nav_panel(
          "General QC",
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
            .qc_plot(ns("plot"))
          )
        ),
        bslib::nav_panel(
          "RLE",
          bslib::layout_sidebar(
            sidebar = bslib::sidebar(
              title = tags$h4("RLE", class = "fs-6 mb-0"), width = 280,
              uiOutput(ns("rle_group_ui")),
              uiOutput(ns("rle_auto_ui")),
              actionButton(ns("rle_render"), "Render", class = "btn-primary")
            ),
            .qc_plot(ns("diag_rle")), .qc_help_note("rle")
          )
        ),
        bslib::nav_panel(
          "Expression density",
          bslib::layout_sidebar(
            sidebar = bslib::sidebar(
              title = tags$h4("Expression density", class = "fs-6 mb-0"), width = 280,
              uiOutput(ns("dens_group_ui")),
              uiOutput(ns("dens_auto_ui")),
              actionButton(ns("dens_render"), "Render", class = "btn-primary")
            ),
            .qc_plot(ns("diag_density")), .qc_help_note("density")
          )
        ),
        bslib::nav_panel(
          "QC Matrix",
          tags$small(class = "text-muted", "Per-sample QC metrics (raw counts):"),
          shinycssloaders::withSpinner(DT::DTOutput(ns("tbl")), proxy.height = "300px")
        )
      )
    ),

    # ---- Sample Correlation --------------------------------------------------
    bslib::nav_panel(
      tags$h4("Sample Correlation", class = "fs-6"),
      bslib::navset_pill(
        bslib::nav_panel(
          "Heatmap",
          bslib::layout_sidebar(
            sidebar = bslib::sidebar(
              title = tags$h4("Correlation heatmap", class = "fs-6 mb-0"), width = 280,
              selectInput(ns("cor_method"), "Correlation method",
                          choices = c("Spearman" = "spearman", "Pearson" = "pearson"),
                          selected = "spearman"),
              uiOutput(ns("cor_anno_ui")),
              actionButton(ns("cor_clear_anno"), "Clear annotation",
                           icon = icon("arrows-rotate"),
                           class = "btn-sm btn-outline-secondary mb-2"),
              uiOutput(ns("cor_auto_ui")),
              actionButton(ns("cor_render"), "Render", class = "btn-primary")
            ),
            .qc_plot(ns("cor_plot")), .qc_help_note("correlation")
          )
        ),
        bslib::nav_panel(
          "Within-group correlation",
          bslib::layout_sidebar(
            sidebar = bslib::sidebar(
              title = tags$h4("Within-group correlation", class = "fs-6 mb-0"), width = 280,
              uiOutput(ns("wg_group_ui")),
              uiOutput(ns("wg_auto_ui")),
              actionButton(ns("wg_render"), "Render", class = "btn-primary")
            ),
            .qc_plot(ns("wg_plot")), .qc_help_note("within_group")
          )
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
    dark <- function() isTRUE(dark_mode())

    # Default grouping column: the design variable / condition, else first col.
    default_group_col <- function() {
      cols <- colnames(SummarizedExperiment::colData(state$working))
      if ("condition" %in% cols) "condition" else cols[1]
    }
    # Per-sample group lookup for the chosen column (named by sample id).
    group_lookup <- function(col) {
      cd <- as.data.frame(SummarizedExperiment::colData(state$working))
      stats::setNames(as.character(cd[[col]]), rownames(cd))
    }

    # --- DRY UI builders (data-dependent, so server-rendered) ---------------
    group_box <- function(input_id, label = "Group / colour by") renderUI({
      req(state$working)
      cols <- colnames(SummarizedExperiment::colData(state$working))
      selectInput(ns(input_id), label, choices = cols, selected = default_group_col())
    })
    auto_box <- function(input_id) renderUI({
      req(state$working)
      checkboxInput(ns(input_id), "Auto-render", value = ncol(state$working) <= 30L)
    })
    # Deferred render: update live when auto is on, else only on the button.
    deferred <- function(auto_id, render_id, spec) {
      rv <- reactiveVal(NULL)
      observe({ if (isTRUE(input[[auto_id]])) rv(spec()) })
      observeEvent(input[[render_id]], rv(spec()))
      rv
    }

    output$dtype_badge <- renderUI(.dtype_badge_ui(state))
    output$diag_badge  <- renderUI(.dtype_badge_ui(state))

    # --- Sample QC: General QC + Matrix -------------------------------------
    qc_tbl <- reactive({
      req(state$working)
      state_derive(state, "qc_metrics", params = list(),
                   expr = function() qc_per_sample_metrics(state$working))
    })

    output$group_ui <- group_box("group")
    output$auto_ui  <- auto_box("auto")

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
    gen_shown <- deferred("auto", "render", current_spec)

    output$plot <- renderPlot({
      validate(need(!is.null(gen_shown()),
                    "Click Render (or enable auto-render) to draw the plot."))
      s <- gen_shown()
      .qc_metric_plot(s$tbl, s$x_axis, s$metric, s$group_lab, sort = s$sort,
                      dark_theme = dark())
    })

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

    # --- Sample QC: RLE + Expression density --------------------------------
    output$rle_group_ui  <- group_box("rle_group")
    output$rle_auto_ui   <- auto_box("rle_auto")
    output$dens_group_ui <- group_box("dens_group")
    output$dens_auto_ui  <- auto_box("dens_auto")

    rle_spec <- reactive({
      req(state$working, input$rle_group)
      dds <- state$working
      rle <- qc_rle_matrix(dds)
      nf <- nrow(rle)
      gmap <- group_lookup(input$rle_group)
      df <- data.frame(
        sample = factor(rep(colnames(rle), each = nf), levels = colnames(rle)),
        value  = as.numeric(rle))
      df$group <- factor(gmap[as.character(df$sample)])
      list(df = df, n = ncol(dds))
    })
    rle_shown <- deferred("rle_auto", "rle_render", rle_spec)
    output$diag_rle <- renderPlot({
      validate(need(!is.null(rle_shown()), "Click Render (or enable auto-render)."))
      .qc_rle_plot(rle_shown()$df, dark(), rle_shown()$n)
    })

    dens_spec <- reactive({
      req(state$working, input$dens_group)
      dds <- state$working
      gmap <- group_lookup(input$dens_group)
      df <- qc_expression_long(dds)
      df$group <- factor(gmap[as.character(df$sample)])
      list(df = df, n = ncol(dds))
    })
    dens_shown <- deferred("dens_auto", "dens_render", dens_spec)
    output$diag_density <- renderPlot({
      validate(need(!is.null(dens_shown()), "Click Render (or enable auto-render)."))
      .qc_density_plot(dens_shown()$df, dark(), dens_shown()$n)
    })

    # --- Dataset diagnostics: Mean-SD ---------------------------------------
    output$diag_auto_ui <- auto_box("diag_auto")
    meansd_spec <- reactive({
      req(state$working)
      dds <- state$working
      vst_mat <- SummarizedExperiment::assay(
        state_derive(state, "vst", params = list(), expr = function() {
          withProgress(message = "Computing variance-stabilizing transform...",
                       value = 1, qc_vst(dds))
        }))
      list(vst = vst_mat)
    })
    meansd_shown <- deferred("diag_auto", "diag_render", meansd_spec)
    output$diag_meansd <- renderPlot({
      validate(need(!is.null(meansd_shown()), "Click Render (or enable auto-render)."))
      .qc_meansd_plot(meansd_shown()$vst, dark_theme = dark())
    })

    # --- Sample Correlation: Heatmap ----------------------------------------
    design_cols <- function() {
      cd_cols <- colnames(SummarizedExperiment::colData(state$working))
      dv <- tryCatch(all.vars(DESeq2::design(state$working)),
                     error = function(e) character(0))
      hit <- intersect(dv, cd_cols)
      if (length(hit)) hit else utils::head(cd_cols, 1)
    }
    output$cor_anno_ui <- renderUI({
      req(state$working)
      cols <- colnames(SummarizedExperiment::colData(state$working))
      selectizeInput(ns("cor_anno"), "Annotate by (one or more)", choices = cols,
                     selected = design_cols(), multiple = TRUE)
    })
    output$cor_auto_ui <- auto_box("cor_auto")
    observeEvent(input$cor_clear_anno, {
      updateSelectizeInput(session, "cor_anno", selected = character(0))
    })

    cor_spec <- reactive({
      req(state$working, input$cor_method)
      method <- input$cor_method
      cm <- state_derive(state, "sample_cor", params = list(method = method),
                         expr = function() {
                           withProgress(message = "Computing sample correlations...",
                                        value = 1,
                                        qc_sample_correlation(state$working, method = method))
                         })
      cd <- as.data.frame(SummarizedExperiment::colData(state$working))
      cols <- intersect(input$cor_anno, colnames(cd))
      anno_df <- if (!length(cols)) NULL else cd[, cols, drop = FALSE]
      list(cm = cm, anno_df = anno_df, n = ncol(state$working), method = method)
    })
    cor_shown <- deferred("cor_auto", "cor_render", cor_spec)
    output$cor_plot <- renderPlot({
      validate(need(!is.null(cor_shown()), "Click Render to draw the correlation heatmap."))
      validate(need(requireNamespace("ComplexHeatmap", quietly = TRUE),
                    "Install 'ComplexHeatmap' to show the correlation heatmap."))
      s <- cor_shown()
      ComplexHeatmap::draw(.qc_correlation_heatmap(
        s$cm, s$anno_df, dark_theme = dark(), n_samples = s$n, method = s$method))
    })

    # --- Sample Correlation: Within-group -----------------------------------
    output$wg_group_ui <- group_box("wg_group", label = "Group by")
    output$wg_auto_ui  <- auto_box("wg_auto")
    wg_spec <- reactive({
      req(state$working, input$wg_group)
      qc_within_group_correlation(state$working, group = input$wg_group)
    })
    wg_shown <- deferred("wg_auto", "wg_render", wg_spec)
    output$wg_plot <- renderPlot({
      validate(need(!is.null(wg_shown()), "Click Render (or enable auto-render)."))
      .qc_within_group_plot(wg_shown(), dark_theme = dark())
    })

    invisible(NULL)
  })
}
