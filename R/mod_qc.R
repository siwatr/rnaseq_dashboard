# Page 2: Quality control & filtering. A sub-tabbed page (navset_card_tab):
#   - "Dataset diagnostics" - VST mean-SD (the one dataset-level variance check).
#   - "Sample QC"           - per-sample pills: General QC / RLE / Expression
#                             density / QC Matrix.
#   - "Sample Correlation"  - pills: Heatmap / Within-group correlation.
#   - "Filtering"           - pills: Samples / Features. The app *suggests*
#                             low-quality items (R/filter_helpers.R); the user
#                             builds a removal pool (DT selection + buttons) and
#                             Applies it via state_mutate (real, undoable removal).
# Diagnostics are pure (R/qc_helpers.R). Removing samples/features invalidates
# downstream derived (rnaseq-bioc); heavy artifacts (metric table, VST,
# correlation, flag tables) are cached via state_derive() keyed on data_version.
# A page-level "Showing:" control hides samples from the sample plots + the
# correlation heatmap without touching the data (view-only; no data_version bump).
# Each diagnostic carries a "How to read this" note below it (.qc_diag_help).

# Semantic 3-colour scheme for the "Removal status" colour-by (reason-aware):
# green = QC pass, yellow = suggested for some other reason, red = suggested for
# the reason of the current plot. A fixed scale (not the qualitative palette);
# the Palette page will later own these values.
.removal_palette <- c(pass = "#2CA02C", suggested_other = "#E6B800",
                      suggested_this = "#D62728")
.removal_labels <- c(pass = "QC pass", suggested_other = "Suggested drop (other)",
                     suggested_this = "Suggested drop (this reason)")
# Sample QC metric -> the flag_samples() reason column it corresponds to (used to
# pick the "this reason" highlight). % spike-in has no drop reason.
.metric_reason <- c(library_size = "low_lib_size", detected = "low_detected",
                    pct_mito = "high_mito")

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
# metric vs metric). Colour/fill by `group`. When `palette` (a named colour
# vector over the group levels) is supplied, a fixed manual scale is used instead
# of the qualitative palette (for the semantic removal-status colouring).
.qc_metric_plot <- function(tbl, x_var, metric, group_lab = NULL,
                            sort = "none", dark_theme = FALSE,
                            palette = NULL, palette_labels = NULL) {
  yy <- .qc_axis(tbl, metric)
  bar <- identical(x_var, "sample")
  if (bar) {
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
  if (!is.null(palette)) {
    labs_arg <- if (is.null(palette_labels)) ggplot2::waiver() else palette_labels
    scale_fn <- if (bar) ggplot2::scale_fill_manual else ggplot2::scale_colour_manual
    p <- p + scale_fn(values = palette, labels = labs_arg, drop = FALSE,
                      name = group_lab %||% "group")
  }
  p + .qc_theme(dark_theme)
}

# Before/after feature-filter density (limma/edgeR RNAseq123 style): the
# distribution of log-expression over endogenous features, before vs after the
# proposed filter. A cleaner after-curve (less low-expression mass) means the
# filter is doing its job. df from qc_filter_density() (sample, value, status).
.qc_filter_density_plot <- function(df, dark_theme = FALSE) {
  ggplot2::ggplot(df, ggplot2::aes(x = .data$value,
                                   group = interaction(.data$sample, .data$status),
                                   colour = .data$status)) +
    ggplot2::geom_density(linewidth = 0.4, alpha = 0.7) +
    ggplot2::scale_colour_manual(values = c(before = "#9aa0a6", after = "#1f77b4"),
                                 name = NULL) +
    ggplot2::labs(x = "log2 expression", y = "density") +
    .qc_theme(dark_theme)
}

# ---- Spike-in (ERCC) dose-response builders --------------------------------

# Per-sample metric labels for the spike-in summary plot/table.
.spike_metric_labels <- c(
  pct_spike        = "% spike-in (of library)",
  n_spike_detected = "Detected spike-in features",
  lod              = "Lowest detected conc (attomoles/uL)",
  slope            = "Dose-response slope (log-log)",
  r_squared        = "Dose-response R-squared")

# Dose-response scatter: known concentration vs observed expression (log-log),
# coloured by group, with a per-sample lm line. `long` carries a `group` column.
# Zeros (undetected) are dropped before the log scales.
.qc_spike_dr_plot <- function(long, dark_theme = FALSE) {
  d <- long[is.finite(long$concentration) & is.finite(long$expression) &
            long$concentration > 0 & long$expression > 0, , drop = FALSE]
  if (!nrow(d)) return(.qc_msg_plot("No detected spike-ins with known concentration."))
  ggplot2::ggplot(d, ggplot2::aes(x = .data$concentration, y = .data$expression,
                                  colour = .data$group)) +
    ggplot2::geom_point(size = 1.4, alpha = 0.7) +
    ggplot2::geom_smooth(ggplot2::aes(group = .data$sample), method = "lm",
                         formula = y ~ x, se = FALSE, linewidth = 0.4) +
    ggplot2::scale_x_log10() + ggplot2::scale_y_log10() +
    ggplot2::labs(x = "known concentration (attomoles/uL)", y = "observed expression",
                  colour = "group") +
    .qc_theme(dark_theme)
}

# Per-sample spike summary: a chosen metric across samples (bar, sorted), filled
# by group. `df` is spike_dose_response()$per_sample with a `group` column.
.qc_spike_summary_plot <- function(df, metric, dark_theme = FALSE) {
  lab <- .spike_metric_labels[[metric]] %||% metric
  d <- df[is.finite(df[[metric]]), , drop = FALSE]
  if (!nrow(d)) return(.qc_msg_plot(paste("No", lab, "to show.")))
  d$sample <- factor(d$sample, levels = d$sample[order(d[[metric]])])
  ggplot2::ggplot(d, ggplot2::aes(x = .data$sample, y = .data[[metric]], fill = .data$group)) +
    ggplot2::geom_col() +
    ggplot2::labs(x = "sample", y = lab, fill = "group") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    .qc_theme(dark_theme)
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
    "candidate outlier. Single-sample groups are omitted."),
  filter_density = paste(
    "Distribution of log2 expression over endogenous features, before vs after",
    "the proposed filter. The tall low-value peak is near-zero genes; a good",
    "filter removes most of it, leaving a cleaner unimodal 'after' curve.",
    "Spike-in / exogenous features are exempt and never removed here."),
  spike_dr = paste(
    "Spike-in titration QC (technical control; not used for normalization here).",
    "Each point is an ERCC control: known input concentration (x) vs observed",
    "expression (y), log-log, with a per-sample fit. A healthy titration tracks a",
    "straight line of slope ~ 1 with high R-squared; samples that fan out, flatten,",
    "or lose low-concentration spike-ins (higher 'lowest detected conc') are",
    "candidates for concern. Undetected (zero) spike-ins are dropped before fitting.",
    "Slope is only meaningful when concentrations span several logs; a vertical",
    "offset between samples reflects a differing spike-in fraction (input mass),",
    "not titration quality. Zero/low spike-ins in only some samples usually just",
    "means those samples were not spiked (mixed designs).")
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

# A blank numericInput reads back as NA; map that to NULL so the flag helpers
# treat that threshold as "rule disabled" (no flag).
.blank_na <- function(x) if (is.null(x) || length(x) != 1L || is.na(x)) NULL else x

# Pool summary badges: suggested / currently-selected / remove-pool / kept-of-total.
.pool_counts <- function(suggested, selected, pooled, total) {
  tags$div(class = "small mb-2 d-flex flex-wrap gap-1 align-items-center",
    tags$span(class = "badge text-bg-warning", paste("Suggested:", suggested)),
    tags$span(class = "badge text-bg-primary", paste("Selected:", selected)),
    tags$span(class = "badge text-bg-danger", paste("Remove pool:", pooled)),
    tags$span(class = "text-muted", sprintf("Keep %d of %d", total - pooled, total)))
}

# Choices shared by the X-axis and metric selectors.
.qc_metric_choices <- c("Library size" = "library_size",
                        "Detected features" = "detected",
                        "% mitochondrial" = "pct_mito",
                        "% spike-in" = "pct_spike")

mod_qc_ui <- function(id) {
  ns <- NS(id)
  # Above-table action row: Select all / Deselect all stage the rows matching the
  # current table search; "Reset … Removal" (warning) un-filters that dimension.
  select_buttons <- function(prefix, reset_label) {
    tags$div(class = "d-flex gap-2 mb-2 align-items-center",
      actionButton(ns(paste0(prefix, "_select_all")), "Select all",
                   class = "btn-sm btn-outline-secondary"),
      actionButton(ns(paste0(prefix, "_deselect_all")), "Deselect all",
                   class = "btn-sm btn-outline-secondary"),
      actionButton(ns(paste0(prefix, "_reset")), reset_label, icon = icon("rotate-left"),
                   class = "btn-sm btn-warning ms-auto"))
  }
  # The removal-pool action buttons (sidebar): move the staged selection into or
  # out of the removal pool.
  pool_buttons <- function(prefix) {
    tags$div(
      class = "mb-2",
      tags$div(class = "fw-semibold small text-body-secondary mb-1", "Removal Pool"),
      tags$div(class = "d-grid gap-1",
        actionButton(ns(paste0(prefix, "_add")), "Add selected to pool",
                     class = "btn-sm btn-outline-primary"),
        actionButton(ns(paste0(prefix, "_remove")), "Remove selected from pool",
                     class = "btn-sm btn-outline-secondary"),
        actionButton(ns(paste0(prefix, "_adopt")), "Adopt remove suggestions",
                     icon = icon("wand-magic-sparkles"), class = "btn-sm btn-outline-primary"),
        actionButton(ns(paste0(prefix, "_clear")), "Clear pool",
                     icon = icon("arrows-rotate"), class = "btn-sm btn-outline-secondary"))
    )
  }
  # View-only "Showing:" control, repeated in every sample-plot sidebar and kept
  # in sync by the server. Default Keep is empty => show all (no subsetting).
  showing_ctrl <- function(s) tagList(
    tags$hr(class = "my-2"),
    selectInput(ns(paste0(s, "_show_by")), "Showing (display only)",
                choices = c("All samples" = "__all__"), selected = "__all__"),
    selectizeInput(ns(paste0(s, "_show_values")), "Keep (blank = show all)",
                   choices = character(0), multiple = TRUE,
                   options = list(placeholder = "(blank = show all)"))
  )
  # A threshold numericInput with its own per-field "Auto" button alongside. The
  # label sits on its own line; the input (wide) and the smaller wand button sit
  # on the next line, vertically centred so they read as aligned. `tip` is the
  # button's hover text (e.g. "Auto threshold: min. library size").
  thr_input <- function(input_id, label, auto_id, tip, ...) {
    tags$div(class = "mb-2",
      tags$label(label, class = "form-label mb-1", `for` = ns(input_id)),
      tags$div(class = "d-flex align-items-center gap-1",
        tags$div(class = "flex-grow-1", numericInput(ns(input_id), label = NULL, ...)),
        bslib::tooltip(
          actionButton(ns(auto_id), NULL, icon = icon("wand-magic-sparkles"),
                       class = "btn-sm btn-outline-primary"),
          tip)))
  }
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
        uiOutput(ns("diag_stale")), .qc_plot(ns("diag_meansd")), .qc_help_note("meansd")
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
              showing_ctrl("gen"),
              uiOutput(ns("auto_ui")),
              actionButton(ns("render"), "Render", class = "btn-primary")
            ),
            uiOutput(ns("gen_stale")), .qc_plot(ns("plot"))
          )
        ),
        bslib::nav_panel(
          "RLE",
          bslib::layout_sidebar(
            sidebar = bslib::sidebar(
              title = tags$h4("RLE", class = "fs-6 mb-0"), width = 280,
              uiOutput(ns("rle_group_ui")),
              showing_ctrl("rle"),
              uiOutput(ns("rle_auto_ui")),
              actionButton(ns("rle_render"), "Render", class = "btn-primary")
            ),
            uiOutput(ns("rle_stale")), .qc_plot(ns("diag_rle")), .qc_help_note("rle")
          )
        ),
        bslib::nav_panel(
          "Expression density",
          bslib::layout_sidebar(
            sidebar = bslib::sidebar(
              title = tags$h4("Expression density", class = "fs-6 mb-0"), width = 280,
              uiOutput(ns("dens_group_ui")),
              showing_ctrl("dens"),
              uiOutput(ns("dens_auto_ui")),
              actionButton(ns("dens_render"), "Render", class = "btn-primary")
            ),
            uiOutput(ns("dens_stale")), .qc_plot(ns("diag_density")), .qc_help_note("density")
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
              showing_ctrl("cor"),
              uiOutput(ns("cor_auto_ui")),
              actionButton(ns("cor_render"), "Render", class = "btn-primary")
            ),
            uiOutput(ns("cor_stale")), .qc_plot(ns("cor_plot")), .qc_help_note("correlation")
          )
        ),
        bslib::nav_panel(
          "Within-group correlation",
          bslib::layout_sidebar(
            sidebar = bslib::sidebar(
              title = tags$h4("Within-group correlation", class = "fs-6 mb-0"), width = 280,
              uiOutput(ns("wg_group_ui")),
              showing_ctrl("wg"),
              uiOutput(ns("wg_auto_ui")),
              actionButton(ns("wg_render"), "Render", class = "btn-primary")
            ),
            uiOutput(ns("wg_stale")), .qc_plot(ns("wg_plot")), .qc_help_note("within_group")
          )
        )
      )
    ),

    # ---- Spike-in (ERCC) -----------------------------------------------------
    bslib::nav_panel(
      tags$h4("Spike-in (ERCC)", class = "fs-6"),
      bslib::layout_sidebar(
        sidebar = bslib::sidebar(
          title = tags$h4("Spike-in titration", class = "fs-6 mb-0"), width = 300,
          helpText("Technical QC of the ERCC titration (not used for normalization)."),
          uiOutput(ns("spike_source_ui")),
          uiOutput(ns("spike_assay_ui")),
          uiOutput(ns("spike_group_ui")),
          uiOutput(ns("spike_auto_ui")),
          actionButton(ns("spike_render"), "Render", class = "btn-primary")
        ),
        uiOutput(ns("spike_msg")),
        uiOutput(ns("spike_stale")),
        bslib::navset_pill(
          bslib::nav_panel("Dose-response",
            .qc_plot(ns("spike_dr_plot")), .qc_help_note("spike_dr")),
          bslib::nav_panel("Per-sample summary",
            selectInput(ns("spike_metric"), "Metric",
                        choices = c("% spike-in" = "pct_spike",
                                    "Detected spike-in features" = "n_spike_detected",
                                    "Lowest detected conc" = "lod",
                                    "Dose-response slope" = "slope",
                                    "Dose-response R-squared" = "r_squared"),
                        selected = "r_squared"),
            uiOutput(ns("spike_cv")),
            .qc_plot(ns("spike_summary_plot")),
            shinycssloaders::withSpinner(DT::DTOutput(ns("spike_table")), proxy.height = "200px"))
        )
      )
    ),

    # ---- Filtering -----------------------------------------------------------
    bslib::nav_panel(
      tags$h4("Filtering", class = "fs-6"),
      bslib::navset_pill(
        bslib::nav_panel(
          "Samples",
          bslib::layout_sidebar(
            sidebar = bslib::sidebar(
              title = tags$h4("Sample filtering", class = "fs-6 mb-0"), width = 300,
              helpText("Flags are advisory. A blank threshold disables that check; the ",
                       tags$strong("Auto"), " buttons fill data-driven thresholds."),
              uiOutput(ns("samp_group_ui")),
              thr_input("samp_lib_min", "Min library size (blank = off)", "samp_lib_auto",
                        tip = "Auto threshold: min. library size", value = NA, min = 0),
              thr_input("samp_detected_min", "Min detected features (blank = off)",
                        "samp_detected_auto", tip = "Auto threshold: min. detected features",
                        value = NA, min = 0),
              thr_input("samp_mito_max", "Max % mitochondrial (blank = off)", "samp_mito_auto",
                        tip = "Auto threshold: max. % mitochondrial", value = NA, min = 0, max = 100),
              numericInput(ns("samp_wg_z"), "Within-group outlier z-cutoff (blank = off)",
                           value = 2, min = 0, step = 0.5),
              actionButton(ns("samp_auto"), "Set auto threshold for all settings",
                           icon = icon("wand-magic-sparkles"),
                           class = "btn-sm btn-outline-primary w-100 mb-2"),
              tags$hr(),
              pool_buttons("samp"),
              uiOutput(ns("samp_counts")),
              actionButton(ns("samp_apply"), "Remove Samples",
                           icon = icon("trash"), class = "btn-danger w-100")
            ),
            tags$small(class = "text-muted",
                       "Select rows to stage them, then use the buttons. Flagged samples are highlighted but never pre-pooled."),
            select_buttons("samp", "Reset Sample Removal"),
            shinycssloaders::withSpinner(DT::DTOutput(ns("samp_tbl")), proxy.height = "300px")
          )
        ),
        bslib::nav_panel(
          "Features",
          bslib::layout_sidebar(
            sidebar = bslib::sidebar(
              title = tags$h4("Feature filtering", class = "fs-6 mb-0"), width = 300,
              helpText("Expression filters run on endogenous features only; spike-in / exogenous are exempt."),
              checkboxInput(ns("feat_use_fbe"), "Use filterByExpr (edgeR)", value = TRUE),
              conditionalPanel(
                "input.feat_use_fbe", ns = ns,
                uiOutput(ns("feat_group_ui"))
              ),
              conditionalPanel(
                "!input.feat_use_fbe", ns = ns,
                numericInput(ns("feat_min_count"), "Min total count", value = 10, min = 0),
                checkboxInput(ns("feat_use_min_samples"), "Also require min detected samples",
                              value = FALSE),
                conditionalPanel(
                  "input.feat_use_min_samples", ns = ns,
                  numericInput(ns("feat_min_samples"), "Min detected samples", value = 2, min = 1)
                )
              ),
              tags$hr(),
              pool_buttons("feat"),
              uiOutput(ns("feat_counts")),
              actionButton(ns("feat_apply"), "Remove Features",
                           icon = icon("trash"), class = "btn-danger w-100"),
              tags$hr(),
              helpText("Drop all spike-in (ERCC) controls if your design has none / they failed."),
              actionButton(ns("feat_drop_spike"), "Remove all spike-in features",
                           icon = icon("flask"), class = "btn-outline-danger w-100")
            ),
            tags$small(class = "text-muted",
                       "The removal pool is pre-seeded with the suggestion; search the table then 'Select all' + 'Add selected to pool' for bulk edits."),
            select_buttons("feat", "Reset Feature Removal"),
            shinycssloaders::withSpinner(DT::DTOutput(ns("feat_tbl")), proxy.height = "300px"),
            .qc_plot(ns("feat_density")), .qc_help_note("filter_density")
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
    # `sig` is a cheap signature of the inputs the plot depends on; after a manual
    # render, a change to `sig` marks the plot "stale" so we can nudge the user to
    # re-render. Returns a list of `value()` and `stale()` reactives.
    deferred <- function(auto_id, render_id, spec, sig) {
      rv <- reactiveVal(NULL)
      last_sig <- reactiveVal(NULL)
      go <- function() { rv(spec()); last_sig(sig()) }
      observe({ if (isTRUE(input[[auto_id]])) go() })
      observeEvent(input[[render_id]], go())
      stale <- reactive({
        if (is.null(rv()) || isTRUE(input[[auto_id]])) return(FALSE)
        !isTRUE(all.equal(last_sig(), sig()))
      })
      list(value = rv, stale = stale)
    }
    # A "settings changed -> re-render" banner, shown above a plot when stale.
    stale_note <- function(d) renderUI({
      if (!isTRUE(d$stale())) return(NULL)
      tags$div(class = "alert alert-warning py-1 px-2 small mb-2 d-flex align-items-center gap-2",
               icon("triangle-exclamation"),
               "Settings changed - click Render to update the plot.")
    })

    output$dtype_badge <- renderUI(.dtype_badge_ui(state))
    output$diag_badge  <- renderUI(.dtype_badge_ui(state))

    # --- "Showing:" display subset (view-only; never mutates the dds) --------
    # One canonical selection (column + kept values) shared across each sample
    # plot's sidebar control. Editing any tab's control updates the canonical
    # state, which is fanned back out to every tab so they stay in sync. Default
    # Keep is empty => show all (no subsetting).
    .show_tabs <- c("gen", "rle", "dens", "wg", "cor")
    show_by_rv     <- reactiveVal("__all__")
    show_values_rv <- reactiveVal(character(0))

    # Value-box choices for the current "show by" column.
    show_val_choices <- function(by) {
      if (is.null(by) || identical(by, "__all__")) return(character(0))
      if (identical(by, "__samples__")) return(colnames(state$working))
      cd <- as.data.frame(SummarizedExperiment::colData(state$working))
      sort(unique(as.character(cd[[by]])))
    }

    # Populate / reset all the per-tab controls when a dataset (re)loads.
    observeEvent(state$working, {
      cols <- colnames(SummarizedExperiment::colData(state$working))
      ch <- c("All samples" = "__all__", stats::setNames(cols, cols),
              "Individual samples" = "__samples__")
      show_by_rv("__all__"); show_values_rv(character(0))
      for (s in .show_tabs) {
        updateSelectInput(session, paste0(s, "_show_by"), choices = ch, selected = "__all__")
        updateSelectizeInput(session, paste0(s, "_show_values"),
                             choices = character(0), selected = character(0))
      }
    })

    # Per-tab edits -> canonical state (guarded so fan-out echoes are no-ops).
    lapply(.show_tabs, function(s) {
      observeEvent(input[[paste0(s, "_show_by")]], {
        v <- input[[paste0(s, "_show_by")]]
        if (is.null(v) || identical(v, show_by_rv())) return()
        show_by_rv(v); show_values_rv(character(0))      # reset values on column switch
      }, ignoreInit = TRUE)
      observeEvent(input[[paste0(s, "_show_values")]], {
        v <- input[[paste0(s, "_show_values")]] %||% character(0)
        if (setequal(v, show_values_rv())) return()
        show_values_rv(v)
      }, ignoreNULL = FALSE, ignoreInit = TRUE)
    })

    # Canonical state -> fan out to every tab's controls (keeps them in sync).
    observeEvent(show_by_rv(), {
      ch <- show_val_choices(show_by_rv())
      for (s in .show_tabs) {
        updateSelectInput(session, paste0(s, "_show_by"), selected = show_by_rv())
        updateSelectizeInput(session, paste0(s, "_show_values"),
                             choices = ch, selected = show_values_rv())
      }
    })
    observeEvent(show_values_rv(), {
      for (s in .show_tabs) {
        updateSelectizeInput(session, paste0(s, "_show_values"), selected = show_values_rv())
      }
    }, ignoreNULL = FALSE)

    # Samples currently shown (always non-empty: a blank Keep box = show all).
    showing_samples <- reactive({
      req(state$working)
      all_s <- colnames(state$working)
      by <- show_by_rv()
      if (identical(by, "__all__")) return(all_s)
      vals <- show_values_rv()
      if (!length(vals)) return(all_s)
      if (identical(by, "__samples__")) return(intersect(all_s, vals))
      cd <- as.data.frame(SummarizedExperiment::colData(state$working))
      all_s[as.character(cd[[by]]) %in% vals]
    })

    # --- Sample QC: General QC + Matrix -------------------------------------
    qc_tbl <- reactive({
      req(state$working)
      state_derive(state, "qc_metrics", params = list(),
                   expr = function() qc_per_sample_metrics(state$working))
    })

    # General QC colour-by gains the removal-status / pool options.
    output$group_ui <- renderUI({
      req(state$working)
      cols <- colnames(SummarizedExperiment::colData(state$working))
      selectInput(ns("group"), "Group / colour by",
                  choices = c(stats::setNames(cols, cols),
                              "Removal status (this metric)" = "__removal__",
                              "In removal pool" = "__pool__"),
                  selected = default_group_col())
    })
    output$auto_ui  <- auto_box("auto")

    # Colour aesthetic for the General QC plot: a colData column, the reason-aware
    # removal status (green/yellow/red), or removal-pool membership.
    sample_aes <- function(col, metric, samples) {
      if (identical(col, "__removal__")) {
        fl <- samp_flags(); i <- match(samples, fl$sample)
        rcol <- .metric_reason[[metric]]
        this <- if (!is.null(rcol)) fl[[rcol]][i] else NULL
        list(values = removal_status(fl$flagged[i], this), lab = "Removal status",
             palette = .removal_palette, labels = .removal_labels)
      } else if (identical(col, "__pool__")) {
        inp <- samples %in% samp_pool()
        list(values = factor(ifelse(inp, "in removal pool", "kept"),
                             levels = c("kept", "in removal pool")),
             lab = "Removal pool",
             palette = c(kept = "#9aa0a6", "in removal pool" = "#D62728"), labels = NULL)
      } else {
        cd <- as.data.frame(SummarizedExperiment::colData(state$working))
        v <- if (!is.null(col) && col %in% colnames(cd)) as.factor(cd[samples, col])
             else factor(rep("all", length(samples)))
        list(values = v, lab = col, palette = NULL, labels = NULL)
      }
    }

    current_spec <- reactive({
      req(state$working, input$x_axis, input$metric)
      list(tbl = qc_tbl(), x_axis = input$x_axis, metric = input$metric,
           sort = input$sort %||% "none", show = showing_samples())
    })
    gen_shown <- deferred("auto", "render", current_spec,
      sig = reactive(list(input$x_axis, input$metric, input$sort,
                          show_by_rv(), show_values_rv(), state$data_version)))
    output$gen_stale <- stale_note(gen_shown)

    output$plot <- renderPlot({
      validate(need(!is.null(gen_shown$value()),
                    "Click Render (or enable auto-render) to draw the plot."))
      s <- gen_shown$value()
      tbl <- s$tbl[s$tbl$sample %in% s$show, , drop = FALSE]
      validate(need(nrow(tbl) > 0, "No samples in the current 'Showing' selection."))
      ae <- sample_aes(input$group %||% default_group_col(), s$metric, tbl$sample)
      tbl$group <- ae$values
      .qc_metric_plot(tbl, s$x_axis, s$metric, ae$lab, sort = s$sort,
                      dark_theme = dark(), palette = ae$palette, palette_labels = ae$labels)
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
      list(df = df, n = ncol(dds), show = showing_samples())
    })
    rle_shown <- deferred("rle_auto", "rle_render", rle_spec,
      sig = reactive(list(input$rle_group, show_by_rv(), show_values_rv(), state$data_version)))
    output$rle_stale <- stale_note(rle_shown)
    output$diag_rle <- renderPlot({
      validate(need(!is.null(rle_shown$value()), "Click Render (or enable auto-render)."))
      s <- rle_shown$value()
      d <- s$df[s$df$sample %in% s$show, , drop = FALSE]
      d$sample <- droplevels(d$sample)
      validate(need(nrow(d) > 0, "No samples in the current 'Showing' selection."))
      .qc_rle_plot(d, dark(), length(levels(d$sample)))
    })

    dens_spec <- reactive({
      req(state$working, input$dens_group)
      dds <- state$working
      gmap <- group_lookup(input$dens_group)
      df <- qc_expression_long(dds)
      df$group <- factor(gmap[as.character(df$sample)])
      list(df = df, n = ncol(dds), show = showing_samples())
    })
    dens_shown <- deferred("dens_auto", "dens_render", dens_spec,
      sig = reactive(list(input$dens_group, show_by_rv(), show_values_rv(), state$data_version)))
    output$dens_stale <- stale_note(dens_shown)
    output$diag_density <- renderPlot({
      validate(need(!is.null(dens_shown$value()), "Click Render (or enable auto-render)."))
      s <- dens_shown$value()
      d <- s$df[s$df$sample %in% s$show, , drop = FALSE]
      d$sample <- droplevels(d$sample)
      validate(need(nrow(d) > 0, "No samples in the current 'Showing' selection."))
      .qc_density_plot(d, dark(), length(levels(d$sample)))
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
    meansd_shown <- deferred("diag_auto", "diag_render", meansd_spec,
      sig = reactive(list(state$data_version)))
    output$diag_stale <- stale_note(meansd_shown)
    output$diag_meansd <- renderPlot({
      validate(need(!is.null(meansd_shown$value()), "Click Render (or enable auto-render)."))
      .qc_meansd_plot(meansd_shown$value()$vst, dark_theme = dark())
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
      list(cm = cm, anno_df = anno_df, n = ncol(state$working), method = method,
           show = showing_samples())
    })
    cor_shown <- deferred("cor_auto", "cor_render", cor_spec,
      sig = reactive(list(input$cor_method, input$cor_anno, show_by_rv(),
                          show_values_rv(), state$data_version)))
    output$cor_stale <- stale_note(cor_shown)
    output$cor_plot <- renderPlot({
      validate(need(!is.null(cor_shown$value()), "Click Render to draw the correlation heatmap."))
      validate(need(requireNamespace("ComplexHeatmap", quietly = TRUE),
                    "Install 'ComplexHeatmap' to show the correlation heatmap."))
      s <- cor_shown$value()
      show <- intersect(colnames(s$cm), s$show)
      validate(need(length(show) > 1,
                    "Need at least two samples in the current 'Showing' selection."))
      cm <- s$cm[show, show, drop = FALSE]
      anno <- if (is.null(s$anno_df)) NULL else s$anno_df[show, , drop = FALSE]
      ComplexHeatmap::draw(.qc_correlation_heatmap(
        cm, anno, dark_theme = dark(), n_samples = length(show), method = s$method))
    })

    # --- Sample Correlation: Within-group -----------------------------------
    output$wg_group_ui <- group_box("wg_group", label = "Group by")
    output$wg_auto_ui  <- auto_box("wg_auto")
    wg_spec <- reactive({
      req(state$working, input$wg_group)
      list(df = qc_within_group_correlation(state$working, group = input$wg_group),
           show = showing_samples())
    })
    wg_shown <- deferred("wg_auto", "wg_render", wg_spec,
      sig = reactive(list(input$wg_group, show_by_rv(), show_values_rv(), state$data_version)))
    output$wg_stale <- stale_note(wg_shown)
    output$wg_plot <- renderPlot({
      validate(need(!is.null(wg_shown$value()), "Click Render (or enable auto-render)."))
      s <- wg_shown$value()
      d <- s$df[s$df$sample %in% s$show, , drop = FALSE]
      .qc_within_group_plot(d, dark_theme = dark())
    })

    # --- Spike-in (ERCC) dose-response --------------------------------------
    spike_ids <- reactive(rownames(state$working)[.detect_spike_features(state$working)])
    # A usable spike_concentration column = present + positive for some spike rows.
    conc_col_ok <- reactive({
      req(state$working); ids <- spike_ids()
      if (!length(ids)) return(FALSE)
      length(spike_features_missing_conc(state$working)) < length(ids)
    })
    output$spike_source_ui <- renderUI({
      req(state$working)
      choices <- c("ERCC Mix 1 (bundled)" = "mix1", "ERCC Mix 2 (bundled)" = "mix2")
      if (conc_col_ok()) choices <- c("Feature metadata (spike_concentration)" = "column", choices)
      selectInput(ns("spike_source"), "Concentration source", choices = choices,
                  selected = if (conc_col_ok()) "column" else "mix1")
    })
    output$spike_assay_ui <- renderUI({
      req(state$working)
      present <- intersect(c("CPM", "TPM", "FPKM"), SummarizedExperiment::assayNames(state$working))
      selectInput(ns("spike_assay"), "Observed assay", choices = union("CPM", present), selected = "CPM")
    })
    output$spike_group_ui <- group_box("spike_group")
    output$spike_auto_ui  <- auto_box("spike_auto")

    spike_spec <- reactive({
      req(state$working, input$spike_source, input$spike_assay)
      source <- input$spike_source; assay <- input$spike_assay
      state_derive(state, "spike_dr", params = list(source = source, assay = assay),
                   expr = function() spike_dose_response(state$working, assay = assay, source = source))
    })
    spike_shown <- deferred("spike_auto", "spike_render", spike_spec,
      sig = reactive(list(input$spike_source, input$spike_assay, state$data_version)))
    output$spike_stale <- stale_note(spike_shown)

    # group column for colouring (applied at render, not part of the cached compute)
    spike_group_col <- function(samples) {
      gmap <- group_lookup(input$spike_group %||% default_group_col())
      factor(gmap[as.character(samples)])
    }

    output$spike_msg <- renderUI({
      req(state$working)
      if (!length(spike_ids())) {
        return(tags$div(class = "alert alert-secondary py-2 small mb-2",
          "No spike-in features detected. Tag ERCC controls on the Feature tab to enable this view."))
      }
      v <- spike_shown$value(); if (is.null(v)) return(NULL)
      notes <- character(0)
      # Match rate when joining the bundled ERCC reference (non-ERCC / custom
      # spike ids resolve to NA and would yield a meaningless fit silently).
      if (input$spike_source %in% c("mix1", "mix2")) {
        cby <- v$long$concentration[!duplicated(v$long$feature)]
        matched <- sum(is.finite(cby))
        if (matched < length(cby)) notes <- c(notes, sprintf(
          "Only %d of %d spike-in id(s) matched the ERCC reference - check the Mix, or designate a concentration column on the Feature tab.",
          matched, length(cby)))
      }
      n_zero <- sum(v$per_sample$n_spike_detected == 0)
      if (n_zero > 0L) notes <- c(notes, sprintf(
        "%d of %d sample(s) have no detected spike-ins - often this just means they were not spiked (mixed designs). Consider removing spike-in features (Filtering tab) if your design has none.",
        n_zero, ncol(state$working)))
      if (!length(notes)) return(NULL)
      tags$div(class = "alert alert-warning py-2 small mb-2", lapply(notes, tags$div))
    })

    output$spike_dr_plot <- renderPlot({
      validate(need(length(spike_ids()) > 0, "No spike-in features in this dataset."))
      validate(need(!is.null(spike_shown$value()), "Click Render (or enable auto-render)."))
      dr <- spike_shown$value()
      validate(need(any(is.finite(dr$long$concentration)),
                    "No known concentrations resolved for these spike-ins - try another source."))
      long <- dr$long; long$group <- spike_group_col(long$sample)
      .qc_spike_dr_plot(long, dark_theme = dark())
    })

    output$spike_summary_plot <- renderPlot({
      validate(need(!is.null(spike_shown$value()), "Click Render (or enable auto-render)."))
      ps <- spike_shown$value()$per_sample
      ps$group <- spike_group_col(ps$sample)
      .qc_spike_summary_plot(ps, input$spike_metric %||% "r_squared", dark_theme = dark())
    })

    output$spike_cv <- renderUI({
      v <- spike_shown$value(); req(v)
      pct <- v$per_sample$pct_spike[is.finite(v$per_sample$pct_spike)]
      if (length(pct) < 2L || mean(pct) == 0) return(NULL)
      tags$div(class = "small text-muted mb-2",
        sprintf("Spike-in fraction across samples: mean %.2f%%, CV %.0f%% (high CV = uneven spike input).",
                mean(pct), 100 * stats::sd(pct) / mean(pct)))
    })

    output$spike_table <- DT::renderDT({
      validate(need(!is.null(spike_shown$value()), "Click Render (or enable auto-render)."))
      ps <- spike_shown$value()$per_sample
      dt_table(data.frame(
        Sample = ps$sample, `% spike` = round(ps$pct_spike, 2),
        Detected = ps$n_spike_detected, `Fit points` = ps$n_points,
        Slope = round(ps$slope, 3), `R-squared` = round(ps$r_squared, 3),
        `Lowest detected conc` = signif(ps$lod, 3), check.names = FALSE))
    })

    # --- Filtering: shared flags + removal pools ----------------------------
    samp_pool <- reactiveVal(character(0))
    feat_pool <- reactiveVal(character(0))

    # Advisory flags, cached on data_version + the rule inputs. Shared by the
    # Filtering tables and the General QC "Removal status" colour-by.
    samp_flags <- reactive({
      req(state$working)
      params <- list(group = input$samp_group, lib = .blank_na(input$samp_lib_min),
                     det = .blank_na(input$samp_detected_min),
                     mito = .blank_na(input$samp_mito_max), z = .blank_na(input$samp_wg_z))
      state_derive(state, "samp_flags", params = params, expr = function() {
        flag_samples(state$working, group = input$samp_group,
                     lib_size_min = .blank_na(input$samp_lib_min),
                     detected_min = .blank_na(input$samp_detected_min),
                     pct_mito_max = .blank_na(input$samp_mito_max),
                     within_group_z = .blank_na(input$samp_wg_z))
      })
    })
    feat_flags <- reactive({
      req(state$working)
      ms <- if (isTRUE(input$feat_use_min_samples)) input$feat_min_samples else NULL
      params <- list(fbe = isTRUE(input$feat_use_fbe), group = input$feat_group,
                     mc = input$feat_min_count %||% 10, ms = ms)
      state_derive(state, "feat_flags", params = params, expr = function() {
        flag_features(state$working, use_filter_by_expr = isTRUE(input$feat_use_fbe),
                      group = input$feat_group, min_count = input$feat_min_count %||% 10,
                      min_samples = ms)
      })
    })

    output$samp_group_ui <- group_box("samp_group", label = "Within-group grouping")
    output$feat_group_ui <- renderUI({
      req(state$working)
      cols <- colnames(SummarizedExperiment::colData(state$working))
      selectInput(ns("feat_group"), "filterByExpr grouping", choices = cols,
                  selected = default_group_col())
    })

    # Sample thresholds: fill data-driven defaults on load and on the Auto button
    # (a degenerate fence -> NA -> blank input -> that rule stays disabled).
    fill_samp_thresholds <- function() {
      req(state$working)
      th <- suggest_sample_thresholds(state$working)
      updateNumericInput(session, "samp_lib_min", value = th$lib_size_min)
      updateNumericInput(session, "samp_detected_min", value = th$detected_min)
      updateNumericInput(session, "samp_mito_max", value = th$pct_mito_max)
    }
    observeEvent(state$working, fill_samp_thresholds())
    observeEvent(input$samp_auto, fill_samp_thresholds())
    # Per-field Auto buttons fill just their own threshold.
    observeEvent(input$samp_lib_auto, updateNumericInput(
      session, "samp_lib_min", value = suggest_sample_thresholds(state$working)$lib_size_min))
    observeEvent(input$samp_detected_auto, updateNumericInput(
      session, "samp_detected_min", value = suggest_sample_thresholds(state$working)$detected_min))
    observeEvent(input$samp_mito_auto, updateNumericInput(
      session, "samp_mito_max", value = suggest_sample_thresholds(state$working)$pct_mito_max))

    # Features default to adopting the suggestion (re-seeded when rules change);
    # samples stay opt-in (highlight only). Both pools reset on any data change.
    observeEvent(feat_flags(), {
      feat_pool(feat_flags()$feature_id[feat_flags()$suggested_drop])
    })
    observeEvent(state$data_version, samp_pool(character(0)))

    # Tables: a read-only boolean "Suggested Removal" column + a separate boolean
    # "In Removal Pool" column; native row-selection is transient staging the
    # pool buttons act on.
    samp_display <- function(fl, pool) data.frame(
      Sample = fl$sample, `Library size (M)` = round(fl$library_size / 1e6, 3),
      Detected = fl$detected, `% mito` = round(fl$pct_mito, 2),
      `Within-grp corr` = round(fl$within_group_corr, 3),
      `Suggested Removal` = fl$flagged, Reason = fl$reason,
      `In Removal Pool` = fl$sample %in% pool,
      check.names = FALSE, stringsAsFactors = FALSE)
    feat_display <- function(fl, pool) data.frame(
      Feature = fl$feature_id, Class = fl$feature_class,
      `Total count` = round(fl$total_count), Detected = fl$n_detected,
      `Mean log` = round(fl$mean_logcounts, 2),
      `Suggested Removal` = fl$suggested_drop, Reason = fl$reason,
      `In Removal Pool` = fl$feature_id %in% pool,
      check.names = FALSE, stringsAsFactors = FALSE)

    output$samp_tbl <- DT::renderDT({
      validate(need(!is.null(state$working), "No dataset loaded."))
      dt_table(samp_display(samp_flags(), isolate(samp_pool())),
               selection = list(mode = "multiple"))
    }, server = TRUE)
    output$feat_tbl <- DT::renderDT({
      validate(need(!is.null(state$working), "No dataset loaded."))
      dt_table(feat_display(feat_flags(), isolate(feat_pool())),
               selection = list(mode = "multiple"))
    }, server = TRUE)
    samp_proxy <- DT::dataTableProxy("samp_tbl")
    feat_proxy <- DT::dataTableProxy("feat_tbl")
    # Refresh only the "In pool" column on pool change (keeps paging + filter).
    observeEvent(samp_pool(), DT::replaceData(
      samp_proxy, samp_display(samp_flags(), samp_pool()),
      resetPaging = FALSE, clearSelection = "none", rownames = FALSE), ignoreInit = TRUE)
    observeEvent(feat_pool(), DT::replaceData(
      feat_proxy, feat_display(feat_flags(), feat_pool()),
      resetPaging = FALSE, clearSelection = "none", rownames = FALSE), ignoreInit = TRUE)

    # Map the (transient) DT row selection back to ids.
    samp_sel  <- reactive(samp_flags()$sample[input$samp_tbl_rows_selected])
    feat_sel  <- reactive(feat_flags()$feature_id[input$feat_tbl_rows_selected])

    # Select all / Deselect all stage the rows matching the current table search.
    observeEvent(input$samp_select_all,   DT::selectRows(samp_proxy, input$samp_tbl_rows_all))
    observeEvent(input$samp_deselect_all, DT::selectRows(samp_proxy, NULL))
    observeEvent(input$feat_select_all,   DT::selectRows(feat_proxy, input$feat_tbl_rows_all))
    observeEvent(input$feat_deselect_all, DT::selectRows(feat_proxy, NULL))

    # Staging selection is consumed once moved to/from the pool -> clear it.
    observeEvent(input$samp_add,    { samp_pool(union(samp_pool(), samp_sel()));   DT::selectRows(samp_proxy, NULL) })
    observeEvent(input$samp_remove, { samp_pool(setdiff(samp_pool(), samp_sel())); DT::selectRows(samp_proxy, NULL) })
    observeEvent(input$samp_adopt,    samp_pool(union(samp_pool(),
                   samp_flags()$sample[samp_flags()$flagged])))
    observeEvent(input$samp_clear,    samp_pool(character(0)))
    observeEvent(input$feat_add,    { feat_pool(union(feat_pool(), feat_sel()));   DT::selectRows(feat_proxy, NULL) })
    observeEvent(input$feat_remove, { feat_pool(setdiff(feat_pool(), feat_sel())); DT::selectRows(feat_proxy, NULL) })
    observeEvent(input$feat_adopt,    feat_pool(union(feat_pool(),
                   feat_flags()$feature_id[feat_flags()$suggested_drop])))
    observeEvent(input$feat_clear,    feat_pool(character(0)))

    output$samp_counts <- renderUI({
      req(state$working)
      .pool_counts(sum(samp_flags()$flagged), length(input$samp_tbl_rows_selected),
                   length(samp_pool()), ncol(state$working))
    })
    output$feat_counts <- renderUI({
      req(state$working)
      .pool_counts(sum(feat_flags()$suggested_drop), length(input$feat_tbl_rows_selected),
                   length(feat_pool()), nrow(state$working))
    })

    # Before/after filtering density: keep = everything not in the removal pool.
    output$feat_density <- renderPlot({
      req(state$working)
      keep <- setdiff(rownames(state$working), feat_pool())
      .qc_filter_density_plot(qc_filter_density(state$working, keep), dark_theme = dark())
    })

    # Apply removal: confirm, then a single state_mutate (undoable, logged).
    confirm_modal <- function(ok_id, msg) showModal(modalDialog(
      title = "Confirm removal", msg, easyClose = TRUE,
      footer = tagList(modalButton("Cancel"),
                       actionButton(ns(ok_id), "Remove", class = "btn-danger"))))
    observeEvent(input$samp_apply, {
      req(length(samp_pool()) > 0)
      confirm_modal("samp_apply_ok",
        sprintf("Remove %d sample(s) from the working dataset? Your original import is kept and can be restored by reloading it.",
                length(samp_pool())))
    })
    observeEvent(input$feat_apply, {
      req(length(feat_pool()) > 0)
      confirm_modal("feat_apply_ok",
        sprintf("Remove %d feature(s) from the working dataset? Your original import is kept and can be restored by reloading it.",
                length(feat_pool())))
    })
    observeEvent(input$samp_apply_ok, {
      ids <- samp_pool(); removeModal(); req(length(ids) > 0)
      state_mutate(state, function(d) drop_samples(d, ids),
                   action = list(action = "filter_samples", n_dropped = length(ids),
                                 dropped = ids))
    })
    # Reset removal: re-add everything removed along that dimension (keeps other
    # edits), itself an undoable state_mutate. No-op notice when nothing removed.
    observeEvent(input$samp_reset, {
      req(state$working, state$original)
      removed <- setdiff(colnames(state$original), colnames(state$working))
      if (!length(removed)) {
        showNotification("No samples have been removed.", type = "message"); return()
      }
      state_mutate(state, function(d) restore_samples(d, state$original),
                   action = list(action = "restore_samples", n_restored = length(removed)))
      showNotification(sprintf("Restored %d removed sample(s).", length(removed)), type = "message")
    })
    observeEvent(input$feat_reset, {
      req(state$working, state$original)
      removed <- setdiff(rownames(state$original), rownames(state$working))
      if (!length(removed)) {
        showNotification("No features have been removed.", type = "message"); return()
      }
      state_mutate(state, function(d) restore_features(d, state$original),
                   action = list(action = "restore_features", n_restored = length(removed)))
      showNotification(sprintf("Restored %d removed feature(s).", length(removed)), type = "message")
    })

    # Remove all spike-in features (sound: size factors use endogenous controlGenes,
    # so this does not affect normalization). Reversible via "Reset Feature Removal".
    observeEvent(input$feat_drop_spike, {
      req(state$working)
      ids <- rownames(state$working)[.detect_spike_features(state$working)]
      if (!length(ids)) {
        showNotification("No spike-in features to remove.", type = "message"); return()
      }
      showModal(modalDialog(
        title = "Remove all spike-in features?",
        sprintf("Drop %d spike-in (ERCC) feature(s) from the working dataset? This removes the spike-in QC; restore via 'Reset Feature Removal'.", length(ids)),
        easyClose = TRUE,
        footer = tagList(modalButton("Cancel"),
                         actionButton(ns("feat_drop_spike_ok"), "Remove", class = "btn-danger"))))
    })
    observeEvent(input$feat_drop_spike_ok, {
      removeModal()
      ids <- rownames(state$working)[.detect_spike_features(state$working)]
      req(length(ids) > 0)
      state_mutate(state, function(d) drop_features(d, ids),
                   action = list(action = "drop_spike_in", n_dropped = length(ids)))
      showNotification(sprintf("Removed %d spike-in feature(s).", length(ids)), type = "message")
    })

    observeEvent(input$feat_apply_ok, {
      ids <- feat_pool(); removeModal(); req(length(ids) > 0)
      # Guard: a removal that zeroes a sample's library would make CPM/logcounts
      # undefined (cpm() errors on a zero-total sample). Block with a clear note.
      keep <- setdiff(rownames(state$working), ids)
      lib <- colSums(as.matrix(SummarizedExperiment::assay(state$working, "counts"))[keep, , drop = FALSE])
      if (any(lib == 0)) {
        showNotification(
          sprintf("Removal would empty %d sample(s) to zero counts (%s). Trim the pool or drop those samples first.",
                  sum(lib == 0), paste(colnames(state$working)[lib == 0], collapse = ", ")),
          type = "error", duration = NULL)
        return(invisible())
      }
      state_mutate(state, function(d) drop_features(d, ids),
                   action = list(action = "filter_features",
                                 method = if (isTRUE(input$feat_use_fbe)) "filterByExpr" else "manual",
                                 min_count = input$feat_min_count %||% 10,
                                 n_dropped = length(ids), dropped = ids))
    })

    invisible(NULL)
  })
}
