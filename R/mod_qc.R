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

# The reason-aware removal-status colour scheme (.removal_palette / .removal_labels)
# and its resolver removal_status_colors() now live in filter_helpers.R, shared
# with the PCA "Suggested removal" aesthetic.
# Sample QC metric -> the flag_samples() reason column it corresponds to (used to
# pick the "this reason" highlight). One column per metric, so % spike-in maps to
# the over-spiked side only; under-spiked samples (low_spike) still show as
# "suggested (other)" on that plot.
.metric_reason <- c(library_size = "low_lib_size", detected = "low_detected",
                    pct_mito = "high_mito", pct_spike = "high_spike")

# Advisory notes for a dose-response result, shared by the Spike-in QC tab
# message and the Filtering sidebar: the ERCC reference match rate (when a
# bundled Mix is the source) and the count of samples with no detected spikes.
.spike_match_notes <- function(dr, source, n_samples) {
  notes <- character(0)
  if (source %in% c("mix1", "mix2")) {
    cby <- dr$long$concentration[!duplicated(dr$long$feature)]
    matched <- sum(is.finite(cby))
    if (matched < length(cby)) notes <- c(notes, sprintf(
      "Only %d of %d spike-in id(s) matched the ERCC reference - check the Mix, or designate a concentration column on the Feature tab.",
      matched, length(cby)))
  }
  n_zero <- sum(dr$per_sample$n_spike_detected == 0)
  if (n_zero > 0L) notes <- c(notes, sprintf(
    "%d of %d sample(s) have no detected spike-ins - often this just means they were not spiked (mixed designs). Consider removing spike-in features (Filtering tab) if your design has none.",
    n_zero, n_samples))
  notes
}

# A threshold numericInput with its own per-field "Auto" wand button alongside.
# The label sits on its own line; the input (wide) + the smaller wand button sit
# on the next line, vertically centred. `tip` is the button's hover text. Used in
# the QC sample-filter accordion (built server-side via renderUI so the general
# and spike-in panels are true siblings in one accordion card).
.qc_thr_input <- function(ns, input_id, label, auto_id, tip, ...) {
  tags$div(class = "mb-2",
    tags$label(label, class = "form-label mb-1", `for` = ns(input_id)),
    tags$div(class = "d-flex align-items-center gap-1",
      tags$div(class = "flex-grow-1", numericInput(ns(input_id), label = NULL, ...)),
      bslib::tooltip(
        actionButton(ns(auto_id), NULL, icon = icon("wand-magic-sparkles"),
                     class = "btn-sm btn-outline-primary"),
        tip)))
}

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

# Plot theme moved to R/mod_plot_engine.R as the shared `.plot_theme()` (reused
# by the PCA page); the QC builders call it below.

# Build the per-sample QC plot. x_var = "sample" (discrete bar of the metric per
# sample, optionally sorted by value) or another metric name (numeric scatter,
# metric vs metric). Colour/fill by `group`. When `palette` (a named colour
# vector over the group levels) is supplied, a fixed manual scale is used instead
# of the qualitative palette (for the semantic removal-status colouring).
# A hover `text` aesthetic, added to a geom only when `interactive` (so static
# ggplot never warns about an unknown aesthetic). Pair with ggplotly(tooltip =
# "text"). The data passed to the geom must carry a `text` column.
.hover_aes <- function(interactive) {
  if (isTRUE(interactive)) ggplot2::aes(text = .data$text) else NULL
}

# Fixed discrete group colour scale from the project palette (a named
# level -> colour vector, via palette_discrete()); NULL keeps thematic's default
# scale. `aes` selects which scales to set: RLE/summary use "fill", density/dose
# use "colour", within-group uses both. Returns a list of scales (or NULL) to add
# to a ggplot with `+`.
.qc_group_scale <- function(palette, aes = "colour") {
  if (is.null(palette)) return(NULL)
  out <- list()
  if ("fill" %in% aes)
    out <- c(out, list(ggplot2::scale_fill_manual(values = palette, drop = FALSE, name = "group")))
  if ("colour" %in% aes)
    out <- c(out, list(ggplot2::scale_colour_manual(values = palette, drop = FALSE, name = "group")))
  out
}

.qc_metric_plot <- function(tbl, x_var, metric, group_lab = NULL,
                            sort = "none", dark_theme = FALSE,
                            palette = NULL, palette_labels = NULL,
                            interactive = FALSE) {
  yy <- .qc_axis(tbl, metric)
  bar <- identical(x_var, "sample")
  if (bar) {
    lvls <- if (identical(sort, "none")) {
      tbl$sample
    } else {
      tbl$sample[order(yy$v, decreasing = identical(sort, "decreasing"))]
    }
    df <- data.frame(x = factor(tbl$sample, levels = lvls), y = yy$v, group = tbl$group,
                     text = sprintf("Sample: %s<br>%s: %s", tbl$sample, yy$lab, signif(yy$v, 4)))
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = .data$y, fill = .data$group)) +
      ggplot2::geom_col(mapping = .hover_aes(interactive)) +
      ggplot2::labs(x = "sample", y = yy$lab, fill = group_lab %||% "group") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  } else {
    xx <- .qc_axis(tbl, x_var)
    df <- data.frame(x = xx$v, y = yy$v, group = tbl$group,
                     text = sprintf("Sample: %s<br>%s: %s<br>%s: %s", tbl$sample,
                                    xx$lab, signif(xx$v, 4), yy$lab, signif(yy$v, 4)))
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = .data$y, colour = .data$group)) +
      ggplot2::geom_point(size = 2, mapping = .hover_aes(interactive)) +
      ggplot2::labs(x = xx$lab, y = yy$lab, colour = group_lab %||% "group")
  }
  if (!is.null(palette)) {
    labs_arg <- if (is.null(palette_labels)) ggplot2::waiver() else palette_labels
    scale_fn <- if (bar) ggplot2::scale_fill_manual else ggplot2::scale_colour_manual
    p <- p + scale_fn(values = palette, labels = labs_arg, drop = FALSE,
                      name = group_lab %||% "group")
  }
  p + .plot_theme(dark_theme)
}

# Before/after feature-filter density (limma/edgeR RNAseq123 style): the
# distribution of log-expression over endogenous features, before vs after the
# proposed filter. A cleaner after-curve (less low-expression mass) means the
# filter is doing its job. df from qc_filter_density() (sample, value, status).
.qc_filter_density_plot <- function(df, dark_theme = FALSE, interactive = FALSE) {
  df$text <- sprintf("Sample: %s (%s)", df$sample, df$status)
  ggplot2::ggplot(df, ggplot2::aes(x = .data$value,
                                   group = interaction(.data$sample, .data$status),
                                   colour = .data$status)) +
    ggplot2::geom_density(linewidth = 0.4, alpha = 0.7, mapping = .hover_aes(interactive)) +
    ggplot2::scale_colour_manual(values = c(before = "#9aa0a6", after = "#1f77b4"),
                                 name = NULL) +
    ggplot2::labs(x = "log2 expression", y = "density") +
    .plot_theme(dark_theme)
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
.qc_spike_dr_plot <- function(long, dark_theme = FALSE, interactive = FALSE,
                              palette = NULL) {
  d <- long[is.finite(long$concentration) & is.finite(long$expression) &
            long$concentration > 0 & long$expression > 0, , drop = FALSE]
  if (!nrow(d)) return(.plot_msg("No detected spike-ins with known concentration."))
  d$text <- sprintf("Feature: %s<br>Sample: %s<br>conc: %s<br>expr: %s",
                    d$feature, d$sample, signif(d$concentration, 4), signif(d$expression, 4))
  ggplot2::ggplot(d, ggplot2::aes(x = .data$concentration, y = .data$expression,
                                  colour = .data$group)) +
    ggplot2::geom_point(size = 1.4, alpha = 0.7, mapping = .hover_aes(interactive)) +
    ggplot2::geom_smooth(ggplot2::aes(group = .data$sample), method = "lm",
                         formula = y ~ x, se = FALSE, linewidth = 0.4) +
    ggplot2::scale_x_log10() + ggplot2::scale_y_log10() +
    ggplot2::labs(x = "known concentration (attomoles/uL)", y = "observed expression",
                  colour = "group") +
    .qc_group_scale(palette, "colour") +
    .plot_theme(dark_theme)
}

# Per-sample spike summary: a chosen metric across samples (bar, sorted), filled
# by group. `df` is spike_dose_response()$per_sample with a `group` column.
.qc_spike_summary_plot <- function(df, metric, dark_theme = FALSE, interactive = FALSE,
                                   palette = NULL, sort = "none") {
  lab <- .spike_metric_labels[[metric]] %||% metric
  d <- df[is.finite(df[[metric]]), , drop = FALSE]
  if (!nrow(d)) return(.plot_msg(paste("No", lab, "to show.")))
  d$text <- sprintf("Sample: %s<br>%s: %s", d$sample, lab, signif(d[[metric]], 4))
  lvls <- if (identical(sort, "none")) d$sample
          else d$sample[order(d[[metric]], decreasing = identical(sort, "decreasing"))]
  d$sample <- factor(d$sample, levels = lvls)
  ggplot2::ggplot(d, ggplot2::aes(x = .data$sample, y = .data[[metric]], fill = .data$group)) +
    ggplot2::geom_col(mapping = .hover_aes(interactive)) +
    ggplot2::labs(x = "sample", y = lab, fill = "group") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    .qc_group_scale(palette, "fill") +
    .plot_theme(dark_theme)
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
    .plot_theme(dark_theme)
}

# Relative log expression boxplots. df: sample (factor), group, value.
.qc_rle_plot <- function(df, dark_theme = FALSE, n_samples = NULL, interactive = FALSE,
                         palette = NULL) {
  n <- n_samples %||% length(unique(df$sample))
  df$text <- sprintf("Sample: %s", df$sample)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$sample, y = .data$value, fill = .data$group)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    ggplot2::geom_boxplot(outlier.size = 0.4, mapping = .hover_aes(interactive)) +
    ggplot2::labs(x = "sample", y = "relative log expression", fill = "group") +
    .qc_group_scale(palette, "fill") +
    .plot_theme(dark_theme)
  if (n > 30) {
    p + ggplot2::theme(axis.text.x = ggplot2::element_blank())
  } else {
    p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  }
}

# Per-sample log-expression density. df: sample (factor), group, value.
.qc_density_plot <- function(df, dark_theme = FALSE, n_samples = NULL, interactive = FALSE,
                             palette = NULL) {
  n <- n_samples %||% length(unique(df$sample))
  df$text <- sprintf("Sample: %s", df$sample)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$value, group = .data$sample,
                                        colour = .data$group)) +
    ggplot2::geom_density(mapping = .hover_aes(interactive)) +
    ggplot2::labs(x = "log2 expression", y = "density", colour = "group") +
    .qc_group_scale(palette, "colour") +
    .plot_theme(dark_theme)
  if (n > 30) p + ggplot2::theme(legend.position = "none") else p
}

# Within-group correlation: boxplot per group + sample points. df from
# qc_within_group_correlation() (sample, group, mean_corr). Singleton groups
# (NA mean_corr) are dropped.
.qc_within_group_plot <- function(df, dark_theme = FALSE, interactive = FALSE,
                                  palette = NULL) {
  d <- df[!is.na(df$mean_corr), , drop = FALSE]
  if (!nrow(d)) return(.plot_msg("No multi-sample groups to summarize."))
  d$text <- sprintf("Sample: %s<br>mean corr: %s", d$sample, signif(d$mean_corr, 4))
  jit_aes <- if (isTRUE(interactive)) {
    ggplot2::aes(colour = .data$group, text = .data$text)
  } else {
    ggplot2::aes(colour = .data$group)
  }
  ggplot2::ggplot(d, ggplot2::aes(x = .data$group, y = .data$mean_corr)) +
    ggplot2::geom_boxplot(ggplot2::aes(fill = .data$group), alpha = 0.35,
                          outlier.shape = NA) +
    ggplot2::geom_jitter(mapping = jit_aes, width = 0.15, height = 0, size = 2) +
    ggplot2::labs(x = "group", y = "mean within-group correlation",
                  fill = "group", colour = "group") +
    .qc_group_scale(palette, c("fill", "colour")) +
    .plot_theme(dark_theme)
}

# Sample-sample correlation heatmap (ComplexHeatmap). Returns a Heatmap, or NULL
# when ComplexHeatmap is unavailable (the module shows a message instead).
# `anno_df` is a samples-by-columns data.frame of metadata to annotate (or NULL);
# colours come from qc_annotation_colors() so they stay stable across re-renders.
# `method` is shown in the legend title + plot title. Placeholder body palette -
# tune via the Themer "Heatmap" sub-tab.
.qc_correlation_heatmap <- function(cor_mat, anno_df = NULL, dark_theme = FALSE,
                                    n_samples = ncol(cor_mat), method = "spearman",
                                    palette_config = NULL, cor_config = NULL,
                                    anno_col = NULL) {
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) return(NULL)
  fg <- if (isTRUE(dark_theme)) "grey90" else "grey10"
  legend_lab <- if (identical(method, "pearson")) "Pearson r" else "Spearman rho"
  title_lab  <- if (identical(method, "pearson")) "Pearson" else "Spearman"
  rng <- range(cor_mat, na.rm = TRUE)
  col_fun <- if (!is.null(cor_config) && !is.null(cor_config$name)) {
    palette_colorramp2(cor_config$name, values = as.numeric(cor_mat),
                       min = cor_config$min, max = cor_config$max,
                       custom = cor_config$custom, reverse = isTRUE(cor_config$reverse))
  } else if (requireNamespace("circlize", quietly = TRUE)) {
    circlize::colorRamp2(c(rng[1], mean(rng), rng[2]),
                         c("#fff7fb", "#74a9cf", "#023858"))
  } else {
    NULL
  }
  top <- NULL
  if (!is.null(anno_df) && ncol(as.data.frame(anno_df)) > 0) {
    adf <- as.data.frame(anno_df)[colnames(cor_mat), , drop = FALSE]
    # Precomputed colours (from the shared resolver) when given; else fall back to
    # the colData-config path for any other caller.
    cols <- anno_col %||% qc_annotation_colors(adf, palette_config)
    top <- ComplexHeatmap::HeatmapAnnotation(df = adf, col = cols)
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

# The ggplot <-> plotly engine toggle (P3f), deferred render, and the plot
# helpers (.plot_msg / .plot_dual / .plotly_max_elements / .muffle_unknown_aes /
# .to_plotly / plot_engine_server) now live in R/mod_plot_engine.R so PCA/DE
# reuse them (P4-pre).

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
  # The view-only "Showing:" control (reusable; see R/mod_plot_subset.R) sits in
  # each plot sidebar via plot_subset_ui(ns, <suffix>), all synced by the server.
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
              plot_subset_ui(ns, "gen"),
              uiOutput(ns("auto_ui")),
              actionButton(ns("render"), "Render", class = "btn-primary")
            ),
            uiOutput(ns("gen_stale")), .plot_dual(ns("plot_container"))
          )
        ),
        bslib::nav_panel(
          "RLE",
          bslib::layout_sidebar(
            sidebar = bslib::sidebar(
              title = tags$h4("RLE", class = "fs-6 mb-0"), width = 280,
              uiOutput(ns("rle_group_ui")),
              plot_subset_ui(ns, "rle"),
              uiOutput(ns("rle_auto_ui")),
              actionButton(ns("rle_render"), "Render", class = "btn-primary")
            ),
            uiOutput(ns("rle_stale")), .plot_dual(ns("diag_rle_container")), .qc_help_note("rle")
          )
        ),
        bslib::nav_panel(
          "Expression density",
          bslib::layout_sidebar(
            sidebar = bslib::sidebar(
              title = tags$h4("Expression density", class = "fs-6 mb-0"), width = 280,
              uiOutput(ns("dens_group_ui")),
              plot_subset_ui(ns, "dens"),
              uiOutput(ns("dens_auto_ui")),
              actionButton(ns("dens_render"), "Render", class = "btn-primary")
            ),
            uiOutput(ns("dens_stale")), .plot_dual(ns("diag_density_container")), .qc_help_note("density")
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
              plot_subset_ui(ns, "cor"),
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
              plot_subset_ui(ns, "wg"),
              uiOutput(ns("wg_auto_ui")),
              actionButton(ns("wg_render"), "Render", class = "btn-primary")
            ),
            uiOutput(ns("wg_stale")), .plot_dual(ns("wg_plot_container")), .qc_help_note("within_group")
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
          plot_subset_ui(ns, "spike"),
          uiOutput(ns("spike_auto_ui")),
          actionButton(ns("spike_render"), "Render", class = "btn-primary")
        ),
        uiOutput(ns("spike_msg")),
        uiOutput(ns("spike_stale")),
        bslib::navset_pill(
          bslib::nav_panel("Dose-response",
            .plot_dual(ns("spike_dr_plot_container")), .qc_help_note("spike_dr")),
          bslib::nav_panel("Per-sample summary",
            selectInput(ns("spike_metric"), "Metric",
                        choices = c("% spike-in" = "pct_spike",
                                    "Detected spike-in features" = "n_spike_detected",
                                    "Lowest detected conc" = "lod",
                                    "Dose-response slope" = "slope",
                                    "Dose-response R-squared" = "r_squared"),
                        selected = "r_squared"),
            selectInput(ns("spike_sort"), "Sort by",
                        choices = c("None" = "none", "Decreasing" = "decreasing",
                                    "Increasing" = "increasing"),
                        selected = "none"),
            uiOutput(ns("spike_cv")),
            .plot_dual(ns("spike_summary_plot_container"))),
          bslib::nav_panel("Spike-in QC Matrix",
            tags$small(class = "text-muted", "Per-sample spike-in QC metrics:"),
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
              # One accordion with the general sample-QC filters and (when the
              # dataset has spike-ins) the spike-in filters as sibling panels.
              # Built server-side (renderUI) so the conditional spike panel is a
              # true sibling in the same accordion card, not a separate card.
              uiOutput(ns("samp_filters_ui")),
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
            .plot_dual(ns("feat_density_container")), .qc_help_note("filter_density")
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

    # Shared plot engine (R/mod_plot_engine.R): the ggplot<->plotly toggle, the
    # deferred-render gate, and the stale banner -- reused by PCA/DE (P4-pre).
    .engine <- plot_engine_server(input, output, session, state)
    use_plotly_base <- .engine$use_plotly_base
    dual_plot       <- .engine$dual_plot
    deferred        <- .engine$deferred
    stale_note      <- .engine$stale_note

    # Default grouping column: the design variable / condition, else first col.
    default_group_col <- function() {
      cols <- colnames(SummarizedExperiment::colData(state$working))
      if ("condition" %in% cols) "condition" else cols[1]
    }
    # Per-sample grouping factor (named by sample) for the chosen field, handling
    # both colData columns and the "This session" removal items (promoted state).
    # Used by the per-sample colour plots (General/RLE/density/spike). The removal
    # status here is the 2-level pass/suggested view (the reason-aware 3-level
    # scheme needs a plot metric, so only General QC's sample_aes uses it).
    # Per-sample grouping factor (named by sample) for a colour plot, via the
    # shared attribute resolver (colData / removal / pool). Falls back to a single
    # "flags pending" level when removal flags are not computed yet.
    group_map <- function(col, samples = colnames(state$working)) {
      res <- aes_resolve(state, col, samples)
      vals <- if (is.null(res)) factor(rep("flags pending", length(samples)))
              else if (is.factor(res$values)) res$values else factor(res$values)
      stats::setNames(vals, samples)
    }
    # Colours for a colour plot's group-field levels (subset to those present), or
    # NULL when a colData column has no project palette config (-> thematic default).
    group_colours <- function(col, levels) {
      res <- aes_resolve(state, col, colnames(state$working))
      if (is.null(res) || is.null(res$colors)) return(NULL)
      res$colors[intersect(as.character(levels), names(res$colors))]
    }

    # --- DRY UI builders (data-dependent, so server-rendered) ---------------
    # `session = TRUE` adds the shared "This session" removal items (for the
    # per-sample colour plots); grouping-semantic selectors keep it FALSE.
    group_box <- function(input_id, label = "Group / colour by", session = FALSE) renderUI({
      req(state$working)
      cols  <- colnames(SummarizedExperiment::colData(state$working))
      items <- if (isTRUE(session)) .session_removal_items else NULL
      selectInput(ns(input_id), label,
                  choices = group_field_choices(cols, items, none = FALSE),
                  selected = default_group_col())
    })
    # A deferred-sig fragment so a plot grouped by a session item re-stales when
    # the promoted pool / flags change (no-op for colData grouping).
    .session_sig <- function(col)
      if (identical(col, "__removal__") || identical(col, "__pool__"))
        list(state$samp_pool, state$samp_flags) else NULL
    auto_box <- function(input_id) renderUI({
      req(state$working)
      checkboxInput(ns(input_id), "Auto-render", value = ncol(state$working) <= 150L)
    })

    output$dtype_badge <- renderUI(.dtype_badge_ui(state))
    output$diag_badge  <- renderUI(.dtype_badge_ui(state))

    # --- "Showing:" display subset (view-only; never mutates the dds) --------
    # Reusable control synced across the sample-plot sidebars (R/mod_plot_subset.R);
    # returns the samples currently shown. Plots filter their data by it; the
    # deferred `sig`s below depend on showing_samples() so a Showing change marks a
    # manual plot stale. Dataset-diagnostic / Spike-in tabs are feature-level (omitted).
    showing_samples <- plot_subset_server(input, output, session, state,
                                          suffixes = c("gen", "rle", "dens", "wg", "cor", "spike"))

    # --- Sample QC: General QC + Matrix -------------------------------------
    qc_tbl <- reactive({
      req(state$working)
      state_derive(state, "qc_metrics", params = list(),
                   expr = function() qc_per_sample_metrics(state$working))
    })

    # General QC colour-by gains the removal-status / pool options (grouped).
    output$group_ui <- renderUI({
      req(state$working)
      cols <- colnames(SummarizedExperiment::colData(state$working))
      selectInput(ns("group"), "Group / colour by",
                  choices = group_field_choices(cols, .session_removal_items, none = FALSE),
                  selected = default_group_col())
    })
    output$auto_ui  <- auto_box("auto")

    # Colour aesthetic for the General QC plot via the shared resolver. The plot
    # metric drives the reason-aware removal highlight (green/amber/red) through
    # `ctx$reason`; colData / pool resolve as elsewhere.
    sample_aes <- function(col, metric, samples) {
      res <- aes_resolve(state, col, samples, ctx = list(reason = metric))
      if (is.null(res))                          # removal before flags computed (rare)
        return(list(values = factor(rep("flags pending", length(samples))),
                    lab = "Suggested removal", palette = NULL, labels = NULL))
      # General QC groups discretely (a numeric colData column is binned to a
      # factor, as the pre-resolver code did); mirrors group_map() for RLE/density.
      vals <- if (is.factor(res$values)) res$values else factor(res$values)
      list(values = vals, lab = res$label, palette = res$colors, labels = res$labels)
    }

    current_spec <- reactive({
      req(state$working, input$x_axis, input$metric)
      list(tbl = qc_tbl(), x_axis = input$x_axis, metric = input$metric,
           sort = input$sort %||% "none", show = showing_samples())
    })
    gen_shown <- deferred("auto", "render", current_spec,
      sig = reactive(list(input$x_axis, input$metric, input$sort,
                          showing_samples(), state$data_version)))
    output$gen_stale <- stale_note(gen_shown)

    plot_gg <- function(interactive) {
      validate(need(!is.null(gen_shown$value()),
                    "Click Render (or enable auto-render) to draw the plot."))
      s <- gen_shown$value()
      tbl <- s$tbl[s$tbl$sample %in% s$show, , drop = FALSE]
      validate(need(nrow(tbl) > 0, "No samples in the current 'Showing' selection."))
      ae <- sample_aes(input$group %||% default_group_col(), s$metric, tbl$sample)
      tbl$group <- ae$values
      .qc_metric_plot(tbl, s$x_axis, s$metric, ae$lab, sort = s$sort,
                      dark_theme = dark(), palette = ae$palette, palette_labels = ae$labels,
                      interactive = interactive)
    }
    dual_plot("plot", plot_gg, n_elements = reactive({
      v <- gen_shown$value(); if (is.null(v)) 0L else sum(v$tbl$sample %in% v$show)
    }))

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
    output$rle_group_ui  <- group_box("rle_group", session = TRUE)
    output$rle_auto_ui   <- auto_box("rle_auto")
    output$dens_group_ui <- group_box("dens_group", session = TRUE)
    output$dens_auto_ui  <- auto_box("dens_auto")

    rle_spec <- reactive({
      req(state$working, input$rle_group)
      dds <- state$working
      rle <- qc_rle_matrix(dds)
      nf <- nrow(rle)
      gmap <- group_map(input$rle_group, colnames(rle))
      df <- data.frame(
        sample = factor(rep(colnames(rle), each = nf), levels = colnames(rle)),
        value  = as.numeric(rle))
      df$group <- gmap[as.character(df$sample)]
      list(df = df, n = ncol(dds), show = showing_samples())
    })
    rle_shown <- deferred("rle_auto", "rle_render", rle_spec,
      sig = reactive(list(input$rle_group, showing_samples(), state$data_version,
                          .session_sig(input$rle_group))))
    output$rle_stale <- stale_note(rle_shown)
    rle_gg <- function(interactive) {
      validate(need(!is.null(rle_shown$value()), "Click Render (or enable auto-render)."))
      s <- rle_shown$value()
      d <- s$df[s$df$sample %in% s$show, , drop = FALSE]
      d$sample <- droplevels(d$sample)
      validate(need(nrow(d) > 0, "No samples in the current 'Showing' selection."))
      .qc_rle_plot(d, dark(), length(levels(d$sample)), interactive = interactive,
                   palette = group_colours(input$rle_group, levels(d$group)))
    }
    dual_plot("diag_rle", rle_gg, n_elements = reactive({
      v <- rle_shown$value(); if (is.null(v)) 0L else sum(v$df$sample %in% v$show)
    }))

    dens_spec <- reactive({
      req(state$working, input$dens_group)
      dds <- state$working
      gmap <- group_map(input$dens_group, colnames(dds))
      df <- qc_expression_long(dds)
      df$group <- gmap[as.character(df$sample)]
      list(df = df, n = ncol(dds), show = showing_samples())
    })
    dens_shown <- deferred("dens_auto", "dens_render", dens_spec,
      sig = reactive(list(input$dens_group, showing_samples(), state$data_version,
                          .session_sig(input$dens_group))))
    output$dens_stale <- stale_note(dens_shown)
    dens_gg <- function(interactive) {
      validate(need(!is.null(dens_shown$value()), "Click Render (or enable auto-render)."))
      s <- dens_shown$value()
      d <- s$df[s$df$sample %in% s$show, , drop = FALSE]
      d$sample <- droplevels(d$sample)
      validate(need(nrow(d) > 0, "No samples in the current 'Showing' selection."))
      .qc_density_plot(d, dark(), length(levels(d$sample)), interactive = interactive,
                       palette = group_colours(input$dens_group, levels(d$group)))
    }
    dual_plot("diag_density", dens_gg, n_elements = reactive({
      v <- dens_shown$value(); if (is.null(v)) 0L else sum(v$df$sample %in% v$show)
    }))

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
    # Annotation choices (multi-select): the shared catalog (no gene, no "(none)";
    # a blank multi-select already means no annotation).
    cor_anno_choices <- function() aes_choices(aes_catalog(state), none = FALSE)
    output$cor_anno_ui <- renderUI({
      req(state$working)
      selectizeInput(ns("cor_anno"), "Annotate by (one or more)",
                     choices = cor_anno_choices(), selected = design_cols(), multiple = TRUE)
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
      samples <- colnames(state$working)
      valid <- vapply(aes_catalog(state), `[[`, "", "key")
      sel <- intersect(input$cor_anno, valid)
      # Snapshot the annotation *values* (deferred); colours resolve live at draw
      # so Palette edits recolour without a re-render.
      anno <- if (!length(sel)) NULL else {
        built <- Filter(Negate(is.null), lapply(sel, function(k) {
          r <- aes_resolve(state, k, samples)
          if (is.null(r)) NULL else list(key = k, name = r$label, values = r$values)
        }))
        if (!length(built)) NULL else list(
          keys = vapply(built, `[[`, "", "key"),
          df = data.frame(stats::setNames(lapply(built, `[[`, "values"),
                                          vapply(built, `[[`, "", "name")),
                          row.names = samples, check.names = FALSE, stringsAsFactors = FALSE))
      }
      list(cm = cm, anno = anno, n = ncol(state$working), method = method,
           show = showing_samples())
    })
    cor_shown <- deferred("cor_auto", "cor_render", cor_spec,
      sig = reactive(list(input$cor_method, input$cor_anno, showing_samples(), state$data_version,
                          if (any(c("__removal__", "__pool__") %in% input$cor_anno))
                            list(state$samp_pool, state$samp_flags) else NULL)))
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
      anno_df <- NULL; anno_col <- NULL
      if (!is.null(s$anno)) {
        anno_df <- s$anno$df[show, , drop = FALSE]
        # Resolve colours live (so a Palette edit recolours without re-render),
        # keyed by track name to match the snapshot value columns.
        res_list <- Filter(Negate(is.null), lapply(s$anno$keys, function(k) aes_resolve(state, k, show)))
        anno_col <- stats::setNames(lapply(res_list, aes_heatmap_col),
                                    vapply(res_list, function(r) r$label, ""))
      }
      ComplexHeatmap::draw(.qc_correlation_heatmap(
        cm, anno_df, dark_theme = dark(), n_samples = length(show), method = s$method,
        anno_col = anno_col, cor_config = state$palette$other$correlation))
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
      sig = reactive(list(input$wg_group, showing_samples(), state$data_version)))
    output$wg_stale <- stale_note(wg_shown)
    wg_gg <- function(interactive) {
      validate(need(!is.null(wg_shown$value()), "Click Render (or enable auto-render)."))
      s <- wg_shown$value()
      d <- s$df[s$df$sample %in% s$show, , drop = FALSE]
      .qc_within_group_plot(d, dark_theme = dark(), interactive = interactive,
                            palette = group_colours(input$wg_group, levels(d$group)))
    }
    dual_plot("wg_plot", wg_gg, n_elements = reactive({
      v <- wg_shown$value(); if (is.null(v)) 0L else sum(v$df$sample %in% v$show)
    }))

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
    # Observed-abundance assay for the dose-response. Read counts scale with molar
    # input *and* transcript length, and ERCCs span ~250-2000 nt, so a
    # length-normalized assay (TPM > FPKM) tracks the known molar concentration
    # more faithfully than CPM. Prefer them when present; fall back to CPM + note.
    spike_default_assay <- reactive({
      present <- SummarizedExperiment::assayNames(state$working)
      if ("TPM" %in% present) "TPM" else if ("FPKM" %in% present) "FPKM" else "CPM"
    })
    output$spike_assay_ui <- renderUI({
      req(state$working)
      present <- intersect(c("TPM", "FPKM", "CPM"), SummarizedExperiment::assayNames(state$working))
      tagList(
        selectInput(ns("spike_assay"), "Observed assay", choices = union(present, "CPM"),
                    selected = spike_default_assay()),
        if (identical(spike_default_assay(), "CPM"))
          helpText("Tip: add a length-normalized assay (TPM/FPKM) on the Assay tab - it tracks ",
                   "molar spike-in input more accurately than CPM for titration.")
      )
    })
    output$spike_group_ui <- group_box("spike_group", session = TRUE)
    output$spike_auto_ui  <- auto_box("spike_auto")

    spike_spec <- reactive({
      req(state$working, input$spike_source, input$spike_assay)
      source <- input$spike_source; assay <- input$spike_assay
      dr <- state_derive(state, "spike_dr", params = list(source = source, assay = assay),
                         expr = function() spike_dose_response(state$working, assay = assay, source = source))
      list(dr = dr, show = showing_samples())   # display subset (plots only)
    })
    # The group is applied at *render* (spike_group_col below), not baked into the
    # cached spec - so the group field (incl. session pool/flags) is deliberately
    # absent from the sig (like General QC's group): the plot recolours live, no
    # stale banner. RLE/density differ: they bake df$group into their spec.
    spike_shown <- deferred("spike_auto", "spike_render", spike_spec,
      sig = reactive(list(input$spike_source, input$spike_assay, showing_samples(),
                          state$data_version)))
    output$spike_stale <- stale_note(spike_shown)

    # group column for colouring (applied at render, not part of the cached compute)
    spike_group_col <- function(samples) {
      group_map(input$spike_group %||% default_group_col(),
                colnames(state$working))[as.character(samples)]
    }

    output$spike_msg <- renderUI({
      req(state$working)
      if (!length(spike_ids())) {
        return(tags$div(class = "alert alert-secondary py-2 small mb-2",
          "No spike-in features detected. Tag ERCC controls on the Feature tab to enable this view."))
      }
      v <- spike_shown$value(); if (is.null(v)) return(NULL)
      notes <- .spike_match_notes(v$dr, input$spike_source, ncol(state$working))
      if (!length(notes)) return(NULL)
      tags$div(class = "alert alert-warning py-2 small mb-2", lapply(notes, tags$div))
    })

    spike_dr_gg <- function(interactive) {
      validate(need(length(spike_ids()) > 0, "No spike-in features in this dataset."))
      validate(need(!is.null(spike_shown$value()), "Click Render (or enable auto-render)."))
      s <- spike_shown$value()
      validate(need(any(is.finite(s$dr$long$concentration)),
                    "No known concentrations resolved for these spike-ins - try another source."))
      long <- s$dr$long[s$dr$long$sample %in% s$show, , drop = FALSE]
      validate(need(nrow(long) > 0, "No samples in the current 'Showing' selection."))
      long$group <- spike_group_col(long$sample)
      .qc_spike_dr_plot(long, dark_theme = dark(), interactive = interactive,
                        palette = group_colours(input$spike_group %||% default_group_col(),
                                                levels(long$group)))
    }
    dual_plot("spike_dr_plot", spike_dr_gg, n_elements = reactive({
      v <- spike_shown$value(); if (is.null(v)) 0L else sum(v$dr$long$sample %in% v$show)
    }))

    spike_summary_gg <- function(interactive) {
      validate(need(!is.null(spike_shown$value()), "Click Render (or enable auto-render)."))
      s <- spike_shown$value()
      ps <- s$dr$per_sample[s$dr$per_sample$sample %in% s$show, , drop = FALSE]
      validate(need(nrow(ps) > 0, "No samples in the current 'Showing' selection."))
      ps$group <- spike_group_col(ps$sample)
      .qc_spike_summary_plot(ps, input$spike_metric %||% "r_squared", dark_theme = dark(),
                             interactive = interactive, sort = input$spike_sort %||% "none",
                             palette = group_colours(input$spike_group %||% default_group_col(),
                                                     levels(ps$group)))
    }
    dual_plot("spike_summary_plot", spike_summary_gg, n_elements = reactive({
      v <- spike_shown$value(); if (is.null(v)) 0L else sum(v$dr$per_sample$sample %in% v$show)
    }))

    output$spike_cv <- renderUI({
      v <- spike_shown$value(); req(v)
      pct <- v$dr$per_sample$pct_spike[is.finite(v$dr$per_sample$pct_spike)]
      if (length(pct) < 2L || mean(pct) == 0) return(NULL)
      tags$div(class = "small text-muted mb-2",
        sprintf("Spike-in fraction across samples: mean %.2f%%, CV %.0f%% (high CV = uneven spike input).",
                mean(pct), 100 * stats::sd(pct) / mean(pct)))
    })

    output$spike_table <- DT::renderDT({
      validate(need(!is.null(spike_shown$value()), "Click Render (or enable auto-render)."))
      ps <- spike_shown$value()$dr$per_sample
      dt_table(data.frame(
        Sample = ps$sample, `% spike` = round(ps$pct_spike, 2),
        Detected = ps$n_spike_detected, `Fit points` = ps$n_points,
        Slope = round(ps$slope, 3), `R-squared` = round(ps$r_squared, 3),
        `Lowest detected conc` = signif(ps$lod, 3), check.names = FALSE))
    })

    # --- Filtering: shared flags + removal pools ----------------------------
    # The sample pool is promoted to shared state so other pages (PCA) can read
    # it; this proxy keeps the reactiveVal-style get/set call sites unchanged.
    samp_pool <- function(x) {
      if (missing(x)) state$samp_pool %||% character(0)
      else { state$samp_pool <- x; invisible(x) }
    }
    feat_pool <- reactiveVal(character(0))

    # Advisory flags, cached on data_version + the rule inputs. Shared by the
    # Filtering tables and the General QC "Removal status" colour-by.
    # Canonical spike source / observed assay - shared with the Spike-in (ERCC)
    # tab so both hit one "spike_dr" cache entry. Defaults match that tab's
    # renderUI selections, so the key is identical even before it is first
    # visited (input$spike_* are then NULL).
    samp_spike_src   <- reactive(input$spike_source %||% (if (isTRUE(conc_col_ok())) "column" else "mix1"))
    samp_spike_assay <- reactive(input$spike_assay %||% spike_default_assay())
    # Per-sample dose-response summary feeding the spike filtering rules; NULL
    # when the dataset has no spike features (rules then stay disabled).
    samp_spike_summary <- reactive({
      req(state$working)
      if (!length(spike_ids())) return(NULL)
      src <- samp_spike_src(); asy <- samp_spike_assay()
      state_derive(state, "spike_dr", params = list(source = src, assay = asy),
                   expr = function() spike_dose_response(state$working, assay = asy, source = src))
    })

    samp_flags <- reactive({
      req(state$working)
      sp_ps <- { sp <- samp_spike_summary(); if (is.null(sp)) NULL else sp$per_sample }
      params <- list(group = input$samp_group, lib = .blank_na(input$samp_lib_min),
                     det = .blank_na(input$samp_detected_min),
                     mito = .blank_na(input$samp_mito_max), z = .blank_na(input$samp_wg_z),
                     src = samp_spike_src(), asy = samp_spike_assay(),
                     spmin = .blank_na(input$samp_spike_min), spmax = .blank_na(input$samp_spike_max),
                     sdet = .blank_na(input$samp_spike_detected_min),
                     r2 = .blank_na(input$samp_dose_r2_min),
                     slmin = .blank_na(input$samp_slope_min), slmax = .blank_na(input$samp_slope_max))
      state_derive(state, "samp_flags", params = params, expr = function() {
        flag_samples(state$working, group = input$samp_group,
                     lib_size_min = .blank_na(input$samp_lib_min),
                     detected_min = .blank_na(input$samp_detected_min),
                     pct_mito_max = .blank_na(input$samp_mito_max),
                     within_group_z = .blank_na(input$samp_wg_z),
                     spike = sp_ps,
                     pct_spike_min = .blank_na(input$samp_spike_min),
                     pct_spike_max = .blank_na(input$samp_spike_max),
                     min_spike_detected = .blank_na(input$samp_spike_detected_min),
                     dose_r2_min = .blank_na(input$samp_dose_r2_min),
                     dose_slope_min = .blank_na(input$samp_slope_min),
                     dose_slope_max = .blank_na(input$samp_slope_max))
      })
    })
    # Mirror the computed sample flags into shared state so PCA (and future plot
    # pages) can colour/shape by "Suggested removal" without re-deriving them.
    # Forces samp_flags() to evaluate even while the Filtering tab is hidden (the
    # threshold inputs stay rendered), so the value is ready when PCA reads it.
    observeEvent(samp_flags(), { state$samp_flags <- samp_flags() }, ignoreNULL = FALSE)

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
    # Does the dataset have spike-ins? Change-detected reactive so the filter
    # accordion below re-renders only when this flips, not on every edit.
    has_spike <- reactive(length(spike_ids()) > 0)
    # Build the sample-filter accordion server-side so the general and (optional)
    # spike-in panels are sibling panels in ONE accordion card. Stays rendered
    # while the tab is hidden so the threshold inputs always exist (the General
    # QC removal colour-by reads them).
    output$samp_filters_ui <- renderUI({
      nsf <- session$ns
      general <- bslib::accordion_panel(
        "Sample QC filters", icon = icon("filter"),
        uiOutput(nsf("samp_group_ui")),
        .qc_thr_input(nsf, "samp_lib_min", "Min library size (blank = off)", "samp_lib_auto",
                      tip = "Auto threshold: min. library size", value = NA, min = 0),
        .qc_thr_input(nsf, "samp_detected_min", "Min detected features (blank = off)",
                      "samp_detected_auto", tip = "Auto threshold: min. detected features",
                      value = NA, min = 0),
        .qc_thr_input(nsf, "samp_mito_max", "Max % mitochondrial (blank = off)", "samp_mito_auto",
                      tip = "Auto threshold: max. % mitochondrial", value = NA, min = 0, max = 100),
        numericInput(nsf("samp_wg_z"), "Within-group outlier z-cutoff (blank = off)",
                     value = 2, min = 0, step = 0.5),
        actionButton(nsf("samp_auto"), "Set auto thresholds",
                     icon = icon("wand-magic-sparkles"), class = "btn-sm btn-outline-primary w-100"))
      panels <- list(general)
      if (isTRUE(has_spike())) {
        spike <- bslib::accordion_panel(
          "Spike-in (ERCC) filters", icon = icon("vial"),
          helpText("Uses the concentration source & observed assay set on the Spike-in (ERCC) tab."),
          uiOutput(nsf("samp_spike_note")),
          .qc_thr_input(nsf, "samp_spike_min", "Min % spike-in (blank = off)", "samp_spike_min_auto",
                        tip = "Auto threshold: min. % spike-in (under-spiked)", value = NA, min = 0, max = 100),
          .qc_thr_input(nsf, "samp_spike_max", "Max % spike-in (blank = off)", "samp_spike_max_auto",
                        tip = "Auto threshold: max. % spike-in (over-spiked)", value = NA, min = 0, max = 100),
          numericInput(nsf("samp_spike_detected_min"), "Min detected spikes (blank = off)",
                       value = NA, min = 0, step = 1),
          .qc_thr_input(nsf, "samp_dose_r2_min", "Min dose-response R2 (blank = off)", "samp_dose_r2_auto",
                        tip = "Auto threshold: min. dose-response R-squared", value = NA, min = 0, max = 1, step = 0.05),
          tags$div(class = "d-flex gap-2",
            tags$div(class = "flex-grow-1",
              numericInput(nsf("samp_slope_min"), "Slope min", value = NA, step = 0.1)),
            tags$div(class = "flex-grow-1",
              numericInput(nsf("samp_slope_max"), "Slope max", value = NA, step = 0.1))),
          actionButton(nsf("samp_spike_auto"), "Set auto spike thresholds",
                       icon = icon("wand-magic-sparkles"), class = "btn-sm btn-outline-primary w-100"))
        panels <- c(panels, list(spike))
      }
      do.call(bslib::accordion, c(list(open = "Sample QC filters"), panels))
    })
    outputOptions(output, "samp_filters_ui", suspendWhenHidden = FALSE)
    # Reuse the QC tab's match-rate / no-detection warning, but only nag when a
    # spike rule is actually enabled (else it is noise on the Filtering tab).
    output$samp_spike_note <- renderUI({
      # Decide whether any spike rule is on first, so merely having spike-ins
      # present (no rule enabled) never forces the dose-response compute here.
      on <- any(vapply(list(input$samp_spike_min, input$samp_spike_max,
                            input$samp_spike_detected_min, input$samp_dose_r2_min,
                            input$samp_slope_min, input$samp_slope_max),
                       function(x) !is.null(.blank_na(x)), logical(1)))
      if (!on) return(NULL)
      sp <- samp_spike_summary(); if (is.null(sp)) return(NULL)
      notes <- .spike_match_notes(sp, samp_spike_src(), ncol(state$working))
      if (!length(notes)) return(NULL)
      tags$div(class = "alert alert-warning py-2 small mb-2", lapply(notes, tags$div))
    })
    output$feat_group_ui <- renderUI({
      req(state$working)
      cols <- colnames(SummarizedExperiment::colData(state$working))
      selectInput(ns("feat_group"), "filterByExpr grouping", choices = cols,
                  selected = default_group_col())
    })

    # Sample thresholds: fill data-driven defaults (a degenerate fence -> NA ->
    # blank input -> that rule stays disabled). The two accordion groups have
    # scoped Auto buttons: the general button fills lib/detected/mito; the spike
    # button fills the spike criteria. On load only the general fields are filled
    # (spike criteria stay blank / opt-in).
    fill_samp_general <- function() {
      req(state$working)
      th <- suggest_sample_thresholds(state$working)
      updateNumericInput(session, "samp_lib_min", value = th$lib_size_min)
      updateNumericInput(session, "samp_detected_min", value = th$detected_min)
      updateNumericInput(session, "samp_mito_max", value = th$pct_mito_max)
    }
    fill_samp_spike <- function() {
      sp <- samp_spike_summary(); if (is.null(sp)) return(invisible())
      th <- suggest_sample_thresholds(state$working, spike = sp$per_sample)
      updateNumericInput(session, "samp_spike_min", value = th$pct_spike_min)
      updateNumericInput(session, "samp_spike_max", value = th$pct_spike_max)
      updateNumericInput(session, "samp_dose_r2_min", value = th$dose_r2_min)
      updateNumericInput(session, "samp_slope_min", value = th$dose_slope_min)
      updateNumericInput(session, "samp_slope_max", value = th$dose_slope_max)
    }
    observeEvent(state$working, fill_samp_general())
    observeEvent(input$samp_auto, fill_samp_general())
    observeEvent(input$samp_spike_auto, fill_samp_spike())
    # Per-field Auto buttons fill just their own threshold.
    observeEvent(input$samp_lib_auto, updateNumericInput(
      session, "samp_lib_min", value = suggest_sample_thresholds(state$working)$lib_size_min))
    observeEvent(input$samp_detected_auto, updateNumericInput(
      session, "samp_detected_min", value = suggest_sample_thresholds(state$working)$detected_min))
    observeEvent(input$samp_mito_auto, updateNumericInput(
      session, "samp_mito_max", value = suggest_sample_thresholds(state$working)$pct_mito_max))
    # Per-field spike Auto buttons (data-driven fence; slope range is fixed).
    spike_auto_val <- function(field) {
      sp <- samp_spike_summary(); if (is.null(sp)) return(NA_real_)
      suggest_sample_thresholds(state$working, spike = sp$per_sample)[[field]]
    }
    observeEvent(input$samp_spike_min_auto, updateNumericInput(
      session, "samp_spike_min", value = spike_auto_val("pct_spike_min")))
    observeEvent(input$samp_spike_max_auto, updateNumericInput(
      session, "samp_spike_max", value = spike_auto_val("pct_spike_max")))
    observeEvent(input$samp_dose_r2_auto, updateNumericInput(
      session, "samp_dose_r2_min", value = spike_auto_val("dose_r2_min")))

    # Features default to adopting the suggestion (re-seeded when rules change);
    # samples stay opt-in (highlight only). Both pools reset on any data change.
    observeEvent(feat_flags(), {
      feat_pool(feat_flags()$feature_id[feat_flags()$suggested_drop])
    })
    observeEvent(state$data_version, samp_pool(character(0)))

    # Tables: a read-only boolean "Suggested Removal" column + a separate boolean
    # "In Removal Pool" column; native row-selection is transient staging the
    # pool buttons act on.
    samp_display <- function(fl, pool) {
      df <- data.frame(
        Sample = fl$sample, `Library size (M)` = round(fl$library_size / 1e6, 3),
        Detected = fl$detected, `% mito` = round(fl$pct_mito, 2),
        `Within-grp corr` = round(fl$within_group_corr, 3),
        check.names = FALSE, stringsAsFactors = FALSE)
      # Spike columns only when the dataset has spike-ins (finite n_spike_detected,
      # which is 0+ when spikes exist and NA when there are none).
      if (any(is.finite(fl$n_spike_detected))) {
        df$`% spike`         <- round(fl$pct_spike, 2)
        df$`Spikes detected` <- fl$n_spike_detected
        df$`Dose R2`         <- round(fl$dose_r2, 3)
        df$`Dose slope`      <- round(fl$dose_slope, 2)
      }
      df$`Suggested Removal` <- fl$flagged
      df$Reason              <- fl$reason
      df$`In Removal Pool`   <- fl$sample %in% pool
      df
    }
    feat_display <- function(fl, pool) data.frame(
      Feature = fl$feature_id, Class = fl$feature_class,
      `Total count` = round(fl$total_count), Detected = fl$n_detected,
      `Mean log` = round(fl$mean_logcounts, 2),
      `Suggested Removal` = fl$suggested_drop, Reason = fl$reason,
      `In Removal Pool` = fl$feature_id %in% pool,
      check.names = FALSE, stringsAsFactors = FALSE)

    # Danger-red TRUEs in the boolean removal columns (survives replaceData since
    # it is baked into the column definitions, not the data).
    color_removal <- function(dt) DT::formatStyle(
      dt, c("Suggested Removal", "In Removal Pool"),
      color      = DT::styleEqual(TRUE, "var(--bs-danger)", default = NULL),
      fontWeight = DT::styleEqual(TRUE, "bold", default = NULL))
    output$samp_tbl <- DT::renderDT({
      validate(need(!is.null(state$working), "No dataset loaded."))
      color_removal(dt_table(samp_display(samp_flags(), isolate(samp_pool())),
                             selection = list(mode = "multiple")))
    }, server = TRUE)
    output$feat_tbl <- DT::renderDT({
      validate(need(!is.null(state$working), "No dataset loaded."))
      color_removal(dt_table(feat_display(feat_flags(), isolate(feat_pool())),
                             selection = list(mode = "multiple")))
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
    feat_density_gg <- function(interactive) {
      req(state$working)
      keep <- setdiff(rownames(state$working), feat_pool())
      .qc_filter_density_plot(qc_filter_density(state$working, keep), dark_theme = dark(),
                              interactive = interactive)
    }
    # Element estimate without recomputing the density frame: the before+after
    # curves span at most 2 x (endogenous features) x samples.
    dual_plot("feat_density", feat_density_gg, n_elements = reactive({
      req(state$working)
      2L * sum(.endogenous_mask(state$working)) * ncol(state$working)
    }))

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
