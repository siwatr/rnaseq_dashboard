# Page 7: Expression. A gene-expression browsing surface (renamed from the old
# "Heatmap" stub). A navset_card_tab with two tabs:
#   Single genes - one feature at a time as a layered violin -> box -> dots overlay,
#                  grouped by a colData variable.
#   Gene sets    - a navset_card_pill of:
#                    Aggregate expression - the SAME layered overlay, but the y value
#                      is a per-sample gene-set score (mean/median across the set's
#                      genes; z-scored by default so highly-expressed genes do not
#                      dominate). Source = a saved set or a quick (uncommitted) search.
#                    Heatmap - a ComplexHeatmap over the set (built in P7c).
#
# The single-gene and aggregate pills share the plot machinery (gene/set source,
# expr-value control, plot engine + deferred render, "Showing" subset, aes_helpers
# colour) via `.expr_dist_server()`: the (expensive) value matrix is cached in
# `derived` behind the Render button, keyed on the assay; the per-gene lookup / the
# set aggregation, grouping, colour, and geom toggles are cheap live re-plots.

# --- Gene-set heatmap (P7c) tuning knobs (option()s, like the distribution
# guards) -- row (genes) and column (samples) label defaults resolve separately,
# and the dendrogram auto-hides above its own threshold.
.hm_row_label_max <- function() getOption("ddsdashboard.heatmap_row_label_max", 50L)
.hm_col_label_max <- function() getOption("ddsdashboard.heatmap_col_label_max", 30L)
.hm_dend_max      <- function() getOption("ddsdashboard.heatmap_dend_max", 100L)

# Build the ComplexHeatmap colour function from a continuous-palette config. In
# z-score mode with both anchors blank the divergent ramp is centred at 0
# (symmetric limits); otherwise the anchors / data range apply. Returns NULL when
# circlize is unavailable (callers validate()).
.hm_col_fun <- function(mat, ramp_cfg, zscored) {
  if (!requireNamespace("circlize", quietly = TRUE)) return(NULL)
  vals <- as.numeric(mat)
  mn <- ramp_cfg$min; mx <- ramp_cfg$max
  blank <- function(a) is.null(a) || !nzchar(trimws(as.character(a)))
  if (isTRUE(zscored) && blank(mn) && blank(mx)) {
    lim <- expr_symmetric_limits(vals); mn <- lim[1]; mx <- lim[2]
  }
  palette_colorramp2(ramp_cfg$name %||% "RColorBrewer: RdBu", values = vals,
                     min = if (blank(mn)) NULL else mn,
                     max = if (blank(mx)) NULL else mx,
                     custom = ramp_cfg$custom, reverse = isTRUE(ramp_cfg$reverse))
}

# Default ramp config by scaling mode: divergent RdBu (centred at 0) for z-scores,
# sequential viridis for raw values.
.hm_ramp_default <- function(zscore = TRUE) {
  if (isTRUE(zscore))
    continuous_palette_default(name = "RColorBrewer: RdBu", reverse = TRUE)
  else continuous_palette_default(name = "viridis: viridis")
}

# Assemble the Heatmap object from a rendered snapshot `s` + live annotation +
# a live colour function. `s` carries the matrix, resolved labels/marks,
# cluster/dend flags, and legend title; `anno` is a HeatmapAnnotation (or NULL)
# and `col_fun` the value ramp -- both resolved live at draw (a value ramp is a
# display aesthetic, so it must not enter the deferred sig / stale the plot).
.hm_build <- function(s, anno, col_fun) {
  args <- list(
    matrix = s$mat, name = s$legend, col = col_fun,
    cluster_rows = isTRUE(s$cluster_rows), cluster_columns = isTRUE(s$cluster_cols),
    show_row_dend = isTRUE(s$show_row_dend), show_column_dend = isTRUE(s$show_col_dend),
    show_row_names = identical(s$row_mode, "all"),
    show_column_names = identical(s$col_mode, "all"),
    row_names_gp = grid::gpar(fontsize = 8),
    column_names_gp = grid::gpar(fontsize = 9),
    top_annotation = anno,
    heatmap_legend_param = list(title = s$legend))
  if (identical(s$row_mode, "all")) args$row_labels <- s$row_labels
  if (identical(s$col_mode, "all")) args$column_labels <- s$col_labels
  # "Selected": names off, mark just the searched rows/columns via anno_mark.
  if (identical(s$row_mode, "selected") && length(s$row_mark_at))
    args$right_annotation <- ComplexHeatmap::rowAnnotation(
      mark = ComplexHeatmap::anno_mark(at = s$row_mark_at, labels = s$row_mark_lab,
                                       labels_gp = grid::gpar(fontsize = 8)))
  if (identical(s$col_mode, "selected") && length(s$col_mark_at))
    args$bottom_annotation <- ComplexHeatmap::columnAnnotation(
      mark = ComplexHeatmap::anno_mark(at = s$col_mark_at, labels = s$col_mark_lab,
                                       which = "column",
                                       labels_gp = grid::gpar(fontsize = 8)))
  do.call(ComplexHeatmap::Heatmap, args)
}

# Colour/fill scale + layered plot for one value vector. `df` carries
# sample/group/value (+ optional `colour`). Distribution geoms fill by a discrete
# colour attribute; dots colour by it. ggbeeswarm draws the dots (geom_jitter
# fallback). Shared by the single-gene and gene-set-aggregate pills.
#' @importFrom ggplot2 .data
.expr_single_plot <- function(df, x_lab, y_lab, title, subtitle = NULL,
                              show_violin = TRUE, show_box = TRUE, show_dots = TRUE,
                              disc_colour = FALSE, colour_lab = NULL,
                              fill_scale = NULL, colour_scale = NULL,
                              violin_width = 0.9, violin_alpha = 0.35,
                              box_width = 0.18, box_alpha = 0.7,
                              dot_method = "quasirandom", dot_size = 1.9,
                              dot_alpha = 0.9, dot_width = 0.4, dot_cex = 1,
                              y_range = NULL, legend_pos = "right",
                              dark_theme = FALSE, interactive = FALSE) {
  # y-axis clamp: pull out-of-range points to the nearest limit (drawn as
  # triangles), and zoom every layer to the range via coord_cartesian (so the
  # violin/box distributions are clipped, not distorted). Points use a clamped y;
  # the distributions keep the true values (coord crops them at draw).
  clamp_on <- !is.null(y_range) && any(is.finite(y_range))
  if (clamp_on) {
    cl <- de_clamp(df$value, y_range[1], y_range[2])
    df$value_dot <- cl$value
    df$oob <- factor(ifelse(cl$clamped, "clamped", "in range"),
                     levels = c("in range", "clamped"))
  }
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$group, y = .data$value))
  fill_map   <- if (disc_colour) ggplot2::aes(fill = .data$colour) else NULL
  colour_map <- if (!is.null(df$colour)) ggplot2::aes(colour = .data$colour) else NULL

  if (isTRUE(show_violin))
    p <- p + ggplot2::geom_violin(mapping = fill_map, alpha = violin_alpha,
                                  colour = NA, scale = "width", width = violin_width,
                                  show.legend = disc_colour)
  if (isTRUE(show_box))
    p <- p + ggplot2::geom_boxplot(mapping = fill_map, width = box_width,
                                   alpha = box_alpha, outlier.shape = NA,
                                   show.legend = FALSE)
  if (isTRUE(show_dots)) {
    dot_aes <- colour_map
    if (clamp_on) {                         # dots use the clamped y + triangle shape
      cmap <- ggplot2::aes(y = .data$value_dot, shape = .data$oob)
      dot_aes <- if (is.null(dot_aes)) cmap else utils::modifyList(dot_aes, cmap)
    }
    if (interactive) {
      dot_aes <- if (is.null(dot_aes)) ggplot2::aes(text = .data$text)
                 else utils::modifyList(dot_aes, ggplot2::aes(text = .data$text))
    }
    dot_args <- list(mapping = dot_aes, size = dot_size, alpha = dot_alpha)
    # cex (point spacing) is native to the beeswarm layout; width (max spread) to
    # quasirandom / jitter. Keeping each control on its own method avoids the
    # ggbeeswarm "duplicated size" warning that cex + a grouping aesthetic trips.
    # Jitter is an explicit choice AND the fallback when ggbeeswarm is absent.
    have_bees <- requireNamespace("ggbeeswarm", quietly = TRUE)
    p <- p + if (have_bees && identical(dot_method, "beeswarm")) {
      do.call(ggbeeswarm::geom_beeswarm, c(dot_args, list(cex = dot_cex)))
    } else if (have_bees && identical(dot_method, "quasirandom")) {
      do.call(ggbeeswarm::geom_quasirandom, c(dot_args, list(width = dot_width)))
    } else {
      do.call(ggplot2::geom_jitter, c(dot_args, list(width = dot_width, height = 0)))
    }
  }

  lab_args <- list(x = x_lab, y = y_lab, title = title, subtitle = subtitle)
  if (disc_colour && (show_violin || show_box)) lab_args$fill <- colour_lab
  if (!is.null(df$colour) && show_dots) lab_args$colour <- colour_lab
  p <- p +
    do.call(ggplot2::labs, lab_args) +
    .plot_theme(dark_theme) +
    ggplot2::theme(legend.position = legend_pos %||% "right",
                   axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
  if (!is.null(fill_scale))   p <- p + fill_scale
  if (!is.null(colour_scale)) p <- p + colour_scale
  if (clamp_on) {
    if (isTRUE(show_dots))
      p <- p + ggplot2::scale_shape_manual(values = c("in range" = 16, "clamped" = 17),
                                           drop = FALSE, guide = "none")
    p <- p + ggplot2::coord_cartesian(ylim = y_range)
  }
  p
}

# --- Shared sidebar UI pieces (parameterized by a prefix so the single-gene and
# gene-set pills get distinct, non-colliding input ids) ----------------------

# The auto-render checkbox + Render button row.
.expr_render_bar <- function(ns, prefix) {
  nid <- function(x) ns(paste0(prefix, x))
  tags$div(class = "d-flex align-items-center gap-2 mb-1",
           tags$div(class = "flex-grow-1", uiOutput(nid("auto_ui"))),
           actionButton(nid("render"), "Render", icon = icon("play"),
                        class = "btn-primary"))
}

# The "Grouping & colour" accordion panel.
.expr_grouping_panel <- function(ns, prefix) {
  nid <- function(x) ns(paste0(prefix, x))
  bslib::accordion_panel(
    "Grouping & colour", icon = icon("layer-group"),
    uiOutput(nid("x_group_ui")),
    uiOutput(nid("colour_ui")))
}

# The "Plot elements" + "Layer styling" accordion panels (returned as a list so
# the caller can splice them into an accordion via do.call).
.expr_style_panels <- function(ns, prefix) {
  nid <- function(x) ns(paste0(prefix, x))
  list(
    bslib::accordion_panel(
      "Plot elements", icon = icon("chart-simple"),
      uiOutput(nid("geom_ui"))),
    bslib::accordion_panel(
      "Layer styling", icon = icon("sliders"),
      tags$h6(class = "fs-4 mb-1", "Violin"),
      sliderInput(nid("violin_width"), "Width", 0.2, 1.5, 0.9, 0.05),
      sliderInput(nid("violin_alpha"), "Opacity", 0, 1, 0.35, 0.05),
      tags$hr(class = "my-2"),
      tags$h6(class = "fs-4 mb-1", "Boxplot"),
      sliderInput(nid("box_width"), "Width", 0.05, 0.8, 0.18, 0.01),
      sliderInput(nid("box_alpha"), "Opacity", 0, 1, 0.7, 0.05),
      tags$hr(class = "my-2"),
      tags$h6(class = "fs-4 mb-1", "Data points"),
      sliderInput(nid("dot_size"), "Size", 0.25, 5, 1.9, 0.25),
      sliderInput(nid("dot_alpha"), "Opacity", 0, 1, 0.9, 0.05),
      uiOutput(nid("dot_method_ui")),
      conditionalPanel(
        sprintf("input['%s'] != 'beeswarm'", nid("dot_method")),
        sliderInput(nid("dot_width"), "Spread (width)", 0, 0.5, 0.4, 0.05)),
      conditionalPanel(
        sprintf("input['%s'] == 'beeswarm'", nid("dot_method")),
        sliderInput(nid("dot_cex"), "Point spacing (cex)", 0.2, 3, 1, 0.1))),
    bslib::accordion_panel(
      "Axis limits", icon = icon("up-right-and-down-left-from-center"),
      helpText(class = "small text-muted",
               "Blank = auto. Points outside the range draw as triangles at the edge."),
      tags$div(class = "small text-muted mb-1", "y-axis (expression value)"),
      bslib::layout_columns(
        col_widths = c(6, 6),
        numericInput(nid("ylim_min"), "min", value = NA),
        numericInput(nid("ylim_max"), "max", value = NA))))
}

mod_expression_ui <- function(id) {
  ns <- NS(id)

  # --- Single genes -------------------------------------------------------
  single_tab <- bslib::nav_panel(
    "Single genes",
    bslib::layout_sidebar(
      sidebar = bslib::sidebar(
        title = tags$h3("Single gene", class = "fs-6 mb-0"), width = 340,
        .expr_render_bar(ns, ""),
        do.call(bslib::accordion, c(
          list(open = c("Gene", "Grouping & colour"),
               bslib::accordion_panel(
                 "Gene", icon = icon("magnifying-glass"),
                 gene_search_ui(ns, "gene", multiple = FALSE, dup_toggle = TRUE,
                                placeholder = "e.g. Gene1"),
                 expr_value_ui(ns, "val")),
               .expr_grouping_panel(ns, "")),
          .expr_style_panels(ns, ""))),
        plot_subset_ui(ns, "expr")
      ),
      bslib::card(
        bslib::card_header(tags$h3("Single-gene expression", class = "fs-6 mb-0")),
        uiOutput(ns("gene_stale")),
        shinycssloaders::withSpinner(.plot_dual(ns("gene_container")),
                                     proxy.height = "400px")
      )
    )
  )

  # --- Gene sets > Aggregate expression -----------------------------------
  aggregate_pill <- bslib::nav_panel(
    "Aggregate expression",
    bslib::layout_sidebar(
      sidebar = bslib::sidebar(
        title = tags$h3("Gene-set score", class = "fs-6 mb-0"), width = 340,
        .expr_render_bar(ns, "set_"),
        do.call(bslib::accordion, c(
          list(open = c("Gene set", "Aggregation"),
               bslib::accordion_panel(
                 "Gene set", icon = icon("magnifying-glass"),
                 radioButtons(ns("set_source"), "Source",
                              c("Saved set" = "saved", "Quick search" = "search"),
                              inline = TRUE),
                 conditionalPanel(
                   sprintf("input['%s'] == 'saved'", ns("set_source")),
                   uiOutput(ns("set_pick_ui"))),
                 conditionalPanel(
                   sprintf("input['%s'] == 'search'", ns("set_source")),
                   gene_search_ui(ns, "setsearch", multiple = TRUE,
                                  search_modes = c("exact", "contains", "regex"),
                                  placeholder = "e.g. Actb, Gapdh, ENSG..."))),
               bslib::accordion_panel(
                 "Aggregation", icon = icon("calculator"),
                 expr_value_ui(ns, "set_val"),
                 radioButtons(ns("set_method"), "Average across genes",
                              c("Mean" = "mean", "Median" = "median"), inline = TRUE),
                 bslib::input_switch(
                   ns("set_zscore"),
                   bslib::tooltip(
                     tags$span("Z-score each gene"),
                     "Standardize every gene across samples before averaging, so highly-expressed genes don't dominate the score. Recommended."),
                   value = TRUE),
                 bslib::input_switch(
                   ns("set_only_expr"),
                   bslib::tooltip(
                     tags$span("Only genes with expression"),
                     "Drop genes with zero counts in every sample before averaging."),
                   value = TRUE)),
               .expr_grouping_panel(ns, "set_")),
          .expr_style_panels(ns, "set_"))),
        plot_subset_ui(ns, "expr_set")
      ),
      bslib::card(
        bslib::card_header(tags$h3("Gene-set aggregate expression", class = "fs-6 mb-0")),
        uiOutput(ns("geneset_stale")),
        shinycssloaders::withSpinner(.plot_dual(ns("geneset_container")),
                                     proxy.height = "400px")
      )
    )
  )

  # --- Gene sets > Heatmap (P7c) ------------------------------------------
  # A ComplexHeatmap over a named set. No auto-render (a heatmap can be slow):
  # Render-only + a spinner. Everything -- matrix, labels, clustering, ramp,
  # annotation -- is snapshotted on Render (a single static plot has no cheap live
  # layer), so a settings change shows a stale banner until the next Render.
  hm_label_panel <- function(axis, title, icon_name) {
    p <- function(x) paste0("hm_", axis, "_", x)
    is_row <- identical(axis, "row")
    bslib::accordion_panel(
      title, icon = icon(icon_name),
      uiOutput(ns(p("source_ui"))),
      radioButtons(ns(p("mode")), "Show labels",
                   c("Auto (by size)" = "auto", "All" = "all",
                     "Selected" = "selected", "None" = "none"),
                   selected = "auto", inline = TRUE),
      conditionalPanel(
        sprintf("input['%s'] == 'selected'", ns(p("mode"))),
        if (is_row)
          gene_search_ui(ns, "hmrowsel", multiple = TRUE,
                         search_modes = c("exact", "contains", "regex"),
                         placeholder = "Genes to label")
        else uiOutput(ns("hm_col_sel_ui"))))
  }
  heatmap_pill <- bslib::nav_panel(
    "Heatmap",
    bslib::layout_sidebar(
      sidebar = bslib::sidebar(
        title = tags$h3("Gene-set heatmap", class = "fs-6 mb-0"), width = 350,
        tags$div(class = "d-flex justify-content-end mb-1",
                 actionButton(ns("hm_render"), "Render", icon = icon("play"),
                              class = "btn-primary")),
        bslib::accordion(
          open = "Gene set",
          bslib::accordion_panel(
            "Gene set", icon = icon("magnifying-glass"),
            radioButtons(ns("hm_source"), "Source",
                         c("Saved set" = "saved", "Quick search" = "search"),
                         inline = TRUE),
            conditionalPanel(
              sprintf("input['%s'] == 'saved'", ns("hm_source")),
              uiOutput(ns("hm_pick_ui"))),
            conditionalPanel(
              sprintf("input['%s'] == 'search'", ns("hm_source")),
              gene_search_ui(ns, "hmsearch", multiple = TRUE,
                             search_modes = c("exact", "contains", "regex"),
                             placeholder = "e.g. Actb, Gapdh, ENSG..."))),
          bslib::accordion_panel(
            "Expression value", icon = icon("calculator"),
            expr_value_ui(ns, "hm_val"),
            bslib::input_switch(
              ns("hm_zscore"),
              bslib::tooltip(tags$span("Z-score each gene"),
                "Standardize every gene across samples (centre 0, sd 1) so genes on different scales are comparable. Recommended."),
              value = TRUE),
            bslib::input_switch(
              ns("hm_only_expr"),
              bslib::tooltip(tags$span("Only genes with expression"),
                "Drop genes with zero counts in every sample."),
              value = TRUE)),
          bslib::accordion_panel(
            "Colour", icon = icon("palette"),
            radioButtons(ns("hm_ramp_src"), "Colour scale",
                         c("Custom" = "custom", "From Palette page (assay)" = "palette"),
                         selected = "custom", inline = TRUE),
            conditionalPanel(
              sprintf("input['%s'] == 'palette'", ns("hm_ramp_src")),
              helpText(class = "small text-muted",
                       "Uses the assay's continuous ramp from the Palette page. Applies to raw values; z-scored data always uses the custom divergent ramp.")),
            conditionalPanel(
              sprintf("input['%s'] == 'custom'", ns("hm_ramp_src")),
              continuous_palette_ui(ns, "hm_ramp", .hm_ramp_default(TRUE)))),
          hm_label_panel("row", "Row labels", "align-left"),
          hm_label_panel("col", "Column labels", "align-left"),
          bslib::accordion_panel(
            "Sample annotation", icon = icon("table-cells"),
            uiOutput(ns("hm_anno_ui")),
            actionButton(ns("hm_clear_anno"), "Clear annotation",
                         icon = icon("xmark"), class = "btn-sm btn-outline-secondary")),
          bslib::accordion_panel(
            "Heatmap elements", icon = icon("sitemap"),
            bslib::input_switch(ns("hm_cluster_rows"), "Cluster rows", value = TRUE),
            conditionalPanel(
              sprintf("input['%s']", ns("hm_cluster_rows")),
              selectInput(ns("hm_row_dend"), "Row dendrogram",
                          c("Auto (by size)" = "auto", "Show" = "show", "Hide" = "hide"),
                          selected = "auto")),
            tags$hr(class = "my-2"),
            bslib::input_switch(ns("hm_cluster_cols"), "Cluster columns", value = TRUE),
            conditionalPanel(
              sprintf("input['%s']", ns("hm_cluster_cols")),
              selectInput(ns("hm_col_dend"), "Column dendrogram",
                          c("Auto (by size)" = "auto", "Show" = "show", "Hide" = "hide"),
                          selected = "auto")))),
        plot_subset_ui(ns, "hm")
      ),
      bslib::card(
        bslib::card_header(tags$h3("Gene-set expression heatmap", class = "fs-6 mb-0")),
        uiOutput(ns("hm_stale")),
        uiOutput(ns("hm_info")),
        shinycssloaders::withSpinner(plotOutput(ns("hm_plot"), height = "560px"),
                                     proxy.height = "560px")
      )
    )
  )

  bslib::navset_card_tab(
    id = ns("tabs"),
    title = tags$h2("Expression", class = "fs-6 mb-0"),
    single_tab,
    bslib::nav_panel(
      "Gene sets",
      bslib::navset_card_pill(aggregate_pill, heatmap_pill)
    )
  )
}

# Shared server for a distribution-overlay pill (single gene OR gene-set score).
# Wires grouping / colour / geom availability + toggles / layer styling, the
# deferred (cached) value matrix, and the dual plot -- all keyed off `cfg$prefix`
# so two instances coexist in one module. `cfg$reduce_factory(val_spec)` returns a
# `function(mm)` turning the cached value matrix into the per-sample vector to plot
# (`list(values, y_lab, title, subtitle, hover)`), or aborting via validate().
# Returns the reactive handles the page (and tests) need.
.expr_dist_server <- function(input, output, session, state, cfg) {
  ns <- session$ns
  P <- cfg$prefix
  pid <- function(x) paste0(P, x)                 # module-relative input/output id
  dark <- cfg$dark; eng <- cfg$engine
  showing_samples <- cfg$showing_samples

  val_spec <- expr_value_server(input, output, session, state, cfg$val_suffix,
                                include_vst = TRUE, include_norm = TRUE,
                                default_fn = expr_default_assay)
  reduce <- cfg$reduce_factory(val_spec)

  coldata <- function() as.data.frame(SummarizedExperiment::colData(state$working))
  discrete_cols <- function() {
    cd <- coldata()
    if (!ncol(cd)) return(character(0))
    names(cd)[vapply(cd, function(x) is.factor(x) || is.character(x) || is.logical(x),
                     logical(1))]
  }
  default_x <- function() {
    disc <- discrete_cols()
    pv <- primary_design_var(state$working)      # variable of interest = last design term
    if (!is.na(pv) && pv %in% disc) return(pv)
    if (length(disc)) disc[1] else "__none__"
  }

  output[[pid("x_group_ui")]] <- renderUI({
    req(state$working)
    selectInput(ns(pid("x_group")), "Group by (x-axis)",
                choices = group_field_choices(discrete_cols(), none = FALSE),
                selected = default_x())
  })
  output[[pid("colour_ui")]] <- renderUI({
    req(state$working)
    selectInput(ns(pid("colour_by")), "Colour by",
                choices = aes_choices(aes_catalog(state), none = TRUE),
                selected = default_x())
  })
  # The grouping / colour selectInputs live in an accordion panel that may be
  # collapsed on load (it is on the aggregate pill). Force these outputs to
  # initialize while hidden so `input$<x_group>` exists on first entry and the
  # plot renders without the user having to expand the panel.
  outputOptions(output, pid("x_group_ui"), suspendWhenHidden = FALSE)
  outputOptions(output, pid("colour_ui"), suspendWhenHidden = FALSE)
  output[[pid("auto_ui")]] <- renderUI({
    req(state$working)
    bslib::input_switch(ns(pid("auto")), "Auto-render", value = TRUE)
  })
  output[[pid("dot_method_ui")]] <- renderUI({
    if (!requireNamespace("ggbeeswarm", quietly = TRUE))
      return(helpText(class = "small text-muted",
                      "ggbeeswarm not installed - using jittered points."))
    selectInput(ns(pid("dot_method")), "Layout",
                c("Quasirandom" = "quasirandom", "Beeswarm" = "beeswarm",
                  "Jitter" = "jitter"),
                selected = "quasirandom")
  })

  group_sizes <- reactive({
    req(state$working, input[[pid("x_group")]])
    cd <- coldata()
    if (!input[[pid("x_group")]] %in% names(cd)) return(integer(0))
    shown <- intersect(rownames(cd), showing_samples())
    g <- cd[shown, input[[pid("x_group")]]]
    tab <- table(g[!is.na(g)])
    as.integer(tab[tab > 0])            # drop empty groups (e.g. filtered by "Showing")
  })
  geom_avail <- reactive(expr_geom_availability(group_sizes()))

  output[[pid("geom_ui")]] <- renderUI({
    av <- geom_avail()
    dist <- isTRUE(av$dist_shown)
    keep <- function(x, default) {
      v <- shiny::isolate(input[[pid(x)]]); if (is.null(v)) default else isTRUE(v)
    }
    tagList(
      if (dist) bslib::input_switch(ns(pid("show_violin")), "Violin",
                                    value = keep("show_violin", TRUE))
      else helpText(class = "small text-muted",
                    "Violin / box hidden: groups are too small to summarize."),
      if (dist) bslib::input_switch(ns(pid("show_box")), "Boxplot",
                                    value = keep("show_box", TRUE)),
      if (isTRUE(av$dots_allowed))
        bslib::input_switch(ns(pid("show_dots")), "Data points",
                            value = keep("show_dots", av$dots_default))
      else helpText(class = "small text-muted",
                    sprintf("Data points hidden: a group has %d+ samples (overplotting).",
                            av$n_max))
    )
  })
  on_flag <- function(x, allowed, default) {
    if (!isTRUE(allowed)) return(FALSE)
    v <- input[[pid(x)]]; if (is.null(v)) default else isTRUE(v)
  }

  # y-axis clamp (blank = auto); a live display input, not part of the render gate.
  y_range <- reactive({
    lo <- suppressWarnings(as.numeric(input[[pid("ylim_min")]] %||% NA))
    hi <- suppressWarnings(as.numeric(input[[pid("ylim_max")]] %||% NA))
    if (is.na(lo) && is.na(hi)) NULL else c(lo, hi)
  })

  colour_resolve <- function(samples) {
    sel <- input[[pid("colour_by")]] %||% "__none__"
    if (identical(sel, "__none__")) return(NULL)
    res <- aes_resolve(state, sel, samples)
    if (is.null(res)) return(NULL)
    disc <- identical(res$kind, "discrete")
    list(values = res$values, discrete = disc, lab = res$label,
         fill_scale = if (disc) aes_ggplot_scale(res, "fill") else NULL,
         colour_scale = aes_ggplot_scale(res, "colour"))
  }

  # The value matrix (cached in `derived`, keyed on the assay) AND the source
  # reduction to per-sample values are gated together behind Render: changing the
  # gene / gene set / value / aggregation waits for a render (or auto-render).
  # Grouping, colour, and geom toggles are cheap live re-plots on top of the gated
  # result. `reduce` returns list(ok, ...) instead of validate()-ing (the deferred
  # spec swallows validate()); build_gg surfaces `msg` when !ok.
  value_spec <- reactive({
    req(state$working)
    assay <- val_spec()$assay
    mm <- state_derive(state, cfg$derive_key, params = list(assay),
                       expr = function() expr_value_matrix(state$working, assay))
    list(mat_cols = colnames(mm$mat), red = reduce(mm))
  })
  out <- eng$deferred(pid("auto"), pid("render"), value_spec,
    sig = reactive(list(val_spec(), cfg$source_sig(), state$data_version)))
  output[[cfg$stale_id]] <- eng$stale_note(out)

  build_gg <- function(interactive) {
    got <- out$value()
    validate(need(!is.null(got),
                  "Click Render (or enable auto-render) to draw the plot."))
    red <- got$red
    validate(need(isTRUE(red$ok), red$msg %||% "Nothing to plot."))
    req(input[[pid("x_group")]])
    cd <- coldata()
    validate(need(input[[pid("x_group")]] %in% names(cd), "Pick a grouping variable."))
    shown <- intersect(names(red$values), showing_samples())
    shown <- intersect(shown, rownames(cd))
    validate(need(length(shown) > 0, "No samples in the current 'Showing' selection."))
    groups <- cd[shown, input[[pid("x_group")]]]
    cr <- colour_resolve(shown)
    df <- expr_long_frame(red$values[shown], groups, shown,
                          colour = if (!is.null(cr)) cr$values else NULL)
    validate(need(nrow(df) > 0, "No samples with a non-missing group."))

    av <- geom_avail()
    sv <- on_flag("show_violin", av$dist_shown, TRUE)
    sb <- on_flag("show_box",    av$dist_shown, TRUE)
    sd <- on_flag("show_dots",   av$dots_allowed, av$dots_default)
    validate(need(sv || sb || sd, "Enable at least one plot element."))
    disc <- !is.null(cr) && isTRUE(cr$discrete)
    if (interactive) df$text <- red$hover(df)

    .expr_single_plot(
      df, x_lab = input[[pid("x_group")]], y_lab = red$y_lab,
      title = red$title, subtitle = red$subtitle,
      show_violin = sv, show_box = sb, show_dots = sd,
      disc_colour = disc, colour_lab = if (!is.null(cr)) cr$lab else NULL,
      fill_scale = if (!is.null(cr)) cr$fill_scale else NULL,
      colour_scale = if (!is.null(cr)) cr$colour_scale else NULL,
      violin_width = input[[pid("violin_width")]] %||% 0.9,
      violin_alpha = input[[pid("violin_alpha")]] %||% 0.35,
      box_width = input[[pid("box_width")]] %||% 0.18,
      box_alpha = input[[pid("box_alpha")]] %||% 0.7,
      dot_method = input[[pid("dot_method")]] %||% "quasirandom",
      dot_size = input[[pid("dot_size")]] %||% 1.9,
      dot_alpha = input[[pid("dot_alpha")]] %||% 0.9,
      dot_width = input[[pid("dot_width")]] %||% 0.4,
      dot_cex = input[[pid("dot_cex")]] %||% 1,
      y_range = y_range(), dark_theme = dark(), interactive = interactive)
  }
  eng$dual_plot(cfg$plot_id, build_gg, n_elements = reactive({
    got <- out$value(); if (is.null(got)) 0L
    else length(intersect(got$mat_cols, showing_samples()))
  }), height = "440px")

  list(out = out, group_sizes = group_sizes, geom_avail = geom_avail,
       build_gg = build_gg, reduce = reduce, val_spec = val_spec)
}

#' @param state the shared app-state object (see [new_app_state()]).
#' @param dark_mode A reactive returning `TRUE` in dark mode (navbar toggle),
#'   driving the plot theme contrast.
#' @return Invisible NULL (consumes the dds; does not mutate it).
mod_expression_server <- function(id, state, dark_mode = reactive(FALSE)) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    dark <- function() isTRUE(dark_mode())
    eng <- plot_engine_server(input, output, session, state)
    # One "Showing:" selection shared (synced) across both plot pills.
    showing_samples <- plot_subset_server(input, output, session, state,
                                          suffixes = c("expr", "expr_set", "hm"))

    # ---- Single genes: gene search -> one row of the value matrix ----------
    gene_search <- gene_search_server(input, output, session, state, "gene",
                                      multiple = FALSE)
    gene_resolved <- reactive({
      r <- gene_search(); if (length(r$records)) r$records[[1]] else NULL
    })
    # Reduce fns return list(ok, ...) rather than calling validate(): they run
    # inside the deferred spec (so the source is gated behind Render too), and a
    # validate() there would be swallowed. build_gg surfaces `msg` when !ok.
    reduce_single_factory <- function(val_spec) function(mm) {
      rf <- gene_resolved()
      if (is.null(rf) || is.na(rf$id))
        return(list(ok = FALSE, msg = "Enter a feature to plot."))
      if (!rf$id %in% rownames(mm$mat))
        return(list(ok = FALSE, msg = sprintf(
          "'%s' is not in the %s matrix (VST is endogenous-only).",
          rf$match, val_spec()$assay)))
      spec <- val_spec()
      raw <- as.numeric(mm$mat[rf$id, ]); names(raw) <- colnames(mm$mat)
      list(ok = TRUE,
           values = expr_transform(raw, spec$transform, spec$pseudocount),
           y_lab = if (identical(spec$transform, "none")) mm$label else spec$label,
           title = sprintf("%s (%s)", rf$match, rf$id), subtitle = NULL,
           gene_id = rf$id, gene_label = rf$match,
           hover = function(df) sprintf("%s\n%s: %.3g", df$sample, rf$match, df$value))
    }
    single <- .expr_dist_server(input, output, session, state, cfg = list(
      prefix = "", val_suffix = "val", derive_key = "expr_value_mat",
      plot_id = "gene", stale_id = "gene_stale",
      dark = dark, engine = eng, showing_samples = showing_samples,
      reduce_factory = reduce_single_factory,
      source_sig = reactive({           # the gene id gates the plot (with the assay)
        rf <- gene_resolved(); if (is.null(rf)) NA_character_ else rf$id %||% NA_character_
      })))
    # Expose single-gene handles at module scope (testServer reads these).
    gene_out      <- single$out         # deferred: value() = list(mat_cols, red)
    group_sizes   <- single$group_sizes
    geom_avail    <- single$geom_avail
    build_gene_gg <- single$build_gg
    gene_values   <- single$reduce

    # ---- Gene sets > Aggregate expression: a set -> one score per sample ---
    set_search <- gene_search_server(input, output, session, state, "setsearch",
                                     multiple = TRUE,
                                     search_modes = c("exact", "contains", "regex"))
    set_ids <- reactive({
      if (identical(input$set_source %||% "saved", "search")) set_search()$ids
      else {
        gs <- state$gene_sets %||% list(); nm <- input$set_pick
        if (is.null(nm) || !nm %in% names(gs)) character(0)
        else gene_set_ids_for(gs[[nm]])                 # full authored membership
      }
    })
    set_title <- reactive({
      if (identical(input$set_source %||% "saved", "search")) {
        n <- length(set_ids())
        sprintf("Custom search (%d gene%s)", n, if (n == 1L) "" else "s")
      } else input$set_pick %||% "Gene set"
    })
    output$set_pick_ui <- renderUI({
      gs <- state$gene_sets %||% list()
      if (!length(gs))
        return(helpText(class = "small text-muted",
          "No saved gene sets yet - define one on Gene Sets > Manage, or use Quick search."))
      cur <- shiny::isolate(input$set_pick)
      sel <- if (!is.null(cur) && cur %in% names(gs)) cur else names(gs)[1]
      selectInput(ns("set_pick"), "Saved gene set", choices = names(gs), selected = sel)
    })

    reduce_set_factory <- function(val_spec) function(mm) {
      ids <- set_ids()
      if (!length(ids))
        return(list(ok = FALSE, msg = "Pick a saved set or search genes to plot."))
      spec <- val_spec()
      zsc <- isTRUE(input$set_zscore %||% TRUE)
      meth <- input$set_method %||% "mean"
      agg <- expr_set_aggregate(
        mm$mat, ids, counts = DESeq2::counts(state$working),
        method = meth, zscore = zsc,
        only_expressed = isTRUE(input$set_only_expr %||% TRUE),
        transform = spec$transform, pseudocount = spec$pseudocount)
      if (is.null(agg$values) || agg$n_used == 0)
        return(list(ok = FALSE, msg = sprintf(
          "No genes to plot: 0 of %d gene%s survive the present / expression filter.",
          agg$n_total, if (agg$n_total == 1L) "" else "s")))
      meth_lab <- if (identical(agg$method, "median")) "Median" else "Mean"
      base_lab <- if (identical(spec$transform, "none")) mm$label else spec$label
      y_lab <- sprintf("%s %s%s", meth_lab, if (zsc) "z-scored " else "", base_lab)
      pct <- if (agg$n_total) round(100 * agg$n_used / agg$n_total) else 0L
      sub <- sprintf("%s expression of %d of %d genes (%d%%) within the set",
                     meth_lab, agg$n_used, agg$n_total, pct)
      if (zsc && agg$n_nonvar > 0L)
        sub <- paste0(sub, sprintf("; %d non-varying (z-score 0)", agg$n_nonvar))
      list(ok = TRUE, values = agg$values, y_lab = y_lab, title = set_title(),
           subtitle = sub, accounting = agg,
           hover = function(df) sprintf("%s\n%s score: %.3g", df$sample, meth_lab, df$value))
    }
    set <- .expr_dist_server(input, output, session, state, cfg = list(
      prefix = "set_", val_suffix = "set_val", derive_key = "expr_set_value_mat",
      plot_id = "geneset", stale_id = "geneset_stale",
      dark = dark, engine = eng, showing_samples = showing_samples,
      reduce_factory = reduce_set_factory,
      source_sig = reactive(list(          # set + aggregation options gate the plot
        set_ids(), input$set_method %||% "mean",
        isTRUE(input$set_zscore %||% TRUE), isTRUE(input$set_only_expr %||% TRUE)))))
    set_values <- set$reduce               # exposed for tests
    build_set_gg <- set$build_gg
    set_out <- set$out

    # ---- Gene sets > Heatmap (P7c) -----------------------------------------
    hm_search <- gene_search_server(input, output, session, state, "hmsearch",
      multiple = TRUE, search_modes = c("exact", "contains", "regex"))
    hm_rowsel <- gene_search_server(input, output, session, state, "hmrowsel",
      multiple = TRUE, search_modes = c("exact", "contains", "regex"))
    hm_val <- expr_value_server(input, output, session, state, "hm_val",
      include_vst = TRUE, include_norm = TRUE, default_fn = expr_default_assay)
    hm_ramp <- continuous_palette_server(input, output, session, "hm_ramp",
                                         .hm_ramp_default(TRUE))

    hm_coldata <- function() as.data.frame(SummarizedExperiment::colData(state$working))
    hm_rowdata <- function() as.data.frame(SummarizedExperiment::rowData(state$working))
    hm_disc <- function(df) names(df)[vapply(df, function(x)
      is.factor(x) || is.character(x) || is.logical(x), logical(1))]

    hm_ids <- reactive({
      if (identical(input$hm_source %||% "saved", "search")) hm_search()$ids
      else {
        gs <- state$gene_sets %||% list(); nm <- input$hm_pick
        if (is.null(nm) || !nm %in% names(gs)) character(0)
        else gene_set_ids_for(gs[[nm]])                  # full authored membership
      }
    })
    output$hm_pick_ui <- renderUI({
      gs <- state$gene_sets %||% list()
      if (!length(gs))
        return(helpText(class = "small text-muted",
          "No saved gene sets yet - define one on Gene Sets > Manage, or use Quick search."))
      cur <- shiny::isolate(input$hm_pick)
      sel <- if (!is.null(cur) && cur %in% names(gs)) cur else names(gs)[1]
      selectInput(ns("hm_pick"), "Saved gene set", choices = names(gs), selected = sel)
    })

    # Label-source selectors (default row = <feature_type>_name -> rownames; col =
    # colnames). Forced to init while hidden so their defaults hold when collapsed.
    output$hm_row_source_ui <- renderUI({
      req(state$working)
      cols <- hm_disc(hm_rowdata())
      guess <- paste0(state$meta$feature_type %||% "gene", "_name")
      def <- if (guess %in% cols) guess else "__id__"
      selectInput(ns("hm_row_source"), "Label by",
                  c("Feature ID (rownames)" = "__id__", stats::setNames(cols, cols)),
                  selected = def)
    })
    output$hm_col_source_ui <- renderUI({
      req(state$working)
      cols <- hm_disc(hm_coldata())
      selectInput(ns("hm_col_source"), "Label by",
                  c("Sample ID (colnames)" = "__id__", stats::setNames(cols, cols)),
                  selected = "__id__")
    })
    output$hm_col_sel_ui <- renderUI({
      req(state$working)
      selectizeInput(ns("hm_col_sel"), NULL, choices = colnames(state$working),
                     multiple = TRUE, options = list(placeholder = "Samples to label"))
    })
    output$hm_anno_ui <- renderUI({
      req(state$working)
      pv <- primary_design_var(state$working)
      selectizeInput(ns("hm_anno"), "Annotate by (one or more)",
                     choices = aes_choices(aes_catalog(state), none = FALSE),
                     selected = if (!is.na(pv)) pv else character(0), multiple = TRUE)
    })
    outputOptions(output, "hm_row_source_ui", suspendWhenHidden = FALSE)
    outputOptions(output, "hm_col_source_ui", suspendWhenHidden = FALSE)
    outputOptions(output, "hm_anno_ui", suspendWhenHidden = FALSE)
    observeEvent(input$hm_clear_anno,
                 updateSelectizeInput(session, "hm_anno", selected = character(0)))

    # The gated snapshot: the value matrix (cached in `derived`), the prepared
    # gene-set matrix, the plotted (column-subset) matrix, resolved labels/marks,
    # cluster/dend flags, ramp config, and annotation *values*. Colours resolve
    # live at draw. Everything here waits for Render (no cheap live layer).
    hm_spec <- reactive({
      req(state$working)
      dds <- state$working
      spec <- hm_val()
      vm <- state_derive(state, "expr_hm_value_mat", params = list(spec$assay),
                         expr = function() expr_value_matrix(dds, spec$assay))
      zsc <- isTRUE(input$hm_zscore %||% TRUE)
      hm <- expr_heatmap_matrix(vm$mat, hm_ids(), counts = DESeq2::counts(dds),
              zscore = zsc, only_expressed = isTRUE(input$hm_only_expr %||% TRUE),
              transform = spec$transform, pseudocount = spec$pseudocount)
      base_lab <- if (identical(spec$transform, "none")) vm$label else spec$label

      out <- list(hm = hm, mat = NULL, zscored = zsc, assay = spec$assay,
                  legend = if (zsc) "z-score" else base_lab,
                  row_mode = "none", col_mode = "none",
                  row_labels = character(0), col_labels = character(0),
                  row_mark_at = integer(0), row_mark_lab = character(0),
                  col_mark_at = integer(0), col_mark_lab = character(0),
                  row_cov = NULL, col_cov = NULL, anno = NULL,
                  cluster_rows = FALSE, cluster_cols = FALSE,
                  show_row_dend = FALSE, show_col_dend = FALSE,
                  empty_msg = "Pick a saved set or search genes, then click Render.")
      if (is.null(hm$mat)) return(out)
      show <- intersect(colnames(hm$mat), showing_samples())
      if (!length(show)) { out$empty_msg <- "No samples in the current 'Showing' selection."; return(out) }
      pm <- hm$mat[, show, drop = FALSE]
      out$mat <- pm
      n_row <- nrow(pm); n_col <- ncol(pm)

      rd <- hm_rowdata(); cd <- hm_coldata()
      out$row_labels <- expr_heatmap_labels(rownames(pm), input$hm_row_source %||% "__id__",
                                            meta = rd, meta_keys = rownames(dds))
      out$col_labels <- expr_heatmap_labels(colnames(pm), input$hm_col_source %||% "__id__",
                                            meta = cd, meta_keys = colnames(dds))
      rmode <- input$hm_row_mode %||% "auto"
      if (identical(rmode, "auto")) rmode <- heatmap_label_default(n_row, .hm_row_label_max())
      cmode <- input$hm_col_mode %||% "auto"
      if (identical(cmode, "auto")) cmode <- heatmap_label_default(n_col, .hm_col_label_max())
      out$row_mode <- rmode; out$col_mode <- cmode
      if (identical(rmode, "selected")) {
        sel <- hm_rowsel()$ids
        out$row_cov <- expr_label_coverage(sel, rownames(pm))
        at <- which(rownames(pm) %in% sel)
        out$row_mark_at <- at; out$row_mark_lab <- out$row_labels[at]
      }
      if (identical(cmode, "selected")) {
        sel <- input$hm_col_sel %||% character(0)
        out$col_cov <- expr_label_coverage(sel, colnames(pm))
        at <- which(colnames(pm) %in% sel)
        out$col_mark_at <- at; out$col_mark_lab <- out$col_labels[at]
      }
      cl_r <- isTRUE(input$hm_cluster_rows %||% TRUE)
      cl_c <- isTRUE(input$hm_cluster_cols %||% TRUE)
      out$cluster_rows <- cl_r && n_row > 1L
      out$cluster_cols <- cl_c && n_col > 1L
      dend <- function(mode, cl, n)
        cl && switch(mode %||% "auto", show = TRUE, hide = FALSE, n <= .hm_dend_max())
      out$show_row_dend <- dend(input$hm_row_dend, out$cluster_rows, n_row)
      out$show_col_dend <- dend(input$hm_col_dend, out$cluster_cols, n_col)

      valid <- vapply(aes_catalog(state), `[[`, "", "key")
      sel_anno <- intersect(input$hm_anno, valid)
      if (length(sel_anno)) {
        built <- Filter(Negate(is.null), lapply(sel_anno, function(k) {
          r <- aes_resolve(state, k, show)
          if (is.null(r)) NULL else list(key = k, name = r$label, values = r$values)
        }))
        if (length(built)) out$anno <- list(
          keys = vapply(built, `[[`, "", "key"),
          df = data.frame(stats::setNames(lapply(built, `[[`, "values"),
                                          vapply(built, `[[`, "", "name")),
                          row.names = show, check.names = FALSE, stringsAsFactors = FALSE))
      }
      out
    })
    # Render-only (no auto id created -> auto never fires); stale banner still works.
    hm_shown <- eng$deferred("hm_no_auto", "hm_render", hm_spec,
      sig = reactive(list(hm_val(), hm_ids(), isTRUE(input$hm_zscore %||% TRUE),
                          isTRUE(input$hm_only_expr %||% TRUE), showing_samples(),
                          input$hm_row_source, input$hm_col_source,
                          input$hm_row_mode, input$hm_col_mode,
                          hm_rowsel()$ids, input$hm_col_sel,
                          isTRUE(input$hm_cluster_rows %||% TRUE),
                          isTRUE(input$hm_cluster_cols %||% TRUE),
                          input$hm_row_dend, input$hm_col_dend, input$hm_anno,
                          state$data_version)))         # ramp is a live aesthetic (not gated)
    output$hm_stale <- eng$stale_note(hm_shown)

    output$hm_info <- renderUI({
      s <- hm_shown$value(); if (is.null(s)) return(NULL)
      hm <- s$hm
      # "not in the value matrix" (not "absent from the dataset"): a spike-in /
      # exogenous gene is in the dds but excluded from an endogenous-only VST
      # matrix, so the same set plots fewer genes under VST than under a stored assay.
      m <- sprintf("%d of %d genes shown%s%s.", hm$n_used, hm$n_total,
                   if (hm$n_absent > 0) sprintf("; %d not in the value matrix", hm$n_absent) else "",
                   if (s$zscored && hm$n_nonvar > 0)
                     sprintf("; %d non-varying (flat row)", hm$n_nonvar) else "")
      lines <- list(tags$div(m))
      if (!is.null(s$row_cov) && s$row_cov$n_hidden > 0)
        lines <- c(lines, list(tags$div(sprintf(
          "Row labels: %d of %d selected not shown (not in the set / current view).",
          s$row_cov$n_hidden, s$row_cov$n_selected))))
      if (!is.null(s$col_cov) && s$col_cov$n_hidden > 0)
        lines <- c(lines, list(tags$div(sprintf(
          "Column labels: %d of %d selected not shown (outside the current 'Showing').",
          s$col_cov$n_hidden, s$col_cov$n_selected))))
      tags$div(class = "small text-muted mb-2", lines)
    })

    output$hm_plot <- renderPlot({
      validate(need(!is.null(hm_shown$value()), "Click Render to draw the heatmap."))
      validate(need(requireNamespace("ComplexHeatmap", quietly = TRUE),
                    "Install 'ComplexHeatmap' to show the gene-set heatmap."))
      s <- hm_shown$value()
      validate(need(!is.null(s$mat) && nrow(s$mat) > 0,
                    s$empty_msg %||% "No genes to plot."))
      validate(need(ncol(s$mat) > 0, "No samples in the current 'Showing' selection."))
      # Colour ramp resolved LIVE (a display aesthetic, not gated): a Custom ramp,
      # or the assay's Palette-page config for a raw (non-z-scored) matrix. So a
      # ramp / Palette edit recolours without a re-render (mirrors the QC heatmap).
      ramp <- if (identical(input$hm_ramp_src %||% "custom", "palette") && !s$zscored)
                (state$palette$assays[[s$assay]] %||% hm_ramp()) else hm_ramp()
      col_fun <- .hm_col_fun(s$mat, ramp, s$zscored)
      validate(need(!is.null(col_fun), "Install 'circlize' for the heatmap colour scale."))
      anno <- NULL
      if (!is.null(s$anno)) {
        show <- colnames(s$mat)
        res_list <- Filter(Negate(is.null),
                           lapply(s$anno$keys, function(k) aes_resolve(state, k, show)))
        if (length(res_list)) {
          anno_col <- stats::setNames(lapply(res_list, aes_heatmap_col),
                                      vapply(res_list, function(r) r$label, ""))
          anno <- ComplexHeatmap::HeatmapAnnotation(
            df = s$anno$df[show, , drop = FALSE], col = anno_col)
        }
      }
      ComplexHeatmap::draw(.hm_build(s, anno, col_fun))
    })
    # Exposed for tests.
    hm_out <- hm_shown

    invisible(NULL)
  })
}
