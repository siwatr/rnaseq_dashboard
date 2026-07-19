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
                              legend_pos = "right", dark_theme = FALSE,
                              interactive = FALSE) {
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
        sliderInput(nid("dot_cex"), "Point spacing (cex)", 0.2, 3, 1, 0.1))))
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
  heatmap_pill <- bslib::nav_panel(
    "Heatmap",
    bslib::card(
      bslib::card_body(
        tags$p(class = "text-muted",
               "The gene-set expression heatmap is coming in the next update (P7c).")
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
    dv <- tryCatch(all.vars(DESeq2::design(state$working)), error = function(e) character(0))
    if (length(dv) && dv[1] %in% disc) return(dv[1])
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
      dark_theme = dark(), interactive = interactive)
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
                                          suffixes = c("expr", "expr_set"))

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

    invisible(NULL)
  })
}
