# Page 7: Expression. A gene-expression browsing surface (renamed from the old
# "Heatmap" stub). A navset_card_tab with two tabs:
#   Single genes - one feature at a time as a layered violin -> box -> dots overlay,
#                  grouped by a colData variable, reusing the shared plot machinery
#                  (gene search, expr-value control, plot engine + deferred render,
#                  "Showing" subset, aes_helpers colour).
#   Gene sets    - a ComplexHeatmap over a named gene set (built in P7b).
#
# Single-gene value extraction (esp. VST) is cached in `derived` and gated behind
# the deferred Render button; geom toggles + colour are cheap live re-plots.

# Colour/fill scale + layered plot for one gene. `df` carries sample/group/value
# (+ optional `colour`). Distribution geoms fill by a discrete colour attribute;
# dots colour by it. ggbeeswarm draws the dots (geom_jitter fallback).
#' @importFrom ggplot2 .data
.expr_single_plot <- function(df, x_lab, y_lab, title, subtitle = NULL,
                              show_violin = TRUE, show_box = TRUE, show_dots = TRUE,
                              disc_colour = FALSE, colour_lab = NULL,
                              fill_scale = NULL, colour_scale = NULL,
                              legend_pos = "right", dark_theme = FALSE,
                              interactive = FALSE) {
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$group, y = .data$value))
  fill_map   <- if (disc_colour) ggplot2::aes(fill = .data$colour) else NULL
  colour_map <- if (!is.null(df$colour)) ggplot2::aes(colour = .data$colour) else NULL

  if (isTRUE(show_violin))
    p <- p + ggplot2::geom_violin(mapping = fill_map, alpha = 0.35,
                                  colour = NA, scale = "width",
                                  show.legend = disc_colour)
  if (isTRUE(show_box))
    p <- p + ggplot2::geom_boxplot(mapping = fill_map, width = 0.18,
                                   alpha = 0.7, outlier.shape = NA,
                                   show.legend = FALSE)
  if (isTRUE(show_dots)) {
    dot_aes <- colour_map
    if (interactive) {
      dot_aes <- if (is.null(dot_aes)) ggplot2::aes(text = .data$text)
                 else utils::modifyList(dot_aes, ggplot2::aes(text = .data$text))
    }
    dot_args <- list(mapping = dot_aes, size = 1.9, alpha = 0.9)
    p <- p + if (requireNamespace("ggbeeswarm", quietly = TRUE)) {
      do.call(ggbeeswarm::geom_quasirandom, dot_args)
    } else {
      do.call(ggplot2::geom_jitter,
              c(dot_args, list(width = 0.18, height = 0)))
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

mod_expression_ui <- function(id) {
  ns <- NS(id)
  bslib::navset_card_tab(
    id = ns("tabs"),
    title = tags$h2("Expression", class = "fs-6 mb-0"),

    # --- Single genes -------------------------------------------------------
    bslib::nav_panel(
      "Single genes",
      bslib::layout_sidebar(
        sidebar = bslib::sidebar(
          title = tags$h3("Single gene", class = "fs-6 mb-0"), width = 340,
          bslib::accordion(
            open = c("Gene", "Grouping & colour"),
            bslib::accordion_panel(
              "Gene", icon = icon("magnifying-glass"),
              gene_search_ui(ns, "gene", multiple = FALSE, dup_toggle = TRUE,
                             placeholder = "e.g. Gene1"),
              expr_value_ui(ns, "val")
            ),
            bslib::accordion_panel(
              "Grouping & colour", icon = icon("layer-group"),
              uiOutput(ns("x_group_ui")),
              uiOutput(ns("colour_ui"))
            ),
            bslib::accordion_panel(
              "Plot elements", icon = icon("chart-simple"),
              uiOutput(ns("geom_ui"))
            ),
            bslib::accordion_panel(
              "Render", icon = icon("play"),
              uiOutput(ns("auto_ui")),
              actionButton(ns("render"), "Render", icon = icon("play"),
                           class = "btn-primary")
            )
          ),
          plot_subset_ui(ns, "expr")
        ),
        bslib::card(
          bslib::card_header(tags$h3("Single-gene expression", class = "fs-6 mb-0")),
          uiOutput(ns("gene_stale")),
          shinycssloaders::withSpinner(.plot_dual(ns("gene_container")),
                                       proxy.height = "400px")
        )
      )
    ),

    # --- Gene sets (P7b) ----------------------------------------------------
    bslib::nav_panel(
      "Gene sets",
      bslib::card(
        bslib::card_body(
          tags$p(class = "text-muted",
                 "The gene-set expression heatmap is coming in the next update (P7b).")
        )
      )
    )
  )
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
    dual_plot <- eng$dual_plot; deferred <- eng$deferred; stale_note <- eng$stale_note
    showing_samples <- plot_subset_server(input, output, session, state, suffixes = "expr")

    coldata <- function() as.data.frame(SummarizedExperiment::colData(state$working))
    # Columns usable as a grouping x-axis: discrete (factor / character / logical).
    discrete_cols <- function() {
      cd <- coldata()
      if (!ncol(cd)) return(character(0))
      names(cd)[vapply(cd, function(x) is.factor(x) || is.character(x) || is.logical(x),
                       logical(1))]
    }
    # Default grouping: first design variable if discrete, else first discrete col.
    default_x <- function() {
      disc <- discrete_cols()
      dv <- tryCatch(all.vars(DESeq2::design(state$working)), error = function(e) character(0))
      if (length(dv) && dv[1] %in% disc) return(dv[1])
      if (length(disc)) disc[1] else "__none__"
    }

    # --- Sidebar controls ---------------------------------------------------
    output$x_group_ui <- renderUI({
      req(state$working)
      selectInput(ns("x_group"), "Group by (x-axis)",
                  choices = group_field_choices(discrete_cols(), none = FALSE),
                  selected = default_x())
    })
    output$colour_ui <- renderUI({
      req(state$working)
      selectInput(ns("colour_by"), "Colour by",
                  choices = aes_choices(aes_catalog(state), none = TRUE),
                  selected = default_x())
    })
    output$auto_ui <- renderUI({
      req(state$working)
      checkboxInput(ns("auto"), "Auto-render", value = TRUE)
    })

    # The shared gene-search + expr-value controls.
    gene_search <- gene_search_server(input, output, session, state, "gene",
                                      multiple = FALSE)
    gene_resolved <- reactive({
      r <- gene_search(); if (length(r$records)) r$records[[1]] else NULL
    })
    val_spec <- expr_value_server(input, output, session, state, "val",
                                  include_vst = TRUE, default_fn = expr_default_assay)

    # Per-group sample counts over the shown samples -> which geoms are offered.
    group_sizes <- reactive({
      req(state$working, input$x_group)
      cd <- coldata()
      if (!input$x_group %in% names(cd)) return(integer(0))
      shown <- intersect(rownames(cd), showing_samples())
      g <- cd[shown, input$x_group]
      tab <- table(g[!is.na(g)])
      as.integer(tab[tab > 0])            # drop empty groups (e.g. filtered out by "Showing")
    })
    geom_avail <- reactive(expr_geom_availability(group_sizes()))

    # Geom toggles: dots offered always (disabled when the group is too large);
    # violin/box only when at least one group is large enough to summarize.
    output$geom_ui <- renderUI({
      av <- geom_avail()
      dist <- isTRUE(av$dist_shown)
      tagList(
        if (dist) bslib::input_switch(ns("show_violin"), "Violin", value = TRUE)
        else helpText(class = "small text-muted",
                      "Violin / box hidden: groups are too small to summarize."),
        if (dist) bslib::input_switch(ns("show_box"), "Boxplot", value = TRUE),
        if (isTRUE(av$dots_allowed))
          bslib::input_switch(ns("show_dots"), "Data points",
                              value = isTRUE(av$dots_default))
        else helpText(class = "small text-muted",
                      sprintf("Data points hidden: a group has %d+ samples (overplotting).",
                              av$n_max))
      )
    })
    # Resolve a toggle whose UI may be absent (dist hidden -> violin/box off;
    # dots not allowed -> off), falling back to the availability default.
    on_flag <- function(id, allowed, default) {
      if (!isTRUE(allowed)) return(FALSE)
      v <- input[[id]]; if (is.null(v)) default else isTRUE(v)
    }

    # --- Deferred value matrix (caches the assay matrix incl. VST) ----------
    # Only the (expensive) value matrix is gated/cached behind Render, keyed on the
    # assay. The per-gene lookup + transform + grouping + colour are cheap live
    # re-plots on top (so changing the gene doesn't need a Render click), and any
    # gene-not-in-matrix message surfaces in the plot builder (a deferred spec would
    # swallow it).
    mat_spec <- reactive({
      req(state$working)
      assay <- val_spec()$assay
      state_derive(state, "expr_value_mat", params = list(assay),
                   expr = function() expr_value_matrix(state$working, assay))
    })
    mat_shown <- deferred("auto", "render", mat_spec,
      sig = reactive(list(val_spec()$assay, state$data_version)))
    output$gene_stale <- stale_note(mat_shown)

    # Colour aesthetic (discrete -> fill + colour; continuous -> dot colour only).
    colour_resolve <- function(samples) {
      sel <- input$colour_by %||% "__none__"
      if (identical(sel, "__none__")) return(NULL)
      res <- aes_resolve(state, sel, samples)
      if (is.null(res)) return(NULL)
      disc <- identical(res$kind, "discrete")
      list(values = res$values, discrete = disc, lab = res$label,
           fill_scale = if (disc) aes_ggplot_scale(res, "fill") else NULL,
           colour_scale = aes_ggplot_scale(res, "colour"))
    }

    # Extract one gene's per-sample values from the cached matrix + spec transform.
    gene_values <- function(mm) {
      rf <- gene_resolved()
      validate(need(!is.null(rf) && !is.na(rf$id), "Enter a feature to plot."))
      validate(need(rf$id %in% rownames(mm$mat),
                    sprintf("'%s' is not in the %s matrix (VST is endogenous-only).",
                            rf$match, val_spec()$assay)))
      spec <- val_spec()
      raw <- as.numeric(mm$mat[rf$id, ]); names(raw) <- colnames(mm$mat)
      list(values = expr_transform(raw, spec$transform, spec$pseudocount),
           y_lab = if (identical(spec$transform, "none")) mm$label else spec$label,
           gene_id = rf$id, gene_label = rf$match)
    }

    # --- The plot (live re-plot from the cached value matrix) ---------------
    build_gene_gg <- function(interactive) {
      mm <- mat_shown$value()
      validate(need(!is.null(mm),
                    "Click Render (or enable auto-render) to plot the gene."))
      req(input$x_group)
      v <- gene_values(mm)
      cd <- coldata()
      validate(need(input$x_group %in% names(cd), "Pick a grouping variable."))
      shown <- intersect(names(v$values), showing_samples())
      shown <- intersect(shown, rownames(cd))
      validate(need(length(shown) > 0, "No samples in the current 'Showing' selection."))
      groups <- cd[shown, input$x_group]
      cr <- colour_resolve(shown)
      df <- expr_long_frame(v$values[shown], groups, shown,
                            colour = if (!is.null(cr)) cr$values else NULL)
      validate(need(nrow(df) > 0, "No samples with a non-missing group."))

      av <- geom_avail()
      sv <- on_flag("show_violin", av$dist_shown, TRUE)
      sb <- on_flag("show_box",    av$dist_shown, TRUE)
      sd <- on_flag("show_dots",   av$dots_allowed, av$dots_default)
      validate(need(sv || sb || sd, "Enable at least one plot element."))
      disc <- !is.null(cr) && isTRUE(cr$discrete)
      if (interactive)
        df$text <- sprintf("%s\n%s: %.3g", df$sample, v$gene_label, df$value)

      .expr_single_plot(
        df, x_lab = input$x_group, y_lab = v$y_lab,
        title = sprintf("%s (%s)", v$gene_label, v$gene_id),
        show_violin = sv, show_box = sb, show_dots = sd,
        disc_colour = disc, colour_lab = if (!is.null(cr)) cr$lab else NULL,
        fill_scale = if (!is.null(cr)) cr$fill_scale else NULL,
        colour_scale = if (!is.null(cr)) cr$colour_scale else NULL,
        dark_theme = dark(), interactive = interactive)
    }
    dual_plot("gene", build_gene_gg, n_elements = reactive({
      mm <- mat_shown$value(); if (is.null(mm)) 0L
      else length(intersect(colnames(mm$mat), showing_samples()))
    }), height = "440px")

    invisible(NULL)
  })
}
