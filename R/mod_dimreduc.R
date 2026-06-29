# Page 4: Dimensionality reduction (PCA). Single-panel PCA scatter + a scree
# (%-variance) plot. The embedding is computed once and cached (state_derive),
# behind the shared deferred render gate; aesthetics (PC axes, colour, shape) are
# cheap re-plots from the cached scores. Colour/shape consume the project Palette
# configs. Reuses the shared plot engine (R/mod_plot_engine.R) and the "Showing"
# subset control (R/mod_plot_subset.R). t-SNE/UMAP are deferred (P4c).

# Recommended PCA inputs (computed) + whatever stored assays exist, de-duplicated.
.pca_assay_choices <- function(dds) {
  stored <- SummarizedExperiment::assayNames(dds)
  base <- c("VST" = "vst", "Normalized log-counts" = "norm_logcounts")
  extra <- setdiff(stored, "vst")
  c(base, stats::setNames(extra, extra))
}

# Build the PCA scatter. `df` carries x/y (the two PCs) + optional `col`/`shp`
# (already factored/numeric) + `text` (hover). `colour_scale`/`shape_scale` are
# pre-built ggplot scales (or NULL). `interactive` adds the plotly hover aes only.
.pca_scatter_plot <- function(df, xlab, ylab, subtitle, colour_lab = NULL,
                              shape_lab = NULL, colour_scale = NULL,
                              shape_scale = NULL, interactive = FALSE) {
  aes_args <- list(x = quote(x), y = quote(y))
  if ("col" %in% names(df)) aes_args$colour <- quote(col)
  if ("shp" %in% names(df)) aes_args$shape  <- quote(shp)
  if (interactive)          aes_args$text   <- quote(text)
  p <- ggplot2::ggplot(df, do.call(ggplot2::aes, aes_args)) +
    ggplot2::geom_point(size = 3, alpha = 0.9) +
    ggplot2::labs(x = xlab, y = ylab, subtitle = subtitle,
                  colour = colour_lab, shape = shape_lab) +
    ggplot2::theme_minimal(base_size = 13)
  if (!is.null(colour_scale)) p <- p + colour_scale
  if (!is.null(shape_scale))  p <- p + shape_scale
  p
}

# Scree bar plot of %variance per PC (top `n_show`); the two plotted PCs highlighted.
.pca_scree_plot <- function(var_pct, highlight = integer(0), n_show = 10L,
                            interactive = FALSE) {
  n <- min(n_show, length(var_pct))
  d <- data.frame(pc = factor(paste0("PC", seq_len(n)), levels = paste0("PC", seq_len(n))),
                  pct = var_pct[seq_len(n)],
                  hi = seq_len(n) %in% highlight)
  d$text <- sprintf("%s: %.1f%%", d$pc, d$pct)
  aes_args <- list(x = quote(pc), y = quote(pct), fill = quote(hi))
  if (interactive) aes_args$text <- quote(text)
  ggplot2::ggplot(d, do.call(ggplot2::aes, aes_args)) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_manual(values = c(`FALSE` = "grey70", `TRUE` = "#1f77b4"),
                               guide = "none") +
    ggplot2::labs(x = NULL, y = "% variance", title = "Variance explained") +
    ggplot2::theme_minimal(base_size = 13)
}

mod_dimreduc_ui <- function(id) {
  ns <- NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      title = tags$h3("Dim. reduction", class = "fs-6 mb-0"), width = 320,
      helpText(class = "small",
               "PCA of the top-variable genes (endogenous only). Pick an input, then Render."),
      uiOutput(ns("assay_ui")),
      uiOutput(ns("assay_note")),
      uiOutput(ns("log_ui")),
      numericInput(ns("n_top"), "Top variable genes", value = 500, min = 2, step = 50),
      uiOutput(ns("pc_ui")),
      uiOutput(ns("colour_ui")),
      uiOutput(ns("gene_ui")),
      uiOutput(ns("shape_ui")),
      plot_subset_ui(ns, "pca"),
      uiOutput(ns("auto_ui")),
      actionButton(ns("render"), "Render", icon = icon("play"), class = "btn-primary")
    ),
    bslib::card(
      bslib::card_header(tags$h3("PCA", class = "fs-6 mb-0")),
      uiOutput(ns("pca_stale")),
      .plot_dual(ns("pca_container"))
    ),
    bslib::card(
      bslib::card_header(tags$h3("Scree", class = "fs-6 mb-0")),
      .plot_dual(ns("scree_container"))
    )
  )
}

#' @param state the shared app-state object (see [new_app_state()]).
#' @return Invisible NULL (consumes the dds; does not mutate it).
mod_dimreduc_server <- function(id, state) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    eng <- plot_engine_server(input, output, session, state)
    dual_plot <- eng$dual_plot; deferred <- eng$deferred; stale_note <- eng$stale_note
    showing_samples <- plot_subset_server(input, output, session, state, suffixes = "pca")

    feature_type <- function() state$meta$feature_type %||% "gene"
    coldata <- function() as.data.frame(SummarizedExperiment::colData(state$working))
    discrete_cols <- function() {
      cd <- coldata(); names(cd)[!vapply(cd, is.numeric, logical(1))]
    }

    # --- Sidebar controls (data-dependent -> server-rendered) ----------------
    output$assay_ui <- renderUI({
      req(state$working)
      selectInput(ns("assay"), "Input (assay / transform)",
                  choices = .pca_assay_choices(state$working), selected = "vst")
    })
    output$assay_note <- renderUI({
      req(input$assay)
      ad <- pca_assay_advice(input$assay)
      if (is.null(ad$msg)) return(NULL)
      cls <- if (identical(ad$tier, "unsuitable")) "alert-danger" else "alert-warning"
      tags$div(class = paste("alert py-1 px-2 small", cls), ad$msg)
    })
    # Log toggle only for non-stabilized stored assays; default per advice.
    output$log_ui <- renderUI({
      req(input$assay)
      if (input$assay %in% c("vst", "norm_logcounts", "logcounts")) return(NULL)
      checkboxInput(ns("log_transform"), "Apply log2(x + 1) before PCA",
                    value = pca_assay_advice(input$assay)$recommend_log)
    })
    output$auto_ui <- renderUI({
      req(state$working)
      checkboxInput(ns("auto"), "Auto-render", value = ncol(state$working) <= 150L)
    })
    output$colour_ui <- renderUI({
      req(state$working)
      ch <- c("(none)" = "__none__", stats::setNames(colnames(coldata()), colnames(coldata())),
              "Gene expression" = "__gene__")
      selectInput(ns("colour_by"), "Colour by", choices = ch, selected = "__none__")
    })
    output$gene_ui <- renderUI({
      if (!identical(input$colour_by, "__gene__")) return(NULL)
      textInput(ns("gene"), sprintf("%s (name or id)", feature_type()),
                placeholder = "e.g. Gene1")
    })
    output$shape_ui <- renderUI({
      req(state$working)
      ch <- c("(none)" = "__none__", stats::setNames(discrete_cols(), discrete_cols()))
      selectInput(ns("shape_by"), "Shape by (discrete)", choices = ch, selected = "__none__")
    })
    # PC-axis selectors, populated from the computed embedding's PC count.
    output$pc_ui <- renderUI({
      np <- pca_n_pc()
      pcs <- paste0("PC", seq_len(np))
      bslib::layout_columns(
        col_widths = c(6, 6),
        selectInput(ns("pc_x"), "X axis", choices = pcs, selected = "PC1"),
        selectInput(ns("pc_y"), "Y axis", choices = pcs,
                    selected = if (np >= 2L) "PC2" else "PC1"))
    })

    # --- Embedding (cached + deferred) --------------------------------------
    log_on <- function() isTRUE(input$log_transform) &&
      !(input$assay %in% c("vst", "norm_logcounts", "logcounts"))
    pca_spec <- reactive({
      req(state$working, input$assay)
      validate(need(ncol(state$working) >= 3L,
                    "PCA needs at least 3 samples to show PC1 vs PC2."))
      n_top <- max(2L, as.integer(input$n_top %||% 500L))
      assay <- input$assay; lg <- log_on()
      res <- state_derive(state, "pca", params = list(assay, n_top, lg), expr = function() {
        dat <- pca_input(state$working, assay, log_transform = lg)
        pc <- compute_pca(dat$mat, n_top = n_top)
        c(pc, list(label = dat$label))
      })
      res
    })
    pca_shown <- deferred("auto", "render", pca_spec,
      sig = reactive(list(input$assay, input$n_top, log_on(), state$data_version)))
    output$pca_stale <- stale_note(pca_shown)
    # PC count for the axis selectors (0 until first render).
    pca_n_pc <- reactive({ v <- pca_shown$value(); if (is.null(v)) 2L else v$n_pc })

    # --- Colour / shape resolution (consumes the Palette configs) -----------
    # Returns list(values, discrete, lab, scale) for the colour aesthetic, or NULL.
    colour_resolve <- function(samples) {
      sel <- input$colour_by %||% "__none__"
      if (identical(sel, "__none__")) return(NULL)
      dds <- state$working
      if (identical(sel, "__gene__")) {
        gq <- input$gene %||% ""
        validate(need(nzchar(gq), "Enter a feature to colour by."))
        id <- lookup_feature(gq, SummarizedExperiment::rowData(dds), feature_type = feature_type())
        validate(need(!is.na(id), sprintf("Feature '%s' not found.", gq)))
        m <- .qc_assay_matrix(dds, "logcounts")
        vals <- as.numeric(m[id, samples])
        g <- .continuous_scale_from(state$palette$assays[["logcounts"]], vals)
        return(list(values = vals, discrete = FALSE, lab = paste0(gq, "\n(logcounts)"),
                    scale = g))
      }
      cd <- coldata()[samples, , drop = FALSE]
      x <- cd[[sel]]
      if (is.numeric(x)) {
        g <- .continuous_scale_from(state$palette$colData[[sel]], x)
        return(list(values = x, discrete = FALSE, lab = sel, scale = g))
      }
      lv <- as.character(x); lv[is.na(lv)] <- "NA"
      f <- factor(lv)
      cfg <- state$palette$colData[[sel]]
      cols <- palette_discrete(levels(f), cfg$colors, cfg$name %||% "Okabe-Ito", cfg$custom)
      list(values = f, discrete = TRUE, lab = sel,
           scale = ggplot2::scale_colour_manual(values = cols, na.value = "grey70"))
    }
    # Continuous colour scale from a Palette continuous config (or a viridis default).
    .continuous_scale_from <- function(cfg, values) {
      cfg <- cfg %||% list(name = "viridis: viridis")
      g <- palette_gradientn(cfg$name %||% "viridis: viridis", values = values,
                             min = cfg$min, max = cfg$max, custom = cfg$custom,
                             reverse = isTRUE(cfg$reverse))
      ggplot2::scale_colour_gradientn(colours = g$colours, values = g$values, limits = g$limits)
    }
    shape_resolve <- function(samples) {
      sel <- input$shape_by %||% "__none__"
      if (identical(sel, "__none__")) return(NULL)
      x <- coldata()[samples, sel]
      lv <- as.character(x); lv[is.na(lv)] <- "NA"
      f <- factor(lv)
      if (nlevels(f) > 6L)
        return(list(skip = TRUE, lab = sel))   # ggplot's shape palette caps at 6
      list(values = f, lab = sel, scale = NULL)
    }

    # --- The two plots (re-plot live from the cached embedding) -------------
    build_pca_gg <- function(interactive) {
      validate(need(!is.null(pca_shown$value()),
                    "Click Render (or enable auto-render) to compute the PCA."))
      v <- pca_shown$value()
      shown <- intersect(rownames(v$scores), showing_samples())
      validate(need(length(shown) > 0, "No samples in the current 'Showing' selection."))
      px <- input$pc_x %||% "PC1"; py <- input$pc_y %||% "PC2"
      validate(need(px %in% colnames(v$scores) && py %in% colnames(v$scores),
                    "Selected PCs are unavailable."))
      df <- data.frame(x = v$scores[shown, px], y = v$scores[shown, py],
                       sample = shown, stringsAsFactors = FALSE)
      cr <- colour_resolve(shown); sr <- shape_resolve(shown)
      colour_lab <- NULL; shape_lab <- NULL; colour_scale <- NULL; shape_scale <- NULL
      if (!is.null(cr)) { df$col <- cr$values; colour_lab <- cr$lab; colour_scale <- cr$scale }
      if (!is.null(sr) && !isTRUE(sr$skip)) { df$shp <- sr$values; shape_lab <- sr$lab }
      pcpct <- function(p) sprintf("%s: %.0f%%", p, v$var_pct[as.integer(sub("PC", "", p))])
      hov <- if (interactive) {
        df$text <- sprintf("%s\n%s, %s", df$sample, pcpct(px), pcpct(py))
        if (!is.null(cr)) df$text <- paste0(df$text, "\n", cr$lab, ": ", df$col)
      }
      .pca_scatter_plot(df, xlab = pcpct(px), ylab = pcpct(py),
                        subtitle = paste0("Input: ", v$label,
                                          " · ", v$n_genes, " top-variable genes"),
                        colour_lab = colour_lab, shape_lab = shape_lab,
                        colour_scale = colour_scale, shape_scale = shape_scale,
                        interactive = interactive)
    }
    dual_plot("pca", build_pca_gg,
              n_elements = reactive({ v <- pca_shown$value(); if (is.null(v)) 0L else nrow(v$scores) }))

    build_scree_gg <- function(interactive) {
      validate(need(!is.null(pca_shown$value()), "Render to see the variance explained."))
      v <- pca_shown$value()
      hi <- as.integer(c(sub("PC", "", input$pc_x %||% "PC1"),
                         sub("PC", "", input$pc_y %||% "PC2")))
      .pca_scree_plot(v$var_pct, highlight = hi, interactive = interactive)
    }
    dual_plot("scree", build_scree_gg,
              n_elements = reactive({ v <- pca_shown$value(); if (is.null(v)) 0L else min(10L, v$n_pc) }),
              height = "220px")

    invisible(NULL)
  })
}
