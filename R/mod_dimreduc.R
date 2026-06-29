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
                              shape_scale = NULL, caption = NULL, point_size = 3,
                              fixed_ratio = FALSE, legend_pos = "right",
                              dark_theme = FALSE, interactive = FALSE) {
  aes_args <- list(x = quote(x), y = quote(y))
  if ("col" %in% names(df)) aes_args$colour <- quote(col)
  if ("shp" %in% names(df)) aes_args$shape  <- quote(shp)
  if (interactive)          aes_args$text   <- quote(text)
  p <- ggplot2::ggplot(df, do.call(ggplot2::aes, aes_args)) +
    ggplot2::geom_point(size = point_size, alpha = 0.9) +
    ggplot2::labs(x = xlab, y = ylab, subtitle = subtitle, caption = caption,
                  colour = colour_lab, shape = shape_lab) +
    .plot_theme(dark_theme) +
    # Override .plot_theme's bottom default (PCA legends can be wide / continuous).
    ggplot2::theme(legend.position = legend_pos %||% "right")
  if (!is.null(colour_scale)) p <- p + colour_scale
  if (!is.null(shape_scale))  p <- p + shape_scale
  if (isTRUE(fixed_ratio))    p <- p + ggplot2::coord_fixed()
  p
}

# Scree bar plot of %variance per PC (top `n_show`); the two plotted PCs highlighted.
.pca_scree_plot <- function(var_pct, highlight = integer(0), n_show = 10L,
                            dark_theme = FALSE, interactive = FALSE) {
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
    .plot_theme(dark_theme)
}

mod_dimreduc_ui <- function(id) {
  ns <- NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      title = tags$h3("Dim. reduction", class = "fs-6 mb-0"), width = 340,
      # Sidebar is grouped into accordions to manage clutter (Embedding controls
      # the cached PCA; Appearance + Plot aesthetics are cheap re-plots).
      bslib::accordion(
        open = c("Embedding", "Appearance"),
        bslib::accordion_panel(
          "Embedding", icon = icon("diagram-project"),
          helpText(class = "small",
                   "PCA of the top-variable genes (endogenous only). Pick an input, then Render."),
          uiOutput(ns("assay_ui")),
          uiOutput(ns("assay_note")),
          uiOutput(ns("log_ui")),
          numericInput(ns("n_top"), "Top variable genes", value = 500, min = 2, step = 50),
          uiOutput(ns("pc_ui")),
          uiOutput(ns("auto_ui")),
          actionButton(ns("render"), "Render", icon = icon("play"), class = "btn-primary")
        ),
        bslib::accordion_panel(
          "Appearance", icon = icon("palette"),
          uiOutput(ns("colour_ui")),
          uiOutput(ns("shape_ui")),
          tags$hr(class = "my-2"),
          tags$div(class = "small fw-semibold mb-1", "Gene expression colour"),
          uiOutput(ns("gene_block"))
        ),
        bslib::accordion_panel(
          "Plot aesthetics", icon = icon("sliders"),
          bslib::input_switch(ns("fixed_ratio"), "Fix 1:1 aspect ratio", value = FALSE),
          sliderInput(ns("point_size"), "Point size", min = 0.25, max = 5,
                      value = 3, step = 0.25),
          selectInput(ns("legend_pos"), "Legend position",
                      choices = c("Right" = "right", "Left" = "left", "Top" = "top",
                                  "Bottom" = "bottom", "None" = "none"),
                      selected = "right")
        )
      ),
      plot_subset_ui(ns, "pca")
    ),
    bslib::card(
      bslib::card_header(tags$h3("PCA", class = "fs-6 mb-0")),
      uiOutput(ns("pca_stale")),
      uiOutput(ns("gene_caption")),
      .plot_dual(ns("pca_container"))
    ),
    bslib::card(
      bslib::card_header(tags$h3("Scree", class = "fs-6 mb-0")),
      .plot_dual(ns("scree_container"))
    )
  )
}

#' @param state the shared app-state object (see [new_app_state()]).
#' @param dark_mode A reactive returning `TRUE` in dark mode (wired from the navbar
#'   `input_dark_mode`); drives the plot theme contrast.
#' @return Invisible NULL (consumes the dds; does not mutate it).
mod_dimreduc_server <- function(id, state, dark_mode = reactive(FALSE)) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    dark <- function() isTRUE(dark_mode())
    eng <- plot_engine_server(input, output, session, state)
    dual_plot <- eng$dual_plot; deferred <- eng$deferred; stale_note <- eng$stale_note
    showing_samples <- plot_subset_server(input, output, session, state, suffixes = "pca")

    feature_type <- function() state$meta$feature_type %||% "gene"
    coldata <- function() as.data.frame(SummarizedExperiment::colData(state$working))
    discrete_cols <- function() {
      cd <- coldata(); names(cd)[!vapply(cd, is.numeric, logical(1))]
    }
    # Discrete colData columns usable for *shape* (ggplot's shape palette caps at
    # 6 levels). NA counts as its own "NA" level (we never silently drop points).
    shape_cols <- function() {
      cd <- coldata()
      Filter(function(c) {
        lv <- unique(ifelse(is.na(cd[[c]]), "NA", as.character(cd[[c]])))
        length(lv) <= 6L
      }, discrete_cols())
    }
    # Per-sample QC metrics offered as continuous colour options ("This session").
    # Share the QC page's derived cache (keyed on data_version) so we don't recompute
    # scater::perCellQCMetrics on every tab visit.
    qc_metrics <- reactive({
      req(state$working)
      state_derive(state, "qc_metrics", params = list(),
                   expr = function() qc_per_sample_metrics(state$working))
    })
    .qc_metric_labels <- c(library_size = "Library size", detected = "Detected features",
                           pct_mito = "% mitochondrial", pct_spike = "% spike-in")
    # Default colour-by: the first DESeq2 design variable (e.g. ~ X + Y -> "X"),
    # then "condition", then none (covers objects with a ~1 / no usable design).
    default_colour_col <- function() {
      cols <- colnames(coldata())
      dv <- tryCatch(all.vars(DESeq2::design(state$working)), error = function(e) character(0))
      if (length(dv) && dv[1] %in% cols) return(dv[1])
      if ("condition" %in% cols) return("condition")
      "__none__"
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
    # Grouped colour selector (<optgroup>s like the Palette selector): dataset
    # metadata vs session-derived options (gene expression + per-sample QC metrics).
    output$colour_ui <- renderUI({
      req(state$working)
      cd <- colnames(coldata())
      groups <- list("General" = c("(none)" = "__none__"))
      if (length(cd)) groups[["Data metadata"]] <- stats::setNames(cd, cd)
      groups[["This session"]] <- c(
        "Gene expression" = "__gene__",
        stats::setNames(paste0("__qc__", names(.qc_metric_labels)),
                        unname(.qc_metric_labels)))
      selectInput(ns("colour_by"), "Colour by", choices = groups,
                  selected = default_colour_col())
    })
    output$shape_ui <- renderUI({
      req(state$working)
      sc <- shape_cols()
      ch <- c("(none)" = "__none__", stats::setNames(sc, sc))
      selectInput(ns("shape_by"), "Shape by (discrete, <=6 values)",
                  choices = ch, selected = "__none__")
    })

    # --- Gene-expression colour block (always embedded; helper text when off) --
    output$gene_block <- renderUI({
      if (!identical(input$colour_by %||% "__none__", "__gene__"))
        return(helpText(class = "small text-muted",
          "Select 'Gene expression' in 'Colour by' to enable feature colouring."))
      # Default the duplicate-columns switch ON for compact rowData (<=10 non-
      # logical columns), where allowing non-unique name columns is convenient.
      rd <- as.data.frame(SummarizedExperiment::rowData(state$working), optional = TRUE)
      n_nonlogical <- if (ncol(rd)) sum(!vapply(rd, is.logical, logical(1))) else 0L
      assays <- SummarizedExperiment::assayNames(state$working)
      assay_sel <- if ("logcounts" %in% assays) "logcounts" else assays[1]
      tagList(
        bslib::input_switch(ns("gene_dup"),
          bslib::tooltip(
            tags$span("Include columns with duplicate values"),
            "Off hides rowData columns whose values are not unique (e.g. gene names that repeat). When such a column is searched, the first matching feature is used."),
          value = n_nonlogical <= 10L),
        uiOutput(ns("gene_searchby_ui")),
        bslib::input_switch(ns("gene_ci"), "Case-insensitive search", value = FALSE),
        textInput(ns("gene"), label = NULL, placeholder = "e.g. Gene1"),
        uiOutput(ns("gene_hint")),
        selectInput(ns("gene_assay"), "Expression assay",
                    choices = stats::setNames(assays, assays), selected = assay_sel),
        selectInput(ns("gene_transform"), "Transformation",
                    choices = c("None" = "none", "log2" = "log2", "log10" = "log10"),
                    selected = "none"),
        helpText(class = "small text-muted",
                 "Use a log transform for linear assays; logcounts / VST are already log-scaled."),
        uiOutput(ns("gene_pseudo_ui"))
      )
    })
    # Search-by field choices depend on the duplicate-columns switch.
    output$gene_searchby_ui <- renderUI({
      req(identical(input$colour_by %||% "", "__gene__"), state$working)
      rd <- SummarizedExperiment::rowData(state$working)
      ch <- feature_search_choices(rd, include_duplicates = isTRUE(input$gene_dup))
      def <- paste0(feature_type(), "_name")
      sel <- if (def %in% ch) def else "__rownames__"
      cur <- shiny::isolate(input$gene_searchby)
      if (!is.null(cur) && cur %in% ch) sel <- cur   # preserve across dup-toggle
      selectInput(ns("gene_searchby"), "Search by", choices = ch, selected = sel)
    })
    # Pseudocount appears only when a log transform is chosen; its default tracks
    # whether the chosen assay is integer-valued (counts) vs continuous.
    output$gene_pseudo_ui <- renderUI({
      req(identical(input$colour_by %||% "", "__gene__"))
      if (identical(input$gene_transform %||% "none", "none")) return(NULL)
      m <- tryCatch(.qc_assay_matrix(state$working, input$gene_assay %||% "logcounts"),
                    error = function(e) NULL)
      is_int <- !is.null(m) && all(m == round(m), na.rm = TRUE)
      numericInput(ns("gene_pseudo"), "Pseudocount",
                   value = if (is_int) 1 else 0.5, min = 0, step = 0.5)
    })
    # PC-axis selectors, populated from the computed embedding's PC count.
    output$pc_ui <- renderUI({
      np <- pca_n_pc()
      pcs <- paste0("PC", seq_len(np))
      # Preserve the current pick across re-renders (a new embedding invalidates
      # this UI even when the PC count is unchanged); isolate to avoid a feedback dep.
      keep <- function(cur, default) if (!is.null(cur) && cur %in% pcs) cur else default
      curx <- shiny::isolate(input$pc_x); cury <- shiny::isolate(input$pc_y)
      bslib::layout_columns(
        col_widths = c(6, 6),
        selectInput(ns("pc_x"), "X axis", choices = pcs, selected = keep(curx, "PC1")),
        selectInput(ns("pc_y"), "Y axis", choices = pcs,
                    selected = keep(cury, if (np >= 2L) "PC2" else "PC1")))
    })

    # --- Embedding (cached + deferred) --------------------------------------
    log_on <- function() {
      a <- input$assay %||% ""
      isTRUE(input$log_transform) && !(a %in% c("vst", "norm_logcounts", "logcounts"))
    }
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

    # --- Gene search (debounced; field-aware) -------------------------------
    # Debounce the free-text gene box so the lookup + suggestion search don't fire
    # on every keystroke (the suggestion regex is O(n features) on a miss).
    gene_q <- debounce(reactive(input$gene %||% ""), 300)
    # The character vector searched for the current "Search by" field (rownames
    # for the feature-id sentinel, else the chosen rowData column as character).
    search_values <- reactive({
      req(state$working)
      field <- input$gene_searchby %||% "__rownames__"
      if (identical(field, "__rownames__")) return(rownames(state$working))
      rd <- as.data.frame(SummarizedExperiment::rowData(state$working), optional = TRUE)
      if (field %in% names(rd)) as.character(rd[[field]]) else rownames(state$working)
    })
    # The resolved feature for the current query/field (id + the actual matched
    # value + match count), or NULL when not colouring by gene / box is empty.
    # Shared by the legend label, the in-card caption, and the hint.
    gene_resolved <- reactive({
      req(identical(input$colour_by %||% "", "__gene__"), state$working)
      q <- trimws(gene_q()); if (!nzchar(q)) return(NULL)
      resolve_feature(q, search_values(), rownames(state$working),
                      case_insensitive = isTRUE(input$gene_ci))
    })
    # Inline feedback under the gene box: duplicate-hit note, "Did you mean ...?"
    # on a near miss, a refine hint when too broad, or a not-found message.
    output$gene_hint <- renderUI({
      req(identical(input$colour_by %||% "", "__gene__"))
      q <- trimws(gene_q()); if (!nzchar(q)) return(NULL)
      rf <- gene_resolved()
      note <- function(cls, ...) tags$div(class = paste("small mb-2", cls), ...)
      if (!is.null(rf) && !is.na(rf$id)) {
        if (rf$n > 1L)
          return(note("text-muted", sprintf("%d features matched '%s'; showing the first (%s).",
                                             rf$n, q, rf$match)))
        return(NULL)
      }
      sg <- suggest_features(q, search_values())  # always case-insensitive (independent of toggle)
      if (isTRUE(sg$over_cap))
        return(note("text-primary", "Too many partial matches - type more to narrow it down."))
      if (length(sg$suggestions)) {
        more <- sg$n_match - length(sg$suggestions)
        txt <- paste(sg$suggestions, collapse = ", ")
        if (more > 0L) txt <- paste0(txt, sprintf(" (+%d more)", more))
        return(note("text-primary", sprintf("Feature '%s' not found. Did you mean: %s?", q, txt)))
      }
      note("text-primary", sprintf("Feature '%s' not found.", q))
    })
    # In-card caption naming the feature + its true unique id (the dds rowname)
    # the expression colour is plotting, shown only when colouring by a gene.
    output$gene_caption <- renderUI({
      if (!identical(input$colour_by %||% "", "__gene__")) return(NULL)
      rf <- gene_resolved()
      if (is.null(rf) || is.na(rf$id)) return(NULL)
      tags$div(class = "small text-muted px-2 pb-1",
               sprintf("Plotting expression of %s (%s)", rf$match, rf$id))
    })

    # --- Colour / shape resolution (consumes the Palette configs) -----------
    # Returns list(values, discrete, lab, scale) for the colour aesthetic, or NULL.
    colour_resolve <- function(samples) {
      sel <- input$colour_by %||% "__none__"
      if (identical(sel, "__none__")) return(NULL)
      dds <- state$working
      if (identical(sel, "__gene__")) {
        q <- trimws(gene_q())
        validate(need(nzchar(q), "Enter a feature to colour by."))
        rf <- gene_resolved()
        validate(need(!is.null(rf) && !is.na(rf$id), sprintf("Feature '%s' not found.", q)))
        assay_name <- input$gene_assay %||% "logcounts"
        # Guard against a stale assay name so the data shown matches the colourbar
        # label (.qc_assay_matrix would silently fall back to log2(CPM+1) otherwise).
        validate(need(assay_name %in% SummarizedExperiment::assayNames(dds),
                      sprintf("Assay '%s' is no longer available.", assay_name)))
        m <- .qc_assay_matrix(dds, assay_name)
        raw <- as.numeric(m[rf$id, samples])
        tf <- input$gene_transform %||% "none"
        pc <- input$gene_pseudo %||% 1
        vals <- expr_transform(raw, tf, pc)
        suffix <- if (identical(tf, "none")) assay_name else paste0(assay_name, ", ", tf)
        # Label with the actual matched value (e.g. "Duxf3"), not the typed query.
        return(list(values = vals, discrete = FALSE, lab = paste0(rf$match, "\n(", suffix, ")"),
                    scale = .continuous_scale_from(state$palette$assays[[assay_name]], vals)))
      }
      if (startsWith(sel, "__qc__")) {
        metric <- sub("^__qc__", "", sel)
        qm <- qc_metrics()
        vals <- as.numeric(qm[samples, metric])
        return(list(values = vals, discrete = FALSE,
                    lab = unname(.qc_metric_labels[metric]) %||% metric, scale = NULL))
      }
      cd <- coldata()[samples, , drop = FALSE]
      x <- cd[[sel]]
      if (is.numeric(x)) {
        return(list(values = x, discrete = FALSE, lab = sel,
                    scale = .continuous_scale_from(state$palette$colData[[sel]], x)))
      }
      lv <- as.character(x); lv[is.na(lv)] <- "NA"
      f <- factor(lv)
      list(values = f, discrete = TRUE, lab = sel,
           scale = .discrete_scale_from(state$palette$colData[[sel]], levels(f)))
    }
    # Colour scales from a Palette config, or NULL when none is set -- mirroring QC's
    # group_palette/.qc_group_scale, so a no-config colouring falls back to thematic's
    # theme palette (not a hard-coded Okabe-Ito / viridis).
    .discrete_scale_from <- function(cfg, levels) {
      if (is.null(cfg)) return(NULL)
      cols <- palette_discrete(levels, cfg$colors, cfg$name %||% "Okabe-Ito", cfg$custom)
      ggplot2::scale_colour_manual(values = cols, na.value = "grey70")
    }
    .continuous_scale_from <- function(cfg, values) {
      if (is.null(cfg)) return(NULL)
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
      colour_lab <- NULL; shape_lab <- NULL; colour_scale <- NULL; shape_scale <- NULL; cap <- NULL
      if (!is.null(cr)) { df$col <- cr$values; colour_lab <- cr$lab; colour_scale <- cr$scale }
      if (!is.null(sr) && !isTRUE(sr$skip)) { df$shp <- sr$values; shape_lab <- sr$lab }
      if (!is.null(sr) && isTRUE(sr$skip))   # ggplot's shape palette caps at 6 levels
        cap <- sprintf("Shape: '%s' has more than 6 groups - not shown.", sr$lab)
      pcpct <- function(p) sprintf("%s: %.1f%%", p, v$var_pct[as.integer(sub("PC", "", p))])
      if (interactive) {
        df$text <- sprintf("%s\n%s, %s", df$sample, pcpct(px), pcpct(py))
        if (!is.null(cr)) df$text <- paste0(df$text, "\n", cr$lab, ": ",
                                            if (cr$discrete) df$col else round(df$col, 3))
      }
      .pca_scatter_plot(df, xlab = pcpct(px), ylab = pcpct(py),
                        subtitle = paste0("Input: ", v$label,
                                          " | ", v$n_genes, " top-variable genes"),
                        colour_lab = colour_lab, shape_lab = shape_lab,
                        colour_scale = colour_scale, shape_scale = shape_scale,
                        caption = cap, point_size = input$point_size %||% 3,
                        fixed_ratio = isTRUE(input$fixed_ratio),
                        legend_pos = input$legend_pos %||% "right", dark_theme = dark(),
                        interactive = interactive)
    }
    dual_plot("pca", build_pca_gg, n_elements = reactive({
      v <- pca_shown$value(); if (is.null(v)) 0L else length(intersect(rownames(v$scores), showing_samples()))
    }))

    build_scree_gg <- function(interactive) {
      validate(need(!is.null(pca_shown$value()), "Render to see the variance explained."))
      v <- pca_shown$value()
      hi <- as.integer(c(sub("PC", "", input$pc_x %||% "PC1"),
                         sub("PC", "", input$pc_y %||% "PC2")))
      .pca_scree_plot(v$var_pct, highlight = hi, dark_theme = dark(), interactive = interactive)
    }
    dual_plot("scree", build_scree_gg,
              n_elements = reactive({ v <- pca_shown$value(); if (is.null(v)) 0L else min(10L, v$n_pc) }),
              height = "220px")

    invisible(NULL)
  })
}
