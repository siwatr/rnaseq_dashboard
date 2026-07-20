# Page 1, "Size factors" tab (Dataset group, after Assay). Size factors are a
# DESeq2 median-of-ratios normalization -- deliberately separate from the
# CPM/TPM/FPKM assays (which are library-size/length normalizations that do NOT
# use size factors). A `navset_card_pill` of three pills:
#   Estimate    - choose what to estimate on (endogenous / spike-in / a custom
#                 gene set / all genes) + the estimator `type`, then estimate.
#                 The config is carried on the dds (metadata) so it survives edits
#                 and structural re-estimations reuse it. DE / PCA / Expression
#                 are consumers.
#   Per-sample  - the current per-sample size factors as a plot (bar per sample,
#                 or points grouped by a colData variable).
#   Compare     - a CONSUMER-ONLY scatter of two size-factor vectors computed
#                 under two different configs (nothing written back to the dds),
#                 to assess how similar two normalizations are (e.g. spike-in vs
#                 endogenous). A structural dds edit re-estimates both (the cache
#                 + deferred sig carry data_version).
# Backed by R/assay_helpers.R (estimate_size_factors + sizefactor_config).

# Shared choice vectors (Estimate + Compare pills).
.sf_control_choices <- c("Endogenous (default)" = "endogenous",
                         "Spike-in" = "spike_in", "Custom set" = "custom",
                         "All genes (discouraged)" = "all_genes")
.sf_type_choices <- c("ratio (median-of-ratios)" = "ratio",
                      "poscounts" = "poscounts", "iterate" = "iterate")

# Build a user-set size-factor config from sidebar inputs. custom_ids only matter
# for control == "custom".
.sf_config_from <- function(control, custom_ids = character(0), type = "ratio") {
  ctrl <- control %||% "endogenous"
  list(control = ctrl,
       custom_ids = if (identical(ctrl, "custom")) custom_ids %||% character(0) else character(0),
       type = type %||% "ratio", provenance = "user")
}

# Human control-set label (status text): "endogenous genes", "spike-in genes", ...
.sf_control_label <- function(control) {
  switch(control %||% "endogenous",
    endogenous = "endogenous genes", spike_in = "spike-in genes",
    custom = "a custom gene set", all_genes = "all genes", control)
}

# Compact "<control> / <type>" label for a Compare axis title.
.sf_config_label <- function(cfg) {
  short <- switch(cfg$control %||% "endogenous",
    endogenous = "Endogenous", spike_in = "Spike-in",
    custom = sprintf("Custom (%d)", length(cfg$custom_ids %||% character(0))),
    all_genes = "All genes", cfg$control)
  sprintf("%s / %s", short, cfg$type %||% "ratio")
}

# The Compare scatter: x/y are two size-factor vectors, a dashed x=y reference
# line + a single linear fit (no per-colour grouping). `colour_scale` is a
# pre-built ggplot scale (or NULL). `interactive` adds the plotly-only hover aes.
.sf_compare_plot <- function(df, xlab, ylab, subtitle, colour_lab = NULL,
                             colour_scale = NULL, fixed_ratio = FALSE,
                             point_size = 3, point_alpha = 0.9,
                             dark_theme = FALSE, interactive = FALSE) {
  aes_args <- list(x = quote(x), y = quote(y))
  if ("col" %in% names(df)) aes_args$colour <- quote(col)
  if (interactive)          aes_args$text   <- quote(text)
  p <- ggplot2::ggplot(df, do.call(ggplot2::aes, aes_args)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
    ggplot2::geom_point(size = point_size, alpha = point_alpha) +
    # One overall fit (inherit.aes = FALSE so a colour aes does not split it).
    ggplot2::geom_smooth(ggplot2::aes(x = x, y = y), data = df, inherit.aes = FALSE,
                         method = "lm", formula = y ~ x, se = FALSE,
                         colour = "#1f77b4", linewidth = 0.8) +
    ggplot2::labs(x = xlab, y = ylab, subtitle = subtitle, colour = colour_lab) +
    .plot_theme(dark_theme) +
    ggplot2::theme(legend.position = "right")
  if (!is.null(colour_scale)) p <- p + colour_scale
  if (isTRUE(fixed_ratio))    p <- p + ggplot2::coord_fixed()
  p
}

mod_sizefactors_ui <- function(id) {
  ns <- NS(id)
  bslib::navset_card_pill(
    # ---- Estimate -----------------------------------------------------------
    bslib::nav_panel(
      "Estimate",
      bslib::layout_sidebar(
        sidebar = bslib::sidebar(
          title = tags$h4("Estimate size factors", class = "fs-6 mb-0"), width = 300,
          radioButtons(ns("sf_control"), "Estimate using:", .sf_control_choices,
                       selected = "endogenous"),
          conditionalPanel(
            sprintf("input['%s'] == 'custom'", ns("sf_control")),
            gene_search_ui(ns, "sf", multiple = TRUE,
                           search_modes = c("exact", "contains", "regex"),
                           label = "Control-gene set",
                           placeholder = "e.g. Actb, Gapdh, ENSG...")),
          # All types honor the control set: estimate_size_factors() estimates on
          # the control-gene row-subset (so DESeq2's "iterate", which ignores
          # controlGenes, still respects it) and the full dds inherits the factors.
          selectInput(ns("sf_type"), "Estimator (type)", .sf_type_choices,
                      selected = "ratio"),
          conditionalPanel(
            sprintf("input['%s'] != 'endogenous'", ns("sf_control")),
            tags$div(class = "small text-warning mt-1",
                     "This normalization is used by DE, PCA, and Expression (not just QC). Spike-in factors can be noisy with few or low-count spikes; 'all genes' ignores the spike-in/exogenous separation and is discouraged.")),
          actionButton(ns("sf_estimate"), "Estimate size factors", class = "btn-primary"),
          helpText(class = "small text-muted mt-2",
                   "Size factors normalize for sequencing depth (used by DE and normalized log-counts). Re-estimating changes only samples' scaling, not the counts.")
        ),
        bslib::card(
          bslib::card_header(tags$h4("Current size factors", class = "fs-6 mb-0")),
          uiOutput(ns("status")),
          DT::DTOutput(ns("table"))
        )
      )
    ),
    # ---- Per-sample ---------------------------------------------------------
    bslib::nav_panel(
      "Per-sample",
      bslib::layout_sidebar(
        sidebar = bslib::sidebar(
          title = tags$h4("Per-sample factors", class = "fs-6 mb-0"), width = 300,
          uiOutput(ns("pers_x_ui")),
          uiOutput(ns("pers_colour_ui")),
          sliderInput(ns("pers_size"), "Point size", min = 0.25, max = 10,
                      value = 2.5, step = 0.25),
          sliderInput(ns("pers_alpha"), "Point opacity", min = 0.1, max = 1,
                      value = 0.9, step = 0.05),
          plot_subset_ui(ns, "pers"),
          uiOutput(ns("pers_auto_ui")),
          actionButton(ns("pers_render"), "Render", class = "btn-primary")
        ),
        bslib::card(
          bslib::card_header(tags$h4("Size factor per sample", class = "fs-6 mb-0")),
          uiOutput(ns("pers_stale")), .plot_dual(ns("pers_container"))
        )
      )
    ),
    # ---- Compare (consumer-only, last pill) ---------------------------------
    bslib::nav_panel(
      "Compare",
      bslib::layout_sidebar(
        sidebar = bslib::sidebar(
          title = tags$h4("Compare methods", class = "fs-6 mb-0"), width = 320,
          helpText(class = "small text-muted",
                   "Compare two size-factor methods. Nothing here is saved to the dataset - this is a view only."),
          bslib::accordion(
            open = c("X-axis method", "Y-axis method"),
            bslib::accordion_panel(
              "X-axis method", icon = icon("arrows-left-right"),
              radioButtons(ns("cmp_x_control"), "Estimate using:", .sf_control_choices,
                           selected = "endogenous"),
              conditionalPanel(
                sprintf("input['%s'] == 'custom'", ns("cmp_x_control")),
                gene_search_ui(ns, "sfx", multiple = TRUE,
                               search_modes = c("exact", "contains", "regex"),
                               label = "Control-gene set (X)")),
              selectInput(ns("cmp_x_type"), "Estimator (type)", .sf_type_choices,
                          selected = "ratio")
            ),
            bslib::accordion_panel(
              "Y-axis method", icon = icon("arrows-up-down"),
              radioButtons(ns("cmp_y_control"), "Estimate using:", .sf_control_choices,
                           selected = "spike_in"),
              conditionalPanel(
                sprintf("input['%s'] == 'custom'", ns("cmp_y_control")),
                gene_search_ui(ns, "sfy", multiple = TRUE,
                               search_modes = c("exact", "contains", "regex"),
                               label = "Control-gene set (Y)")),
              selectInput(ns("cmp_y_type"), "Estimator (type)", .sf_type_choices,
                          selected = "ratio")
            ),
            bslib::accordion_panel(
              "Plot settings", icon = icon("sliders"),
              bslib::input_switch(ns("cmp_fixed"), "Fix 1:1 aspect ratio", value = TRUE),
              uiOutput(ns("cmp_colour_ui")),
              sliderInput(ns("cmp_size"), "Point size", min = 0.25, max = 10,
                          value = 3, step = 0.25),
              sliderInput(ns("cmp_alpha"), "Point opacity", min = 0.1, max = 1,
                          value = 0.9, step = 0.05)
            )
          ),
          plot_subset_ui(ns, "cmp"),
          uiOutput(ns("cmp_auto_ui")),
          actionButton(ns("cmp_render"), "Render", icon = icon("play"), class = "btn-primary")
        ),
        bslib::card(
          bslib::card_header(tags$h4("Size-factor comparison", class = "fs-6 mb-0")),
          uiOutput(ns("cmp_stale")), .plot_dual(ns("cmp_container"))
        )
      )
    )
  )
}

#' @param state the shared app-state object (see [new_app_state()]).
#' @param dark_mode A reactive returning `TRUE` in dark mode (drives plot contrast).
#' @return Invisible NULL.
mod_sizefactors_server <- function(id, state, dark_mode = reactive(FALSE)) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    dark <- function() isTRUE(dark_mode())
    eng <- plot_engine_server(input, output, session, state)
    dual_plot <- eng$dual_plot; deferred <- eng$deferred; stale_note <- eng$stale_note
    showing_samples <- plot_subset_server(input, output, session, state,
                                          suffixes = c("pers", "cmp"))

    # Gene-search instances: the Estimate custom set + the two Compare custom sets.
    sf_search  <- gene_search_server(input, output, session, state, "sf",
                                     multiple = TRUE,
                                     search_modes = c("exact", "contains", "regex"))
    sfx_search <- gene_search_server(input, output, session, state, "sfx",
                                     multiple = TRUE,
                                     search_modes = c("exact", "contains", "regex"))
    sfy_search <- gene_search_server(input, output, session, state, "sfy",
                                     multiple = TRUE,
                                     search_modes = c("exact", "contains", "regex"))

    coldata <- function() as.data.frame(SummarizedExperiment::colData(state$working))

    # ==== Estimate pill =====================================================
    # Reflect the dataset's stored config in the inputs (on load / any data change).
    observeEvent(state$working, {
      req(state$working)
      cfg <- sizefactor_config(state$working)
      updateRadioButtons(session, "sf_control", selected = cfg$control)
      updateSelectInput(session, "sf_type", selected = cfg$type)
    })

    pending_config <- reactive({
      .sf_config_from(input$sf_control, sf_search()$ids, input$sf_type)
    })

    observeEvent(input$sf_estimate, {
      req(state$working)
      cfg <- pending_config()
      if (identical(cfg$control, "custom") && !length(cfg$custom_ids)) {
        showNotification("Enter at least one control gene for a custom set.", type = "error")
        return()
      }
      # Always re-estimate (cheap) so the user can re-run at will; commit only when
      # values or the recorded config actually change, so a repeat of the same
      # known config does not needlessly invalidate downstream caches. A loaded
      # (unknown) config always commits -- the fresh estimate differs from the
      # inherited vector, and this is how the user replaces it with a known one.
      cand <- tryCatch(estimate_size_factors(state$working, cfg),
                       error = function(e) { showNotification(conditionMessage(e), type = "error"); NULL })
      if (is.null(cand)) return()
      cur <- sizefactor_config(state$working)
      cur_sf <- tryCatch(DESeq2::sizeFactors(state$working), error = function(e) NULL)
      unchanged <- identical(cur$provenance, "user") &&
                   identical(cur$control, cfg$control) && setequal(cur$custom_ids, cfg$custom_ids) &&
                   identical(cur$type, cfg$type) && !is.null(cur_sf) &&
                   isTRUE(all.equal(unname(cur_sf), unname(DESeq2::sizeFactors(cand))))
      if (unchanged) {
        showNotification("Size factors re-estimated - values unchanged.", type = "message")
        return()
      }
      state_mutate(state, function(d) cand,
                   action = list(action = "size_factors", control = cfg$control, type = cfg$type,
                                 n_control = length(.sf_control_index(state$working, cfg))))
      showNotification("Size factors updated.", type = "message")
    })

    output$status <- renderUI({
      req(state$working)
      dds <- state$working
      sf <- tryCatch(DESeq2::sizeFactors(dds), error = function(e) NULL)
      if (is.null(sf))
        return(tags$p(class = "text-warning",
                      "Size factors are not set for this dataset. Estimate them below."))
      cfg <- sizefactor_config(dds)
      rng <- sprintf("Range %.3f-%.3f.", min(sf), max(sf))
      # For a loaded object we don't know the upstream control set / estimator, so
      # don't invent them -- only report what WE computed (auto default / user set).
      txt <- switch(cfg$provenance %||% "auto",
        loaded = sprintf(
          "Inherited from the loaded object; the control genes and estimator used upstream are unknown. %s Re-estimate below to set a known configuration.", rng),
        user = sprintf("Set on this tab from %s (%d control genes), estimator '%s'. %s",
                       .sf_control_label(cfg$control), length(.sf_control_index(dds, cfg)),
                       cfg$type %||% "ratio", rng),
        sprintf("Endogenous default (median-of-ratios), computed on load from %d endogenous genes. %s",
                length(.sf_control_index(dds, cfg)), rng))
      tags$p(class = "text-muted mb-2", txt)
    })

    output$table <- DT::renderDT({
      req(state$working)
      sf <- tryCatch(DESeq2::sizeFactors(state$working), error = function(e) NULL)
      validate(need(!is.null(sf), "No size factors to show."))
      df <- data.frame(sample = names(sf) %||% as.character(seq_along(sf)),
                       size_factor = round(as.numeric(sf), 4),
                       stringsAsFactors = FALSE)
      dt_table(df)
    })

    # ==== Per-sample pill ===================================================
    output$pers_x_ui <- renderUI({
      req(state$working)
      selectInput(ns("pers_x"), "X-axis",
                  choices = c("Sample" = "__sample__",
                              group_field_choices(colnames(coldata()), none = FALSE)),
                  selected = "__sample__")
    })
    output$pers_colour_ui <- renderUI({
      req(state$working)
      selectInput(ns("pers_colour"), "Colour by",
                  choices = aes_choices(aes_catalog(state), none = TRUE),
                  selected = "__none__")
    })
    output$pers_auto_ui <- renderUI({
      req(state$working)
      bslib::input_switch(ns("pers_auto"), "Auto-render", value = ncol(state$working) <= 200L)
    })

    # Colour aesthetic (optional): resolve via the shared attribute resolver; a
    # NULL result (e.g. removal flags not ready) simply means "no colour".
    pers_colour_resolve <- function(samples) {
      sel <- input$pers_colour %||% "__none__"
      if (identical(sel, "__none__")) return(NULL)
      res <- aes_resolve(state, sel, samples)
      if (is.null(res)) return(NULL)
      list(values = res$values, discrete = identical(res$kind, "discrete"),
           lab = res$label, scale = res)
    }

    pers_shown <- deferred("pers_auto", "pers_render",
      reactive({ req(state$working); input$pers_x %||% "__sample__" }),
      sig = reactive(list(input$pers_x, showing_samples(), state$data_version)))
    output$pers_stale <- stale_note(pers_shown)

    build_pers_gg <- function(interactive) {
      validate(need(!is.null(pers_shown$value()),
                    "Click Render (or enable auto-render) to plot size factors."))
      req(state$working)
      sf <- tryCatch(DESeq2::sizeFactors(state$working), error = function(e) NULL)
      validate(need(!is.null(sf), "No size factors yet - estimate them on the Estimate pill."))
      shown <- intersect(colnames(state$working), showing_samples())
      validate(need(length(shown) > 0, "No samples in the current 'Showing' selection."))
      xsel <- input$pers_x %||% "__sample__"
      psize <- input$pers_size %||% 2.5; palpha <- input$pers_alpha %||% 0.9
      cr <- pers_colour_resolve(shown)
      df <- data.frame(sample = shown, sf = as.numeric(sf[shown]), stringsAsFactors = FALSE)

      if (identical(xsel, "__sample__")) {
        df$x <- factor(df$sample, levels = df$sample)
        aes_args <- list(x = quote(x), y = quote(sf))
        if (!is.null(cr) && cr$discrete) aes_args$fill <- quote(col)
        if (interactive) aes_args$text <- quote(text)
        if (!is.null(cr) && cr$discrete) df$col <- cr$values
        if (interactive) df$text <- sprintf("Sample: %s<br>Size factor: %s",
                                            df$sample, signif(df$sf, 4))
        p <- ggplot2::ggplot(df, do.call(ggplot2::aes, aes_args)) +
          ggplot2::geom_col(alpha = palpha) +
          ggplot2::labs(x = "sample", y = "Size factor",
                        fill = if (!is.null(cr) && cr$discrete) cr$lab else NULL) +
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
          .plot_theme(dark())
        if (!is.null(cr) && cr$discrete) {
          sc <- aes_ggplot_scale(cr$scale, "fill")
          if (!is.null(sc)) p <- p + sc
        }
        return(p)
      }

      # Grouped by a colData variable: points (jittered) of size factor per group.
      g <- coldata()[shown, xsel]
      df$grp <- if (is.numeric(g)) g else { gg <- as.character(g); gg[is.na(gg)] <- "NA"; factor(gg) }
      aes_args <- list(x = quote(grp), y = quote(sf))
      if (!is.null(cr)) aes_args$colour <- quote(col)
      if (interactive)  aes_args$text   <- quote(text)
      if (!is.null(cr)) df$col <- cr$values
      if (interactive)  df$text <- sprintf("Sample: %s<br>%s: %s<br>Size factor: %s",
                                           df$sample, xsel, as.character(df$grp), signif(df$sf, 4))
      pos <- if (is.factor(df$grp)) ggplot2::position_jitter(width = 0.15, height = 0)
             else ggplot2::position_identity()
      p <- ggplot2::ggplot(df, do.call(ggplot2::aes, aes_args)) +
        ggplot2::geom_point(size = psize, alpha = palpha, position = pos) +
        ggplot2::labs(x = xsel, y = "Size factor",
                      colour = if (!is.null(cr)) cr$lab else NULL) +
        .plot_theme(dark())
      if (!is.null(cr)) {
        sc <- aes_ggplot_scale(cr$scale, "colour")
        if (!is.null(sc)) p <- p + sc
      }
      p
    }
    dual_plot("pers", build_pers_gg, n_elements = reactive({
      length(intersect(colnames(state$working %||% list()), showing_samples()))
    }))

    # ==== Compare pill (consumer-only) ======================================
    output$cmp_colour_ui <- renderUI({
      req(state$working)
      selectInput(ns("cmp_colour"), "Colour by",
                  choices = aes_choices(aes_catalog(state), none = TRUE),
                  selected = "__none__")
    })
    output$cmp_auto_ui <- renderUI({
      req(state$working)
      bslib::input_switch(ns("cmp_auto"), "Auto-render", value = FALSE)
    })

    cfg_x <- reactive(.sf_config_from(input$cmp_x_control, sfx_search()$ids, input$cmp_x_type))
    cfg_y <- reactive(.sf_config_from(input$cmp_y_control, sfy_search()$ids, input$cmp_y_type))

    # Compute the two size-factor vectors READ-ONLY (state is never mutated). Each
    # is estimated on the working dds under its config; an empty spike-in / custom
    # set is carried as ok = FALSE (a validate() inside a deferred spec is
    # swallowed, so error state is data, not a condition). Cached in `derived`
    # keyed on the two configs + data_version, so a structural dds edit
    # (data_version bump) re-estimates both, while a mere re-view reuses the cache.
    compare_value <- reactive({
      req(state$working)
      cx <- cfg_x(); cy <- cfg_y()
      state_derive(state, "sf_compare", params = list(cx, cy), expr = function() {
        one <- function(cfg) tryCatch(
          DESeq2::sizeFactors(estimate_size_factors(state$working, cfg)),
          error = function(e) structure(NA_real_, err = conditionMessage(e)))
        sx <- one(cx); sy <- one(cy)
        if (!is.null(attr(sx, "err"))) return(list(ok = FALSE, msg = paste0("X axis: ", attr(sx, "err"))))
        if (!is.null(attr(sy, "err"))) return(list(ok = FALSE, msg = paste0("Y axis: ", attr(sy, "err"))))
        list(ok = TRUE, sf_x = sx, sf_y = sy,
             lab_x = .sf_config_label(cx), lab_y = .sf_config_label(cy))
      })
    })
    compare_shown <- deferred("cmp_auto", "cmp_render", compare_value,
      sig = reactive(list(cfg_x(), cfg_y(), state$data_version, showing_samples())))
    output$cmp_stale <- stale_note(compare_shown)

    cmp_colour_resolve <- function(samples) {
      sel <- input$cmp_colour %||% "__none__"
      if (identical(sel, "__none__")) return(NULL)
      res <- aes_resolve(state, sel, samples)
      if (is.null(res)) return(NULL)
      list(values = res$values, discrete = identical(res$kind, "discrete"),
           lab = res$label, scale = aes_ggplot_scale(res, "colour"))
    }

    build_cmp_gg <- function(interactive) {
      v <- compare_shown$value()
      validate(need(!is.null(v),
                    "Click Render (or enable auto-render) to compare size factors."))
      validate(need(isTRUE(v$ok), v$msg %||% "Cannot compute size factors."))
      shown <- intersect(names(v$sf_x), showing_samples())
      validate(need(length(shown) > 0, "No samples in the current 'Showing' selection."))
      df <- data.frame(sample = shown, x = as.numeric(v$sf_x[shown]),
                       y = as.numeric(v$sf_y[shown]), stringsAsFactors = FALSE)
      cr <- cmp_colour_resolve(shown)
      colour_lab <- NULL; colour_scale <- NULL
      if (!is.null(cr)) { df$col <- cr$values; colour_lab <- cr$lab; colour_scale <- cr$scale }
      r2 <- tryCatch(summary(stats::lm(y ~ x, df))$r.squared, error = function(e) NA_real_)
      subtitle <- if (is.finite(r2)) sprintf("Linear fit R2 = %.3f (dashed line = x=y)", r2)
                  else "Dashed line = x=y"
      if (interactive) {
        df$text <- sprintf("Sample: %s<br>%s: %s<br>%s: %s", df$sample,
                           v$lab_x, signif(df$x, 4), v$lab_y, signif(df$y, 4))
        if (!is.null(cr)) df$text <- paste0(df$text, "<br>", cr$lab, ": ",
                                            if (cr$discrete) as.character(df$col) else round(df$col, 3))
      }
      .sf_compare_plot(df, xlab = v$lab_x, ylab = v$lab_y, subtitle = subtitle,
                       colour_lab = colour_lab, colour_scale = colour_scale,
                       fixed_ratio = isTRUE(input$cmp_fixed),
                       point_size = input$cmp_size %||% 3, point_alpha = input$cmp_alpha %||% 0.9,
                       dark_theme = dark(), interactive = interactive)
    }
    dual_plot("cmp", build_cmp_gg, n_elements = reactive({
      v <- compare_shown$value()
      if (is.null(v) || !isTRUE(v$ok)) 0L else length(intersect(names(v$sf_x), showing_samples()))
    }))

    invisible(NULL)
  })
}
