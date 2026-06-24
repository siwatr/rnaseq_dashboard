# Project-wide colour palette configuration.
#
# The single source of truth for discrete colour mappings: the user opts in per
# metadata column, assigns a base palette, and pins specific levels. The config
# lives in `state$palette` (a UI preference, untouched by load/reset) and is read
# by every plot via palette_helpers.R -- the QC ggplots (group scales) and the
# ComplexHeatmap annotations (qc_annotation_colors()) now agree.
#
# P3g-a wires the **colData** (sample annotation) group only, discrete columns
# only. rowData / assays / Other groups and continuous palettes + config
# import/export land in P3g-b.

# Distinct levels of a column, in plot order (factor levels, else sorted unique).
.pal_levels <- function(x) {
  if (is.factor(x)) levels(x) else sort(unique(as.character(x[!is.na(x)])))
}
# A column is configurable here when it is discrete (not numeric) with >= 1 level.
.pal_is_discrete <- function(x) !is.numeric(x) && length(.pal_levels(x)) >= 1L
# Sanitize a column name into an input-id-safe token.
.pal_safe <- function(col) gsub("[^A-Za-z0-9]", "_", col)
# Show a level-count warning above this many pickers.
.pal_many_levels <- 12L

mod_palette_ui <- function(id) {
  ns <- NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      title = tags$h4("Palette", class = "fs-6 mb-0"), width = 320,
      helpText(paste("Set colour conventions per metadata attribute. Configured",
                     "colours feed both the ggplot plots and the ComplexHeatmap",
                     "annotations, so colours stay consistent across the app.")),
      uiOutput(ns("add_ui")),
      helpText(class = "text-muted small mt-2",
               "Feature (rowData), assay, and other palettes - plus continuous",
               "scales and config import/export - arrive in a later slice.")
    ),
    bslib::navset_card_pill(
      bslib::nav_panel(
        tags$h4("Sample (colData)", class = "fs-6 mb-0"),
        uiOutput(ns("panels"))
      ),
      bslib::nav_panel(
        tags$h4("Feature (rowData)", class = "fs-6 mb-0"),
        tags$p(class = "text-muted p-2", "Coming in a later slice (P3g-b).")
      ),
      bslib::nav_panel(
        tags$h4("Assays", class = "fs-6 mb-0"),
        tags$p(class = "text-muted p-2", "Coming in a later slice (P3g-b).")
      ),
      bslib::nav_panel(
        tags$h4("Other", class = "fs-6 mb-0"),
        tags$p(class = "text-muted p-2", "Coming in a later slice (P3g-b).")
      )
    )
  )
}

#' @param state the shared app-state object (see [new_app_state()]).
#' @return Invisible NULL.
mod_palette_server <- function(id, state) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Structure trigger: bumped on add / remove / dataset-change so the panels
    # rebuild on those structural events only. Per-level picker edits and base-
    # palette changes write state and push values back to the pickers via
    # update*() (panels read the config via isolate()), so pickers never reset
    # mid-interaction and no observe loop forms.
    struct <- reactiveVal(0L)
    bump <- function() struct(shiny::isolate(struct()) + 1L)
    observeEvent(state$data_version, bump())   # dataset change -> columns/levels change

    has_picker <- requireNamespace("shinyWidgets", quietly = TRUE)
    pin_id <- function(col, i) paste0("pin_", .pal_safe(col), "_", i)
    update_picker <- function(i_id, value) {
      if (has_picker) shinyWidgets::updateColorPickr(session, i_id, value = value)
      else updateTextInput(session, i_id, value = value)
    }

    # colData columns + their discreteness, for the current dataset.
    coldata_df <- reactive({
      req(state$working)
      as.data.frame(SummarizedExperiment::colData(state$working))
    })

    # --- Sidebar: add-a-mapping control -------------------------------------
    output$add_ui <- renderUI({
      if (is.null(state$working)) {
        return(helpText(class = "text-muted", "Load a dataset to configure colours."))
      }
      cd <- coldata_df()
      discrete <- names(cd)[vapply(cd, .pal_is_discrete, logical(1))]
      configured <- names(state$palette$colData)
      choices <- setdiff(discrete, configured)
      tagList(
        selectInput(ns("add_col"), "Add colour mapping for",
                    choices = if (length(choices)) choices else character(0)),
        actionButton(ns("add_btn"), "Add", icon = icon("plus"),
                     class = "btn-sm btn-primary",
                     disabled = if (!length(choices)) NA else NULL),
        if (!length(choices) && length(discrete))
          helpText(class = "text-muted small mt-1", "All discrete columns are configured.")
      )
    })

    # --- Add: seed a default discrete config, then register its observers ----
    observeEvent(input$add_btn, {
      col <- input$add_col
      req(col, is.null(state$palette$colData[[col]]))
      p <- state$palette
      p$colData[[col]] <- list(type = "discrete", palette = "Okabe-Ito",
                               pins = stats::setNames(character(0), character(0)))
      state$palette <- p
      register_col(col)
      bump()
    })

    # Per-column observers (per-level pins + base palette + reset + remove +
    # preview), registered exactly once per column. A "pin" is recorded only when
    # a picker differs from the base palette's auto colour for that level, so the
    # base-palette dropdown keeps meaning. reset/palette-change push values to the
    # pickers via update*(), which is also what reverts them in the live UI.
    registered <- new.env(parent = emptyenv())
    register_col <- function(col) {
      key <- .pal_safe(col)
      if (!is.null(registered[[key]])) return(invisible())
      registered[[key]] <- TRUE
      lvls <- .pal_levels(shiny::isolate(coldata_df())[[col]])

      cur_palette <- function() state$palette$colData[[col]]$palette %||% "Okabe-Ito"

      # One observer per level: classify the picker value as pin vs auto-fill.
      lapply(seq_along(lvls), function(i) {
        observeEvent(input[[pin_id(col, i)]], {
          if (is.null(state$palette$colData[[col]])) return()
          raw <- input[[pin_id(col, i)]]
          auto <- palette_discrete(lvls, NULL, cur_palette())[[lvls[i]]]
          hex <- if (is.null(raw) || !nzchar(raw)) NA_character_ else norm_color(raw)
          p <- state$palette
          pins <- p$colData[[col]]$pins
          if (!is.na(hex) && !identical(unname(hex), unname(auto))) {
            pins[[lvls[i]]] <- unname(hex)               # user override -> pin
          } else {
            pins <- pins[setdiff(names(pins), lvls[i])]  # equals auto -> not a pin
          }
          p$colData[[col]]$pins <- pins
          state$palette <- p
        }, ignoreInit = TRUE)
      })

      # Base palette: persist + move every non-pinned picker to the new auto colour.
      observeEvent(input[[paste0("pal_", key)]], {
        val <- input[[paste0("pal_", key)]]
        if (is.null(val) || identical(val, cur_palette())) return()
        p <- state$palette; p$colData[[col]]$palette <- val; state$palette <- p
        pins <- names(p$colData[[col]]$pins)
        auto <- palette_discrete(lvls, NULL, val)
        for (i in seq_along(lvls)) if (!lvls[i] %in% pins)
          update_picker(pin_id(col, i), unname(auto[[lvls[i]]]))
      }, ignoreInit = TRUE)

      # Reset: clear pins + push every picker back to the base-palette auto colour.
      observeEvent(input[[paste0("reset_", key)]], {
        p <- state$palette
        p$colData[[col]]$pins <- stats::setNames(character(0), character(0))
        state$palette <- p
        auto <- palette_discrete(lvls, NULL, cur_palette())
        for (i in seq_along(lvls)) update_picker(pin_id(col, i), unname(auto[[lvls[i]]]))
      })

      observeEvent(input[[paste0("remove_", key)]], {
        p <- state$palette; p$colData[[col]] <- NULL; state$palette <- p
        bump()
      })

      # Live preview swatch row (reads state$palette, so it updates on every pick).
      output[[paste0("preview_", key)]] <- renderPlot({
        cd <- coldata_df()
        if (!col %in% colnames(cd)) return(NULL)
        cfg <- state$palette$colData[[col]]; req(cfg)
        lv <- .pal_levels(cd[[col]])
        cols <- palette_discrete(lv, cfg$pins, cfg$palette)
        df <- data.frame(level = factor(names(cols), levels = names(cols)), y = 1L)
        ggplot2::ggplot(df, ggplot2::aes(x = .data$level, y = .data$y, fill = .data$level)) +
          ggplot2::geom_col(width = 1) +
          ggplot2::scale_fill_manual(values = cols, guide = "none") +
          ggplot2::labs(x = NULL, y = NULL) +
          ggplot2::theme_minimal(base_size = 11) +
          ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                         axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                         panel.grid = ggplot2::element_blank())
      }, height = 90)
    }

    # --- Panels: one accordion item per configured colData column ------------
    output$panels <- renderUI({
      struct()                                  # rebuild on add/remove/palette change
      if (is.null(state$working)) {
        return(tags$p(class = "text-muted p-2", "Load a dataset to configure colours."))
      }
      cd <- shiny::isolate(coldata_df())
      cfg <- shiny::isolate(state$palette$colData)
      cols <- intersect(names(cfg), colnames(cd))
      if (!length(cols)) {
        return(tags$p(class = "text-muted p-2",
                      'No colour mappings yet. Use "Add colour mapping" in the sidebar.'))
      }
      panels <- lapply(cols, function(col) .palette_panel(ns, col, cd[[col]], cfg[[col]],
                                                          has_picker, pin_id))
      do.call(bslib::accordion, c(list(open = TRUE, multiple = TRUE), panels))
    })

    invisible(NULL)
  })
}

# Build the accordion panel for one discrete colData column. Pickers are seeded
# from the resolved colour (pin if present, else the base-palette auto colour).
.palette_panel <- function(ns, col, x, cfg, has_picker, pin_id) {
  lvls <- .pal_levels(x)
  resolved <- palette_discrete(lvls, cfg$pins, cfg$palette %||% "Okabe-Ito")
  key <- .pal_safe(col)

  one_picker <- function(i) {
    lev <- lvls[i]; val <- unname(resolved[[lev]]); pid <- ns(pin_id(col, i))
    if (has_picker) {
      shinyWidgets::colorPickr(pid, label = lev, selected = val, update = "save",
                               width = "100%", useAsButton = TRUE,
                               interaction = list(input = TRUE, save = TRUE, clear = FALSE))
    } else {
      tags$div(class = "d-flex align-items-center gap-2 mb-1",
        tags$span(style = sprintf(
          "display:inline-block;width:1.1rem;height:1.1rem;border-radius:3px;border:1px solid #888;background:%s;",
          val)),
        textInput(pid, label = NULL, value = val, width = "120px"),
        tags$span(class = "small text-muted", lev))
    }
  }

  body <- tagList(
    selectInput(ns(paste0("pal_", key)), "Base palette",
                choices = palette_qualitative_names(),
                selected = cfg$palette %||% "Okabe-Ito"),
    if (length(lvls) > .pal_many_levels)
      tags$div(class = "alert alert-warning py-1 px-2 small",
               sprintf("%d levels - colours are interpolated/recycled; pinning each is tedious.",
                       length(lvls))),
    tags$div(class = "small fw-semibold mb-1", "Levels"),
    tags$div(class = if (has_picker) "d-flex flex-wrap gap-2 mb-2" else "mb-2",
             lapply(seq_along(lvls), one_picker)),
    plotOutput(ns(paste0("preview_", key)), height = "90px"),
    tags$div(class = "d-flex gap-2 mt-2",
      bslib::tooltip(
        actionButton(ns(paste0("reset_", key)), "Reset to palette",
                     icon = icon("rotate-left"), class = "btn-sm btn-outline-secondary"),
        "Clear all pinned colours; every level follows the base palette."),
      actionButton(ns(paste0("remove_", key)), "Remove mapping",
                   icon = icon("trash"), class = "btn-sm btn-outline-danger"))
  )
  bslib::accordion_panel(title = col, value = col, body)
}
