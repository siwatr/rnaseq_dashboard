# Project-wide colour palette configuration.
#
# The single source of truth for discrete colour mappings: the user opts in per
# metadata column, chooses a base palette (a two-level type -> name selector),
# and may hand-edit individual level colours. The config lives in `state$palette`
# (a UI preference, untouched by load/reset) and is read by every plot via
# palette_helpers.R -- the QC ggplots (group scales) and the ComplexHeatmap
# annotations (qc_annotation_colors()) now agree.
#
# Config shape per colData column:
#   list(type, name, colors)  where `colors` is the FULL level -> hex map.
# Choosing a (type, name) regenerates `colors` from that palette; hand-editing a
# level updates one entry and flips type -> "Custom" (so the selector no longer
# implies the named palette is being modified). Resolution + the palette catalogue
# live in palette_helpers.R.
#
# P3g-a wires the colData group only. rowData / assays / Other groups, continuous
# palettes, and config import/export land in P3g-b.

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

# Reference swatch: a small gradient bar for a palette (discrete bands for
# qualitative, a smooth ramp for sequential/divergent). Static / pure.
.pal_ref_swatch <- function(name, type) {
  cols <- palette_colors(type, name, if (identical(type, "Qualitative")) 8L else 9L)
  if (!length(cols)) return(NULL)
  tags$div(class = "mb-2",
    tags$div(class = "small mb-1", name),
    tags$div(style = sprintf(
      "height:22px;border-radius:4px;border:1px solid var(--bs-border-color);background:linear-gradient(to right, %s);",
      paste(cols, collapse = ", "))))
}

# The Reference tab: every catalogue palette grouped by type. Pure/static.
.pal_reference_ui <- function() {
  panel <- function(type) {
    nm <- palette_names(type)
    bslib::accordion_panel(type,
      if (length(nm)) lapply(nm, .pal_ref_swatch, type = type)
      else tags$p(class = "text-muted small", "Install the source package to preview these."))
  }
  tagList(
    tags$p(class = "text-muted small",
           "Reference swatches for the built-in palettes (qualitative shown as 8 bands; sequential/divergent as a ramp)."),
    bslib::accordion(panel("Qualitative"), panel("Sequential"), panel("Divergent"),
                     open = "Qualitative", multiple = TRUE)
  )
}

mod_palette_ui <- function(id) {
  ns <- NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      title = tags$h4("Palette", class = "fs-6 mb-0"), width = 320,
      helpText(paste("Set colour conventions per metadata attribute. Configured",
                     "colours feed both the ggplot plots and the ComplexHeatmap",
                     "annotations, so colours stay consistent across the app.")),
      uiOutput(ns("add_ui")),
      actionButton(ns("collapse_all"), "Collapse all", icon = icon("compress"),
                   class = "btn-sm btn-outline-secondary"),
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
        tags$h4("Reference", class = "fs-6 mb-0"),
        .pal_reference_ui()
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
    # rebuild on those structural events only. Per-level edits and palette
    # changes write state and push values back to the controls via update*()
    # (panels read the config via isolate()), so controls never reset mid-edit.
    struct <- reactiveVal(0L)
    bump <- function() struct(shiny::isolate(struct()) + 1L)

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
      choices <- setdiff(discrete, names(state$palette$colData))
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

    # Per-column observers (per-level pickers + palette type/name + reset +
    # remove + preview), registered once per column. Observers capture the
    # column's levels, so they are destroyed + re-registered whenever the dataset
    # changes (a persisted config may now map a column with different levels). All
    # observers use ignoreInit so a stale input value never re-fires on (re)registration.
    registered  <- new.env(parent = emptyenv())   # key -> TRUE once registered
    obs_handles <- new.env(parent = emptyenv())    # key -> list of observer handles

    unregister_col <- function(key) {
      h <- obs_handles[[key]]
      if (!is.null(h)) lapply(h, function(o) if (!is.null(o)) o$destroy())
      obs_handles[[key]] <- NULL
      registered[[key]] <- NULL
    }

    # Write a single field of a column's config (wholesale list assign).
    set_cfg <- function(col, field, value) {
      p <- state$palette; p$colData[[col]][[field]] <- value; state$palette <- p
    }

    register_col <- function(col) {
      key <- .pal_safe(col)
      if (!is.null(registered[[key]])) return(invisible())
      registered[[key]] <- TRUE
      # Levels are read *live* (cur_lvls) inside every observer, so the observers
      # stay correct if a later dataset reload changes this column's levels -- no
      # destroy/re-register dance (which races with the deferred data_version
      # flush). Only the per-level observer *count* is fixed at registration: a
      # reload that adds levels beyond it leaves the extra levels palette-filled
      # (re-add the column to hand-edit them).
      cur_lvls <- function() .pal_levels(shiny::isolate(coldata_df())[[col]])
      n0 <- length(cur_lvls())

      # Apply a named (non-Custom) palette: regenerate the full colour map and
      # push every picker to its new colour.
      apply_named <- function(type, name) {
        lvls <- cur_lvls()
        cols <- palette_discrete(lvls, NULL, type, name)
        p <- state$palette
        p$colData[[col]]$type <- type
        p$colData[[col]]$name <- name
        p$colData[[col]]$colors <- cols
        state$palette <- p
        for (i in seq_along(lvls)) update_picker(pin_id(col, i), unname(cols[[lvls[i]]]))
      }

      # Palette type: repopulate the name choices; for a real palette also apply
      # its first entry. "Custom" keeps the current colours (user-defined).
      type_obs <- observeEvent(input[[paste0("type_", key)]], {
        type <- input[[paste0("type_", key)]]
        if (identical(type, "Custom")) {
          set_cfg(col, "type", "Custom"); set_cfg(col, "name", "Custom palette")
          updateSelectInput(session, paste0("name_", key),
                            choices = "Custom palette", selected = "Custom palette")
        } else {
          nm <- palette_names(type)
          updateSelectInput(session, paste0("name_", key), choices = nm, selected = nm[1])
          apply_named(type, nm[1])
        }
      }, ignoreInit = TRUE)

      # Palette name within the current (non-Custom) type.
      name_obs <- observeEvent(input[[paste0("name_", key)]], {
        if (identical(state$palette$colData[[col]]$type, "Custom")) return()
        apply_named(state$palette$colData[[col]]$type %||% "Qualitative",
                    input[[paste0("name_", key)]])
      }, ignoreInit = TRUE)

      # One observer per level: a hand-edit updates that colour and flips the
      # palette to "Custom" (the colour map already holds every level).
      pin_obs <- lapply(seq_len(n0), function(i) {
        observeEvent(input[[pin_id(col, i)]], {
          if (is.null(state$palette$colData[[col]])) return()
          lvls <- cur_lvls(); if (i > length(lvls)) return()
          lev <- lvls[i]
          hex <- norm_color(input[[pin_id(col, i)]])
          if (is.na(hex)) return()
          p <- state$palette
          cols <- p$colData[[col]]$colors
          if (identical(unname(cols[[lev]]), unname(hex))) return()       # no change
          cols[[lev]] <- unname(hex)
          p$colData[[col]]$colors <- cols
          was_named <- !identical(p$colData[[col]]$type, "Custom")
          if (was_named) { p$colData[[col]]$type <- "Custom"
                           p$colData[[col]]$name <- "Custom palette" }
          state$palette <- p
          if (was_named) {                                  # reflect in the selectors
            updateSelectInput(session, paste0("type_", key), selected = "Custom")
            updateSelectInput(session, paste0("name_", key),
                              choices = "Custom palette", selected = "Custom palette")
          }
        }, ignoreInit = TRUE)
      })

      # Reset: back to the default palette (Qualitative / Okabe-Ito).
      reset_obs <- observeEvent(input[[paste0("reset_", key)]], {
        updateSelectInput(session, paste0("type_", key), selected = "Qualitative")
        updateSelectInput(session, paste0("name_", key),
                          choices = palette_names("Qualitative"), selected = "Okabe-Ito")
        apply_named("Qualitative", "Okabe-Ito")
      }, ignoreInit = TRUE)

      remove_obs <- observeEvent(input[[paste0("remove_", key)]], {
        p <- state$palette; p$colData[[col]] <- NULL; state$palette <- p
        unregister_col(key)                              # drop this column's observers
        bump()
      }, ignoreInit = TRUE)

      obs_handles[[key]] <- c(pin_obs, list(type_obs, name_obs, reset_obs, remove_obs))

      # Live preview swatch row (reads state$palette, so it updates on every edit).
      output[[paste0("preview_", key)]] <- renderPlot({
        cd <- coldata_df()
        if (!col %in% colnames(cd)) return(NULL)
        cfg <- state$palette$colData[[col]]; req(cfg)
        lv <- .pal_levels(cd[[col]])
        cols <- palette_discrete(lv, cfg$colors, cfg$type %||% "Qualitative",
                                 cfg$name %||% "Okabe-Ito", cfg$custom)
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
      invisible()
    }

    # --- Add: seed a default discrete config, then register its observers ----
    observeEvent(input$add_btn, {
      col <- input$add_col
      req(col, is.null(state$palette$colData[[col]]))
      lvls <- .pal_levels(coldata_df()[[col]])
      p <- state$palette
      p$colData[[col]] <- list(type = "Qualitative", name = "Okabe-Ito",
                               colors = palette_discrete(lvls, NULL, "Qualitative", "Okabe-Ito"))
      state$palette <- p
      register_col(col)
      bump()
    })

    # Collapse all open accordion panels.
    observeEvent(input$collapse_all, bslib::accordion_panel_close(ns("acc"), values = TRUE))

    # Dataset change: rebuild the panels so a surviving config's pickers reflect
    # the new levels. Observers read levels live (see register_col), so they need
    # no re-registration; configs for columns absent from the new dataset simply
    # don't render (the config persists, harmless).
    observeEvent(state$data_version, bump(), ignoreNULL = FALSE)

    # --- Panels: one accordion item per configured colData column ------------
    output$panels <- renderUI({
      struct()                                  # rebuild on add / remove / dataset change
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
      do.call(bslib::accordion,
              c(list(id = ns("acc"), open = TRUE, multiple = TRUE), panels))
    })

    invisible(NULL)
  })
}

# Build the accordion panel for one discrete colData column. Controls are seeded
# from the resolved colours (the config's `colors`, filled from its palette).
.palette_panel <- function(ns, col, x, cfg, has_picker, pin_id) {
  lvls <- .pal_levels(x)
  type <- cfg$type %||% "Qualitative"
  name <- cfg$name %||% "Okabe-Ito"
  resolved <- palette_discrete(lvls, cfg$colors, type, name, cfg$custom)
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
    bslib::layout_columns(
      col_widths = c(5, 7),
      selectInput(ns(paste0("type_", key)), "Palette type",
                  choices = palette_type_names(), selected = type),
      selectInput(ns(paste0("name_", key)), "Palette",
                  choices = palette_names(type), selected = name)),
    if (length(lvls) > .pal_many_levels)
      tags$div(class = "alert alert-warning py-1 px-2 small",
               sprintf("%d levels - colours are interpolated; pinning each is tedious.",
                       length(lvls))),
    tags$div(class = "small fw-semibold mb-1", "Levels"),
    tags$div(class = if (has_picker) "d-flex flex-wrap gap-2 mb-2" else "mb-2",
             lapply(seq_along(lvls), one_picker)),
    plotOutput(ns(paste0("preview_", key)), height = "90px"),
    tags$div(class = "d-flex gap-2 mt-2",
      bslib::tooltip(
        actionButton(ns(paste0("reset_", key)), "Reset to palette",
                     icon = icon("rotate-left"), class = "btn-sm btn-outline-secondary"),
        "Reset to the default palette (Qualitative / Okabe-Ito)."),
      actionButton(ns(paste0("remove_", key)), "Remove mapping",
                   icon = icon("trash"), class = "btn-sm btn-outline-danger"))
  )
  bslib::accordion_panel(title = col, value = col, body)
}
