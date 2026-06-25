# Project-wide colour palette configuration.
#
# The single source of truth for colour mappings, organized by the aspect of the
# dds object it colours: colData (sample annotations), rowData (feature
# annotations), assays (expression ramps), and "other" (app-internal maps:
# removal status, sample-correlation score). The config lives in `state$palette`
# (a UI preference, untouched by load/reset) and is read by every plot via
# palette_helpers.R -- the QC ggplots, the ComplexHeatmap annotations
# (qc_annotation_colors()), and the correlation-heatmap ramp now agree.
#
# Config shape: state$palette[[domain]][[item]] is either
#   discrete:   list(name, colors)          (full level -> hex map)
#   continuous: list(name, min, max, custom) (palette + anchors; min/max are a
#               number, a "p<pct>" percentile, or "" = the data range)
# Resolution + the palette catalogue live in palette_helpers.R. P3g-b wires the
# colData/rowData/assays/other groups + continuous palettes; JSON config
# import/export (P3g-c) and factor management (P3g-d) follow.

# Distinct levels of a column, in plot order (factor levels, else sorted unique).
.pal_levels <- function(x) {
  if (is.factor(x)) levels(x) else sort(unique(as.character(x[!is.na(x)])))
}
# Sanitize an item name into an input-id-safe token.
.pal_safe <- function(x) gsub("[^A-Za-z0-9]", "_", x)
# Show a level-count warning above this many pickers.
.pal_many_levels <- 12L

# Note listing discrete attributes hidden from the add control (over the cap).
.pal_hidden_note <- function(items, max_n) {
  more <- length(items) > 5L
  shown <- if (more) utils::head(items, 5L) else items
  sprintf("%d attribute(s) hidden (> %d unique values; raise option ddsdashboard.palette_max_levels to allow). %s%s%s",
          length(items), max_n,
          if (more) "Affected fields (first 5): " else "Affected fields: ",
          paste(shown, collapse = ", "), if (more) ", ..." else "")
}

# The four Setting pills: domain -> display label.
.pal_domains <- c(colData = "Sample", rowData = "Feature",
                  assays = "Assay", other = "Other")
# "other" maps are fixed, app-internal.
.pal_other_items <- c("removal_status", "correlation")
.pal_other_kind  <- c(removal_status = "discrete", correlation = "continuous")
.pal_removal_levels <- c("pass", "suggested_other", "suggested_this")

# Preview swatch for a palette. `discrete` -> equal-width solid blocks (no
# interpolation); otherwise a smooth gradient ramp. `name` is the resolvable
# palette name; the visible label is its clean form. Static / pure.
.pal_ref_swatch <- function(name, discrete) {
  cols <- palette_colors(name, if (discrete) 8L else 9L)
  if (!length(cols)) return(NULL)
  bar <- if (discrete) {
    tags$div(class = "d-flex",
      style = "height:22px;border-radius:4px;overflow:hidden;border:1px solid var(--bs-border-color);",
      lapply(cols, function(cl) tags$div(style = sprintf("flex:1;background:%s;", cl))))
  } else {
    tags$div(style = sprintf(
      "height:22px;border-radius:4px;border:1px solid var(--bs-border-color);background:linear-gradient(to right, %s);",
      paste(cols, collapse = ", ")))
  }
  tags$div(class = "mb-2", tags$div(class = "small mb-1", .pal_label(name)), bar)
}

# A gradient bar for a continuous palette name (or NULL). Pure.
.pal_gradient_bar <- function(name, custom = NULL, reverse = FALSE) {
  cols <- .continuous_stops(name, custom, reverse = reverse)
  tags$div(style = sprintf(
    "height:26px;border-radius:4px;border:1px solid var(--bs-border-color);background:linear-gradient(to right, %s);",
    paste(cols, collapse = ", ")))
}

# The Preview tab: every catalogue palette grouped, qualitative groups as
# discrete blocks. Pure/static; the accordion id (ns("ref_acc")) drives Collapse.
.pal_reference_ui <- function(ns) {
  panel <- function(type) {
    nm <- palette_names(type); disc <- .pal_type_discrete(type)
    bslib::accordion_panel(type,
      if (length(nm)) lapply(nm, .pal_ref_swatch, discrete = disc)
      else tags$p(class = "text-muted small", "Install the source package to preview these."))
  }
  groups <- setdiff(palette_type_names(), "Custom")   # Custom isn't a real preset
  tagList(
    tags$div(class = "d-flex justify-content-between align-items-center mb-2",
      tags$p(class = "text-muted small mb-0",
             "Swatches for the built-in palettes (qualitative shown as discrete blocks; sequential/divergent as a ramp)."),
      actionButton(ns("collapse_ref"), "Collapse all", icon = icon("compress"),
                   class = "btn-sm btn-outline-secondary")),
    do.call(bslib::accordion,
            c(list(id = ns("ref_acc"), multiple = TRUE, open = groups[1]),
              lapply(groups, panel)))
  )
}

mod_palette_ui <- function(id) {
  ns <- NS(id)
  blurb <- c(colData = "Set colours per sample-metadata column (colData). They feed the QC plots and the ComplexHeatmap annotations.",
             rowData = "Set colours per feature-metadata column (rowData). Used for heatmap row annotations (P5).",
             assays  = "Set the expression colour ramp per assay. Used for the expression heatmap / PCA gene colouring (P4/P5).",
             other   = "Recolour app-internal maps: the QC removal-status colours and the sample-correlation ramp.")
  pill <- function(dom) {
    label <- .pal_domains[[dom]]
    bslib::nav_panel(label,
      bslib::layout_sidebar(
        sidebar = bslib::sidebar(
          title = tags$h4(paste(label, "colours"), class = "fs-6 mb-0"), width = 310,
          helpText(blurb[[dom]]),
          uiOutput(ns(paste0("addui_", dom))),
          actionButton(ns(paste0("collapse_", dom)), "Collapse all",
                       icon = icon("compress"), class = "btn-sm btn-outline-secondary")),
        uiOutput(ns(paste0("panels_", dom)))))
  }
  bslib::navset_card_tab(
    title = tags$h3("Palette", class = "fs-6 mb-0 pe-3"),
    bslib::nav_panel(
      tags$h4("Setting", class = "fs-6"),
      do.call(bslib::navset_pill, lapply(names(.pal_domains), pill))
    ),
    bslib::nav_panel(
      tags$h4("Preview", class = "fs-6"),
      .pal_reference_ui(ns)
    )
  )
}

#' @param state the shared app-state object (see [new_app_state()]).
#' @return Invisible NULL.
mod_palette_server <- function(id, state) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    struct <- reactiveVal(0L)
    bump <- function() struct(shiny::isolate(struct()) + 1L)
    observeEvent(state$data_version, bump(), ignoreNULL = FALSE)

    has_picker <- requireNamespace("shinyWidgets", quietly = TRUE)
    update_picker <- function(i_id, value) {
      if (has_picker) shinyWidgets::updateColorPickr(session, i_id, value = value)
      else updateTextInput(session, i_id, value = value)
    }

    # --- Domain accessors (read state$working live; isolate at call sites) ----
    dom_frame <- function(dom) {
      if (dom == "colData") as.data.frame(SummarizedExperiment::colData(state$working))
      else as.data.frame(SummarizedExperiment::rowData(state$working))
    }
    dom_items <- function(dom) {
      switch(dom,
        colData = colnames(dom_frame("colData")),
        rowData = colnames(dom_frame("rowData")),
        assays  = SummarizedExperiment::assayNames(state$working),
        other   = .pal_other_items)
    }
    dom_needs_data <- function(dom) dom %in% c("colData", "rowData", "assays")
    dom_kind <- function(dom, item) {
      if (dom == "assays") return("continuous")
      if (dom == "other")  return(.pal_other_kind[[item]])
      if (is.numeric(dom_frame(dom)[[item]])) "continuous" else "discrete"
    }
    dom_levels <- function(dom, item) {
      if (dom == "other") return(.pal_removal_levels)
      .pal_levels(dom_frame(dom)[[item]])
    }
    # Underlying data class for the accordion badge (helps the future factor PR).
    dom_class <- function(dom, item) {
      if (dom == "assays") return("numeric")
      if (dom == "other")  return(if (item == "correlation") "numeric" else "factor")
      class(dom_frame(dom)[[item]])[1]
    }
    # Default config for a freshly added item. A few items carry presets: the
    # removal-status map keeps its QC green/amber/red, and the correlation ramp
    # defaults to a reversed RdBu (high correlation = red) anchored to [-1, 1].
    default_cfg <- function(dom, item) {
      if (dom == "other" && item == "removal_status")
        return(list(name = "Custom palette", colors = .removal_palette))
      if (dom == "other" && item == "correlation")
        return(list(name = "RColorBrewer: RdBu", min = "-1", max = "1",
                    reverse = TRUE, custom = NULL))
      if (dom_kind(dom, item) == "continuous")
        list(name = "viridis: viridis", min = "", max = "", reverse = FALSE, custom = NULL)
      else
        list(name = "Okabe-Ito",
             colors = palette_discrete(dom_levels(dom, item), NULL, "Okabe-Ito"))
    }

    # --- Per-item id helpers + observer registry ----------------------------
    key_of  <- function(dom, item) paste0(dom, "__", .pal_safe(item))
    pin_id  <- function(dom, item, i) paste0("pin_", key_of(dom, item), "_", i)
    iid     <- function(prefix, dom, item) paste0(prefix, "_", key_of(dom, item))

    registered  <- new.env(parent = emptyenv())
    obs_handles <- new.env(parent = emptyenv())
    unregister_item <- function(key) {
      h <- obs_handles[[key]]
      if (!is.null(h)) lapply(h, function(o) if (!is.null(o)) o$destroy())
      obs_handles[[key]] <- NULL
      registered[[key]] <- NULL
    }
    set_cfg <- function(dom, item, field, value) {
      p <- state$palette; p[[dom]][[item]][[field]] <- value; state$palette <- p
    }

    register_item <- function(dom, item) {
      key <- key_of(dom, item)
      if (!is.null(registered[[key]])) return(invisible())
      registered[[key]] <- TRUE
      cur <- function() state$palette[[dom]][[item]]

      remove_obs <- observeEvent(input[[iid("remove", dom, item)]], {
        p <- state$palette; p[[dom]][[item]] <- NULL; state$palette <- p
        unregister_item(key); bump()
      }, ignoreInit = TRUE)

      if (dom_kind(dom, item) == "continuous") {
        # Continuous: palette name + min/max anchors. Preview is the palette
        # gradient (range is resolved against real data at the consumer).
        name_obs <- observeEvent(input[[iid("cname", dom, item)]], {
          if (is.null(cur())) return()
          set_cfg(dom, item, "name", input[[iid("cname", dom, item)]] %||% "viridis: viridis")
        }, ignoreInit = TRUE)
        min_obs <- observeEvent(input[[iid("cmin", dom, item)]],
          if (!is.null(cur())) set_cfg(dom, item, "min", input[[iid("cmin", dom, item)]]),
          ignoreInit = TRUE)
        max_obs <- observeEvent(input[[iid("cmax", dom, item)]],
          if (!is.null(cur())) set_cfg(dom, item, "max", input[[iid("cmax", dom, item)]]),
          ignoreInit = TRUE)
        rev_obs <- observeEvent(input[[iid("crev", dom, item)]],
          if (!is.null(cur())) set_cfg(dom, item, "reverse", isTRUE(input[[iid("crev", dom, item)]])),
          ignoreInit = TRUE)
        # Custom-ramp colours: 5 pickers; store the filled (non-empty) ones.
        ccol_obs <- lapply(1:5, function(j) {
          observeEvent(input[[iid(paste0("ccol", j), dom, item)]], {
            if (is.null(cur())) return()
            vals <- vapply(1:5, function(k) {
              v <- input[[iid(paste0("ccol", k), dom, item)]]
              if (is.null(v) || !nzchar(v)) NA_character_ else norm_color(v)
            }, character(1))
            set_cfg(dom, item, "custom", unname(vals[!is.na(vals)]))
          }, ignoreInit = TRUE)
        })
        reset_obs <- observeEvent(input[[iid("creset", dom, item)]], {
          p <- state$palette
          p[[dom]][[item]] <- list(name = "viridis: viridis", min = "", max = "",
                                   reverse = FALSE, custom = NULL)
          state$palette <- p
          updateSelectInput(session, iid("cname", dom, item), selected = "viridis: viridis")
          updateTextInput(session, iid("cmin", dom, item), value = "")
          updateTextInput(session, iid("cmax", dom, item), value = "")
          updateCheckboxInput(session, iid("crev", dom, item), value = FALSE)
        }, ignoreInit = TRUE)
        obs_handles[[key]] <- c(ccol_obs,
          list(name_obs, min_obs, max_obs, rev_obs, reset_obs, remove_obs))
        output[[iid("cpreview", dom, item)]] <- renderUI({
          cfg <- cur(); req(cfg)
          .pal_gradient_bar(cfg$name %||% "viridis: viridis", cfg$custom, isTRUE(cfg$reverse))
        })
        return(invisible())
      }

      # Discrete: base palette + per-level hand-edits (flip to Custom).
      n0 <- length(dom_levels(dom, item))
      apply_named <- function(name) {
        lvls <- dom_levels(dom, item)
        cols <- palette_discrete(lvls, NULL, name)
        p <- state$palette
        p[[dom]][[item]]$name <- name; p[[dom]][[item]]$colors <- cols
        state$palette <- p
        for (i in seq_along(lvls)) update_picker(pin_id(dom, item, i), unname(cols[[lvls[i]]]))
      }
      name_obs <- observeEvent(input[[iid("name", dom, item)]], {
        if (is.null(cur())) return()
        name <- input[[iid("name", dom, item)]]; if (is.null(name)) return()
        if (identical(name, "Custom palette")) set_cfg(dom, item, "name", "Custom palette")
        else apply_named(name)
      }, ignoreInit = TRUE)
      pin_obs <- lapply(seq_len(n0), function(i) {
        observeEvent(input[[pin_id(dom, item, i)]], {
          if (is.null(cur())) return()
          lvls <- dom_levels(dom, item); if (i > length(lvls)) return()
          lev <- lvls[i]; hex <- norm_color(input[[pin_id(dom, item, i)]])
          if (is.na(hex)) return()
          p <- state$palette; cols <- p[[dom]][[item]]$colors
          if (identical(unname(cols[[lev]]), unname(hex))) return()
          cols[[lev]] <- unname(hex); p[[dom]][[item]]$colors <- cols
          was_named <- !identical(p[[dom]][[item]]$name, "Custom palette")
          if (was_named) p[[dom]][[item]]$name <- "Custom palette"
          state$palette <- p
          if (was_named)
            updateSelectInput(session, iid("name", dom, item), selected = "Custom palette")
        }, ignoreInit = TRUE)
      })
      reset_obs <- observeEvent(input[[iid("reset", dom, item)]], {
        apply_named("Okabe-Ito")
        updateSelectInput(session, iid("name", dom, item), selected = "Okabe-Ito")
      }, ignoreInit = TRUE)
      obs_handles[[key]] <- c(pin_obs, list(name_obs, reset_obs, remove_obs))
      invisible()
    }

    # High-cardinality guard: hide attributes above the hard cap from the add
    # control; warn (confirm) between the warn threshold and the cap.
    max_levels  <- function() getOption("ddsdashboard.palette_max_levels", 50L)
    warn_levels <- function() getOption("ddsdashboard.palette_warn_levels", 10L)
    just_added <- reactiveVal(NULL)        # so a rebuild opens only the new panel
    do_add <- function(dom, item) {
      p <- state$palette; p[[dom]][[item]] <- default_cfg(dom, item); state$palette <- p
      register_item(dom, item); just_added(list(dom = dom, item = item)); bump()
    }
    pending_add <- reactiveVal(NULL)
    observeEvent(input$confirm_add, {
      pa <- pending_add(); req(pa)
      do_add(pa$dom, pa$item); pending_add(NULL); removeModal()
    })

    # --- Wire each domain's add control, panels, and collapse ----------------
    for (dom in names(.pal_domains)) local({
      d <- dom
      output[[paste0("addui_", d)]] <- renderUI({
        if (dom_needs_data(d) && is.null(state$working))
          return(helpText(class = "text-muted", "Load a dataset to configure colours."))
        unconfigured <- setdiff(dom_items(d), names(state$palette[[d]]))
        over_cap <- Filter(function(it)
          dom_kind(d, it) == "discrete" && length(dom_levels(d, it)) > max_levels(),
          unconfigured)
        choices <- setdiff(unconfigured, over_cap)
        tagList(
          selectInput(ns(paste0("addsel_", d)), "Add colour mapping for",
                      choices = if (length(choices)) choices else character(0)),
          actionButton(ns(paste0("addbtn_", d)), "Add", icon = icon("plus"),
                       class = "btn-sm btn-primary",
                       disabled = if (!length(choices)) NA else NULL),
          if (!length(choices) && !length(over_cap) && length(dom_items(d)))
            helpText(class = "text-muted small mt-1", "All available items are configured."),
          if (length(over_cap))
            helpText(class = "text-muted small mt-1", .pal_hidden_note(over_cap, max_levels())))
      })
      observeEvent(input[[paste0("addbtn_", d)]], {
        item <- input[[paste0("addsel_", d)]]
        req(item, is.null(state$palette[[d]][[item]]))
        if (dom_kind(d, item) == "discrete") {
          n <- length(dom_levels(d, item))
          if (n > warn_levels()) {           # over-cap items are already filtered out
            pending_add(list(dom = d, item = item))
            showModal(modalDialog(title = "Many unique values",
              sprintf("'%s' has %d unique values; adding it creates %d colour pickers and may be slow.",
                      item, n, n),
              footer = tagList(modalButton("Cancel"),
                               actionButton(ns("confirm_add"), "Proceed", class = "btn-warning"))))
            return()
          }
        }
        do_add(d, item)
      })
      observeEvent(input[[paste0("collapse_", d)]],
                   bslib::accordion_panel_close(paste0("acc_", d), values = TRUE))
      output[[paste0("panels_", d)]] <- renderUI({
        struct()
        if (dom_needs_data(d) && is.null(state$working))
          return(tags$p(class = "text-muted p-2", "Load a dataset to configure colours."))
        cfg <- shiny::isolate(state$palette[[d]])
        items <- intersect(names(cfg), shiny::isolate(dom_items(d)))
        if (!length(items))
          return(tags$p(class = "text-muted p-2",
                        'No colour mappings yet. Use "Add colour mapping" in the sidebar.'))
        panels <- lapply(items, function(it) .palette_item_panel(ns, d, it, cfg[[it]],
                                                                 shiny::isolate(dom_kind(d, it)),
                                                                 shiny::isolate(dom_class(d, it)),
                                                                 shiny::isolate(dom_levels(d, it)),
                                                                 has_picker))
        # Preserve which panels are open across rebuilds (a rebuild otherwise
        # re-opens all): keep the currently-open set + the just-added item.
        open_now <- shiny::isolate(input[[paste0("acc_", d)]])
        ja <- shiny::isolate(just_added())
        open_set <- unique(c(open_now, if (!is.null(ja) && ja$dom == d) ja$item))
        do.call(bslib::accordion,
                c(list(id = ns(paste0("acc_", d)), multiple = TRUE,
                       open = if (length(open_set)) open_set else FALSE), panels))
      })
    })

    observeEvent(input$collapse_ref,
                 bslib::accordion_panel_close("ref_acc", values = TRUE))

    invisible(NULL)
  })
}

# Type (Discrete/Continuous) + data-class badges shown in each accordion title.
.pal_badges <- function(kind, klass) {
  tagList(
    tags$span(class = "badge rounded-pill text-bg-secondary ms-2 fw-normal",
              if (identical(kind, "continuous")) "Continuous" else "Discrete"),
    if (!is.null(klass) && nzchar(klass))
      tags$span(class = "badge rounded-pill text-bg-light ms-1 fw-normal", klass))
}

# Build the accordion panel for one item (discrete or continuous).
.palette_item_panel <- function(ns, dom, item, cfg, kind, klass, levels, has_picker) {
  key <- paste0(dom, "__", .pal_safe(item))
  id  <- function(prefix) ns(paste0(prefix, "_", key))
  remove_btn <- actionButton(id("remove"), "Remove mapping", icon = icon("trash"),
                             class = "btn-sm btn-outline-danger")

  body <- if (kind == "continuous") {
    name <- cfg$name %||% "viridis: viridis"
    # Up to 5 custom-ramp pickers, shown only when "Custom ramp" is selected
    # (client-side conditionalPanel keyed on the namespaced cname input). Slots
    # default to a valid viridis spread (never NULL -- a null initial colour makes
    # the colour picker throw and breaks every picker on the page); clear a slot to
    # drop it from the ramp.
    defaults5 <- norm_color(palette_colors("viridis: viridis", 5))
    ccol_picker <- function(j) {
      cv <- if (length(cfg$custom) >= j) cfg$custom[j] else defaults5[j]
      if (has_picker) {
        shinyWidgets::colorPickr(id(paste0("ccol", j)), label = NULL, selected = cv,
                                 update = "save", useAsButton = TRUE,
                                 interaction = list(input = TRUE, save = TRUE, clear = TRUE))
      } else {
        textInput(id(paste0("ccol", j)), label = NULL, value = cv, width = "90px")
      }
    }
    custom_ui <- conditionalPanel(
      condition = sprintf("input['%s'] == 'Custom ramp'", id("cname")),
      tags$div(class = "small fw-semibold mb-1", "Custom ramp colours (fill 2-5, low -> high)"),
      tags$div(class = "d-flex flex-wrap gap-2 mb-2 pal-level-row", lapply(1:5, ccol_picker)))
    tagList(
      selectInput(id("cname"), "Continuous palette",
                  choices = palette_continuous_choices()[
                    c("viridis", "Brewer: Sequential", "Brewer: Divergent", "Custom")],
                  selected = name),
      checkboxInput(id("crev"), "Reverse direction", value = isTRUE(cfg$reverse)),
      bslib::layout_columns(col_widths = c(6, 6),
        textInput(id("cmin"), "Min anchor", value = cfg$min %||% "",
                  placeholder = "e.g. 0  or  p5"),
        textInput(id("cmax"), "Max anchor", value = cfg$max %||% "",
                  placeholder = "e.g. 100  or  p95")),
      helpText(class = "text-muted small",
               "Leave blank to use the data min/max, or enter a number, or a percentile like p5 / p95."),
      custom_ui,
      uiOutput(id("cpreview")),
      tags$div(class = "d-flex gap-2 mt-2",
        bslib::tooltip(
          actionButton(id("creset"), "Reset", icon = icon("rotate-left"),
                       class = "btn-sm btn-outline-secondary"),
          "Reset to viridis with data-range anchors."),
        remove_btn))
  } else {
    name <- cfg$name %||% "Okabe-Ito"
    resolved <- palette_discrete(levels, cfg$colors, name, cfg$custom)
    # One compact row per level: the colour swatch/picker, then the level label to
    # its right (larger font). Rows stack vertically (kept for future factor
    # reordering).
    one_picker <- function(i) {
      lev <- levels[i]; val <- unname(resolved[[lev]]); pid <- ns(pin_id_str(key, i))
      picker <- if (has_picker) {
        shinyWidgets::colorPickr(pid, label = NULL, selected = val, update = "save",
                                 useAsButton = TRUE,
                                 interaction = list(input = TRUE, save = TRUE, clear = FALSE))
      } else {
        tagList(
          tags$span(style = sprintf(
            "display:inline-block;width:1.4rem;height:1.4rem;border-radius:3px;border:1px solid #888;background:%s;", val)),
          textInput(pid, label = NULL, value = val, width = "110px"))
      }
      tags$div(class = "d-flex align-items-center gap-2 mb-1 pal-level-row",
               picker, tags$span(class = "fw-medium", lev))
    }
    tagList(
      selectInput(id("name"), "Base palette", choices = palette_choices(), selected = name),
      if (length(levels) > .pal_many_levels)
        tags$div(class = "alert alert-warning py-1 px-2 small",
                 sprintf("%d levels - colours are interpolated; pinning each is tedious.",
                         length(levels))),
      tags$div(class = "small fw-semibold mb-1", "Levels"),
      tags$div(class = "mb-2", lapply(seq_along(levels), one_picker)),
      tags$div(class = "d-flex gap-2 mt-2",
        bslib::tooltip(
          actionButton(id("reset"), "Reset to palette", icon = icon("rotate-left"),
                       class = "btn-sm btn-outline-secondary"),
          "Reset to the default palette (Okabe-Ito)."),
        remove_btn))
  }
  bslib::accordion_panel(title = tags$span(item, .pal_badges(kind, klass)),
                         value = item, body)
}

# pin id string (without ns); mirrors the server's pin_id().
pin_id_str <- function(key, i) paste0("pin_", key, "_", i)
