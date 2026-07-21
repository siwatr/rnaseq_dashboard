# Reusable "continuous palette" control -- a base ramp (viridis / Brewer
# sequential / divergent / Custom ramp) + min/max anchors + reverse, plus an
# "Edit palette" copy-to-custom and a Reset. Extracted from the Palette page's
# per-item continuous panel so any plot that needs its own colour ramp (the P7
# Expression heatmap now; the Phase-8 shared heatmap controller later) reuses the
# exact same widgets + config shape instead of re-building them.
#
# Like mod_expr_value / mod_plot_subset, it operates in the *host* module's
# namespace: a page calls continuous_palette_ui(ns, suffix, default) in its
# sidebar and continuous_palette_server(input, output, session, suffix, default)
# once. `suffix` distinguishes instances. The server returns a reactive config
# `list(name, min, max, custom, reverse)` -- the same shape stored in
# `state$palette` -- which the caller feeds to palette_colorramp2() /
# palette_gradientn() to build a ComplexHeatmap col_fun or a ggplot scale.

# A sensible starting config (a white -> black custom ramp is never a good
# default, so seed a named ramp). Callers usually pass their own `default`.
continuous_palette_default <- function(name = "viridis: viridis", min = "", max = "",
                                       reverse = FALSE, custom = NULL) {
  list(name = name, min = min, max = max, reverse = isTRUE(reverse),
       custom = custom %||% c("#FFFFFF", "#000000"))
}

# The widget cluster. `default` seeds the initial values (the server keeps the
# authoritative config afterwards). Ids are prefixed by `<suffix>_`.
continuous_palette_ui <- function(ns, suffix, default = continuous_palette_default()) {
  id <- function(p) ns(paste0(suffix, "_", p))
  cfg <- default
  name <- cfg$name %||% "viridis: viridis"
  has_picker <- requireNamespace("shinyWidgets", quietly = TRUE)
  cur_custom <- if (length(cfg$custom) >= 2L) cfg$custom else c("#FFFFFF", "#000000")
  seed5 <- grDevices::colorRampPalette(cur_custom)(5)

  ccol_picker <- function(j) {
    if (has_picker)
      shinyWidgets::colorPickr(id(paste0("ccol", j)), label = NULL, selected = seed5[j],
                               update = "save", useAsButton = TRUE,
                               interaction = list(input = TRUE, save = TRUE, clear = FALSE))
    else
      textInput(id(paste0("ccol", j)), label = NULL, value = seed5[j], width = "90px")
  }
  custom_ui <- conditionalPanel(
    condition = sprintf("input['%s'] == 'Custom ramp'", id("cname")),
    selectInput(id("cnstops"), "Number of colours", choices = 2:5,
                selected = length(cur_custom)),
    tags$div(class = "small fw-semibold mb-1", "Ramp colours (low -> high)"),
    tags$div(class = "d-flex flex-wrap gap-2 mb-2 pal-level-row",
      lapply(1:5, function(j)
        conditionalPanel(
          condition = sprintf("parseInt(input['%s']) >= %d", id("cnstops"), j),
          ccol_picker(j)))))

  tagList(
    selectInput(id("cname"), "Colour palette",
                choices = palette_continuous_choices()[
                  c("viridis", "Brewer: Sequential", "Brewer: Divergent", "Custom")],
                selected = name),
    conditionalPanel(
      condition = sprintf("input['%s'] != 'Custom ramp'", id("cname")),
      bslib::tooltip(
        actionButton(id("cedit"), "Edit palette", icon = icon("sliders"),
                     class = "btn-sm btn-outline-secondary mb-2"),
        "Copy this palette's colours into an editable 5-stop Custom ramp.")),
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
    bslib::tooltip(
      actionButton(id("creset"), "Reset ramp", icon = icon("rotate-left"),
                   class = "btn-sm btn-outline-secondary mt-2"),
      "Clear the anchors, reverse, and any custom-ramp colours."))
}

# Wire the observers + preview and return a reactive config. `default` is the
# initial config (must match the UI's `default`). The returned reactive yields
# `list(name, min, max, custom, reverse)`.
continuous_palette_server <- function(input, output, session, suffix,
                                      default = continuous_palette_default()) {
  iid <- function(p) paste0(suffix, "_", p)          # module-relative input id
  cfg <- reactiveVal(default)
  set_field <- function(field, value) {
    p <- cfg(); p[[field]] <- value; cfg(p)
  }
  update_picker <- function(pid, value) {
    if (requireNamespace("shinyWidgets", quietly = TRUE))
      shinyWidgets::updateColorPickr(session, pid, value = value)
    else updateTextInput(session, pid, value = value)
  }
  cur_custom <- function() {
    cc <- cfg()$custom
    if (length(cc) >= 2L) cc else c("#FFFFFF", "#000000")
  }
  n_stops <- function() {
    n <- suppressWarnings(as.integer(input[[iid("cnstops")]]))
    if (length(n) != 1L || is.na(n)) length(cur_custom()) else n
  }

  observeEvent(input[[iid("cname")]],
    set_field("name", input[[iid("cname")]] %||% "viridis: viridis"), ignoreInit = TRUE)
  observeEvent(input[[iid("cmin")]], set_field("min", input[[iid("cmin")]]), ignoreInit = TRUE)
  observeEvent(input[[iid("cmax")]], set_field("max", input[[iid("cmax")]]), ignoreInit = TRUE)
  observeEvent(input[[iid("crev")]],
    set_field("reverse", isTRUE(input[[iid("crev")]])), ignoreInit = TRUE)

  # Stop count resamples the current ramp (preserving endpoints/shape) + repaints.
  observeEvent(input[[iid("cnstops")]], {
    newN <- suppressWarnings(as.integer(input[[iid("cnstops")]]))
    if (is.na(newN) || newN < 2L || newN > 5L) return()
    if (length(cfg()$custom) == newN) return()
    new <- grDevices::colorRampPalette(cur_custom())(newN)
    set_field("custom", new)
    for (j in seq_len(newN)) update_picker(iid(paste0("ccol", j)), new[j])
  }, ignoreInit = TRUE)

  # Each visible picker (1..N) contributes one anchor colour.
  lapply(1:5, function(j) {
    observeEvent(input[[iid(paste0("ccol", j))]], {
      n <- n_stops(); if (j > n) return()
      vals <- vapply(seq_len(n), function(k) {
        v <- norm_color(input[[iid(paste0("ccol", k))]])
        if (is.na(v)) "#000000" else v
      }, character(1))
      set_field("custom", unname(vals))
    }, ignoreInit = TRUE)
  })

  # Edit palette: copy a named ramp into an editable 5-stop Custom ramp, current
  # reverse baked into the order then zeroed (on-screen gradient unchanged).
  observeEvent(input[[iid("cedit")]], {
    nm <- cfg()$name %||% "viridis: viridis"
    if (identical(nm, "Custom ramp")) return()
    stops <- .continuous_stops(nm, NULL, n = 5L, reverse = isTRUE(cfg()$reverse))
    p <- cfg(); p$name <- "Custom ramp"; p$custom <- unname(stops); p$reverse <- FALSE
    cfg(p)
    updateSelectInput(session, iid("cname"), selected = "Custom ramp")
    updateSelectInput(session, iid("cnstops"), selected = 5L)
    updateCheckboxInput(session, iid("crev"), value = FALSE)
    for (j in seq_len(5L)) update_picker(iid(paste0("ccol", j)), stops[j])
  }, ignoreInit = TRUE)

  # Reset = keep the chosen palette, clear anchors/reverse (Custom -> white/black).
  observeEvent(input[[iid("creset")]], {
    nm <- cfg()$name %||% "viridis: viridis"
    new <- list(name = nm, min = "", max = "", reverse = FALSE,
                custom = if (identical(nm, "Custom ramp")) c("#FFFFFF", "#000000")
                         else cfg()$custom)
    cfg(new)
    updateSelectInput(session, iid("cname"), selected = new$name)
    updateTextInput(session, iid("cmin"), value = "")
    updateTextInput(session, iid("cmax"), value = "")
    updateCheckboxInput(session, iid("crev"), value = FALSE)
    cc <- if (length(new$custom) >= 2L) new$custom else c("#FFFFFF", "#000000")
    updateSelectInput(session, iid("cnstops"), selected = length(cc))
    for (j in seq_along(cc)) update_picker(iid(paste0("ccol", j)), cc[j])
  }, ignoreInit = TRUE)

  output[[iid("cpreview")]] <- renderUI({
    cf <- cfg(); req(cf)
    .pal_gradient_bar(cf$name, cf$custom, reverse = isTRUE(cf$reverse))
  })

  reactive(cfg())
}
