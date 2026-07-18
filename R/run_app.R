# Default ceiling for browser uploads (.rds objects and tables). Shiny's own
# default is only 5 MB, far too small for real DESeqDataSet .rds files.
.default_max_upload_mb <- 500L

# The app theme: Bootstrap 5 with a light custom palette + a system-font stack
# (no runtime web-font fetch, so it works offline), plus the compiled custom
# rules (denser DT tables + dark-mode fixes) in inst/www/custom.scss. The numeric
# starting values here are meant to be tuned interactively with
# `bslib::run_with_themer(run_app())`; dark mode is handled by Bootstrap's
# [data-bs-theme] and the scoped rules in custom.scss.
.app_theme <- function() {
  theme <- bslib::bs_theme(
    version = 5,
    # primary = "#2C7FB8",
    bg = "#fffefb", fg = "#111111", 
    primary = "#14D1AF", secondary = "#41ebcc", info = "#C4A381",
    success = "#14D1AF", warning = "#eaaa3c", danger = "#e34125",
    dark = "#3a353f", light = "#FAFAF8",
    "border-radius" = "0.5rem", 
    heading_font = bslib::font_collection(
      sass::font_google("Unbounded"), 
      sass::font_google("Zain"), 
      # sass::font_google("Julius Sans One"), 
      # sass::font_google("Major Mono Display"), 
      sass::font_google("Roboto"), 
      sass::font_google("Noto Sans"), 
      "system-ui", "-apple-system"
    ),
    base_font = bslib::font_collection(
      sass::font_google("SUSE"), 
      # sass::font_google("Vend Sans"),
      sass::font_google("Nunito"), 
      sass::font_google("Noto Sans"), 
      sass::font_google("Roboto"),
      "Helvetica Neue", "Arial",
      "system-ui", "-apple-system", "Segoe UI", "sans-serif"
    ),
    code_font = bslib::font_collection(
      sass::font_google("SUSE Mono"), 
      "Menlo", "Monaco", "Consolas", "monospace"),
    preset = "bootstrap"
  )
  scss <- system.file("www", "custom.scss", package = "ddsdashboard")
  if (nzchar(scss)) theme <- bslib::bs_add_rules(theme, readLines(scss))
  theme
}

# Candidate ComplexHeatmap palettes (low -> high) shown as gradient swatches on
# the Themer "Heatmap" sub-tab, so a palette can be chosen for the QC sample-
# correlation heatmap (and later the expression heatmap). Pure reference data.
.heatmap_palettes <- list(
  "Blues (sequential)"        = c("#fff7fb", "#74a9cf", "#023858"),
  "Viridis"                   = c("#440154", "#21908C", "#FDE725"),
  "Magma"                     = c("#000004", "#B63679", "#FCFDBF"),
  "Blue-White-Red (diverging)" = c("#2166AC", "#F7F7F7", "#B2182B")
)

# A static component gallery shown as the last navbar tab when
# run_app(themer_mode = TRUE). Pair with bslib::run_with_themer() to watch theme
# / colour choices apply to headings, body text, the Bootstrap role colours
# (primary/secondary/success/info/warning/danger/light/dark) across buttons /
# badges / alerts, and candidate heatmap palettes. Static UI only -- no server.
.themer_ui <- function() {
  roles <- c("primary", "secondary", "success", "info",
             "warning", "danger", "light", "dark")
  each <- function(fn) lapply(roles, fn)
  cap <- function(x) paste0(toupper(substring(x, 1, 1)), substring(x, 2))
  swatch <- function(name, cols) {
    grad <- paste(cols, collapse = ", ")
    tags$div(class = "mb-3",
      tags$div(class = "small fw-semibold mb-1", name),
      tags$div(style = sprintf(
        "height:28px;border-radius:0.25rem;background:linear-gradient(to right, %s);",
        grad)))
  }
  bslib::nav_panel(
    tags$h2("Themer", class = "fs-6 mb-0"),
    bslib::navset_pill(
      bslib::nav_panel(
        "Components",
        bslib::layout_columns(
          col_widths = c(6, 6, 6, 6),
          bslib::card(
            bslib::card_header("Typography"),
            tags$h1("Heading 1"), tags$h2("Heading 2"), tags$h3("Heading 3"),
            tags$h4("Heading 4"), tags$h5("Heading 5"), tags$h6("Heading 6"),
            tags$p("Body text - the quick brown fox jumps over the lazy dog."),
            tags$p(class = "mb-0",
                   tags$a(href = "#", "a link"), " | ", tags$code("inline code"),
                   " | ", tags$small(class = "text-muted", "small muted text"))
          ),
          bslib::card(
            bslib::card_header("Buttons"),
            div(class = "d-flex flex-wrap gap-2 mb-2",
                each(function(r) tags$button(type = "button",
                                             class = paste0("btn btn-", r), cap(r)))),
            div(class = "d-flex flex-wrap gap-2",
                each(function(r) tags$button(type = "button",
                                             class = paste0("btn btn-outline-", r), cap(r))))
          ),
          bslib::card(
            bslib::card_header("Badges"),
            div(class = "d-flex flex-wrap gap-2",
                each(function(r) tags$span(class = paste("badge rounded-pill",
                                                         paste0("text-bg-", r)), cap(r))))
          ),
          bslib::card(
            bslib::card_header("Alerts"),
            each(function(r) div(class = paste0("alert alert-", r), role = "alert",
                                 paste(cap(r), "alert")))
          )
        )
      ),
      bslib::nav_panel(
        "Heatmap",
        bslib::card(
          bslib::card_header("Candidate heatmap palettes"),
          tags$p(class = "text-muted small",
                 "Reference swatches for the ComplexHeatmap colormap (QC sample correlation, later the expression heatmap)."),
          lapply(names(.heatmap_palettes),
                 function(nm) swatch(nm, .heatmap_palettes[[nm]]))
        )
      )
    )
  )
}

#' Launch the dds dashboard
#'
#' @param max_upload_mb Maximum size (in MB) for a single browser upload — covers
#'   the `.rds` object and the counts/sample-sheet tables. Sets the global
#'   `shiny.maxRequestSize` option. Raise it for very large datasets.
#' @param themer_mode If `TRUE`, add a "Themer" tab (a component gallery for
#'   tuning the theme); pair with [bslib::run_with_themer()]. Default `FALSE`.
#' @param ... Passed to [shiny::shinyApp()].
#' @return A Shiny app object (invisibly when run interactively).
#' @export
run_app <- function(max_upload_mb = .default_max_upload_mb, themer_mode = FALSE, ...) {
  options(shiny.maxRequestSize = max_upload_mb * 1024^2)
  shiny::shinyApp(ui = app_ui(themer_mode = themer_mode), server = app_server, ...)
}

#' Dashboard UI
#'
#' The multi-page shell. Each page is a Shiny module (see `R/mod_*.R`).
#' @param themer_mode If `TRUE`, append a "Themer" component-gallery tab for
#'   theme tuning (see [run_app()]). Default `FALSE`.
#' @return A bslib `page_navbar` UI definition.
#' @export
app_ui <- function(themer_mode = FALSE) {
  bslib::page_navbar(
    # Real heading element so it reads as the app title (distinct from the tab
    # labels) and picks up the theme's heading font; sized down to fit the navbar.
    title = tags$h1("DDS Dashboard", class = "fs-3 fw-bold mb-0 pb-0 pe-3",
                    style = "color:#8b58db;"),
    theme = .app_theme(),
    bslib::nav_panel(tags$h2("Dataset", class = "fs-6 mb-0"),
      bslib::navset_card_tab(
        bslib::nav_panel(tags$h3("Import", class = "fs-6"),        mod_input_ui("input")),
        bslib::nav_panel(tags$h3("Sample", class = "fs-6"),  mod_metadata_ui("metadata")),
        bslib::nav_panel(tags$h3("Feature", class = "fs-6"), mod_feature_ui("feature")),
        bslib::nav_panel(tags$h3("Assay", class = "fs-6"),        mod_assay_ui("assay")),
        bslib::nav_panel(tags$h3("Design", class = "fs-6"),
          bslib::card(bslib::card_body(
            tags$p(class = "text-muted",
                   "Set the model design and reference (control) levels early. This only sets the stage - running differential expression happens on the DE page (the design there stays in sync with this)."),
            tags$div(style = "max-width: 480px;",
                     mod_design_builder_ui("design_input")))))
      )
    ),
    bslib::nav_panel(tags$h2("QC", class = "fs-6 mb-0"),        mod_qc_ui("qc")),
    bslib::nav_panel(tags$h2("DimReduc", class = "fs-6 mb-0"),  mod_dimreduc_ui("dimreduc")),
    bslib::nav_panel(tags$h2("DE", class = "fs-6 mb-0"),        mod_de_ui("de")),
    bslib::nav_panel(tags$h2("Gene Sets", class = "fs-6 mb-0"), mod_geneset_ui("geneset")),
    bslib::nav_panel(tags$h2("Expression", class = "fs-6 mb-0"), mod_expression_ui("expression")),
    bslib::nav_panel(tags$h2("Palette", class = "fs-6 mb-0"),   mod_palette_ui("palette")),
    bslib::nav_panel(tags$h2("Export", class = "fs-6 mb-0"),    mod_export_ui("export")),
    if (isTRUE(themer_mode)) .themer_ui(),
    bslib::nav_spacer(),
    # Status bar + global view controls share ONE nav_item so they sit together
    # (separate nav_items get navbar spacing between them). The dark-mode switch
    # keeps its global, un-namespaced id so app_server reads input$dark_mode.
    # No `mode` arg => follows the OS prefers-color-scheme, persisted per browser.
    bslib::nav_item(
      tags$div(class = "d-flex align-items-center gap-2",
               mod_statusbar_ui("statusbar"),
               bslib::input_dark_mode(id = "dark_mode")))
  )
}

#' Dashboard server
#'
#' Creates the shared app-state object and threads it into every module. Modules
#' read it via `state_dds()`, edit via `state_mutate()`, and cache via
#' `state_derive()`. See the `shiny-module` skill for the contract.
#' @param input,output,session Shiny server arguments.
#' @return Invisible `NULL`.
#' @export
app_server <- function(input, output, session) {
  # Ensure a sane upload ceiling even if launched directly (not via run_app()).
  if (is.null(getOption("shiny.maxRequestSize"))) {
    options(shiny.maxRequestSize = .default_max_upload_mb * 1024^2)
  }
  # Make plots (ggplot2/base/lattice) follow the active bslib theme, including
  # the dark-mode switch. Takes effect for the Phase-3 QC plots onward.
  thematic::thematic_shiny(font = "auto")
  state <- new_app_state()
  mod_input_server("input", state)
  mod_metadata_server("metadata", state)
  mod_feature_server("feature", state)
  mod_assay_server("assay", state)
  mod_design_builder_server("design_input", state)   # Dataset > Design tab (synced with the DE page)
  mod_qc_server("qc", state, dark_mode = reactive(identical(input$dark_mode, "dark")))
  mod_dimreduc_server("dimreduc", state, dark_mode = reactive(identical(input$dark_mode, "dark")))
  mod_de_server("de", state, dark_mode = reactive(identical(input$dark_mode, "dark")))
  mod_geneset_server("geneset", state, dark_mode = reactive(identical(input$dark_mode, "dark")))
  mod_expression_server("expression", state, dark_mode = reactive(identical(input$dark_mode, "dark")))
  mod_palette_server("palette", state)
  mod_export_server("export", state)
  mod_statusbar_server("statusbar", state)
  invisible(NULL)
}
