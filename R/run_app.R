#' Launch the dds dashboard
#'
#' @param ... Passed to [shiny::shinyApp()].
#' @return A Shiny app object (invisibly when run interactively).
#' @export
run_app <- function(...) {
  shiny::shinyApp(ui = app_ui(), server = app_server, ...)
}

#' Dashboard UI
#'
#' The multi-page shell. Each page is a Shiny module (see `R/mod_*.R`).
#' @return A bslib `page_navbar` UI definition.
#' @export
app_ui <- function() {
  bslib::page_navbar(
    title = "dds dashboard",
    theme = bslib::bs_theme(version = 5),
    bslib::nav_panel("Input",     mod_input_ui("input")),
    bslib::nav_panel("QC",        mod_qc_ui("qc")),
    bslib::nav_panel("Process",   mod_process_ui("process")),
    bslib::nav_panel("DimReduc",  mod_dimreduc_ui("dimreduc")),
    bslib::nav_panel("DE",        mod_de_ui("de")),
    bslib::nav_panel("Heatmap",   mod_heatmap_ui("heatmap")),
    bslib::nav_panel("Export",    mod_export_ui("export"))
  )
}

#' Dashboard server
#'
#' Owns the canonical `dds`. Each module receives the current `dds` as a
#' reactive and returns a (possibly updated) `dds`; the top-level server
#' rewires the canonical copy. See the `shiny-module` skill for the contract.
#' @param input,output,session Shiny server arguments.
#' @return Invisible `NULL`.
#' @export
app_server <- function(input, output, session) {
  # Modules that can mutate the object return it; consumers just read it.
  dds <- mod_input_server("input")
  dds <- mod_qc_server("qc", dds = dds)
  dds <- mod_process_server("process", dds = dds)
  mod_dimreduc_server("dimreduc", dds = dds)
  dds <- mod_de_server("de", dds = dds)
  mod_heatmap_server("heatmap", dds = dds)
  mod_export_server("export", dds = dds)
  invisible(NULL)
}
