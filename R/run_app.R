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
    bslib::nav_panel("Export",    mod_export_ui("export")),
    bslib::nav_spacer(),
    bslib::nav_item(mod_statusbar_ui("statusbar"))
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
  state <- new_app_state()
  mod_input_server("input", state)
  mod_qc_server("qc", state)
  mod_process_server("process", state)
  mod_dimreduc_server("dimreduc", state)
  mod_de_server("de", state)
  mod_heatmap_server("heatmap", state)
  mod_export_server("export", state)
  mod_statusbar_server("statusbar", state)
  invisible(NULL)
}
