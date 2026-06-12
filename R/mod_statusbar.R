# Persistent dataset-status bar shown in the navbar on every page. Reads the
# app state via state_meta() and renders compact badges: data type, dimensions,
# assays, design, edit count, and a warning when single-cell was coerced per-cell.

mod_statusbar_ui <- function(id) {
  ns <- NS(id)
  uiOutput(ns("status"), inline = TRUE)
}

.badge <- function(text, class = "text-bg-secondary") {
  tags$span(class = paste("badge rounded-pill", class), text)
}

mod_statusbar_server <- function(id, state) {
  moduleServer(id, function(input, output, session) {
    output$status <- renderUI({
      m <- state_meta(state)
      if (!isTRUE(m$loaded)) return(.badge("no dataset loaded", "text-bg-light"))
      badges <- list(
        .badge(m$data_type, "text-bg-info"),
        .badge(sprintf("%s features x %s samples", m$n_features, m$n_samples)),
        .badge(paste("assays:", paste(m$assays, collapse = ", "))),
        .badge(paste("design:", m$design))
      )
      if (m$n_edits > 0L) badges <- c(badges, list(.badge(sprintf("%d edits", m$n_edits))))
      if (isTRUE(m$sce_per_cell)) {
        badges <- c(badges, list(.badge("per-cell: stats unreliable", "text-bg-warning")))
      }
      tags$span(class = "d-flex gap-1 align-items-center", badges)
    })
  })
}
