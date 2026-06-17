# Page 2: Quality control & filtering
# scater QC metric plots, variance-stabilization plot, sample correlation
# heatmap; suggest + apply sample/feature filtering. Removing samples or
# features must invalidate downstream assays and the DESeq fit (rnaseq-bioc).

mod_qc_ui <- function(id) {
  ns <- NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      title = "QC & filtering",
      actionButton(ns("render"), "Render QC", class = "btn-primary")
      # TODO: metric selector, group-by column, filter thresholds
    ),
    bslib::card(
      bslib::card_header(tags$h3("QC", class = "fs-6 mb-0")),
      plotOutput(ns("qc_plot"))
    )
  )
}

#' @param state the shared app-state object (see [new_app_state()]).
#' @return Invisible NULL.
mod_qc_server <- function(id, state) {
  moduleServer(id, function(input, output, session) {
    # TODO: compute QC metrics (per-sample/per-cell) behind the render button;
    # filter via state_mutate(); cache via state_derive().
    invisible(NULL)
  })
}
