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
      bslib::card_header("QC"),
      plotOutput(ns("qc_plot"))
    )
  )
}

#' @param dds reactive() yielding the current DESeqDataSet.
#' @return reactive() yielding the (possibly filtered) DESeqDataSet.
mod_qc_server <- function(id, dds) {
  moduleServer(id, function(input, output, session) {
    # TODO: compute scater per-cell QC metrics behind the render button.
    reactive(dds())
  })
}
