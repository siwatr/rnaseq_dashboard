# Page 6: Expression heatmap
# ComplexHeatmap over a gene set (default = DEGs). Expression default
# log10(TPM + 0.01), fall back to log10(CPM + 0.01). Top annotation from the
# design column(s). A sub-panel controls palette / range / annotation colors.

mod_heatmap_ui <- function(id) {
  ns <- NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      title = "Heatmap",
      textAreaInput(ns("genes"), "Genes of interest (one per line)"),
      actionButton(ns("render"), "Render", class = "btn-primary")
      # TODO: palette / range / annotation controls
    ),
    bslib::card(
      bslib::card_header("Expression heatmap"),
      plotOutput(ns("heatmap"))
    )
  )
}

#' @param dds reactive() yielding the current DESeqDataSet.
#' @return Invisible NULL (consumes dds, does not mutate it).
mod_heatmap_server <- function(id, dds) {
  moduleServer(id, function(input, output, session) {
    # TODO: build ComplexHeatmap behind the render button.
    invisible(NULL)
  })
}
