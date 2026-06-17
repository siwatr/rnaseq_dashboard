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
      bslib::card_header(tags$h3("Expression heatmap", class = "fs-6 mb-0")),
      plotOutput(ns("heatmap"))
    )
  )
}

#' @param state the shared app-state object (see [new_app_state()]).
#' @return Invisible NULL (consumes dds, does not mutate it).
mod_heatmap_server <- function(id, state) {
  moduleServer(id, function(input, output, session) {
    # TODO: build ComplexHeatmap behind the render button.
    invisible(NULL)
  })
}
