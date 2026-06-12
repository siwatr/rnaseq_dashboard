# Page 4: Dimensionality reduction & gene expression
# PCA / t-SNE / UMAP from top-variable genes; up to 4 panels (layouts 1/2/4),
# color/shape by metadata or gene expression. Cache VST/PCA results keyed on
# dds state + params; invalidate on sample/feature change (see CLAUDE.md).

mod_dimreduc_ui <- function(id) {
  ns <- NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      title = "Dim. reduction",
      selectInput(ns("method"), "Method", c("PCA", "t-SNE", "UMAP")),
      numericInput(ns("n_top"), "Top variable genes", value = 500, min = 2),
      actionButton(ns("render"), "Render", class = "btn-primary")
    ),
    bslib::card(
      bslib::card_header("Reduced dimensions"),
      plotOutput(ns("plot"))
    )
  )
}

#' @param dds reactive() yielding the current DESeqDataSet.
#' @return Invisible NULL (consumes dds, does not mutate it).
mod_dimreduc_server <- function(id, dds) {
  moduleServer(id, function(input, output, session) {
    # TODO: compute embedding behind the render button; cache it.
    invisible(NULL)
  })
}
