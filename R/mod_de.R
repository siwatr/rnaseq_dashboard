# Page 5: Differential expression
# Design the formula, run DESeq2::DESeq() (rerun after any sample/feature
# change), set thresholds, add sig/DEG columns, build MA / volcano / direct
# comparison plots. See rnaseq-bioc for the sig/DEG rules and axis mappings.

mod_de_ui <- function(id) {
  ns <- NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      title = "Differential expression",
      numericInput(ns("padj"), "padj threshold", value = 0.05,
                   min = 0, max = 1, step = 0.01),
      numericInput(ns("lfc"), "abs(log2FC) threshold", value = log2(2),
                   min = 0, step = 0.1),
      selectInput(ns("plot_type"), "Plot",
                  c("MA", "Volcano", "Direct comparison")),
      actionButton(ns("run"), "Run / render", class = "btn-primary")
    ),
    bslib::card(
      bslib::card_header("DE plot"),
      plotOutput(ns("plot"))
    )
  )
}

#' @param dds reactive() yielding the current DESeqDataSet.
#' @return reactive() yielding the DESeqDataSet (with DESeq fit attached).
mod_de_server <- function(id, dds) {
  moduleServer(id, function(input, output, session) {
    # TODO: run DESeq() behind the button; cache fit; build the chosen plot
    # with the sig/DEG columns and axis clamping (triangle for clamped points).
    reactive(dds())
  })
}
