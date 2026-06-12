# Page 3: Processing
# Add normalized assays (TPM, FPKM, CPM) to the dds. TPM/FPKM require an
# effective feature length; degrade to CPM when absent (see rnaseq-bioc and
# R/utils_normalization.R).

mod_process_ui <- function(id) {
  ns <- NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      title = "Add assays",
      checkboxGroupInput(ns("assays"), "Assays to compute",
                         choices = c("CPM", "TPM", "FPKM")),
      actionButton(ns("apply"), "Add assays", class = "btn-primary")
    ),
    bslib::card(
      bslib::card_header("Assays present"),
      verbatimTextOutput(ns("assay_list"))
    )
  )
}

#' @param dds reactive() yielding the current DESeqDataSet.
#' @return reactive() yielding the DESeqDataSet with added assays.
mod_process_server <- function(id, dds) {
  moduleServer(id, function(input, output, session) {
    # TODO: on $apply, compute selected assays via cpm()/tpm()/fpkm() and
    # attach them, then return the updated object.
    reactive(dds())
  })
}
