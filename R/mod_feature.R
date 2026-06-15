# Page 1, "Feature info" tab: view/tag feature metadata. For now this is
# feature_class tagging (endogenous / spike_in / exogenous) with a server-side
# feature search, plus a class summary. Annotation (OrgDb/GTF) and a fuller
# rowData view arrive in a later PR. Backed by set_feature_class() in
# R/metadata_helpers.R; tagging applies immediately via state_mutate().

mod_feature_ui <- function(id) {
  ns <- NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      title = "Tag feature class", width = 340,
      selectizeInput(ns("feat"), "Find features", choices = NULL, multiple = TRUE,
                     options = list(placeholder = "gene name or id")),
      radioButtons(ns("fclass"), "Set class to",
                   c("endogenous", "spike_in", "exogenous")),
      actionButton(ns("tag"), "Apply tag", class = "btn-primary")
    ),
    bslib::card(
      bslib::card_header("Feature classes"),
      verbatimTextOutput(ns("summary"))
    )
  )
}

#' @param state the shared app-state object (see [new_app_state()]).
#' @return Invisible NULL.
mod_feature_server <- function(id, state) {
  moduleServer(id, function(input, output, session) {

    # Server-side search, labelled by <feature_type>_name when present.
    observeEvent(state$working, {
      req(state$working)
      ids <- rownames(state$working)
      name_col <- paste0(state_meta(state)$feature_type, "_name")
      rd <- SummarizedExperiment::rowData(state$working)
      labels <- if (name_col %in% colnames(rd)) paste0(rd[[name_col]], " (", ids, ")") else ids
      updateSelectizeInput(session, "feat",
                           choices = stats::setNames(ids, labels), server = TRUE)
    })

    observeEvent(input$tag, {
      req(state$working, input$feat)
      cls <- input$fclass
      ids <- input$feat
      state_mutate(state, function(d) set_feature_class(d, ids, cls),
                   action = list(action = "set_feature_class", class = cls, n = length(ids)))
      showNotification(sprintf("Tagged %d feature(s) as %s.", length(ids), cls), type = "message")
    })

    output$summary <- renderPrint({
      req(state$working)
      print(table(SummarizedExperiment::rowData(state$working)$feature_class))
    })

    invisible(NULL)
  })
}
