# Page 1, "Feature info" tab: view + annotate + tag feature metadata. OrgDb
# annotation populates gene name / chromosome / description; feature_class
# tagging marks endogenous / spike_in / exogenous. Both apply via state_mutate.
# GTF annotation and a fuller rowData editor come in later PRs.

mod_feature_ui <- function(id) {
  ns <- NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      title = "Annotate & tag", width = 340,
      tags$strong("Annotate from OrgDb"),
      selectInput(ns("organism"), "Organism",
                  c("Mouse (org.Mm.eg.db)" = "mouse", "Human (org.Hs.eg.db)" = "human")),
      selectInput(ns("id_type"), "Feature id type",
                  c("Auto-detect" = "auto", "Ensembl" = "ensembl",
                    "Entrez" = "entrez", "Symbol" = "symbol")),
      helpText(textOutput(ns("detected_id"), inline = TRUE)),
      actionButton(ns("annotate"), "Annotate from OrgDb", class = "btn-primary"),
      hr(),
      tags$strong("Tag feature class"),
      selectizeInput(ns("feat"), "Find features", choices = NULL, multiple = TRUE,
                     options = list(placeholder = "gene name or id")),
      radioButtons(ns("fclass"), "Set class to",
                   c("endogenous", "spike_in", "exogenous")),
      actionButton(ns("tag"), "Apply tag")
    ),
    bslib::card(
      bslib::card_header("Feature metadata"),
      textOutput(ns("coverage")),
      DT::DTOutput(ns("rowdata"))
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

    output$detected_id <- renderText({
      req(state$working)
      sprintf("Detected: %s", detect_id_type(rownames(state$working)))
    })

    # Server-side feature search, labelled by <feature_type>_name when present.
    observeEvent(state$working, {
      req(state$working)
      ids <- rownames(state$working)
      name_col <- paste0(state_meta(state)$feature_type, "_name")
      rd <- SummarizedExperiment::rowData(state$working)
      labels <- if (name_col %in% colnames(rd)) paste0(rd[[name_col]], " (", ids, ")") else ids
      updateSelectizeInput(session, "feat",
                           choices = stats::setNames(ids, labels), server = TRUE)
    })

    observeEvent(input$annotate, {
      req(state$working)
      organism <- input$organism
      id_type <- if (identical(input$id_type, "auto")) NULL else input$id_type
      ft <- state_meta(state)$feature_type
      ok <- tryCatch({
        state_mutate(state,
          function(d) annotate_with_orgdb(d, organism = organism, id_type = id_type, feature_type = ft),
          action = list(action = "annotate_orgdb", organism = organism))
        TRUE
      }, error = function(e) { showNotification(conditionMessage(e), type = "error", duration = NULL); FALSE })
      req(ok)
      # Re-detect the feature unit now that a name column may exist.
      newft <- detect_feature_type(state$working)
      state$meta <- utils::modifyList(state$meta, list(feature_type = newft$feature_type))
      cov <- annotation_coverage(state$working, paste0(newft$feature_type, "_name"))
      showNotification(sprintf("Annotated %d of %d endogenous features.", cov$matched, cov$total),
                       type = "message")
    })

    observeEvent(input$tag, {
      req(state$working, input$feat)
      cls <- input$fclass; ids <- input$feat
      state_mutate(state, function(d) set_feature_class(d, ids, cls),
                   action = list(action = "set_feature_class", class = cls, n = length(ids)))
      showNotification(sprintf("Tagged %d feature(s) as %s.", length(ids), cls), type = "message")
    })

    output$coverage <- renderText({
      req(state$working)
      name_col <- paste0(state_meta(state)$feature_type, "_name")
      cov <- annotation_coverage(state$working, name_col)
      if (cov$matched == 0L) "Not annotated yet."
      else sprintf("Annotated %d of %d endogenous features%s.", cov$matched, cov$total,
                   if (cov$matched < cov$total) " (a GTF upload can fill gaps; coming soon)" else "")
    })

    output$rowdata <- DT::renderDT({
      req(state$working)
      rd <- as.data.frame(SummarizedExperiment::rowData(state$working))
      keep <- intersect(c(paste0(state_meta(state)$feature_type, "_name"),
                          "chromosome", "description", "feature_class", "feature_length"),
                        colnames(rd))
      DT::datatable(cbind(id = rownames(state$working), rd[, keep, drop = FALSE]),
                    rownames = FALSE, selection = "none",
                    options = list(pageLength = 10, scrollX = TRUE))
    }, server = TRUE)

    output$summary <- renderPrint({
      req(state$working)
      print(table(SummarizedExperiment::rowData(state$working)$feature_class))
    })

    invisible(NULL)
  })
}
