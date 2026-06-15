# Page 1, "Feature info" tab: OrgDb annotation (applies immediately) plus the
# shared metadata draft editor configured for rowData (features). The editable
# table replaces the old manual feature_class tagging: feature_class is edited
# directly (validated to its allowed values) or set in bulk on the filtered
# rows. Editing logic lives in R/mod_meta_editor.R.

.feature_editor_opts <- list(
  slot = "rowData", title = "Feature information", row_noun = "feature",
  allow_merge = FALSE, allow_row_rename = FALSE, bulk_class = TRUE
)

mod_feature_ui <- function(id) {
  ns <- NS(id)
  annotation_ui <- tagList(
    tags$strong("Annotate from OrgDb"),
    selectInput(ns("organism"), "Organism",
                c("Mouse (org.Mm.eg.db)" = "mouse", "Human (org.Hs.eg.db)" = "human")),
    selectInput(ns("id_type"), "Feature id type",
                c("Auto-detect" = "auto", "Ensembl" = "ensembl",
                  "Entrez" = "entrez", "Symbol" = "symbol")),
    helpText(textOutput(ns("detected_id"), inline = TRUE)),
    actionButton(ns("annotate"), "Annotate from OrgDb", class = "btn-primary"),
    hr()
  )
  meta_editor_ui(ns("editor"), .feature_editor_opts,
                 extra_sidebar = annotation_ui,
                 extra_main = div(class = "mb-2", textOutput(ns("coverage"))))
}

#' @param state the shared app-state object (see [new_app_state()]).
#' @return Invisible NULL.
mod_feature_server <- function(id, state) {
  moduleServer(id, function(input, output, session) {
    meta_editor_server("editor", state, .feature_editor_opts)

    output$detected_id <- renderText({
      req(state$working)
      sprintf("Detected: %s", detect_id_type(rownames(state$working)))
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
      newft <- detect_feature_type(state$working)
      state$meta <- utils::modifyList(state$meta, list(feature_type = newft$feature_type))
      cov <- annotation_coverage(state$working, paste0(newft$feature_type, "_name"))
      showNotification(sprintf("Annotated %d of %d endogenous features.", cov$matched, cov$total),
                       type = "message")
    })

    output$coverage <- renderText({
      req(state$working)
      name_col <- paste0(state_meta(state)$feature_type, "_name")
      cov <- annotation_coverage(state$working, name_col)
      if (cov$matched == 0L) "Not annotated yet."
      else sprintf("Annotated %d of %d endogenous features%s.", cov$matched, cov$total,
                   if (cov$matched < cov$total) " (a GTF upload can fill gaps; coming soon)" else "")
    })

    invisible(NULL)
  })
}
