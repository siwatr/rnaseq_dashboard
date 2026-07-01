# A reusable sub-module: the guided DESeq2 design builder. Embedded in BOTH the
# Input/Design tab (sets the stage early) and the DE page's Design & Contrasts
# tab. Both instances read/write design(working) via state_set_design(), so they
# stay in sync (a design applied in one shows in the other). It only sets the
# design + reference levels -- it never runs DE. Returns a reactive giving the
# committed design + its full-rank status, so a host (the DE page) can gate Run.
# See the shiny-de-analysis + shiny-module skills.

mod_design_builder_ui <- function(id, title = NULL) {
  ns <- NS(id)
  tagList(
    if (!is.null(title)) tags$h4(title, class = "fs-6"),
    selectInput(ns("primary"), "Variable of interest", choices = NULL),
    selectizeInput(ns("covariates"), "Covariates (optional)", choices = NULL,
                   multiple = TRUE,
                   options = list(placeholder = "none")),
    uiOutput(ns("reflevels")),
    tags$div(class = "mt-2 mb-1 small text-muted", "Model formula"),
    tags$pre(class = "small mb-2", textOutput(ns("formula"), inline = TRUE)),
    uiOutput(ns("rank_badge")),
    actionButton(ns("apply"), "Apply design", class = "btn-primary btn-sm mt-2")
  )
}

#' @param state the shared app-state object (see [new_app_state()]).
#' @return A reactive yielding `list(design, ok)` for the committed design.
#' @noRd
mod_design_builder_server <- function(id, state) {
  moduleServer(id, function(input, output, session) {

    # Re-sync the controls to the committed design whenever the dataset or the
    # design changes (so the two embedded instances mirror each other). Keyed on
    # both versions; design_version bumps on Apply, data_version on load/edit.
    sync_key <- reactive(list(state$data_version, state$design_version))
    observeEvent(sync_key(), {
      dds <- state$working
      if (is.null(dds)) return()
      fac <- de_design_factors(dds)
      dv  <- tryCatch(all.vars(DESeq2::design(dds)), error = function(e) character(0))
      dv  <- intersect(dv, fac)                       # only discrete design vars
      prim <- if (length(dv)) dv[1] else if (length(fac)) fac[1] else character(0)
      cov  <- setdiff(dv, prim)
      updateSelectInput(session, "primary", choices = fac, selected = prim)
      updateSelectizeInput(session, "covariates", choices = fac, selected = cov)
    }, ignoreNULL = FALSE)

    # The terms the user has selected (primary first; primary excluded from covs).
    terms <- reactive({
      p <- input$primary %||% character(0)
      p <- p[nzchar(p)]
      cov <- setdiff(input$covariates %||% character(0), p)
      unique(c(p, cov))
    })

    pending_design <- reactive({
      tm <- terms()
      if (!length(tm)) return(stats::as.formula("~ 1"))
      stats::as.formula(paste("~", paste(tm, collapse = " + ")))
    })

    output$formula <- renderText(paste(deparse(pending_design()), collapse = " "))

    # Reference (control) level selector per discrete term.
    output$reflevels <- renderUI({
      dds <- state$working
      if (is.null(dds)) return(NULL)
      cd <- SummarizedExperiment::colData(dds)
      sels <- lapply(terms(), function(col) {
        x <- cd[[col]]
        if (is.null(x) || !(is.factor(x) || is.character(x))) return(NULL)
        lv <- if (is.factor(x)) levels(droplevels(x)) else sort(unique(as.character(x)))
        cur <- if (is.factor(x)) levels(droplevels(x))[1] else lv[1]
        selectInput(session$ns(paste0("ref_", col)),
                    sprintf("Reference level: %s", col), choices = lv, selected = cur)
      })
      sels <- Filter(Negate(is.null), sels)
      if (!length(sels)) return(NULL)
      tagList(tags$div(class = "small text-muted", "Reference (control) levels"), sels)
    })

    rank_now <- reactive({
      dds <- state$working
      if (is.null(dds)) return(list(ok = FALSE, msg = "No dataset loaded."))
      de_full_rank(pending_design(), SummarizedExperiment::colData(dds))
    })

    output$rank_badge <- renderUI({
      rk <- rank_now()
      if (isTRUE(rk$ok)) {
        tags$span(class = "badge rounded-pill text-bg-success", "design is full rank")
      } else {
        tags$div(class = "text-danger small", rk$msg %||% "Design is not full rank.")
      }
    })

    observeEvent(input$apply, {
      dds <- state$working
      if (is.null(dds)) return()
      rk <- rank_now()
      if (!isTRUE(rk$ok)) {
        showNotification(rk$msg %||% "Design is not full rank.", type = "error")
        return()
      }
      cd <- SummarizedExperiment::colData(dds)
      refs <- list()
      for (col in terms()) {
        x <- cd[[col]]
        if (!is.null(x) && (is.factor(x) || is.character(x))) {
          v <- input[[paste0("ref_", col)]]
          if (!is.null(v) && nzchar(v)) refs[[col]] <- v
        }
      }
      state_set_design(state, pending_design(), relevel = refs, action = list(source = id))
      showNotification("Design applied.", type = "message", duration = 3)
    })

    # The committed design + rank status (for a host to gate on).
    reactive({
      dds <- state$working
      if (is.null(dds)) return(list(design = NULL, ok = FALSE))
      d <- tryCatch(DESeq2::design(dds), error = function(e) NULL)
      ok <- !is.null(d) &&
        isTRUE(de_full_rank(d, SummarizedExperiment::colData(dds))$ok)
      list(design = d, ok = ok)
    })
  })
}
