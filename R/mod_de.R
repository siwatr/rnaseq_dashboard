# Page 5: Differential expression. A navset_card_tab of Design & Contrasts (this
# PR, P5b) / DE Plots / Results Table (P5c). The design is committed to the dds via
# the shared mod_design_builder (synced with the Input/Design tab); contrasts are
# stored in state$de and the fit is cached via state_derive keyed on
# (data_version, design_version). See the shiny-de-analysis skill.

# --- Design & Contrasts tab UI ---------------------------------------------
.de_design_ui <- function(ns) {
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      title = "Design & fit", width = 340,
      mod_design_builder_ui(ns("design")),
      tags$hr(),
      selectInput(ns("shrink"), "LFC shrinkage",
                  c("apeglm (needs coefficient)" = "apeglm",
                    "ashr" = "ashr", "none" = "none"), selected = "apeglm"),
      bslib::tooltip(
        actionButton(ns("run"), "Run DESeq2", class = "btn-primary"),
        "Fits DESeq2 on the current design and extracts every stored contrast. Re-run after a data or design change."),
      uiOutput(ns("run_note"))
    ),
    bslib::card(
      bslib::card_header(tags$h4("Contrasts", class = "fs-6 mb-0")),
      tags$p(class = "text-muted small",
             "A contrast compares a test level against a control (reference) level of a factor."),
      bslib::layout_columns(
        col_widths = c(4, 3, 3, 2),
        selectInput(ns("c_var"), "Factor", choices = NULL),
        selectInput(ns("c_test"), "Test", choices = NULL),
        selectInput(ns("c_control"), "Control", choices = NULL),
        tags$div(class = "d-flex align-items-end",
                 actionButton(ns("c_add"), "Add", class = "btn-secondary btn-sm w-100"))
      ),
      uiOutput(ns("contrast_list")),
      bslib::layout_columns(
        col_widths = c(8, 4),
        selectInput(ns("active"), "Active contrast", choices = NULL),
        tags$div(class = "d-flex align-items-end",
                 actionButton(ns("remove_active"), "Remove",
                              class = "btn-outline-danger btn-sm w-100"))
      ),
      uiOutput(ns("summary"))
    )
  )
}

.de_placeholder <- function(text) {
  bslib::card(bslib::card_body(tags$p(class = "text-muted", text)))
}

mod_de_ui <- function(id) {
  ns <- NS(id)
  bslib::navset_card_tab(
    id = ns("tabs"),
    title = tags$h3("Differential expression", class = "fs-6"),
    bslib::nav_panel("Design & Contrasts", .de_design_ui(ns)),
    bslib::nav_panel("DE Plots",
                     .de_placeholder("MA / volcano / direct-comparison plots arrive in P5c.")),
    bslib::nav_panel("Results Table",
                     .de_placeholder("The results table arrives in P5c."))
  )
}

# Contrast-spec labels currently stored.
.de_labels <- function(state) {
  vapply((state$de %||% list())$contrasts %||% list(),
         function(s) s$label, character(1))
}

#' @param state the shared app-state object (see [new_app_state()]).
#' @return Invisible `NULL`.
#' @noRd
mod_de_server <- function(id, state) {
  moduleServer(id, function(input, output, session) {

    # The shared design builder (synced with the Input/Design tab via state).
    design_status <- mod_design_builder_server("design", state)

    # --- contrast pickers -------------------------------------------------
    # A contrast can only be on a factor that is IN the model (else DESeq2::results
    # errors), so restrict the Factor to the discrete design terms; preserve the
    # user's selection across unrelated data edits when it's still valid.
    observeEvent(list(state$data_version, state$design_version), {
      dds <- state$working
      if (is.null(dds)) return()
      dv <- intersect(tryCatch(all.vars(DESeq2::design(dds)), error = function(e) character(0)),
                      de_design_factors(dds))
      sel <- if (isTRUE(input$c_var %in% dv)) input$c_var else if (length(dv)) dv[1] else NULL
      updateSelectInput(session, "c_var", choices = dv, selected = sel)
    }, ignoreNULL = FALSE)

    observeEvent(input$c_var, {
      dds <- state$working
      if (is.null(dds) || is.null(input$c_var) || !nzchar(input$c_var)) return()
      lv <- de_contrast_levels(dds, input$c_var)
      if (length(lv) < 2) return()
      updateSelectInput(session, "c_test", choices = lv, selected = lv[length(lv)])
      updateSelectInput(session, "c_control", choices = lv, selected = lv[1])
    })

    observeEvent(input$c_add, {
      req(input$c_var, input$c_test, input$c_control)
      if (identical(input$c_test, input$c_control)) {
        showNotification("Test and control levels must differ.", type = "error"); return()
      }
      label <- sprintf("%s: %s vs %s", input$c_var, input$c_test, input$c_control)
      if (label %in% .de_labels(state)) {
        showNotification("That contrast is already added.", type = "warning"); return()
      }
      de <- state$de %||% list()
      de$contrasts <- c(de$contrasts %||% list(),
                        list(list(var = input$c_var, test = input$c_test,
                                  control = input$c_control, label = label)))
      state$de <- de
    })

    observeEvent(input$remove_active, {
      label <- input$active
      if (is.null(label) || !nzchar(label)) return()
      de <- state$de %||% list()
      keep <- .de_labels(state) != label
      de$contrasts <- (de$contrasts %||% list())[keep]
      de$results <- de$results[setdiff(names(de$results %||% list()), label)]
      if (identical(de$active, label)) {
        de$active <- if (length(de$results)) names(de$results)[1] else NULL
      }
      state$de <- de
    })

    # Keep the active-contrast selector in sync with the stored list.
    observeEvent(.de_labels(state), {
      labs <- .de_labels(state)
      sel <- (state$de %||% list())$active %||% (if (length(labs)) labs[1] else NULL)
      updateSelectInput(session, "active", choices = labs, selected = sel)
    }, ignoreNULL = FALSE)
    observeEvent(input$active, {
      if (is.null(input$active) || !nzchar(input$active)) return()
      de <- state$de %||% list(); de$active <- input$active; state$de <- de
    })

    output$contrast_list <- renderUI({
      specs <- (state$de %||% list())$contrasts %||% list()
      if (!length(specs)) return(tags$p(class = "text-muted small", "No contrasts added yet."))
      tags$div(class = "d-flex flex-wrap gap-1 mb-2",
        lapply(specs, function(s)
          tags$span(class = "badge rounded-pill text-bg-secondary", s$label)))
    })

    # --- run DESeq2 -------------------------------------------------------
    observeEvent(input$run, {
      dds <- state$working
      if (is.null(dds)) { showNotification("Load a dataset first.", type = "error"); return() }
      if (!isTRUE(design_status()$ok)) {
        showNotification("Set a valid (full-rank) design first (Apply design).", type = "error"); return()
      }
      specs <- (state$de %||% list())$contrasts %||% list()
      if (!length(specs)) { showNotification("Add at least one contrast.", type = "error"); return() }
      shrink <- input$shrink %||% "apeglm"
      # Only contrasts whose factor is in the model and whose levels still exist
      # can be extracted; a stale spec (factor dropped from the design, or a level
      # removed by a later edit) is skipped + reported rather than sinking the batch.
      dv <- intersect(tryCatch(all.vars(DESeq2::design(dds)), error = function(e) character(0)),
                      de_design_factors(dds))
      valid <- vapply(specs, function(s)
        s$var %in% dv && all(c(s$test, s$control) %in% de_contrast_levels(dds, s$var)),
        logical(1))
      if (!any(valid)) {
        showNotification("No stored contrast is valid for the current design.", type = "error")
        return()
      }
      ok <- tryCatch({
        shiny::withProgress(message = "Running DESeq2", value = 0.2, {
          # data_version is carried by state_derive's version stamp; only the
          # design_version needs to be an explicit param (state_derive has no
          # concept of it) so the fit refits on a design-only change.
          fit <- state_derive(state, "de_fit",
            params = list(desv = state$design_version %||% 0L),
            expr = function() de_run(dds))
          shiny::incProgress(0.5, message = "Extracting contrasts")
          results <- list(); methods <- list()
          for (s in specs[valid]) {
            df <- de_results(fit, c(s$var, s$test, s$control), shrink_type = shrink)
            results[[s$label]] <- df
            methods[[s$label]] <- attr(df, "shrink_method") %||% "none"
          }
          de <- state$de %||% list()
          de$results <- results
          de$methods <- methods
          de$active  <- if (!is.null(de$active) && de$active %in% names(results)) de$active else names(results)[1]
          de$shrink  <- shrink
          de$stamp   <- list(dv = state$data_version, desv = state$design_version %||% 0L)
          state$de <- de
        })
        TRUE
      }, error = function(e) {
        showNotification(paste("DESeq2 failed:", conditionMessage(e)),
                         type = "error", duration = NULL); FALSE
      })
      if (ok) {
        n_skip <- sum(!valid)
        showNotification(
          if (n_skip) sprintf("DE run complete (%d contrast(s) skipped as invalid for the current design).", n_skip)
          else "DE run complete.",
          type = if (n_skip) "warning" else "message", duration = 4)
      }
    })

    output$run_note <- renderUI({
      st <- de_status(state)
      if (identical(st, "none")) return(NULL)
      if (identical(st, "stale"))
        tags$div(class = "text-warning small mt-2",
                 "Results are stale (data/design changed) - re-run.")
      else tags$div(class = "text-success small mt-2", "Results are current.")
    })

    # --- per-contrast DEG summary (default thresholds; P5c makes them adjustable) ---
    output$summary <- renderUI({
      res <- (state$de %||% list())$results %||% list()
      if (!length(res)) return(NULL)
      methods <- (state$de %||% list())$methods %||% list()
      rows <- lapply(names(res), function(lab) {
        s <- de_summary(de_classify_table(res[[lab]]))
        tags$tr(
          tags$td(lab),
          tags$td(class = "text-end", s[["up"]]),
          tags$td(class = "text-end", s[["down"]]),
          tags$td(class = "text-end", s[["total"]]),
          tags$td(class = "text-muted small", methods[[lab]] %||% ""))
      })
      tags$div(class = "mt-3",
        tags$h5("DEG summary", class = "fs-6"),
        tags$table(class = "table table-sm",
          tags$thead(tags$tr(tags$th("Contrast"), tags$th(class = "text-end", "Up"),
                             tags$th(class = "text-end", "Down"),
                             tags$th(class = "text-end", "Total"), tags$th("Shrinkage"))),
          tags$tbody(rows)),
        tags$p(class = "text-muted small",
               "Counts at padj < 0.05 and |log2FC| >= 1 (standard LFC). Adjustable thresholds arrive with the plots (P5c)."))
    })

    invisible(NULL)
  })
}
