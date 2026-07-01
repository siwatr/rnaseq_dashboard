# Page 5: Differential expression. A navset_card_tab of Design & Contrasts (P5b) /
# DE Plots / Results Table (P5c). The design is committed to the dds via the shared
# mod_design_builder (synced with the Input/Design tab). The DESeq2 FIT and the
# per-contrast RESULT EXTRACTION are separate: Run DESeq2 fits (contrast-free);
# extraction is reactive (Auto-update) or on demand, cached per contrast. Contrast
# specs live in state$de and carry a validity tier (extractable/not_in_design/
# invalid). See the shiny-de-analysis skill.

# Badge classes + tips for the three contrast-validity tiers.
.de_tier_class <- c(extractable = "text-bg-primary",
                    not_in_design = "text-bg-warning",
                    invalid = "text-bg-danger")
.de_tier_tip <- c(
  extractable   = "In the current design - produces results.",
  not_in_design = "Factor is valid but not in the current design - add it to the design to extract (kept for later).",
  invalid       = "Factor or a level no longer exists in the data - kept, and recoverable if you restore it.")

# A small colour legend for the contrast badges, derived from the tier map so it
# can't drift from the badges themselves.
.de_tier_label <- c(extractable = "in design", not_in_design = "not in design",
                    invalid = "invalid")
.de_legend <- function() {
  item <- function(cls, txt) tags$span(class = "d-inline-flex align-items-center me-2",
    tags$span(class = paste("badge rounded-pill me-1", cls), " "), txt)
  tags$div(class = "small text-muted mt-2",
    lapply(names(.de_tier_label),
           function(t) item(.de_tier_class[[t]], .de_tier_label[[t]])))
}

# --- Design & Contrasts tab UI ---------------------------------------------
.de_design_ui <- function(ns) {
  bslib::layout_sidebar(
    # This tab is a FORM whose height grows with the number of contrasts, so opt
    # out of the fill/flex container (which would compress + overlap the content):
    # let it flow at natural height and scroll with the page.
    fillable = FALSE,
    sidebar = bslib::sidebar(
      title = "Design & fit", width = 360,
      mod_design_builder_ui(ns("design")),
      tags$hr(),
      tags$div(class = "fw-semibold mb-2", "Differential expression fit"),
      bslib::tooltip(
        selectInput(ns("shrink"), "LFC shrinkage",
                    c("apeglm (recommended)" = "apeglm", "ashr" = "ashr", "none" = "none"),
                    selected = "apeglm"),
        "Shrinks noisy log2 fold-changes toward zero for low-count genes. apeglm (recommended) is stable but needs the contrast as a model coefficient (control = reference level); ashr works for any contrast; none = raw LFCs."),
      uiOutput(ns("run_note")),                       # fit status, above the button
      tags$div(class = "d-grid mt-2",
        bslib::tooltip(
          actionButton(ns("run"), "Run DESeq2", class = "btn fw-semibold",
                       style = "background-color:#8b58db;border-color:#8b58db;color:#fff;"),
          "Fits DESeq2 on the current design (this can take a while). Re-run after a data or design change."))
    ),
    bslib::card(
      fill = FALSE,
      bslib::card_header(tags$h4("Contrasts", class = "fs-6 mb-0")),

      ## Row 1 -- Add contrast
      tags$h5("Add contrast", class = "fs-6"),
      tags$p(class = "text-muted small",
             "A contrast compares a test level against a control (reference) level of a design factor."),
      bslib::layout_columns(
        col_widths = c(4, 4, 4),
        selectInput(ns("c_var"), "Factor", choices = NULL),
        selectInput(ns("c_test"), "Test", choices = NULL),
        selectInput(ns("c_control"), "Control", choices = NULL)
      ),
      actionButton(ns("c_add"), "Add contrast", class = "btn-secondary btn-sm"),
      tags$hr(),

      ## Row 2 -- Defined contrasts + Remove
      bslib::layout_columns(
        col_widths = c(7, 5),
        tags$div(
          tags$h5("Defined contrasts", class = "fs-6"),
          .de_legend(),                               # legend right under the header
          uiOutput(ns("contrast_list"))),             # badges below the legend
        tags$div(
          tags$h5("Remove contrasts", class = "fs-6"),
          selectizeInput(ns("remove_multi"), NULL, choices = NULL, multiple = TRUE,
                         options = list(placeholder = "select contrasts to remove")),
          tags$div(class = "d-flex gap-2 flex-wrap",
            actionButton(ns("remove_sel"), "Remove selected", class = "btn-secondary btn-sm"),
            actionButton(ns("remove_invalid"), "Remove invalid", class = "btn-danger btn-sm"),
            actionButton(ns("remove_all"), "Remove all", class = "btn-danger btn-sm")))
      ),
      tags$hr(),

      ## Row 3 -- DEG summary (extraction controls under the header, table below)
      tags$h5("DEG summary", class = "fs-6"),
      bslib::tooltip(
        bslib::input_switch(ns("auto_update"), "Auto-update results", value = TRUE),
        "On: results refresh automatically as you add/remove contrasts (no re-fit needed). Off: use the Update results button."),
      uiOutput(ns("update_btn")),
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

    # The fitted DESeqDataSet lives in the derived env (non-reactive) under a
    # stamp; extraction reads it without triggering a re-fit. It bypasses the
    # state_derive accessor deliberately: state_derive keys only on data_version,
    # but the fit must also invalidate on a design change (design_version) -- hence
    # the manual (dv, desv) stamp. NULL if no current fit.
    current_fit <- function() {
      if (!exists("de_fit", envir = state$derived, inherits = FALSE)) return(NULL)
      ent <- get("de_fit", envir = state$derived)
      cur <- list(dv = state$data_version, desv = state$design_version %||% 0L)
      if (!identical(ent$stamp, cur)) return(NULL)
      ent$value
    }

    # Extract results for every EXTRACTABLE contrast not already cached; prune
    # results whose contrast is gone / no longer extractable. Cheap + idempotent,
    # so it never loops (never re-triggers its own dependencies).
    do_extract <- function() {
      fit <- current_fit(); if (is.null(fit)) return(invisible())
      dds <- state$working; if (is.null(dds)) return(invisible())
      de <- state$de %||% list()
      # Extraction always uses the CURRENT shrinkage; if it changed since the last
      # extraction, drop the cached results so they recompute under the new method.
      shrink <- input$shrink %||% de$shrink %||% "apeglm"
      if (!identical(de$shrink, shrink)) { de$results <- list(); de$methods <- list() }
      de$shrink <- shrink
      specs <- de$contrasts %||% list()
      results <- de$results %||% list(); methods <- de$methods %||% list()
      extractable <- Filter(function(s) de_contrast_validity(dds, s) == "extractable", specs)
      ex_labels <- vapply(extractable, function(s) s$label, character(1))
      # Only the not-yet-cached contrasts cost anything (results()/lfcShrink can be
      # slow); show a progress spinner while they compute (auto or manual path).
      todo <- Filter(function(s) is.null(results[[s$label]]), extractable)
      if (length(todo)) {
        shiny::withProgress(message = "Extracting DE results", value = 0, {
          for (s in todo) {
            shiny::incProgress(1 / length(todo), detail = s$label)
            df <- de_results(fit, c(s$var, s$test, s$control), shrink_type = shrink)
            results[[s$label]] <- df
            methods[[s$label]] <- attr(df, "shrink_method") %||% "none"
          }
        })
      }
      results <- results[intersect(names(results), ex_labels)]
      methods <- methods[intersect(names(methods), ex_labels)]
      de$results <- results; de$methods <- methods
      if (is.null(de$active) || !(de$active %in% names(results))) {
        de$active <- if (length(results)) names(results)[1] else NULL
      }
      state$de <- de
    }
    auto_on <- function() isTRUE(input$auto_update %||% TRUE)
    maybe_extract <- function() if (auto_on()) do_extract()

    # Results out of date vs the controls (auto-off only): an extractable contrast
    # without a cached result, or a shrinkage change not yet applied, given a fit.
    needs_update <- function() {
      if (auto_on() || is.null(current_fit())) return(FALSE)
      de <- state$de %||% list(); dds <- state$working
      specs <- de$contrasts %||% list()
      ex <- vapply(specs, function(s) de_contrast_validity(dds, s) == "extractable", logical(1))
      ex_labels <- vapply(specs[ex], function(s) s$label, character(1))
      missing <- length(setdiff(ex_labels, names(de$results %||% list())))
      shrink_mismatch <- length(de$results %||% list()) > 0 &&
        !identical(de$shrink %||% "apeglm", input$shrink %||% "apeglm")
      missing > 0 || shrink_mismatch
    }

    # --- contrast pickers (reactive to data + design + factor) ------------
    observeEvent(list(state$data_version, state$design_version), {
      dds <- state$working
      if (is.null(dds)) return()
      dv <- intersect(tryCatch(all.vars(DESeq2::design(dds)), error = function(e) character(0)),
                      de_design_factors(dds))
      sel <- if (isTRUE(input$c_var %in% dv)) input$c_var else if (length(dv)) dv[1] else NULL
      updateSelectInput(session, "c_var", choices = dv, selected = sel)
    }, ignoreNULL = FALSE)

    # Levels refresh on factor OR data/design change (a dataset switch or a value
    # relabel keeps c_var == "condition" but changes its levels).
    observeEvent(list(input$c_var, state$data_version, state$design_version), {
      dds <- state$working
      if (is.null(dds) || is.null(input$c_var) || !nzchar(input$c_var)) return()
      lv <- de_contrast_levels(dds, input$c_var)
      test_sel <- if (isTRUE(input$c_test %in% lv)) input$c_test else if (length(lv)) lv[length(lv)] else NULL
      ctrl_sel <- if (isTRUE(input$c_control %in% lv)) input$c_control else if (length(lv)) lv[1] else NULL
      updateSelectInput(session, "c_test", choices = lv, selected = test_sel)
      updateSelectInput(session, "c_control", choices = lv, selected = ctrl_sel)
    })

    # --- add / remove contrasts -------------------------------------------
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
      maybe_extract()
    })

    drop_labels <- function(labels) {
      de <- state$de %||% list()
      keep <- !(.de_labels(state) %in% labels)
      de$contrasts <- (de$contrasts %||% list())[keep]
      de$results <- de$results[setdiff(names(de$results %||% list()), labels)]
      de$methods <- de$methods[setdiff(names(de$methods %||% list()), labels)]
      if (!is.null(de$active) && de$active %in% labels) {
        de$active <- if (length(de$results)) names(de$results)[1] else NULL
      }
      state$de <- de
    }
    observeEvent(input$remove_sel, {
      if (length(input$remove_multi)) drop_labels(input$remove_multi)
    })
    observeEvent(input$remove_all, {
      de <- state$de %||% list()
      de$contrasts <- list(); de$results <- list(); de$methods <- list(); de$active <- NULL
      state$de <- de
    })
    observeEvent(input$remove_invalid, {
      dds <- state$working; if (is.null(dds)) return()
      specs <- (state$de %||% list())$contrasts %||% list()
      bad <- vapply(specs, function(s) de_contrast_validity(dds, s) == "invalid", logical(1))
      if (any(bad)) drop_labels(vapply(specs[bad], function(s) s$label, character(1)))
    })

    # Keep the Remove multi-select choices synced with the stored contrasts.
    observeEvent(.de_labels(state), {
      updateSelectizeInput(session, "remove_multi", choices = .de_labels(state), server = FALSE)
    }, ignoreNULL = FALSE)

    # Notify (keep, don't auto-remove) when a data/design change leaves contrasts
    # non-extractable. Load resets state$de, so no false alarm on load.
    observeEvent(list(state$data_version, state$design_version), {
      dds <- state$working; if (is.null(dds)) return()
      specs <- (state$de %||% list())$contrasts %||% list()
      if (!length(specs)) return()
      nbad <- sum(vapply(specs, function(s) de_contrast_validity(dds, s) != "extractable", logical(1)))
      if (nbad > 0) {
        showNotification(sprintf(
          "%d contrast(s) are not extractable under the current design/data (kept - see the Contrasts card; remove them there if you don't need them).", nbad),
          type = "warning", duration = 6)
      }
    }, ignoreInit = TRUE)

    output$contrast_list <- renderUI({
      dds <- state$working
      specs <- (state$de %||% list())$contrasts %||% list()
      if (!length(specs)) return(tags$p(class = "text-muted small", "No active contrast."))
      badges <- lapply(specs, function(s) {
        v <- if (is.null(dds)) "invalid" else de_contrast_validity(dds, s)
        bslib::tooltip(
          tags$span(class = paste("badge rounded-pill me-1 mb-1", .de_tier_class[[v]]),
                    style = "font-size:1rem;", s$label),
          .de_tier_tip[[v]])
      })
      tags$div(class = "d-flex flex-wrap mt-3", badges)
    })

    # --- run DESeq2 (fit only) --------------------------------------------
    observeEvent(input$run, {
      dds <- state$working
      if (is.null(dds)) { showNotification("Load a dataset first.", type = "error"); return() }
      if (!isTRUE(design_status()$ok)) {
        showNotification("Set a valid (full-rank) design first (Apply design).", type = "error"); return()
      }
      shrink <- input$shrink %||% "apeglm"
      ok <- tryCatch({
        shiny::withProgress(message = "Running DESeq2", value = 0.4, {
          fit <- de_run(dds)
          stamp <- list(dv = state$data_version, desv = state$design_version %||% 0L)
          assign("de_fit", list(value = fit, stamp = stamp), envir = state$derived)
          de <- state$de %||% list()
          de$stamp <- stamp; de$shrink <- shrink
          de$results <- list(); de$methods <- list()   # fresh fit -> re-extract
          state$de <- de
        })
        TRUE
      }, error = function(e) {
        showNotification(paste("DESeq2 failed:", conditionMessage(e)),
                         type = "error", duration = NULL); FALSE
      })
      if (ok) {
        do_extract()                                    # always extract after a fresh fit
        showNotification("DESeq2 fit complete.", type = "message", duration = 3)
      }
    })

    # A shrinkage change needs re-extraction (new shrunk LFC), not a re-fit.
    # do_extract() clears + recomputes when it sees the mismatch; with auto off the
    # results stay (under the old shrink) until Update results, and needs_update()
    # flags them out of date.
    observeEvent(input$shrink, if (auto_on()) do_extract(), ignoreInit = TRUE)

    # Flipping Auto-update on catches results up.
    observeEvent(input$auto_update, if (auto_on()) do_extract(), ignoreInit = TRUE)
    observeEvent(input$update_results, do_extract(), ignoreInit = TRUE)

    output$update_btn <- renderUI({
      if (auto_on()) return(NULL)
      # Opaque; primary when an update is pending, secondary otherwise.
      cls <- if (needs_update()) "btn-primary btn-sm mt-1" else "btn-secondary btn-sm mt-1"
      actionButton(session$ns("update_results"), "Update results",
                   icon = icon("arrows-rotate"), class = cls)
    })

    output$run_note <- renderUI({
      switch(de_status(state),
        none    = tags$div(class = "text-muted small mb-2", "No DESeq2 results available."),
        stale   = tags$div(class = "text-warning small mb-2",
                           icon("triangle-exclamation"), " Re-run DESeq2 needed (data or design changed)."),
        current = tags$div(class = "text-success small mb-2", "Results are current."))
    })

    # --- per-contrast DEG summary -----------------------------------------
    output$summary <- renderUI({
      res <- (state$de %||% list())$results %||% list()
      # Echo the staleness the sidebar/status bar show, so the DEG table below is
      # never read as current when the data/design/shrink moved on.
      note <- if (identical(de_status(state), "stale"))
                tags$p(class = "text-warning small mb-1",
                       "DEG counts below are stale (data or design changed) - re-run DESeq2.")
              else if (needs_update())
                tags$p(class = "text-warning small mb-1",
                       "Results are out of date - click Update results.")
              else NULL
      if (!length(res)) {
        msg <- if (!is.null(current_fit())) "Add a contrast to extract results."
               else "Run DESeq2 to see the DEG summary."
        return(tagList(note, tags$p(class = "text-muted small", msg)))
      }
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
      tags$div(
        note,
        tags$table(class = "table table-sm",
          tags$thead(tags$tr(tags$th("Contrast"), tags$th(class = "text-end", "Up"),
                             tags$th(class = "text-end", "Down"),
                             tags$th(class = "text-end", "Total"), tags$th("Shrinkage"))),
          tags$tbody(rows)),
        tags$p(class = "text-muted small",
               "Counts at padj < 0.05 and |log2FC| >= 1 (standard LFC); adjustable thresholds arrive with the plots (P5c). apeglm shrinkage needs the control level to be the factor's reference; other contrasts fall back to ashr (see the Shrinkage column)."))
    })

    invisible(NULL)
  })
}
