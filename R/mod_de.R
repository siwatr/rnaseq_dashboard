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

# The shared thresholds (+ optional contrast-to-view) accordion group. Rendered in
# every DE tab's sidebar with a per-tab `suffix` (design/plots/table); the server
# keeps the copies in sync (contrast-to-view via state$de$active; padj/lfc/shrunk
# via a canonical thr store), so the controls are consistent on all tabs. The ids
# are suffixed because a Shiny input id can't appear twice. `include_view = FALSE`
# drops the "Contrast to view" selector on the Design & Contrasts tab, where a
# single "view" contrast is meaningless (that tab lists all contrasts) -- the
# thresholds still apply, driving the DEG summary there.
.de_view_controls <- function(ns, suffix, include_view = TRUE) {
  bslib::accordion_panel(
    if (include_view) "Contrast & thresholds" else "DEG thresholds",
    icon = icon("circle-half-stroke"),
    if (include_view)
      selectInput(ns(paste0("view_", suffix)), "Contrast to view", choices = NULL),
    numericInput(ns(paste0("padj_", suffix)), "padj threshold",
                 value = 0.05, min = 0, max = 1, step = 0.01),
    numericInput(ns(paste0("lfc_", suffix)), "abs(log2FC) threshold",
                 value = log2(2), min = 0, step = 0.1),
    bslib::input_switch(ns(paste0("shrunk_", suffix)), "Use shrunk LFC", value = FALSE))
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
          "Fits DESeq2 on the current design (this can take a while). Re-run after a data or design change.")),
      tags$hr(),
      bslib::accordion(open = FALSE, .de_view_controls(ns, "design", include_view = FALSE))
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
      actionButton(ns("c_add"), "Add contrast", class = "btn-secondary btn-sm", width="10%"),
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
            actionButton(ns("remove_invalid"), "Remove invalid", class = "btn-danger btn-sm", style = "color:#fff;"),
            actionButton(ns("remove_all"), "Remove all", class = "btn-danger btn-sm", style = "color:#fff;")))
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

# Segmented plot-type control (pills when shinyWidgets is present, else inline radios).
.de_plot_type_control <- function(id) {
  choices <- c("MA", "Volcano", "Direct comparison")
  if (requireNamespace("shinyWidgets", quietly = TRUE)) {
    shinyWidgets::radioGroupButtons(id, label = NULL, choices = choices,
                                    selected = "MA", size = "sm", status = "primary")
  } else {
    radioButtons(id, label = NULL, choices = choices, selected = "MA", inline = TRUE)
  }
}

# One field-based axis-limit control (min/max), labelled with the axis it drives.
.de_clamp_field <- function(ns, prefix, label) {
  tagList(
    tags$div(class = "small text-muted mb-1", label),
    bslib::layout_columns(
      col_widths = c(6, 6),
      numericInput(ns(paste0(prefix, "_min")), "min", value = NA),
      numericInput(ns(paste0(prefix, "_max")), "max", value = NA)))
}

# --- DE Plots tab UI -------------------------------------------------------
.de_plots_ui <- function(ns) {
  pt <- ns("plot_type")
  shows <- function(types) paste(sprintf("input['%s'] == '%s'", pt, types), collapse = " || ")
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      title = "Plot controls", width = 340,
      bslib::accordion(
        open = c("Contrast & thresholds", "Appearance"),
        .de_view_controls(ns, "plots"),
        bslib::accordion_panel("Appearance", icon = icon("palette"),
          selectInput(ns("colour_by"), "Colour by",
                      c("DEG status" = "DEG", "None" = "__none__",
                        "log2FC" = "__lfc__", "baseMean" = "baseMean")),
          sliderInput(ns("point_size"), "Point size", min = 0.25, max = 10, value = 3, step = 0.25),
          sliderInput(ns("point_alpha"), "Point opacity", min = 0.1, max = 1, value = 0.85, step = 0.05),
          selectInput(ns("legend_pos"), "Legend position",
                      c(Right = "right", Left = "left", Top = "top", Bottom = "bottom", None = "none"),
                      selected = "right"),
          conditionalPanel(shows("Direct comparison"),
            tags$hr(),
            tags$div(class = "small text-muted mb-1", "Direct comparison"),
            bslib::input_switch(ns("fixed_ratio"), "Fix 1:1 aspect ratio", value = FALSE),
            expr_value_ui(ns, "direct"))),
        bslib::accordion_panel("Labels", icon = icon("tag"),
          bslib::input_switch(ns("show_labels"), "Show gene labels", value = FALSE),
          numericInput(ns("top_n"), "Label top N by padj", value = 0, min = 0, max = 100, step = 1),
          gene_search_ui(ns, "label", multiple = TRUE,
                         search_modes = c("exact", "contains", "regex"),
                         label = "Genes of interest")),
        bslib::accordion_panel("Axis limits", icon = icon("up-right-and-down-left-from-center"),
          tags$p(class = "small text-muted",
                 "Blank = auto; out-of-range points draw as triangles. A limit is shared across plot types wherever its field recurs."),
          conditionalPanel(shows(c("MA", "Volcano")),
            .de_clamp_field(ns, "clamp_lfc", "log2FC (MA y-axis / Volcano x-axis)")),
          conditionalPanel(shows("Volcano"),
            .de_clamp_field(ns, "clamp_neglogp", "-log10(padj) (Volcano y-axis)")),
          conditionalPanel(shows("MA"),
            .de_clamp_field(ns, "clamp_bm", "log10(baseMean) (MA x-axis)")),
          conditionalPanel(shows("Direct comparison"),
            .de_clamp_field(ns, "clamp_expr", "expression (Direct x & y)")))
      )
    ),
    bslib::card(
      full_screen = TRUE,
      bslib::card_header(
        tags$div(class = "d-flex align-items-center justify-content-between flex-wrap gap-2",
          tags$h4("DE plot", class = "fs-6 mb-0"),
          .de_plot_type_control(ns("plot_type")))),
      uiOutput(ns("plots_deg_summary")),
      # Render controls above the plot; button on its own line under the toggle.
      tags$div(class = "mb-2",
        checkboxInput(ns("plot_auto"), "Auto-render", value = TRUE),
        actionButton(ns("plot_render"), "Render", icon = icon("play"),
                     class = "btn-sm btn-primary")),
      uiOutput(ns("de_plot_stale")),
      .plot_dual(ns("de_plot_container"))
    )
  )
}

# --- Results Table tab UI --------------------------------------------------
.de_table_ui <- function(ns) {
  bslib::layout_sidebar(
    fillable = FALSE,
    sidebar = bslib::sidebar(
      title = "Table", width = 320,
      bslib::accordion(open = "Contrast & thresholds", .de_view_controls(ns, "table")),
      bslib::input_switch(ns("sig_only"), "Significant only", value = FALSE)),
    bslib::card(
      bslib::card_header(tags$h4("DE results", class = "fs-6 mb-0")),
      uiOutput(ns("table_deg_summary")),
      tags$h5("DESeq2 Results", class = "fs-6 mt-1"),
      DT::DTOutput(ns("de_table")))
  )
}

mod_de_ui <- function(id) {
  ns <- NS(id)
  bslib::navset_card_tab(
    id = ns("tabs"),
    title = tags$h3("Differential expression", class = "fs-6"),
    bslib::nav_panel("Design & Contrasts", .de_design_ui(ns)),
    bslib::nav_panel("DE Plots", .de_plots_ui(ns)),
    bslib::nav_panel("Results Table", .de_table_ui(ns))
  )
}

# Contrast-spec labels currently stored.
.de_labels <- function(state) {
  vapply((state$de %||% list())$contrasts %||% list(),
         function(s) s$label, character(1))
}

#' @param state the shared app-state object (see [new_app_state()]).
#' @param dark_mode reactive; `TRUE` in dark mode (drives plot contrast).
#' @return Invisible `NULL`.
#' @noRd
mod_de_server <- function(id, state, dark_mode = reactive(FALSE)) {
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
                   icon = icon("arrows-rotate"), class = cls, width = "10%")
    })

    # Reflects the FIT (has DESeq() been run and is it up to date?), not the
    # extracted results -- so a fit with no contrasts yet still reads "up to date".
    output$run_note <- renderUI({
      switch(de_fit_status(state),
        none    = tags$div(class = "text-muted small mb-2",
                           "No DESeq2 fit yet - run DESeq2 to enable results."),
        stale   = tags$div(class = "text-warning small mb-2",
                           icon("triangle-exclamation"),
                           " Fit is out of date (data or design changed) - re-run DESeq2."),
        current = tags$div(class = "text-success small mb-2",
                           icon("check"), " Fit is up to date."))
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
      deg_col <- if (isTRUE(thr$shrunk)) "DEG_shrunk" else "DEG"
      rows <- lapply(names(res), function(lab) {
        s <- de_summary(
          de_classify_table(res[[lab]], thr$padj %||% 0.05, thr$lfc %||% log2(2)), deg_col)
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
               sprintf("Counts at padj < %s and |log2FC| >= %s (%s LFC), using the thresholds shared across all DE tabs. apeglm shrinkage needs the control level to be the factor's reference; other contrasts fall back to ashr (see the Shrinkage column).",
                       thr$padj %||% 0.05, round(thr$lfc %||% log2(2), 3),
                       if (isTRUE(thr$shrunk)) "shrunk" else "standard")))
    })

    # ================= DE Plots + Results Table (P5c) =================
    dark <- function() isTRUE(dark_mode())
    eng <- plot_engine_server(input, output, session, state)
    dual_plot <- eng$dual_plot; deferred <- eng$deferred; stale_note <- eng$stale_note
    feature_type <- function() (state$meta %||% list())$feature_type %||% "feature"

    # --- shared "Contrast & thresholds" controls (synced across all 3 tabs) ---
    # One copy of the group lives in each tab's sidebar (design/plots/table). The
    # contrast-to-view is synced through state$de$active; the thresholds through a
    # canonical `thr` store. Guarded fan-out so the echoes never loop.
    de_view_suffixes <- c("design", "plots", "table")
    thr <- reactiveValues(padj = 0.05, lfc = log2(2), shrunk = FALSE)

    result_labels <- reactive(names((state$de %||% list())$results %||% list()))
    observeEvent(list(result_labels(), (state$de %||% list())$active), {
      labs <- result_labels()
      a <- (state$de %||% list())$active
      if (is.null(a) || !(a %in% labs)) a <- if (length(labs)) labs[1] else NULL
      for (sfx in de_view_suffixes) {
        updateSelectInput(session, paste0("view_", sfx), choices = labs, selected = a)
      }
    }, ignoreNULL = FALSE)
    set_active <- function(lab) {
      if (is.null(lab) || !nzchar(lab)) return()
      de <- state$de %||% list()
      if (!identical(de$active, lab)) { de$active <- lab; state$de <- de }
    }

    # Per-tab edits -> canonical (guarded so a fan-out echo is a no-op).
    lapply(de_view_suffixes, function(sfx) {
      observeEvent(input[[paste0("view_", sfx)]], set_active(input[[paste0("view_", sfx)]]))
      observeEvent(input[[paste0("padj_", sfx)]], {
        v <- input[[paste0("padj_", sfx)]]; if (!identical(thr$padj, v)) thr$padj <- v
      }, ignoreInit = TRUE)
      observeEvent(input[[paste0("lfc_", sfx)]], {
        v <- input[[paste0("lfc_", sfx)]]; if (!identical(thr$lfc, v)) thr$lfc <- v
      }, ignoreInit = TRUE)
      observeEvent(input[[paste0("shrunk_", sfx)]], {
        v <- input[[paste0("shrunk_", sfx)]]; if (!identical(thr$shrunk, v)) thr$shrunk <- v
      }, ignoreInit = TRUE)
    })
    # Canonical -> every copy (handler is isolated, so reading inputs here can't loop).
    observeEvent(thr$padj, for (sfx in de_view_suffixes) {
      id <- paste0("padj_", sfx)
      if (!identical(input[[id]], thr$padj)) updateNumericInput(session, id, value = thr$padj)
    })
    observeEvent(thr$lfc, for (sfx in de_view_suffixes) {
      id <- paste0("lfc_", sfx)
      if (!identical(input[[id]], thr$lfc)) updateNumericInput(session, id, value = thr$lfc)
    })
    observeEvent(thr$shrunk, for (sfx in de_view_suffixes) {
      id <- paste0("shrunk_", sfx)
      if (!identical(input[[id]], thr$shrunk))
        bslib::update_switch(id, value = thr$shrunk, session = session)
    })

    active_spec <- reactive({
      a <- (state$de %||% list())$active
      Find(function(s) identical(s$label, a), (state$de %||% list())$contrasts %||% list())
    })
    active_raw <- reactive({
      a <- (state$de %||% list())$active
      if (is.null(a)) return(NULL)
      (state$de$results %||% list())[[a]]
    })

    # --- DEG stats for the active contrast (title/subtitle + summary tables) ---
    # down/up/no_change counts from the chosen DEG column (shrunk toggle).
    deg_counts <- function(d) {
      col <- if (isTRUE(thr$shrunk)) "DEG_shrunk" else "DEG"
      v <- d[[col]]
      list(down = sum(v == "down", na.rm = TRUE),
           up   = sum(v == "up", na.rm = TRUE),
           nc   = sum(v == "no_change", na.rm = TRUE))
    }
    # The shrinkage method actually shown: the extraction method when viewing shrunk
    # LFC, else "None" (standard LFC is unshrunk). Title-cased "None".
    active_method <- reactive({
      if (!isTRUE(thr$shrunk)) return("None")
      m <- ((state$de %||% list())$methods %||% list())[[(state$de %||% list())$active]]
      if (is.null(m) || identical(m, "none")) "None" else m
    })
    # A compact one-row DEG summary table for the active contrast (auto-updating).
    # Shown above both the plot and the results table for visual consistency.
    deg_summary_tags <- function() {
      res <- active_raw()
      if (is.null(res)) {
        return(tags$p(class = "text-muted small mb-2", "No results for the selected contrast."))
      }
      d <- de_classify_table(res, thr$padj %||% 0.05, thr$lfc %||% log2(2))
      cnt <- deg_counts(d)
      cell <- function(x) tags$td(class = "text-end", x)
      tags$table(class = "table table-sm mb-2 w-auto",
        tags$thead(tags$tr(
          tags$th("Contrast"), tags$th(class = "text-end", "Down"),
          tags$th(class = "text-end", "Up"), tags$th(class = "text-end", "No change"),
          tags$th(class = "text-end", "Total"), tags$th("Shrinkage"))),
        tags$tbody(tags$tr(
          tags$td((state$de %||% list())$active), cell(cnt$down), cell(cnt$up),
          cell(cnt$nc), cell(cnt$down + cnt$up),
          tags$td(class = "text-muted small", active_method()))))
    }
    output$plots_deg_summary <- renderUI(deg_summary_tags())
    output$table_deg_summary <- renderUI(deg_summary_tags())

    # DEG palette (from the Palette page's Other -> DEG; NULL -> Pink-Blue default).
    deg_colors <- reactive({
      cfg <- (state$palette$other %||% list())$DEG
      palette_discrete(c("up", "down", "no_change"),
                       cfg$colors, cfg$name %||% "DEG: Pink-Blue", cfg$custom)
    })

    # --- deferred classified results (the DATA; display aesthetics stay live) ---
    de_spec <- reactive({
      res <- req(active_raw())
      de_classify_table(res, thr$padj %||% 0.05, thr$lfc %||% log2(2))
    })
    # The sig includes a results-identity token (the fit stamp + shrink method) so a
    # re-fit / re-extraction of the SAME contrast marks the plot stale even though
    # extraction never bumps data_version (auto-render off would otherwise show old
    # values). Display aesthetics stay OUT of the sig (cheap live re-plots).
    de_shown <- deferred("plot_auto", "plot_render", de_spec,
      sig = reactive(list((state$de %||% list())$active, thr$padj, thr$lfc,
                          state$data_version, (state$de %||% list())$stamp,
                          (state$de %||% list())$shrink)))
    output$de_plot_stale <- stale_note(de_shown)

    # Field-based axis clamps (a limit follows its field across plot types).
    clamp_ranges <- reactive({
      rng <- function(p) {
        lo <- suppressWarnings(as.numeric(input[[paste0(p, "_min")]] %||% NA))
        hi <- suppressWarnings(as.numeric(input[[paste0(p, "_max")]] %||% NA))
        if (is.na(lo) && is.na(hi)) NULL else c(lo, hi)
      }
      list(lfc = rng("clamp_lfc"), neglogp = rng("clamp_neglogp"),
           bm = rng("clamp_bm"), expr = rng("clamp_expr"))
    })

    # Direct-comparison: the shared assay/transform/pseudocount control, then group
    # means of that (transformed) assay over the active contrast's control/test.
    direct_value <- expr_value_server(input, output, session, state, "direct")
    direct_means <- reactive({
      s <- active_spec(); dds <- state$working; req(s, dds)
      x <- SummarizedExperiment::colData(dds)[[s$var]]
      ctrl <- colnames(dds)[!is.na(x) & x == s$control]
      test <- colnames(dds)[!is.na(x) & x == s$test]
      validate(need(length(ctrl) && length(test), "The contrast's groups are empty in the current data."))
      spec <- direct_value()
      present <- SummarizedExperiment::assayNames(dds)
      assay <- if (spec$assay %in% present) spec$assay
               else if ("logcounts" %in% present) "logcounts" else present[1]
      de_group_means(dds, assay, ctrl, test,
                     transform = spec$transform, pseudocount = spec$pseudocount)
    })

    # --- gene labels (top-N by padj + searched genes via the shared module) ---
    # The searched "Genes of interest" box is the shared gene-search module (multi
    # mode; DE now gets an explicit "Search by" column picker + case toggle). Its
    # `label_hint` output reports unmatched terms.
    label_search <- gene_search_server(input, output, session, state, "label",
                                       multiple = TRUE,
                                       search_modes = c("exact", "contains", "regex"))
    display_names <- function(ids) {
      rd <- SummarizedExperiment::rowData(state$working)
      fn <- paste0(feature_type(), "_name")
      if (fn %in% colnames(rd)) {
        nm <- as.character(rd[ids, fn]); ifelse(is.na(nm) | !nzchar(nm), ids, nm)
      } else ids
    }
    label_ids <- function(d) {
      ids <- character(0)
      n <- suppressWarnings(as.integer(input$top_n %||% 0)); if (is.na(n)) n <- 0L
      if (n > 0L) ids <- rownames(d)[utils::head(order(d$padj, na.last = NA), n)]
      unique(intersect(c(ids, label_search()$ids), rownames(d)))
    }

    # --- the plot builder + dual_plot ------------------------------------
    de_colour_for <- function(d, ids) {
      ckey <- input$colour_by %||% "DEG"
      if (identical(ckey, "__none__")) return(NULL)
      d2 <- d[match(ids, rownames(d)), , drop = FALSE]
      shr <- isTRUE(thr$shrunk)
      if (identical(ckey, "DEG")) {
        de_colour_resolve(d2, if (shr) "DEG_shrunk" else "DEG", deg_colors())
      } else if (identical(ckey, "__lfc__")) {
        de_colour_resolve(d2, if (shr) "log2FoldChange_shrunk" else "log2FoldChange")
      } else {
        de_colour_resolve(d2, ckey)                 # baseMean
      }
    }
    labels_df <- function(d, ptype, lfc_col) {
      if (!isTRUE(input$show_labels)) return(NULL)
      ids <- label_ids(d); if (!length(ids)) return(NULL)
      sub <- d[ids, , drop = FALSE]
      xy <- switch(ptype,
        "MA"      = list(x = log10(sub$baseMean), y = sub[[lfc_col]]),
        "Volcano" = list(x = sub[[lfc_col]], y = -log10(sub$padj)),
        "Direct comparison" = {
          gm <- direct_means(); m <- gm[match(ids, gm$id), ]; list(x = m$control, y = m$test)
        })
      data.frame(x = xy$x, y = xy$y, label = display_names(ids), stringsAsFactors = FALSE)
    }
    build_de_gg <- function(interactive) {
      d <- de_shown$value()
      validate(need(!is.null(d), "Click Render to draw the plot (or enable auto-render)."))
      ptype <- input$plot_type %||% "MA"
      shr <- isTRUE(thr$shrunk)
      lfc_col <- if (shr) "log2FoldChange_shrunk" else "log2FoldChange"
      validate(need(lfc_col %in% names(d) && any(is.finite(d[[lfc_col]])),
                    "Shrunk LFC is unavailable for this contrast - switch off 'Use shrunk LFC'."))
      rg <- clamp_ranges()
      ps <- input$point_size %||% 1.4; pa <- input$point_alpha %||% 0.85
      labs <- labels_df(d, ptype, lfc_col)
      gg <- switch(ptype,
        "MA" = de_ma_gg(d, lfc_col, de_colour_for(d, rownames(d)),
                        x_range = rg$bm, y_range = rg$lfc, labels = labs,
                        point_size = ps, point_alpha = pa, dark = dark(),
                        interactive = interactive),
        "Volcano" = de_volcano_gg(d, lfc_col, de_colour_for(d, rownames(d)),
                        x_range = rg$lfc, y_range = rg$neglogp, labels = labs,
                        point_size = ps, point_alpha = pa, dark = dark(),
                        interactive = interactive),
        "Direct comparison" = {
          gm <- direct_means(); s <- active_spec()
          de_direct_gg(gm, value_label = direct_value()$label,
                       control_label = s$control %||% "control",
                       test_label = s$test %||% "test",
                       colour = de_colour_for(d, gm$id),
                       x_range = rg$expr, y_range = rg$expr, labels = labs,
                       point_size = ps, point_alpha = pa, dark = dark(),
                       fixed_ratio = isTRUE(input$fixed_ratio), interactive = interactive)
        })
      # Embed the contrast + DEG stats in the plot itself (for figure export).
      s <- active_spec(); cnt <- deg_counts(d)
      ttl <- if (!is.null(s)) sprintf("%s: %s vs %s", s$var, s$test, s$control) else NULL
      sub <- sprintf("DEG: %d Down | %d Up | %d No Change\n(LFC Shrinkage: %s)",
                     cnt$down, cnt$up, cnt$nc, active_method())
      gg + .plot_theme(dark()) +
        ggplot2::labs(title = ttl, subtitle = sub) +
        ggplot2::theme(legend.position = input$legend_pos %||% "right")
    }
    dual_plot("de_plot", build_de_gg,
              n_elements = reactive({ v <- de_shown$value(); if (is.null(v)) 0L else nrow(v) }),
              height = "460px")

    # --- Results Table ----------------------------------------------------
    output$de_table <- DT::renderDT({
      res <- active_raw()
      validate(need(!is.null(res), "No results - run DESeq2 on the Design & Contrasts tab."))
      shr <- isTRUE(thr$shrunk)
      d <- de_classify_table(res, thr$padj %||% 0.05, thr$lfc %||% log2(2))
      lfc_col <- if (shr) "log2FoldChange_shrunk" else "log2FoldChange"
      deg_col <- if (shr) "DEG_shrunk" else "DEG"
      sig_col <- if (shr) "sig_shrunk" else "sig"
      rd <- SummarizedExperiment::rowData(state$working)
      fn <- paste0(feature_type(), "_name")
      out <- data.frame(id = rownames(d), stringsAsFactors = FALSE)
      if (fn %in% colnames(rd)) out[[fn]] <- as.character(rd[rownames(d), fn])
      se_col <- if (shr && "lfcSE_shrunk" %in% names(d)) "lfcSE_shrunk" else "lfcSE"
      out$baseMean  <- round(d$baseMean, 1)
      out[[lfc_col]] <- round(d[[lfc_col]], 3)
      out$lfcSE     <- round(d[[se_col]], 3)   # matches the chosen LFC (shrunk SE when shrunk)
      out$pvalue    <- signif(d$pvalue, 3)
      out$padj      <- signif(d$padj, 3)
      out$DEG       <- as.character(d[[deg_col]])
      if (isTRUE(input$sig_only)) out <- out[which(d[[sig_col]]), , drop = FALSE]
      cols <- deg_colors()
      DT::formatStyle(dt_table(out), "DEG",
        color = DT::styleEqual(names(cols), unname(cols)),
        fontWeight = DT::styleEqual(c("up", "down"), c("bold", "bold"), default = NULL))
    })

    invisible(NULL)
  })
}
