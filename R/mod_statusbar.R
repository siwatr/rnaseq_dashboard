# Persistent dataset-status bar shown in the navbar on every page. Reads the
# app state via state_meta() and renders compact badges: data type, dimensions,
# assays, design, edit count, and a warning when single-cell was coerced per-cell.
# A separate, periodically-refreshed badge reports the R process memory (RSS).

mod_statusbar_ui <- function(id) {
  ns <- NS(id)
  tags$span(class = "d-flex gap-2 align-items-center",
            uiOutput(ns("status"), inline = TRUE),
            uiOutput(ns("actions"), inline = TRUE),
            uiOutput(ns("mem"), inline = TRUE),
            # Global plot-engine toggle (static, so it survives status re-renders).
            # Default off = static ggplot. Read app-wide via state$plot_interactive.
            bslib::tooltip(
              bslib::input_switch(ns("interactive"), "Interactive plots", value = FALSE),
              "Render exploratory plots with plotly (hover, zoom). Off = static ggplot. Heavy plots stay static with a per-plot 'render anyway' option."))
}

# A small icon button, disabled (native HTML, no shinyjs) when `enabled` is FALSE.
.sb_btn <- function(id, ic, cls, enabled, tip) {
  b <- actionButton(id, NULL, icon = icon(ic), class = paste("btn-sm py-0", cls))
  if (!isTRUE(enabled)) b <- tagAppendAttributes(b, disabled = NA, `aria-disabled` = "true")
  bslib::tooltip(b, tip)
}

.badge <- function(text, class = "text-bg-secondary") {
  tags$span(class = paste("badge rounded-pill", class), text)
}

# Human-readable bytes (binary units). Pure + testable.
.format_bytes <- function(bytes) {
  if (is.null(bytes) || is.na(bytes)) return(NA_character_)
  units <- c("B", "KB", "MB", "GB", "TB")
  i <- if (bytes <= 0) 1L else min(length(units), floor(log(bytes, 1024)) + 1L)
  sprintf("%.1f %s", bytes / 1024^(i - 1L), units[i])
}

# Current R process RSS in bytes, or NA when `ps` is unavailable.
.process_rss <- function() {
  if (!requireNamespace("ps", quietly = TRUE)) return(NA_real_)
  tryCatch(as.numeric(ps::ps_memory_info(ps::ps_handle())[["rss"]]),
           error = function(e) NA_real_)
}

mod_statusbar_server <- function(id, state) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    output$status <- renderUI({
      m <- state_meta(state)
      if (!isTRUE(m$loaded)) return(.badge("no dataset loaded", "text-bg-light"))
      badges <- list(
        .badge(m$data_type, "text-bg-info"),
        bslib::tooltip(
          .badge(sprintf("%d features | %d samples | %d assays",
                         m$n_features, m$n_samples, length(m$assays))),
          paste("Assays:", paste(m$assays, collapse = ", "))),
        .badge(paste("design:", m$design))
      )
      # Flag non-endogenous features (spike-in / exogenous) when present.
      if (isTRUE(m$n_spike_in > 0L) || isTRUE(m$n_exogenous > 0L)) {
        parts <- c(if (m$n_spike_in > 0L) sprintf("%d spike-in", m$n_spike_in),
                   if (m$n_exogenous > 0L) sprintf("%d exogenous", m$n_exogenous))
        badges <- c(badges, list(bslib::tooltip(
          .badge(paste(parts, collapse = " / "), "text-bg-light"),
          "Non-endogenous features present (excluded from normalization, variable-gene selection, and expression filtering).")))
      }
      # `N edits` = net distance from the original (drives Reset); `N undo limit`
      # = how many of those edits can still be stepped back (capped at the
      # snapshot depth), which explains why Undo stops before reaching 0 edits.
      if (m$n_edits > 0L) {
        badges <- c(badges, list(
          bslib::tooltip(.badge(sprintf("%d edits", m$n_edits), "text-bg-light"),
                         "Edits applied since the dataset was loaded; Reset reverts all of them."),
          bslib::tooltip(.badge(sprintf("%d undo limit", m$n_undo), "text-bg-light"),
                         sprintf("Undo steps available now (at most %d are kept).", .undo_depth))))
      }
      if (isTRUE(m$sce_per_cell)) {
        badges <- c(badges, list(.badge("per-cell: stats unreliable", "text-bg-warning")))
      }
      tags$span(class = "d-flex gap-1 align-items-center", badges)
    })

    # Global Undo / Reset (the originally-designed state affordances). Undo steps
    # back the last committed edit (any tab); Reset restores the loaded object.
    # Disabled when unavailable; hidden until a dataset is loaded.
    output$actions <- renderUI({
      m <- state_meta(state)
      if (!isTRUE(m$loaded)) return(NULL)
      tags$span(class = "d-flex gap-1 align-items-center",
        .sb_btn(ns("undo"), "rotate-left", "btn-outline-secondary",
                length(state$undo_stack) > 0L, "Undo last edit"),
        .sb_btn(ns("reset"), "arrows-rotate", "btn-outline-danger",
                m$n_edits > 0L, "Reset to original (undo all edits)"))
    })
    # Global plot-engine toggle -> shared state (read by plot modules, e.g. QC).
    observeEvent(input$interactive, state$plot_interactive <- isTRUE(input$interactive),
                 ignoreNULL = FALSE)
    observeEvent(input$undo, state_undo(state))
    observeEvent(input$reset, {
      showModal(modalDialog(
        title = "Reset to original?",
        "Discard ALL edits and restore the dataset exactly as loaded? The undo history is cleared.",
        easyClose = TRUE,
        footer = tagList(modalButton("Cancel"),
                         actionButton(session$ns("reset_confirm"), "Reset", class = "btn-danger"))))
    })
    observeEvent(input$reset_confirm, { removeModal(); state_reset(state) })

    # Session memory footprint (whole R process), refreshed on its own timer so it
    # does not re-render the dataset badges. Hidden when `ps` is unavailable.
    output$mem <- renderUI({
      invalidateLater(5000, session)
      rss <- .process_rss()
      if (is.na(rss)) return(NULL)
      .badge(paste("session memory:", .format_bytes(rss)), "text-bg-light")
    })
  })
}
