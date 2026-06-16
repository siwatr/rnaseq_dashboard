# Persistent dataset-status bar shown in the navbar on every page. Reads the
# app state via state_meta() and renders compact badges: data type, dimensions,
# assays, design, edit count, and a warning when single-cell was coerced per-cell.
# A separate, periodically-refreshed badge reports the R process memory (RSS).

mod_statusbar_ui <- function(id) {
  ns <- NS(id)
  tags$span(class = "d-flex gap-1 align-items-center",
            uiOutput(ns("status"), inline = TRUE),
            uiOutput(ns("mem"), inline = TRUE))
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
    output$status <- renderUI({
      m <- state_meta(state)
      if (!isTRUE(m$loaded)) return(.badge("no dataset loaded", "text-bg-light"))
      badges <- list(
        .badge(m$data_type, "text-bg-info"),
        .badge(sprintf("%s features x %s samples", m$n_features, m$n_samples)),
        .badge(paste("assays:", paste(m$assays, collapse = ", "))),
        .badge(paste("design:", m$design))
      )
      if (m$n_edits > 0L) badges <- c(badges, list(.badge(sprintf("%d edits", m$n_edits))))
      if (isTRUE(m$sce_per_cell)) {
        badges <- c(badges, list(.badge("per-cell: stats unreliable", "text-bg-warning")))
      }
      tags$span(class = "d-flex gap-1 align-items-center", badges)
    })

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
