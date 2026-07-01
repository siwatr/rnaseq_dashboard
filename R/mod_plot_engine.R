# Shared plot engine: the ggplot <-> plotly toggle, deferred (render-button)
# rendering, and the small plot helpers used by any plot-visualization page (QC
# today; PCA/DE next). Extracted from mod_qc.R so the pages agree on behaviour.
#
# Like mod_plot_subset, the server piece operates in the *host* module's
# namespace (it wires outputs/inputs the host's UI declares), so a page calls
# `eng <- plot_engine_server(input, output, session, state)` once and uses
# `eng$dual_plot()` / `eng$deferred()` / `eng$stale_note()` / `eng$use_plotly_base()`.
# The UI placeholders (`.plot_dual()`) and pure helpers are plain functions.

# A blank plot carrying a centered message (graceful degradation / empty state).
.plot_msg <- function(msg) {
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0, y = 0, label = msg) +
    ggplot2::theme_void()
}

# UI placeholder for a toggleable plot: a container the server fills with either
# a plotOutput (static) or a plotlyOutput (interactive). The spinner moves inside
# the container (built server-side) so it still fires during the actual render.
.plot_dual <- function(id) uiOutput(id)

# Per-plot interactivity is gated on an *element budget* (rows of the plotted
# data ~ rendered glyphs, the real ggplotly cost driver), not a sample cap - so a
# light per-sample plot goes interactive while a heavy features x samples plot
# falls back to static (with a per-plot override). The budget is an option() so
# deployments / power users can tune it without a settings page.
.plotly_max_elements <- function() {
  getOption("ddsdashboard.plotly_max_elements", 5000L)
}

# Muffle ggplot2's "Ignoring unknown aesthetics: text" - emitted at construction
# for the interactive builders' plotly-only hover aes; expected and noise-only.
.muffle_unknown_aes <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    if (grepl("unknown aesthetic", conditionMessage(w), ignore.case = TRUE))
      invokeRestart("muffleWarning")
  })
}

# Convert a ggplot to an interactive plotly figure. Tooltip is driven by the
# `text` aesthetic the builders add when interactive. Theming (thematic / dark
# mode) is intentionally NOT restyled here yet - a full theme pass is planned.
.to_plotly <- function(p) plotly::ggplotly(p, tooltip = "text")

# Shared plot theme. thematic recolors fg/bg/accent to follow the live bslib theme
# (incl. dark mode) at draw time; `dark_theme` is the explicit lever for element
# choices thematic does not manage (gridline + text contrast). Used by QC + PCA.
.plot_theme <- function(dark_theme = FALSE) {
  grid <- if (isTRUE(dark_theme)) "grey35" else "grey85"
  text <- if (isTRUE(dark_theme)) "grey90" else "gray5"
  ggplot2::theme() +
    ggplot2::theme(
      panel.background = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(linewidth = 1, colour = grid),
      panel.grid.minor = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_line(linewidth = 1, colour = grid),
      text = ggplot2::element_text(size = 20, color = text),
      legend.position = "bottom"
    )
}

#' Host-namespace plot-engine helpers (interactivity toggle + deferred render)
#'
#' Call once in a page's `moduleServer` body. Returns a list of closures that
#' wire toggleable, deferred plots in the host module's namespace.
#'
#' @param input,output,session The host module server arguments.
#' @param state The shared app-state object (see [new_app_state()]).
#' @return A list with `use_plotly_base` (reactive), `dual_plot(id, gg, n_elements,
#'   height)`, `deferred(auto_id, render_id, spec, sig)`, and `stale_note(d)`.
#' @export
plot_engine_server <- function(input, output, session, state) {
  ns <- session$ns

  # The global engine toggle resolves to "use plotly as the base engine" only if
  # the package is installed; whether a given plot actually renders interactive is
  # decided per plot against the element budget (see dual_plot).
  use_plotly_base <- reactive({
    isTRUE(state$plot_interactive) && requireNamespace("plotly", quietly = TRUE) &&
      !is.null(state$working)
  })

  # Wire a toggleable plot. `gg(interactive)` returns the ggplot (or validate()):
  # the static output always builds it with interactive = FALSE (no plotly-only
  # `text` aes, so no warning); the plotly output with interactive = TRUE.
  # `n_elements()` is a cheap reactive estimate of rendered glyphs (rows of the
  # plotted data). A plot renders interactive when the engine is on AND it is
  # within the element budget OR the user clicked its per-plot "render anyway"
  # (a sticky override, reset when the data changes or the global toggle flips -
  # so flipping the toggle off/on is the way back to the default gate). The
  # static fallback shows a one-line note + the override button above the plot.
  dual_plot <- function(id, gg, n_elements, height = "300px") {
    forced <- reactiveVal(FALSE)
    observeEvent(input[[paste0(id, "_force")]], forced(TRUE))
    observeEvent(state$data_version, forced(FALSE))     # new data -> re-gate
    observeEvent(state$plot_interactive, forced(FALSE)) # toggle flip -> re-gate

    over_budget <- reactive(isTRUE(use_plotly_base()) && n_elements() > .plotly_max_elements())
    interactive_here <- reactive(isTRUE(use_plotly_base()) &&
                                 (n_elements() <= .plotly_max_elements() || isTRUE(forced())))

    # `gg` is intentionally a plain function, not memoized: the expensive work
    # already sits behind `deferred()`/`state_derive`, and only one of the two
    # outputs is ever in the DOM, so there is no shared double-compute to cache.
    output[[paste0(id, "_static")]] <- renderPlot(gg(FALSE))
    if (requireNamespace("plotly", quietly = TRUE)) {
      # renderPlot natively turns a validate()/req() condition into its grey
      # message; renderPlotly does NOT, so catch it and render the message as a
      # figure. The plotly-only `text` aes warning is muffled at construction.
      output[[paste0(id, "_plotly")]] <- plotly::renderPlotly(.muffle_unknown_aes({
        p <- tryCatch(gg(TRUE), shiny.silent.error = function(e) {
          msg <- conditionMessage(e)
          .plot_msg(if (nzchar(msg)) msg else "Click Render (or enable auto-render).")
        })
        .to_plotly(p)
      }))
    }
    output[[paste0(id, "_container")]] <- renderUI({
      inner <- if (isTRUE(interactive_here())) {
        plotly::plotlyOutput(ns(paste0(id, "_plotly")), height = height)
      } else {
        plotOutput(ns(paste0(id, "_static")), height = height)
      }
      tagList(
        if (isTRUE(over_budget()) && !isTRUE(forced()) && !is.null(state$working))
          tags$div(class = "alert alert-secondary py-2 px-2 small mb-2",
            tags$div(sprintf(
              "This plot would draw ~%s plotting elements (%d samples) - interactive rendering may be slow, so it is shown as a static image.",
              format(n_elements(), big.mark = ","), ncol(state$working))),
            actionButton(ns(paste0(id, "_force")), "Render interactive anyway",
                         icon = icon("bolt"), class = "btn-sm btn-outline-primary mt-1")),
        shinycssloaders::withSpinner(inner, proxy.height = height))
    })
  }

  # Gate an expensive plot spec behind an explicit Render button (+ an auto-render
  # toggle). Returns `value` (a reactiveVal of the last computed spec) and `stale`
  # (TRUE when inputs changed since the last manual render). `sig()` is a cheap
  # signature of the inputs the spec depends on.
  deferred <- function(auto_id, render_id, spec, sig) {
    rv <- reactiveVal(NULL)
    last_sig <- reactiveVal(NULL)
    go <- function() { rv(spec()); last_sig(sig()) }
    observe({ if (isTRUE(input[[auto_id]])) go() })
    observeEvent(input[[render_id]], go())
    stale <- reactive({
      if (is.null(rv()) || isTRUE(input[[auto_id]])) return(FALSE)
      !isTRUE(all.equal(last_sig(), sig()))
    })
    list(value = rv, stale = stale)
  }
  # A "settings changed -> re-render" banner, shown above a plot when stale.
  stale_note <- function(d) renderUI({
    if (!isTRUE(d$stale())) return(NULL)
    tags$div(class = "alert alert-warning py-1 px-2 small mb-2 d-flex align-items-center gap-2",
             icon("triangle-exclamation"),
             "Settings changed - click Render to update the plot.")
  })

  list(use_plotly_base = use_plotly_base, dual_plot = dual_plot,
       deferred = deferred, stale_note = stale_note)
}
