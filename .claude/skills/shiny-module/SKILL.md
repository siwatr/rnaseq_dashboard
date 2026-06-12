---
name: shiny-module
description: Scaffold a new page/feature in this RNA-seq dashboard as a Shiny module following the project's conventions — paired mod_<name>_ui/server functions, the shared app-state object (original/working/derived/data_version/history), versioned cached derivations, bslib layout, and deferred (button-gated) rendering. Use when adding a new page, splitting a page into sub-modules, or refactoring inline server logic into a module.
metadata:
  type: project
  version: "2.0"
---

# Scaffolding a dashboard module

Every page in this app is a Shiny module so state and namespacing stay isolated. Pair with the `shiny-bslib` skill for layout and `rnaseq-bioc` for any data logic.

## The shared app-state object

The canonical store is **not** a bare `dds` reactive — it's one `reactiveValues` (`state`) threaded into every module (see `CLAUDE.md` → State model). Modules interact with it through a small helper API (lives in `R/state.R`):

- `state$working` — current `dds`; `state$original` — immutable load; `state$derived` — cached heavy artifacts; `state$data_version` — bumped on any data edit; `state$history` — action log.
- `state_dds(state)` — reactive read of `state$working` (use this where the old code read `dds()`).
- `state_mutate(state, fn, action)` — apply `fn(working)` → new `dds`, **bump `data_version`**, append `action` (name + params) to `history`. The single entry point for edits (filters, metadata, annotation).
- `state_derive(state, key, params, expr)` — get-or-compute a cached artifact (VST, PCA, DESeq fit, DE table) stamped with the current `data_version`; recomputes only when stale. Wrap heavy work behind the render button **and** this helper.
- `state_reset(state)` / `state_undo(state)` — restore `original` / pop the snapshot stack.

This keeps invalidation in one place: a module never hand-rolls staleness — mutating via `state_mutate` invalidates everything downstream automatically.

## The contract

- **Naming:** `mod_<page>_ui(id, ...)` and `mod_<page>_server(id, state, ...)`. File `R/mod_<page>.R`.
- **Namespacing:** all input IDs go through `ns <- NS(id)` in the UI; bare IDs inside `moduleServer()`.
- **Edits go through `state_mutate`**; derived results go through `state_derive`. Modules return `invisible(NULL)` — they mutate shared `state`, they don't pass a `dds` back.
- **Deferred rendering:** expensive work is gated behind an "apply/render" button (`bindEvent()`/`eventReactive()`), never live on every control.
- **Stale-aware UI:** show a stale badge / disable "use" when a needed `derived` artifact is stale; gate on prerequisites (e.g. DE needs a design) with a banner.

## Template

```r
# R/mod_<page>.R

mod_<page>_ui <- function(id) {
  ns <- NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      # controls -> ns("...")
      actionButton(ns("render"), "Render", class = "btn-primary")
    ),
    bslib::card(
      bslib::card_header("..."),
      shinycssloaders::withSpinner(plotOutput(ns("plot")))
    )
  )
}

#' @param state the shared app-state reactiveValues
mod_<page>_server <- function(id, state) {
  moduleServer(id, function(input, output, session) {
    # Read current dds via the helper; compute heavy artifacts via state_derive.
    embedding <- eventReactive(input$render, {
      dds <- req(state_dds(state)())
      state_derive(state, key = "pca",
                   params = list(n_top = input$n_top, assay = input$assay),
                   expr = function() compute_pca(dds, input$n_top, input$assay))
    })

    output$plot <- renderPlot({ req(embedding()); ggplot2::ggplot(...) })

    # An edit (e.g. filtering) goes through state_mutate — bumps version + logs.
    observeEvent(input$apply_filter, {
      state_mutate(state,
        fn = function(dds) dds[keep_rows(dds, input$min_count), ],
        action = list(name = "filter_features", min_count = input$min_count))
    })

    invisible(NULL)
  })
}
```

## Top-level wiring

```r
# R/run_app.R
server <- function(input, output, session) {
  state <- new_app_state()           # reactiveValues: original/working/derived/data_version/history
  mod_input_server("input",       state)   # populates original + working
  mod_qc_server("qc",             state)
  mod_process_server("process",   state)
  mod_dimreduc_server("dimreduc", state)
  mod_de_server("de",             state)
  mod_heatmap_server("heatmap",   state)
  mod_export_server("export",     state)   # reads working + history for the report
  mod_statusbar_server("statusbar", state) # persistent data-type / n / stale badges
}
```

## Checklist for a new module

- [ ] `mod_<page>_ui` namespaces every input with `ns()`; uses bslib (`layout_sidebar`, `card`, `layout_column_wrap`), not `fluidRow`/`column`.
- [ ] `mod_<page>_server(id, state)` reads via `state_dds()`, edits via `state_mutate()`, caches via `state_derive()`.
- [ ] Heavy work is gated behind an explicit button (`eventReactive`/`bindEvent`) and `req()`-guarded.
- [ ] No hand-rolled staleness — invalidation is implicit through `state_mutate` bumping `data_version`.
- [ ] Prerequisite banner + stale badge where a derived artifact is required.
- [ ] Any new mutating action logs to `history` (powers the reproducibility export).
