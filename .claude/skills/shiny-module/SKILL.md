---
name: shiny-module
description: Scaffold a new page/feature in this RNA-seq dashboard as a Shiny module following the project's conventions — paired mod_<name>_ui/server functions, the shared reactive dds threaded in and returned out, bslib layout, and deferred (button-gated) plot rendering. Use when adding a new page, splitting a page into sub-modules, or refactoring inline server logic into a module.
metadata:
  type: project
  version: "1.0"
---

# Scaffolding a dashboard module

Every page in this app is a Shiny module so state and namespacing stay isolated. Use this pattern for new pages and sub-panels. Pair with the `shiny-bslib` skill for layout and `rnaseq-bioc` for any data logic.

## The contract

- **One reactive `dds` flows through.** A module *reads* the current `dds` (passed in as a reactive) and *returns* an updated `dds` reactive. The top-level server owns the canonical copy and rewires it: `dds <- mod_qc_server("qc", dds = dds)`.
- **Naming:** `mod_<page>_ui(id, ...)` and `mod_<page>_server(id, dds, ...)`. File `R/mod_<page>.R`.
- **Namespacing:** all input IDs go through `ns <- NS(id)` in the UI; the server uses bare IDs inside `moduleServer()`.
- **Deferred rendering:** expensive plots/computations are gated behind an "apply/render" button via `bindEvent()` / `eventReactive()` — never recompute live on every control.

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
      plotOutput(ns("plot"))
    )
  )
}

#' @param dds reactive() returning the current DESeqDataSet
#' @return reactive() returning the (possibly updated) DESeqDataSet
mod_<page>_server <- function(id, dds) {
  moduleServer(id, function(input, output, session) {
    # Gate heavy work behind the button + current dds.
    result <- eventReactive(input$render, {
      req(dds())
      # ... compute on dds() ...
    })

    output$plot <- renderPlot({
      req(result())
      # ggplot2 ...
    })

    # Return the (updated) object so the parent can rewire canonical state.
    reactive(dds())
  })
}
```

## Top-level wiring

```r
# app.R / R/run_app.R
ui <- bslib::page_navbar(
  title = "dds dashboard",
  theme = bslib::bs_theme(version = 5),
  bslib::nav_panel("Input",      mod_input_ui("input")),
  bslib::nav_panel("QC",         mod_qc_ui("qc")),
  bslib::nav_panel("Process",    mod_process_ui("process")),
  bslib::nav_panel("DimReduc",   mod_dimreduc_ui("dimreduc")),
  bslib::nav_panel("DE",         mod_de_ui("de")),
  bslib::nav_panel("Heatmap",    mod_heatmap_ui("heatmap")),
  bslib::nav_panel("Export",     mod_export_ui("export"))
)

server <- function(input, output, session) {
  dds <- mod_input_server("input")          # produces the initial dds
  dds <- mod_qc_server("qc", dds = dds)      # each page may refine it
  dds <- mod_process_server("process", dds = dds)
  # dimreduc / de / heatmap read dds; export consumes it
}
```

## Checklist for a new module

- [ ] `mod_<page>_ui` namespaces every input with `ns()`.
- [ ] `mod_<page>_server` takes `dds` reactive, returns a `dds` reactive (even if unchanged).
- [ ] Heavy work is `eventReactive`/`bindEvent` on an explicit button.
- [ ] If the module edits samples/features, downstream assays + DESeq fit are flagged stale (see `rnaseq-bioc`).
- [ ] `req()` guards against missing `dds`/empty inputs.
- [ ] Layout uses bslib (`layout_sidebar`, `card`, `layout_column_wrap`), not legacy `fluidRow`/`column`.
