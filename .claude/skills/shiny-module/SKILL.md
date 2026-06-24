---
name: shiny-module
description: Scaffold a new page/feature in this RNA-seq dashboard as a Shiny module following the project's conventions — paired mod_<name>_ui/server functions, the shared app-state object (original/working/data_version/history/undo_stack/meta + a derived cache), versioned cached derivations, global + scoped undo/reset, bslib layout, the deferred-render + stale-note pattern, the reusable "Plot Showing" subset control, the removal-pool select-rows→commit idiom, and the composable sub-module / draft-editor patterns. Use when adding a new page, splitting a page into sub-modules, or refactoring inline server logic into a module.
metadata:
  type: project
  version: "2.2"
---

# Scaffolding a dashboard module

Every page in this app is a Shiny module so state and namespacing stay isolated. Pair with the `shiny-bslib` skill for layout and `rnaseq-bioc` for any data logic.

## The shared app-state object

The canonical store is **not** a bare `dds` reactive — it's one `reactiveValues` (`state`) threaded into every module (see `CLAUDE.md` → State model). Modules interact with it through a small helper API (lives in `R/state.R`):

- `state$working` — current `dds`; `state$original` — immutable load; `state$data_version` — bumped on any data edit; `state$history` — action log; `state$undo_stack` — short snapshot stack (depth 5); `state$meta` — app-level flags (`data_type`, `feature_type`, `sce_per_cell`) not stored on the `dds`. `state$derived` — cached heavy artifacts, a plain **environment** (not a reactive field), staleness keyed on `data_version`.
- `state_meta(state)` — summary list for the status bar (`loaded`, `data_type`, `feature_type`, `n_features`/`n_samples`, `assays`, `design`, `n_edits`, `n_undo`, feature-class counts `n_endogenous`/`n_spike_in`/`n_exogenous`, `data_version`). Read `meta$feature_type` here for adaptive `<unit>_name` labels.
- `state_dds(state)` — reactive read of `state$working` (use this where the old code read `dds()`).
- `state_mutate(state, fn, action)` — apply `fn(working)` → new `dds`, **bump `data_version`**, append `action` (name + params) to `history`. The single entry point for edits (filters, metadata, annotation).
- `state_derive(state, key, params, expr)` — get-or-compute a cached artifact (VST, PCA, DESeq fit, DE table) stamped with the current `data_version`; recomputes only when stale. Wrap heavy work behind the render button **and** this helper.
- `state_reset(state)` / `state_undo(state)` — restore `original` / pop the snapshot stack. **Wired globally** as the status-bar Undo / Reset buttons (`mod_statusbar`); every committed edit is undoable. **Scoped resets** complement them: `reset_metadata_slot(working, original, slot)` (the metadata editor's "Reset to original" — reverts one slot, auto-commits) and `restore_samples()`/`restore_features()` (QC "Reset Sample/Feature Removal" — re-add removed items from `original`, keeping other edits). All are themselves `state_mutate`s (undoable + logged).

This keeps invalidation in one place: a module never hand-rolls staleness — mutating via `state_mutate` invalidates everything downstream automatically.

## The contract

- **Naming:** `mod_<page>_ui(id, ...)` and `mod_<page>_server(id, state, ...)`. File `R/mod_<page>.R`.
- **Namespacing:** all input IDs go through `ns <- NS(id)` in the UI; bare IDs inside `moduleServer()`.
- **Edits go through `state_mutate`**; derived results go through `state_derive`. Modules return `invisible(NULL)` — they mutate shared `state`, they don't pass a `dds` back.
- **Deferred rendering:** expensive work is gated behind an "apply/render" button (`bindEvent()`/`eventReactive()`), never live on every control. The QC pages use a reusable helper `deferred(auto_id, render_id, spec, sig)` returning `list(value, stale)`: `spec()` is the heavy data reactive; `sig` is a *cheap* reactive signature of the inputs the plot depends on; the plot re-pulls on the Render button (or live when an `auto_id` checkbox is on), and `stale()` is TRUE when `sig` changed since the last render. Pair it with `stale_note()` to render a "Settings changed — click Render" banner above the plot.
- **Stale-aware UI:** show a stale badge / disable "use" when a needed `derived` artifact is stale; gate on prerequisites (e.g. DE needs a design) with a banner.
- **ggplot↔plotly engine toggle (dual-output pattern):** a global status-bar switch (`mod_statusbar`) sets `state$plot_interactive` (default `FALSE` = static ggplot). A toggleable plot is **a `gg(interactive)` builder fn feeding `dual_plot("<id>", gg, n_elements, height)`** (in `mod_qc.R`): it registers a static `plotOutput` (always `gg(FALSE)` — no plotly-only aes) + an interactive `plotlyOutput` (`plotly::ggplotly(gg(TRUE), tooltip = "text")`) behind a `uiOutput` container that shows whichever the per-plot decision selects (spinner inside). UI placeholder: `.qc_dual_plot(ns("<id>_container"))`. Interactivity needs `use_plotly_base()` (toggle ON **and** `requireNamespace("plotly")` — a `Suggests`, graceful static fallback) and is then **gated per plot on an element budget**: `n_elements()` (a cheap estimate ≈ rows of the plotted data, the real `ggplotly` cost driver — *not* sample count) vs `.plotly_max_elements()` = `getOption("ddsdashboard.plotly_max_elements", 5000L)` (an `option()`, no settings page). Over budget → static fallback + a **sticky per-plot "Render interactive anyway"** override (a `forced` `reactiveVal` set by a namespaced button, reset on `state$data_version` change or any toggle flip — so flipping the global switch off/on is the way back to the default gate; never raise the cap). Builders add a hover **`text` aesthetic only when `interactive = TRUE`** (via `.hover_aes`, so static never warns; plotly path muffles the residual warning with `.muffle_unknown_aes`). `validate(need(...))` lives in `gg`; `renderPlot` shows it natively, the plotly path **catches** the condition into a message figure (`renderPlotly` would otherwise leak a widget error). **Excludes** ComplexHeatmap / non-ggplot outputs. Reuse for any future plot page (PCA/DE).
- **Read-only tables go through `dt_table()`** (`R/utils_table.R`) — the standard `DT::datatable()` wrapper (per-column filters, rows-per-page selector, search, horizontal scroll, `rownames = FALSE`). It defaults to **`selection = "none"`** (read-only); actionable tables opt in with `selection = list(mode = "multiple")`. Extend via its `options=`/`selection=`/`...` passthrough rather than hand-rolling. Per-value cell colour via `DT::formatStyle()` + `DT::styleEqual()` (e.g. danger-red `TRUE`s; `feature_class` orange/purple). The editable metadata editor (`mod_meta_editor`) builds its own `DT::datatable()` (editable cells + opt-in row selection) and is the exception.
- **Tooltips use `bslib::tooltip()`**, never the native `title=` attribute — consistent hover styling/delay app-wide.

## Composable sub-modules (established patterns)

A page can embed a reusable sub-module (`mod_<shared>.R`) that, unlike a page, **returns a
value** instead of `invisible(NULL)`:

- **Reader / picker → a settable `reactiveVal`.** `mod_gtf_reader_server()` returns its `gtf`
  `reactiveVal`; the host reads it (`gtf_obj()`) and tests inject by *setting* it
  (`gtf_obj(<GRanges>)`) instead of driving the file input. Keep it settable for testability.

- **Draft editor → `list(draft, set)` (compose-then-commit).** `mod_meta_editor_server()` keeps
  edits in a local `draft` `reactiveVal` and commits them in **one** `state_mutate()` on Save.
  It returns `list(draft = <reactive>, set = function(dds))`, so a host page composes extra edits
  (OrgDb/GTF annotation, sheet merge) onto the *same* draft via `editor$set(fn(editor$draft()))`
  — no `data_version` bump until Save. **Never write those composed edits straight to
  `state$working`:** the editor re-seeds its draft from `state$working`, so doing so silently
  drops the user's unsaved edits (a real bug we hit). Edits flow draft → Save → `state_mutate`.

- **"Plot Showing" subset → `showing_samples` reactive (the standard for plot tabs).**
  `R/mod_plot_subset.R`: put `plot_subset_ui(ns, <suffix>)` (a collapsible "Plot Showing" accordion:
  show-by column → keep-values) in *each* plot sidebar, and call `plot_subset_server(input, output,
  session, state, suffixes)` **once** — it syncs all instances to one canonical selection and returns a
  `showing_samples` reactive. It's **display-only**: plots filter their data by `showing_samples()`
  (and fold it into the deferred `sig` so a change marks the plot stale) — it never mutates the `dds`
  or bumps `data_version`. Any page with per-sample plots should embed it.

- **Removal-pool select-rows → commit (actionable tables).** QC Filtering stages a DT row selection,
  moves it into a *pool* `reactiveVal` via buttons (Add/Remove selected, Select all/Deselect all over
  the search-filtered rows), and **Applies** the pool in one `state_mutate` (a confirm modal for
  destructive removals). The same idiom drives the feature editor's bulk `feature_class` assign
  ("Set on selected rows" / "Set on filtered rows" via `input$<tbl>_rows_selected` / `_rows_all`).

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
  state <- new_app_state()           # original/working/data_version/history/undo_stack/meta + derived env
  mod_input_server("input",         state)   # Load: populates original + working
  mod_metadata_server("metadata",   state)   # Sample info (uses mod_meta_editor)
  mod_feature_server("feature",     state)   # Feature info: annotation (mod_gtf_reader + mod_meta_editor)
  mod_assay_server("assay",         state)   # normalized assays / size factors / feature_length
  mod_qc_server("qc",               state)
  mod_dimreduc_server("dimreduc",   state)
  mod_de_server("de",               state)
  mod_heatmap_server("heatmap",     state)
  mod_export_server("export",       state)   # reads working + history for the report
  mod_statusbar_server("statusbar", state)   # persistent data-type / n / stale badges
}
```

## Dynamic per-row UI with per-element observers (e.g. `mod_palette`)

When a module renders a variable set of rows/columns each with its own inputs + observers (the Palette page: one panel per metadata column, one colour picker per level), two hard-won rules apply — learned debugging `mod_palette` (P3g-a):

- **Register observers once; read the changing data *live* inside them. Do NOT destroy + re-register on every `data_version` change.** Tearing observers down and rebuilding them inside an `observeEvent(state$data_version, …)` *races with the deferred reactive flush*: the rebuild can land in the same flush as a user input event, so the freshly-created observer (with `ignoreInit`) skips the event while the old one is already destroyed — the write silently vanishes (and it's non-deterministic). Instead capture only the column *name*, fix the per-element observer *count* at registration, and read levels live: `cur_lvls <- function() .pal_levels(isolate(coldata_df())[[col]]); … lvls <- cur_lvls(); lev <- lvls[i]`. This keeps the index→level mapping correct after a dataset reload with no re-registration. (Documented limitation: a reload that *adds* levels beyond the original count leaves the extras palette-filled until the column is re-added.)
- **Per-row action observers need `ignoreInit = TRUE`.** A Reset/Remove `actionButton` carries a click counter; if you ever re-register the row's observers, the fresh observer fires immediately on the *stale* counter value and re-triggers the action (the "remove → can't re-add" bug). All the row's pin/select/reset/remove observers use `ignoreInit = TRUE`.
- Drive the panels off a `struct` `reactiveVal` bumped on add/remove/`data_version`, and read the config via `isolate()` inside the `renderUI` so per-element edits don't rebuild the panel mid-interaction; push programmatic value changes back to the widgets with `update*()` (e.g. `updateColorPickr`), not by re-rendering.
- **Normalize colours to 6-digit hex before comparing/storing them.** `viridisLite::viridis()` (and other ramp functions) return **8-digit hex with alpha** (`#440154FF`), but a colour picker (`shinyWidgets::colorPickr`) — and most "standard" colour inputs — echo back **6-digit** (`#440154`). So a "did the user edit this picker?" check of *stored-value vs picker-echo* mis-fires for those palettes (and only those — Brewer/Okabe-Ito are already 6-digit, so the bug looks intermittent). Run every resolved colour through a `col2rgb`→`rgb` normalizer (`norm_color()` here) at the resolver boundary so stored and echoed values compare equal. Same trap applies anywhere you diff a programmatic colour against a widget's reported colour.
- **`bslib::accordion_panel_close()` / `_open()` / `nav_select()` etc. take the *unnamespaced* id inside a module.** They send via `session$sendInputMessage`, which already applies the module namespace — so pass `"acc"` (the bare id you'd give `accordion(id = ns("acc"))` in the UI), not `ns("acc")`, or it double-namespaces and silently no-ops.

**Testing dynamic-UI modules with `shiny::testServer`:** the harness *batches* reactive flushes, so a `renderUI` rebuild can collide with a following `setInputs` in a way that never happens in the running app (where a dataset load and a widget click are separate event-loop cycles). Put `session$flushReact()` between a simulated load and a simulated interaction. And don't assert a fragile UI *round-trip* (input → observer → state → `update*()` → input) that only the real browser completes — assert the **deterministic data guarantee** instead (e.g. that the resolver fills a reloaded dataset's levels with no stale-name corruption). The core same-cycle behaviours (add, edit, select, reset, remove, re-add) *are* deterministic in `testServer` and should be covered directly.

## Checklist for a new module

- [ ] `mod_<page>_ui` namespaces every input with `ns()`; uses bslib (`layout_sidebar`, `card`, `layout_column_wrap`), not `fluidRow`/`column`.
- [ ] `mod_<page>_server(id, state)` reads via `state_dds()`, edits via `state_mutate()`, caches via `state_derive()`.
- [ ] Heavy work is gated behind an explicit button (`eventReactive`/`bindEvent`) and `req()`-guarded.
- [ ] No hand-rolled staleness — invalidation is implicit through `state_mutate` bumping `data_version`.
- [ ] Prerequisite banner + stale badge where a derived artifact is required.
- [ ] Per-sample plot tabs embed the "Plot Showing" control (`plot_subset_ui` per sidebar + one `plot_subset_server`); plots filter by `showing_samples()`.
- [ ] Any new mutating action logs to `history` (powers the reproducibility export).
