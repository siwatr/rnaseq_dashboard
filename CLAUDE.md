# CLAUDE.md

Guidance for working in this repository — the distilled, load-bearing version. The full design narrative (the *why*) is in [dev_ref/rough_design.md](dev_ref/rough_design.md); the detailed, living phase plan + status is in [dev_ref/roadmap.md](dev_ref/roadmap.md) (the source of truth for what's done / next).

## What this is

An R **Shiny** dashboard that lets people with little to no bioinformatics background **explore, manipulate, and visualize RNA-seq results**. The input is an R object — a `DESeqDataSet` (`dds`) or a `SingleCellExperiment` (`sce`) — or a **counts matrix + sample sheet** (CSV/TSV/XLSX) that the app assembles into a `dds`.

**Scope is bulk-first** (decided in review). Build the bulk path end-to-end first; keep the design general enough to add single-cell later. An `sce` is supported via two routes: **pseudobulk** aggregation by a sample grouping, or **per-cell coercion to `dds`** behind a *"statistically inaccurate and slow"* warning. **Force pseudobulk above ~1k cells.** (Why per-cell DESeq2 is unsound: [dev_ref/normalization-scran-vs-deseq.md](dev_ref/normalization-scran-vs-deseq.md).)

The minimum viable input is a `dds` carrying only a raw count matrix; everything else (sample metadata, gene metadata, feature lengths, normalized assays, DE results) can be added inside the app.

## Tech stack

- **R + Shiny**, UI built with **bslib** (Bootstrap 5). Use the `shiny-bslib` and `shiny-bslib-theming` skills for layout/theming. `page_navbar()` for the multi-page structure; `DT` for editable tables; `shinycssloaders` for spinners; `thematic` so ggplots match the theme.
- **Bioconductor:** `DESeq2`, `SummarizedExperiment`, `SingleCellExperiment`, `scater` (QC), `rtracklayer` (GTF), `edgeR` (`filterByExpr`), `apeglm`/`ashr` (LFC shrinkage), `AnnotationDbi` + `org.*.eg.db` (annotation), `scran` (single-cell, later).
- **Plots:** `ggplot2` is the default for all plots. `ComplexHeatmap` is the heatmap engine (do not substitute pheatmap/heatmaply).
- **Tables / IO:** `readr`/`readxl` in, `writexl` out; `rmarkdown` + `sessioninfo` for the reproducibility report.

Domain conventions (object shape, assays, normalization math, QC metrics, DE result columns, `feature_class`/`feature_length`) live in the **`rnaseq-bioc`** skill — consult it before touching data logic.

## Project layout

This is a **plain R package** named `ddsdashboard` (not golem). Reusable helpers live in `R/` and are exported so other lab projects can `library(ddsdashboard)` and call them; the Shiny app is just `run_app()`.

```
DESCRIPTION / NAMESPACE   # package metadata + exports (regenerate NAMESPACE with devtools::document())
environment.yml           # conda/mamba env (source of truth for deps)
R/
  run_app.R               # run_app(), app_ui(), app_server() — the navbar shell + module wiring
  state.R                 # new_app_state() + the state_*() helper API (the canonical store)
  mod_<page>.R            # one Shiny module per page (input, metadata, feature, assay, qc, dimreduc, de, heatmap, export)
  mod_<shared>.R          # reusable sub-modules embedded in pages (mod_meta_editor, mod_gtf_reader, mod_statusbar)
  <topic>_helpers.R       # pure, exported, tested helpers (annotation_, gtf_, metadata_, assay_, load_)
  utils_normalization.R   # cpm(), tpm(), fpkm(), logcounts_from_counts() — pure, dependency-free, tested
  utils_lookup.R          # lookup_feature() — feature id <-> name resolution
  utils_table.R           # dt_table() — the standard read-only DT wrapper (filters + rows-per-page)
  mock_data.R             # make_mock_dds() etc. — fixtures used by tests + the demo
data-raw/                 # scripts that build mock dds/sce fixtures for tests + demo
dev_ref/                  # developer reference notes (rough_design.md, roadmap.md, normalization-scran-vs-deseq.md, object-states-and-resets.md) — all build-ignored
tests/testthat/           # unit tests for helpers + module servers (shiny::testServer)
inst/extdata/             # bundled reference data (ERCC Mix1/Mix2 concentrations: ercc_concentrations.csv, see ERCC_SOURCE.md)
```

The file conventions: **`mod_<page>.R`** = a navbar page; **`mod_<shared>.R`** = a reusable sub-module embedded inside pages (e.g. the draft metadata editor, the GTF reader); **`<topic>_helpers.R`** = pure, Shiny-free, exported functions that the modules call (so they stay reusable and unit-testable).

Module contract and the shared-`state` wiring are documented in the **`shiny-module`** skill (which also covers the reusable sub-module + draft-editor patterns). New pages follow that pattern. Domain/statistical conventions (object shape, assays, normalization, QC, DE schema, `feature_class`/`feature_length`, ERCC) live in **`rnaseq-bioc`**; the **`annotation`** skill covers the OrgDb + GTF row-annotation write conventions (consult it before touching mapping/annotation logic). A **`de-analysis`** skill (design/contrast/shrinkage) may be added when Phase 5 lands. The **`bioc-reviewer`** agent reviews any data-logic change; the **`shiny-feature-builder`** agent implements a page/helper end-to-end.

**Dependency policy:** `DESCRIPTION` `Imports` holds the always-on framework deps the app loads unconditionally — currently `AnnotationDbi`, `bslib`, `DESeq2`, `DT`, `methods`, `readr`, `readxl`, `shiny`, `SummarizedExperiment`, `thematic`. Heavy or optional packages are loaded **on demand via `requireNamespace()`** (with graceful degradation / a `skip_if_not_installed()` in tests) and stay in `Suggests` *even though code calls them* — currently `ComplexHeatmap`, `GenomicRanges`, `IRanges`, `S4Vectors`, `rtracklayer`, `org.*.eg.db`, `ps`, `scater`, `SingleCellExperiment`. They're also pinned in `environment.yml`. Promote `Suggests → Imports` only when a feature loads the package unconditionally. The pure base-R utilities (`utils_normalization.R`) deliberately take no dependencies so they stay reusable and fast to test.

## Environment & running

Dependencies are managed with **mamba** (channels: `conda-forge`, `bioconda` only — never `defaults`/`anaconda`). The environment is `rnaseq_dashboard`.

```bash
mamba env create -f environment.yml      # first time
mamba env update -f environment.yml --prune   # after editing deps
conda activate rnaseq_dashboard
```

Dev loop (from the repo root, env active):

```r
devtools::load_all()      # source the package
run_app()                 # launch the dashboard
devtools::document()      # regenerate NAMESPACE/man from roxygen
devtools::test()          # run testthat
devtools::check()         # full R CMD check
```

Clone-and-run for end users: create the env, then `R -e 'devtools::load_all(); run_app()'` (or install once with `devtools::install_local(".")` and `library(ddsdashboard); run_app()`).

## Architecture

Multi-page dashboard. Each page is a **Shiny module** (`mod_<page>_ui()` / `mod_<page>_server()`) so pages stay isolated and testable. Build a new page with the **`shiny-module`** skill.

### State model (the load-bearing design)

App state is one `reactiveValues` (`new_app_state()` in `R/state.R`) interacted with through the `state_*()` helper API (`state_load/mutate/derive/reset/undo`, `state_dds`, `state_meta`):

- **`original`** — the object exactly as loaded; never mutated; the reset target.
- **`working`** — the current `dds` after edits/filters/normalization (always a `dds`, even for `sce` input via pseudobulk/coercion).
- **`data_version`** — bumped on **any** edit to samples, features, or design-relevant metadata.
- **`history`** — an ordered action log (one entry per operation, with parameters); renders the reproducibility export.
- **`undo_stack`** — a short snapshot stack (depth 5) of prior `working` objects for `state_undo()`.
- **`meta`** — app-level flags not stored on the `dds` itself: `data_type` (bulk/single-cell), `feature_type` (drives `<feature_type>_name` labels and the status bar), `sce_per_cell`.
- **`derived`** — cached heavy artifacts (size factors, VST, PCA, DESeq fit, DE result tables), keyed by name. **A plain environment, not a reactive field** (writing cache entries must not trigger reactivity); staleness is keyed on `data_version` via `state_derive()`.

**Invalidation (conservative v1):** each `derived` artifact records the `data_version` it was computed under and is **stale** when that differs from the current `data_version`. The UI shows a stale badge + re-run button and **never silently reuses stale results**. Dependency order: `raw counts → filter → normalized assays → variable genes → PCA → DESeq fit → DE results → heatmap`. **Reset-to-original** restores `working <- original` + clears `derived`; **undo** keeps a short snapshot stack. These two are surfaced as **global Undo / Reset buttons in the status bar** (`state_undo()` / `state_reset()`); every committed edit (metadata Save, assay add, sample/feature removal) is undoable. **Scoped, local resets** complement them: the metadata editor's "Reset to original" reverts just that slot (`colData`/`rowData`) and **auto-commits** (no Save; it also shows an "unsaved changes" badge by Save), and the QC Filtering pills' "Reset Sample/Feature Removal" re-add removed items from `original` while keeping other edits (`reset_metadata_slot()`, `restore_samples()`/`restore_features()`). See [dev_ref/object-states-and-resets.md](dev_ref/object-states-and-resets.md). A **persistent dataset-status bar** on every page surfaces data type, n samples/features, assays, design-set?, and DE current/stale. The `history` log renders the **reproducibility export** (R script / Quarto report). Modules still follow the `dds`-in/`dds`-out contract (see `shiny-module`), but the canonical store is this state object.

Pages (navbar order). The first navbar entry, **Input**, is itself a sub-tabbed group (`navset_card_tab`) holding Dataset / Sample / Feature / Assay; the rest are top-level pages, plus a persistent status-bar `nav_item`:

1. **Input** (sub-tabbed group):
   - **Dataset** — load `dds`/`sce` **or** counts matrix + sample sheet; show `show(dds)`.
   - **Sample** — a `navset_card_pill` of *Sample Metadata* (edit `colData` via editable `DT`) and *Additional Metadata* (upload + bind a sample sheet). Built on the shared draft editor (`mod_meta_editor`).
   - **Feature** — a `navset_card_pill` of *Feature Metadata* (edit `rowData`, set `feature_class`/feature unit), *OrgDb Annotation*, and *GTF Annotation*. Annotation (OrgDb organism-first, then GTF which overrides) composes onto the editor draft and commits on Save; detect/confirm feature type; tag exogenous / ERCC spike-in. See the **`annotation`** skill.
   - **Assay** — add normalized assays (CPM/TPM/FPKM); size factors on **endogenous** via `controlGenes`; store `feature_length`.
2. **QC & filtering** — one **unified** page (a `navset_card_tab`) carrying a **data-type badge**. Sub-tabs: **Dataset diagnostics** — the VST mean–SD plot (rank-mean vs sd, points coloured by local density à la `geom_pointdensity` but computed with base `MASS::kde2d` to avoid an extra dep; VST via `qc_vst()`, which falls back from `DESeq2::vst()` to `varianceStabilizingTransformation()` on small/low-count data). **Sample QC** (pills): *General QC* (per-sample metrics — library size, detected, % mito, % spike-in; `qc_helpers.R`, `scater::perCellQCMetrics` + base-R fallback), *RLE*, *Expression density*, *QC Matrix* (the `dt_table`); each ggplot pill has its own group/colour control. **Sample Correlation** (pills): *Heatmap* (`ComplexHeatmap`, log-counts, Spearman/Pearson shown in legend+title, multi-column annotation with a Clear-annotation button) and *Within-group correlation* (`qc_within_group_correlation()` — per-sample mean correlation to same-group samples; candidate-outlier read). **Spike-in (ERCC)** (`ercc_helpers.R`) — **titration QC** (*not* normalization; size factors stay endogenous-`controlGenes`): a *Dose-response* log–log scatter of known concentration vs observed expression with a per-sample fit, and a *Per-sample summary* (% spike-in, detected-spike count, slope, R², lowest detected conc + spike-fraction CV). Concentration source is selectable — a `spike_concentration` rowData column (when present) or the bundled **ERCC Mix 1/2** (`ercc_concentrations()`, join by id); observed assay is a linear depth-normalized one (**TPM > FPKM > CPM**-preferred since length-normalized abundance tracks molar spike-in input; never counts/logcounts); zeros dropped before the log, fit NA-safe below 3 points (`spike_dose_response()`). All sample diagnostics compute on **endogenous** features only. **Diagnostic plots carry a short "How to read this" guideline below them** (`.qc_diag_help`). Heavy outputs show `shinycssloaders` spinners; VST/correlation/load show `withProgress`. **Filtering** (pills *Samples* / *Features*, `filter_helpers.R`): the app *suggests* low-quality items with per-reason detail (`flag_samples()` — opt-in library-size/detected/%-mito + within-group-correlation-outlier checks, each **blank threshold = disabled**, with an **Auto** button that fills robust median±3·MAD fences via `suggest_sample_thresholds()`; `flag_features()` — all-zero + `edgeR::filterByExpr()` or a manual `rowSums`/detected rule); **expression filters run on `endogenous` only — `spike_in`/`exogenous` exempt** (consistent with `controlGenes`); **avoid `HTSFilter`**. The user builds a **removal pool** (read-only *Suggested* column + separate *In pool* column; DT row-selection stages, buttons commit — Add/Remove selected, **Add all in current view**, **Adopt suggestions** as a union, Clear) and **Applies** it via `state_mutate` (a real, undoable removal → `drop_features()`/`drop_samples()`, which `refresh_assays()` and re-estimate size factors). **Features pre-seed the pool with the suggestion; samples are highlight-only / opt-in** (rnaseq-bioc: small-n outliers are flagged, not auto-dropped). The Features pill shows the **before/after log-CPM filtering density** (`qc_filter_density()`); Sample QC plots offer a **"Removal status" colour-by** (reason-aware green/yellow/red via `removal_status()`). A **"Showing:"** control lives in every sample-plot sidebar (Show-by column → *Keep* values, blank = show all) — all instances stay **synced** to one canonical selection (`show_by_rv`/`show_values_rv`); it hides samples from the sample plots + correlation heatmap **without** mutating the `dds` (view-only — diagnostics still compute on all samples; no `data_version` bump) and re-renders only on the tab's Render/auto (folded into each deferred spec). DT selected-row highlight uses the theme primary colour app-wide (`inst/www/custom.scss`).
3. **Dimensionality reduction** — **PCA-focused** (t-SNE/UMAP gated by sample count). Top-variable genes (default 500, assay default `logcounts`), excluding non-endogenous. Up to 4 panels (1/2/4); compute the embedding once, vary aesthetics per panel; color/shape by metadata or gene expression (via `lookup_feature()`).
4. **Differential expression** — guided design builder (+ reference level, full-rank check) and contrast picker; multiple stored contrasts; `DESeq2::DESeq()` + `lfcShrink`. Results carry **both** `log2FoldChange`/`log2FoldChange_shrunk` and matching `sig`/`DEG`/`sig_shrunk`/`DEG_shrunk`; shrinkage toggle = column selection. MA / volcano / direct-comparison plots with axis clamping (triangle markers).
5. **Expression heatmap** — `ComplexHeatmap` over a gene set (default = DEGs, fallback top-variable). **Default per-gene z-score** (toggle to raw `log10(TPM/CPM + 0.01)`). Row names hidden by default with `anno_mark()` for genes of interest; column names auto-hidden when > 30 samples.
6. **Palette** (`mod_palette.R`) — project-wide colour configuration; the **single source of truth** for discrete colour mappings. The user opts in per metadata column (empty by default), chooses a **base palette** from one `selectInput` with `<optgroup>`s (Custom / Qualitative / Brewer: Qualitative / Brewer: Sequential / Brewer: Divergent / viridis; option labels drop the `<pkg>:` prefix, values stay resolvable), and may hand-edit individual level colours (which flips the palette to **Custom palette**, so the selector never implies the named palette was modified). Per-column config is `list(name, colors)` where `colors` is the full level→hex map (everything else is inferred from `name`); it lives in `state$palette` (a UI preference, untouched by load/reset, no `data_version`) and is resolved by **`palette_helpers.R`**: `palette_discrete(levels, colors, name, custom)` is the one resolver feeding **both** the ggplot group scales (via the QC builders' `palette =` arg + `.qc_group_scale()`) **and** `qc_annotation_colors(df, config)` for ComplexHeatmap — so plots and heatmap annotations agree. `palette_colors(name, n)` infers the source from the name and enumerates from `viridisLite`/`RColorBrewer` (optional Suggests; graceful Okabe-Ito fallback); `palette_choices()` builds the grouped `selectInput` choices. **All resolved colours are normalized to 6-digit hex** (`norm_color()`, which also accepts R names like `gray50` / CSS names) — important because `viridis()` returns 8-digit alpha hex and a colour-picker echoes 6-digit, so without normalization the equality check would mis-fire and spuriously flip the selector to Custom. The page follows the **3-level layout** (like Input/QC): a `navset_card_tab` of **Setting** / **Preview**. Setting holds the customizable pills (Sample / Feature / Assay / Other), each with its own `layout_sidebar` (the QC pattern). The panel machinery is **generalized over four domains** (`colData`/`rowData`/`assays`/`other`), each item **discrete *or* continuous** — colData/rowData by column type; assays continuous; `other` = a discrete removal-status map + a continuous correlation ramp. **Continuous** mappings are `list(name, min, max, custom)`: a sequential/divergent/viridis palette + anchors (a number, a `p<pct>` percentile, or blank = data range), resolved against the real data **at the consumer** via `palette_resolve_range()` → `palette_colorramp2()` (ComplexHeatmap) / `palette_gradientn()` (ggplot `scale_*_gradientn` pieces). Wired consumers: numeric `colData` annotations + the correlation-heatmap ramp (`qc_annotation_colors`/`.qc_correlation_heatmap` `cor_config`) and the QC removal-status map (`removal_palette()`); rowData/assay configs are stored for P4/P5. **Preview** is sidebar-free: catalogue swatches grouped (qualitative as discrete blocks, sequential/divergent/viridis as ramps) + its own Collapse all. (Collapse calls `accordion_panel_close()` with the *unnamespaced* id — `sendInputMessage` already namespaces.) Per-item observers read levels **live** and use `ignoreInit`; ids are `<prefix>_<domain>__<item>`. `thematic` still owns fonts/fg/bg and the *default* scale when an item has no config. **P3g-a/b** are done; **P3g-c** adds JSON config import/export + a custom continuous-ramp picker UI; **P3g-d** adds factor management (coerce/reorder levels). The legacy `themer_mode` gallery's Heatmap sub-tab is superseded by the Preview tab.
7. **Export** — `dds` object, DE tables (XLSX), plots, and the reproducibility script/report.

## Conventions that must hold

- **`rowData(dds)$feature_class`** is **always present** (factor: `endogenous` default / `spike_in` / `exogenous`). Spike-in/exogenous are **flagged, never dropped** — excluded from size-factor estimation (`controlGenes = which(feature_class == "endogenous")`) and variable-gene selection, but still in the object so their expression stays plottable.
- **`rowData(dds)$feature_length`** holds effective length (bases) when known; it may come from a *different* feature type than the rows (e.g. exonic length for gene rows) — never assume it equals the row's genomic span. GTF-derived length uses the direct `rtracklayer` union-exon route (not a `TxDb`). TPM/FPKM degrade to CPM when absent.
- **`logcounts`** assay is added automatically on load; **CPM** may be auto-added too. All count-derived assays are **recomputed whenever features change** (drives `data_version`).
- **DE result columns** (every results table carries both LFC variants):
  - `log2FoldChange` (standard) and `log2FoldChange_shrunk` (`lfcShrink`).
  - `sig`/`DEG` and `sig_shrunk`/`DEG_shrunk` — `sig` rule: `!is.na(padj) & padj < 0.05 & abs(<chosen LFC>) >= log2(2)`; `padj` is shared. `DEG` factor levels `up`/`down`/`no_change`. Thresholds user-adjustable; the shrinkage toggle selects which column set drives plots.
- **Feature lookup & adaptive labels:** `lookup_feature()` resolves a query against `rowData`, defaulting to `<feature_type>_name` (e.g. `gene_name`), falling back to `rownames`. The feature unit (gene/transcript/…) is detected and relabeled dynamically; prompt when unknown; let the user correct it.
- **Heading hierarchy in the UI:** `h1` = app title (navbar); `h2` = top-level `nav_panel` labels ("Input", "QC", …); `h3` = Input sub-tab labels ("Dataset", "Sample", "Feature", "Assay") and `card_header` titles on non-Input pages; `h4` = `card_header` titles inside Input sub-tabs and `navset_card_pill` titles. All are visually scaled with `fs-*` to fit their context. Tab/pill navigation labels (plain strings in `nav_panel`) are not headings.
- **Read-only tables render via `dt_table()`** (`R/utils_table.R`) — the standard DT wrapper giving per-column filters, a rows-per-page selector, a search box, and horizontal scroll. Don't hand-roll `DT::datatable()` options for display tables; extend via its `options=`/`...` passthrough. The editable metadata editor (`mod_meta_editor`) is the deliberate exception (it has its own editable/selection config).
- **Tooltips use `bslib::tooltip()`** (Bootstrap 5), never the native `title=` HTML attribute — so hover styling and delay stay consistent app-wide (e.g. the `mod_meta_editor` reset buttons, the Filtering Auto buttons).
- **Plot-visualization tabs embed the view-only "Showing:" subset control** (`R/mod_plot_subset.R`): `plot_subset_ui(ns, <suffix>)` in each plot sidebar + one `plot_subset_server(input, output, session, state, suffixes)` call returning a `showing_samples` reactive the plots filter their data by (and the deferred `sig`s depend on). It's display-only — never mutates the `dds` or bumps `data_version`. This is the standard for any page with per-sample plots (QC today; PCA/DE later).
- **Deferred rendering:** plots can be slow. Gate re-rendering behind an explicit "apply / render" button (`bindEvent()` / `eventReactive()`), not live reactivity on every input.
- **ggplot↔plotly engine toggle** (`mod_qc.R`): a global status-bar switch sets `state$plot_interactive` (default off = static ggplot). Toggleable ggplot plots use the **dual-output** pattern — a `gg(interactive)` builder fn feeds `dual_plot("<id>", gg, n_elements)`, which wires a static `plotOutput` (always `gg(FALSE)`) + an interactive `plotlyOutput` (`ggplotly`, `gg(TRUE)`) behind a `uiOutput` container; the UI placeholder is `.qc_dual_plot(ns("<id>_container"))`. Interactivity needs `use_plotly_base()` (toggle ON **and** `requireNamespace("plotly")` — a `Suggests`, graceful fallback) and is then gated **per plot on an element budget**: `n_elements()` (≈ rows of the plotted data, the real `ggplotly` cost driver) vs `.plotly_max_elements()` = `getOption("ddsdashboard.plotly_max_elements", 5000L)` — tunable without a settings page. Over budget → static + a **sticky per-plot "Render interactive anyway"** override (`forced` reactiveVal, reset on `data_version` change or toggle flip). Hover comes from a `text` aes the builders add only when `interactive = TRUE`; the plotly path muffles the residual aesthetic warning (`.muffle_unknown_aes`) and catches `validate()`/`req()` into a message figure. **Not** for ComplexHeatmap / non-ggplot outputs. Reuse for future plot pages (PCA/DE).
- **Caching is the `derived` store** (see State model): VST/PCA/DESeq fit are version-stamped and invalidated when `data_version` changes. `bindCache()` or keyed entries both work.

## Workflow

- Build the app in **phases** (bulk-first). **The detailed phase plan + current status lives in [dev_ref/roadmap.md](dev_ref/roadmap.md) — consult/update it as the source of truth** (don't duplicate the full list here). At a glance: P1 (skeleton), P2 (annotation + normalization + theming), and **P3a–f** (QC, filtering, edit-history controls, ERCC dose-response, spike-in QC filtering, ggplot↔plotly engine toggle) plus **P3g-a/b** (Palette page — discrete + continuous, all four domain pills) are **done**; **P3g-c (Palette config JSON import/export + custom-ramp UI) is next**, then P3g-d (factor management), P4 (PCA), P5 (DESeq2 + DE plots + heatmap + add the `de-analysis` skill), P6 (Export & reproducibility: R script/Quarto report, plot export, DE tables), P7 (single-cell, lowest priority).
- **Build mock-`dds` fixtures early** (`data-raw/`, used by `tests/`) so every phase has data to test against.
- **Commit after each meaningful, self-contained change** with a descriptive message.
- Heavy Bioconductor objects: don't print full matrices to logs; prefer `show()` / dimensions / `head()`.
