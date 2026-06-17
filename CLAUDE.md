# CLAUDE.md

Guidance for working in this repository. See [rough_design.md](rough_design.md) for the full design narrative — this file is the distilled, load-bearing version.

## What this is

An R **Shiny** dashboard that lets people with little to no bioinformatics background **explore, manipulate, and visualize RNA-seq results**. The input is an R object — a `DESeqDataSet` (`dds`) or a `SingleCellExperiment` (`sce`) — or a **counts matrix + sample sheet** (CSV/TSV/XLSX) that the app assembles into a `dds`.

**Scope is bulk-first** (decided in review). Build the bulk path end-to-end first; keep the design general enough to add single-cell later. An `sce` is supported via two routes: **pseudobulk** aggregation by a sample grouping, or **per-cell coercion to `dds`** behind a *"statistically inaccurate and slow"* warning. **Force pseudobulk above ~1k cells.** (Why per-cell DESeq2 is unsound: [dev_ref/normalization-scran-vs-deseq.md](dev_ref/normalization-scran-vs-deseq.md).)

The minimum viable input is a `dds` carrying only a raw count matrix; everything else (sample metadata, gene metadata, feature lengths, normalized assays, DE results) can be added inside the app.

## Tech stack

- **R + Shiny**, UI built with **bslib** (Bootstrap 5). Use the `shiny-bslib` and `shiny-bslib-theming` skills for layout/theming. `page_navbar()` for the multi-page structure; `DT` for editable tables; `shinycssloaders` for spinners; `thematic` so ggplots match the theme.
- **Bioconductor:** `DESeq2`, `SummarizedExperiment`, `SingleCellExperiment`, `scater` (QC), `rtracklayer` (GTF), `edgeR` (`filterByExpr`), `apeglm`/`ashr` (LFC shrinkage), `vsn` (`meanSdPlot`), `AnnotationDbi` + `org.*.eg.db` (annotation), `scran` (single-cell, later).
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
  mock_data.R             # make_mock_dds() etc. — fixtures used by tests + the demo
data-raw/                 # scripts that build mock dds/sce fixtures for tests + demo
dev_ref/                  # developer reference notes (e.g. normalization-scran-vs-deseq.md)
tests/testthat/           # unit tests for helpers + module servers (shiny::testServer)
inst/extdata/             # bundled reference data (e.g. ERCC ids + concentrations)
```

The file conventions: **`mod_<page>.R`** = a navbar page; **`mod_<shared>.R`** = a reusable sub-module embedded inside pages (e.g. the draft metadata editor, the GTF reader); **`<topic>_helpers.R`** = pure, Shiny-free, exported functions that the modules call (so they stay reusable and unit-testable).

Module contract and the shared-`state` wiring are documented in the **`shiny-module`** skill (which also covers the reusable sub-module + draft-editor patterns). New pages follow that pattern. Domain/statistical conventions (object shape, assays, normalization, QC, DE schema, `feature_class`/`feature_length`, ERCC) live in **`rnaseq-bioc`**; the **`annotation`** skill covers the OrgDb + GTF row-annotation write conventions (consult it before touching mapping/annotation logic). A **`de-analysis`** skill (design/contrast/shrinkage) may be added when Phase 5 lands. The **`bioc-reviewer`** agent reviews any data-logic change; the **`shiny-feature-builder`** agent implements a page/helper end-to-end.

**Dependency policy:** `DESCRIPTION` `Imports` holds the always-on framework deps the app loads unconditionally — currently `AnnotationDbi`, `bslib`, `DESeq2`, `DT`, `methods`, `readr`, `readxl`, `shiny`, `SummarizedExperiment`. Heavy or optional packages are loaded **on demand via `requireNamespace()`** (with graceful degradation / a `skip_if_not_installed()` in tests) and stay in `Suggests` *even though code calls them* — currently `ComplexHeatmap`, `GenomicRanges`, `IRanges`, `S4Vectors`, `rtracklayer`, `org.*.eg.db`, `ps`, `scater`, `SingleCellExperiment`. They're also pinned in `environment.yml`. Promote `Suggests → Imports` only when a feature loads the package unconditionally. The pure base-R utilities (`utils_normalization.R`) deliberately take no dependencies so they stay reusable and fast to test.

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

**Invalidation (conservative v1):** each `derived` artifact records the `data_version` it was computed under and is **stale** when that differs from the current `data_version`. The UI shows a stale badge + re-run button and **never silently reuses stale results**. Dependency order: `raw counts → filter → normalized assays → variable genes → PCA → DESeq fit → DE results → heatmap`. **Reset-to-original** restores `working <- original` + clears `derived`; **undo** keeps a short snapshot stack. A **persistent dataset-status bar** on every page surfaces data type, n samples/features, assays, design-set?, and DE current/stale. The `history` log renders the **reproducibility export** (R script / Quarto report). Modules still follow the `dds`-in/`dds`-out contract (see `shiny-module`), but the canonical store is this state object.

Pages (navbar order). The first navbar entry, **Input**, is itself a sub-tabbed group (`navset_card_tab`) holding Load / Sample info / Feature info / Assay; the rest are top-level pages, plus a persistent status-bar `nav_item`:

1. **Input** (sub-tabbed group):
   - **Load** — load `dds`/`sce` **or** counts matrix + sample sheet; show `show(dds)`.
   - **Sample info** — a `navset_card_pill` of *Sample Metadata* (edit `colData` via editable `DT`) and *Additional Metadata* (upload + bind a sample sheet). Built on the shared draft editor (`mod_meta_editor`).
   - **Feature info** — a `navset_card_pill` of *Feature Metadata* (edit `rowData`, set `feature_class`/feature unit), *OrgDb Annotation*, and *GTF Annotation*. Annotation (OrgDb organism-first, then GTF which overrides) composes onto the editor draft and commits on Save; detect/confirm feature type; tag exogenous / ERCC spike-in. See the **`annotation`** skill.
   - **Assay** — add normalized assays (CPM/TPM/FPKM); size factors on **endogenous** via `controlGenes`; store `feature_length`.
2. **QC & filtering** — one **unified** page with a **data-type badge**; per-sample (bulk) vs per-cell (single-cell) metrics (library size, detected features, % mito, % spike-in); `vsn::meanSdPlot` VST plot; sample–sample correlation heatmap (log-counts). Gentle outlier flags (not aggressive auto-drop) for small-n bulk. Feature filtering: always drop all-zero; `edgeR::filterByExpr()` default + manual `rowSums` threshold. **Avoid `HTSFilter`** (too slow).
3. **Dimensionality reduction** — **PCA-focused** (t-SNE/UMAP gated by sample count). Top-variable genes (default 500, assay default `logcounts`), excluding non-endogenous. Up to 4 panels (1/2/4); compute the embedding once, vary aesthetics per panel; color/shape by metadata or gene expression (via `lookup_feature()`).
4. **Differential expression** — guided design builder (+ reference level, full-rank check) and contrast picker; multiple stored contrasts; `DESeq2::DESeq()` + `lfcShrink`. Results carry **both** `log2FoldChange`/`log2FoldChange_shrunk` and matching `sig`/`DEG`/`sig_shrunk`/`DEG_shrunk`; shrinkage toggle = column selection. MA / volcano / direct-comparison plots with axis clamping (triangle markers).
5. **Expression heatmap** — `ComplexHeatmap` over a gene set (default = DEGs, fallback top-variable). **Default per-gene z-score** (toggle to raw `log10(TPM/CPM + 0.01)`). Row names hidden by default with `anno_mark()` for genes of interest; column names auto-hidden when > 30 samples.
6. **Export** — `dds` object, DE tables (XLSX), plots, and the reproducibility script/report.

## Conventions that must hold

- **`rowData(dds)$feature_class`** is **always present** (factor: `endogenous` default / `spike_in` / `exogenous`). Spike-in/exogenous are **flagged, never dropped** — excluded from size-factor estimation (`controlGenes = which(feature_class == "endogenous")`) and variable-gene selection, but still in the object so their expression stays plottable.
- **`rowData(dds)$feature_length`** holds effective length (bases) when known; it may come from a *different* feature type than the rows (e.g. exonic length for gene rows) — never assume it equals the row's genomic span. GTF-derived length uses the direct `rtracklayer` union-exon route (not a `TxDb`). TPM/FPKM degrade to CPM when absent.
- **`logcounts`** assay is added automatically on load; **CPM** may be auto-added too. All count-derived assays are **recomputed whenever features change** (drives `data_version`).
- **DE result columns** (every results table carries both LFC variants):
  - `log2FoldChange` (standard) and `log2FoldChange_shrunk` (`lfcShrink`).
  - `sig`/`DEG` and `sig_shrunk`/`DEG_shrunk` — `sig` rule: `!is.na(padj) & padj < 0.05 & abs(<chosen LFC>) >= log2(2)`; `padj` is shared. `DEG` factor levels `up`/`down`/`no_change`. Thresholds user-adjustable; the shrinkage toggle selects which column set drives plots.
- **Feature lookup & adaptive labels:** `lookup_feature()` resolves a query against `rowData`, defaulting to `<feature_type>_name` (e.g. `gene_name`), falling back to `rownames`. The feature unit (gene/transcript/…) is detected and relabeled dynamically; prompt when unknown; let the user correct it.
- **Deferred rendering:** plots can be slow. Gate re-rendering behind an explicit "apply / render" button (`bindEvent()` / `eventReactive()`), not live reactivity on every input.
- **Caching is the `derived` store** (see State model): VST/PCA/DESeq fit are version-stamped and invalidated when `data_version` changes. `bindCache()` or keyed entries both work.

## Workflow

- Build the app in **phases** (bulk-first; see roadmap in [rough_design.md](rough_design.md)): P1 = layout + status bar + import (rds + tabular) + export shell *(done)*; P2 = metadata edit, annotation (OrgDb + GTF), normalization, GTF memory minimization + session memory monitor, feature-info annotation UI redesign *(done)*; **P3 = QC + sample/feature filtering (+ ERCC) ← next**; P4 = PCA dim-reduction; P5 = DESeq2 + DE plots + heatmap; P6 (later) = single-cell + pseudobulk. A theming pass (dark mode + `thematic`) is queued before the Phase-4 plots.
- **Build mock-`dds` fixtures early** (`data-raw/`, used by `tests/`) so every phase has data to test against.
- **Commit after each meaningful, self-contained change** with a descriptive message.
- Heavy Bioconductor objects: don't print full matrices to logs; prefer `show()` / dimensions / `head()`.
