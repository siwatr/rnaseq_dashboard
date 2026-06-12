# CLAUDE.md

Guidance for working in this repository. See [rough_design.md](rough_design.md) for the full design narrative — this file is the distilled, load-bearing version.

## What this is

An R **Shiny** dashboard that lets people with little to no bioinformatics background **explore, manipulate, and visualize bulk RNA-seq results**. The input is an R object — a `DESeqDataSet` (`dds`) or a `SingleCellExperiment` (`sce`). An `sce` is converted to a `dds` on load, because differential expression (DESeq2) is a core feature.

The minimum viable input is a single `dds` carrying only a raw count matrix; everything else (sample metadata, gene metadata, feature lengths, normalized assays, DE results) can be added inside the app.

## Tech stack

- **R + Shiny**, UI built with **bslib** (Bootstrap 5). Use the `shiny-bslib` and `shiny-bslib-theming` skills for layout/theming. Prefer `page_navbar()` for the multi-page structure.
- **Bioconductor core:** `DESeq2`, `SummarizedExperiment`, `SingleCellExperiment`, `scater` (QC), `rtracklayer` (GTF import).
- **Plots:** `ggplot2` is the default for all plots. `ComplexHeatmap` is the heatmap engine (do not substitute pheatmap/heatmaply).
- **Tables / IO:** reading user metadata from `.xlsx`, `.tsv`, `.csv`.

Domain conventions (object shape, assays, normalization math, QC metrics, DE result columns) live in the **`rnaseq-bioc`** skill — consult it before touching data logic.

## Project layout

This is a **plain R package** named `ddsdashboard` (not golem). Reusable helpers live in `R/` and are exported so other lab projects can `library(ddsdashboard)` and call them; the Shiny app is just `run_app()`.

```
DESCRIPTION / NAMESPACE   # package metadata + exports (regenerate NAMESPACE with devtools::document())
environment.yml           # conda/mamba env (source of truth for deps)
R/
  run_app.R               # run_app(), app_ui(), app_server() — the navbar shell + module wiring
  mod_<page>.R            # one Shiny module per page (input, qc, process, dimreduc, de, heatmap, export)
  utils_normalization.R   # cpm(), tpm(), fpkm(), logcounts_from_counts() — pure, exported, tested
  utils_lookup.R          # lookup_feature() — gene id <-> name resolution
tests/testthat/           # unit tests for the utils
```

Module contract and the `dds`-in/`dds`-out wiring are documented in the **`shiny-module`** skill. New pages follow that pattern.

**Dependency policy:** only what the code actually calls today sits in `DESCRIPTION` `Imports` (`shiny`, `bslib`). Bioconductor packages (`DESeq2`, `scater`, `ComplexHeatmap`, …) are in `Suggests` and `environment.yml`; promote a package to `Imports` when a feature starts calling it. The pure base-R utilities (`utils_normalization.R`) deliberately take no dependencies so they stay reusable and fast to test.

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

The shared state is **one reactive `dds` object** threaded through every module (a `reactiveVal` or single-field `reactiveValues`). Modules read it and return an updated `dds`; the top-level server owns the canonical copy. Any edit to samples or features **invalidates downstream results** — most importantly, removing a sample/feature requires re-running `DESeq2::DESeq()`.

Pages, in order:

1. **Input data** — load `dds`/`sce`; show `show(dds)`; edit `colData` (sample metadata) directly or via uploaded table; map `rowData`/`rowRanges` from an uploaded GTF (feature type, name column, selected `mcols`, effective feature length); tag exogenous / ERCC spike-in genes.
2. **QC & filtering** — `scater` QC metric bar/box plots (library size, detected features, % mito, % spike-in), variance-stabilization plot, sample–sample correlation heatmap (log-counts). Suggest samples to drop (with reasons) and let the user override. Feature filtering: drop all-zero genes always; drop low-count genes by a user threshold on `rowSums()`. **Avoid `HTSFilter`** (too slow).
3. **Processing** — add normalized assays (TPM, FPKM, CPM).
4. **Dimensionality reduction** — PCA (DESeq2 or manual), t-SNE, UMAP from top-variable genes (default 500, assay default `logcounts`). Up to 4 plot panels (layouts: 1 / 2 / 4), color/shape by metadata (discrete or continuous) or by gene expression. One panel may be a box/violin plot of a gene.
5. **Differential expression** — design the formula, run `DESeq2::DESeq()`, set thresholds; build MA / volcano / direct-comparison scatter plots.
6. **Expression heatmap** — `ComplexHeatmap` over a gene set (default = DEGs); expression default `log10(TPM + 0.01)`, fall back to `log10(CPM + 0.01)`.
7. **Export** — `dds` object, DE result tables, plots.

## Conventions that must hold

- **`logcounts`** assay is added automatically on load. **`CPM`** may be auto-added too, but both must be **recomputed whenever features are added or removed**.
- **DE result columns** added to every results table:
  - `sig` (logical): `TRUE` when significant. Default rule: `!is.na(padj) & padj < 0.05 & abs(log2FoldChange) >= log2(2)`.
  - `DEG` (factor, levels `up` / `down` / `no_change`).
  - Default thresholds: `padj < 0.05`, `abs(log2FoldChange) >= log2(2)` — both user-adjustable.
- **Gene ID → human-readable name lookup:** feature IDs (e.g. Ensembl gene IDs) are not readable. Provide a helper that resolves a query against `rowData(dds)`, defaulting to the `<feature_type>_name` column (e.g. `gene_name` when feature type is `gene`), falling back to `rownames(dds)`.
- **Effective feature length** for TPM/FPKM may come from a *different* feature type than the rows (e.g. exonic length for gene-level rows). Keep length as an explicit, separately-sourced value — never assume it equals the row feature's span.
- **Deferred rendering:** plots can be slow. Gate re-rendering behind an explicit "apply / render" button (`bindEvent()` / `eventReactive()`), not live reactivity on every input.
- **Cache expensive computations** — VST and PCA (and the `DESeq()` fit) are reused across panels and re-renders. Cache them keyed on the current `dds` state plus the relevant parameters (assay, n top-variable genes, formula), and **invalidate the cache whenever samples or features change** (same trigger as the stale-state rule above). `bindCache()` or a keyed entry in the shared state both work.

## Workflow

- Build the app in **phases** (see roadmap in [rough_design.md](rough_design.md)): Phase 1 = layout + import + export shell; Phase 2 = metadata edit, QC, filtering, normalization; Phase 3 = dim-reduction page; Phase 4 = DESeq2 + DE plots.
- **Commit after each meaningful, self-contained change** with a descriptive message.
- Heavy Bioconductor objects: don't print full matrices to logs; prefer `show()` / dimensions / `head()`.
