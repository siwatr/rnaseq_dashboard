---
name: shiny-de-analysis
description: How the DESeq2 differential-expression page in this RNA-seq dashboard runs and plots DE. Covers the guided design/contrast builder (reference level + full-rank check), running DESeq2 + LFC shrinkage (apeglm/ashr/normal fallback), the dual-LFC sig/DEG(_shrunk) result schema, the size-factor edge cases (poscounts fallback), and the DE-plot conventions (MA / volcano / direct-comparison; per-feature colour resolved locally but configured via the Palette other/DEG slot; axis clamping -> triangle markers; layered structure for the P6 gene-set overlay; ggrepel labels). Engine implemented as pure helpers in R/de_helpers.R; the page is mod_de + the shared mod_design_builder. Use when building or changing the DE page, the shared design builder, or any DE plot/result logic. Pair with rnaseq-bioc (the statistical conventions), shiny-module, and shiny-plot-aesthetics.
metadata:
  type: project
  version: "1.0"
---

# Differential expression (DESeq2) — the DE module

The DE page (`mod_de`) fits DESeq2, stores contrasts, and draws MA / volcano /
direct-comparison plots + a results table. The **statistical rules live in
`rnaseq-bioc`**; this skill covers **how the shiny app implements and plots them**
via the pure helpers in [R/de_helpers.R](../../../R/de_helpers.R). Design is
committed to the dds (shared with the Input **Design** tab via `mod_design_builder`).
The **fit** and the per-contrast **result extraction are separate**: Run DESeq2 fits
(contrast-free), stored in the `derived` env under a `(data_version, design_version)`
stamp; extraction is reactive (Auto-update) or on demand, cached per contrast in
`state$de$results`. The page is a `navset_card_tab`: **Design & Contrasts** / **DE
Plots** / **Results Table** (all wired; P5a–c complete).

## Design + contrast (guided, not free-text)

- **Candidate factors:** `de_design_factors(dds)` (discrete colData columns, >= 2 levels).
- **Reference level:** `de_relevel(dds, col, ref)` — the DE *primitive* for the control
  level (Phase 8 factor management generalizes reordering and reuses it). Committed to the
  dds so `apeglm` can shrink by coefficient.
- **Full-rank guard before fitting:** `de_full_rank(design, colData)` → `list(ok, rank,
  ncoef, msg)`; block the Run button + show `msg` when `!ok` (confounded/collinear terms).
- **Contrast picker:** `de_contrast_levels(dds, var)` → level choices; a contrast is
  `c(variable, test, control)`; support **multiple stored contrasts** (a module `reactiveVal`).
- **Shared builder:** `mod_design_builder` is embedded in *both* the Input **Design** tab and the
  DE *Design & Contrasts* tab; both read/write `design(working)` so they stay in sync. Design +
  relevel are **design-scoped** state edits — they invalidate the DE fit but **not** `data_version`
  (PCA/VST/QC don't depend on the model).

## Fit + shrinkage + the dual-LFC schema

- `de_run(dds)` — ensures robust size factors (see below) then `DESeq2::DESeq()`. **Rerun from
  raw counts whenever the data or design changes** (cache keyed on `data_version` + design signature).
- `de_results(fit, contrast, shrink_type)` → a `data.frame` with the DESeq2 columns **plus
  `log2FoldChange_shrunk`** (`de_shrink()`).
- **`de_shrink()` fallback chain:** `apeglm` needs a **coefficient** (`de_coef_name(contrast)`,
  available when control is the reference) → else **`ashr`** (takes the contrast) → else built-in
  **`type="normal"`** → else `NA`. All probed with `requireNamespace()` (apeglm/ashr are `Suggests`).
- **Classification is cheap + threshold-driven:** `de_classify_table(df, padj_cut, lfc_cut)` adds
  `sig`/`DEG` (+ `sig_shrunk`/`DEG_shrunk`). `classify()` rule: `!is.na(padj) & padj < cut &
  abs(lfc) >= lfc_cut`; `DEG` factor levels **`up`/`down`/`no_change`**. **Precompute both LFC
  variants; the shrinkage toggle is column selection (no refit); thresholds re-derive sig/DEG (no refit)**
  — so padj/|LFC| live in the plot/table controls, *not* the fit cache key.
- `de_summary(df, deg_col)` → `up`/`down`/`total` for the per-contrast summary.

## Size-factor edge cases (keep an eye on these)

`estimate_size_factors_endogenous()` ([R/assay_helpers.R](../../../R/assay_helpers.R)) uses
median-of-ratios on endogenous `controlGenes`. Known failure modes:

- **Every (control) gene has >= 1 zero** → all geometric means are 0 → *"every gene contains at
  least one zero, cannot compute log geometric means"*. **Handled:** the helper catches this and
  retries **`type="poscounts"`** (geometric mean over positive counts).
- **An all-zero / near-empty sample** → degenerate size factor; **a tiny gene set** after
  aggressive filtering → unstable. Surface, don't silently proceed.
- `DESeq()` estimates size factors internally when unset — so pre-setting robust ones (via the
  hardened helper, which `de_run()` calls) injects the poscounts fallback into the whole pipeline.

## DE Plots + Results Table (the P5c tabs)

- **Plot engine:** the DE Plots tab reuses `plot_engine_server` (`dual_plot` / `deferred` /
  `stale_note`, [R/mod_plot_engine.R]) exactly like PCA. A **segmented control** (MA | Volcano |
  Direct) drives **one** deferred plot; `n_elements = nrow(results)` gates the plotly fallback.
  **No `mod_plot_subset`** — DE plots are per-feature, not per-sample.
- **Deferred sig = the DATA only** (contrast + padj + |LFC| + `data_version`); the spec is the
  classified df (`de_classify_table`, carrying both `DEG`/`DEG_shrunk`). **Display aesthetics stay
  live** (plot type, Use-shrunk, colour-by, labels, clamps, point/alpha/legend) — cheap re-plots,
  never in the sig (the shiny-plot-aesthetics rule).
- **Colour is per-feature, resolved locally** — `de_colour_resolve(df, key, colors)` (the shared
  `aes_helpers` resolver is per-*sample*). `key` = `"__none__"` / `"DEG"`/`"DEG_shrunk"` / a numeric
  DE column. **Config is centralized:** the DEG palette comes from the Palette **`other/DEG`** slot
  — a curated **"DEG palette" set** (Pink-Blue default / Orange-Purple / Red-Blue / Coral-Teal),
  led by `deg_palette_choices()` for the 3-level DEG item, which **also offers the generic discrete
  palettes** (Okabe-Ito / Brewer / viridis) so DEG can be coloured with any scheme; the DEG set is
  resolved by a `"DEG:"` branch in `palette_colors()` (kept out of `palette_type_names()` so other
  items never see it; still shown in the Palette Preview). Read it with
  `palette_discrete(c("up","down","no_change"), cfg$colors, cfg$name %||% "DEG: Pink-Blue", cfg$custom)`.
- **Point draw order:** the shared scatter reorders rows so a discrete DEG colour draws
  **`no_change` first** (up/down land on top); continuous colours draw low-first. So DEGs are never
  hidden behind the grey background — the standard "arrange before plotting" idiom, done in `.de_scatter`.
- **Axis limits are field-based AND set real coord limits.** Clamp/limit by *field* (log2FC /
  -log10(padj) / log10(baseMean) / expression), not raw X/Y, so a limit follows its field across plot
  types (log2FC is MA-y and volcano-x). `.de_scatter` (a) `de_clamp`s out-of-range points → **shape**
  (circle vs triangle), never dropped, **and** (b) applies `coord_cartesian(xlim/ylim)` (or
  `coord_fixed` when squared) so a range **wider than the data extends the axis** (e.g. a symmetric
  log2FC window for a figure), not only pulls outliers in. Gene labels are clamped to the same range
  so they sit at the boundary too. The Direct plot's **1:1 aspect** toggle passes `fixed_ratio` into
  `de_direct_gg` (one coord, no double-`coord_*`).
- **Builders:** `de_ma_gg` / `de_volcano_gg` / `de_direct_gg(de_group_means(...))` — each takes
  `interactive` (hover `text` aes for plotly only), `point_size`/`point_alpha`, `dark` (high-contrast
  bold ggrepel labels with a `bg.color` halo), and an optional ggrepel `labels` frame (top-N by padj +
  ad-hoc searched genes via `resolve_feature`); `de_direct_gg` also takes `fixed_ratio`.
- **Direct plot expression value = the shared `expr_value_*` control** ([R/mod_expr_value.R] —
  `expr_value_ui(ns, suffix)` + `expr_value_server(input, output, session, state, suffix)`, the
  host-namespace idiom like `mod_plot_subset`): **assay + transform (none/log2/log10) + pseudocount**,
  applied via `de_group_means(dds, assay, ctrl, test, transform, pseudocount)` (transform → average).
  This is the **standard expression-value control** the P7 Expression page reuses; the axis label
  comes from `expr_value_label()`.
- **Shared "Contrast & thresholds" group in *all three* tabs** (Design & Contrasts / DE Plots /
  Results Table): the `.de_view_controls(ns, suffix)` accordion panel (contrast-to-view + padj + |LFC|
  + Use-shrunk), rendered once per tab with a per-tab `suffix` (ids can't repeat). The server syncs
  the copies — contrast via `state$de$active`, thresholds via a canonical `thr` reactiveValues — with
  **guarded fan-out** (no loop). Downstream reads `thr$*`, not the raw inputs. The **Render / Auto-render**
  controls sit **above the plot** in the card (not buried in the sidebar). The **Results Table** is a
  `dt_table()` of the active contrast with DEG `formatStyle`/`styleEqual` colouring (same `other/DEG`
  palette) + a significant-only filter.
- **Fit vs results status:** the note above **Run DESeq2** reflects the **fit** (`de_fit_status(state)`:
  none / stale / current) — a fit with no contrasts yet still reads "up to date". `de_status(state)`
  (results-level) drives the DEG-summary staleness note.
- **Layered for the P6 gene-set overlay (channel separation):** reserve **shape** for the clamp
  (triangle), **fill/stroke** for gene-set DEG-membership (filled vs hollow "donut", `shape 21`
  `fill=NA`), **colour** for DEG-status-or-set. So the clamp and the donut never collide. The
  gene-set overlay itself (select a saved set to highlight) **depends on `state$gene_sets`
  (Phase 6)** — a P6 follow-on; P5's layered builders accept it without rework.

## Checklist

- [ ] Design goes through `de_full_rank` before Run; reference set via `de_relevel` (committed to the dds).
- [ ] Fit cached on `data_version` + design signature; **rerun from raw counts** on data/design change.
- [ ] Results carry both `log2FoldChange(_shrunk)`; shrinkage toggle = column selection; thresholds re-derive sig/DEG (no refit).
- [ ] `de_shrink` degrades apeglm → ashr → normal → NA via `requireNamespace()`.
- [ ] Size factors go through the hardened endogenous estimator (poscounts fallback).
- [ ] DE colour via `de_colour_resolve` (local), DEG palette from the curated `other/DEG` set; clamp = shape, never drop.
- [ ] DE Plots deferred sig = data-only; display aesthetics (plot type / colour / labels / clamps) stay live.
- [ ] Shared "Contrast & thresholds" group in all 3 tabs; synced via `state$de$active` + a canonical `thr` (guarded fan-out, no loop).
- [ ] Field-based axis limits set real coord limits (extend the axis); labels clamped; Direct plot `fixed_ratio` 1:1.
- [ ] DEG points draw `no_change` first (DEGs on top); labels are high-contrast (`dark` + halo).
- [ ] Direct plot uses the shared `expr_value_*` control (assay/transform/pseudocount) → `de_group_means`.
- [ ] Run-note uses `de_fit_status` (fit), summary uses `de_status` (results).
- [ ] Plots are layered so the P6 gene-set overlay drops in without rework.
