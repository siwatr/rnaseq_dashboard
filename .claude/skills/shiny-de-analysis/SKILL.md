---
name: shiny-de-analysis
description: How the DESeq2 differential-expression page in this RNA-seq dashboard runs and plots DE. Covers the guided design/contrast builder (reference level + full-rank check), running DESeq2 + LFC shrinkage (apeglm/ashr/normal fallback), the dual-LFC sig/DEG(_shrunk) result schema, the size-factor edge cases (poscounts fallback), and the DE-plot conventions (MA / volcano / direct-comparison; per-feature colour resolved locally but configured via the Palette other/DEG slot; axis clamping -> triangle markers; layered structure for the P6 gene-set overlay; ggrepel labels). Engine implemented as pure helpers in R/de_helpers.R; the page is mod_de + the shared mod_design_builder. Use when building or changing the DE page, the shared design builder, or any DE plot/result logic. Pair with rnaseq-bioc (the statistical conventions), shiny-module, and shiny-plot-aesthetics.
metadata:
  type: project
  version: "0.1"
---

# Differential expression (DESeq2) — the DE module

The DE page (`mod_de`) fits DESeq2, stores contrasts, and draws MA / volcano /
direct-comparison plots + a results table. The **statistical rules live in
`rnaseq-bioc`**; this skill covers **how the shiny app implements and plots them**
via the pure helpers in [R/de_helpers.R](../../../R/de_helpers.R). Design is
committed to the dds (shared with the Input **Design** tab via `mod_design_builder`);
the fit + results are cached via `state_derive`.

> Scaffold (P5a). The engine + these conventions land first; the module wiring
> (P5b design/contrast builder + fit; P5c plots/table) fills in against this spec.
> Finalize this skill at P5c once the UI reality is settled.

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

## DE plots (per-feature; layered; customizable)

- **Colour is per-feature, so it's resolved locally** — `de_colour_resolve(df, key, colors)` (the
  shared `aes_helpers` resolver is per-*sample*). `key` = `"__none__"` / `"DEG"`/`"DEG_shrunk"`
  (discrete) / a numeric DE column (continuous). **The colour config is still centralized:** the
  DEG palette comes from the Palette **`other/DEG`** slot (P5c) — pass its colours in as `colors`.
- **Axis clamping → triangles:** `de_clamp(x, lo, hi)` pulls out-of-range points to the limit and
  flags them; the builders map the flag to **shape** (circle vs triangle), never dropping points.
- **Builders:** `de_ma_gg()` (x = `log10(baseMean)`, y = LFC), `de_volcano_gg()` (x = LFC, y =
  `-log10(padj)`), `de_direct_gg(mean_df)` (x = control mean, y = test mean; `de_group_means()`
  computes the means + a y=x guide). Each takes `interactive` (adds a hover `text` aes for the
  plotly path only) and an optional ggrepel `labels` frame — reuse the shared `dual_plot` engine.
- **Layered for the P6 gene-set overlay (channel separation):** reserve **shape** for the clamp
  (triangle), **fill/stroke** for gene-set DEG-membership (filled vs hollow "donut", `shape 21`
  `fill=NA`), **colour** for DEG-status-or-set. So the clamp and the donut never collide. The
  gene-set overlay itself (select a saved set to highlight) **depends on `state$gene_sets`
  (Phase 6)** — P5 ships colour + labels (top-N + ad-hoc searched genes) and the layer structure.

## Checklist

- [ ] Design goes through `de_full_rank` before Run; reference set via `de_relevel` (committed to the dds).
- [ ] Fit cached on `data_version` + design signature; **rerun from raw counts** on data/design change.
- [ ] Results carry both `log2FoldChange(_shrunk)`; shrinkage toggle = column selection; thresholds re-derive sig/DEG (no refit).
- [ ] `de_shrink` degrades apeglm → ashr → normal → NA via `requireNamespace()`.
- [ ] Size factors go through the hardened endogenous estimator (poscounts fallback).
- [ ] DE colour via `de_colour_resolve` (local), DEG palette from `other/DEG`; clamp = shape, never drop.
- [ ] Plots are layered so the P6 gene-set overlay drops in without rework.
