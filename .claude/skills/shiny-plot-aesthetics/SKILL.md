---
name: shiny-plot-aesthetics
description: How any plot page in this RNA-seq dashboard turns a per-sample "colour / group / annotation by" attribute into values + colours, via the one shared catalog + resolver (R/aes_helpers.R). Covers the attribute model (colData / gene / QC metrics / removal / pool / spike-in), the grouped selector choices (group_field_choices / aes_choices), the resolver API (aes_resolve / aes_ggplot_scale / aes_heatmap_col / aes_annotation), the project standard that every attribute is palette-customizable (Palette "Other" + aes_other_palette_items), and the deferred values-snapshot / colours-live rule for heatmap annotations. Use when adding/colouring/annotating a plot (PCA, QC, t-SNE/UMAP, DE/expression heatmap) or adding a new colour/annotation attribute. Pair with shiny-module (module contract, dual_plot engine, mod_plot_subset, deferred render) and shiny-bslib-theming.
metadata:
  type: project
  version: "1.0"
---

# Plot colour / annotation aesthetics

Every plot that colours, groups, or annotates by a **per-sample attribute** (PCA, the
QC ggplot plots, the QC correlation heatmap, and the future t-SNE/UMAP and DE/expression
heatmap pages) goes through **one shared catalog + resolver** in `R/aes_helpers.R`. Do
**not** hand-roll colour scales or per-attribute `if/else` in a module — that duplication
is exactly what this machinery removed (it used to live three times: PCA, QC plots, the
heatmap annotation).

> Mechanics this skill does **not** cover (see the `shiny-module` skill): the
> `dual_plot()` ggplot↔plotly engine, the `deferred()`/`stale_note()` render gate, and
> the `mod_plot_subset` "Plot Showing" control. This skill is only about **what attribute
> can be picked and how it resolves to values + colours**.

## The attribute model

An attribute is identified by a `key`:

| key | source | kind |
|---|---|---|
| `"<colData col>"` | a sample-metadata column | discrete or continuous (by class; NA → an explicit `"NA"` level) |
| `"__gene__"` | gene expression (PCA only — that page owns the gene-search UI) | continuous |
| `"__qc__<metric>"` | a per-sample QC metric (`library_size`/`detected`/`pct_mito`) | continuous |
| `"__removal__"` | the QC suggested-removal status | discrete |
| `"__pool__"` | removal-pool membership (`Kept`/`In removal pool`) | discrete |
| `"__spike__<m>"` | ERCC dose-response (`slope`/`r_squared`/`lod`/`n_spike_detected`) | continuous, **only when the dataset has spike-ins** |
| `"__qc__pct_spike"` | `% spike-in` — grouped under **Spike-in** (also spike-gated) | continuous |

`aes_catalog(state, gene = FALSE)` returns the available descriptors —
`list(key, label, group, kind, loc)` — gated to the current dataset (gene only when
asked; spike only when `aes:::.aes_has_spike(state)`). `group` is one of **General /
This session / Spike-in / Data metadata** (fixed selector order).

## Resolver API (`R/aes_helpers.R`)

- **`aes_resolve(state, key, samples, ctx = list())`** → `list(values, kind, label, colors,
  labels, ramp_config)`, or **`NULL`** when not resolvable yet (e.g. removal before
  `state$samp_flags` exists, or a spike key on a spike-free dataset). `ctx$reason` (a QC
  metric) switches removal to the reason-aware 3-level highlight (General QC only);
  `"__gene__"` takes caller-precomputed `ctx$values`/`ctx$label`/`ctx$assay`.
- **`aes_ggplot_scale(res, aes_name = "colour")`** → a ready `scale_*_manual`/`scale_*_gradientn`,
  or `NULL` to keep thematic's default (no project palette configured).
- **`aes_heatmap_col(res)`** → a named colour vector (discrete) or a `circlize::colorRamp2`
  (continuous) for `ComplexHeatmap::HeatmapAnnotation(col = …)`.
- **`aes_annotation(state, keys, samples)`** → `list(df, col)` for a multi-track heatmap
  annotation (skips keys that resolve `NULL`).
- **Choices:** `aes_choices(catalog, kinds = c("discrete","continuous"), none, max_levels, state)`
  → grouped `selectInput` choices (filtered by kind; `max_levels` caps discrete columns for
  a *shape* selector). Built on `group_field_choices(coldata_cols, session_items, none,
  spike_items)` in `R/utils_group_choices.R` (the General → This session → Spike-in → Data
  metadata optgroup builder; flattens to a bare vector when only Data metadata remains).

**Every consumer reads `state` (palette configs, promoted `samp_flags`/`samp_pool`, the
`qc_metrics`/`spike_dr` derived caches), so all pages resolve an attribute identically.**
`aes_helpers` never reads module inputs — it stays page-independent.

## Consumer patterns

- **ggplot (PCA, QC General/RLE/density/spike):** `res <- aes_resolve(...)`; plot
  `res$values` with `aes_ggplot_scale(res)` (or `NULL` → thematic). A discrete-grouping
  ggplot factorizes continuous values (`if (is.factor(res$values)) res$values else factor(...)`)
  — see `group_map`/`sample_aes` in `mod_qc.R`.
- **ComplexHeatmap annotation (QC correlation heatmap):** **snapshot the annotation
  *values* in the deferred `spec` (so they update only on Render), but resolve the
  *colours* live in the `renderPlot` draw** (`aes_resolve` + `aes_heatmap_col`, keyed by
  `res$label`) so a Palette edit recolours without a re-render. The draw passes a
  precomputed `anno_col` to `.qc_correlation_heatmap()`. Preserve this split.
- **Display-only:** colour/shape/annotation are cheap re-plots — they must **not** enter
  the `state_derive` cache key or the `deferred()` `sig` (e.g. the PCA embedding is cached
  on `assay/n_top/log`, never on `colour_by`).

## The palette standard (load-bearing)

**Every attribute selectable in a colour/annotation control has a palette-config slot and
is editable on the Palette page.** The slot is the descriptor's `loc`:
`state$palette$colData[[col]]` / `state$palette$assays[[assay]]` /
`state$palette$other[[item]]`. The resolver reads it; `aes_ggplot_scale`/`aes_heatmap_col`
build the scale/ramp from it; **no config → thematic/viridis default** (never a
hard-coded Okabe-Ito).

- The customizable **"Other"** items are **`aes_other_palette_items()`** (removal_status,
  `__pool__`, the `__qc__*` + `__spike__*` ramps) **plus** the app-internal `correlation`
  ramp. `mod_palette.R`'s `.pal_other_meta()` = that function + `correlation`; it drives
  the "Other" pill's items, kinds, levels, classes, and friendly labels (`dom_item_label`).
- `removal_status_colors(config)` + `.removal_labels`/`.removal_labels_2` (in
  `filter_helpers.R`) own the reason-aware green/amber/red. (`removal_palette()` was
  removed — use the resolver.)
- **Feature-axis note (P7e):** the Palette page has a **`geneset`** discrete domain
  (`state$palette$geneset[[<annotated-set name>]]` → level→hex), driven by the session
  gene-set store, **not** part of the sample-axis `aes_helpers` catalog. It colours the
  Gene Sets Annotation tab's composition bar (`gene_set_anno_colors()`) and — in P7e-2 —
  the Expression heatmap's **row (feature) annotation** via a separate feature-axis
  resolver (`expr_row_annotation()`), the row-axis analog of `aes_annotation()`. Don't
  route feature/gene annotation through the sample-only `aes_helpers`.

## Adding a new colour/annotation attribute (recipe)

1. **Catalog:** add a descriptor in `aes_catalog()` — `key`, `label`, `group`, `kind`,
   and a `loc` palette slot. Gate it (`if (.aes_has_spike(state))` etc.) if it isn't
   universally available.
2. **Resolve:** add a branch in `aes_resolve()` returning per-sample `values` (+ `colors`
   for a fixed discrete map, or `ramp_config = state$palette$<loc>` for continuous).
   Return `NULL` if it can't resolve yet so consumers degrade (`validate`/skip).
3. **Palette slot:** list it in `aes_other_palette_items()` (so it's customizable) — unless
   it lives in `colData`/`assays` (already customizable). Add a `preset_for()` default in
   `mod_palette.R` if it needs one.
4. **Choices:** nothing — `aes_choices`/`group_field_choices` pick it up by `group`. Only
   extend `group_field_choices` for a *new optgroup*.
5. **Tests:** catalog membership/gating, `aes_resolve` per-sample alignment, the choices
   optgroup, and (if customizable) a config round-trips through `palette_to_json`.

The PCA / QC / heatmap consumers then offer + colour by it with **no edits** (they read
the catalog/resolver generically).

## Shared state fields (cross-page)

- `state$samp_pool` — staged removal pool (sample ids); `state$samp_flags` — latest
  `flag_samples()` frame (or `NULL`). Promoted from the QC page (`mod_qc` proxies the pool
  and mirrors the flags) so other pages can colour/annotate by removal/pool. Session UI
  state — no `data_version` impact; cleared on data load.

## Checklist

- [ ] Colour/group/annotation resolves via `aes_resolve` (+ `aes_ggplot_scale` /
      `aes_heatmap_col` / `aes_annotation`) — no bespoke per-attribute scales.
- [ ] Selector choices come from `aes_choices(aes_catalog(state, …))` (or
      `group_field_choices` directly for a fixed session-item set).
- [ ] New attribute → `aes_catalog` + a `loc` palette slot (+ `aes_other_palette_items`
      when "Other") + availability gate; `aes_resolve` returns `NULL` when not ready.
- [ ] Display aesthetics never enter the `state_derive` key / `deferred` `sig`.
- [ ] Heatmap annotation keeps the values-snapshot-in-spec / colours-live-in-draw split.
- [ ] No-config colouring falls back to thematic/viridis (not a hard-coded palette).
