# Roadmap — ddsdashboard

Living, detailed phase plan + status. **This file is the source of truth for "what's
done / next".** `CLAUDE.md` links here; the design narrative (the *why*) lives in
[rough_design.md](rough_design.md). Keep this updated as PRs merge.

Legend: ✅ done · ⬅️ next · ⬜ later

Build is **bulk-first**. Cross-cutting concerns (state/invalidation model, caching,
progress indicators, reproducibility export, mock-`dds` fixtures) are threaded throughout.

---

## Phase 1 — skeleton ✅
- bslib `page_navbar` layout + persistent dataset-status bar.
- Import: `.rds` `dds`/`sce` **and** counts-matrix + sample-sheet (CSV/TSV/XLSX).
- Export page shell.

## Phase 2 — annotation & normalization ✅
- Editable metadata (shared draft editor `mod_meta_editor` over colData/rowData).
- Annotation: OrgDb-first, then authoritative GTF; `feature_class` flags.
- Normalization assays (CPM/TPM/FPKM); size factors on endogenous `controlGenes`.
- GTF memory minimization (`mod_gtf_reader`, optional repeatable filtering + Remove) +
  session-memory (RSS) status-bar badge.
- Feature-info UI redesign (one `navset_card_pill` per page; per-tab sidebar + preview +
  match-coverage banner). Theming pass (custom palette + dark mode + `thematic`).

## Phase 3 — QC, filtering, ERCC (sub-PR'd)
- **P3a ✅** per-sample QC metrics + Sample QC tab (`qc_helpers.R`, scater + base-R fallback). [PR #11]
- **P3b ✅** dataset diagnostics (VST mean–SD), RLE, expression density, Sample Correlation
  heatmap + within-group correlation; mock Palette page. [PR #12]
- **P3c ✅** filtering (Samples/Features pills: auto-flag → removal-pool → `state_mutate`;
  view-only "Showing:" subset; removal-status colouring). [PR #13]
- **(edit-history) ✅** global Undo/Reset in status bar; scoped metadata "Reset to original"
  (auto-commit) + unsaved badge; QC "Reset Sample/Feature Removal"; net-edit + undo-limit
  badges. [PR #14] — see [object-states-and-resets.md](object-states-and-resets.md).
- **P3d ✅** ERCC spike-in dose-response (titration QC) + spike-in visibility (Feature class
  counts, status badge, "Spike-in Conc." designation, remove-all-spikes). Bundled ERCC
  Mix1/2 ref (`inst/extdata/ercc_concentrations.csv`). Plus polish: reusable "Plot Showing"
  module (`mod_plot_subset.R`), `dt_table` read-only by default, selection-based
  feature_class assignment, per-value DT colours, auto-render threshold → 150. [PR #15]
- **P3e ✅** **Filtering by spike-in QC metrics** — QC → Filtering → Samples flags/selects samples by
  spike content & fit (two-sided % spike-in fence, detected-spike count, dose-response R²/slope,
  +`<3-point` reason), feeding the existing removal-pool flow. `flag_samples()`/`suggest_sample_thresholds()`
  gain opt-in spike criteria fed by the shared `spike_dr` cache; reuses the Spike-in QC tab's
  source/assay (observed assay now **TPM > FPKM > CPM**-preferred, the former deferred fix).
  Filtering sidebar reorganized into collapsible "Sample QC filters" / "Spike-in (ERCC) filters"
  accordions with scoped Auto buttons. New knowledge note [ercc-spike-in.md](ercc-spike-in.md). [PR #17]
- **P3f ✅** **ggplot ↔ plotly engine toggle** — global status-bar "Interactive plots" switch
  (default off = static ggplot) writing `state$plot_interactive`; QC ggplot plots render via a
  dual-output container (`dual_plot` → static `plotOutput` / interactive `plotlyOutput` via
  `ggplotly`). Gated **per plot on an element budget** (estimated rendered glyphs ≈ rows of the
  plotted data, the true `ggplotly` cost driver) rather than a sample cap; budget is an `option()`
  (`ddsdashboard.plotly_max_elements`, default 5000 — tunable without a settings page). Over budget,
  a plot falls back to static with a **sticky per-plot "Render interactive anyway"** override (reset
  when data changes or the toggle flips off/on). Sample-name (+ value) hover `text` aes on the 7
  toggled plots. **Excluded:** VST mean–SD and the ComplexHeatmap correlation heatmap. `plotly` is a
  `Suggests` (on-demand, graceful fallback). [PR #18]
- **P3g-a ✅** project-wide **Palette page** — discrete **colData** wiring. New `palette_helpers.R`
  engine (`palette_qualitative*`, `norm_color`, `palette_discrete`) is the single resolver for discrete
  level→colour mappings: explicit pins layered on a named base palette (ggplot default / Okabe-Ito /
  Viridis-d / Set2 / Dark2). `qc_annotation_colors(df, config)` and the QC ggplot group scales (via a
  shared `palette =` plumbing + `.qc_group_scale()`) both read it, so the QC plots and the correlation
  heatmap agree; per-level pins are honoured. Palette page is empty-by-default, opt-in per column
  (accordion panel, `shinyWidgets::colorPickr` with `col2rgb` R-name/CSS normalization + textInput
  fallback, live preview). `state$palette` is a UI pref (untouched by load/reset, no `data_version`).
  Single optgroup palette selector (Custom / Qualitative / Brewer: Qual/Seq/Div / viridis), the full
  `viridisLite`/`RColorBrewer` catalogue, a Preview tab, Collapse-all. Config shape `list(name, colors)`.
  [PR #19]
- **P3g-b ✅** project-wide **Palette page** — **continuous palettes + the Feature/Assay/Other pills**.
  The page is a `navset_card_tab` (Setting / Preview); Setting holds Sample / Feature / Assay / Other
  pills, each a `layout_sidebar`. The panel machinery is generalized over four domains
  (`colData`/`rowData`/`assays`/`other`), each item discrete *or* continuous (by column type; assays
  continuous; `other` = a discrete removal-status map + a continuous correlation ramp). New continuous
  engine: `palette_continuous_choices()`, `palette_resolve_range()` (number / `p<pct>` / data-range
  anchors), `palette_colorramp2()` (→ `circlize::colorRamp2` for heatmaps) and `palette_gradientn()`
  (→ ggplot `scale_*_gradientn` pieces). Wired consumers: numeric `colData` annotations + the
  correlation-heatmap ramp (`qc_annotation_colors` / `.qc_correlation_heatmap` `cor_config`) and the
  QC removal-status map (`removal_palette()`); rowData/assay configs are stored for the P4/P5
  consumers. Continuous config shape `list(name, min, max, reverse, custom)`. Custom continuous ramp
  is an **N-stops (2-5) selector** whose pickers resample via `colorRampPalette` when the count changes
  (default white -> black); reverse-direction checkbox; high-cardinality guard (warn/cap options);
  Type/Class accordion badges; presets (`removal_status`, `correlation`). [PR #TBD]
- **P3g-c ✅** Palette **config import/export** (JSON via `jsonlite`) + **Config tab**. New
  `palette_to_json()`/`palette_from_json()` pure helpers round-trip `state$palette` as a versioned,
  faithful mirror (kind inferred from keys; discrete `colors` serialized as a `{level: hex}` object,
  `custom` as an array). The Config tab (3rd, beside Setting/Preview) is a `layout_sidebar`: **selective
  export** (one `checkboxGroupInput` per non-empty domain, all-checked default, Select/Deselect all; an
  `export_palette()` reactive feeds both a live JSON preview and the download, so the preview is exactly
  what downloads — `ddsdashboard-palette-<date>.json`); **import** via one Replace/Merge/Cancel modal that
  classifies each dataset-present item (kind mismatch + over-cap → skip+report; a known palette whose
  colours diverge → a conflict radio Keep-colours/Force-palette) and re-wires items via the new
  **register-on-visible** observe. Plus an **"Edit palette" button** on non-custom continuous palettes:
  extracts the ramp to 5 editable Custom-ramp anchors (current `reverse` baked into the order, then
  zeroed — gradient unchanged). Robustness: `dom_levels()` is `other`-item-aware. [PR #TBD]
- **P3g-d ⬅️ next** Palette **factor management** — coerce a `colData`/`rowData` column to factor + reorder its
  levels (drives both plot order and the palette mapping). (The legacy Themer gallery can be retired —
  its Heatmap sub-tab is already superseded by the Preview tab.)

## Phase 4 — Dimensionality reduction (PCA) (sub-PR'd)
- **P4-pre ✅** Extracted the shared plot engine (`R/mod_plot_engine.R`): the ggplot↔plotly
  toggle (`dual_plot`/`use_plotly_base`), the deferred render gate (`deferred`/`stale_note`),
  and the plot helpers (`.plot_msg`/`.plot_dual`/`.plotly_max_elements`/`.muffle_unknown_aes`/
  `.to_plotly`), as a host-namespace submodule (`plot_engine_server()`) reused by QC + PCA.
  Behaviour-preserving refactor of `mod_qc.R`. [PR #22]
- **P4a ✅** Single-panel **PCA**. New `R/dimreduc_helpers.R`: `pca_assay_advice()`
  (tiers inputs — recommended VST/logcounts/normalized-log-counts, log-first CPM/TPM/FPKM,
  unsuitable raw counts), `pca_input()` (endogenous-only matrix + honest subtitle label; VST
  default with broad logcounts fallback; `norm_logcounts` estimates endogenous size factors),
  `top_variable_features()`, `compute_pca()` (`prcomp(t(top), center, !scale)`, %var from
  `sdev^2`, deterministic PC sign). The page (`R/mod_dimreduc.R`) caches the embedding via
  `state_derive` behind the deferred Render gate; PCA scatter with PC-axis selectors, colour
  by metadata (discrete/continuous Palette configs) or gene expression (logcounts gradient),
  shape by discrete metadata (NA→"NA"), the input-used subtitle, a **scree** %-variance bar,
  and the reusable "Showing" subset (display-only — embedding/axes stay stable). Default
  colour-by = first design var → condition → none; no-config colouring uses the thematic
  default (not Okabe-Ito); dark-mode aware (shared `.plot_theme`); plot-aesthetics accordion
  (1:1 ratio + point size). [PR #23]
- **P4a-2 ✅ done** PCA **colour & feature-search enhancements** (sidebar accordions; grouped
  colour/shape selectize with a "This session" group = gene expression + per-sample QC metrics;
  always-embedded gene-colour block with field-aware Search-by, dup-column toggle, case-insensitive
  lookup, debounced "Did you mean?" hint, expression assay + transform + pseudocount; shape ≤6-level
  guard; legend-position select). Pure helpers `resolve_feature()`/`suggest_features()`/
  `feature_search_choices()`/`expr_transform()` + `lookup_feature(case_insensitive=)`. QC metrics
  share the `derived` cache; gene-colour guards a stale assay name. [PR #24]
- **P4b** multi-panel (1/2/4 layouts), embedding computed once, aesthetics per panel.
- **P4c (later)** t-SNE/UMAP, hard-gated by sample count (~≥30); Rtsne/uwot optional Suggests.
- **P4-removal ✅ done** Promoted the removal-pool + sample-flag state from `mod_qc.R` into shared
  `state` (`samp_pool` + `samp_flags`; `mod_qc` proxies the pool to state and mirrors the computed
  flags). PCA exposes **"Suggested removal"** / **"In removal pool"** under the colour-by + shape-by
  "This session" group, resolved from shared state. Extracted `group_field_choices()` (optgroups
  **General → This session → Data metadata**) and applied it across the **QC per-sample colour
  selectors** (General, RLE, Expression density, Spike-in) with the same session items, routed
  through session-aware `group_map()`/`group_colours()`; shared removal colour scheme +
  `removal_status_colors()` moved to `filter_helpers.R`. **Decision (revised after implementation):**
  the **grouping-semantic** selectors (Within-group correlation "Group by", filtering "Within-group
  grouping") + the Feature selector stay **colData-only** — grouping a correlation/flagging
  computation by removal-status/pool isn't a meaningful biological grouping, and "Suggested removal"
  on the flagging selector is *circular* (the flags would group by their own output).
  The **sample-correlation heatmap "Annotate by"** multi-select also adopts the grouped layout and
  offers the session items (Suggested removal / Removal pool, discrete) **and** the per-sample QC
  metrics (library size / detected / % mito / % spike-in, continuous) as top-annotation tracks,
  via a mixed-type `anno_df` + merged session-track colours through `qc_annotation_colors()`. [PR #25]
  Still-open follow-ups it sets up: **promote the "Showing" subset to app-state** so other plot
  pages reuse it; and P4b/P5 plot pages can now reuse the session removal aesthetics.

### P4a-2 — detailed plan (PCA colour & feature-search)

Builds on P4a's `mod_dimreduc.R`; pure helpers in `dimreduc_helpers.R`/`utils_lookup.R` (tested).
Sidebar reorganised into **bslib accordions** to manage clutter: **Embedding** (input assay · log ·
n_top · Render/auto, open) · **Appearance** (colour/shape + an always-visible Gene-expression block) ·
**Plot aesthetics** · **Showing**.

- **Grouped colour/shape selectize** (`selectInput` optgroups, like the Palette selector): **Data
  metadata** = `colData` columns; **This session** = Gene expression + per-sample QC metrics
  (`qc_per_sample_metrics()` → library_size/detected/pct_mito/pct_spike, continuous). (Removal-status
  / in-pool join this group in P4-removal.) Shape-by: discrete `colData` columns **filtered to ≤6
  levels** (hide the rest), label "Shape by (discrete, ≤6 values)"; keep the runtime guard defensively.
- **Gene-expression block — always embedded** (not conditionally revealed): when Colour by ≠ Gene
  expression, show helper text "Select 'Gene expression' in 'Colour by' to enable this." Controls:
  - **Search by** `selectInput`: candidate `rowData` fields **excluding logical columns** (keep
    integer/numeric — some IDs are numeric, e.g. Entrez) + **"Feature ID (rownames)"** (always; unique
    by construction). Default = `<feature_type>_name` if present, else Feature ID.
  - **"Include columns with duplicate values"** `input_switch` (tooltip): OFF hides columns where
    `anyDuplicated(col) != 0`; default **ON when the logical-filtered `rowData` has ≤10 columns**, else OFF.
  - **Case-insensitive** `input_switch` for the lookup.
  - **Duplicate handling** (case-folding or genuine dups): take the **first** match + warn
    "N features matched '<q>'; showing the first."
  - **Search hint** (on an exact miss, query length **≥2**): search the chosen field **always
    case-insensitively (independent of the case-insensitive toggle**, so case-sensitive "duxf3" still
    suggests "Duxf3"); **raw-match cap = 100** → over cap show "too many partial matches; type more to
    narrow it down" (no list); else rank by match position (prefix `^q` > word-start > substring;
    tiebreak shorter length then alpha), show **top 5** + "(+N more)". Pure tested helper
    `suggest_features(query, values, n=5, cap=100)`; **debounce the gene input ~300 ms**; pre-lowercase
    the search vector cached per `data_version`. Message tiers: found→colour; miss+hits→"Did you
    mean …?"; miss+too-broad→refine; miss+none→"not found".
  - **Expression assay** `selectInput` (any assay; default logcounts) + **Transformation** select
    (none [default] / log2 / log10); when log, a **Pseudocount** numeric appears — default **1 for
    integer assays** (detect via `all(m == round(m), na.rm=TRUE)`), else **0.5**. Already-log assays
    (logcounts/VST) default to none. Colourbar label reflects assay + transform.
- **`lookup_feature()`/`suggest_features()`**: add case-insensitive matching + a match **count**
  (for the warning); the "Feature ID" choice searches rownames; `suggest_features` does the ranked
  regex hint. All unit-tested (prefix ranks first, cap triggers, min-length gate, duplicate count).
- **Plot aesthetics**: swap the 1:1-ratio checkbox for **`bslib::input_switch()`**; add a
  **legend-position** `selectInput` (right [default]/left/top/bottom/none → `theme(legend.position=)`,
  overriding `.plot_theme`'s bottom for PCA).

## Phase 5 — Differential expression + heatmap ⬜
- Guided design builder (reference level + full-rank check) + contrast picker; multiple
  stored contrasts; `DESeq2::DESeq()` + `lfcShrink`. Dual LFC columns
  (`log2FoldChange`(_shrunk)) + `sig`/`DEG`(_shrunk); shrinkage toggle = column selection.
- MA / volcano / direct-comparison plots with axis clamping (triangle markers).
- Expression heatmap (`ComplexHeatmap`; default per-gene z-score; `anno_mark` for genes of
  interest; column names auto-hidden > 30 samples).
- **ComplexHeatmap plot-controls PR** — the heatmap controllers (annotation/colour/clustering
  options) for the DE expression heatmap and the QC sample-correlation heatmap get a dedicated PR
  here, so we can factor which control sub-modules are **shared** vs **QC-only** (the P3f engine
  toggle deliberately left both heatmaps static — they are not ggplot).
- **Create the `de-analysis` skill** (guided design/contrast/shrinkage, dual-LFC schema,
  MA/volcano conventions) as part of this phase, so it matches the real DESeq2 implementation.

## Phase 6 — Export & reproducibility ⬜
Fill in the Export page (shell since P1; currently downloads the processed `dds` only):
- **Reproducibility R script / Quarto report** generated from the `history` action log (loaded data,
  edits, filters, design, thresholds, `sessioninfo`) — the provenance trail + publishable record.
- **Plot export** via explicit device capture (`png()`/`pdf()` + `draw()`/`print()`), not `ggsave()`.
- **DE result tables → XLSX** (`writexl`); processed `dds` export (already works).

## Phase 7 — single-cell (later, lowest priority) ⬜
- SCE ingestion, per-cell QC (`scran`/`scater`), pseudobulk aggregation
  (`scuttle::aggregateAcrossCells`), warned per-cell DESeq2 (< ~1k cells), t-SNE/UMAP.

---

## Deferred / wishlist (revisit when relevant)
- **Full bright/dark theme customization (post-P5)** — an (almost) fully customized bslib light/dark
  theme (brand colours, typography, component styling) coordinated with `thematic` and a `plotly`
  layout theme so interactive plots match the app. Plan after the DESeq2 phase (P5). This is also
  when P3f's plotly figures get themed (currently `ggplotly` does not inherit `thematic` — accepted).
- **Two-sided spike highlight on the General QC % spike-in plot** — the "Suggested removal" colour-by
  maps `pct_spike` to the over-spiked side only; an under-spiked sample shows as
  "suggested (other)". Symmetrise if users ask.
- **Within-group-correlation auto-flag tuning** (threshold/z-score UX beyond v1).
- **Advanced filter builder** for metadata tables (rows of criteria + AND/OR gate) — deferred
  in P2; built-in per-column DT filters cover AND-across-columns; revisit if OR is needed.
- **Global "Showing:" subset** promoted to an app-state field once P4 adds more sample plots
  (currently QC-page-local via `mod_plot_subset`).
- **Undo/Reset depth** (`.undo_depth = 5`) — revisit if large/single-cell data makes the
  snapshot stack heavy.

## Merged-PR log (recent)
| PR | Title |
|----|-------|
| #11 | P3a: per-sample QC metrics + Sample QC tab |
| #12 | P3b: QC dataset diagnostics, sample correlation & within-group QC |
| #13 | P3c: sample/feature filtering, auto-flagging & display subset |
| #14 | Edit-history controls: global Undo/Reset, scoped metadata reset, QC reset-removal |
| #15 | P3d: ERCC spike-in dose-response & spike-in QC (+ polish bundle) |
| #17 | P3e: filtering by spike-in QC metrics (+ ERCC reference note) |
| #18 | P3f: ggplot↔plotly engine toggle (global switch + sample cap + hover labels) |
| #19 | P3g-a: project palette — discrete colData wiring |
| #20 | P3g-b: continuous palettes + Feature/Assay/Other pills |
| #21 | P3g-c: palette config JSON import/export + Edit palette |
| #22 | P4-pre: extract shared plot engine (mod_plot_engine) |
| #23 | P4a: dimensionality reduction — single-panel PCA |
| #24 | P4a-2: PCA colour & feature-search enhancements |
| #25 | P4-removal: promote removal state; PCA removal/pool; QC colour-selector unification; heatmap session/QC annotations |
| #27 | anno-refactor-A: shared color/annotation attribute catalog + resolver (parity) |

**v0.1.0** released after #21 (end of P3).

### Color/annotation standardization (epic, post-P4-removal)
One shared attribute **catalog + resolver** (`R/aes_helpers.R`) so QC + PCA + future pages (t-SNE/UMAP, DE heatmap) agree on the selectable colour/annotation attributes and every attribute's palette is customizable.
- **PR A ✅ (#27)** extract the shared catalog/resolver; migrate PCA + QC + heatmap annotation onto it (parity).
- **PR B ✅** generalize the Palette **"Other"** domain so removal-pool + per-sample QC metrics are user-customizable (`aes_other_palette_items()` is the single source; `.aes_pool` reads its config; pool seeds a grey/red preset, QC metrics viridis). **Project standard:** every colour/annotation attribute has a palette slot and is editable on the Palette page.
- **PR C ⬅️ next** spike-in QC metrics (`slope`/`r²`/`lod`/`n_spike_detected`) as attributes + a dedicated **"Spike-in"** optgroup, availability-gated on spikes, reusing the `spike_dr` cache.
