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

## Phase 4 — Dimensionality reduction ⬜
- PCA-focused (t-SNE/UMAP gated by sample count). Top-variable genes (default 500, assay
  default `logcounts`, endogenous only). Up to 4 panels (1/2/4); compute embedding once,
  vary aesthetics per panel; colour/shape by metadata or gene expression (`lookup_feature()`).
- Embed the reusable "Plot Showing" control (per the convention in CLAUDE.md).

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
