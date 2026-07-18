# Roadmap — ddsdashboard

Living, detailed phase plan + status. **This file is the source of truth for "what's
done / next".** `CLAUDE.md` links here; the design narrative (the *why*) lives in
[rough_design.md](rough_design.md). Keep this updated as PRs merge. **When a sub-phase
merges, collapse its pre-build detail to the status bullet + a Merged-PR log row** — the
shipped code is the truth, so verbose plan blocks shouldn't accumulate here.

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
- **P3g-d → moved to Phase 8** Palette **factor management** — coerce a `colData`/`rowData` column to factor +
  reorder its levels (drives both plot order and the palette mapping). (The legacy Themer gallery can be
  retired — its Heatmap sub-tab is already superseded by the Preview tab.) **Folded into Phase 8
  (Visualization enhancements)** so the path to DE is shortest; it is no longer the immediate next item.

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
- **P4b → moved to Phase 8** multi-panel (1/2/4 layouts), **ReduceDim panels only**, embedding computed
  once, aesthetics per panel. Real value: same embedding **side-by-side with different colour-by** (e.g.
  condition vs. biological replicate). Deferred for the complex multi-panel UI, not for lack of merit.
  (The original "mixed gene box/violin panel" idea is **dropped** — superseded by Phase 7 Expression >
  Single genes.)
- **P4c → moved to Phase 8** t-SNE/UMAP, hard-gated by sample count (~≥30); Rtsne/uwot optional Suggests.
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

## Phase 5 — Differential expression ✅
DE statistics + DEG plots (the heatmap split out to Phase 7). Sub-PR'd P5a→P5c.
- **P5a ✅** pure engine `de_helpers.R` — `classify`/dual-LFC schema, `de_full_rank`/`de_relevel`/
  contrast levels, `de_run`/`de_results`/`de_shrink` (apeglm→ashr→normal fallback), MA/volcano/direct
  builders + `de_clamp`/`de_colour_resolve`; hardened `estimate_size_factors_endogenous` (poscounts
  fallback); `apeglm`/`ashr`/`ggrepel` Suggests; the **`shiny-de-analysis` skill**. [PR #32]
- **P5b/P5b-2 ✅** shared `mod_design_builder` (in the Dataset **Design** sub-tab + the DE tab);
  **design-scoped** `state_set_design` (bumps `design_version`, not `data_version`); the DE Design &
  Contrasts tab — contrast picker with 3 validity tiers, **Run = fit only** (stored in `derived` under
  a stamp) + **reactive/auto per-contrast extraction** (cached in `state$de`); status-bar DE badge.
  [PR #33, #34]
- **P5c ✅** DE Plots (segmented MA/volcano/direct on the `dual_plot` engine; deferred sig = data-only;
  local per-feature colour via `de_colour_resolve` + the curated Palette **`other/DEG`** set;
  field-based axis limits → triangles; 1:1 Direct toggle; ggrepel labels) + Results Table (`dt_table`
  DEG-coloured + significant-only) + a **shared "Contrast to view"** selector. [PR #35]

## Phase 6 — Gene Sets ✅ (sub-PR'd P6a→P6e)
A dedicated page (between DE and Expression, a `navset_card_tab`: **Manage** | **Compare**) to
*define, record, manage, compare, and share* **named gene sets of interest** — DE *seeds* them,
the Phase 7 Expression heatmap *consumes* them. Full design + rationale in the approved plan:
`~/.claude/plans/vast-giggling-ripple.md` (mirrored by the `gene-sets-phase-plan` project memory).

**Load-bearing decisions:** a shared **`mod_gene_search`** sub-module (extracted from the inline
PCA/DE gene-search, retrofitted into both — DE gains PCA's explicit "Search by" column picker);
**non-destructive storage** — `state$gene_sets[[name]] = list(ids, kind, annotation, source)` keeps
the **full authored membership**, "present in dataset" is a live derived view (powers Compare's
"Within this dataset" toggle); sets are **snapshots** (source controls are live previews, Add
freezes ids); a per-add **New/Append** toggle; the **annotated** layer (id→label sets, combine/
overlap UI, a `Palette > Gene Set` sub-tab) is **deferred to P7** with its heatmap consumer
(storage is forward-compatible). Recorded in the reproducibility export (Phase 9).

- **P6a ✅** shared `mod_gene_search.R` (single/multi, host-namespace; Exact/Contains/Regex match
  modes + tiered miss hints) + retrofit PCA + DE; this roadmap expansion + the plan-pointer
  memory. [PR #36]
- **P6b ✅** `state$gene_sets` (structured, non-destructive) + `gene_set_helpers.R`
  (`new_gene_set`/`gene_set_commit`/`gene_set_present`/`gene_set_absent`) + `mod_geneset.R` and
  its **staging** Manage tab — "Build a gene set" (source pills → live Preview → New/Add Save,
  New rejects a name clash) beside "Your gene sets"; `staged()` = a named list following the
  active source, so the P6c import reuses the layout. + nav wiring. [PR #37]
- **P6c ✅** rich tabular import (CSV/TSV/XLSX `dt_table` view/filter/select, ID-column pick +
  a **1:many keep-all** resolver — *not* `lookup_feature`, see the ID-match convention —
  add-filtered / add-selected, annotation-split → N staged sets, live stats) **+ the multi-set
  Save UI** (staged > 1: pick sets, name prefix, conflict resolver / auto-rename, or union-add
  into existing). `.read_user_table` promoted to `load_helpers.R`. Doc-sync changed to per-PR. [PR #38]
- **P6d ✅** file round-trip — `gene_sets_to_json`/`_gmt`/`_tsv` + inverses + a
  `gene_sets_from_file()` extension/sniff dispatcher (JSON faithful — `source`/`kind`, 1-gene
  array via `I()`; GMT `source`→description; long TSV `set`/`id`/`annotation`) + a **"Gene set
  file" source pill** (named sets stage → the existing multi-set Save, non-destructive; reuses
  the table-import ID-match scheme — match field auto-detect, keep-all/first, keep-unmatched — so
  a foreign id scheme like `gene_name` resolves) + a standalone **Export** section (selective set
  multi-selector + Select/Deselect all + format selectize + a live capped Preview).
- **P6e ✅** Compare tab (a `navset_card_pill` Stats | Overlap, under a shared **"Sets to
  visualize"** multiselect — extracted `.gs_set_multiselect_ui/_server`, reused by Export — and a
  **"Within this dataset only"** toggle governing both, `gene_set_present` vs full `ids`). **Stats**
  = a present/absent set-size bar on the shared `dual_plot` engine (horizontal default + a Vertical
  toggle; colours from a new static **`other/gene_set_presence`** Palette item via
  `gene_set_presence_colors()`). **Overlap** = Euler / Venn / UpSet with a **per-type-cap
  auto-switch to UpSet** (`.gs_overlap_type`: Euler ≤3, Venn ≤4). Pure helpers `gene_set_size_frame`
  / `gene_set_overlap_list` / `gene_set_ids_for`. **Dep note:** `eulerr` (area-proportional Euler +
  Venn) has **no conda-forge osx-arm64 build** and `eulerr >= 7` needs Rust, so it is a CRAN
  install (documented exception in `environment.yml`); the conda-clean **`ggVennDiagram`** is the
  Venn fallback when eulerr is absent (Euler option hidden), UpSet via `ComplexHeatmap`. + doc-sync
  (per-PR) + propose `v0.4.0`.

## Phase 7 — Expression ⬅️ next ⬜
Renamed from "Heatmap": a gene-expression **browsing surface** (more than a heatmap). A
`navset_card_tab` with two tabs.
- **Single genes** — one feature at a time; a layered overlay (back→front: violin → boxplot →
  dots) reusing the shared plot machinery (`dual_plot` engine, deferred render, `mod_plot_subset`
  "Showing:", `aes_helpers` colour). Controls: feature search (reuse the PCA gene block); x-axis =
  a metadata grouping var (default first design var; `group_field_choices()`, colData-only); y-axis
  assay (default log-norm-counts when size factors exist, else TPM→FPKM→CPM→logcounts→counts via a
  new `expr_default_assay()`); group/colour-by (independent, `aes_resolve()`); transform+pseudocount
  (`expr_transform()`); plot-element toggles (`geom_violin`/`geom_boxplot`/`ggbeeswarm::geom_quasirandom`
  with `geom_jitter` fallback). **Sample guards** (`N` = max samples/group; `G1`/`G2` as `option()`s,
  proposed 10/50): dots default ON when `N ≤ G1`, never allowed when `N ≥ G2`; violin/box hidden when
  all groups `N < G1`. (Subsumes the single-gene box/violin panel the old design attached to the PCA
  multi-panel; *deferred to wishlist:* facet-by-2nd-variable, mean±error overlay.)
- **Gene sets** — a `ComplexHeatmap` over a named set from Phase 6: default per-gene z-score (toggle
  raw `log10(TPM/CPM + 0.01)`); row names hidden + `anno_mark()` for genes of interest; column names
  auto-hidden > 30 samples; top annotation via `aes_annotation()` (default = design column). Default
  set when none chosen = DEGs (fallback top-variable). **Basic controls in v1**; the shared
  heatmap-controller refactor + advanced options (row k-means) are Phase 8.

## Phase 8 — Visualization enhancements ⬜
The gathered deferred plot/UX items (no home in the build order until now):
- **P3g-d** Palette factor management — coerce a `colData`/`rowData` column to factor + reorder its
  levels (drives plot order + the palette mapping). Retire the legacy Themer gallery.
- **P4b** multi-panel PCA — ReduceDim panels only (1/2/4 layouts; embedding computed once, aesthetics
  per panel). Same embedding side-by-side with different colour-by (condition vs. biological replicate);
  deferred for the complex multi-panel UI.
- **P4c** t-SNE/UMAP — sample-count-gated (~≥30); Rtsne/uwot optional Suggests.
- **ComplexHeatmap shared plot-controls PR** — factor the heatmap controllers (annotation/colour/
  clustering) **shared** vs **QC-only** across the QC sample-correlation heatmap and the Phase 7
  Expression heatmap; advanced options (row k-means). (The P3f engine toggle left both heatmaps
  static — they are not ggplot.)
- Relevant **wishlist** items — full light/dark theme + plotly theming (post-DE), promote the
  "Showing:" subset to app-state, two-sided spike highlight, within-group-corr auto-flag tuning.

## Phase 9 — Export & reproducibility ⬜
Fill in the Export page (shell since P1; currently downloads the processed `dds` only):
- **Reproducibility R script / Quarto report** generated from the `history` action log (loaded data,
  edits, filters, design, thresholds, gene sets, `sessioninfo`) — the provenance trail + publishable record.
- **Plot export** via explicit device capture (`png()`/`pdf()` + `draw()`/`print()`), not `ggsave()`.
- **DE result tables → XLSX** (`writexl`); processed `dds` export (already works).

## Phase 10 — single-cell (later, lowest priority) ⬜
- SCE ingestion, per-cell QC (`scran`/`scater`), pseudobulk aggregation
  (`scuttle::aggregateAcrossCells`), warned per-cell DESeq2 (< ~1k cells), t-SNE/UMAP.

---

## Deferred / wishlist (revisit when relevant)
> Several of these are now **scheduled in Phase 8 (Visualization enhancements)** — see that phase
> for the theme/plotly-theming, "Showing:"-subset promotion, two-sided spike highlight, and
> within-group-corr tuning. They remain listed here for context until built.
- **Full bright/dark theme customization (post-P5 → Phase 8)** — an (almost) fully customized bslib
  light/dark theme (brand colours, typography, component styling) coordinated with `thematic` and a
  `plotly` layout theme so interactive plots match the app. This is also when P3f's plotly figures get
  themed (currently `ggplotly` does not inherit `thematic` — accepted).
- **Two-sided spike highlight on the General QC % spike-in plot (→ Phase 8)** — the "Suggested removal"
  colour-by maps `pct_spike` to the over-spiked side only; an under-spiked sample shows as
  "suggested (other)". Symmetrise if users ask.
- **Within-group-correlation auto-flag tuning (→ Phase 8)** (threshold/z-score UX beyond v1).
- **Advanced filter builder** for metadata tables (rows of criteria + AND/OR gate) — deferred
  in P2; built-in per-column DT filters cover AND-across-columns; revisit if OR is needed.
- **Global "Showing:" subset (→ Phase 8)** promoted to an app-state field now that more sample plots
  exist (PCA; soon the Expression single-gene plot) — currently QC-page-local via `mod_plot_subset`.
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
| #27 | anno-refactor-A: shared colour/annotation catalog + resolver (`aes_helpers.R`, parity) |
| #28 | anno-refactor-B: generalize Palette "Other" (removal-pool + QC metrics customizable) |
| #29 | anno-refactor-C: spike-in QC metrics as attributes + "Spike-in" optgroup |
| #30 | docs: shared colour/annotation resolver convention (`shiny-plot-aesthetics` skill) |
| #31 | docs: split DE/Expression phases, add Gene Sets page, gather deferred items |
| #32 | P5a: DESeq2 DE engine (`de_helpers.R`) + size-factor hardening + `shiny-de-analysis` skill |
| #33 | P5b: DE design/contrast builder + Dataset/Design tab + DESeq2 fit (design-scoped) |
| #34 | P5b-2: DE contrast validity + reactive extraction + UI polish |

**v0.1.0** released after #21 (end of P3).

> **Colour/annotation standardization epic ✅ (#27–#29):** one shared attribute catalog +
> resolver (`R/aes_helpers.R`) so QC + PCA + future pages agree on selectable attributes and
> every attribute's palette is customizable (Palette "Other"). See the `shiny-plot-aesthetics`
> skill for the full reference.
