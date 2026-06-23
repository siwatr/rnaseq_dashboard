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
- **P3e ⬅️ next** **Filtering by spike-in QC metrics** — let QC → Filtering → Samples flag/select
  samples by spike content & fit (% spike-in, detected-spike count, dose-response R²/slope,
  lowest-detected conc), feeding the existing removal-pool flow (`flag_samples()` gains spike
  criteria).
- **P3f ⬜** UI-polish bundle — primarily a **ggplot ↔ plotly engine toggle** for the
  interactive plots; sweep up small UX items.
- **P3g ⬜** project-wide **Palette page** full wiring (mock landed in P3b): feed `thematic`
  qualitative/sequential palettes + `qc_annotation_colors()` for ComplexHeatmap; pin
  metadata levels to fixed colours.

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
- **ERCC dose-response default assay** — should prefer **TPM/FPKM when available**, falling back to
  **CPM with a warning** (the current code defaults to CPM unconditionally). Fix in a future PR; the
  `rnaseq-bioc` skill already documents the desired behaviour.
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
