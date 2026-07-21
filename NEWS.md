# ddsdashboard (development version)

# ddsdashboard 0.4.3

* **Gene sets > Heatmap** (P7c): a `ComplexHeatmap` over a named gene set (blank until you pick
  one). Render-only with a spinner — a heatmap is slow, and a single static plot has no cheap live
  layer, so every control is snapshotted on Render and a settings change shows a stale banner. The
  value matrix (any assay / VST / normalized log-counts) is cached behind the gate; per-gene
  z-score is the default (toggle raw `log10`), with constant rows drawn flat rather than breaking
  the plot. Row/column labels are duplicate- and NA-safe (a gene symbol shared by several ids shows
  for each; a missing label falls back to the id) and can show all / a searched selection / none,
  defaulting by size. Sample annotation, separate cluster + 3-state (Off/Auto/On) dendrogram
  controls, and a divergent ramp centred at 0 for z-scores (the legend keeps the underlying value,
  e.g. "z-score / log10(TPM + 1)"). Because a redraw is expensive, colour and annotation changes are
  gated too — they mark the plot stale rather than auto-redrawing. The Render button and stale banner
  sit above the plot; the sidebar adds Collapse/Expand-all and an in-app Plot-size panel (height px /
  width %), plus the "Showing:" sample subset.
* Extracted a reusable **continuous-palette control** (`continuous_palette_ui`/`_server`) from the
  Palette page so any plot can carry its own colour ramp (base palette + anchors + reverse + custom
  stops); the heatmap uses it now and the planned shared heatmap controller will too.
* **Heatmap k-means** (P7d): a "Clustering (k-means)" panel splits the heatmap into row and/or
  column clusters (k = 1 turns it off), computed on the displayed values with a visible seed and a
  Redo button for a fresh clustering; slices are labelled `C1 / (count)`. Each **row cluster can be
  saved as a gene set** in one click (prefix defaults to the source set's name, with kmeans
  provenance recorded), and the **column clustering can be saved to a colData column** for reuse as
  a sample annotation (samples hidden by the "Showing" subset are marked `unclustered`). So clusters
  flow into the rest of the app, the gene-set export, and the annotation controls. Leave the seed
  blank for a fresh (non-reproducible) clustering each render; toggle whether the k-means slices are
  ordered by clustering. Saving is conflict-guarded — a name clash prompts Overwrite / Abort rather
  than silently overwriting, the name field clears on success, and saving a column-cluster column no
  longer resets the annotation selectors.
# ddsdashboard 0.4.2

* Size-factor normalization is decoupled from assay assignment and gets its own **Size factors**
  Input sub-tab — a pill group of **Estimate / Per-sample / Compare**. *Estimate*: choose what to
  estimate on (endogenous / spike-in / a custom set via gene search / all genes) and the
  `estimateSizeFactors` type; every type honors the control set via a control-gene row-subset
  inherit (so `iterate`, which ignores `controlGenes`, respects it too). The config rides on the
  dds so structural edits re-estimate under it. Size factors are materialized by default at load
  (an object's own factors are respected), so DE / PCA / Expression read one consistent, visible
  normalization instead of re-deriving it; DE never writes them back. Re-estimating the same
  config is a no-op.
* New **Per-sample** and **Compare** size-factor plots. Per-sample shows the current factors as a
  bar (or points grouped by a `colData` variable). Compare is a consumer-only scatter of two
  size-factor vectors estimated under two different methods (x=y line + linear-fit R², both
  centered to geometric mean 1) to gauge how similar two normalizations are — it never mutates the
  dataset. Size factor is also a plottable **General QC** metric.
* **Design convention:** the variable of interest is now the *last* term of the design formula
  (covariates first, DESeq2 results-name style). A single resolver (`primary_design_var()`) drives
  the default group/colour/contrast pickers across DE / PCA / Expression / QC / Size factors so
  they agree. Term order does not affect DE correctness (results are extracted via `contrast=`).
* Adding a normalized assay (CPM/TPM/FPKM) no longer invalidates the DESeq2 fit or the cached
  VST. Those are now *content-addressed* — keyed on a fingerprint of the count matrix, sample /
  feature sets, and size factors (`dds_content_fingerprint()`) rather than the coarse
  `data_version` — so a structure-preserving edit preserves them, while a real sample/feature
  drop or a size-factor change still invalidates them.

# ddsdashboard 0.4.1

Phase 7 (Expression) begins — the **Single genes** tab (P7a) and the **Gene sets > Aggregate
expression** pill (P7b). The old "Heatmap" nav entry is renamed **Expression** (a
`navset_card_tab`: *Single genes* + *Gene sets*, the latter a pill group of *Aggregate
expression* + a *Heatmap* stub).

* **Single genes**: one feature at a time as a layered violin → boxplot → dots overlay, grouped
  by a colData variable and coloured by any attribute (the shared `aes_helpers` resolver). The
  expression value control gains opt-in **VST** and **normalized log-counts** choices; the value
  matrix is cached behind a deferred Render gate. Sample-count guards decide which geoms are
  offered; a dot-layout selector (`ggbeeswarm` beeswarm/quasirandom, `geom_jitter` fallback) plus
  per-layer width/opacity/size controls and a **y-axis limit** (out-of-range points draw as
  boundary triangles).
* **Gene sets > Aggregate expression**: the same overlay, but the y value is a per-sample
  **gene-set score** (`expr_set_aggregate()`) — drop absent + optionally all-zero-count genes,
  transform, **per-gene z-score by default**, then mean/median down to one vector. Source = a
  saved set or a quick uncommitted search; the subtitle reports how many of the set's genes were
  used. Default value = TPM when feature lengths exist, else VST.
* **Auto-render** controls are now toggles (`bslib::input_switch`) app-wide; the Expression pills
  gate the gene/set source behind the Render button when auto-render is off.

# ddsdashboard 0.4.0

Phase 6 (Gene Sets) complete — adds the **Compare tab** (P6e) alongside the Manage tab.

* Compare tab (a `navset_card_pill` **Stats** | **Overlap**) laid out like the DE Plots tab:
  the two shared controls — a "Sets to visualize" multiselect and a "Within this dataset only"
  toggle (`gene_set_present` vs the full authored membership) — live in the tab sidebar and
  govern both sub-pills.
* Stats: a present/absent set-size bar on the shared plot engine, with a horizontal/vertical
  toggle, an "Order by" control (none / size / name), and colours from a new static
  `other/gene_set_presence` Palette item.
* Overlap: Euler / Venn / UpSet, where the diagram-type default follows the set count (Euler
  ≤3, Venn 4, UpSet ≥5) and a type that cannot draw the current count shows a message rather
  than silently substituting another plot.
* `eulerr` draws the area-proportional Euler + Venn diagrams (a CRAN install — it has no
  conda-forge Apple-Silicon build; see `environment.yml`), degrading to `ggVennDiagram` (Venn
  only) when absent; UpSet uses `ComplexHeatmap`.

# ddsdashboard 0.3.1

Phase 6 (Gene Sets) — the Manage tab (P6a–P6d).

* Gene Sets page: define, record, and manage named gene sets of interest (DE seeds them;
  the Phase 7 Expression heatmap will consume them). Non-destructive storage — a set keeps
  its full authored membership; "present / absent in the dataset" is a live derived view,
  so a data edit only notifies about newly-absent ids, never trims the set.
* Shared free-text gene-search module (`mod_gene_search`) with exact / contains / regex
  match modes and tiered miss hints, retrofitted into PCA and DE (DE gained the explicit
  "Search by" column picker).
* Staging Manage tab: a "Build a gene set" card (Paste / From DE DEGs / Top-variable /
  Import table / Gene set file → a live Preview → New or Add-to-existing Save) beside a
  "Your gene sets" store.
* Tabular import (CSV / TSV / XLSX) with view/filter/select, ID-column + match-field pick
  (auto-detected), 1:many keep-all/first name matching, and annotation-split into N sets.
* File round-trip: JSON / GMT / long-TSV import (with the same ID-match scheme) and a
  selective, previewable export.

# ddsdashboard 0.3.0

Phase 5 — differential expression (DESeq2).

* Design & contrast builder shared between the Dataset "Design" sub-tab and the DE
  page, with a full-rank check and guided reference-level selection. Design changes are
  design-scoped: they invalidate the DE fit but not the PCA/VST/QC caches.
* `Run DESeq2` fits the model (contrast-free, cached under a data + design stamp),
  separate from per-contrast result extraction (reactive/auto or on demand), with
  apeglm -> ashr -> normal LFC shrinkage fallback.
* Results carry both standard and shrunk log2 fold-changes plus matching significance /
  DEG classifications; thresholds re-classify without re-fitting.
* DE Plots tab: MA / volcano / direct-comparison on the shared plot engine, per-feature
  DEG colouring, field-based axis limits (real coord limits + triangle clamping), gene
  labels, and the contrast + DEG stats embedded in the plot title/subtitle.
* Results Table tab with DEG colouring and a significant-only filter; a shared
  "Contrast & thresholds" control group across all DE tabs, plus per-contrast DEG
  summary tables.
* Shared expression-value control (assay / transform / pseudocount) for the
  direct-comparison plot, reused later by the Expression page.
* Curated "DEG palette" set (Pink-Blue / Orange-Purple / Red-Blue / Coral-Teal) alongside
  the generic discrete palettes, configured on the Palette page's Other -> DEG.

# ddsdashboard 0.2.0

Phase 4 — dimensionality reduction (PCA) and the shared colour/annotation system.

* Single-panel PCA (`mod_dimreduc`): VST-by-default input with assay/log advice,
  top-variable-gene selection (endogenous only), PC-axis selectors, colour by metadata
  or gene expression, shape by metadata, a scree %-variance bar, and deferred rendering.
* Shared plot engine extracted from the QC page (`mod_plot_engine`): the ggplot<->plotly
  `dual_plot`, the deferred render gate, and the staleness note.
* Shared colour/annotation attribute catalog + resolver (`aes_helpers`): the single place
  PCA, the QC plots, and heatmap annotations turn a per-sample attribute into values +
  colours, reading the Palette configs.
* Sample-removal state promoted to app level and exposed on PCA; QC colour selectors
  unified through the resolver.

# ddsdashboard 0.1.0

Initial release — Phases 1–3.

* App skeleton: `page_navbar` shell, one Shiny module per page, and the `new_app_state()`
  store with the `state_*()` API (working/original objects, `data_version`, history,
  undo/reset).
* Input: load a `dds`/`sce` or assemble one from a counts matrix + sample sheet; editable
  sample/feature metadata; OrgDb + GTF annotation; normalized assays (CPM/TPM/FPKM) with
  endogenous-only size factors.
* QC & filtering: per-sample QC metrics, dataset diagnostics (VST mean–SD), sample
  correlation, ERCC spike-in dose-response, and a removal-pool filtering workflow with
  auto-flagging.
* Edit-history controls (global undo/reset, scoped metadata reset) and the ggplot<->plotly
  engine toggle.
* Palette page: project-wide discrete + continuous colour configuration with JSON
  import/export.
