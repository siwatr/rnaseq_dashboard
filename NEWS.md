# ddsdashboard (development version)

* Phase 7 (Expression) — Single genes + Gene sets heatmap next.

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
