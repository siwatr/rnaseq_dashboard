# Shiny App for exploring and basic processing of RNA-seq output

Hi, I'm planning to make a R shiny app for exploring the results obtained from an RNA-seq analysis pipeline. The input for this app can be an R object (supporting DESeqDataSet (dds) and SingleCellExperiment (sce)).

> **Reviewed scope (2026-06):** the app is **bulk-first**. The first version targets bulk RNA-seq (and low-input single-cell such as SMART-seq, < 1k cells); full single-cell support is a later phase but the design must not preclude it. See "Reviewed design decisions" below.


## Goal 

The aim is to allow people with little to no bioinformatic knowledge to explore their datasets, as well as manipulating and performing basic analysis on these input object. 

***
## Reviewed design decisions (2026-06)

These were settled in a design review and refine the rough plan below. Where a section conflicts with this list, this list wins.

1. **Scope: bulk-first.** Build the bulk path end-to-end first. An `sce` is supported via two routes: **pseudobulk** aggregation by a sample-metadata grouping (`scuttle::aggregateAcrossCells()`), or **per-cell coercion to `dds`** behind a prominent *"statistically inaccurate and slow"* warning. **Force pseudobulk above ~1k cells** (configurable). Why per-cell DESeq2 is unsound: see [dev_ref/normalization-scran-vs-deseq.md](dev_ref/normalization-scran-vs-deseq.md).
2. **Input formats:** accept a **counts matrix + sample sheet** (CSV/TSV/XLSX) and build the `dds` in-app, in addition to a saved `.rds` `dds`/`sce`.
3. **Annotation (three-tier, organism-first):** (a) **OrgDb** by organism + ID type (auto-detected; build optional) with ERCC auto-detect, (b) **full GTF upload**, (c) **hybrid** — OrgDb base + manual fill for unmatched features, prompting a GTF when too many are missing. **A supplied GTF is authoritative** and overrides OrgDb for matching features.
4. **`rowData` conventions (always present):** a `feature_class` factor (`endogenous` default / `spike_in` / `exogenous`) on every feature, and `feature_length` when determined (feeds TPM/FPKM). Spike-in/exogenous features are **flagged, never dropped** — excluded from size-factor estimation (`controlGenes`) and variable-gene selection, but still plottable.
5. **Feature unit is adaptive:** features may be genes, transcripts, etc. Detect and label the UI dynamically; prompt when unknown; let the user correct the guess.
6. **Dimensionality reduction is PCA-focused** for bulk; t-SNE/UMAP are gated behind a low-sample warning (they suit the single-cell path).
7. **DE results carry both standard and shrunken LFC** (`log2FoldChange` + `log2FoldChange_shrunk`), each with its own `sig`/`DEG` (`sig_shrunk`/`DEG_shrunk`); the shrinkage toggle just selects which column set drives the plot. Significance uses `padj` (shared) + the chosen LFC.
8. **State, invalidation & reproducibility** follow a single model (see "State management & reproducibility" before the roadmap): an immutable `original`, a mutable `working` `dds`, version-stamped cached `derived` results, conservative downstream invalidation, reset/undo, and an action-log that exports a runnable R script / Quarto report.

***
<!--
#region MARK: App Usage
-->
## What this can be used for:

* QC of the datasets -- e.g., detecting samples with low library size, small detected features, high mitochondrial reads, or samples with poor correlation with counts of the same condition
* Modification of the object -- modification of sample and feature metadata (e.g., via `SummarizedExperiment::colData()` and `SummarizedExperiment::rowData()`, respectively). 
* Filtering out low quality samples
* Filtering out uninformative features (e.g., low counts genes)
* Transforming quantification (e.g., raw counts -> logcounts, or raw counts -> TPM, etc.)
* Visualizing PCA, tSNE, or UMAP build from hight variant genes.
  * Modify shapes or color by categorical or continuous value from sample metadata. 
  * Or change color based on expression of certain genes.
* Perform DESeq2 analysis.
* Visualization of differential analysis results:
  * make MA or volcano plots.
* Building gene expression heatmaps of gene set of interest.
* Exporting plots, DESeq2 results table, or R object (e.g., processed dds object)


***
<!--
#region MARK: Design
-->
## design:

The app can be multi-page dashboard design, where:

<!--
MARK: 1st page 
-->
### The first page: Input data

The first page is for input loading and metadata manipulation. User can look at the basic info of its datasets. If `sce` object is given, it should be converted to `dds` object as we should be able to perform differential expression analysis of this dataset at some point.


#### Manipulating sample metadata

The most minimum form of input is a single `dds` object with only raw count matrix. 
Ideally, the user can modify sample metadata as they want inside the app. Alternatively, they can upload a sample metadata table (.xlsx, .tsv, or .csv format) in which the row names or value in certain column match sample names of the `dds` object. 

> **Reviewed (2026-06):**
> - Sample metadata is edited on a dedicated **"Sample info"** tab (feature/`rowData` editing lives on a separate **"Feature info"** tab). The read-only colData preview on the Load tab is removed as redundant.
> - Editing uses a **draft model**: in-cell edits, add/remove/rename column, rename sample, and sheet-merge all accumulate in a local draft and are committed with an explicit **Save** (one history entry); **Reset to last save** and **Reset to original** discard. (Earlier auto-apply-on-edit gave too little feedback.)
> - **Protected columns** = design variables (`all.vars(design(dds))`): cannot be removed; renaming one rewrites the design formula to follow. Sample renames must stay unique/non-empty.
> - A merged sample sheet **overwrites** same-named columns (new wins) and the UI reports which were overwritten.
> - The viewing table uses DT's **built-in per-column filters** (`filter = "top"`: categorical dropdowns + numeric ranges; AND across columns). **Wishlist (deferred):** an *advanced filter builder* — rows of criteria with an AND/OR logic gate + column selector + value/range. Deferred because it is heavy and bug-prone, and a custom server-side filtered view conflicts with in-table cell editing (row-index mapping); built-in filters already cover the common AND-across-columns case. Revisit if OR-across-columns becomes necessary.


#### Manipulating feature (gene) metadata

We expected the input to have many genes (more than 50,000 genes in the case of mouse genome), thus manual manipulation of `rowData` can be challenging. User can upload the GTF annotation of their datasets (read via `rtracklayer` package). They can select what feature type should be used for their dataset, which column corresponding to their feature names (`rownames(dds)`), then they can select some columns from `mcols(gtf)` that will be incorporated in the feature metadata (`rowData(dds)` or `rowRanges(dds)`). 
The feature length (important for TPM or FPKM normalization) can optionally be determined here based on certain column of `rowData(dds)` (some quantification tools report these values and it might be coming with the input already), alternatively it can be determined from the GTF that is read in the current session. 

Certain group of genes (explicitly specified or narrow down using `dplyr::filter(as.data.frame(mcols(gtf)))`) can be marked as exogenous genes (e.g., being over-expressed) or spike-in genes (check ERCC-spike-in). Note that I have little experience with UMI, so I don't know if we should also handle it here as well?

__Note__ on feature length, the effective feature length may not comming from the same type of feature represented in the datasets! For example, if gene ID or gene names is used as the feature, but quantification of read count was done only on exon, then the effective feature length used for TPM and FPKM normalization should come from the exonic part of that gene.

> **Reviewed:**
> - **Annotation is three-tier and organism-first** (decision 3): default to OrgDb by organism + auto-detected ID type; GTF upload stays for custom/exogenous features and, when provided, **overrides** OrgDb for matching features. Bundle an ERCC reference in `inst/extdata/` for spike-in auto-detection and the dose–response plot.
> - **Effective feature length from a GTF** uses the direct `rtracklayer` route (robust to custom GTFs), not a `TxDb`: `import(gtf)` → keep `type == "exon"` → split by gene id → `reduce()` → `sum(width(...))`. Store as `rowData(dds)$feature_length`.
> - **Always add `rowData(dds)$feature_class`** (`endogenous` default / `spike_in` / `exogenous`). Mark, don't drop — these features stay in the object (so exogenous expression is still viewable in PCA/boxplots) but are excluded from size-factor estimation and variable-gene selection.
> - **Feature type is adaptive** (decision 5): detect gene/transcript/… and relabel the UI; prompt if unknown; allow correction. Drives the `<feature_type>_name` lookup column.

> **Implemented (2026-06):** Sample info and Feature info share one **draft metadata editor** (`mod_meta_editor`) over `colData`/`rowData`: editable cells, per-column filters, add / **multi-remove** / rename columns, Save / Reset. `feature_class` is edited **inline** (validated to its allowed values) or set in bulk on the filtered rows — replacing the old manual tagging. Sample info adds sheet-merge + sample rename; protected columns (design vars; `feature_class`) can't be removed.
>
> **Implemented (2026-06):** The Input page is organized by the three `SummarizedExperiment` components as tabs — **Load | Sample info | Feature info | Assay** (replacing the old standalone "Process" page). **OrgDb annotation** (organism + auto-detected ID type → `<feature_type>_name` and `description`, filling only where mapped) lives in the **Feature info** tab (chromosome is left to the GTF, since OrgDb's `CHR` accessor is deprecated); **normalization assays** (CPM always; TPM/FPKM when a complete `feature_length` exists) + endogenous size factors live in the **Assay** tab.
>
> **Implemented (2026-06):** **GTF annotation** in the **Feature info** tab (`gtf_helpers.R`): upload a GTF/GFF, match dds rows by an auto-resolved column (Ensembl `gene_id` else `gene_name`, version-stripped; overridable), and import selected `mcols` columns into `rowData` — **authoritative over OrgDb** (`gene_name`→`<feature_type>_name`, `seqnames`→`chromosome`), filling where matched and never wiping unmatched features. **`feature_length` is an explicit, optional action** with two sources: adopt an existing numeric `rowData` column, or compute the union length over a **user-chosen feature `type`** (default `exon`; `gene`/`transcript` for whole-body/nascent quantification) via `import()`→subset `type`→split→`reduce()`→`sum(width())`. Partial GTF coverage leaves `feature_length` incomplete (NA for unmatched), so `has_feature_length()` stays FALSE and TPM/FPKM remain off until complete. Still to come: **ERCC concentration**/dose–response (with the QC page, using real ERCC92 values).

<!--
MARK: 2nd page 
-->
### The second page: Quality control and filtering

The next step will be QC of the datasets. This can include: 
* A bar and/or box plots of QC matrix (e.g., by `scater` package). Bar plot could work when there are not many samples within the datasets. The box plot (x axis groupped by certain column of the `colData()`) can be used when input contain many samples. The QC for this category can include:
  * Library size (in millions).
  * Number of detected features.
  * Percent mitochondrial reads
  * (when applicable) percent spike-in reads
  * (when applicable) Scatter plot of known spike-in concentration and detected expression (TPM or FPKM)
  * Is there any common QC for dataset with UMI? (we don't expect to deal with this type of dataset, but it's worth thinking about it).
* A Variant stabilization plot
* Correlation heatmap between samples, preferably by log-counts values. 

Feel free to suggest any other common QC matrix/plots!

> **Reviewed:**
> - **One unified QC page** with a visible **data-type badge** (bulk vs single-cell). The plot types are shared; only the unit of observation changes (per-sample for bulk, per-cell for single-cell). `scater::perCellQCMetrics()` for single-cell; direct per-sample computation (or coercion) for bulk.
> - The **"Variant stabilization plot"** = `vsn::meanSdPlot()` on VST-transformed data (and/or `DESeq2::plotDispEsts()`).
> - **% mitochondrial / spike-in** depend on annotation — if no annotation is loaded, disable the metric and say why rather than erroring.
> - **UMI:** for bulk, UMI counts are just deduplicated counts; no special handling needed downstream.


#### Sample filtering

The app can suggest which samples should be discarded, and for what reason (e.g., those suggested by `scater`). Users can also decide to adopt our suggestion, and/or manually select additional samples to be excluded from the dataset. 

#### Feature filtering

Feature filtering can be done in semi-automating manner. features with zero counts in all samples should be removed. Features with low counts (either imposed by hard number of `rowSums(count_matrix) < x`; where `x` is a hard threshold set up user, or some other automatic way (please suggest)) can also be removed.

__Note__ that the function from the package `HTSFilter` can take long time to run, so we will avoid using it.

> **Reviewed:**
> - **Sample filtering:** for typical small-n bulk designs (3–12 samples), MAD-based auto-outlier detection (`scater::isOutlier`) is unreliable — *show* the metrics with gentle flags and let the user decide, rather than auto-suggesting aggressive drops.
> - **Feature filtering:** use `edgeR::filterByExpr()` (fast, design-aware) as the smart automatic default, alongside the manual `rowSums` threshold. Always drop all-zero features.


<!--
MARK: 3rd page 
-->
### The third page: Processing

The process in this step can be somewhat optional. User should be able to add new assays, such as TPM, FPKM, or CPM, to the dataset.

> **Reviewed:** assays are computed by the package's `cpm()/tpm()/fpkm()` helpers. **Size factors are estimated on endogenous genes only** (`DESeq2::estimateSizeFactors(dds, controlGenes = which(feature_class == "endogenous"))`). TPM/FPKM require `rowData(dds)$feature_length`; **degrade to CPM** when length is unavailable. All count-derived assays are recomputed when features are added/removed (see state model).


<!--
MARK: 4th page 
-->
### The fourth page: Dimensional Reduction and gene expression:

The fourth page can be used for exploring dimensional reduction plots (hereafter _ReduceDim_ plot for short), such as PCA, tSNE and UMAP, that constructed from one of the following methods:
* PCA extracted from DESeq2 function
* Manual plot of PCA, tSNE, and UMAP.

> **Reviewed:** the bulk path is **PCA-focused** (PCA, optionally MDS). t-SNE/UMAP need many points to be meaningful, so on small bulk designs they are **gated behind a low-sample warning** (they become first-class on the single-cell path). When showing multiple panels, **compute the embedding once** and let panels differ only by aesthetic mapping (e.g. colour by condition vs. replicate). Non-endogenous features are excluded from variable-gene selection.

User can specify which assay can be (default: `logcounts`), as well as how many top high variant genes (all genes or any specified number, default: 500) to be used for making ReduceDim plot.

#### Page panel: 

Ideally, this page should be able to show up to 4 plots panels at the same time. This can either be: 
(a) 4 ReduceDim plots, or 
(b) 3 ReduceDim plots with 1 box plot (prioritized being the last plot (right, or bottom right panel, if possible)) or violin plot of expression of specific gene of interest.

The page can start by 1 RedeceDim plot (e.g., PCA), and only adding more plot panel if user decided to add more panel -- e.g., if they want to look at plot colored by condition and look side-by-side with the same plot colored by biological replicate. Thus, the possible layout for this page is 1, 2, and 4 panels.


#### Plot customization: 

The color and/or shape of each data point (sample) can be modified based on sample metadata (discrete or continuous) or expression of certain genes (continuous values). 

Note that some feature ID, such as gene ID are not very human readable, and thus we would want to have a helper function to look up for the ID (i.e., which row of the `dds`) from values within certain column of the `rowData(dds)`. By default, this should be `<feature_type>_name` column, for example, if feature type was specified to be `gene`, then we are looking for a `gene_name`. Fall back to `rownames(dds)` if the default column doesn't exist.

<!--
MARK: 5th page 
-->
### The fifth page: Differential expression analysis

User can choose to perform differential expression analysis by `DESeq2::DESeq()` function at this step here. Since some datasets will already have DESeq2 results attached to it, so this step is optional. However, it will be required to re-run `DESeq2::DESeq()` function every time that a sample or feature is removed from the current instance of `dds` object. 

DESeq2 results will be extracted as data frame for making plots.

User can also design the formula for the analysis. Threshold for considering which feature are differentially expressed can be set by user:
  * log2-fold-change: default log2(2) -> `abs(log2FoldChange) >= log2(2)`
  * adjusted p-value: default: 0.05 -> `!is.na(padj) & (padj < 0.05)`

A two new columns will be added to the DESeq2 results:
* `sig` (logical): `TRUE` if the feature is significantly differentially expressed.
* `DEG` (factor - `c("up", "down", "no_change")`): Noting what type of differential expressed genes (DEG) each feature is.

If the DESeq2 result is available, then user can pick the results names (constructed by 3 values: column, test, control). Then making one of the 3 types of scatter plots:
* MA plot (starting point) -- 
  * x: `log10(baseMean)`
  * y: `log2FoldChange`
* Volcano plot
  * x: `log2FoldChange`
  * y: `-log10(padj)`
* direct expression comparison:
  * x: mean expression of the control condition
  * y: mean expression of the test condition

Each feature (point) is colored based on value from one of the columns in DESeq2 result table (Default: DEG column).

User can have some freedom to modify the range of x- and y-axis. Data point outside the range will be (virtually) set to the value of the limit, with triangle shape.

> **Reviewed:**
> - **Guided design, not free-text.** Build the formula via a UI: pick the variable of interest + optional covariates, set the **reference (control) level**, and validate the model is **full rank** before running. Replace free-text "results names" with a **contrast picker** listing factor levels. Support **multiple stored contrasts**.
> - **LFC shrinkage is precomputed, toggled at view time.** At DE time, store both `log2FoldChange` (standard `results()`) and `log2FoldChange_shrunk` (`DESeq2::lfcShrink()`, apeglm/ashr), each with its own `sig`/`DEG` (`sig_shrunk`/`DEG_shrunk`). `padj` is shared. The shrinkage toggle selects which column set drives the plot/colour — no recompute. Significance default: `!is.na(padj) & padj < 0.05 & abs(log2FoldChange) >= log2(2)`.
> - DE must **rerun from raw counts** on the current sample/feature subset whenever the data changes (state model).

<!--
MARK: 6th page 
-->
### The sixth page: Expression heatmap

User can import set of genes of interest, and this can be used to make a gene expression heatmap plot. By default, plots will be made using differentially expressed genes. 
Default value for expression is log10(TPM + 0.01). If TPM is not specified, call back to log10(CPM + 0.01) instead. 
Heatmap top annotation can indicate value in sample metadata (default: values in column used as design for `dds`)

Additional subpanel can be added as a space to change the color mapping of heatmap color (palette, color range, etc), column/row annotation color mapping. This can include other heatmap setup: e.g., 
* boolean value setting for `ComplexHeatmap::Heatmap()` function -- allowed only show column/row name and dendrogram at the first version.
* Other noon-boolean setting for `ComplexHeatmap::Heatmap()`, e.g., row kmean.

> **Reviewed:**
> - **Default scaling = per-gene z-score (row scaling)**, with a toggle to show raw `log10(TPM/CPM + 0.01)`. Raw values let high-expressors dominate.
> - **Row (gene) names hidden by default** (heatmaps usually carry many genes); let the user mark genes of interest, labelled via `ComplexHeatmap::anno_mark()` (leader lines) so "where is my gene" works without thousands of labels.
> - **Column (sample) names** are user-toggleable and **auto-hidden when > 30 samples** (default, configurable).
> - Report any imported genes-of-interest that weren't found. Default gene set = DEGs, falling back to top-variable genes when no DE has been run. Export via explicit device capture (`png()/pdf()` + `draw()`), not `ggsave()`.

<!--
MARK: 6th page 
-->
### The seventh page: Data export

This is for exporting the dds object, DESeq2 result tables, and plots.

> **Reviewed:** also export a **reproducibility artifact** — a runnable R script and/or Quarto/R Markdown report — generated from the action log (see state model). It captures the loaded data, filters, design, thresholds, and package versions, doubling as the provenance trail and the publishable record. Result tables export to XLSX.

***
## State management & reproducibility

A single model serves invalidation, reset/undo, and reproducible export. App state (a `reactiveValues`) holds:

* **`original`** — the object exactly as loaded; never mutated; the reset target.
* **`working`** — the current `dds` after edits/filters/normalization.
* **`derived`** — cached heavy artifacts (size factors, VST, PCA, DESeq fit, DE result tables), keyed by name.
* **`data_version`** — bumped on **any** edit to samples, features, or design-relevant metadata.
* **`history`** — an ordered **action log** (one entry per user operation, with parameters).

**Invalidation (v1, conservative):** each derived artifact records the `data_version` it was computed under and is **stale** when that differs from the current `data_version`. Any edit bumps the version, marking everything downstream stale; the UI shows a **stale badge** + re-run button and never silently reuses stale results. Dependency order: `raw counts → filter → normalized assays → variable genes → PCA → DESeq fit → DE results → heatmap`.

**Reset / undo:** *reset-to-original* sets `working <- original` and clears `derived`; *undo* keeps a short stack of previous `working` snapshots (start depth 1).

**Reproducibility:** the `history` log renders to the exportable R script / Quarto report — same record that drives invalidation. A **persistent dataset-status bar** on every page surfaces this state (data type, n samples/features, assays, design set?, DE current vs stale).


<!-- 
#region MARK: Roadmap
 -->

## Road map

The project is implemented in phases, **bulk-first** (revised in review). Cross-cutting concerns — the state/invalidation model, caching, progress indicators, reproducibility export, and **mock-`dds` test fixtures** — are built early and threaded throughout.

### Phase 1 (skeleton — done)
* Basic Shiny app layout (bslib `page_navbar`) + persistent dataset-status bar
* Data import (1st page): `.rds` **and** counts-matrix + sample-sheet tabular input; show `show(dds)`
* Data export page shell

### Phase 2
* Metadata manipulation (editable tables); annotation (OrgDb-first, then GTF); `feature_class` flags
* QC (unified bulk/single-cell page) + sample/feature filtering (`filterByExpr`)
* Normalization assays (CPM/TPM/FPKM; size factors on endogenous)

### Phase 3
* PCA-focused dimensional-reduction page

### Phase 4
* DESeq2 steps: guided design + contrast builder + reference level; LFC shrinkage (dual columns)
* DE plots (MA / volcano / direct comparison); expression heatmap (row z-score default)

### Phase 5 (later — single-cell)
* Single-cell ingestion, per-cell QC (`scran`/`scater`), pseudobulk aggregation, warned per-cell DESeq2; t-SNE/UMAP for single-cell

<!-- 
#region MARK: Roadmap
 -->

## misc

* Consider using `ggplot2` package as the main plot package. `ComplexHeatmap` is preferred as main package for heatmap.
* Since showing the plot can take a while, we can add a button that will trigger re-making of the plot when user is satisfied with their setting.
* logcounts assay can be added automatically to the dds object. Consider doing this with CPM as well, but this has to be updated each time the new features is added or removed.
* Consider making git commit each time you add a considerable, meaningful change to the repo.

