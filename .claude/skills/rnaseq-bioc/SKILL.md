---
name: rnaseq-bioc
description: Domain reference for working with RNA-seq Bioconductor objects in this app — DESeqDataSet/SummarizedExperiment/SingleCellExperiment structure, assays, normalization math (CPM/TPM/FPKM) and size factors on endogenous genes, feature_class/feature_length conventions, scater QC metrics, filterByExpr, and the DESeq2 dual-LFC (standard + shrunken) result schema. Use when writing or reviewing any code that touches the dds object, its assays, colData/rowData, normalization, QC, or DESeq2 results.
metadata:
  type: reference
  version: "1.1"
---

# RNA-seq Bioconductor reference

Encodes the domain knowledge this app depends on so it doesn't get re-derived (or gotten subtly wrong) each time. **Bulk-first:** the bulk path is built first; single-cell is a later phase via pseudobulk (see end). For bulk, UMI counts are just deduplicated counts — no special handling.

## Object model

A `DESeqDataSet` *is a* `RangedSummarizedExperiment`. Accessors you will use constantly:

| Need | Accessor |
|---|---|
| Count / value matrices | `assay(dds, "counts")`, `assays(dds)`, `assayNames(dds)` |
| Sample (column) metadata | `colData(dds)` — a `DataFrame`, one row per sample |
| Feature (row) metadata | `rowData(dds)` / `rowRanges(dds)` (a `GRanges` when ranges exist) |
| Dimensions | `nrow` = features, `ncol` = samples; `rownames` = feature IDs, `colnames` = sample IDs |
| Design / model | `design(dds)` |
| Size factors | `sizeFactors(dds)` |

`sce` → `dds`: build raw counts with `DESeq2::DESeqDataSet(sce, design = ~ 1)` (or `DESeqDataSetFromMatrix`), carrying over `colData` and `rowData`. A placeholder design (`~ 1`) is fine until the user sets one on the DE page. **For real single-cell input, prefer pseudobulk** (`scuttle::aggregateAcrossCells()` by a sample grouping) → a bulk-like `dds`; **force it above ~1k cells**. Per-cell coercion is allowed below that only behind a *"statistically inaccurate and slow"* warning — see [normalization-scran-vs-deseq.md](../../../dev_ref/normalization-scran-vs-deseq.md).

**Edits invalidate downstream state.** Dropping a sample or feature subsets the object (`dds[keep_rows, keep_cols]`) and means any stored `DESeq()` fit, normalized assays, and DE results are stale and must be recomputed. The app tracks this with a `data_version` stamp (see the `shiny-module` skill's state model) — never silently reuse a stale `derived` artifact.

## Assays and normalization

Keep raw integer counts in `"counts"` untouched; add normalized values as new assays.

- **logcounts** (auto-added on load): `log2(cpm + 1)` is a reasonable default, or `assay(dds, "logcounts") <- log2(counts + 1)`. Recompute on feature add/remove.
- **CPM** (counts per million): `cpm = counts / colSums(counts) * 1e6`. Library size is per-sample `colSums(counts)`; recompute when the feature set changes.
- **TPM** (needs effective feature length `L`, in bases): `rate = counts / L`; `tpm = rate / colSums(rate) * 1e6`. Length-normalize *first*, then per-sample scale.
- **FPKM**: `fpkm = counts / (L/1000) / (colSums(counts)/1e6)`.

**Effective length caveat:** `L` may come from a different feature type than the rows (e.g. exonic length for gene-level rows). Store it as **`rowData(dds)$feature_length`**, aligned to `rownames(dds)`. Source it from an upstream `rowData` column, or compute from a GTF via the **direct `rtracklayer` route** (robust to custom GTFs): `import(gtf)` → keep `type == "exon"` → split by gene id → `reduce()` per group → `sum(width(...))`. Avoid building a `TxDb`. Never assume `L` = the row's genomic span. TPM/FPKM are unavailable until `L` exists — degrade gracefully to CPM.

## Feature classes (always present)

Add **`rowData(dds)$feature_class`** unconditionally — a factor with levels `endogenous` (default for every feature) / `spike_in` / `exogenous`. Set during annotation (ERCC pattern `^ERCC-`, GTF/manual marking). **Flag, never drop** — these features stay in the object so their expression remains plottable, but their *role* changes:

- **Size factors on endogenous only:** `DESeq2::estimateSizeFactors(dds, controlGenes = which(rowData(dds)$feature_class == "endogenous"))`. (DESeq2's `controlGenes` selects which genes *form* the size-factor estimate.)
- **Variable-gene selection** (PCA, heatmap default): subset to `feature_class == "endogenous"` before ranking by variance (`matrixStats::rowVars` / `MatrixGenerics::rowVars`).
- **DE:** keep all genes in the fit; the column lets the UI optionally grey/hide spike-in/exogenous in MA/volcano.

## QC metrics (scater)

`scater::addPerCellQCMetrics()` / `perCellQCMetrics()` compute per-sample:

- **library size** (`sum`, total counts) — report in millions.
- **detected features** (`detected`, genes with count > 0).
- **% mitochondrial** — pass mito genes via `subsets=` (identify by chromosome `MT`/`chrM` in `rowRanges`, or a `gene_name` `^mt-`/`^MT-` pattern).
- **% spike-in** — for ERCC, the spike-in subset (`^ERCC-` names, or rows tagged on the input page).
- **spike-in dose–response** — **titration QC, not normalization** (size factors stay endogenous-`controlGenes`; spike-ins are QC-only here). Per-sample log–log `lm` of observed expression vs. known concentration: slope ≈ 1 + high R² = a healthy titration. Observed assay = a linear depth-normalized one — **prefer TPM/FPKM when available, else CPM (with a warning)**; **never counts or logcounts**. Drop zeros before the log; return NA below 3 usable points; report the lowest *detected* spiked concentration (LOD). Concentrations come from a `spike_concentration` rowData column or the bundled ERCC Mix 1/2 reference (`inst/extdata/ercc_concentrations.csv`). Helpers: `ercc_concentrations()`, `resolve_spike_concentration()`, `spike_dose_response()` in `R/ercc_helpers.R`.

Useful QC plots: per-metric bar plot (few samples) or box plot grouped by a `colData` column (many samples); a sample–sample correlation heatmap on log-counts (`ComplexHeatmap`, correlation of `assay(dds, "logcounts")`); and the **mean–SD variance-stabilization plot** — VST via `DESeq2::vst()` (fallback `varianceStabilizingTransformation()` on small/low-count data), plotting `rank(rowMeans)` vs per-gene `sd`, points coloured by local density via base `MASS::kde2d` (no `vsn`/`hexbin` dependency) with a red running-median trend; a flat trend ≈ variance stabilized.

Filtering: always drop all-zero features (`rowSums(counts) == 0`); offer **`edgeR::filterByExpr()`** (fast, design-aware) as the smart automatic default, plus a manual `rowSums(counts)` threshold. **Do not use `HTSFilter`** — too slow for an interactive app. For small-n bulk (3–12 samples), MAD-based outlier detection is unreliable — *flag* samples, don't auto-drop.

**QC page is unified** across bulk/single-cell (same plot types, per-sample vs per-cell unit) with a visible **data-type badge**. (The mean–SD variance-stabilization plot is the base-`MASS::kde2d` reimplementation described above — `vsn`/`hexbin` were intentionally dropped.)

## DESeq2 differential expression

Use a **guided design builder**, not free-text: pick the variable of interest + optional covariates, set the **reference (control) level** (`relevel()`), and validate the model is **full rank** before fitting. Replace free-text "results names" with a **contrast picker** over factor levels; support **multiple stored contrasts**.

```r
dds <- DESeq2::DESeq(dds)              # rerun from raw counts after any sample/feature change
res <- DESeq2::results(dds, contrast = c(column, test, control))
shr <- DESeq2::lfcShrink(dds, coef = ..., type = "apeglm")   # precompute once
df  <- as.data.frame(res)
df$log2FoldChange_shrunk <- shr$log2FoldChange[match(rownames(df), rownames(shr))]
```

**Precompute both LFC variants and classify each** (the shrinkage toggle then just selects columns — no recompute). `padj` is shared (shrinkage changes only LFC/lfcSE):

```r
classify <- function(lfc, padj, padj_cut, lfc_cut) {
  sig <- !is.na(padj) & padj < padj_cut & abs(lfc) >= lfc_cut
  list(
    sig = sig,
    DEG = factor(ifelse(!sig, "no_change", ifelse(lfc > 0, "up", "down")),
                 levels = c("up", "down", "no_change"))
  )
}
std <- classify(df$log2FoldChange,        df$padj, 0.05, log2(2))
shk <- classify(df$log2FoldChange_shrunk, df$padj, 0.05, log2(2))
df$sig        <- std$sig;  df$DEG        <- std$DEG
df$sig_shrunk <- shk$sig;  df$DEG_shrunk <- shk$DEG
```

Defaults: `padj_cut = 0.05`, `lfc_cut = log2(2)` (user-adjustable). Memory is a non-issue (~2–3 MB/contrast at 50k genes).

DE scatter plots (ggplot2; color by `DEG` by default):

- **MA:** x = `log10(baseMean)`, y = `log2FoldChange`.
- **Volcano:** x = `log2FoldChange`, y = `-log10(padj)`.
- **Direct comparison:** x = mean expression in control, y = mean expression in test.

Axis clamping: points outside a user-set range are drawn at the limit with a **triangle** shape (clamp the value for plotting, flag with a shape aesthetic) rather than dropped.

## Gene ID → name lookup

Feature IDs (Ensembl etc.) aren't readable. Resolve queries against `rowData(dds)`, defaulting to `<feature_type>_name` (e.g. `gene_name` when feature type is `gene`), falling back to `rownames(dds)` when that column is absent. Return the matching row index/ID.
