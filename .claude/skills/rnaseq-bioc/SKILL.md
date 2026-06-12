---
name: rnaseq-bioc
description: Domain reference for working with bulk RNA-seq Bioconductor objects in this app — DESeqDataSet/SummarizedExperiment/SingleCellExperiment structure, assay slots, normalization math (CPM/TPM/FPKM), scater QC metrics, and the DESeq2 differential-expression result conventions used here. Use when writing or reviewing any code that touches the dds object, its assays, colData/rowData, normalization, QC, or DESeq2 results.
metadata:
  type: reference
  version: "1.0"
---

# RNA-seq Bioconductor reference

Encodes the domain knowledge this app depends on so it doesn't get re-derived (or gotten subtly wrong) each time. Bulk RNA-seq only — UMI/droplet single-cell QC is out of scope.

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

`sce` → `dds`: build raw counts with `DESeq2::DESeqDataSet(sce, design = ~ 1)` (or `DESeqDataSetFromMatrix`), carrying over `colData` and `rowData`. A placeholder design (`~ 1`) is fine until the user sets one on the DE page.

**Edits invalidate downstream state.** Dropping a sample or feature subsets the object (`dds[keep_rows, keep_cols]`) and means any stored `DESeq()` fit, normalized assays, and DE results are stale and must be recomputed.

## Assays and normalization

Keep raw integer counts in `"counts"` untouched; add normalized values as new assays.

- **logcounts** (auto-added on load): `log2(cpm + 1)` is a reasonable default, or `assay(dds, "logcounts") <- log2(counts + 1)`. Recompute on feature add/remove.
- **CPM** (counts per million): `cpm = counts / colSums(counts) * 1e6`. Library size is per-sample `colSums(counts)`; recompute when the feature set changes.
- **TPM** (needs effective feature length `L`, in bases): `rate = counts / L`; `tpm = rate / colSums(rate) * 1e6`. Length-normalize *first*, then per-sample scale.
- **FPKM**: `fpkm = counts / (L/1000) / (colSums(counts)/1e6)`.

**Effective length caveat:** `L` may come from a different feature type than the rows (e.g. exonic length for gene-level rows). Treat it as a separately-sourced vector aligned to `rownames(dds)`, sourced either from a `rowData` column the upstream tool provided or computed from a GTF in-session. Never assume `L` = the row feature's genomic span. TPM/FPKM are unavailable until `L` exists — degrade gracefully to CPM.

## QC metrics (scater)

`scater::addPerCellQCMetrics()` / `perCellQCMetrics()` compute per-sample:

- **library size** (`sum`, total counts) — report in millions.
- **detected features** (`detected`, genes with count > 0).
- **% mitochondrial** — pass mito genes via `subsets=` (identify by chromosome `MT`/`chrM` in `rowRanges`, or a `gene_name` `^mt-`/`^MT-` pattern).
- **% spike-in** — for ERCC, the spike-in subset (`^ERCC-` names, or rows tagged on the input page).
- **spike-in dose–response** — scatter of known ERCC concentration vs. observed TPM/FPKM; a good linear fit is a quality signal.

Useful QC plots: per-metric bar plot (few samples) or box plot grouped by a `colData` column (many samples); a variance-stabilization mean–variance plot (`DESeq2::vst()`); and a sample–sample correlation heatmap on log-counts (`ComplexHeatmap`, correlation of `assay(dds, "logcounts")`).

Filtering: always drop all-zero features (`rowSums(counts) == 0`); drop low-count features by a user threshold on `rowSums(counts)`. **Do not use `HTSFilter`** — too slow for an interactive app.

## DESeq2 differential expression

```r
dds <- DESeq2::DESeq(dds)              # rerun after any sample/feature change
res <- DESeq2::results(dds, contrast = c(column, test, control))
df  <- as.data.frame(res)
```

Augment every results frame with the app's two standard columns:

```r
df$sig <- !is.na(df$padj) & df$padj < padj_cut & abs(df$log2FoldChange) >= lfc_cut
df$DEG <- factor(
  ifelse(!df$sig, "no_change", ifelse(df$log2FoldChange > 0, "up", "down")),
  levels = c("up", "down", "no_change")
)
```

Defaults: `padj_cut = 0.05`, `lfc_cut = log2(2)`. A contrast/result name is `c(column, test, control)`.

DE scatter plots (ggplot2; color by `DEG` by default):

- **MA:** x = `log10(baseMean)`, y = `log2FoldChange`.
- **Volcano:** x = `log2FoldChange`, y = `-log10(padj)`.
- **Direct comparison:** x = mean expression in control, y = mean expression in test.

Axis clamping: points outside a user-set range are drawn at the limit with a **triangle** shape (clamp the value for plotting, flag with a shape aesthetic) rather than dropped.

## Gene ID → name lookup

Feature IDs (Ensembl etc.) aren't readable. Resolve queries against `rowData(dds)`, defaulting to `<feature_type>_name` (e.g. `gene_name` when feature type is `gene`), falling back to `rownames(dds)` when that column is absent. Return the matching row index/ID.
