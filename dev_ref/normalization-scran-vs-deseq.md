# Normalization: scran vs. DESeq2 — and why it shapes the bulk-first decision

A developer reference for the `rnaseq_dashboard` project. The question: when we normalize counts, why does the app use **DESeq2 median-of-ratios size factors** for bulk, and **scran pooling/deconvolution** for single-cell — and why can't we just run DESeq2 size factors on a single-cell object?

## TL;DR

- **DESeq2 (median-of-ratios)** estimates **one size factor per sample** from the geometric-mean reference across genes. It assumes most genes are not differentially expressed and that each library has **enough non-zero counts** to form a stable ratio. Ideal for **bulk** (tens of samples, deep per-sample coverage).
- **scran (pooling + deconvolution)** estimates **one size factor per cell** by pooling many cells, normalizing the summed profile, then deconvolving back to per-cell factors. It was designed precisely because per-cell counts are **shallow and full of zeros**, which breaks the median-of-ratios geometric mean.
- Running DESeq2 size factors per-cell on a real single-cell matrix **fails or degrades**: the geometric-mean reference needs genes detected in *all* columns, and with sparse single-cell data that set is tiny or empty → many `NA`/zero factors. This is a core reason the app is **bulk-first** and **forces pseudobulk above ~1k cells**.

## DESeq2 median-of-ratios (the bulk default)

For sample *j*, DESeq2 computes a size factor `s_j` as the **median across genes** of the ratio of that sample's count to a **per-gene geometric-mean reference**:

```
ref_i = geometric_mean_j( K_ij )            # reference "pseudo-sample" for gene i
s_j   = median_i ( K_ij / ref_i )           # over genes with ref_i > 0
```

Counts are then divided by `s_j`. Key properties and assumptions:

- **Robust to composition / a few very-high genes** because it uses the *median* ratio, not the total — better than naive total-count (library-size) scaling.
- The **geometric-mean reference is only defined for genes with a non-zero count in *every* sample** (a single zero makes the geometric mean zero, and that gene is dropped from the median). In bulk this still leaves thousands of usable genes.
- Assumes **most genes are not DE** (the median ratio reflects technical depth, not biology).
- One factor per **sample**; the unit is a biological library with deep coverage.

DESeq2 offers `type = "poscounts"` (geometric mean over positive counts only) and `controlGenes =` to estimate factors on a chosen subset — the latter is exactly how this app keeps spike-in/exogenous features out of normalization (`controlGenes = which(feature_class == "endogenous")`).

## scran pooling + deconvolution (the single-cell default)

`scran::computeSumFactors()` / `pooledSizeFactors()` (often via `scran::quickCluster()` first) was built for the single-cell regime where each column is **one cell** with **low total counts and many zeros**:

1. **Pool** many cells together and sum their counts → a much denser pseudo-profile.
2. Normalize the **pooled** profile against a reference (the median-of-ratios idea now works, because pooling filled in the zeros).
3. Use **many overlapping pools** and solve a linear system to **deconvolve** the pool-level factors back to **per-cell** size factors.
4. Pre-clustering (`quickCluster`) keeps pools roughly homogeneous so the "most genes not DE" assumption holds within a pool.

The result is per-cell factors that are stable even when any individual cell is too sparse to normalize on its own. (After factors are set, `scuttle::logNormCounts()` produces the `logcounts` assay used for single-cell PCA/UMAP.)

## Why per-cell DESeq2 size factors break

- **The shared-detection set collapses.** Median-of-ratios needs genes detected across the columns it compares. With thousands of sparse cells, the set of genes non-zero in *all* cells approaches empty → unstable or `NA` size factors, or DESeq2 erroring out.
- **`poscounts` helps but doesn't fix the statistics.** Even with positive-count geometric means, a cell with a handful of detected genes gives a noisy single-cell factor; there's no pooling to borrow strength.
- **The DE model is wrong too.** DESeq2's negative-binomial dispersion model assumes replicate libraries; treating cells as independent replicates ignores within-sample correlation and inflates significance — hence the app's **"statistically inaccurate"** warning on per-cell DESeq2.
- **Cost.** Dispersion estimation and `nbinomWaldTest` over `genes × 10k+ cells` is slow and memory-heavy.

## How this informs the app

- **Bulk path:** DESeq2 median-of-ratios size factors, estimated on endogenous genes via `controlGenes`. This is the canonical, correct route and is what the current `cpm()/tpm()/fpkm()` utilities and the planned `logcounts` assay sit alongside.
- **Single-cell path (later phase):** scran pooling/deconvolution for per-cell factors + `scuttle::logNormCounts()`; per-cell QC via `scater`. Reuse the same QC *page* (per-cell unit) per the unified-QC plan.
- **The bridge — pseudobulk:** `scuttle::aggregateAcrossCells()` sums counts within a sample grouping, turning a single-cell object into a **bulk-like `dds`** where median-of-ratios and the DESeq2 NB model are valid again. The app **forces pseudobulk above ~1k cells** and recommends it below that, because pseudobulk restores the statistical footing that per-cell DESeq2 lacks.

## References (for implementers)

- Anders & Huber 2010, *Genome Biology* — median-of-ratios normalization (DESeq).
- Love, Huber & Anders 2014, *Genome Biology* — DESeq2.
- Lun, Bach & Marioni 2016, *Genome Biology* — pooling/deconvolution normalization for single-cell (scran).
- Amezquita et al. 2020, *Nature Methods* — "Orchestrating single-cell analysis with Bioconductor" (OSCA); normalization and pseudobulk chapters.
