# Shiny App for exploring and basic processing of RNA-seq output

Hi, I'm planning to make a R shiny app for exploring the results obtained from single-cell RNA-seq analysis pipeline. The input for this app can be an R object (supporting DESeqDataSet (dds) and SingleCellExperiment (sce)).


## Goal 

The aim is to allow people with little to no bioinformatic knowledge to explore their datasets, as well as manipulating and performing basic analysis on these input object. 

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


#### Manipulating feature (gene) metadata

We expected the input to have many genes (more than 50,000 genes in the case of mouse genome), thus manual manipulation of `rowData` can be challenging. User can upload the GTF annotation of their datasets (read via `rtracklayer` package). They can select what feature type should be used for their dataset, which column corresponding to their feature names (`rownames(dds)`), then they can select some columns from `mcols(gtf)` that will be incorporated in the feature metadata (`rowData(dds)` or `rowRanges(dds)`). 
The feature length (important for TPM or FPKM normalization) can optionally be determined here based on certain column of `rowData(dds)` (some quantification tools report these values and it might be coming with the input already), alternatively it can be determined from the GTF that is read in the current session. 

Certain group of genes (explicitly specified or narrow down using `dplyr::filter(as.data.frame(mcols(gtf)))`) can be marked as exogenous genes (e.g., being over-expressed) or spike-in genes (check ERCC-spike-in). Note that I have little experience with UMI, so I don't know if we should also handle it here as well?

__Note__ on feature length, the effective feature length may not comming from the same type of feature represented in the datasets! For example, if gene ID or gene names is used as the feature, but quantification of read count was done only on exon, then the effective feature length used for TPM and FPKM normalization should come from the exonic part of that gene.

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


#### Sample filtering

The app can suggest which samples should be discarded, and for what reason (e.g., those suggested by `scater`). Users can also decide to adopt our suggestion, and/or manually select additional samples to be excluded from the dataset. 

#### Feature filtering

Feature filtering can be done in semi-automating manner. features with zero counts in all samples should be removed. Features with low counts (either imposed by hard number of `rowSums(count_matrix) < x`; where `x` is a hard threshold set up user, or some other automatic way (please suggest)) can also be removed.

__Note__ that the function from the package `HTSFilter` can take long time to run, so we will avoid using it.


<!--
MARK: 3rd page 
-->
### The third page: Processing

The process in this step can be somewhat optional. User should be able to add new assays, such as TPM, FPKM, or CPM, to the dataset.


<!--
MARK: 4th page 
-->
### The fourth page: Dimensional Reduction and gene expression:

The fourth page can be used for exploring dimensional reduction plots (hereafter _ReduceDim_ plot for short), such as PCA, tSNE and UMAP, that constructed from one of the following methods:
* PCA extracted from DESeq2 function
* Manual plot of PCA, tSNE, and UMAP.

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
  * log2-fold-change: default log2(2) -> `abs(log2FoldChange) <= log2(2)`
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

<!--
MARK: 6th page 
-->
### The sixth page: Data export

This is for exporting the dds object, DESeq2 result tables, and plots.


<!-- 
#region MARK: Roadmap
 -->

## Road map

The project can be implemented in different phases. Here is my suggestion, you are free to comment or suggest on how to do this differently.

### Phase 1

* Implement basic shiny app layout
* Allow data import (1st page), show stats of the input (`show(dds)`)
* Set up data export page

### Phase 2
* Metadata manipulation
* Implement QC
* Implement sample/feature filtering
* Adding normalization assays


### Phase 3:
* Making dimensional reduction visualization page

### Phase 4:
* Implement DESeq2-related steps
* Implement DESeq2 plots

<!-- 
#region MARK: Roadmap
 -->

## misc

* Consider using `ggplot2` package as the main plot package. `ComplexHeatmap` is preferred as main package for heatmap.
* Since showing the plot can take a while, we can add a button that will trigger re-making of the plot when user is satisfied with their setting.
* logcounts assay can be added automatically to the dds object. Consider doing this with CPM as well, but this has to be updated each time the new features is added or removed.
* Consider making git commit each time you add a considerable, meaningful change to the repo.

