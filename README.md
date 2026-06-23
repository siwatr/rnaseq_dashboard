# dds_dashboard

An R Shiny app for exploring, manipulating, and visualizing bulk RNA-seq results
stored in `DESeqDataSet` (`dds`) — or convertible `SingleCellExperiment` (`sce`) —
objects. Built for users with little to no bioinformatics background: QC, metadata
editing, filtering, normalization, dimensionality reduction, DESeq2 differential
expression, heatmaps, and export.

Packaged as a plain R package (`ddsdashboard`); the reusable helpers in `R/` can be
used from other projects via `library(ddsdashboard)`.

## Setup

Dependencies are managed with [mamba](https://github.com/conda-forge/miniforge)
(channels: `conda-forge` + `bioconda`).

```bash
mamba env create -f environment.yml
conda activate rnaseq_dashboard
```

## Run

```bash
R -e 'devtools::load_all("."); run_app()'
```

Or install once and launch:

```r
devtools::install_local(".")
library(ddsdashboard)
run_app()
```

## Develop

```r
devtools::load_all()    # source the package
devtools::document()    # regenerate NAMESPACE/man from roxygen
devtools::test()        # run unit tests
devtools::check()       # full R CMD check
```

See [CLAUDE.md](CLAUDE.md) for conventions, [dev_ref/rough_design.md](dev_ref/rough_design.md)
for the design narrative, and [dev_ref/roadmap.md](dev_ref/roadmap.md) for the phase plan.
