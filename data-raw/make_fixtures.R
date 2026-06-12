# Materialize demo/test fixtures from the deterministic generators in
# R/mock_data.R. Run from the repo root in the project env:
#   mamba run -n rnaseq_dashboard Rscript data-raw/make_fixtures.R
#
# Tests call make_mock_dds() directly (seeded), so they don't need these files;
# the .rds here backs a future "load demo data" button in the app.

devtools::load_all(".")

dir.create("inst/extdata", recursive = TRUE, showWarnings = FALSE)

dds <- make_mock_dds(n_genes = 200, n_per_group = 4, n_spike = 10, seed = 1)
saveRDS(dds, "inst/extdata/mock_dds.rds")
message("Wrote inst/extdata/mock_dds.rds  (",
        nrow(dds), " features x ", ncol(dds), " samples)")
