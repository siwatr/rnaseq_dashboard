test_that("Feature-info OrgDb annotation populates a name column", {
  skip_if_not_installed("DESeq2")
  skip_if_not_installed("org.Mm.eg.db")
  state <- new_app_state()
  shiny::testServer(mod_feature_server, args = list(state = state), {
    counts <- data.frame(id = c("ENSMUSG00000000001", "ENSMUSG00000000003"),
                         S1 = c(5L, 2L), S2 = c(3L, 7L))
    samples <- data.frame(sample = c("S1", "S2"), cond = c("a", "b"))
    state_load(state, ensure_feature_class(tabular_to_dds(counts, samples)),
               source = "tabular", meta = list(feature_type = "gene"))
    session$flushReact()
    v0 <- state$data_version
    session$setInputs(organism = "mouse", id_type = "ensembl")
    session$setInputs(annotate = 1)
    # Annotation lands in the draft, not state$working, until Save.
    expect_equal(state$data_version, v0)
    expect_false("gene_name" %in% colnames(SummarizedExperiment::rowData(state$working)))
    rd <- SummarizedExperiment::rowData(editor$draft())
    expect_true("gene_name" %in% colnames(rd))
    expect_equal(as.character(rd$gene_name[1]), "Gnai3")
  })
})

# Reading/trimming a GTF now lives in mod_gtf_reader (see test-mod_gtf_reader.R);
# the tests below inject the confirmed GRanges via gtf_obj() and exercise the
# draft-based annotation/length wiring on the Feature page.

test_that("Feature-info feature-unit selector updates meta without bumping version", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_feature_server, args = list(state = state), {
    state_load(state, make_mock_dds(n_genes = 10, n_per_group = 2, n_spike = 1, seed = 1),
               source = "demo", meta = list(feature_type = "gene"))
    session$flushReact()
    v0 <- state$data_version
    session$setInputs(feature_type = "transcript")
    expect_equal(state_meta(state)$feature_type, "transcript")
    expect_equal(state$data_version, v0)   # labeling change only, no data edit
  })
})

test_that("Feature-info GTF annotation and length apply to the draft until Save", {
  skip_if_not_installed("DESeq2")
  skip_if_not_installed("rtracklayer")
  skip_if_not_installed("GenomicRanges")
  state <- new_app_state()
  shiny::testServer(mod_feature_server, args = list(state = state), {
    ids <- c("ENSG00000000001", "ENSG00000000002", "ENSG00000000003")
    counts <- matrix(c(5L, 2L, 9L, 3L, 7L, 1L), nrow = 3, dimnames = list(ids, c("S1", "S2")))
    samples <- data.frame(cond = c("a", "b"), row.names = c("S1", "S2"))
    dds <- ensure_feature_class(DESeq2::DESeqDataSetFromMatrix(counts, samples, ~ 1))
    state_load(state, dds, source = "demo", meta = list(feature_type = "gene"))
    session$flushReact()
    # Inject the parsed GTF directly (skip the fileInput upload path).
    gtf_obj(import_gtf(system.file("extdata", "demo_annotation.gtf", package = "ddsdashboard")))
    v0 <- state$data_version

    # Attribute import (its own button), then length compute (its own button).
    session$setInputs(gtf_match = "auto", gtf_import = c("gene_name", "seqnames"), apply_gtf = 1)
    session$setInputs(gtf_len_type = "exon", compute_len = 1)

    # Both land in the draft; state$working is untouched until Save.
    expect_equal(state$data_version, v0)
    expect_false("feature_length" %in% colnames(SummarizedExperiment::rowData(state$working)))
    rd <- SummarizedExperiment::rowData(editor$draft())
    expect_equal(as.character(rd$gene_name), c("GeneA", "GeneB", "GeneC"))
    expect_equal(as.numeric(rd$feature_length), c(202, 201, 501))
    expect_true(has_feature_length(editor$draft()))
  })
})

test_that("Feature-info: GTF apply warns before overwrite, applies on Proceed, flags matches", {
  skip_if_not_installed("DESeq2")
  skip_if_not_installed("rtracklayer")
  skip_if_not_installed("GenomicRanges")
  state <- new_app_state()
  shiny::testServer(mod_feature_server, args = list(state = state), {
    ids <- c("ENSG00000000001", "ENSG00000000002", "ENSG00000000003")
    counts <- matrix(c(5L, 2L, 9L, 3L, 7L, 1L), nrow = 3, dimnames = list(ids, c("S1", "S2")))
    samples <- data.frame(cond = c("a", "b"), row.names = c("S1", "S2"))
    dds <- ensure_feature_class(DESeq2::DESeqDataSetFromMatrix(counts, samples, ~ 1))
    SummarizedExperiment::rowData(dds)$gene_name <- c("OLD1", "OLD2", "OLD3")  # populated -> overwrite
    state_load(state, dds, source = "demo", meta = list(feature_type = "gene"))
    session$flushReact()
    gtf_obj(import_gtf(system.file("extdata", "demo_annotation.gtf", package = "ddsdashboard")))

    session$setInputs(gtf_match = "auto", gtf_import = "gene_name", gtf_flag = TRUE, apply_gtf = 1)
    # Overwrite of an existing populated column -> deferred behind the modal.
    expect_equal(as.character(SummarizedExperiment::rowData(editor$draft())$gene_name),
                 c("OLD1", "OLD2", "OLD3"))

    session$setInputs(ow_proceed = 1)            # confirm
    rd <- SummarizedExperiment::rowData(editor$draft())
    expect_equal(as.character(rd$gene_name), c("GeneA", "GeneB", "GeneC"))  # now applied
    expect_true("in_gtf" %in% colnames(rd))      # match-flag column added
    expect_true(all(rd$in_gtf))
  })
})

test_that("Feature-info: GTF edits compose with unsaved draft edits", {
  skip_if_not_installed("DESeq2")
  skip_if_not_installed("rtracklayer")
  skip_if_not_installed("GenomicRanges")
  state <- new_app_state()
  shiny::testServer(mod_feature_server, args = list(state = state), {
    ids <- c("ENSG00000000001", "ENSG00000000002", "ENSG00000000003")
    counts <- matrix(c(5L, 2L, 9L, 3L, 7L, 1L), nrow = 3, dimnames = list(ids, c("S1", "S2")))
    samples <- data.frame(cond = c("a", "b"), row.names = c("S1", "S2"))
    dds <- ensure_feature_class(DESeq2::DESeqDataSetFromMatrix(counts, samples, ~ 1))
    state_load(state, dds, source = "demo", meta = list(feature_type = "gene"))
    session$flushReact()
    gtf_obj(import_gtf(system.file("extdata", "demo_annotation.gtf", package = "ddsdashboard")))
    # An unsaved draft edit: tag feature 1 as spike_in (mutate the draft directly,
    # as a cell edit would; state$working stays unchanged).
    editor$set(set_feature_class(editor$draft(), ids[1], "spike_in"))
    # GTF compute should preserve that pending edit.
    session$setInputs(gtf_len_type = "exon", compute_len = 1)
    rd <- SummarizedExperiment::rowData(editor$draft())
    expect_equal(as.character(rd$feature_class[1]), "spike_in")   # draft edit kept
    expect_equal(as.numeric(rd$feature_length), c(202, 201, 501)) # plus GTF length
  })
})

test_that("Feature tab: 'Set as spike-in concentration' writes spike_concentration on Save", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_feature_server, args = list(state = state), {
    dds <- make_mock_dds(n_genes = 30, n_per_group = 2, n_spike = 3, seed = 1)
    # Add a numeric column carrying a (different) concentration for spike rows.
    dds <- add_meta_column(dds, "rowData", "my_conc", "numeric", NA)
    spike <- rownames(dds)[ddsdashboard:::.detect_spike_features(dds)]
    rd <- SummarizedExperiment::rowData(dds)
    rd$my_conc[match(spike, rownames(dds))] <- c(5, 50, 500)
    SummarizedExperiment::rowData(dds) <- rd
    state_load(state, dds, source = "demo")
    session$flushReact()
    session$setInputs(spike_col = "my_conc", set_spike_conc = 1)
    # mock already has a spike_concentration column -> overwrite guard fires.
    session$setInputs(ow_proceed = 1)
    session$setInputs(`editor-save` = 1)            # Save lives in the editor submodule
    sc <- SummarizedExperiment::rowData(state$working)$spike_concentration
    expect_equal(sc[match(spike, rownames(state$working))], c(5, 50, 500))
  })
})
