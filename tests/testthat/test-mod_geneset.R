# Gene Sets page server (R/mod_geneset.R) -- the "Build a set" staging workflow (P6b-2).

test_that("mod_geneset UI mounts the build + your-sets cards", {
  ui <- as.character(mod_geneset_ui("gs"))
  expect_match(ui, "gs-source")            # source pills
  expect_match(ui, "gs-paste_q")           # shared gene-search box
  expect_match(ui, "gs-preview_table")     # staged preview
  expect_match(ui, "gs-save_mode")         # New / Add pills
  expect_match(ui, "gs-new_name")
  expect_match(ui, "gs-sets_table")        # your gene sets
  expect_match(ui, "gs-members_table")
})

# Stage a paste of rownames + create/append via the Save section.
stage_paste <- function(session, q) {
  session$setInputs(source = "paste", paste_searchby = "__rownames__",
                    paste_mode = "exact", paste_ci = FALSE, paste_q = q)
  session$elapse(300); session$flushReact()
}

test_that("Paste -> New creates a set; a name clash is rejected (no auto-suffix)", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_geneset_server, args = list(state = state), {
    dds <- ensure_logcounts(make_mock_dds(n_genes = 20, n_per_group = 3, n_spike = 0, seed = 1))
    state_load(state, dds, source = "demo", meta = list(feature_type = "gene"))
    session$flushReact()
    dv <- state$data_version
    r1 <- rownames(state$working)[1]; r2 <- rownames(state$working)[2]; r3 <- rownames(state$working)[3]

    stage_paste(session, paste(r1, r2, sep = ", "))
    expect_setequal(names(staged()), "Pasted genes")
    expect_setequal(staged()[["Pasted genes"]], c(r1, r2))

    session$setInputs(save_mode = "new", new_name = "mySet")
    session$setInputs(create = 1); session$flushReact()
    expect_named(state$gene_sets, "mySet")
    expect_setequal(state$gene_sets$mySet$ids, c(r1, r2))
    expect_equal(state$data_version, dv)                 # not a data edit

    # A clash is REJECTED (no mySet_2), unlike the old auto-suffix.
    stage_paste(session, r3)
    session$setInputs(new_name = "mySet"); session$setInputs(create = 2); session$flushReact()
    expect_named(state$gene_sets, "mySet")               # still just one
    expect_false("mySet_2" %in% names(state$gene_sets))
    expect_setequal(state$gene_sets$mySet$ids, c(r1, r2))
  })
})

test_that("Add-to-existing unions the staged genes into selected sets", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_geneset_server, args = list(state = state), {
    dds <- ensure_logcounts(make_mock_dds(n_genes = 20, n_per_group = 3, n_spike = 0, seed = 2))
    state_load(state, dds, source = "demo", meta = list(feature_type = "gene"))
    session$flushReact()
    r1 <- rownames(state$working)[1]; r2 <- rownames(state$working)[2]; r3 <- rownames(state$working)[3]

    stage_paste(session, r1)
    session$setInputs(save_mode = "new", new_name = "A"); session$setInputs(create = 1); session$flushReact()
    expect_setequal(state$gene_sets$A$ids, r1)

    stage_paste(session, paste(r2, r3, sep = ", "))
    session$setInputs(save_mode = "add", add_targets = "A"); session$setInputs(add_existing = 1)
    session$flushReact()
    expect_setequal(state$gene_sets$A$ids, c(r1, r2, r3))   # unioned
  })
})

test_that("literal-add only applies when searching by feature ID", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_geneset_server, args = list(state = state), {
    dds <- ensure_logcounts(make_mock_dds(n_genes = 20, n_per_group = 3, n_spike = 0, seed = 3))
    state_load(state, dds, source = "demo", meta = list(feature_type = "gene"))
    session$flushReact()
    r1 <- rownames(state$working)[1]

    # Search by gene_name with literal ON: an unmatched token must NOT be added
    # literally (it would be a name, not an id).
    session$setInputs(source = "paste", paste_searchby = "gene_name", paste_mode = "exact",
                      paste_ci = FALSE, paste_literal = TRUE, paste_q = "GHOST")
    session$elapse(300); session$flushReact()
    expect_length(staged(), 0L)                          # GHOST not added

    # Search by feature ID with literal ON: the literal id IS added.
    session$setInputs(paste_searchby = "__rownames__", paste_q = paste(r1, "GHOST_ID", sep = ", "))
    session$elapse(300); session$flushReact()
    expect_setequal(staged()[["Pasted genes"]], c(r1, "GHOST_ID"))

    # ...but only in EXACT mode: in contains/regex an unmatched entry is a
    # PATTERN, not an id, so it must never be committed literally.
    session$setInputs(paste_mode = "regex", paste_q = "^ZZZ_no_match.*")
    session$elapse(300); session$flushReact()
    expect_length(staged(), 0L)
    session$setInputs(paste_mode = "contains", paste_q = "ZZZ_no_match")
    session$elapse(300); session$flushReact()
    expect_length(staged(), 0L)
  })
})

test_that("'Use shrunk LFC' falls back to standard LFCs when shrinkage never ran", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_geneset_server, args = list(state = state), {
    dds <- ensure_logcounts(make_mock_dds(n_genes = 20, n_per_group = 3, n_spike = 0, seed = 7))
    state_load(state, dds, source = "demo", meta = list(feature_type = "gene"))
    session$flushReact()
    rn <- rownames(state$working)
    # A results frame as de_results() builds it with shrink = "none": the
    # log2FoldChange_shrunk column EXISTS but is all-NA.
    df <- data.frame(baseMean = 100,
                     padj = c(0.001, 0.001, rep(0.9, length(rn) - 2)),
                     log2FoldChange = c(3, -3, rep(0, length(rn) - 2)),
                     log2FoldChange_shrunk = NA_real_, row.names = rn)
    de <- state$de; de$results <- list("C1" = df); de$active <- "C1"; state$de <- de
    session$flushReact()

    session$setInputs(source = "deg", deg_contrast = "C1", deg_dir = "both",
                      deg_padj = 0.05, deg_lfc = 1, deg_shrunk = TRUE)
    session$flushReact()
    expect_false(deg_shrunk_ok())          # column present but no real values
    expect_equal(deg_col(), "DEG")         # falls back instead of an empty set
    expect_setequal(staged()[[1]], rn[1:2])
  })
})

test_that("DE DEGs + top-variable sources stage the right ids", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_geneset_server, args = list(state = state), {
    dds <- ensure_logcounts(make_mock_dds(n_genes = 30, n_per_group = 3, n_spike = 4, seed = 4))
    state_load(state, dds, source = "demo", meta = list(feature_type = "gene"))
    session$flushReact()
    rn <- rownames(state$working)
    df <- data.frame(baseMean = 100,
                     padj = c(0.001, 0.001, 0.001, 0.001, rep(0.9, length(rn) - 4)),
                     log2FoldChange = c(3, 2, -3, -2, rep(0, length(rn) - 4)), row.names = rn)
    de <- state$de; de$results <- list("C1" = df); de$active <- "C1"; state$de <- de
    session$flushReact()

    session$setInputs(source = "deg", deg_contrast = "C1", deg_dir = "up",
                      deg_padj = 0.05, deg_lfc = 1, deg_shrunk = FALSE); session$flushReact()
    expect_setequal(staged()[[1]], rn[1:2])
    session$setInputs(save_mode = "new", new_name = "up"); session$setInputs(create = 1); session$flushReact()
    expect_setequal(state$gene_sets$up$ids, rn[1:2])
    expect_match(state$gene_sets$up$source, "^DE: C1")

    session$setInputs(source = "topvar", topvar_n = 5); session$flushReact()
    expect_length(staged()[[1]], 5L)
    expect_false(any(grepl("^ERCC", staged()[[1]])))     # endogenous only
  })
})

test_that("membership is non-destructive across a feature drop; load clears sets", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_geneset_server, args = list(state = state), {
    dds <- ensure_logcounts(make_mock_dds(n_genes = 20, n_per_group = 3, n_spike = 0, seed = 5))
    state_load(state, dds, source = "demo", meta = list(feature_type = "gene"))
    session$flushReact()
    r1 <- rownames(state$working)[1]; r2 <- rownames(state$working)[2]

    stage_paste(session, paste(r1, r2, sep = ", "))
    session$setInputs(save_mode = "new", new_name = "s"); session$setInputs(create = 1); session$flushReact()
    expect_setequal(state$gene_sets$s$ids, c(r1, r2))

    state_mutate(state, function(d) d[-1, ], action = list(action = "drop_feature"))
    session$flushReact()
    expect_setequal(state$gene_sets$s$ids, c(r1, r2))               # intact
    expect_equal(gene_set_absent(state$gene_sets$s, rownames(state$working)), r1)

    state_load(state, dds, source = "demo", meta = list(feature_type = "gene"))
    session$flushReact()
    expect_length(state$gene_sets, 0L)
  })
})

test_that("rename / delete act on the selected row; members follow selection", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_geneset_server, args = list(state = state), {
    dds <- ensure_logcounts(make_mock_dds(n_genes = 20, n_per_group = 3, n_spike = 0, seed = 6))
    state_load(state, dds, source = "demo", meta = list(feature_type = "gene"))
    session$flushReact()
    r1 <- rownames(state$working)[1]
    stage_paste(session, r1)
    session$setInputs(save_mode = "new", new_name = "orig"); session$setInputs(create = 1); session$flushReact()

    # No selection -> members placeholder (selected_name NULL).
    expect_null(selected_name())
    session$setInputs(sets_table_rows_selected = 1); session$flushReact()
    expect_equal(selected_name(), "orig")

    session$setInputs(rename_to = "renamed"); session$setInputs(rename_ok = 1); session$flushReact()
    expect_named(state$gene_sets, "renamed")

    session$setInputs(sets_table_rows_selected = 1); session$flushReact()
    session$setInputs(delete = 1); session$flushReact()
    expect_length(state$gene_sets, 0L)
  })
})
