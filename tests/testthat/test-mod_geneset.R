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

# --- Table import (P6c) -----------------------------------------------------

# Write a DESeq2-result-shaped CSV and hand it to the module's fileInput.
.write_import_csv <- function(rn) {
  path <- withr::local_tempfile(fileext = ".csv", .local_envir = parent.frame())
  utils::write.csv(
    data.frame(gene = c(rn[1:2], rn[3:4], "GHOST_ID"),
               direction = c("up", "up", "down", "down", "up"),
               padj = 0.001, stringsAsFactors = FALSE),
    path, row.names = FALSE)
  path
}

test_that("table import stages one set; ID column resolves; unmatched excluded", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_geneset_server, args = list(state = state), {
    dds <- ensure_logcounts(make_mock_dds(n_genes = 20, n_per_group = 3, n_spike = 0, seed = 8))
    state_load(state, dds, source = "demo", meta = list(feature_type = "gene"))
    session$flushReact()
    rn <- rownames(state$working)
    path <- .write_import_csv(rn)

    session$setInputs(source = "file",
                      tbl_file = list(datapath = path, name = "res.csv"),
                      tbl_header = TRUE, tbl_id_col = "gene", tbl_match_by = "__rownames__",
                      tbl_rows = "view", tbl_anno_cols = character(0),
                      tbl_sep = ".", tbl_literal = FALSE)
    session$flushReact()
    expect_null(tbl_raw())                          # nothing loaded until the Load click
    session$setInputs(tbl_load = 1); session$flushReact()
    expect_equal(nrow(tbl_raw()), 5L)
    # No split -> ONE staged set; the unmatched GHOST_ID is dropped.
    expect_length(staged(), 1L)
    expect_setequal(staged()[[1]], rn[1:4])

    # Literal ON (matching against feature ids) keeps the unmatched id.
    session$setInputs(tbl_literal = TRUE); session$flushReact()
    expect_setequal(staged()[[1]], c(rn[1:4], "GHOST_ID"))
  })
})

test_that("annotation-split stages N sets; multi-save creates one set each", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_geneset_server, args = list(state = state), {
    dds <- ensure_logcounts(make_mock_dds(n_genes = 20, n_per_group = 3, n_spike = 0, seed = 9))
    state_load(state, dds, source = "demo", meta = list(feature_type = "gene"))
    session$flushReact()
    rn <- rownames(state$working)
    path <- .write_import_csv(rn)

    session$setInputs(source = "file",
                      tbl_file = list(datapath = path, name = "res.csv"),
                      tbl_header = TRUE, tbl_id_col = "gene", tbl_match_by = "__rownames__",
                      tbl_rows = "view", tbl_literal = FALSE,
                      tbl_anno_cols = "direction", tbl_sep = ".")
    session$setInputs(tbl_load = 1); session$flushReact()
    # Split by direction -> two staged sets.
    expect_setequal(names(staged()), c("up", "down"))
    expect_setequal(staged()$up, rn[1:2])
    expect_setequal(staged()$down, rn[3:4])
    expect_equal(output$multi_staged, "1")        # the Save section swaps to multi

    # Multi-save: create one session set per staged set, with a prefix.
    session$setInputs(multi_pick = c("up", "down"), multi_mode = "new",
                      multi_prefix = "res_", multi_autoname = FALSE)
    session$flushReact()
    session$setInputs(multi_create = 1); session$flushReact()
    expect_setequal(names(state$gene_sets), c("res_up", "res_down"))
    expect_setequal(state$gene_sets$res_up$ids, rn[1:2])
    expect_match(state$gene_sets$res_up$source, "^import: res.csv")
  })
})

test_that("multi-save blocks name clashes unless auto-rename is on", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_geneset_server, args = list(state = state), {
    dds <- ensure_logcounts(make_mock_dds(n_genes = 20, n_per_group = 3, n_spike = 0, seed = 10))
    state_load(state, dds, source = "demo", meta = list(feature_type = "gene"))
    session$flushReact()
    rn <- rownames(state$working)
    path <- .write_import_csv(rn)
    session$setInputs(source = "file",
                      tbl_file = list(datapath = path, name = "res.csv"),
                      tbl_header = TRUE, tbl_id_col = "gene", tbl_match_by = "__rownames__",
                      tbl_rows = "view", tbl_literal = FALSE,
                      tbl_anno_cols = "direction", tbl_sep = ".",
                      multi_pick = c("up", "down"), multi_mode = "new",
                      multi_prefix = "", multi_autoname = FALSE)
    session$setInputs(tbl_load = 1); session$flushReact()
    session$setInputs(multi_create = 1); session$flushReact()
    expect_setequal(names(state$gene_sets), c("up", "down"))

    # Re-creating the same names is BLOCKED while auto-rename is off.
    session$setInputs(multi_create = 2); session$flushReact()
    expect_setequal(names(state$gene_sets), c("up", "down"))

    # With auto-rename on, the clash resolves via the suffix path.
    session$setInputs(multi_autoname = TRUE); session$flushReact()
    session$setInputs(multi_create = 3); session$flushReact()
    expect_setequal(names(state$gene_sets), c("up", "down", "up_2", "down_2"))

    # Multi 'Add to existing' unions ALL picked staged sets into the targets.
    session$setInputs(multi_mode = "add", multi_targets = "up"); session$flushReact()
    session$setInputs(multi_add = 1); session$flushReact()
    expect_setequal(state$gene_sets$up$ids, rn[1:4])
  })
})

test_that("table import honours the row scope (filtered view vs selected)", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_geneset_server, args = list(state = state), {
    dds <- ensure_logcounts(make_mock_dds(n_genes = 20, n_per_group = 3, n_spike = 0, seed = 11))
    state_load(state, dds, source = "demo", meta = list(feature_type = "gene"))
    session$flushReact()
    rn <- rownames(state$working)
    path <- .write_import_csv(rn)
    session$setInputs(source = "file",
                      tbl_file = list(datapath = path, name = "res.csv"),
                      tbl_header = TRUE, tbl_id_col = "gene", tbl_match_by = "__rownames__",
                      tbl_anno_cols = character(0), tbl_sep = ".", tbl_literal = FALSE)
    session$setInputs(tbl_load = 1); session$flushReact()

    # A filtered view (DT reports the surviving rows) drives the staged set.
    session$setInputs(tbl_rows = "view", tbl_table_rows_all = c(1L, 2L)); session$flushReact()
    expect_setequal(staged()[[1]], rn[1:2])

    # Explicit row selection instead.
    session$setInputs(tbl_rows = "selected", tbl_table_rows_selected = c(3L)); session$flushReact()
    expect_setequal(staged()[[1]], rn[3])
  })
})

test_that("header toggle: no-header files get Column_<i> names + auto-detected match field", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_geneset_server, args = list(state = state), {
    dds <- ensure_logcounts(make_mock_dds(n_genes = 20, n_per_group = 3, n_spike = 0, seed = 12))
    state_load(state, dds, source = "demo", meta = list(feature_type = "gene"))
    session$flushReact()
    rn <- rownames(state$working)
    # A header-less single-column list of feature ids.
    path <- withr::local_tempfile(fileext = ".txt")
    writeLines(rn[1:5], path)

    session$setInputs(source = "file", tbl_file = list(datapath = path, name = "ids.txt"),
                      tbl_header = FALSE, tbl_load = 1)
    session$flushReact()
    df <- tbl_raw()
    expect_equal(names(df), "Column_1")            # friendly name, not readr's X1
    expect_equal(nrow(df), 5L)
    # Auto-detect picks feature ids (rownames) for this column.
    expect_equal(.gs_best_match_field(df$Column_1, state$working), "__rownames__")

    session$setInputs(tbl_id_col = "Column_1", tbl_match_by = "__rownames__",
                      tbl_rows = "view"); session$flushReact()
    expect_setequal(staged()[[1]], rn[1:5])
  })
})

test_that("auto-detect matches a gene_name column when the ID column holds names", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_geneset_server, args = list(state = state), {
    dds <- ensure_logcounts(make_mock_dds(n_genes = 20, n_per_group = 3, n_spike = 0, seed = 13))
    state_load(state, dds, source = "demo", meta = list(feature_type = "gene"))
    session$flushReact()
    nm <- as.character(SummarizedExperiment::rowData(state$working)$gene_name)[1:5]
    # A column of gene NAMES should auto-detect the gene_name field, not rownames.
    expect_equal(.gs_best_match_field(nm, state$working), "gene_name")
    # And an id column still auto-detects rownames.
    expect_equal(.gs_best_match_field(rownames(state$working)[1:5], state$working), "__rownames__")
  })
})

test_that("a name matching several features keeps ALL by default, or the first on request", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_geneset_server, args = list(state = state), {
    dds <- make_mock_dds(n_genes = 20, n_per_group = 3, n_spike = 0, seed = 14)
    rd <- SummarizedExperiment::rowData(dds)
    rd$gene_name[2] <- rd$gene_name[1]                    # a duplicated gene name
    SummarizedExperiment::rowData(dds) <- rd
    state_load(state, ensure_logcounts(dds), source = "demo", meta = list(feature_type = "gene"))
    session$flushReact()
    rn <- rownames(state$working); dup <- as.character(rd$gene_name[1])

    path <- withr::local_tempfile(fileext = ".csv")
    utils::write.csv(data.frame(name = dup, grp = "a"), path, row.names = FALSE)
    session$setInputs(source = "file", tbl_file = list(datapath = path, name = "n.csv"),
                      tbl_header = TRUE, tbl_load = 1, tbl_id_col = "name",
                      tbl_match_by = "gene_name", tbl_rows = "view",
                      tbl_anno_cols = character(0))
    session$flushReact()
    # Keep-all (default): the ambiguous name expands to BOTH feature ids.
    expect_setequal(staged()[[1]], c(rn[1], rn[2]))

    # First-only: just the first matching feature.
    session$setInputs(tbl_multi = "first"); session$flushReact()
    expect_equal(staged()[[1]], rn[1])
  })
})

test_that("the header toggle re-reads reactively after the first Load (no re-click)", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_geneset_server, args = list(state = state), {
    dds <- ensure_logcounts(make_mock_dds(n_genes = 20, n_per_group = 3, n_spike = 0, seed = 15))
    state_load(state, dds, source = "demo", meta = list(feature_type = "gene"))
    session$flushReact()
    rn <- rownames(state$working)
    path <- withr::local_tempfile(fileext = ".csv")
    # No real header: the first data row would be swallowed as one if header = TRUE.
    writeLines(c(rn[1], rn[2], rn[3]), path)

    session$setInputs(source = "file", tbl_file = list(datapath = path, name = "ids.csv"),
                      tbl_header = TRUE, tbl_load = 1)
    session$flushReact()
    expect_equal(nrow(tbl_raw()), 2L)                    # first line taken as the header
    # Flip the toggle -- NO second Load click -- and it re-parses.
    session$setInputs(tbl_header = FALSE); session$flushReact()
    expect_equal(nrow(tbl_raw()), 3L)
    expect_equal(names(tbl_raw()), "Column_1")
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
