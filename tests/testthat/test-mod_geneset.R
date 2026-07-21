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
  expect_match(ui, "gs-anno_members_ui")   # Annotation tab builder
  expect_match(ui, "gs-anno_build")
  expect_match(ui, "gs-anno_panels")
})

test_that("Annotation tab builds a kind='annotated' record and deletes it", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_geneset_server, args = list(state = state), {
    dds <- ensure_logcounts(make_mock_dds(n_genes = 20, n_per_group = 3, n_spike = 0, seed = 5))
    state_load(state, dds, source = "demo", meta = list(feature_type = "gene"))
    session$flushReact()
    rn <- rownames(state$working)
    # Two simple member sets with a one-gene overlap (rn[2]).
    state$gene_sets <- list(up   = new_gene_set(rn[1:2], source = "paste"),
                            down = new_gene_set(rn[2:3], source = "paste"))
    session$flushReact()
    expect_setequal(simple_names(), c("up", "down"))

    session$setInputs(anno_members = c("up", "down"), anno_shared = "concat",
                      anno_sep = ";", anno_name = "DE dir")
    session$flushReact()
    expect_true(anno_combined()$ok)
    session$setInputs(anno_build = 1); session$flushReact()
    rec <- state$gene_sets[["DE dir"]]
    expect_equal(rec$kind, "annotated")
    expect_equal(unname(rec$annotation[rn[2]]), "up;down")     # overlap concatenated
    expect_setequal(annotated_names(), "DE dir")
    expect_match(rec$source, "^combine: up, down")

    # A name clash is rejected (no second record, store unchanged).
    session$setInputs(anno_name = "DE dir", anno_build = 2); session$flushReact()
    expect_length(annotated_names(), 1L)

    # Delete via the sidebar selected-delete path.
    session$setInputs(anno_del_pick = "DE dir", anno_del_sel = 1); session$flushReact()
    session$setInputs(anno_del_sel_ok = 1); session$flushReact()
    expect_length(annotated_names(), 0L)
    # The simple member sets are untouched.
    expect_setequal(names(state$gene_sets), c("up", "down"))
  })
})

test_that("annotation rename (set + level) propagates to the Palette geneset config", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_geneset_server, args = list(state = state), {
    dds <- ensure_logcounts(make_mock_dds(n_genes = 12, n_per_group = 3, n_spike = 0, seed = 6))
    state_load(state, dds, source = "demo", meta = list(feature_type = "gene"))
    session$flushReact()
    rn <- rownames(state$working)
    state$gene_sets <- list(A = new_gene_set(rn[1:4], kind = "annotated",
      annotation = stats::setNames(c("up", "up", "down", "down"), rn[1:4])))
    state$palette <- list(geneset = list(A = list(name = "Custom palette",
      colors = c(up = "#FF0000", down = "#0000FF"))))
    session$flushReact()

    # Rename set A -> B via the panel modal (k = 1): store + palette key follow.
    session$setInputs(anno_rename_to_1 = "B", anno_rename_ok_1 = 1); session$flushReact()
    expect_null(state$gene_sets[["A"]])
    expect_false(is.null(state$gene_sets[["B"]]))
    expect_equal(state$palette$geneset[["B"]]$colors[["up"]], "#FF0000")
    expect_null(state$palette$geneset[["A"]])

    # Rename level up -> UP: annotation map + palette colour key follow.
    session$setInputs(anno_level_from_1 = "up", anno_level_to_1 = "UP", anno_level_ok_1 = 1)
    session$flushReact()
    expect_equal(unname(state$gene_sets$B$annotation[rn[1]]), "UP")
    expect_equal(state$palette$geneset$B$colors[["UP"]], "#FF0000")
    expect_false("up" %in% names(state$palette$geneset$B$colors))

    # Merge guard: renaming UP -> down (an existing level) is refused.
    session$setInputs(anno_level_from_1 = "UP", anno_level_to_1 = "down", anno_level_ok_1 = 2)
    session$flushReact()
    expect_true(all(c("UP", "down") %in% unique(unname(state$gene_sets$B$annotation))))
  })
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

    session$setInputs(source = "deg", deg_contrast = "C1", deg_dir = c("up", "down"),
                      deg_padj = 0.05, deg_lfc = 1, deg_shrunk = TRUE)
    session$flushReact()
    df1 <- deg_classify_one("C1")
    expect_false(deg_shrunk_avail(df1))    # column present but no real values
    expect_equal(deg_col_for(df1), "DEG")  # falls back instead of an empty set
    # Group build: one set per direction, named "<contrast> <dir>".
    st <- staged()
    expect_setequal(names(st), c("C1 up", "C1 down"))
    expect_setequal(unlist(st, use.names = FALSE), rn[1:2])
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

test_that("DE group build stages one set per (contrast x direction) and multi-saves", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_geneset_server, args = list(state = state), {
    dds <- ensure_logcounts(make_mock_dds(n_genes = 30, n_per_group = 3, n_spike = 0, seed = 9))
    state_load(state, dds, source = "demo", meta = list(feature_type = "gene"))
    session$flushReact()
    rn <- rownames(state$working)
    mk <- function(up, down) data.frame(
      baseMean = 100,
      padj = c(rep(0.001, length(up) + length(down)), rep(0.9, length(rn) - length(up) - length(down))),
      log2FoldChange = c(rep(3, length(up)), rep(-3, length(down)), rep(0, length(rn) - length(up) - length(down))),
      row.names = c(up, down, setdiff(rn, c(up, down))))
    de <- state$de
    de$results <- list(C1 = mk(rn[1:2], rn[3:4]), C2 = mk(rn[5:7], rn[8]))
    de$active <- "C1"; state$de <- de
    session$flushReact()

    session$setInputs(source = "deg", deg_contrast = c("C1", "C2"),
                      deg_dir = c("up", "down"), deg_padj = 0.05, deg_lfc = 1,
                      deg_shrunk = FALSE)
    session$flushReact()
    st <- staged()
    expect_setequal(names(st), c("C1 up", "C1 down", "C2 up", "C2 down"))
    expect_setequal(st[["C1 up"]], rn[1:2])
    expect_setequal(st[["C2 up"]], rn[5:7])
    expect_setequal(st[["C2 down"]], rn[8])

    # >1 staged set -> the multi-set Save path. Create all four under a prefix.
    session$setInputs(multi_pick = names(st), multi_prefix = "grp_",
                      multi_autoname = FALSE)
    session$flushReact()
    session$setInputs(multi_create = 1); session$flushReact()
    expect_true(all(c("grp_C1 up", "grp_C1 down", "grp_C2 up", "grp_C2 down") %in%
                      names(state$gene_sets)))
    expect_setequal(state$gene_sets[["grp_C2 up"]]$ids, rn[5:7])
    expect_match(state$gene_sets[["grp_C1 up"]]$source, "^DE: C1, C2")
  })
})

test_that("DE group build resolves the shrunk column PER contrast (mixed availability)", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_geneset_server, args = list(state = state), {
    dds <- ensure_logcounts(make_mock_dds(n_genes = 30, n_per_group = 3, n_spike = 0, seed = 12))
    state_load(state, dds, source = "demo", meta = list(feature_type = "gene"))
    session$flushReact()
    rn <- rownames(state$working)
    mk <- function(shrunk_na) data.frame(
      baseMean = 100,
      padj = c(0.001, 0.001, rep(0.9, length(rn) - 2)),
      log2FoldChange = c(3, -3, rep(0, length(rn) - 2)),
      log2FoldChange_shrunk = if (shrunk_na) NA_real_ else c(3, -3, rep(0, length(rn) - 2)),
      row.names = rn)
    de <- state$de
    de$results <- list(C1 = mk(FALSE), C2 = mk(TRUE))   # C2's shrunk col is all-NA
    de$active <- "C1"; state$de <- de
    session$flushReact()
    session$setInputs(source = "deg", deg_contrast = c("C1", "C2"),
                      deg_dir = c("up", "down"), deg_padj = 0.05, deg_lfc = 1,
                      deg_shrunk = TRUE)
    session$flushReact()
    st <- staged()
    # C2 must NOT be silently emptied by a global shrunk choice off C1 -- it falls
    # back to standard DEG per contrast, so its up/down sets are populated.
    expect_setequal(names(st), c("C1 up", "C1 down", "C2 up", "C2 down"))
    expect_setequal(st[["C2 up"]], rn[1]); expect_setequal(st[["C2 down"]], rn[2])
    expect_setequal(st[["C1 up"]], rn[1]); expect_setequal(st[["C1 down"]], rn[2])
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

test_that("import fast-track combines split groups into one annotated set", {
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
                      tbl_rows = "view", tbl_literal = FALSE,
                      tbl_anno_cols = "direction", tbl_sep = ".")
    session$setInputs(tbl_load = 1); session$flushReact()
    # Fast-track: one annotated set, annotation = the split-group value.
    session$setInputs(tbl_anno_name = "DE dir", tbl_anno_build = 1); session$flushReact()
    rec <- state$gene_sets[["DE dir"]]
    expect_equal(rec$kind, "annotated")
    expect_equal(unname(rec$annotation[rn[1]]), "up")
    expect_equal(unname(rec$annotation[rn[3]]), "down")
    expect_match(rec$source, "^import annotation: res.csv")
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

test_that("gene-set file import stages the file's named sets, then commits them", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_geneset_server, args = list(state = state), {
    dds <- ensure_logcounts(make_mock_dds(n_genes = 20, n_per_group = 3, n_spike = 0, seed = 16))
    state_load(state, dds, source = "demo", meta = list(feature_type = "gene"))
    session$flushReact()
    rn <- rownames(state$working)

    # Export a two-set store to a GMT file, then import it back.
    path <- withr::local_tempfile(fileext = ".gmt")
    writeLines(gene_sets_to_gmt(list(SetA = new_gene_set(rn[1:3], source = "paste"),
                                     SetB = new_gene_set(rn[4:5], source = "paste"))), path)

    session$setInputs(source = "gsfile",
                      gsfile_file = list(datapath = path, name = "mysets.gmt"))
    session$setInputs(gsfile_load = 1); session$flushReact()
    expect_setequal(names(staged()), c("SetA", "SetB"))
    expect_setequal(staged()$SetA, rn[1:3])
    expect_equal(output$multi_staged, "1")            # file import -> multi-set Save

    session$setInputs(multi_pick = c("SetA", "SetB"), multi_mode = "new",
                      multi_prefix = "", multi_autoname = FALSE)
    session$flushReact()
    session$setInputs(multi_create = 1); session$flushReact()
    expect_setequal(names(state$gene_sets), c("SetA", "SetB"))
    expect_setequal(state$gene_sets$SetA$ids, rn[1:3])
    expect_match(state$gene_sets$SetA$source, "^import: mysets.gmt")

    # A single-set file still routes to the multi-set Save (names come from it).
    p2 <- withr::local_tempfile(fileext = ".json")
    writeLines(gene_sets_to_json(list(Solo = new_gene_set(rn[6], source = "paste"))), p2)
    session$setInputs(gsfile_file = list(datapath = p2, name = "solo.json"))
    session$setInputs(gsfile_load = 1); session$flushReact()
    expect_equal(names(staged()), "Solo")
    expect_equal(output$multi_staged, "1")
  })
})

test_that("gene-set file import matches by a name column and honours keep-unmatched", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_geneset_server, args = list(state = state), {
    dds <- ensure_logcounts(make_mock_dds(n_genes = 20, n_per_group = 3, n_spike = 0, seed = 18))
    state_load(state, dds, source = "demo", meta = list(feature_type = "gene"))
    session$flushReact()
    rn <- rownames(state$working)
    gn <- as.character(SummarizedExperiment::rowData(state$working)$gene_name)

    # A file authored with gene NAMES + one name absent from the dataset.
    path <- withr::local_tempfile(fileext = ".gmt")
    writeLines(gene_sets_to_gmt(list(
      NameSet = new_gene_set(c(gn[1:3], "MISSING_NAME"), source = "paste"))), path)
    session$setInputs(source = "gsfile",
                      gsfile_file = list(datapath = path, name = "names.gmt"))
    session$setInputs(gsfile_load = 1); session$flushReact()

    # Matching against gene_name resolves the names to feature ids; the absent
    # name is dropped by default (literal off).
    session$setInputs(gsfile_match_by = "gene_name", gsfile_multi = "all",
                      gsfile_literal = FALSE)
    session$flushReact()
    expect_setequal(staged()$NameSet, rn[1:3])
    expect_equal(gsfile_resolved()$unmatched, "MISSING_NAME")

    # Keep-unmatched adds the absent id verbatim (non-destructive authorship).
    session$setInputs(gsfile_literal = TRUE); session$flushReact()
    expect_setequal(staged()$NameSet, c(rn[1:3], "MISSING_NAME"))
  })
})

test_that("the keep-unmatched toggle renders in the Preview accordion per source", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_geneset_server, args = list(state = state), {
    dds <- ensure_logcounts(make_mock_dds(n_genes = 20, n_per_group = 3, n_spike = 0, seed = 19))
    state_load(state, dds, source = "demo", meta = list(feature_type = "gene"))
    session$flushReact()
    html <- function() as.character(output$preview_literal_ui$html)

    # Paste + exact search-by-id -> the paste toggle.
    session$setInputs(source = "paste", paste_searchby = "__rownames__", paste_mode = "exact")
    session$flushReact()
    expect_match(html(), "paste_literal")
    # A contains search has no literal escape hatch -> no toggle.
    session$setInputs(paste_mode = "contains"); session$flushReact()
    expect_true(is.null(output$preview_literal_ui) || !nzchar(html()))
    # gsfile always offers it.
    session$setInputs(source = "gsfile"); session$flushReact()
    expect_match(html(), "gsfile_literal")
  })
})

test_that("Export enables only with sets, and exports the selected subset", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_geneset_server, args = list(state = state), {
    dds <- ensure_logcounts(make_mock_dds(n_genes = 20, n_per_group = 3, n_spike = 0, seed = 17))
    state_load(state, dds, source = "demo", meta = list(feature_type = "gene"))
    session$flushReact()
    expect_match(as.character(output$export_ui$html), "enable export")

    rn <- rownames(state$working)
    stage_paste(session, rn[1]);  session$setInputs(save_mode = "new", new_name = "A")
    session$setInputs(create = 1); session$flushReact()
    stage_paste(session, rn[2]);  session$setInputs(new_name = "B")
    session$setInputs(create = 1); session$flushReact()
    expect_setequal(names(state$gene_sets), c("A", "B"))

    # The render default selects all (the browser reflects it into input); mirror
    # that here -> both round-trip through the chosen serializer.
    session$setInputs(export_fmt = "json", export_which = c("A", "B")); session$flushReact()
    expect_setequal(names(export_sets()), c("A", "B"))
    expect_setequal(names(gene_sets_from_json(export_txt())), c("A", "B"))
    # Selecting a subset limits the export (and thus the download).
    session$setInputs(export_which = "A"); session$flushReact()
    expect_equal(names(export_sets()), "A")
    expect_equal(names(gene_sets_from_json(export_txt())), "A")
    # An empty selection exports nothing -- an emptied selectize reports NULL, so
    # the absent-selection path must resolve to empty, NOT to "all".
    session$setInputs(export_which = character(0)); session$flushReact()
    expect_length(export_sets(), 0L)
    session$setInputs(export_which = NULL); session$flushReact()
    expect_length(export_sets(), 0L)
  })
})

# ---- Compare tab (P6e) -----------------------------------------------------

test_that("mod_geneset UI mounts the Compare tab (Stats + Overlap)", {
  ui <- as.character(mod_geneset_ui("gs"))
  expect_match(ui, "gs-cmp_sets_ui")       # shared sets-to-visualize control
  expect_match(ui, "gs-cmp_within")        # within-dataset toggle
  expect_match(ui, "gs-stats_container")   # dual_plot Stats bar
  expect_match(ui, "gs-overlap_plot")      # Overlap plot
  expect_match(ui, "gs-overlap_type_ui")
})

test_that(".gs_first_suitable_type picks Euler<=3 / Venn=4 / UpSet>4 (Venn<=4 sans eulerr)", {
  expect_equal(.gs_first_suitable_type(1, TRUE), "euler")
  expect_equal(.gs_first_suitable_type(3, TRUE), "euler")
  expect_equal(.gs_first_suitable_type(4, TRUE), "venn")
  expect_equal(.gs_first_suitable_type(5, TRUE), "upset")
  # Without eulerr, Euler is unavailable: <=4 -> Venn.
  expect_equal(.gs_first_suitable_type(2, FALSE), "venn")
  expect_equal(.gs_first_suitable_type(5, FALSE), "upset")
})

test_that(".gs_type_valid enforces per-type caps (no silent substitution)", {
  expect_true(.gs_type_valid("euler", 3, TRUE))
  expect_false(.gs_type_valid("euler", 4, TRUE))         # over cap -> invalid (message, no plot)
  expect_false(.gs_type_valid("euler", 2, FALSE))        # needs eulerr
  expect_true(.gs_type_valid("venn", 4, TRUE))
  expect_false(.gs_type_valid("venn", 5, TRUE))
  expect_true(.gs_type_valid("upset", 99, TRUE))
  expect_false(.gs_type_valid("upset", 0, TRUE))
})

test_that("Compare: shared selection + Stats frame follow the within-dataset toggle", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_geneset_server, args = list(state = state), {
    dds <- ensure_logcounts(make_mock_dds(n_genes = 30, n_per_group = 3, n_spike = 0, seed = 21))
    state_load(state, dds, source = "demo", meta = list(feature_type = "gene"))
    rn <- rownames(state$working)
    state$gene_sets <- list(
      A = new_gene_set(rn[1:10]),
      B = new_gene_set(c(rn[6:12], "ghostA", "ghostB")))   # 7 present, 2 absent
    session$setInputs(cmp_sets = c("A", "B"), cmp_within = FALSE,
                      stats_auto = TRUE, stats_vertical = FALSE)
    session$flushReact()
    expect_setequal(cmp_selected(), c("A", "B"))

    # within = FALSE: present + absent rows; B carries 2 absent.
    fr <- stats_frame()
    expect_equal(levels(fr$status), c("present", "absent"))
    expect_equal(fr$n[fr$set == "B" & fr$status == "absent"], 2L)
    expect_false(is.null(stats_shown$value()))            # auto-render populated it

    # within = TRUE: present-only, one row per set.
    session$setInputs(cmp_within = TRUE); session$flushReact()
    fw <- stats_frame()
    expect_true(all(fw$status == "present"))
    expect_equal(fw$n[fw$set == "B"], 7L)
  })
})

test_that("Compare: Overlap drops empty sets, defaults type by count, flags invalid picks", {
  skip_if_not_installed("DESeq2")
  skip_if_not_installed("ComplexHeatmap")
  state <- new_app_state()
  shiny::testServer(mod_geneset_server, args = list(state = state), {
    dds <- ensure_logcounts(make_mock_dds(n_genes = 40, n_per_group = 3, n_spike = 0, seed = 22))
    state_load(state, dds, source = "demo", meta = list(feature_type = "gene"))
    rn <- rownames(state$working)
    state$gene_sets <- list(
      A = new_gene_set(rn[1:12]), B = new_gene_set(rn[8:20]),
      C = new_gene_set(rn[15:25]), D = new_gene_set(rn[1:6]),
      E = new_gene_set(rn[3:9]),
      Ghost = new_gene_set(c("nope1", "nope2")))          # all-absent -> empty within
    # 3 sets, type not picked -> default follows the count (first suitable).
    session$setInputs(cmp_sets = c("A", "B", "C"), cmp_within = FALSE, overlap_auto = TRUE)
    session$flushReact()
    od <- overlap_data()
    expect_equal(length(od$lst), 3L)
    expect_equal(od$type, .gs_first_suitable_type(3, have_eulerr()))
    expect_true(.gs_type_valid(od$type, 3L, have_eulerr()))

    # A within-dataset view drops the all-absent Ghost set from the overlap.
    session$setInputs(cmp_sets = c("A", "Ghost"), cmp_within = TRUE); session$flushReact()
    expect_equal(length(overlap_lst()), 1L)               # Ghost contributes nothing

    # User picks Venn, then selects 5 sets -> invalid, so the plot must show the
    # message (not a substituted UpSet): the validity check is FALSE.
    session$setInputs(cmp_sets = c("A", "B", "C", "D", "E"), cmp_within = FALSE,
                      overlap_type = "venn"); session$flushReact()
    expect_equal(overlap_data()$type, "venn")             # respected, not switched
    expect_false(.gs_type_valid(overlap_data()$type, length(overlap_lst()), have_eulerr()))

    # UpSet handles any n and builds.
    session$setInputs(overlap_type = "upset"); session$flushReact()
    expect_equal(overlap_data()$type, "upset")
    p <- .gs_overlap_plot(overlap_data()$lst, "upset")
    expect_true(inherits(p, c("Heatmap", "HeatmapList")))
  })
})
