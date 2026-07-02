# Shared gene-search sub-module (R/mod_gene_search.R). Runs in the host module's
# namespace, so exercise it through a tiny wrapper moduleServer that returns the
# reactive it produces.

test_that("gene_search_ui builds the suffixed controls (single vs multi)", {
  ns <- shiny::NS("host")
  single <- paste(as.character(gene_search_ui(ns, "g", multiple = FALSE)), collapse = " ")
  expect_match(single, "host-g_searchby_ui")
  expect_match(single, "host-g_ci")
  expect_match(single, "host-g_q")
  expect_match(single, "host-g_hint")
  expect_no_match(single, "host-g_dup")                 # dup switch off by default
  expect_no_match(single, "host-g_mode")                # no selector for one mode
  expect_no_match(single, "host-g_split_hint")          # no split hint in single mode

  dup <- paste(as.character(gene_search_ui(ns, "g", dup_toggle = TRUE)), collapse = " ")
  expect_match(dup, "host-g_dup")

  # Multi mode: textarea + a split hint; the mode selector appears only when more
  # than one search mode is offered.
  multi <- paste(as.character(gene_search_ui(ns, "g", multiple = TRUE,
                    search_modes = c("exact", "contains", "regex"))), collapse = " ")
  expect_match(multi, "textarea")
  expect_match(multi, "host-g_mode")
  expect_match(multi, "host-g_split_hint")
  expect_match(multi, "Contains")                        # mode label rendered
  # A single-mode multi box shows no mode selector.
  multi1 <- paste(as.character(gene_search_ui(ns, "g", multiple = TRUE)), collapse = " ")
  expect_no_match(multi1, "host-g_mode")
})

test_that(".gene_search_split: exact/contains split on comma+newline, regex on newline", {
  expect_equal(ddsdashboard:::.gene_search_split("A, b\nC", "exact"), c("A", "b", "C"))
  expect_equal(ddsdashboard:::.gene_search_split("A, b\nC", "contains"), c("A", "b", "C"))
  # A comma is valid regex, so regex keeps 'A, b' as one pattern (splits newlines only).
  expect_equal(ddsdashboard:::.gene_search_split("A, b\nC", "regex"), c("A, b", "C"))
  expect_equal(ddsdashboard:::.gene_search_split("x, x, y", "exact"), c("x", "y"))  # deduped
})

# A wrapper module returning the resolver reactive, for testServer.
.gs_wrapper <- function(multiple, search_modes = "exact") {
  function(id, state) shiny::moduleServer(id, function(input, output, session)
    gene_search_server(input, output, session, state, "g", multiple = multiple,
                       search_modes = search_modes))
}

test_that("single mode resolves one id, reports a miss, honours 'Search by'", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(.gs_wrapper(FALSE), args = list(state = state), {
    dds <- ensure_logcounts(make_mock_dds(n_genes = 30, n_per_group = 3, n_spike = 0, seed = 1))
    state_load(state, dds, source = "demo", meta = list(feature_type = "gene"))
    session$flushReact()
    rid <- rownames(state$working)[1]
    nm  <- as.character(SummarizedExperiment::rowData(state$working)$gene_name[1])

    # Search by feature id (rownames).
    session$setInputs(g_searchby = "__rownames__", g_ci = FALSE, g_q = rid)
    session$elapse(300); session$flushReact()
    r <- session$returned()
    expect_equal(r$ids, rid)
    expect_length(r$records, 1L)
    expect_equal(r$records[[1]]$id, rid)
    expect_length(r$unmatched, 0L)

    # Search by the gene_name column resolves to the same id.
    session$setInputs(g_searchby = "gene_name", g_q = nm)
    session$elapse(300); session$flushReact()
    expect_equal(session$returned()$ids, rid)

    # A miss -> no ids, the query reported as unmatched, and a not-found hint.
    session$setInputs(g_q = "NoSuchGene")
    session$elapse(300); session$flushReact()
    r2 <- session$returned()
    expect_length(r2$ids, 0L)
    expect_equal(r2$unmatched, "NoSuchGene")
    expect_match(as.character(output$g_hint$html), "not found")
  })
})

test_that("multi mode splits on comma/newline, dedups, reports unmatched", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(.gs_wrapper(TRUE), args = list(state = state), {
    dds <- ensure_logcounts(make_mock_dds(n_genes = 30, n_per_group = 3, n_spike = 0, seed = 2))
    state_load(state, dds, source = "demo", meta = list(feature_type = "gene"))
    session$flushReact()
    r1 <- rownames(state$working)[1]; r2 <- rownames(state$working)[2]

    session$setInputs(g_searchby = "__rownames__", g_ci = FALSE,
                      g_q = paste0(r1, ", ", r2, "\n", r1, ", ZZZ"))
    session$elapse(300); session$flushReact()
    out <- session$returned()
    expect_setequal(out$ids, c(r1, r2))          # deduped (r1 appears twice)
    expect_equal(out$unmatched, "ZZZ")
    expect_match(as.character(output$g_hint$html), "Not found: ZZZ")
  })
})

test_that("case-insensitive toggle controls matching", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(.gs_wrapper(FALSE), args = list(state = state), {
    dds <- ensure_logcounts(make_mock_dds(n_genes = 30, n_per_group = 3, n_spike = 0, seed = 3))
    state_load(state, dds, source = "demo", meta = list(feature_type = "gene"))
    session$flushReact()
    nm <- as.character(SummarizedExperiment::rowData(state$working)$gene_name[5])

    session$setInputs(g_searchby = "gene_name", g_ci = FALSE, g_q = tolower(nm))
    session$elapse(300); session$flushReact()
    expect_length(session$returned()$ids, 0L)     # case-sensitive miss

    session$setInputs(g_ci = TRUE)
    session$elapse(300); session$flushReact()
    r <- session$returned()
    expect_length(r$ids, 1L)
    expect_equal(r$records[[1]]$match, nm)         # stored (correct-case) value
  })
})

test_that("contains + regex modes match many; regex does not comma-split", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  wrap <- .gs_wrapper(TRUE, c("exact", "contains", "regex"))
  shiny::testServer(wrap, args = list(state = state), {
    # mock gene_name: mt-Gene1..mt-Gene5 then Gene1, Gene2, ...
    dds <- ensure_logcounts(make_mock_dds(n_genes = 12, n_per_group = 3, n_spike = 0, seed = 4))
    state_load(state, dds, source = "demo", meta = list(feature_type = "gene"))
    session$flushReact()
    mt_ids <- rownames(state$working)[grepl("^mt-",
                 as.character(SummarizedExperiment::rowData(state$working)$gene_name))]
    expect_true(length(mt_ids) >= 3L)

    # Contains "mt-" on gene_name -> all the mt- features.
    session$setInputs(g_searchby = "gene_name", g_ci = FALSE, g_mode = "contains", g_q = "mt-")
    session$elapse(300); session$flushReact()
    expect_setequal(session$returned()$ids, mt_ids)

    # Regex "^mt-" -> the same set.
    session$setInputs(g_mode = "regex", g_q = "^mt-")
    session$elapse(300); session$flushReact()
    expect_setequal(session$returned()$ids, mt_ids)

    # Regex treats "mt-Gene1, mt-Gene2" as ONE pattern (comma not split) -> no hit.
    session$setInputs(g_q = "mt-Gene1, mt-Gene2")
    session$elapse(300); session$flushReact()
    r_rx <- session$returned()
    expect_length(r_rx$ids, 0L)
    expect_equal(r_rx$n_query, 1L)
    # Exact mode DOES comma-split the same text -> both genes hit.
    session$setInputs(g_mode = "exact", g_q = "mt-Gene1, mt-Gene2")
    session$elapse(300); session$flushReact()
    expect_length(session$returned()$ids, 2L)

    # An invalid regex is reported, not thrown.
    session$setInputs(g_mode = "regex", g_q = "mt-Gene[")
    session$elapse(300); session$flushReact()
    r_bad <- session$returned()
    expect_equal(r_bad$invalid, "mt-Gene[")
    expect_match(as.character(output$g_hint$html), "Invalid regex")
  })
})

test_that("multi-mode hint tiers: coverage, per-term suggestions, wrong-column", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  wrap <- .gs_wrapper(TRUE, c("exact", "contains", "regex"))
  shiny::testServer(wrap, args = list(state = state), {
    dds <- ensure_logcounts(make_mock_dds(n_genes = 12, n_per_group = 3, n_spike = 0, seed = 5))
    state_load(state, dds, source = "demo", meta = list(feature_type = "gene"))
    session$flushReact()

    # Exact miss with a small unmatched count -> a "did you mean" suggestion.
    # "mt-Gene" is not an exact name but is a substring of mt-Gene1..5.
    session$setInputs(g_searchby = "gene_name", g_mode = "exact", g_ci = FALSE, g_q = "mt-Gene")
    session$elapse(300); session$flushReact()
    h1 <- as.character(output$g_hint$html)
    expect_match(h1, "of 1 not found")
    expect_match(h1, "did you mean")

    # Many unmatched (>=20 and >=50%) -> the wrong-column hint.
    fakes <- paste0("FAKE", seq_len(25))
    session$setInputs(g_q = paste(fakes, collapse = "\n"))
    session$elapse(300); session$flushReact()
    h2 <- as.character(output$g_hint$html)
    expect_match(h2, "25 of 25 not found")
    expect_match(h2, "different ID column")
  })
})
