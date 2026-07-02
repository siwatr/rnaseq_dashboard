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

  dup <- paste(as.character(gene_search_ui(ns, "g", dup_toggle = TRUE)), collapse = " ")
  expect_match(dup, "host-g_dup")

  # Multi mode swaps the single-line input for a textarea.
  multi <- paste(as.character(gene_search_ui(ns, "g", multiple = TRUE)), collapse = " ")
  expect_match(multi, "textarea")
})

# A wrapper module returning the resolver reactive, for testServer.
.gs_wrapper <- function(multiple) {
  function(id, state) shiny::moduleServer(id, function(input, output, session)
    gene_search_server(input, output, session, state, "g", multiple = multiple))
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
