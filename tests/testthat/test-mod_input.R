test_that("input module loads demo data into shared state", {
  skip_if_not_installed("DESeq2")

  st <- new_app_state()
  shiny::testServer(mod_input_server, args = list(state = st), {
    session$setInputs(source = "demo")
    session$setInputs(load = 1)

    expect_true(!is.null(st$working))
    expect_s4_class(st$working, "DESeqDataSet")
    m <- state_meta(st)
    expect_true(m$loaded)
    expect_equal(m$data_type, "bulk")
    expect_equal(m$feature_type, "gene")
    # on-load conventions applied
    expect_true("logcounts" %in% SummarizedExperiment::assayNames(st$working))
    expect_true("feature_class" %in% colnames(SummarizedExperiment::rowData(st$working)))
  })
})
