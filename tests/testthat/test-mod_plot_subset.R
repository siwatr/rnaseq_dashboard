test_that("plot_subset_ui renders the show-by + keep controls for a suffix", {
  ns <- shiny::NS("qc")
  html <- paste(as.character(ddsdashboard:::plot_subset_ui(ns, "gen")), collapse = " ")
  expect_match(html, "qc-gen_show_by")
  expect_match(html, "qc-gen_show_values")
  expect_match(html, "Showing (display only)", fixed = TRUE)
})
# The synced-state + showing_samples() behaviour is exercised end-to-end via the
# mod_qc testServer tests (Showing subset / cross-tab sync).
