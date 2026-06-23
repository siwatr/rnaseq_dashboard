test_that("plot_subset_ui renders a 'Plot Showing' accordion with the controls", {
  ns <- shiny::NS("qc")
  html <- paste(as.character(ddsdashboard:::plot_subset_ui(ns, "gen")), collapse = " ")
  expect_match(html, "qc-gen_show_by")
  expect_match(html, "qc-gen_show_values")
  expect_match(html, "Plot Showing", fixed = TRUE)
})
# The synced-state + showing_samples() behaviour is exercised end-to-end via the
# mod_qc testServer tests (Showing subset / cross-tab sync).
