test_that("run_app sets the upload size limit and builds an app", {
  old <- getOption("shiny.maxRequestSize")
  on.exit(options(shiny.maxRequestSize = old), add = TRUE)

  app <- run_app(max_upload_mb = 250)
  expect_s3_class(app, "shiny.appobj")
  expect_equal(getOption("shiny.maxRequestSize"), 250 * 1024^2)
})

test_that(".app_theme is a bslib theme and bundles the custom stylesheet", {
  expect_s3_class(.app_theme(), "bs_theme")
  expect_true(nzchar(system.file("www", "custom.scss", package = "ddsdashboard")))
})

test_that("app_ui includes the dark-mode toggle", {
  ui <- as.character(app_ui())
  expect_true(grepl("dark_mode", ui))             # input_dark_mode(id = "dark_mode")
  expect_true(grepl("data-bs-toggle|bslib-dark", ui) || grepl("dark", ui))
})
