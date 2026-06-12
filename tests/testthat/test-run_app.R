test_that("run_app sets the upload size limit and builds an app", {
  old <- getOption("shiny.maxRequestSize")
  on.exit(options(shiny.maxRequestSize = old), add = TRUE)

  app <- run_app(max_upload_mb = 250)
  expect_s3_class(app, "shiny.appobj")
  expect_equal(getOption("shiny.maxRequestSize"), 250 * 1024^2)
})
