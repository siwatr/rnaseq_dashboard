test_that(".format_bytes renders binary units", {
  expect_equal(.format_bytes(0), "0.0 B")
  expect_equal(.format_bytes(1023), "1023.0 B")
  expect_equal(.format_bytes(1024), "1.0 KB")
  expect_equal(.format_bytes(1024^2 * 2), "2.0 MB")
  expect_equal(.format_bytes(1024^3), "1.0 GB")
  expect_true(is.na(.format_bytes(NA)))
})

test_that("statusbar renders dataset and memory badges", {
  state <- new_app_state()
  shiny::testServer(mod_statusbar_server, args = list(state = state), {
    expect_true(any(grepl("no dataset loaded", as.character(output$status))))
    if (requireNamespace("ps", quietly = TRUE)) {
      expect_true(any(grepl("mem:", as.character(output$mem))))
    }
  })
})
