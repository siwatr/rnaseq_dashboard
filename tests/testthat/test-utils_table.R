test_that("dt_table returns a datatable widget with the standard options", {
  w <- dt_table(head(mtcars))
  expect_s3_class(w, "datatables")
  opts <- w$x$options
  expect_equal(opts$pageLength, 10L)
  expect_equal(opts$lengthMenu, c(10, 25, 50, 100))
  expect_true(isTRUE(opts$scrollX))
  expect_match(opts$dom, "l")   # rows-per-page selector
  expect_match(opts$dom, "f")   # search box
  # filter = "top" adds a header filter row to the widget.
  expect_true(!is.null(w$x$filter) && !identical(w$x$filter, "none"))
})

test_that("dt_table lets callers override and extend options", {
  w <- dt_table(head(mtcars), page_length = 25L, options = list(dom = "tp"))
  expect_equal(w$x$options$pageLength, 25L)
  expect_equal(w$x$options$dom, "tp")        # override wins
  expect_equal(w$x$options$lengthMenu, c(10, 25, 50, 100))  # default retained
})
