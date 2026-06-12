test_that("state_load / mutate / derive / undo / reset behave correctly", {
  skip_if_not_installed("DESeq2")
  shiny::reactiveConsole(TRUE)
  on.exit(shiny::reactiveConsole(FALSE), add = TRUE)

  st <- new_app_state()
  expect_false(isTRUE(state_meta(st)$loaded))

  dds <- make_mock_dds(n_genes = 30, n_per_group = 2, n_spike = 2, seed = 3)
  state_load(st, dds, source = "demo", meta = list(data_type = "bulk"))
  expect_equal(st$data_version, 1L)
  expect_equal(length(st$history), 1L)
  expect_identical(st$working, st$original)
  expect_equal(state_meta(st)$n_features, nrow(dds))
  expect_equal(state_meta(st)$n_edits, 0L)

  # state_derive caches under the current data_version + params.
  calls <- 0L
  f <- function() { calls <<- calls + 1L; 42 }
  expect_equal(state_derive(st, "x", list(a = 1), f), 42)
  expect_equal(state_derive(st, "x", list(a = 1), f), 42)
  expect_equal(calls, 1L)                                  # second call cached
  state_derive(st, "x", list(a = 2), f)
  expect_equal(calls, 2L)                                  # params changed -> recompute

  # state_mutate bumps the version, invalidates the cache, logs, supports undo.
  v0 <- st$data_version
  state_mutate(st, function(d) d[seq_len(10), ], action = list(action = "subset", n = 10))
  expect_equal(st$data_version, v0 + 1L)
  expect_equal(nrow(st$working), 10L)
  expect_equal(state_meta(st)$n_edits, 1L)
  state_derive(st, "x", list(a = 1), f)
  expect_equal(calls, 3L)                                  # cache invalidated by the edit

  state_undo(st)
  expect_equal(nrow(st$working), nrow(dds))

  state_reset(st)
  expect_identical(st$working, st$original)
  expect_equal(length(st$undo_stack), 0L)
})

test_that("state_mutate errors when nothing is loaded", {
  shiny::reactiveConsole(TRUE)
  on.exit(shiny::reactiveConsole(FALSE), add = TRUE)
  st <- new_app_state()
  expect_error(state_mutate(st, identity), "No dataset loaded")
})
