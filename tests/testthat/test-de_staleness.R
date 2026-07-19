# Content-addressed DESeq-fit / VST caching: an assay-add must not invalidate the
# fit or VST (counts / samples / features / size factors are unchanged), while a
# real structural or size-factor change must.

test_that("dds_content_fingerprint tracks samples, features, and (effective) size factors", {
  skip_if_not_installed("DESeq2")
  dds <- ensure_logcounts(make_mock_dds(n_genes = 40, n_per_group = 3, n_spike = 4, seed = 1))
  fp <- dds_content_fingerprint(dds)
  expect_equal(fp$rn, rownames(dds))
  expect_equal(fp$cn, colnames(dds))
  expect_null(dds_content_fingerprint(NULL))
  # NULL size factors == the endogenous estimate they would become: materializing
  # that same default is NOT a change (this is what fixes the assay-add bug).
  dds_sf <- estimate_size_factors_endogenous(dds)
  expect_identical(fp$sf, dds_content_fingerprint(dds_sf)$sf)
  # ...but a genuinely different size-factor vector IS a change.
  dds_other <- dds_sf
  DESeq2::sizeFactors(dds_other) <- DESeq2::sizeFactors(dds_sf) * 2
  expect_false(identical(dds_content_fingerprint(dds_sf), dds_content_fingerprint(dds_other)))
  # dropping a sample / feature changes it
  expect_false(identical(fp, dds_content_fingerprint(dds[, -1])))
  expect_false(identical(fp, dds_content_fingerprint(dds[-1, ])))
})

test_that(".clear_derived keeps content-addressed keys only when asked", {
  shiny::reactiveConsole(TRUE); on.exit(shiny::reactiveConsole(FALSE), add = TRUE)
  state <- new_app_state()
  assign("de_fit", list(value = 1), envir = state$derived)
  assign("vst", list(value = 2), envir = state$derived)
  assign("pca", list(value = 3), envir = state$derived)
  .clear_derived(state, keep = .content_derived_keys)
  expect_true(exists("de_fit", envir = state$derived, inherits = FALSE))
  expect_true(exists("vst", envir = state$derived, inherits = FALSE))
  expect_false(exists("pca", envir = state$derived, inherits = FALSE))  # version-addressed -> wiped
  .clear_derived(state)                                                  # full wipe
  expect_false(exists("de_fit", envir = state$derived, inherits = FALSE))
})

test_that("state_derive validates on a supplied version token", {
  shiny::reactiveConsole(TRUE); on.exit(shiny::reactiveConsole(FALSE), add = TRUE)
  state <- new_app_state()
  calls <- 0L
  gen <- function(v) state_derive(state, "vst", params = list(),
                                  expr = function() { calls <<- calls + 1L; v },
                                  version = list(tag = v))
  expect_equal(gen("a"), "a"); expect_equal(calls, 1L)
  expect_equal(gen("a"), "a"); expect_equal(calls, 1L)   # same version -> cached
  expect_equal(gen("b"), "b"); expect_equal(calls, 2L)   # new version -> recompute
})

test_that("adding an assay preserves the DESeq fit + VST; a structural edit stales them", {
  skip_if_not_installed("DESeq2")
  shiny::reactiveConsole(TRUE); on.exit(shiny::reactiveConsole(FALSE), add = TRUE)
  state <- new_app_state()
  # Load WITHOUT pre-estimating size factors -- exactly the app's state when a user
  # runs DE before ever adding an assay (de_run estimates internally; working keeps
  # NULL size factors). This is the scenario that regressed.
  dds <- ensure_logcounts(make_mock_dds(n_genes = 60, n_per_group = 3, n_spike = 6, seed = 2))
  state_load(state, dds, source = "demo", meta = list(feature_type = "gene"))
  expect_null(DESeq2::sizeFactors(state$working))

  # Simulate a completed fit + extracted results, and a content-addressed VST.
  assign("de_fit", list(value = "FIT", stamp = .de_stamp(state)), envir = state$derived)
  state$de <- list(stamp = .de_stamp(state), results = list(c1 = data.frame(x = 1)))
  vst_calls <- 0L
  get_vst <- function() state_derive(state, "vst", params = list(),
    version = dds_content_fingerprint(state$working),
    expr = function() { vst_calls <<- vst_calls + 1L; "VST" })
  get_vst()
  expect_equal(vst_calls, 1L)
  expect_equal(de_fit_status(state), "current")
  expect_equal(de_status(state), "current")

  # Add a normalized assay (re-estimates size factors to the SAME values).
  state_mutate(state, function(d) estimate_size_factors_endogenous(add_normalized_assays(d, "CPM")),
               action = list(action = "add_assays", assays = "CPM"))
  expect_true(exists("de_fit", envir = state$derived, inherits = FALSE))  # not wiped
  expect_equal(de_fit_status(state), "current")                          # fit still valid
  expect_equal(de_status(state), "current")                              # results still valid
  get_vst(); expect_equal(vst_calls, 1L)                                 # VST reused, not recomputed

  # A real structural edit (drop a sample) invalidates fit + VST.
  state_mutate(state, function(d) d[, -1], action = list(action = "drop_sample"))
  expect_equal(de_fit_status(state), "stale")
  get_vst(); expect_equal(vst_calls, 2L)                                 # VST recomputed
})
