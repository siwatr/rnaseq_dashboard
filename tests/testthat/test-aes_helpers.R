# Shared color/annotation attribute catalog + resolver (R/aes_helpers.R).

mk_state <- function(seed = 1, gene_extra = TRUE) {
  skip_if_not_installed("DESeq2")
  # Allow reactiveValues reads/writes at the console (these helpers read state
  # outside a Shiny reactive context in unit tests).
  shiny::reactiveConsole(TRUE)
  withr::defer(shiny::reactiveConsole(FALSE), envir = parent.frame())
  state <- new_app_state()
  dds <- ensure_logcounts(make_mock_dds(n_genes = 80, n_per_group = 3, n_spike = 4, seed = seed))
  state_load(state, dds, source = "demo", meta = list(feature_type = "gene"))
  s <- colnames(state$working)
  state$samp_flags <- data.frame(sample = s, flagged = c(TRUE, rep(FALSE, length(s) - 1L)),
                                 stringsAsFactors = FALSE)
  state$samp_pool <- s[2]
  state
}

test_that("aes_catalog lists colData + QC metrics + removal/pool (+ gene when asked)", {
  state <- mk_state()
  cat0 <- aes_catalog(state, gene = FALSE)
  keys <- vapply(cat0, `[[`, "", "key")
  expect_true(all(c("condition", "__qc__library_size", "__removal__", "__pool__") %in% keys))
  expect_false("__gene__" %in% keys)
  expect_true("__gene__" %in% vapply(aes_catalog(state, gene = TRUE), `[[`, "", "key"))
  # Kinds + groups are tagged correctly.
  byk <- function(k) Filter(function(d) d$key == k, cat0)[[1]]
  expect_equal(byk("condition")$kind, "discrete")
  expect_equal(byk("__qc__library_size")$kind, "continuous")
  expect_equal(byk("__removal__")$group, "This session")
  expect_equal(byk("condition")$group, "Data metadata")
})

test_that("aes_choices builds grouped optgroups and respects kind + level cap", {
  state <- mk_state()
  ch <- aes_choices(aes_catalog(state, gene = TRUE), none = TRUE)
  expect_equal(names(ch)[1], "General")
  expect_true("This session" %in% names(ch) && "Data metadata" %in% names(ch))
  expect_true("__removal__" %in% unlist(ch))
  # Discrete-only (shape) drops continuous attributes.
  disc <- aes_choices(aes_catalog(state), kinds = "discrete", none = TRUE, state = state)
  expect_false("__qc__library_size" %in% unlist(disc))
  expect_true("__removal__" %in% unlist(disc))
})

test_that("aes_resolve aligns colData/QC/removal/pool to the given samples", {
  state <- mk_state()
  s <- colnames(state$working)
  # colData discrete -> factor + (no config) NULL colours.
  r <- aes_resolve(state, "condition", s)
  expect_equal(r$kind, "discrete"); expect_s3_class(r$values, "factor"); expect_null(r$colors)
  # QC metric -> numeric aligned by sample id.
  q <- aes_resolve(state, "__qc__library_size", s)
  expect_equal(q$kind, "continuous"); expect_length(q$values, length(s))
  # removal: sample 1 flagged -> "Suggested removal"; others "QC pass".
  rm <- aes_resolve(state, "__removal__", s)
  expect_equal(as.character(rm$values[1]), "Suggested removal")
  expect_true("QC pass" %in% levels(rm$values))
  expect_true(all(levels(rm$values) %in% names(rm$colors)))
  # pool: sample 2 staged.
  pl <- aes_resolve(state, "__pool__", s)
  expect_equal(as.character(pl$values[2]), "In removal pool")
  expect_named(pl$colors, c("Kept", "In removal pool"))
})

test_that("aes_resolve removal is NULL before flags exist; reason-aware adds a third level", {
  state <- mk_state(); s <- colnames(state$working)
  state$samp_flags <- NULL
  expect_null(aes_resolve(state, "__removal__", s))
  # Reason-aware path uses the 3-level labels when a metric reason is given.
  state$samp_flags <- data.frame(sample = s, flagged = TRUE, high_mito = c(TRUE, rep(FALSE, length(s) - 1L)),
                                 stringsAsFactors = FALSE)
  rr <- aes_resolve(state, "__removal__", s, ctx = list(reason = "pct_mito"))
  expect_true(any(grepl("this reason", levels(rr$values))))
})

test_that("aes_resolve honours a colData discrete palette config", {
  state <- mk_state()
  lv <- levels(factor(as.character(SummarizedExperiment::colData(state$working)$condition)))
  state$palette <- list(colData = list(condition = list(name = "Okabe-Ito",
    colors = palette_discrete(lv, NULL, "Okabe-Ito"))))
  r <- aes_resolve(state, "condition", colnames(state$working))
  expect_true(is.character(r$colors)); expect_named(r$colors, lv, ignore.order = TRUE)
})

test_that("aes_ggplot_scale: NULL without config, manual/gradient with one", {
  state <- mk_state(); s <- colnames(state$working)
  expect_null(aes_ggplot_scale(aes_resolve(state, "condition", s)))            # no config
  expect_null(aes_ggplot_scale(aes_resolve(state, "__qc__detected", s)))       # continuous, no config
  expect_s3_class(aes_ggplot_scale(aes_resolve(state, "__removal__", s)), "Scale")  # fixed colours
  expect_s3_class(aes_ggplot_scale(aes_resolve(state, "__pool__", s), "fill"), "Scale")
})

test_that("aes_annotation builds a mixed-type frame + colour list", {
  state <- mk_state(); s <- colnames(state$working)
  a <- aes_annotation(state, c("condition", "__removal__", "__pool__", "__qc__library_size"), s)
  expect_true(all(c("condition", "Suggested removal", "Removal pool", "Library size")
                  %in% colnames(a$df)))
  expect_s3_class(a$df[["Removal pool"]], "factor")
  expect_true(is.numeric(a$df[["Library size"]]))
  expect_true("In removal pool" %in% names(a$col[["Removal pool"]]))
  skip_if_not_installed("circlize")
  expect_true(is.function(a$col[["Library size"]]))           # colorRamp2
})
