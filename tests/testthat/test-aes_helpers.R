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

test_that("spike-in metrics are catalog attributes only when the dataset has spikes", {
  with_spike <- mk_state(seed = 5)                       # n_spike = 4
  keys <- vapply(aes_catalog(with_spike), `[[`, "", "key")
  expect_true("__spike__slope" %in% keys)
  byk <- function(cat, k) Filter(function(d) d$key == k, cat)[[1]]
  cat_s <- aes_catalog(with_spike)
  d <- byk(cat_s, "__spike__slope")
  expect_equal(d$group, "Spike-in"); expect_equal(d$kind, "continuous")
  # % spike-in (pct_spike) is grouped under "Spike-in", not "This session".
  expect_equal(byk(cat_s, "__qc__pct_spike")$group, "Spike-in")
  expect_equal(byk(cat_s, "__qc__library_size")$group, "This session")
  # No spikes -> no spike attributes offered.
  st0 <- new_app_state()
  shiny::reactiveConsole(TRUE); withr::defer(shiny::reactiveConsole(FALSE))
  state_load(st0, ensure_logcounts(make_mock_dds(n_genes = 50, n_per_group = 3, n_spike = 0, seed = 6)),
             source = "demo")
  keys0 <- vapply(aes_catalog(st0), `[[`, "", "key")
  expect_false(any(grepl("^__spike__", keys0)))
  expect_false("__qc__pct_spike" %in% keys0)             # % spike-in is spike-gated too
})

test_that("aes_resolve resolves a spike metric per sample (cached via spike_dr)", {
  state <- mk_state(seed = 5); s <- colnames(state$working)
  r <- aes_resolve(state, "__spike__n_spike_detected", s)
  expect_equal(r$kind, "continuous"); expect_length(r$values, length(s))
  expect_true(all(is.finite(r$values)))                  # detected-spike counts
  expect_true(exists("spike_dr", envir = state$derived, inherits = FALSE))
  expect_equal(r$label, "Detected spike features")
})

test_that("aes_choices builds a 'Spike-in' optgroup when spikes exist", {
  state <- mk_state(seed = 5)
  ch <- aes_choices(aes_catalog(state, gene = TRUE), none = TRUE)
  expect_true("Spike-in" %in% names(ch))
  expect_true("__spike__r_squared" %in% unlist(ch))
  # Shape (discrete only) excludes the continuous spike metrics.
  disc <- aes_choices(aes_catalog(state), kinds = "discrete", none = TRUE, state = state)
  expect_false("Spike-in" %in% names(disc))
})

test_that("aes_other_palette_items lists removal/pool/QC-metric + spike customizable attrs", {
  it <- aes_other_palette_items()
  expect_true(all(c("removal_status", "__pool__", "__qc__library_size",
                    "__spike__slope") %in% names(it)))
  expect_equal(it[["__pool__"]]$kind, "discrete")
  expect_equal(it[["__pool__"]]$levels, c("Kept", "In removal pool"))
  expect_equal(it[["__qc__library_size"]]$kind, "continuous")
  expect_equal(it[["__qc__pct_mito"]]$label, "% mitochondrial")
})

test_that("aes_resolve pool honours an 'other' palette config (else the default)", {
  state <- mk_state(); s <- colnames(state$working)
  d <- aes_resolve(state, "__pool__", s)
  expect_equal(unname(d$colors["In removal pool"]), "#D62728")   # built-in default
  # A configured custom map overrides the pool colours.
  state$palette <- list(other = list(`__pool__` = list(name = "Custom palette",
    colors = c(Kept = "#101010", "In removal pool" = "#FEFEFE"))))
  d2 <- aes_resolve(state, "__pool__", s)
  expect_equal(unname(d2$colors["Kept"]), "#101010")
  expect_equal(unname(d2$colors["In removal pool"]), "#FEFEFE")
})

test_that("a configured QC-metric ramp flows through to ggplot + heatmap scales", {
  state <- mk_state(); s <- colnames(state$working)
  # No config -> thematic default (NULL ggplot scale); heatmap default ramp.
  r0 <- aes_resolve(state, "__qc__library_size", s)
  expect_null(r0$ramp_config)
  expect_null(aes_ggplot_scale(r0))
  # Configure a continuous ramp in the Palette "Other" slot for this metric.
  state$palette <- list(other = list(`__qc__library_size` = list(
    name = "viridis: viridis", min = "", max = "", custom = NULL, reverse = FALSE)))
  r <- aes_resolve(state, "__qc__library_size", s)
  expect_equal(r$ramp_config$name, "viridis: viridis")
  expect_s3_class(aes_ggplot_scale(r), "Scale")              # gradient scale built
  skip_if_not_installed("circlize")
  expect_true(is.function(aes_heatmap_col(r)))               # colorRamp2 from the config
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
