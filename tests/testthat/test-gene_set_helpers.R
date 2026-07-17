# Pure gene-set helpers (R/gene_set_helpers.R).

test_that("new_gene_set dedups, drops NA/blank, stamps kind/source", {
  s <- new_gene_set(c("a", "b", "a", NA, "", "c"), source = "paste")
  expect_equal(s$ids, c("a", "b", "c"))
  expect_equal(s$kind, "simple")
  expect_null(s$annotation)
  expect_equal(s$source, "paste")
  expect_error(new_gene_set("a", kind = "bogus"))
})

test_that("gene_set_commit New auto-suffixes on a name clash", {
  sets <- list()
  r1 <- gene_set_commit(sets, "myset", c("a", "b"), "new", source = "paste")
  expect_named(r1$sets, "myset")
  expect_equal(r1$name, "myset")
  expect_equal(r1$sets$myset$ids, c("a", "b"))

  r2 <- gene_set_commit(r1$sets, "myset", c("c"), "new")
  expect_equal(r2$name, "myset_2")               # collision -> suffix
  expect_setequal(names(r2$sets), c("myset", "myset_2"))
  expect_equal(r2$sets$myset$ids, c("a", "b"))   # original untouched

  r3 <- gene_set_commit(r2$sets, "myset", c("d"), "new")
  expect_equal(r3$name, "myset_3")               # next free suffix
})

test_that("gene_set_commit Append unions + dedups; Append to a new name creates", {
  sets <- gene_set_commit(list(), "s", c("a", "b"), "new")$sets
  r <- gene_set_commit(sets, "s", c("b", "c"), "append")
  expect_equal(r$name, "s")
  expect_setequal(r$sets$s$ids, c("a", "b", "c"))     # union, deduped

  # Append to a non-existent name just creates it.
  r2 <- gene_set_commit(sets, "fresh", c("x"), "append")
  expect_equal(r2$name, "fresh")
  expect_equal(r2$sets$fresh$ids, "x")
})

test_that("gene_set_commit rejects an empty name", {
  expect_error(gene_set_commit(list(), "  ", c("a"), "new"), "non-empty")
})

test_that("gene_set_present / gene_set_absent derive the live view (order-preserving)", {
  set <- new_gene_set(c("g3", "g1", "gX"))
  feats <- c("g1", "g2", "g3")
  expect_equal(gene_set_present(set, feats), c("g3", "g1"))   # set order preserved
  expect_equal(gene_set_absent(set, feats), "gX")
  # Accepts a bare id vector too.
  expect_equal(gene_set_present(c("g1", "gX"), feats), "g1")
  expect_equal(gene_set_absent(c("g1", "gX"), feats), "gX")
})
