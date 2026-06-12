test_that("cpm columns sum to 1e6", {
  m <- matrix(c(1, 1, 2, 0, 3, 7), nrow = 3,
              dimnames = list(c("g1", "g2", "g3"), c("s1", "s2")))
  out <- cpm(m)
  expect_equal(colSums(out), c(s1 = 1e6, s2 = 1e6))
  expect_equal(dim(out), dim(m))
})

test_that("tpm length-normalizes then scales to 1e6", {
  m <- matrix(c(10, 10, 10, 10), nrow = 2,
              dimnames = list(c("g1", "g2"), c("s1", "s2")))
  len <- c(1000, 4000)
  out <- tpm(m, len)
  expect_equal(colSums(out), c(s1 = 1e6, s2 = 1e6))
  # g1 is 4x shorter than g2 with equal counts -> 4x the TPM.
  expect_equal(unname(out["g1", "s1"] / out["g2", "s1"]), 4)
})

test_that("fpkm matches the explicit formula", {
  m <- matrix(c(1000, 3000), nrow = 2,
              dimnames = list(c("g1", "g2"), "s1"))
  len <- c(2000, 1000)
  lib <- sum(m)  # 4000
  expected <- m[, 1] / (len / 1e3) / (lib / 1e6)
  expect_equal(fpkm(m, len)[, 1], expected)
})

test_that("logcounts_from_counts is log2(cpm + 1)", {
  m <- matrix(c(1, 0, 4, 5), nrow = 2,
              dimnames = list(c("g1", "g2"), c("s1", "s2")))
  expect_equal(logcounts_from_counts(m), log2(cpm(m) + 1))
})

test_that("input validation rejects bad matrices and lengths", {
  m <- matrix(c(1, 2, 3, 4), nrow = 2)
  expect_error(cpm(matrix(c(-1, 2, 3, 4), nrow = 2)), "non-negative")
  expect_error(cpm("not a matrix"), "numeric matrix")
  expect_error(tpm(m, c(1000)), "one value per feature")
  expect_error(tpm(m, c(0, 1000)), "positive")
})

test_that("zero-library samples are flagged by name", {
  m <- matrix(c(0, 0, 5, 7), nrow = 2,
              dimnames = list(c("g1", "g2"), c("empty", "ok")))
  expect_error(cpm(m), "empty")
})
