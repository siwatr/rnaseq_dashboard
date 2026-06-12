test_that("lookup uses <feature_type>_name by default", {
  rd <- data.frame(gene_name = c("Actb", "Gapdh", "Sox2"),
                   row.names = c("ENSG1", "ENSG2", "ENSG3"))
  expect_equal(lookup_feature(c("Sox2", "Actb"), rd, feature_type = "gene"),
               c("ENSG3", "ENSG1"))
})

test_that("lookup falls back to ids when the name column is absent", {
  rd <- data.frame(biotype = c("protein_coding", "lincRNA"),
                   row.names = c("ENSG1", "ENSG2"))
  expect_equal(lookup_feature(c("ENSG2"), rd, feature_type = "gene"), "ENSG2")
})

test_that("an explicit column overrides feature_type", {
  rd <- data.frame(symbol = c("Actb", "Gapdh"),
                   row.names = c("ENSG1", "ENSG2"))
  expect_equal(lookup_feature("Gapdh", rd, column = "symbol"), "ENSG2")
})

test_that("unmatched queries yield NA", {
  rd <- data.frame(gene_name = "Actb", row.names = "ENSG1")
  expect_true(is.na(lookup_feature("Nope", rd, feature_type = "gene")))
})
