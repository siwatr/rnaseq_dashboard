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

test_that("lookup_feature can match case-insensitively", {
  rd <- data.frame(gene_name = c("Actb", "Gapdh"), row.names = c("ENSG1", "ENSG2"))
  expect_true(is.na(lookup_feature("actb", rd, feature_type = "gene")))
  expect_equal(lookup_feature("actb", rd, feature_type = "gene",
                              case_insensitive = TRUE), "ENSG1")
})

test_that("resolve_feature returns the first id and a match count", {
  vals <- c("Actb", "Gapdh", "Actb")
  ids  <- c("ENSG1", "ENSG2", "ENSG3")
  r <- resolve_feature("Actb", vals, ids)
  expect_equal(r$id, "ENSG1")     # first match
  expect_equal(r$n, 2L)           # duplicated name
  expect_equal(resolve_feature("Gapdh", vals, ids)$n, 1L)
  miss <- resolve_feature("Nope", vals, ids)
  expect_true(is.na(miss$id)); expect_equal(miss$n, 0L)
})

test_that("resolve_feature case-insensitive folds query and values", {
  r <- resolve_feature("duxf3", c("Duxf3", "Actb"), c("ID1", "ID2"),
                       case_insensitive = TRUE)
  expect_equal(r$id, "ID1"); expect_equal(r$n, 1L)
  expect_true(is.na(resolve_feature("duxf3", c("Duxf3"), "ID1")$id))  # sensitive miss
})

test_that("suggest_features ranks prefix before word-start before substring", {
  vals <- c("Duxf3", "Cux1", "Reduxin", "Aldux")  # all contain 'dux' (ci)
  sg <- suggest_features("dux", vals)
  expect_false(sg$over_cap)
  # Duxf3 (prefix), then Aldux/Reduxin (substring, ranked by length then alpha)
  expect_equal(sg$suggestions[1], "Duxf3")
  expect_true(all(c("Duxf3", "Reduxin", "Aldux") %in% sg$suggestions))
  expect_false("Cux1" %in% sg$suggestions)         # no 'dux'
})

test_that("suggest_features is always case-insensitive and honours min length", {
  expect_equal(suggest_features("d", c("Duxf3"))$suggestions, character(0))  # < 2 chars
  expect_equal(suggest_features("DUX", c("Duxf3"))$suggestions, "Duxf3")
})

test_that("suggest_features caps too-broad searches and de-dups display values", {
  many <- paste0("Gene", sprintf("%03d", 1:150))
  sg <- suggest_features("gene", many, cap = 100)
  expect_true(sg$over_cap)
  expect_equal(sg$suggestions, character(0))
  dup <- suggest_features("act", c("Actb", "Actb", "Acta1"))
  expect_equal(dup$n_match, 2L)                     # distinct only
  expect_equal(length(dup$suggestions), 2L)
})

test_that("feature_search_choices drops logical cols, keeps numeric, adds Feature ID", {
  rd <- data.frame(gene_name = c("Actb", "Gapdh"), entrez = c(11461L, 14433L),
                   is_mito = c(FALSE, TRUE), row.names = c("ENSG1", "ENSG2"))
  ch <- feature_search_choices(rd)
  expect_true("gene_name" %in% ch)
  expect_true("entrez" %in% ch)                     # numeric id kept
  expect_false("is_mito" %in% ch)                   # logical dropped
  expect_equal(unname(ch["Feature ID (row names)"]), "__rownames__")
})

test_that("feature_search_choices can hide duplicate-valued columns", {
  rd <- data.frame(gene_name = c("Actb", "Actb"), uniq = c("a", "b"),
                   row.names = c("ENSG1", "ENSG2"))
  with_dup <- feature_search_choices(rd, include_duplicates = TRUE)
  no_dup   <- feature_search_choices(rd, include_duplicates = FALSE)
  expect_true("gene_name" %in% with_dup)
  expect_false("gene_name" %in% no_dup)             # duplicated values hidden
  expect_true("uniq" %in% no_dup)
  expect_true("__rownames__" %in% no_dup)           # Feature ID always present
})
