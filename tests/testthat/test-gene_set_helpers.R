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

test_that("split_ids_by_group keys groups by the annotation column(s)", {
  df <- data.frame(id = c("a", "b", "c", "d"),
                   dir = c("up", "up", "down", "down"),
                   grp = c("x", "y", "x", "y"), stringsAsFactors = FALSE)
  # No annotation columns -> one "All" group.
  expect_equal(split_ids_by_group(df, "id"), list(All = c("a", "b", "c", "d")))
  # One column -> one set per level.
  s1 <- split_ids_by_group(df, "id", "dir")
  expect_setequal(names(s1), c("up", "down"))
  expect_equal(s1$up, c("a", "b")); expect_equal(s1$down, c("c", "d"))
  # Two columns -> combined key with the separator.
  s2 <- split_ids_by_group(df, "id", c("dir", "grp"), sep = ".")
  expect_setequal(names(s2), c("up.x", "up.y", "down.x", "down.y"))
  expect_equal(s2[["up.x"]], "a")
  expect_equal(names(split_ids_by_group(df, "id", c("dir", "grp"), sep = "_"))[1], "down_x")
})

test_that("split_ids_by_group drops NA/blank ids + empty groups; keys NA annotations", {
  df <- data.frame(id = c("a", NA, "", "b", "a"),
                   dir = c("up", "up", "up", NA, "up"), stringsAsFactors = FALSE)
  s <- split_ids_by_group(df, "id", "dir")
  expect_equal(s$up, c("a"))              # NA/blank dropped, deduped
  expect_equal(s[["NA"]], "b")            # NA annotation keys as "NA"
  # An all-blank table yields nothing.
  expect_length(split_ids_by_group(data.frame(id = c(NA, "")), "id"), 0L)
  expect_error(split_ids_by_group(df, "nope"), "not in the table")
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

# ---- File round-trip (P6d): JSON / GMT / TSV serializers ------------------

test_that("JSON round-trips ids, names, kind, source (faithful)", {
  sets <- list(Up   = new_gene_set(c("g1", "g2", "g3"), source = "DE: A_vs_B"),
               Solo = new_gene_set("g9", source = "paste"))
  rt <- gene_sets_from_json(gene_sets_to_json(sets))
  expect_equal(names(rt), names(sets))
  expect_equal(lapply(rt, `[[`, "ids"), lapply(sets, `[[`, "ids"))
  expect_equal(rt$Up$source, "DE: A_vs_B")            # source survives
  expect_equal(rt$Solo$ids, "g9")                      # a 1-gene set stays an array
  expect_equal(rt$Up$kind, "simple")
})

test_that("gene_sets_to_json emits an empty object for an empty store", {
  js <- gene_sets_to_json(list())
  expect_match(js, "\"gene_sets\": \\{\\}")
  expect_length(gene_sets_from_json(js), 0L)
})

test_that("gene_sets_from_json tolerates a bare {name: [ids]} object", {
  rt <- gene_sets_from_json('{"S1": ["a", "b"], "S2": ["c"]}')
  expect_equal(names(rt), c("S1", "S2"))
  expect_equal(rt$S1$ids, c("a", "b"))
  expect_equal(rt$S1$kind, "simple")
})

test_that("GMT round-trips ids/names; source goes to the description field", {
  sets <- list(Up = new_gene_set(c("g1", "g2"), source = "DE: A_vs_B"),
               Dn = new_gene_set("g3", source = "paste"))
  gmt <- gene_sets_to_gmt(sets)
  expect_match(strsplit(gmt, "\n")[[1]][1], "^Up\tDE: A_vs_B\tg1\tg2$")
  rt <- gene_sets_from_gmt(gmt)
  expect_equal(names(rt), names(sets))
  expect_equal(lapply(rt, `[[`, "ids"), lapply(sets, `[[`, "ids"))
  expect_equal(rt$Up$source, "DE: A_vs_B")             # description -> source
})

test_that("gene_sets_from_gmt reads a description-less line and skips blanks", {
  rt <- gene_sets_from_gmt(c("SetA\t\tg1\tg2", "", "SetB\tmy desc\tg3"))
  expect_equal(names(rt), c("SetA", "SetB"))
  expect_equal(rt$SetA$ids, c("g1", "g2"))
  expect_equal(rt$SetA$source, "import: gmt")           # blank desc -> default label
  expect_equal(rt$SetB$source, "my desc")
})

test_that("TSV round-trips ids/names via the long set/id form", {
  sets <- list(Up = new_gene_set(c("g1", "g2")), Dn = new_gene_set("g3"))
  tsv <- gene_sets_to_tsv(sets)
  expect_match(strsplit(tsv, "\n")[[1]][1], "^set\tid\tannotation$")
  rt <- gene_sets_from_tsv(tsv)
  expect_equal(names(rt), names(sets))
  expect_equal(lapply(rt, `[[`, "ids"), lapply(sets, `[[`, "ids"))
})

test_that("gene_sets_from_tsv requires set + id columns", {
  expect_error(gene_sets_from_tsv("foo\tbar\nx\ty"), "set.*id")
})

test_that("gene_sets_from_file auto-detects format by extension", {
  sets <- list(Up = new_gene_set(c("g1", "g2"), source = "paste"))
  for (fmt in c("json", "gmt", "tsv")) {
    tf <- tempfile(fileext = paste0(".", fmt))
    writeLines(get(paste0("gene_sets_to_", fmt))(sets), tf)
    rt <- gene_sets_from_file(tf)
    expect_equal(names(rt), "Up", info = fmt)
    expect_equal(rt$Up$ids, c("g1", "g2"), info = fmt)
  }
})

test_that("gene_sets_from_file sniffs content when the extension is unknown", {
  sets <- list(Up = new_gene_set(c("g1", "g2")))
  tf <- tempfile(fileext = ".dat"); writeLines(gene_sets_to_json(sets), tf)
  expect_equal(gene_sets_from_file(tf)$Up$ids, c("g1", "g2"))
})
