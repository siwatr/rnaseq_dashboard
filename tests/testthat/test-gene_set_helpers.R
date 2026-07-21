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

# ---- combine_gene_set_annotation (P7e) ------------------------------------

test_that("combine_gene_set_annotation: single-set membership + level order", {
  sets <- list(up = new_gene_set(c("g1", "g2")), down = new_gene_set(c("g3", "g4")))
  r <- combine_gene_set_annotation(sets)
  expect_equal(unname(r$annotation[c("g1", "g2", "g3", "g4")]),
               c("up", "up", "down", "down"))
  expect_equal(r$levels, c("up", "down"))       # member-set input order
  expect_length(r$shared_ids, 0L)
  expect_identical(r$note, "")
})

test_that("combine_gene_set_annotation: concat is the default for shared genes", {
  sets <- list(up = new_gene_set(c("g1", "g2")), hyp = new_gene_set(c("g2", "g3")))
  r <- combine_gene_set_annotation(sets)               # sep = ";"
  expect_equal(unname(r$annotation["g2"]), "up;hyp")   # input order, joined
  expect_setequal(r$shared_ids, "g2")
  expect_match(r$note, "1 gene")
  expect_equal(r$levels, c("up", "hyp", "up;hyp"))     # singles first, combo last
  # Custom separator.
  expect_equal(unname(combine_gene_set_annotation(sets, sep = "|")$annotation["g2"]),
               "up|hyp")
})

test_that("combine_gene_set_annotation: label + first strategies", {
  sets <- list(A = new_gene_set(c("g1", "g2")), B = new_gene_set(c("g2", "g3")))
  lab <- combine_gene_set_annotation(sets, shared = "label")
  expect_equal(unname(lab$annotation["g2"]), "multiple")
  expect_true("multiple" %in% lab$levels)
  expect_equal(unname(combine_gene_set_annotation(sets, shared = "label",
                                                  label = "both")$annotation["g2"]), "both")
  fst <- combine_gene_set_annotation(sets, shared = "first")
  expect_equal(unname(fst$annotation["g2"]), "A")      # first set in input order
})

test_that("combine_gene_set_annotation: 'label' guards against a member-set-name clash", {
  sets <- list(multiple = new_gene_set("g1"), B = new_gene_set(c("g1", "g2")))
  expect_error(combine_gene_set_annotation(sets, shared = "label"), "collides")
  expect_error(combine_gene_set_annotation(sets, shared = "label", label = ""), "non-empty")
})

test_that("combine_gene_set_annotation: empty / one-set / NA-blank edges", {
  expect_length(combine_gene_set_annotation(list())$annotation, 0L)
  expect_length(combine_gene_set_annotation(NULL)$levels, 0L)
  one <- combine_gene_set_annotation(list(S = new_gene_set(c("g1", "g2"))))
  expect_equal(one$levels, "S")
  expect_length(one$shared_ids, 0L)
  # Blank ids dropped; a set with no ids contributes nothing.
  edge <- combine_gene_set_annotation(list(A = c("g1", "", NA), B = character(0)))
  expect_equal(names(edge$annotation), "g1")
  expect_equal(edge$levels, "A")
})

test_that("gene_set_annotation_composition counts per level, within-aware", {
  set <- new_gene_set(c("g1", "g2", "g3", "g4"), kind = "annotated",
                      annotation = c(g1 = "up", g2 = "up", g3 = "down", g4 = "up"))
  comp <- gene_set_annotation_composition(set)
  expect_setequal(comp$level, c("up", "down"))
  expect_equal(comp$n[comp$level == "up"], 3L)
  comp2 <- gene_set_annotation_composition(set, feature_ids = c("g1", "g3"), within = TRUE)
  expect_equal(comp2$n[comp2$level == "up"], 1L)
  expect_equal(comp2$n[comp2$level == "down"], 1L)
  expect_equal(nrow(gene_set_annotation_composition(new_gene_set("g1"))), 0L)  # no annotation
})

test_that("gene_set_anno_colors: default scheme + palette override (normalized hex)", {
  d <- gene_set_anno_colors(c("up", "down"))
  expect_setequal(names(d), c("up", "down"))
  expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", d)))
  cfg <- list(name = "Custom palette", colors = c(up = "#ff0000", down = "#0000ff"))
  o <- gene_set_anno_colors(c("up", "down"), cfg)
  expect_equal(unname(o[["up"]]), "#FF0000")
  expect_equal(unname(o[["down"]]), "#0000FF")
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

test_that("a bare object with a set named 'gene_sets' is not mistaken for the wrapper", {
  # The wrapper is detected by the version KEY, not the 'gene_sets' name.
  rt <- gene_sets_from_json('{"gene_sets": ["a", "b"], "other": ["c"]}')
  expect_setequal(names(rt), c("gene_sets", "other"))
  expect_equal(rt$gene_sets$ids, c("a", "b"))
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

test_that("TSV round-trip preserves an id literally equal to 'NA'", {
  sets <- list(S = new_gene_set(c("g1", "NA", "g2")))
  rt <- gene_sets_from_tsv(gene_sets_to_tsv(sets))
  expect_equal(rt$S$ids, c("g1", "NA", "g2"))          # the string "NA" is not dropped
})

test_that("gene_sets_from_file sniffs a GMT saved with a .txt extension", {
  sets <- list(Up = new_gene_set(c("g1", "g2"), source = "paste"))
  tf <- tempfile(fileext = ".txt"); writeLines(gene_sets_to_gmt(sets), tf)
  rt <- gene_sets_from_file(tf)
  expect_equal(names(rt), "Up")
  expect_equal(rt$Up$ids, c("g1", "g2"))
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

# ---- Compare tab (P6e): size frame / overlap list / presence colours ------

test_that("gene_set_ids_for honours the within-dataset toggle", {
  s <- new_gene_set(c("g1", "g2", "gX"))
  feats <- c("g1", "g2", "g3")
  expect_setequal(gene_set_ids_for(s, feats, within = FALSE), c("g1", "g2", "gX"))
  expect_setequal(gene_set_ids_for(s, feats, within = TRUE), c("g1", "g2"))
  # Bare id vector accepted too.
  expect_setequal(gene_set_ids_for(c("g1", "gX"), feats, within = TRUE), "g1")
})

test_that("gene_set_size_frame: present+absent (within=FALSE) vs present-only (within=TRUE)", {
  sets <- list(A = new_gene_set(c("g1", "g2", "gX")),   # 2 present, 1 absent
               B = new_gene_set(c("g2", "g3")))          # 2 present, 0 absent
  feats <- c("g1", "g2", "g3")
  fr <- gene_set_size_frame(sets, feats, within = FALSE)
  expect_equal(levels(fr$set), c("A", "B"))             # input order preserved
  expect_equal(levels(fr$status), c("present", "absent"))
  expect_equal(fr$n[fr$set == "A" & fr$status == "present"], 2L)
  expect_equal(fr$n[fr$set == "A" & fr$status == "absent"], 1L)
  expect_equal(fr$n[fr$set == "B" & fr$status == "absent"], 0L)
  # within=TRUE -> one present row per set (no absent rows).
  fw <- gene_set_size_frame(sets, feats, within = TRUE)
  expect_equal(nrow(fw), 2L)
  expect_true(all(fw$status == "present"))
  expect_equal(fw$n, c(2L, 2L))
  # Empty store -> empty (but typed) frame.
  expect_equal(nrow(gene_set_size_frame(list(), feats)), 0L)
})

test_that("gene_set_size_frame orders the set factor by size or name", {
  sets <- list(mid = new_gene_set(c("g1", "g2")),         # size 2
               big = new_gene_set(c("g1", "g2", "g3")),   # size 3
               sm  = new_gene_set("g1"))                   # size 1
  feats <- c("g1", "g2", "g3")
  expect_equal(levels(gene_set_size_frame(sets, feats, order = "none")$set),
               c("mid", "big", "sm"))                      # input order
  expect_equal(levels(gene_set_size_frame(sets, feats, order = "inc")$set),
               c("sm", "mid", "big"))
  expect_equal(levels(gene_set_size_frame(sets, feats, order = "dec")$set),
               c("big", "mid", "sm"))
  expect_equal(levels(gene_set_size_frame(sets, feats, order = "az")$set),
               c("big", "mid", "sm"))
  expect_equal(levels(gene_set_size_frame(sets, feats, order = "za")$set),
               c("sm", "mid", "big"))
})

test_that("gene_set_overlap_list keeps unique ids per set + follows the toggle", {
  sets <- list(A = new_gene_set(c("g1", "g2", "gX")), B = new_gene_set(c("g2", "g3")))
  feats <- c("g1", "g2", "g3")
  full <- gene_set_overlap_list(sets, feats, within = FALSE)
  expect_equal(names(full), c("A", "B"))
  expect_setequal(full$A, c("g1", "g2", "gX"))
  within <- gene_set_overlap_list(sets, feats, within = TRUE)
  expect_setequal(within$A, c("g1", "g2"))              # gX dropped
  expect_length(gene_set_overlap_list(list()), 0L)
})

test_that("gene_set_presence_colors: default scheme + config override", {
  d <- gene_set_presence_colors()
  expect_named(d, c("present", "absent"))
  expect_match(d[["present"]], "^#[0-9A-Fa-f]{6}$")
  ov <- gene_set_presence_colors(list(name = "Custom palette",
                                      colors = c(present = "#111111", absent = "#eeeeee")))
  expect_equal(unname(ov[["present"]]), "#111111")
  expect_equal(unname(toupper(ov[["absent"]])), "#EEEEEE")
})
