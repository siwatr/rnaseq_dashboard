test_that("norm_color accepts hex, R names, and CSS names; rejects garbage", {
  expect_equal(norm_color("#1f77b4"), "#1F77B4")
  expect_equal(norm_color("gray50"), "#7F7F7F")     # R-specific numbered grey
  expect_equal(norm_color("steelblue"), norm_color("#4682B4"))
  expect_true(is.na(norm_color("not-a-colour")))
  expect_true(is.na(norm_color("")))
  expect_true(is.na(norm_color(NA_character_)))
  # vectorized, length-preserving
  out <- norm_color(c("red", "bogus", "#00FF00"))
  expect_length(out, 3L)
  expect_equal(out[1], "#FF0000")
  expect_true(is.na(out[2]))
})

test_that("palette_qualitative_names lists the built-ins with ggplot default first", {
  nm <- palette_qualitative_names()
  expect_identical(nm[1], "ggplot default")
  expect_true(all(c("Okabe-Ito", "Viridis (discrete)", "Set2", "Dark2") %in% nm))
})

test_that("palette_qualitative returns n colours and interpolates past the stops", {
  expect_length(palette_qualitative("Okabe-Ito", 3), 3L)
  expect_equal(palette_qualitative("Okabe-Ito", 2), c("#E69F00", "#56B4E9"))
  # 10 levels > 8 Okabe-Ito stops -> interpolated, still 10 distinct-ish colours
  expect_length(palette_qualitative("Okabe-Ito", 10), 10L)
  # ggplot default reproduces scales::hue_pal()(n) without a scales dependency
  expect_length(palette_qualitative("ggplot default", 4), 4L)
  expect_true(all(grepl("^#", palette_qualitative("ggplot default", 4))))
  # unknown palette falls back to Okabe-Ito
  expect_equal(palette_qualitative("nope", 2), palette_qualitative("Okabe-Ito", 2))
  expect_length(palette_qualitative("Okabe-Ito", 0), 0L)
})

test_that("palette_discrete fills from the palette and preserves level order", {
  cols <- palette_discrete(c("b", "a", "c"), palette = "Okabe-Ito")
  expect_named(cols, c("b", "a", "c"))            # order preserved (not sorted)
  expect_equal(unname(cols), palette_qualitative("Okabe-Ito", 3))
  expect_true(all(grepl("^#", cols)))
})

test_that("palette_discrete honours pins (normalized) over the auto fill", {
  cols <- palette_discrete(c("control", "treated"),
                           mapping = c(treated = "gray50"), palette = "Okabe-Ito")
  expect_equal(unname(cols["treated"]), "#7F7F7F")   # pinned + normalized
  expect_equal(unname(cols["control"]), "#E69F00")   # auto-filled
})

test_that("palette_discrete ignores invalid and non-matching pins", {
  cols <- palette_discrete(c("a", "b"),
                           mapping = c(a = "garbage", z = "#000000"), palette = "Okabe-Ito")
  expect_equal(unname(cols["a"]), "#E69F00")          # invalid pin -> auto fill
  expect_false("z" %in% names(cols))                  # non-matching pin dropped
})

test_that("palette_discrete handles empty / factor input", {
  expect_length(palette_discrete(character(0)), 0L)
  f <- factor(c("x", "y", "x"), levels = c("y", "x"))
  expect_named(palette_discrete(f), c("y", "x"))      # factor level order
})
