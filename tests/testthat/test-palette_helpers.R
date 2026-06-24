test_that("norm_color accepts hex, R names, and CSS names; rejects garbage", {
  expect_equal(norm_color("#1f77b4"), "#1F77B4")
  expect_equal(norm_color("gray50"), "#7F7F7F")     # R-specific numbered grey
  expect_equal(norm_color("steelblue"), norm_color("#4682B4"))
  expect_true(is.na(norm_color("not-a-colour")))
  expect_true(is.na(norm_color("")))
  expect_true(is.na(norm_color(NA_character_)))
  out <- norm_color(c("red", "bogus", "#00FF00"))
  expect_length(out, 3L)
  expect_equal(out[1], "#FF0000")
  expect_true(is.na(out[2]))
})

test_that("palette types + names enumerate the catalogue", {
  expect_setequal(palette_type_names(),
                  c("Qualitative", "Sequential", "Divergent", "Custom"))
  q <- palette_names("Qualitative")
  expect_identical(q[1], "ggplot default")
  expect_true("Okabe-Ito" %in% q)
  expect_identical(palette_names("Custom"), "Custom palette")
  # Package palettes are formatted "<pkg>: <name>" when the package is present.
  if (requireNamespace("viridisLite", quietly = TRUE))
    expect_true("viridis: magma" %in% palette_names("Sequential"))
  if (requireNamespace("RColorBrewer", quietly = TRUE)) {
    expect_true(any(grepl("^RColorBrewer: ", palette_names("Qualitative"))))
    expect_true(any(grepl("^RColorBrewer: ", palette_names("Divergent"))))
  }
})

test_that("palette_colors returns n colours across types and interpolates", {
  expect_length(palette_colors("Qualitative", "Okabe-Ito", 3), 3L)
  expect_equal(palette_colors("Qualitative", "Okabe-Ito", 2), c("#E69F00", "#56B4E9"))
  expect_length(palette_colors("Qualitative", "Okabe-Ito", 10), 10L)   # >8 -> interp
  expect_length(palette_colors("Qualitative", "ggplot default", 4), 4L)
  expect_equal(palette_colors("Qualitative", "nope", 2),
               palette_colors("Qualitative", "Okabe-Ito", 2))          # unknown -> fallback
  expect_length(palette_colors("Qualitative", "Okabe-Ito", 0), 0L)
  # Sequential / divergent sampled to n discrete colours.
  if (requireNamespace("viridisLite", quietly = TRUE)) {
    v <- palette_colors("Sequential", "viridis: viridis", 5)
    expect_length(v, 5L); expect_true(all(grepl("^#", v)))
  }
  if (requireNamespace("RColorBrewer", quietly = TRUE)) {
    b <- palette_colors("Divergent", "RColorBrewer: RdBu", 4)
    expect_length(b, 4L); expect_true(all(grepl("^#", b)))
    expect_length(palette_colors("Divergent", "RColorBrewer: RdBu", 2), 2L)  # <3 handled
  }
  # Custom uses the supplied anchors verbatim (normalization happens downstream
  # in palette_discrete), else Okabe-Ito.
  expect_equal(palette_colors("Custom", "Custom palette", 2, custom = c("#000000", "#ffffff")),
               c("#000000", "#ffffff"))
  expect_equal(palette_colors("Custom", "Custom palette", 2),
               palette_colors("Qualitative", "Okabe-Ito", 2))
})

test_that("palette_discrete fills from the palette and preserves level order", {
  cols <- palette_discrete(c("b", "a", "c"), type = "Qualitative", name = "Okabe-Ito")
  expect_named(cols, c("b", "a", "c"))            # order preserved (not sorted)
  expect_equal(unname(cols), palette_colors("Qualitative", "Okabe-Ito", 3))
  expect_true(all(grepl("^#", cols)))
})

test_that("palette_discrete honours explicit colours (normalized) over the fill", {
  cols <- palette_discrete(c("control", "treated"), colors = c(treated = "gray50"),
                           type = "Qualitative", name = "Okabe-Ito")
  expect_equal(unname(cols["treated"]), "#7F7F7F")   # explicit + normalized
  expect_equal(unname(cols["control"]), "#E69F00")   # filled from palette
})

test_that("palette_discrete ignores invalid and non-matching colours", {
  cols <- palette_discrete(c("a", "b"), colors = c(a = "garbage", z = "#000000"),
                           type = "Qualitative", name = "Okabe-Ito")
  expect_equal(unname(cols["a"]), "#E69F00")          # invalid -> palette fill
  expect_false("z" %in% names(cols))                  # non-matching dropped
})

test_that("palette_discrete handles empty / factor input", {
  expect_length(palette_discrete(character(0)), 0L)
  f <- factor(c("x", "y", "x"), levels = c("y", "x"))
  expect_named(palette_discrete(f), c("y", "x"))      # factor level order
})
