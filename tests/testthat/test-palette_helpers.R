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

test_that("palette groups + names + grouped choices enumerate the catalogue", {
  expect_setequal(palette_type_names(),
                  c("Custom", "Qualitative", "Brewer: Qualitative",
                    "Brewer: Sequential", "Brewer: Divergent", "viridis"))
  expect_identical(palette_type_names()[1], "Custom")        # Custom first
  expect_identical(utils::tail(palette_type_names(), 1), "viridis")  # viridis last
  # Qualitative is the dependency-free pair only (Brewer is its own group now).
  expect_setequal(palette_names("Qualitative"), c("ggplot default", "Okabe-Ito"))
  expect_identical(palette_names("Custom"), "Custom palette")
  if (requireNamespace("viridisLite", quietly = TRUE))
    expect_true("viridis: magma" %in% palette_names("viridis"))
  if (requireNamespace("RColorBrewer", quietly = TRUE)) {
    expect_true(all(grepl("^RColorBrewer: ", palette_names("Brewer: Qualitative"))))
    expect_true(length(palette_names("Brewer: Divergent")) > 0)
  }
  # palette_choices(): grouped, with clean labels (no "<pkg>: " prefix) but
  # resolvable values.
  ch <- palette_choices()
  expect_named(ch, palette_type_names())
  expect_true("Okabe-Ito" %in% ch[["Qualitative"]])
  if (requireNamespace("viridisLite", quietly = TRUE)) {
    vir <- ch[["viridis"]]
    expect_true("viridis: magma" %in% unname(vir))           # value resolvable
    expect_true("magma" %in% names(vir))                     # label clean
  }
})

test_that("palette_colors resolves a palette name to n colours, inferring type", {
  expect_length(palette_colors("Okabe-Ito", 3), 3L)
  expect_equal(palette_colors("Okabe-Ito", 2), c("#E69F00", "#56B4E9"))
  expect_length(palette_colors("Okabe-Ito", 10), 10L)        # >8 -> interpolated
  expect_length(palette_colors("ggplot default", 4), 4L)
  expect_equal(palette_colors("nope", 2), palette_colors("Okabe-Ito", 2))  # unknown -> fallback
  expect_length(palette_colors("Okabe-Ito", 0), 0L)
  if (requireNamespace("viridisLite", quietly = TRUE)) {
    v <- palette_colors("viridis: viridis", 5)
    expect_length(v, 5L); expect_true(all(grepl("^#", v)))
  }
  if (requireNamespace("RColorBrewer", quietly = TRUE)) {
    expect_length(palette_colors("RColorBrewer: RdBu", 4), 4L)
    expect_length(palette_colors("RColorBrewer: RdBu", 2), 2L)   # <3 handled
  }
  # "Custom palette" ramps through the supplied anchors (verbatim), else Okabe-Ito.
  expect_equal(palette_colors("Custom palette", 2, custom = c("#000000", "#ffffff")),
               c("#000000", "#ffffff"))
  expect_equal(palette_colors("Custom palette", 2), palette_colors("Okabe-Ito", 2))
})

test_that("palette_discrete fills from the palette, preserves order, normalizes", {
  cols <- palette_discrete(c("b", "a", "c"), name = "Okabe-Ito")
  expect_named(cols, c("b", "a", "c"))            # order preserved (not sorted)
  expect_equal(unname(cols), norm_color(palette_colors("Okabe-Ito", 3)))
  expect_true(all(grepl("^#[0-9A-F]{6}$", cols)))  # 6-digit, normalized
  # Sequential palettes (8-digit alpha hex) are normalized to 6 digits, so a
  # picker echo (#440154) equals the stored value -> no spurious "edit".
  if (requireNamespace("viridisLite", quietly = TRUE)) {
    v <- palette_discrete(c("a", "b"), name = "viridis: viridis")
    expect_true(all(grepl("^#[0-9A-F]{6}$", v)))
  }
})

test_that("palette_discrete honours explicit colours (normalized) over the fill", {
  cols <- palette_discrete(c("control", "treated"), colors = c(treated = "gray50"),
                           name = "Okabe-Ito")
  expect_equal(unname(cols["treated"]), "#7F7F7F")   # explicit + normalized
  expect_equal(unname(cols["control"]), "#E69F00")   # filled from palette
})

test_that("palette_discrete ignores invalid and non-matching colours", {
  cols <- palette_discrete(c("a", "b"), colors = c(a = "garbage", z = "#000000"),
                           name = "Okabe-Ito")
  expect_equal(unname(cols["a"]), "#E69F00")          # invalid -> palette fill
  expect_false("z" %in% names(cols))                  # non-matching dropped
})

test_that("palette_discrete handles empty / factor input", {
  expect_length(palette_discrete(character(0)), 0L)
  f <- factor(c("x", "y", "x"), levels = c("y", "x"))
  expect_named(palette_discrete(f), c("y", "x"))      # factor level order
})
