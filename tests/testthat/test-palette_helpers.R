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

test_that("palette_continuous_choices offers continuous groups + a custom ramp", {
  ch <- palette_continuous_choices()
  expect_true("Custom" %in% names(ch))
  expect_identical(unname(ch[["Custom"]]), "Custom ramp")
  expect_false("Qualitative" %in% names(ch))     # qualitative excluded
  if (requireNamespace("viridisLite", quietly = TRUE))
    expect_true("viridis: magma" %in% unname(ch[["viridis"]]))
})

test_that("palette_resolve_range handles numbers, p<pct> percentiles, and fallbacks", {
  v <- 1:100
  expect_equal(palette_resolve_range(v, 0, 100), c(0, 100))
  expect_equal(palette_resolve_range(v), c(1, 100))                 # data min/max
  pr <- palette_resolve_range(v, "p10", "p90")
  expect_equal(round(pr), c(11, 90))                                # ~10th/90th pct
  expect_equal(palette_resolve_range(v, NULL, 50), c(1, 50))        # mixed
  # Degenerate / empty inputs stay valid (hi > lo).
  r0 <- palette_resolve_range(numeric(0)); expect_true(r0[2] > r0[1])
  rc <- palette_resolve_range(c(5, 5, 5)); expect_true(rc[2] > rc[1])
  expect_true(is.na(ddsdashboard:::.resolve_anchor("p200", 1:10)))  # invalid pct
})

test_that("palette_colorramp2 + palette_gradientn resolve named + custom ramps", {
  skip_if_not_installed("circlize")
  f <- palette_colorramp2("viridis: viridis", values = 0:10, min = 0, max = 10)
  expect_type(f, "closure")
  expect_match(f(0), "^#")                                          # maps a value to hex
  # Custom ramp uses the supplied anchor colours.
  g <- palette_gradientn("Custom ramp", values = 0:10, custom = c("#000000", "#ffffff"))
  expect_named(g, c("colours", "values", "limits"))
  expect_equal(g$limits, c(0, 10))
  expect_equal(g$values, seq(0, 1, length.out = length(g$colours)))
  expect_true(all(grepl("^#[0-9A-F]{6}$", g$colours)))
})

test_that("continuous reverse flips the ramp; custom ramp uses supplied colours", {
  g  <- palette_gradientn("viridis: viridis", 0:10)
  gr <- palette_gradientn("viridis: viridis", 0:10, reverse = TRUE)
  expect_equal(gr$colours, rev(g$colours))
  c1 <- palette_gradientn("Custom ramp", 0:10, custom = c("#000000", "#ffffff"))
  c2 <- palette_gradientn("Custom ramp", 0:10, custom = c("#000000", "#ffffff"), reverse = TRUE)
  expect_equal(c2$colours, rev(c1$colours))
  expect_equal(c1$colours[1], "#000000")           # ramps low -> high through custom
})

test_that("palette_to_json / palette_from_json round-trip a mixed config", {
  p <- list(
    colData = list(condition = list(name = "Okabe-Ito",
                                    colors = c(control = "#E69F00", treated = "#56B4E9"))),
    assays  = list(logcounts = list(name = "viridis: viridis", min = "", max = "",
                                    reverse = FALSE, custom = c("#FFFFFF", "#000000"))),
    other   = list(correlation = list(name = "RColorBrewer: RdBu", min = "-1", max = "1",
                                      reverse = TRUE, custom = NULL)))
  j <- palette_to_json(p)
  # The discrete colours serialize to a {level: hex} object (named, not an array).
  expect_match(j, '"control"')
  expect_match(j, '"ddsdashboard_palette_version"')
  rt <- palette_from_json(j)
  expect_equal(rt$colData$condition$colors, c(control = "#E69F00", treated = "#56B4E9"))
  expect_equal(rt$assays$logcounts$custom, c("#FFFFFF", "#000000"))
  expect_true(isTRUE(rt$other$correlation$reverse))
  expect_equal(rt$other$correlation$min, "-1")
})

test_that("palette_to_json drops empty domains; empty config -> {}", {
  expect_match(palette_to_json(list()), "\"palette\": {}", fixed = TRUE)
  # An all-empty / unknown-only config also yields {}.
  expect_match(palette_to_json(list(colData = list(), bogus = list(x = 1))),
               "\"palette\": {}", fixed = TRUE)
})

test_that("palette_from_json accepts a bare palette, drops unknown domains, normalizes", {
  bare <- '{"colData":{"x":{"name":"Okabe-Ito","colors":{"a":"gray50"}}},"bogus":{"y":{}}}'
  out <- palette_from_json(bare)
  expect_named(out, "colData")                       # bogus domain dropped
  expect_equal(out$colData$x$colors[["a"]], "#7F7F7F")  # gray50 normalized to hex
})

test_that("palette_from_json drops empty items and errors on invalid JSON", {
  out <- palette_from_json('{"palette":{"colData":{"empty":{}}}}')
  expect_length(out, 0L)                             # nothing usable
  expect_error(palette_from_json("{not json"))
})

# --- DEG palette set (P5c) --------------------------------------------------

test_that("DEG palette set resolves each named scheme in level order", {
  expect_setequal(palette_names("DEG palette"),
                  c("DEG: Pink-Blue", "DEG: Orange-Purple", "DEG: Red-Blue", "DEG: Coral-Teal"))
  expect_true(.pal_type_discrete("DEG palette"))
  expect_equal(palette_colors("DEG: Pink-Blue", 3), c("#B54661", "#235675", "gray80"))
  # palette_discrete maps positionally onto the DEG factor levels + normalizes hex
  cols <- palette_discrete(c("up", "down", "no_change"), NULL, "DEG: Red-Blue")
  expect_equal(names(cols), c("up", "down", "no_change"))
  expect_equal(unname(cols[c("up", "down")]), c("#D62728", "#1F77B4"))
})

test_that("deg_palette_choices offers only the DEG schemes + Custom", {
  ch <- deg_palette_choices()
  expect_equal(names(ch), c("DEG palette", "Custom"))
  expect_setequal(unname(ch[["DEG palette"]]), palette_names("DEG palette"))
  expect_equal(unname(ch[["Custom"]]), "Custom palette")
  # the generic per-item catalogue must NOT include DEG palettes
  expect_false("DEG palette" %in% palette_type_names())
})
