# Project-wide colour palette engine. Pure, dependency-light helpers that are
# the single source of truth for resolving a *discrete* colour mapping: an
# explicit per-level colour map layered on top of a named base palette. Consumed
# by the Palette page (mod_palette.R), the QC plots + correlation heatmap
# (mod_qc.R via qc_annotation_colors()), and later the PCA/DE/heatmap pages.
#
# A palette is identified by a (type, name) pair:
#   - type: "Qualitative" / "Sequential" / "Divergent" / "Custom"
#   - name: a palette within that type. Package palettes are "<pkg>: <name>"
#     (e.g. "viridis: magma", "RColorBrewer: Set2"); "Qualitative" also has the
#     dependency-free "ggplot default" and "Okabe-Ito"; "Custom" has the single
#     "Custom palette" (a user-defined ramp, default Okabe-Ito).
# Sequential/divergent palettes are usable on discrete data (n colours sampled
# from the ramp). viridisLite / RColorBrewer are optional (Suggests): when
# absent, only the built-in qualitative palettes are offered and everything
# falls back to Okabe-Ito.

# Okabe-Ito: an 8-colour, colour-blind-safe qualitative palette (no dependency).
.okabe_ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
                "#0072B2", "#D55E00", "#CC79A7", "#999999")

# viridisLite options offered under "Sequential".
.viridis_options <- c("viridis", "magma", "plasma", "inferno",
                      "cividis", "mako", "rocket", "turbo")

# ggplot2's default discrete colours == scales::hue_pal()(n). Reproduced with
# grDevices::hcl so we take no scales dependency.
.gg_hue <- function(n) {
  if (n < 1L) return(character(0))
  h <- seq(15, 375, length.out = n + 1L)[seq_len(n)]
  grDevices::hcl(h = h, c = 100, l = 65)
}

# RColorBrewer palette names of a given category ("qual"/"seq"/"div"), prefixed
# "RColorBrewer: "; empty when the package is unavailable.
.brewer_names <- function(category) {
  if (!requireNamespace("RColorBrewer", quietly = TRUE)) return(character(0))
  ci <- RColorBrewer::brewer.pal.info
  paste0("RColorBrewer: ", rownames(ci)[ci$category == category])
}

#' Palette types and the palette names within each
#'
#' @return `palette_type_names()`: the available types. `palette_names(type)`:
#'   the palette names selectable for that type (package palettes formatted
#'   `"<pkg>: <name>"`). Names depend on which optional packages are installed.
#' @export
palette_type_names <- function() c("Qualitative", "Sequential", "Divergent", "Custom")

#' @rdname palette_type_names
#' @param type One of [palette_type_names()].
#' @export
palette_names <- function(type = "Qualitative") {
  viridis <- if (requireNamespace("viridisLite", quietly = TRUE))
    paste0("viridis: ", .viridis_options) else character(0)
  switch(type,
    Qualitative = c("ggplot default", "Okabe-Ito", .brewer_names("qual")),
    Sequential  = c(viridis, .brewer_names("seq")),
    Divergent   = .brewer_names("div"),
    Custom      = "Custom palette",
    character(0))
}

#' Grouped palette choices for a `selectInput`
#'
#' A named list (type -> palette names) suitable as `selectInput(choices = ...)`,
#' rendering each type as an `<optgroup>`. The selected value is the palette
#' *name* (unique across types), which is all [palette_colors()] /
#' [palette_discrete()] need.
#' @return A named list of character vectors.
#' @export
palette_choices <- function() {
  types <- palette_type_names()
  stats::setNames(lapply(types, palette_names), types)
}

#' Names of the built-in qualitative palettes (back-compat shim)
#' @return Character vector of qualitative palette names.
#' @export
palette_qualitative_names <- function() palette_names("Qualitative")

#' Generate `n` colours from the Okabe-Ito palette (or `"ggplot default"`)
#'
#' A thin helper retained for internal use and the Reference swatches; general
#' palette generation goes through [palette_colors()].
#' @param name `"Okabe-Ito"` or `"ggplot default"`.
#' @param n Number of colours.
#' @return `n` hex colours (interpolated past 8 for Okabe-Ito).
#' @export
palette_qualitative <- function(name = "Okabe-Ito", n) {
  n <- as.integer(n)
  if (is.na(n) || n < 1L) return(character(0))
  if (identical(name, "ggplot default")) return(.gg_hue(n))
  if (n <= length(.okabe_ito)) .okabe_ito[seq_len(n)]
  else grDevices::colorRampPalette(.okabe_ito)(n)
}

#' Generate `n` colours from a palette name
#'
#' Everything needed is inferred from `name` (so a single selector value drives
#' it): `"Custom palette"` ramps through `custom` (default Okabe-Ito); package
#' palettes are `"viridis: <opt>"` / `"RColorBrewer: <pal>"`; plus the
#' dependency-free `"ggplot default"` and `"Okabe-Ito"`. Sequential/divergent
#' ramps are sampled to `n` discrete colours; qualitative palettes are taken
#' verbatim and interpolated only past their stops. Unknown palettes / missing
#' optional packages fall back to Okabe-Ito.
#'
#' @param name A palette name from [palette_names()].
#' @param n Number of colours required.
#' @param custom For `"Custom palette"`, anchor colours to ramp through (default
#'   Okabe-Ito).
#' @return A character vector of `n` hex colours.
#' @export
palette_colors <- function(name = "Okabe-Ito", n, custom = NULL) {
  n <- as.integer(n)
  if (is.na(n) || n < 1L) return(character(0))
  ramp <- function(stops) {
    stops <- stops[!is.na(stops)]
    if (!length(stops)) stops <- .okabe_ito
    if (n <= length(stops)) stops[seq_len(n)] else grDevices::colorRampPalette(stops)(n)
  }
  if (identical(name, "Custom palette")) return(ramp(if (length(custom)) custom else .okabe_ito))
  if (identical(name, "ggplot default")) return(.gg_hue(n))
  if (identical(name, "Okabe-Ito")) return(ramp(.okabe_ito))
  parts <- strsplit(name, ": ", fixed = TRUE)[[1]]
  pkg <- parts[1L]; pal <- if (length(parts) > 1L) parts[2L] else parts[1L]
  if (identical(pkg, "viridis") && requireNamespace("viridisLite", quietly = TRUE)) {
    return(viridisLite::viridis(n, option = pal))
  }
  if (identical(pkg, "RColorBrewer") && requireNamespace("RColorBrewer", quietly = TRUE)) {
    ci <- RColorBrewer::brewer.pal.info
    if (pal %in% rownames(ci)) {
      mx <- ci[pal, "maxcolors"]
      base <- RColorBrewer::brewer.pal(max(3L, min(mx, n)), pal)
      return(if (n <= length(base)) base[seq_len(n)] else grDevices::colorRampPalette(base)(n))
    }
  }
  ramp(.okabe_ito)                                   # fallback
}

#' Normalize a colour string to hex
#'
#' Accepts any colour R understands — hex (`"#1f77b4"`), R names (`"gray50"`),
#' or CSS names (`"steelblue"`) — and returns a `#RRGGBB` string, so the same
#' value can feed both ggplot scales and a JS colour picker. Invalid/empty input
#' yields `NA`.
#'
#' @param x Character vector of colour strings.
#' @return Character vector of `#RRGGBB` hex colours (or `NA`), same length as `x`.
#' @export
norm_color <- function(x) {
  vapply(as.character(x), function(v) {
    if (is.na(v) || !nzchar(trimws(v))) return(NA_character_)
    rgb <- tryCatch(grDevices::col2rgb(trimws(v)), error = function(e) NULL)
    if (is.null(rgb)) NA_character_
    else grDevices::rgb(rgb[1L], rgb[2L], rgb[3L], maxColorValue = 255)
  }, character(1), USE.NAMES = FALSE)
}

#' Resolve a discrete level -> colour mapping
#'
#' The single resolver behind every discrete colour scale in the app: an explicit
#' per-level colour map (`colors`) layered over an auto-filled base palette, so
#' the QC ggplots and the ComplexHeatmap annotations agree. Deterministic and
#' stable across re-renders (level order preserved). Levels without an explicit
#' colour are filled from the `(type, name)` palette — so a new dataset's levels
#' still get sensible colours.
#'
#' @param levels Character/factor levels to colour.
#' @param colors Named character vector (`level = colour`) of explicit colours;
#'   may be hex / R names / CSS names. Invalid or non-matching entries are
#'   ignored. `NULL` for none.
#' @param name The base palette name (see [palette_names()]) for unmapped levels.
#' @param custom Custom-ramp anchors when `name = "Custom palette"`.
#' @return A named character vector (`level = #RRGGBB`) over all `levels`, in
#'   their original order. All colours are normalized to 6-digit hex (via
#'   [norm_color()]) so equality checks against picker echoes are stable.
#' @export
palette_discrete <- function(levels, colors = NULL, name = "Okabe-Ito",
                             custom = NULL) {
  levels <- if (is.factor(levels)) levels(levels) else as.character(levels)
  levels <- levels[!is.na(levels)]
  if (!length(levels)) return(stats::setNames(character(0), character(0)))
  out <- stats::setNames(norm_color(palette_colors(name, length(levels), custom)), levels)
  if (!is.null(colors) && length(colors)) {
    m <- norm_color(colors)
    names(m) <- names(colors)
    hit <- intersect(levels, names(m)[!is.na(m)])
    out[hit] <- m[hit]
  }
  out
}
