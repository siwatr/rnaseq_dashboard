# Project-wide colour palette engine. Pure, dependency-light helpers that are
# the single source of truth for resolving a *discrete* colour mapping: explicit
# per-level pins layered on top of a named base palette. Consumed by the Palette
# page (mod_palette.R), the QC plots + correlation heatmap (mod_qc.R via
# qc_annotation_colors()), and later the PCA/DE/heatmap pages. The continuous
# resolver + config import/export land in P3g-b.

# Okabe-Ito: an 8-colour, colour-blind-safe qualitative palette (no dependency).
.okabe_ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
                "#0072B2", "#D55E00", "#CC79A7", "#999999")

# Built-in qualitative palettes (fixed hex stops; interpolated when a column has
# more levels than stops). "ggplot default" is generated on the fly (see below)
# so it always yields exactly n evenly-spaced hues.
.palette_registry <- list(
  "Okabe-Ito"          = .okabe_ito,
  "Viridis (discrete)" = c("#440154", "#414487", "#2A788E",
                           "#22A884", "#7AD151", "#FDE725"),
  "Set2"  = c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3",
              "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3"),
  "Dark2" = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A",
              "#66A61E", "#E6AB02", "#A6761D", "#666666")
)

# ggplot2's default discrete colours == scales::hue_pal()(n). Reproduced with
# grDevices::hcl so we take no scales dependency.
.gg_hue <- function(n) {
  if (n < 1L) return(character(0))
  h <- seq(15, 375, length.out = n + 1L)[seq_len(n)]
  grDevices::hcl(h = h, c = 100, l = 65)
}

#' Names of the built-in qualitative palettes
#'
#' @return Character vector of palette names accepted by [palette_qualitative()]
#'   and [palette_discrete()].
#' @export
palette_qualitative_names <- function() c("ggplot default", names(.palette_registry))

#' Generate `n` colours from a named qualitative palette
#'
#' @param name One of [palette_qualitative_names()]. Unknown names fall back to
#'   `"Okabe-Ito"`.
#' @param n Number of colours required.
#' @return A character vector of `n` hex colours. Fixed palettes are interpolated
#'   (via [grDevices::colorRampPalette()]) when `n` exceeds the number of stops;
#'   `"ggplot default"` always yields `n` evenly-spaced hues.
#' @export
palette_qualitative <- function(name = "ggplot default", n) {
  n <- as.integer(n)
  if (is.na(n) || n < 1L) return(character(0))
  if (identical(name, "ggplot default")) return(.gg_hue(n))
  stops <- .palette_registry[[name]]
  if (is.null(stops)) stops <- .okabe_ito
  if (n <= length(stops)) stops[seq_len(n)]
  else grDevices::colorRampPalette(stops)(n)
}

# Whether `n` levels exceed a palette's distinct stops (so colours get
# interpolated rather than taken verbatim). "ggplot default" never overflows.
.palette_overflows <- function(name, n) {
  if (identical(name, "ggplot default")) return(FALSE)
  stops <- .palette_registry[[name]]
  if (is.null(stops)) stops <- .okabe_ito
  as.integer(n) > length(stops)
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
#' The single resolver behind every discrete colour scale in the app: explicit
#' per-level pins (`mapping`) override an auto-filled base `palette`, so the QC
#' ggplots and the ComplexHeatmap annotations agree. Deterministic and stable
#' across re-renders (level order preserved).
#'
#' @param levels Character/factor levels to colour.
#' @param mapping Named character vector of pins (`level = colour`); colours may
#'   be hex / R names / CSS names (normalized via [norm_color()]). Invalid or
#'   non-matching pins are ignored. `NULL` for none.
#' @param palette Base palette name (see [palette_qualitative_names()]) used to
#'   fill levels without a pin.
#' @return A named character vector (`level = #RRGGBB`) over all `levels`, in
#'   their original order.
#' @export
palette_discrete <- function(levels, mapping = NULL, palette = "ggplot default") {
  levels <- if (is.factor(levels)) levels(levels) else as.character(levels)
  levels <- levels[!is.na(levels)]
  if (!length(levels)) return(stats::setNames(character(0), character(0)))
  out <- stats::setNames(palette_qualitative(palette, length(levels)), levels)
  if (!is.null(mapping) && length(mapping)) {
    m <- norm_color(mapping)
    names(m) <- names(mapping)
    hit <- intersect(levels, names(m)[!is.na(m)])
    out[hit] <- m[hit]
  }
  out
}
