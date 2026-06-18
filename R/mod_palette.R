# Project-wide colour palette configuration (mock / placeholder).
#
# Intent (not yet wired): one page where the user sets colour conventions for
# the whole project, so every plot stays consistent. It will be the single
# source of truth that (a) configures `thematic`'s qualitative/sequential
# palettes for the ggplot-based plots (QC, PCA, DE), and (b) supplies the `col`
# mappings for ComplexHeatmap annotations (see qc_annotation_colors()). Designed
# to cooperate with thematic, not fight it: thematic keeps theming ggplot colour
# scales, this page just chooses which palette/level-colours it uses. Pinning a
# specific metadata level to a specific colour (a lab convention) is the part
# that goes beyond thematic and will live here.

mod_palette_ui <- function(id) {
  ns <- NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      title = tags$h4("Palette", class = "fs-6 mb-0"), width = 280,
      helpText("Project-wide colour configuration. Coming soon.")
    ),
    bslib::card(
      bslib::card_header(tags$h3("Project palette", class = "fs-6 mb-0")),
      bslib::card_body(
        tags$p(class = "text-muted",
               paste("Set colour conventions for discrete metadata (e.g. condition ->",
                     "fixed colours) and for continuous scales here. These will feed",
                     "both the ggplot plots (via thematic) and the ComplexHeatmap",
                     "annotations, so colours stay consistent across every plot.")),
        tags$p(class = "text-muted mb-0",
               "This page is a placeholder; configuration controls land in a later slice.")
      )
    )
  )
}

#' @param state the shared app-state object (see [new_app_state()]).
#' @return Invisible NULL.
mod_palette_server <- function(id, state) {
  moduleServer(id, function(input, output, session) {
    # TODO: project palette controls; persist into state$meta and apply via
    # thematic (qualitative/sequential) + qc_annotation_colors().
    invisible(NULL)
  })
}
