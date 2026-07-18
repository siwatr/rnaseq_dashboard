# Reusable "expression value" control -- assay + transform + pseudocount -- for
# any plot that turns an assay into a displayed value. Used by the DE direct-
# comparison plot now and intended for the P7 Expression page (single-gene +
# gene-set), so the three controls stay consistent across pages. Like
# mod_plot_subset, it operates in the *host* module's namespace (ui builds the
# inputs, server wires them + returns a reactive spec): a page calls
# expr_value_ui(ns, suffix) in a sidebar and expr_value_server(input, output,
# session, state, suffix) once. `suffix` distinguishes instances on a page.

# The three inputs. Embed inside an accordion panel / conditionalPanel as needed.
# Pseudocount only shows for a log transform (undefined for "none").
expr_value_ui <- function(ns, suffix, assay_label = "Expression value (assay)") {
  tp <- ns(paste0(suffix, "_transform"))
  tagList(
    selectInput(ns(paste0(suffix, "_assay")), assay_label, choices = NULL),
    selectInput(ns(paste0(suffix, "_transform")), "Transform",
                c("None" = "none", "log2" = "log2", "log10" = "log10"),
                selected = "none"),
    conditionalPanel(
      sprintf("input['%s'] != 'none'", tp),
      numericInput(ns(paste0(suffix, "_pc")), "Pseudocount", value = 1, min = 0, step = 0.5))
  )
}

# Populate the assay choices from the working dds (default logcounts), auto-pick a
# sensible transform when the assay changes (none for a log assay, log10 else),
# and return a reactive spec the caller feeds to de_group_means() / de_transform_matrix().
#
# `include_vst` prepends a synthetic "VST" choice (resolved by the caller via
# expr_value_matrix()/qc_vst()); `default_fn` (a function(dds) -> value key, e.g.
# expr_default_assay) picks the initial selection. Both default off so DE's use is
# unchanged.
expr_value_server <- function(input, output, session, state, suffix,
                              include_vst = FALSE, default_fn = NULL) {
  a_id <- paste0(suffix, "_assay")
  t_id <- paste0(suffix, "_transform")
  pc_id <- paste0(suffix, "_pc")

  observeEvent(state$working, {
    dds <- state$working; if (is.null(dds)) return()
    an <- SummarizedExperiment::assayNames(dds)
    if (!length(an)) an <- "counts"
    choices <- an
    if (isTRUE(include_vst)) choices <- c("VST" = "vst", stats::setNames(an, an))
    keys <- unname(choices)
    sel <- if (!is.null(input[[a_id]]) && input[[a_id]] %in% keys) input[[a_id]]
           else if (!is.null(default_fn)) {
             d <- tryCatch(default_fn(dds), error = function(e) NULL)
             if (!is.null(d) && d %in% keys) d
             else if ("logcounts" %in% keys) "logcounts" else keys[1]
           }
           else if ("logcounts" %in% keys) "logcounts" else keys[1]
    updateSelectInput(session, a_id, choices = choices, selected = sel)
  }, ignoreNULL = FALSE)

  # An assay that already looks log-scale (logcounts, VST, normalized log-counts)
  # defaults to no transform; linear abundance assays to log10.
  observeEvent(input[[a_id]], {
    a <- input[[a_id]]; if (is.null(a) || !nzchar(a)) return()
    is_log <- grepl("log|vst", a, ignore.case = TRUE)
    updateSelectInput(session, t_id, selected = if (is_log) "none" else "log10")
  }, ignoreInit = TRUE)

  reactive({
    a <- input[[a_id]] %||% "logcounts"
    tf <- input[[t_id]] %||% "none"
    pc <- suppressWarnings(as.numeric(input[[pc_id]] %||% 1)); if (is.na(pc)) pc <- 1
    list(assay = a, transform = tf, pseudocount = pc,
         label = expr_value_label(a, tf, pc))
  })
}

#' Axis / legend label for an expression-value spec
#'
#' `"none"` shows the bare assay name; a log transform shows `"<transform>(<assay>
#' + <pseudocount>)"` so a log axis is never mistaken for a linear one.
#' @param assay Assay name.
#' @param transform `"none"`, `"log2"`, or `"log10"`.
#' @param pseudocount Pseudocount used for the log transform.
#' @return A single label string.
#' @export
expr_value_label <- function(assay, transform = "none", pseudocount = 1) {
  if (identical(transform, "none")) return(assay)
  sprintf("%s(%s + %s)", transform, assay, format(pseudocount))
}
