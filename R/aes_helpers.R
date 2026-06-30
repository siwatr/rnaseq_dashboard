# Shared color/annotation attribute catalog + resolver. The single source of
# truth for "what per-sample attributes can colour/group/annotate a plot, and how
# each resolves to values + colours", used by PCA, the QC plots, and the QC
# correlation-heatmap annotation (and future plot pages). Consolidates logic that
# used to be duplicated as PCA's colour_resolve/session_removal, the QC
# group_map/group_colours/sample_aes, and the heatmap .cor_anno_* helpers.
#
# An attribute is identified by a `key`:
#   "<colData col>"     - a sample-metadata column (discrete or continuous)
#   "__qc__<metric>"    - a per-sample QC metric (continuous)
#   "__removal__"       - the QC suggested-removal status (discrete)
#   "__pool__"          - removal-pool membership (discrete)
#   "__gene__"          - gene expression (continuous; PCA-only, resolved by the
#                         caller which owns the gene-search UI)
# These read the shared app-state (palette configs, promoted samp_flags/samp_pool,
# the qc_metrics derived cache) so every page resolves them identically.

# Per-sample QC metric labels (shared with the catalog + the selectors).
.aes_qc_labels <- c(library_size = "Library size", detected = "Detected features",
                    pct_mito = "% mitochondrial", pct_spike = "% spike-in")

# The qc_per_sample_metrics() frame via the shared derived cache (keyed on
# data_version), so colour-by-QC-metric never recomputes scater metrics per page.
.aes_qc_metrics <- function(state) {
  state_derive(state, "qc_metrics", params = list(),
               expr = function() qc_per_sample_metrics(state$working))
}

# Fixed pool-membership colours (the default when no Palette "Other" config).
.aes_pool_colors <- c(Kept = "#9aa0a6", "In removal pool" = "#D62728")

# Per-sample spike-in (ERCC) dose-response metrics offered as continuous "Spike-in"
# attributes (from spike_dose_response()$per_sample). Available only when the
# dataset has spike-in features.
.aes_spike_labels <- c(slope = "Dose-response slope", r_squared = "Dose-response R-squared",
                       lod = "Lowest detected conc", n_spike_detected = "Detected spike features")

# Does the dataset have spike-in features? (gates the spike attributes.)
.aes_has_spike <- function(state) {
  !is.null(state$working) && any(.detect_spike_features(state$working))
}
# Canonical dose-response source/assay - the same *defaults* the QC Spike-in tab
# uses, so at those defaults colouring/annotating by a spike metric reuses the
# QC tab's `spike_dr` cache entry. (If the QC user picks a non-default source/
# assay there, this path keeps the canonical key and simply forgoes the shared
# entry - aes_helpers stays page-independent and never reads module inputs.)
.aes_spike_assay <- function(state) {
  present <- SummarizedExperiment::assayNames(state$working)
  if ("TPM" %in% present) "TPM" else if ("FPKM" %in% present) "FPKM" else "CPM"
}
.aes_spike_src <- function(state) {
  ids <- rownames(state$working)[.detect_spike_features(state$working)]
  ok <- length(ids) && length(spike_features_missing_conc(state$working)) < length(ids)
  if (ok) "column" else "mix1"
}
# The per-sample dose-response frame via the shared `spike_dr` derived cache.
.aes_spike_metrics <- function(state) {
  src <- .aes_spike_src(state); asy <- .aes_spike_assay(state)
  state_derive(state, "spike_dr", params = list(source = src, assay = asy),
               expr = function() spike_dose_response(state$working, assay = asy, source = src))$per_sample
}

#' Customizable "Other"-domain attribute palette items
#'
#' The session/derived attributes whose palette config lives in the Palette page
#' "Other" domain (`state$palette$other[[item]]`): the removal status, removal
#' pool, and the per-sample QC metrics. The single source of truth shared by the
#' resolver (which reads the configs) and the Palette page (which sets them), so
#' the two stay in sync. (The sample-correlation ramp is a separate app-internal
#' "other" item owned by the Palette page, not an attribute, and not listed here.)
#'
#' @return A named list keyed by palette item id, each
#'   `list(kind, levels, class, label)`.
#' @export
aes_other_palette_items <- function() {
  items <- list(
    removal_status = list(kind = "discrete",
                          levels = c("pass", "suggested_other", "suggested_this"),
                          class = "factor", label = "Suggested removal"),
    `__pool__` = list(kind = "discrete", levels = c("Kept", "In removal pool"),
                      class = "factor", label = "Removal pool"))
  for (m in names(.aes_qc_labels))
    items[[paste0("__qc__", m)]] <- list(kind = "continuous", levels = character(0),
                                         class = "numeric", label = unname(.aes_qc_labels[m]))
  # Spike-in dose-response metrics are customizable too (dataset-independent here;
  # the *selectors* gate them on the dataset having spikes via aes_catalog()).
  for (m in names(.aes_spike_labels))
    items[[paste0("__spike__", m)]] <- list(kind = "continuous", levels = character(0),
                                            class = "numeric", label = unname(.aes_spike_labels[m]))
  items
}

#' Catalog of color/annotation attributes for the current dataset
#'
#' @param state The shared app-state (see [new_app_state()]).
#' @param gene Include the gene-expression attribute (PCA only). Default `FALSE`.
#' @return A list of descriptors, each `list(key, label, group, kind, loc)` where
#'   `group` is one of `"General"`/`"This session"`/`"Data metadata"`, `kind` is
#'   `"discrete"`/`"continuous"`, and `loc` names the palette-config slot
#'   (`list(domain, item)`) or `NULL`.
#' @export
aes_catalog <- function(state, gene = FALSE) {
  out <- list()
  cd <- as.data.frame(SummarizedExperiment::colData(state$working))
  # This session: gene (optional) + QC metrics (continuous) + removal/pool
  # (discrete). `% spike-in` is grouped under "Spike-in" instead (see below).
  if (isTRUE(gene))
    out <- c(out, list(list(key = "__gene__", label = "Gene expression",
                            group = "This session", kind = "continuous", loc = NULL)))
  for (m in setdiff(names(.aes_qc_labels), "pct_spike"))
    out <- c(out, list(list(key = paste0("__qc__", m), label = unname(.aes_qc_labels[m]),
                            group = "This session", kind = "continuous",
                            loc = list("other", paste0("__qc__", m)))))
  out <- c(out, list(
    list(key = "__removal__", label = "Suggested removal", group = "This session",
         kind = "discrete", loc = list("other", "removal_status")),
    list(key = "__pool__", label = "Removal pool", group = "This session",
         kind = "discrete", loc = list("other", "__pool__"))))
  # Spike-in: % spike-in (a per-sample QC metric) + the dose-response metrics,
  # only when the dataset has spikes.
  if (.aes_has_spike(state)) {
    out <- c(out, list(list(key = "__qc__pct_spike", label = unname(.aes_qc_labels["pct_spike"]),
                            group = "Spike-in", kind = "continuous",
                            loc = list("other", "__qc__pct_spike"))))
    for (m in names(.aes_spike_labels))
      out <- c(out, list(list(key = paste0("__spike__", m), label = unname(.aes_spike_labels[m]),
                              group = "Spike-in", kind = "continuous",
                              loc = list("other", paste0("__spike__", m)))))
  }
  # Data metadata: colData columns, typed by class.
  for (col in colnames(cd))
    out <- c(out, list(list(key = col, label = col, group = "Data metadata",
                            kind = if (is.numeric(cd[[col]])) "continuous" else "discrete",
                            loc = list("colData", col))))
  out
}

# Distinct-level count of a colData column (NA -> its own level), for the shape
# cap. Pure helper shared with the catalog filters.
.aes_n_levels <- function(x) length(unique(ifelse(is.na(x), "NA", as.character(x))))

#' Grouped selectInput choices from a catalog
#'
#' @param catalog Output of [aes_catalog()].
#' @param kinds Keep attributes whose `kind` is in this set (default both).
#' @param none Prepend a "General" group with a single `(none)` entry.
#' @param max_levels For discrete attributes, drop those with more than this many
#'   levels (uses the live data; `NULL` = no cap). Used for the shape selector.
#' @param state Needed only when `max_levels` is set (to read level counts).
#' @return A named list of `<optgroup>`s for `selectInput(choices = )`, ordered
#'   General -> This session -> Data metadata (see [group_field_choices()]).
#' @export
aes_choices <- function(catalog, kinds = c("discrete", "continuous"), none = FALSE,
                        max_levels = NULL, state = NULL) {
  keep <- Filter(function(d) d$kind %in% kinds, catalog)
  if (!is.null(max_levels) && !is.null(state)) {
    cd <- as.data.frame(SummarizedExperiment::colData(state$working))
    keep <- Filter(function(d) {
      if (d$kind != "discrete") return(TRUE)
      if (identical(d$group, "Data metadata") && d$key %in% colnames(cd))
        return(.aes_n_levels(cd[[d$key]]) <= max_levels)
      TRUE                                  # session discrete items are small
    }, keep)
  }
  in_group <- function(g) keep[vapply(keep, function(d) identical(d$group, g), logical(1))]
  items_of <- function(ds) stats::setNames(vapply(ds, `[[`, "", "key"),
                                           vapply(ds, `[[`, "", "label"))
  session <- in_group("This session")
  spike   <- in_group("Spike-in")
  coldata <- in_group("Data metadata")
  group_field_choices(vapply(coldata, `[[`, "", "key"), items_of(session), none = none,
                      spike_items = items_of(spike))
}

# --- Value + colour resolution ---------------------------------------------

# Removal-status bundle. `reason` (a QC metric name) enables the reason-aware
# 3-level highlight (General QC); NULL gives the 2-level pass/suggested view.
.aes_removal <- function(state, samples, reason = NULL) {
  fl <- state$samp_flags
  if (is.null(fl)) return(NULL)
  i <- match(samples, fl$sample)
  flagged <- fl$flagged[i]
  if (is.null(reason)) {
    st <- removal_status(flagged); labs <- .removal_labels_2
  } else {
    rcol <- .metric_reason[[reason]]
    this <- if (!is.null(rcol)) fl[[rcol]][i] else NULL
    st <- removal_status(flagged, this); labs <- .removal_labels
  }
  pal <- removal_status_colors(state$palette$other$removal_status)
  present <- levels(droplevels(st))
  values <- factor(unname(labs[as.character(st)]), levels = unname(labs[present]))
  colors <- stats::setNames(unname(pal[present]), unname(labs[present]))
  list(values = values, kind = "discrete", label = "Suggested removal",
       colors = colors, labels = NULL, ramp_config = NULL)
}

.aes_pool <- function(state, samples) {
  pool <- state$samp_pool %||% character(0)
  lv <- c("Kept", "In removal pool")
  cfg <- state$palette$other[["__pool__"]]
  colors <- if (is.null(cfg)) .aes_pool_colors
            else palette_discrete(lv, cfg$colors, cfg$name %||% "Okabe-Ito", cfg$custom)
  list(values = factor(ifelse(samples %in% pool, "In removal pool", "Kept"), levels = lv),
       kind = "discrete", label = "Removal pool", colors = colors,
       labels = NULL, ramp_config = NULL)
}

#' Resolve an attribute to per-sample values + colour mapping
#'
#' @param state The shared app-state.
#' @param key Attribute key (see file header).
#' @param samples Character vector of sample ids to resolve over.
#' @param ctx Optional context: `reason` (QC metric, reason-aware removal); for
#'   `"__gene__"`, `values`/`label` precomputed by the caller and `assay` (the
#'   palette-config assay).
#' @return `list(values, kind, label, colors, labels, ramp_config)`, or `NULL`
#'   when the attribute cannot be resolved yet (e.g. removal before flags exist).
#' @export
aes_resolve <- function(state, key, samples, ctx = list()) {
  if (identical(key, "__removal__")) return(.aes_removal(state, samples, ctx$reason))
  if (identical(key, "__pool__"))    return(.aes_pool(state, samples))
  if (identical(key, "__gene__")) {
    return(list(values = ctx$values, kind = "continuous", label = ctx$label,
                colors = NULL, labels = NULL,
                ramp_config = state$palette$assays[[ctx$assay %||% ""]]))
  }
  if (startsWith(key, "__qc__")) {
    m <- sub("^__qc__", "", key)
    qm <- .aes_qc_metrics(state)
    return(list(values = as.numeric(qm[samples, m]), kind = "continuous",
                label = unname(.aes_qc_labels[m]) %||% m, colors = NULL, labels = NULL,
                ramp_config = state$palette$other[[key]]))
  }
  if (startsWith(key, "__spike__")) {
    if (!.aes_has_spike(state)) return(NULL)
    m <- sub("^__spike__", "", key)
    ps <- .aes_spike_metrics(state)
    return(list(values = as.numeric(ps[match(samples, ps$sample), m]), kind = "continuous",
                label = unname(.aes_spike_labels[m]) %||% m, colors = NULL, labels = NULL,
                ramp_config = state$palette$other[[key]]))
  }
  # colData column.
  cd <- as.data.frame(SummarizedExperiment::colData(state$working))
  if (!key %in% colnames(cd)) return(NULL)
  x <- cd[samples, key]
  if (is.numeric(x))
    return(list(values = x, kind = "continuous", label = key, colors = NULL,
                labels = NULL, ramp_config = state$palette$colData[[key]]))
  lv <- as.character(x); lv[is.na(lv)] <- "NA"; f <- factor(lv)
  cfg <- state$palette$colData[[key]]
  colors <- if (is.null(cfg)) NULL
            else palette_discrete(levels(f), cfg$colors, cfg$name %||% "Okabe-Ito", cfg$custom)
  list(values = f, kind = "discrete", label = key, colors = colors,
       labels = NULL, ramp_config = NULL)
}

#' Build a ggplot scale from a resolved attribute
#'
#' @param res An [aes_resolve()] result (non-NULL).
#' @param aes_name `"colour"` (default) or `"fill"`.
#' @return A ggplot scale to add with `+`, or `NULL` to keep thematic's default
#'   (no project palette configured).
#' @export
aes_ggplot_scale <- function(res, aes_name = "colour") {
  manual <- if (aes_name == "fill") ggplot2::scale_fill_manual else ggplot2::scale_colour_manual
  gradn  <- if (aes_name == "fill") ggplot2::scale_fill_gradientn else ggplot2::scale_colour_gradientn
  if (identical(res$kind, "discrete")) {
    if (is.null(res$colors)) return(NULL)
    labs <- if (is.null(res$labels)) ggplot2::waiver() else res$labels
    return(manual(values = res$colors, labels = labs, na.value = "grey70"))
  }
  cfg <- res$ramp_config
  if (is.null(cfg)) return(NULL)                             # no config -> thematic default
  g <- palette_gradientn(cfg$name %||% "viridis: viridis", values = res$values,
                         min = cfg$min, max = cfg$max, custom = cfg$custom,
                         reverse = isTRUE(cfg$reverse))
  gradn(colours = g$colours, values = g$values, limits = g$limits)
}

#' Build a ComplexHeatmap annotation colour entry from a resolved attribute
#'
#' @param res An [aes_resolve()] result (non-NULL).
#' @return A named colour vector (discrete) or a [circlize::colorRamp2()]
#'   function (continuous; `NULL` if circlize is unavailable).
#' @export
aes_heatmap_col <- function(res) {
  if (identical(res$kind, "discrete"))
    return(res$colors %||% palette_discrete(levels(res$values), NULL, "Okabe-Ito"))
  cfg <- res$ramp_config
  if (!is.null(cfg) && !is.null(cfg$name))
    return(palette_colorramp2(cfg$name, values = res$values, min = cfg$min, max = cfg$max,
                              custom = cfg$custom, reverse = isTRUE(cfg$reverse)))
  if (!requireNamespace("circlize", quietly = TRUE)) return(NULL)
  rng <- range(res$values, na.rm = TRUE)                    # default viridis-like ramp
  if (!is.finite(diff(rng)) || diff(rng) == 0) rng <- c(rng[1] - 0.5, rng[1] + 0.5)
  circlize::colorRamp2(seq(rng[1], rng[2], length.out = 3),
                       c("#440154", "#21908C", "#FDE725"))
}

#' Build a HeatmapAnnotation data.frame + colour list for several attributes
#'
#' @param state The shared app-state.
#' @param keys Character vector of attribute keys (annotation tracks).
#' @param samples Sample ids (rows of the annotation frame).
#' @return `list(df, col)` for [ComplexHeatmap::HeatmapAnnotation()], or `NULL`
#'   when nothing resolves. Tracks that cannot resolve yet are skipped.
#' @export
aes_annotation <- function(state, keys, samples) {
  built <- Filter(Negate(is.null), lapply(keys, function(k) {
    r <- aes_resolve(state, k, samples)
    if (is.null(r)) NULL else list(name = r$label, res = r)
  }))
  if (!length(built)) return(NULL)
  names_ <- vapply(built, `[[`, "", "name")
  df <- data.frame(stats::setNames(lapply(built, function(b) b$res$values), names_),
                   row.names = samples, check.names = FALSE, stringsAsFactors = FALSE)
  col <- stats::setNames(lapply(built, function(b) aes_heatmap_col(b$res)), names_)
  list(df = df, col = col)
}
