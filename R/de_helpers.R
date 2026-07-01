# Phase 5: differential expression. Pure, exported helpers the mod_de page calls
# -- validating a design (full rank), releveling the reference, enumerating
# contrast levels, running DESeq2 + shrinking LFCs into the dual-LFC schema
# (log2FoldChange(_shrunk) + sig/DEG(_shrunk)), and building the MA / volcano /
# direct-comparison plots (axis clamping -> triangle markers). The statistical
# conventions live in the rnaseq-bioc / de-analysis skills; this file implements
# them. Shiny-free and unit-tested; the module wires + caches (state_derive).

# --- classification (the dual-LFC sig/DEG schema) --------------------------

#' Classify features into significance + DEG direction
#'
#' The rnaseq-bioc rule: significant when `padj` is present and below `padj_cut`
#' **and** the (chosen) `abs(lfc)` is at least `lfc_cut`. `DEG` is a factor with
#' fixed levels `up`/`down`/`no_change`.
#'
#' @param lfc Numeric log2 fold-changes (standard or shrunk).
#' @param padj Numeric adjusted p-values (shared across LFC variants).
#' @param padj_cut Adjusted-p cutoff (default 0.05).
#' @param lfc_cut Absolute-LFC cutoff (default `log2(2)`).
#' @return A list with `sig` (logical) and `DEG` (factor `up`/`down`/`no_change`).
#' @export
classify <- function(lfc, padj, padj_cut = 0.05, lfc_cut = log2(2)) {
  sig <- !is.na(padj) & padj < padj_cut & !is.na(lfc) & abs(lfc) >= lfc_cut
  DEG <- factor(ifelse(!sig, "no_change", ifelse(lfc > 0, "up", "down")),
                levels = c("up", "down", "no_change"))
  list(sig = sig, DEG = DEG)
}

#' Add sig/DEG(_shrunk) columns to a results frame from thresholds
#'
#' Cheap post-processing (no refit) so the padj/LFC thresholds can drive the
#' plots/table live. Adds the shrunk pair only when `log2FoldChange_shrunk` is
#' present.
#' @param df A results `data.frame` (from [de_results()]).
#' @param padj_cut,lfc_cut Thresholds (see [classify()]).
#' @return `df` with `sig`/`DEG` (+ `sig_shrunk`/`DEG_shrunk`) columns.
#' @export
de_classify_table <- function(df, padj_cut = 0.05, lfc_cut = log2(2)) {
  std <- classify(df$log2FoldChange, df$padj, padj_cut, lfc_cut)
  df$sig <- std$sig
  df$DEG <- std$DEG
  if ("log2FoldChange_shrunk" %in% names(df)) {
    shk <- classify(df$log2FoldChange_shrunk, df$padj, padj_cut, lfc_cut)
    df$sig_shrunk <- shk$sig
    df$DEG_shrunk <- shk$DEG
  }
  df
}

#' Count DEGs by direction
#' @param df A classified results frame.
#' @param deg_col Which DEG column (`"DEG"` or `"DEG_shrunk"`).
#' @return A named integer vector `up`/`down`/`total`.
#' @export
de_summary <- function(df, deg_col = "DEG") {
  d <- df[[deg_col]]
  c(up    = sum(d == "up", na.rm = TRUE),
    down  = sum(d == "down", na.rm = TRUE),
    total = sum(d %in% c("up", "down"), na.rm = TRUE))
}

# --- design + contrast building --------------------------------------------

#' Is a design full rank on this sample table?
#'
#' Builds the model matrix and compares its numeric rank to the number of
#' coefficients; a deficient rank means confounded/collinear terms that DESeq2
#' will reject. Unused factor levels are dropped first (matching a filtered dds).
#' @param design A one-sided model formula (e.g. `~ condition`).
#' @param col_data A `DataFrame`/`data.frame` of sample metadata.
#' @return A list: `ok` (logical), `rank`, `ncoef`, and `msg` (`NULL` when ok).
#' @export
de_full_rank <- function(design, col_data) {
  cd <- droplevels(as.data.frame(col_data))
  mm <- tryCatch(stats::model.matrix(design, data = cd), error = function(e) e)
  if (inherits(mm, "error")) {
    return(list(ok = FALSE, rank = NA_integer_, ncoef = NA_integer_,
                msg = paste0("Could not build the model matrix: ", conditionMessage(mm))))
  }
  rank <- qr(mm)$rank
  ncoef <- ncol(mm)
  ok <- isTRUE(rank == ncoef)
  list(ok = ok, rank = rank, ncoef = ncoef,
       msg = if (ok) NULL else sprintf(paste0(
         "Design is not full rank (rank %d < %d model coefficients) - some terms ",
         "are collinear or confounded. Drop a covariate or add replication."),
         rank, ncoef))
}

#' Candidate design factors: discrete colData columns with >= 2 levels
#' @param dds A `DESeqDataSet`.
#' @return Character vector of column names.
#' @export
de_design_factors <- function(dds) {
  cd <- SummarizedExperiment::colData(dds)
  keep <- vapply(names(cd), function(nm) {
    x <- cd[[nm]]
    (is.factor(x) || is.character(x) || is.logical(x)) &&
      length(unique(x[!is.na(x)])) >= 2L
  }, logical(1))
  names(cd)[keep]
}

#' Levels of a design variable (for contrast test/control pickers)
#' @param dds A `DESeqDataSet`.
#' @param var A colData column name.
#' @return Character vector of levels (factor order preserved; else sorted uniques).
#' @export
de_contrast_levels <- function(dds, var) {
  x <- SummarizedExperiment::colData(dds)[[var]]
  if (is.null(x)) return(character(0))
  if (is.factor(x)) levels(droplevels(x)) else sort(unique(as.character(x[!is.na(x)])))
}

#' Relevel a colData factor so `ref` is the reference (first) level
#'
#' The DE primitive for setting the control level; Phase 8 factor management
#' generalizes reordering and reuses this.
#' @param dds A `DESeqDataSet`.
#' @param col A colData column name.
#' @param ref The level to make the reference.
#' @return `dds` with the column releveled.
#' @export
de_relevel <- function(dds, col, ref) {
  cd <- SummarizedExperiment::colData(dds)
  x <- cd[[col]]
  if (is.null(x)) stop("Column '", col, "' is not in colData.", call. = FALSE)
  if (!is.factor(x)) x <- factor(x)
  if (!ref %in% levels(x)) stop("Level '", ref, "' is not in '", col, "'.", call. = FALSE)
  cd[[col]] <- stats::relevel(x, ref = ref)
  SummarizedExperiment::colData(dds) <- cd
  dds
}

#' The coefficient name DESeq2 assigns to a two-level factor contrast
#'
#' Matches `DESeq2::resultsNames()` when `control` is the reference level; used to
#' route `apeglm` shrinkage (which needs a coefficient, not a contrast).
#' @param contrast `c(variable, test, control)`.
#' @return The coefficient name, e.g. `"condition_treated_vs_control"`.
#' @export
de_coef_name <- function(contrast) {
  paste0(contrast[1], "_", contrast[2], "_vs_", contrast[3])
}

# --- fitting + results ------------------------------------------------------

# Drop unused levels of the design factors before fitting. Critical: a metadata
# relabel can leave an EMPTY factor level (an all-zero model-matrix column ->
# rank-deficient), so `DESeq()` fails "full model matrix is less than full rank"
# even though the design-builder badge (which droplevels) reports full rank. This
# keeps the fit and the badge in agreement; also rejects NA-in-design early with a
# clear message rather than DESeq2's opaque one.
.de_droplevels_design <- function(dds) {
  vars <- tryCatch(all.vars(DESeq2::design(dds)), error = function(e) character(0))
  cd <- SummarizedExperiment::colData(dds)
  for (v in intersect(vars, colnames(cd))) {
    x <- cd[[v]]
    if (is.factor(x) || is.character(x)) {
      if (anyNA(x)) {
        stop("Design factor '", v, "' has missing values; fix the metadata before running DE.",
             call. = FALSE)
      }
      cd[[v]] <- droplevels(if (is.factor(x)) x else factor(x))
    }
  }
  SummarizedExperiment::colData(dds) <- cd
  dds
}

#' Run DESeq2 with a robust size-factor step
#'
#' Drops unused design-factor levels (so the fit agrees with [de_full_rank()] and
#' an empty level doesn't trip "full model matrix is less than full rank"), then
#' ensures size factors exist (via [estimate_size_factors_endogenous()], which
#' falls back to the positive-counts estimator on sparse data) before fitting.
#' @param dds A `DESeqDataSet` (design already set).
#' @param quiet Passed to `DESeq2::DESeq()`.
#' @return The fitted `DESeqDataSet`.
#' @export
de_run <- function(dds, quiet = TRUE) {
  dds <- .de_droplevels_design(dds)
  sf <- tryCatch(DESeq2::sizeFactors(dds), error = function(e) NULL)
  if (is.null(sf)) dds <- estimate_size_factors_endogenous(dds)
  DESeq2::DESeq(dds, quiet = quiet)
}

#' Validity of a stored contrast against the current dds
#'
#' Distinguishes three states so the UI can label + decide extraction:
#' `"extractable"` (factor is in the design and both levels exist -> can produce
#' results now), `"not_in_design"` (factor + levels exist but the factor isn't a
#' current design term -> kept; add it to the design to extract), and `"invalid"`
#' (the factor column or a level no longer exists in the data -> kept unless
#' removed; recoverable if restored).
#' @param dds A `DESeqDataSet`.
#' @param spec A contrast spec `list(var, test, control, ...)`.
#' @return One of `"extractable"`, `"not_in_design"`, `"invalid"`.
#' @export
de_contrast_validity <- function(dds, spec) {
  cd <- SummarizedExperiment::colData(dds)
  if (is.null(spec$var) || is.null(cd[[spec$var]])) return("invalid")
  lv <- de_contrast_levels(dds, spec$var)
  if (!all(c(spec$test, spec$control) %in% lv)) return("invalid")
  design_vars <- tryCatch(all.vars(DESeq2::design(dds)), error = function(e) character(0))
  if (!(spec$var %in% design_vars)) return("not_in_design")
  "extractable"
}

#' Shrunken log2 fold-changes for a contrast (apeglm -> ashr -> normal fallback)
#'
#' `apeglm` needs a coefficient (available when `control` is the reference level);
#' `ashr` takes the contrast directly; `type="normal"` is the built-in last
#' resort. All optional packages are probed with `requireNamespace()` and every
#' call is wrapped so a failure degrades to `NA` rather than erroring.
#' @param fit A fitted `DESeqDataSet` (from [de_run()]).
#' @param contrast `c(variable, test, control)`.
#' @param type `"apeglm"` (default), `"ashr"`, or `"none"`.
#' @return A numeric vector aligned to `rownames(fit)` (`NA` where unavailable),
#'   carrying an `attr(., "method")` naming the estimator that actually ran
#'   (`"apeglm"`/`"ashr"`/`"normal"`/`"none"`) so the UI can label honestly.
#' @export
de_shrink <- function(fit, contrast, type = c("apeglm", "ashr", "none")) {
  type <- match.arg(type)
  ids <- rownames(fit)
  na_vec <- rep(NA_real_, length(ids))
  if (identical(type, "none")) return(structure(na_vec, method = "none"))

  coef_name <- de_coef_name(contrast)
  has_coef <- coef_name %in% DESeq2::resultsNames(fit)
  try_shrink <- function(...) tryCatch(DESeq2::lfcShrink(fit, ..., quiet = TRUE),
                                       error = function(e) NULL)
  sh <- NULL; used <- "none"
  if (identical(type, "apeglm") && has_coef && requireNamespace("apeglm", quietly = TRUE)) {
    sh <- try_shrink(coef = coef_name, type = "apeglm"); if (!is.null(sh)) used <- "apeglm"
  }
  if (is.null(sh) && requireNamespace("ashr", quietly = TRUE)) {
    sh <- try_shrink(contrast = contrast, type = "ashr"); if (!is.null(sh)) used <- "ashr"
  }
  if (is.null(sh) && has_coef) {
    sh <- try_shrink(coef = coef_name, type = "normal"); if (!is.null(sh)) used <- "normal"
  }
  if (is.null(sh)) return(structure(na_vec, method = "none"))
  structure(as.numeric(sh$log2FoldChange)[match(ids, rownames(sh))], method = used)
}

#' A results table carrying both LFC variants
#'
#' `DESeq2::results()` for the standard `log2FoldChange` + shared `padj`, plus a
#' `log2FoldChange_shrunk` column from [de_shrink()]. Classification (`sig`/`DEG`)
#' is added separately by [de_classify_table()] so thresholds stay cheap to tweak.
#' @param fit A fitted `DESeqDataSet`.
#' @param contrast `c(variable, test, control)`.
#' @param shrink_type `"apeglm"` (default), `"ashr"`, or `"none"`.
#' @return A `data.frame` (rownames = feature ids) with the DESeq2 result columns
#'   plus `log2FoldChange_shrunk`.
#' @export
de_results <- function(fit, contrast, shrink_type = c("apeglm", "ashr", "none")) {
  shrink_type <- match.arg(shrink_type)
  df <- as.data.frame(DESeq2::results(fit, contrast = contrast))
  sh <- de_shrink(fit, contrast, shrink_type)
  df$log2FoldChange_shrunk <- as.numeric(sh)
  attr(df, "shrink_method") <- attr(sh, "method")   # what actually ran (for the UI label)
  df
}

#' Per-group mean expression for the direct-comparison plot
#'
#' Averages the chosen `assay` **on its native scale** for each group. Pass a
#' depth-normalized / log assay (e.g. `logcounts`, the default) — averaging raw
#' `counts` is depth-confounded and not comparable across the axes. Both groups
#' must be non-empty and their ids present in the dds.
#' @param dds A `DESeqDataSet`.
#' @param assay Assay name to average (default `"logcounts"`).
#' @param control_samples,test_samples Sample-id vectors for the two groups.
#' @return A `data.frame(id, control, test)`.
#' @export
de_group_means <- function(dds, assay = "logcounts", control_samples, test_samples) {
  m <- as.matrix(SummarizedExperiment::assay(dds, assay))
  if (!length(control_samples) || !length(test_samples)) {
    stop("Both groups need at least one sample.", call. = FALSE)
  }
  miss <- setdiff(c(control_samples, test_samples), colnames(m))
  if (length(miss)) {
    stop("Unknown sample(s): ", paste(miss, collapse = ", "), call. = FALSE)
  }
  data.frame(
    id      = rownames(m),
    control = rowMeans(m[, control_samples, drop = FALSE], na.rm = TRUE),
    test    = rowMeans(m[, test_samples, drop = FALSE], na.rm = TRUE),
    row.names = NULL, stringsAsFactors = FALSE)
}

# --- plot building (MA / volcano / direct comparison) ----------------------

#' Clamp values to a range, flagging what was moved
#'
#' Out-of-range points are pulled to the nearest limit (and flagged) so they draw
#' as triangles rather than being dropped.
#' @param x Numeric vector.
#' @param lo,hi Optional finite limits (`NULL`/`NA` = no clamp on that side).
#' @return A list `value` (clamped) + `clamped` (logical).
#' @export
de_clamp <- function(x, lo = NULL, hi = NULL) {
  v <- x
  clamped <- rep(FALSE, length(x))
  if (!is.null(lo) && is.finite(lo)) {
    below <- !is.na(v) & v < lo; v[below] <- lo; clamped <- clamped | below
  }
  if (!is.null(hi) && is.finite(hi)) {
    above <- !is.na(v) & v > hi; v[above] <- hi; clamped <- clamped | above
  }
  list(value = v, clamped = clamped)
}

# Default DEG colours (overridden by the Palette other/DEG config at the call site).
.de_deg_colors <- c(up = "#D62728", down = "#1F77B4", no_change = "#B0B0B0")

#' Resolve a per-feature colour aesthetic for a DE plot
#'
#' Local to the DE module (the shared `aes_helpers` resolver is per-sample). DEG
#' columns get a discrete manual scale (colours from the Palette `other/DEG`
#' config when supplied); numeric DE columns get a continuous scale; `"__none__"`
#' means no colour mapping.
#' @param df A classified results frame.
#' @param key `"__none__"`, `"DEG"`/`"DEG_shrunk"`, or a numeric column name.
#' @param colors Optional named `up`/`down`/`no_change` colours (Palette config).
#' @return `NULL`, or a list `values` / `kind` / `label` / `scale`.
#' @export
de_colour_resolve <- function(df, key, colors = NULL) {
  if (is.null(key) || identical(key, "__none__")) return(NULL)
  if (key %in% c("DEG", "DEG_shrunk")) {
    cols <- colors %||% .de_deg_colors
    list(values = df[[key]], kind = "discrete", label = "DEG",
         scale = ggplot2::scale_colour_manual(values = cols, drop = FALSE,
                                               na.value = "grey80"))
  } else {
    list(values = suppressWarnings(as.numeric(df[[key]])), kind = "continuous",
         label = key, scale = ggplot2::scale_colour_viridis_c())
  }
}

.de_require <- function(df, cols) {
  miss <- setdiff(cols, names(df))
  if (length(miss)) stop("Missing column(s): ", paste(miss, collapse = ", "), call. = FALSE)
  invisible(TRUE)
}

# Shared scatter: clamp -> shape(triangle) for out-of-range, optional colour, an
# optional ggrepel label layer, and a hover `text` aes only when interactive.
.de_scatter <- function(d, xlab, ylab, colour = NULL, x_range = NULL, y_range = NULL,
                        labels = NULL, point_size = 1.4, interactive = FALSE) {
  cx <- de_clamp(d$x, x_range[1], x_range[2])
  cy <- de_clamp(d$y, y_range[1], y_range[2])
  d$x <- cx$value
  d$y <- cy$value
  d$clamped <- factor(ifelse(cx$clamped | cy$clamped, "clamped", "in range"),
                      levels = c("in range", "clamped"))
  has_col <- !is.null(colour)
  if (has_col) d$colour <- colour$values
  if (interactive) {
    d$text <- sprintf("%s\n%s: %.3g\n%s: %.3g", d$id, xlab, d$x, ylab, d$y)
  }

  base_map <- if (has_col) {
    ggplot2::aes(x = .data$x, y = .data$y, shape = .data$clamped, colour = .data$colour)
  } else {
    ggplot2::aes(x = .data$x, y = .data$y, shape = .data$clamped)
  }
  hov <- if (interactive) ggplot2::aes(text = .data$text) else NULL

  p <- ggplot2::ggplot(d, base_map) +
    ggplot2::geom_point(mapping = hov, size = point_size, alpha = 0.85) +
    ggplot2::scale_shape_manual(values = c("in range" = 16, "clamped" = 17),
                                drop = FALSE, guide = "none") +
    ggplot2::labs(x = xlab, y = ylab, colour = if (has_col) colour$label else NULL)
  if (has_col && !is.null(colour$scale)) p <- p + colour$scale
  if (!is.null(labels) && nrow(labels)) {
    geom_lab <- if (requireNamespace("ggrepel", quietly = TRUE)) {
      ggrepel::geom_text_repel
    } else {
      ggplot2::geom_text
    }
    p <- p + geom_lab(data = labels,
                      mapping = ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
                      inherit.aes = FALSE, size = 3, max.overlaps = Inf)
  }
  p
}

#' MA plot: x = log10(baseMean), y = (chosen) log2FoldChange
#'
#' @param df A classified results frame.
#' @param lfc_col Which LFC column (`"log2FoldChange"` or `"log2FoldChange_shrunk"`).
#' @param colour A [de_colour_resolve()] result, or `NULL`.
#' @param x_range,y_range Optional `c(lo, hi)` axis clamps (triangle markers).
#' @param labels Optional `data.frame(x, y, label)` for a ggrepel layer.
#' @param point_size Point size.
#' @param interactive Add a hover `text` aes (for the plotly path).
#' @return A ggplot object.
#' @export
de_ma_gg <- function(df, lfc_col = "log2FoldChange", colour = NULL,
                     x_range = NULL, y_range = NULL, labels = NULL,
                     point_size = 1.4, interactive = FALSE) {
  .de_require(df, c("baseMean", lfc_col))
  d <- data.frame(id = rownames(df), x = log10(df$baseMean), y = df[[lfc_col]])
  .de_scatter(d, "log10(baseMean)", lfc_col, colour, x_range, y_range,
              labels, point_size, interactive)
}

#' Volcano plot: x = (chosen) log2FoldChange, y = -log10(padj)
#' @inheritParams de_ma_gg
#' @return A ggplot object.
#' @export
de_volcano_gg <- function(df, lfc_col = "log2FoldChange", colour = NULL,
                          x_range = NULL, y_range = NULL, labels = NULL,
                          point_size = 1.4, interactive = FALSE) {
  .de_require(df, c("padj", lfc_col))
  d <- data.frame(id = rownames(df), x = df[[lfc_col]], y = -log10(df$padj))
  .de_scatter(d, lfc_col, "-log10(padj)", colour, x_range, y_range,
              labels, point_size, interactive)
}

#' Direct-comparison plot: x = control mean, y = test mean expression
#' @param mean_df A `data.frame(id, control, test)` (from [de_group_means()]).
#' @param value_label Scale/assay shown in the axis labels (e.g. `"logcounts"`),
#'   so the y=x guide isn't read as a linear comparison when the assay is log-scale.
#' @param colour,x_range,y_range,labels,point_size,interactive As in [de_ma_gg()].
#' @return A ggplot object.
#' @export
de_direct_gg <- function(mean_df, value_label = "mean", colour = NULL,
                         x_range = NULL, y_range = NULL,
                         labels = NULL, point_size = 1.4, interactive = FALSE) {
  .de_require(mean_df, c("control", "test"))
  d <- data.frame(id = mean_df$id, x = mean_df$control, y = mean_df$test)
  .de_scatter(d, sprintf("control (%s)", value_label), sprintf("test (%s)", value_label),
              colour, x_range, y_range, labels, point_size, interactive) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                         colour = "grey60", linewidth = 0.3)
}
