# Phase 7: Expression page. Pure, exported helpers the mod_expression page calls --
# choosing/resolving the value matrix (incl. VST and normalized log-counts), the
# per-gene z-score used by the gene-set heatmap, the single-gene long frame, and
# the sample-count guards that decide which geoms (violin/box/dots) are offered.
# Unlike the PCA helpers these do NOT restrict to endogenous rows (a gene set or a
# searched single gene may be a spike-in/exogenous feature) -- except VST, which is
# inherently computed on the endogenous diagnostic rows via qc_vst().

#' Default expression value key for a dataset
#'
#' Picks the most interpretable value for the Expression page: size-factor
#' normalized log-counts when size factors exist, otherwise the variance-
#' stabilized transform. Both are computed value keys (resolved by
#' [expr_value_matrix()]); VST is always attemptable with a log-counts fallback
#' at resolve time, so no stored-assay fallthrough is needed here.
#'
#' @param dds A `DESeqDataSet`.
#' @return A single value key: `"norm_logcounts"` or `"vst"`.
#' @export
expr_default_assay <- function(dds) {
  if (!is.null(DESeq2::sizeFactors(dds))) "norm_logcounts" else "vst"
}

#' Resolve an Expression value matrix (genes x samples) + an honest label
#'
#' Handles the two computed inputs (`"vst"`, `"norm_logcounts"`) and any stored
#' assay. VST is computed on the endogenous diagnostic rows (via [qc_vst()]) and
#' falls back to log-counts on any failure; `"norm_logcounts"` estimates
#' endogenous size factors when absent. Stored assays are returned as-is (the
#' caller's transform/pseudocount then applies on top).
#'
#' @param dds A `DESeqDataSet`.
#' @param assay A value key (`"vst"`, `"norm_logcounts"`, or a stored assay name).
#' @return A list with `mat` (numeric genes x samples) and `label` (for the axis /
#'   subtitle, reflecting any fallback).
#' @export
expr_value_matrix <- function(dds, assay = "norm_logcounts") {
  fallback_logcounts <- function(why)
    list(mat = .qc_assay_matrix(dds, "logcounts"),
         label = paste0("logcounts (", why, ")"))

  if (identical(assay, "vst")) {
    res <- tryCatch(
      list(mat = as.matrix(SummarizedExperiment::assay(qc_vst(dds, blind = TRUE))),
           label = "VST (blind)"),
      error = function(e) NULL)
    return(res %||% fallback_logcounts("VST unavailable"))
  }

  if (identical(assay, "norm_logcounts")) {
    res <- tryCatch({
      d <- if (is.null(DESeq2::sizeFactors(dds))) estimate_size_factors_endogenous(dds) else dds
      nc <- as.matrix(DESeq2::counts(d, normalized = TRUE))
      list(mat = log2(nc + 1), label = "normalized log-counts (log2)")
    }, error = function(e) NULL)
    return(res %||% fallback_logcounts("normalized counts unavailable"))
  }

  list(mat = .qc_assay_matrix(dds, assay), label = assay)
}

#' Per-gene (row-wise) z-score of a value matrix
#'
#' Centers and scales each row to mean 0 / sd 1 -- the heatmap default so genes on
#' very different expression scales are comparable. Zero-variance rows (constant
#' expression) would divide by zero; they are returned as all-zero rows rather than
#' `NaN`. Non-finite inputs are treated as the row mean before scaling.
#'
#' @param mat A numeric genes x samples matrix.
#' @return A matrix the same shape, each row z-scored (constant rows -> 0).
#' @export
row_zscore <- function(mat) {
  mat <- as.matrix(mat)
  if (!nrow(mat) || !ncol(mat)) return(mat)
  mu <- rowMeans(mat, na.rm = TRUE)
  centered <- mat - mu
  n_row <- rowSums(!is.na(mat))                          # per-row non-NA count
  sdv <- sqrt(rowSums(centered^2, na.rm = TRUE) / pmax(n_row - 1L, 1L))
  sdv[!is.finite(sdv) | sdv == 0] <- 1  # constant row -> centered values are 0
  z <- centered / sdv
  z[!is.finite(z)] <- 0
  dimnames(z) <- dimnames(mat)
  z
}

#' Which distribution geoms the single-gene plot should offer
#'
#' Decides whether data points are allowed / on by default and whether the
#' violin/box distributions are shown, from the per-group sample counts. The two
#' concerns use separate thresholds: the distribution geoms appear once any group
#' is large enough to summarize (`dist_min`), while data points default on up to a
#' larger group size (`dots_max`) and are disallowed only for genuinely overplotted
#' groups (`dots_hard`).
#'
#' @param group_sizes Integer vector of per-group sample counts (>= 1 each).
#' @param dist_min Minimum group size for the violin/box geoms to be offered
#'   (default `getOption("ddsdashboard.expr_dist_min", 10)`).
#' @param dots_max Data points default ON when the largest group is below this
#'   (default `getOption("ddsdashboard.expr_dots_max", 100)`).
#' @param dots_hard Hard cap above which data points are never drawn
#'   (default `getOption("ddsdashboard.expr_dots_hard", 500)`).
#' @return A list: `n_max` (largest group), `dots_allowed`, `dots_default`,
#'   `dist_shown` (whether violin/box are offered).
#' @export
expr_geom_availability <- function(group_sizes,
                                   dist_min = getOption("ddsdashboard.expr_dist_min", 10L),
                                   dots_max = getOption("ddsdashboard.expr_dots_max", 100L),
                                   dots_hard = getOption("ddsdashboard.expr_dots_hard", 500L)) {
  sizes <- as.integer(group_sizes[is.finite(group_sizes)])
  n_max <- if (length(sizes)) max(sizes) else 0L
  dots_allowed <- n_max < dots_hard
  list(
    n_max = n_max,
    dots_allowed = dots_allowed,
    dots_default = dots_allowed && n_max < dots_max,
    dist_shown = length(sizes) > 0 && any(sizes >= dist_min)
  )
}

#' Assemble the single-gene long data frame
#'
#' Joins one gene's per-sample values to a grouping vector (and optional colour
#' vector), dropping samples with a missing group. The grouping keeps its factor
#' level order when supplied as a factor.
#'
#' @param values Named (or same-order) numeric vector of the gene's values.
#' @param groups Per-sample grouping (character/factor), same length/order as
#'   `values`.
#' @param samples Sample ids (same length/order); defaults to `names(values)`.
#' @param colour Optional per-sample colour attribute (same length/order).
#' @return A data.frame with `sample`, `group` (factor), `value`, and (if given)
#'   `colour`.
#' @export
expr_long_frame <- function(values, groups, samples = names(values), colour = NULL) {
  n <- length(values)
  if (length(groups) != n) stop("groups must match values in length", call. = FALSE)
  if (is.null(samples)) samples <- as.character(seq_len(n))
  df <- data.frame(
    sample = as.character(samples),
    group  = if (is.factor(groups)) groups else factor(groups),
    value  = as.numeric(values),
    stringsAsFactors = FALSE)
  if (!is.null(colour)) df$colour <- colour
  df <- df[!is.na(df$group), , drop = FALSE]
  df$group <- droplevels(df$group)                       # drop groups with no shown samples
  df
}

#' Aggregate a gene set's expression into one value per sample
#'
#' Collapses the value-matrix rows of a gene set to a single per-sample score
#' (mean or median) -- the basis for the Gene Sets "Aggregate expression" plot.
#' Genes absent from the dataset are dropped; with `only_expressed = TRUE` genes
#' with zero counts in every sample are dropped too (needs `counts`). The kept
#' rows are transformed (`transform`/`pseudocount`), optionally per-gene z-scored
#' across samples (the default, so highly-expressed genes do not dominate the
#' average), then averaged down to one vector. All the accounting is returned so
#' the caller can caption "<n> of <N> genes" and warn about z-scored constant
#' genes.
#'
#' @param mat A numeric genes x samples value matrix (e.g. from
#'   [expr_value_matrix()]).
#' @param ids Authored gene-set ids (matched against `rownames(mat)`).
#' @param counts Optional raw counts matrix for the `only_expressed` filter;
#'   `NULL` skips that filter (all present genes are used).
#' @param method `"mean"` (default) or `"median"` across genes, per sample.
#' @param zscore Per-gene z-score across samples before averaging (default TRUE).
#' @param only_expressed Drop genes with all-zero counts (default TRUE).
#' @param transform Passed to [expr_transform()] (`"none"`/`"log2"`/`"log10"`),
#'   applied to the kept rows before z-scoring / averaging.
#' @param pseudocount Added before a log transform.
#' @return A list: `values` (named per-sample numeric, or `NULL` when no gene
#'   survives), `n_total` (authored ids), `n_present` (in the dataset),
#'   `n_absent`, `n_used` (genes averaged), `n_nonvar` (kept genes constant
#'   across samples -> z-score 0), `ids_used`, `method`.
#' @export
expr_set_aggregate <- function(mat, ids, counts = NULL,
                               method = c("mean", "median"),
                               zscore = TRUE, only_expressed = TRUE,
                               transform = "none", pseudocount = 1) {
  method <- match.arg(method)
  mat <- as.matrix(mat)
  sel <- .expr_set_rows(rownames(mat), ids, counts = counts,
                        only_expressed = only_expressed)
  n_total <- sel$n_total; n_present <- length(sel$present)
  used <- sel$used; n_used <- length(used)

  out <- list(values = NULL, n_total = n_total, n_present = n_present,
              n_absent = n_total - n_present, n_used = n_used,
              n_nonvar = 0L, ids_used = used, method = method)
  if (!n_used) return(out)

  sub <- expr_transform(mat[used, , drop = FALSE], transform, pseudocount)
  # Genes with no spread across samples become 0 under z-scoring; count them so
  # the caller can warn (a constant gene contributes nothing to a z-scored score).
  rng <- apply(sub, 1, function(r) { r <- r[is.finite(r)]; if (!length(r)) 0 else diff(range(r)) })
  out$n_nonvar <- sum(rng == 0)
  if (isTRUE(zscore)) sub <- row_zscore(sub)
  agg <- if (identical(method, "median")) apply(sub, 2, stats::median, na.rm = TRUE)
         else colMeans(sub, na.rm = TRUE)
  names(agg) <- colnames(mat)
  out$values <- agg
  out
}

# Shared gene-set row selection (used by expr_set_aggregate + expr_heatmap_matrix):
# authored ids -> those present in the matrix -> optionally drop genes with all-zero
# counts. Order-preserving. `only_expressed` needs `counts`; genes absent from
# `counts` can't be judged and are kept (conservative). Returns the deduped `ids`,
# `n_total`, the `present` ids, and the `used` ids (after the expression filter).
.expr_set_rows <- function(mat_rownames, ids, counts = NULL, only_expressed = TRUE) {
  ids <- unique(as.character(ids)); ids <- ids[!is.na(ids)]
  present <- ids[ids %in% mat_rownames]
  used <- present
  if (isTRUE(only_expressed) && length(present) && !is.null(counts)) {
    counts <- as.matrix(counts)
    in_counts <- present[present %in% rownames(counts)]
    expressed <- in_counts[rowSums(counts[in_counts, , drop = FALSE], na.rm = TRUE) > 0]
    keep <- union(expressed, setdiff(present, rownames(counts)))
    used <- present[present %in% keep]                   # preserve the set's order
  }
  list(ids = ids, n_total = length(ids), present = present, used = used)
}

#' Prepare a gene-set expression matrix for the heatmap
#'
#' The heatmap counterpart of [expr_set_aggregate()]: it selects the same gene-set
#' rows (present in the dataset, optionally dropping all-zero-count genes),
#' transforms them, and -- by default -- per-gene z-scores across samples so genes
#' on very different expression scales are comparable. Unlike the aggregate helper
#' it keeps the full genes x samples matrix (no averaging). Constant (zero-variance)
#' rows z-score to 0 (via [row_zscore()]) rather than `NaN`, so they still draw as a
#' flat mid-colour row instead of breaking the heatmap; their count is returned so
#' the caller can note it.
#'
#' @param mat A numeric genes x samples value matrix (e.g. from
#'   [expr_value_matrix()]).
#' @param ids Authored gene-set ids (matched against `rownames(mat)`).
#' @param counts Optional raw counts matrix for the `only_expressed` filter.
#' @param zscore Per-gene z-score across samples (default TRUE).
#' @param only_expressed Drop genes with all-zero counts (default TRUE).
#' @param transform Passed to [expr_transform()], applied before z-scoring.
#' @param pseudocount Added before a log transform.
#' @return A list: `mat` (the prepared genes x samples matrix, or `NULL` when no
#'   gene survives), `n_total`, `n_present`, `n_absent`, `n_used`, `n_nonvar`
#'   (constant kept rows), `ids_used`, `zscored`.
#' @export
expr_heatmap_matrix <- function(mat, ids, counts = NULL, zscore = TRUE,
                                only_expressed = TRUE, transform = "none",
                                pseudocount = 1) {
  mat <- as.matrix(mat)
  sel <- .expr_set_rows(rownames(mat), ids, counts = counts,
                        only_expressed = only_expressed)
  n_total <- sel$n_total; n_present <- length(sel$present)
  used <- sel$used; n_used <- length(used)
  out <- list(mat = NULL, n_total = n_total, n_present = n_present,
              n_absent = n_total - n_present, n_used = n_used, n_nonvar = 0L,
              ids_used = used, zscored = isTRUE(zscore))
  if (!n_used) return(out)

  sub <- expr_transform(mat[used, , drop = FALSE], transform, pseudocount)
  rng <- apply(sub, 1, function(r) { r <- r[is.finite(r)]; if (!length(r)) 0 else diff(range(r)) })
  out$n_nonvar <- sum(rng == 0)
  if (isTRUE(zscore)) sub <- row_zscore(sub)
  out$mat <- sub
  out
}

#' Display labels for heatmap rows or columns (duplicate- and NA-safe)
#'
#' Builds the vector of *display* labels for a set of matrix dimnames without
#' touching the matrix's actual `dimnames` -- the matrix keeps the unique dds ids
#' as its row/column names (the stable key for subsetting, `anno_mark` indexing,
#' and clustering), while ComplexHeatmap's `row_labels=`/`column_labels=` receive
#' this vector. Because it is display-only, duplicate labels (many gene ids share
#' one `gene_name`) are fine. When `source` is a metadata column, a row whose value
#' is `NA`/empty falls back to its id, so every element still gets a label.
#'
#' @param keys The matrix dimnames to label (gene ids for rows, sample ids for
#'   columns).
#' @param source `"__id__"` to use `keys` verbatim, or a `meta` column name.
#' @param meta The `rowData`/`colData` data frame (or `NULL`).
#' @param meta_keys Row keys of `meta` aligned to the dds (default
#'   `rownames(meta)`); `keys` are matched against these.
#' @return A character vector the same length/order as `keys`.
#' @export
expr_heatmap_labels <- function(keys, source = "__id__", meta = NULL,
                                meta_keys = rownames(meta)) {
  keys <- as.character(keys)
  if (identical(source, "__id__") || is.null(meta) || !source %in% colnames(meta))
    return(keys)
  lab <- as.character(meta[[source]])[match(keys, meta_keys)]
  bad <- is.na(lab) | !nzchar(trimws(lab))
  lab[bad] <- keys[bad]
  lab
}

#' Default label-display mode for a heatmap axis
#'
#' Show all labels when the axis has few elements, otherwise none -- the sensible
#' default before the user overrides it. Rows (genes) and columns (samples) use
#' their own thresholds at the call site.
#'
#' @param n Number of elements on the axis.
#' @param max_n Threshold at/below which all labels are shown.
#' @return `"all"` or `"none"`.
#' @export
heatmap_label_default <- function(n, max_n) if (n <= max_n) "all" else "none"

#' How many searched labels can actually be shown
#'
#' For the "show selected" label mode: the user searches ids from the whole dataset,
#' but only those present on the plotted axis (a gene in the chosen set, a sample in
#' the current "Showing" subset) can be marked. Reports the split so the UI can note
#' "X of Y labels cannot be shown".
#'
#' @param selected Searched ids.
#' @param present_keys The axis dimnames actually plotted.
#' @return A list: `n_selected`, `n_shown`, `n_hidden`.
#' @export
expr_label_coverage <- function(selected, present_keys) {
  selected <- unique(as.character(selected)); selected <- selected[!is.na(selected)]
  n_shown <- sum(selected %in% present_keys)
  list(n_selected = length(selected), n_shown = n_shown,
       n_hidden = length(selected) - n_shown)
}

#' k-means cluster the rows of a value matrix (reproducibly)
#'
#' Computed *outside* `ComplexHeatmap::Heatmap()` so the caller controls the seed,
#' keeps the membership (to save clusters as gene sets), and can label slices with
#' member counts. Clusters are relabelled by decreasing size (cluster 1 = largest)
#' for a stable, meaningful order. The global RNG state is saved and restored, so
#' seeding here never perturbs the rest of the session.
#'
#' @param mat A numeric rows x cols matrix (e.g. the displayed heatmap matrix).
#'   Non-finite values are treated as 0 before clustering.
#' @param k Number of clusters. `< 2` (or fewer than 2 rows) returns `NULL`
#'   (no split). Clamped to the number of rows.
#' @param seed Integer seed for reproducibility, or `NA`/`NULL` for **no seed** --
#'   a fresh (non-reproducible) clustering each call, which advances the global
#'   RNG. A concrete seed is RNG-safe (state saved/restored).
#' @param nstart Passed to [stats::kmeans()] (random restarts).
#' @return A named integer vector (names = `rownames(mat)`) of cluster ids, or
#'   `NULL` when no split applies / clustering fails.
#' @export
expr_kmeans <- function(mat, k, seed = 1L, nstart = 10L) {
  mat <- as.matrix(mat)
  k <- suppressWarnings(as.integer(k)[1]); n <- nrow(mat)
  if (is.na(k) || k < 2L || n < 2L) return(NULL)
  if (!all(is.finite(mat))) mat[!is.finite(mat)] <- 0
  k <- min(k, n)
  relabel <- function(cl) {                      # cluster 1 = largest, stable order
    tab <- sort(table(cl), decreasing = TRUE)
    remap <- stats::setNames(seq_along(tab), names(tab))
    out <- unname(remap[as.character(cl)]); names(out) <- rownames(mat); out
  }
  if (k >= n) return(relabel(seq_len(n)))        # each row its own cluster
  seed <- suppressWarnings(as.integer(seed)[1])
  if (!is.na(seed)) {                            # fixed seed -> reproducible + RNG-safe
    old <- if (exists(".Random.seed", envir = .GlobalEnv))
      get(".Random.seed", envir = .GlobalEnv) else NULL
    on.exit(if (is.null(old)) suppressWarnings(rm(".Random.seed", envir = .GlobalEnv))
            else assign(".Random.seed", old, envir = .GlobalEnv), add = TRUE)
    set.seed(seed)
  }                                              # else: no seed -> random, RNG advances
  km <- tryCatch(stats::kmeans(mat, centers = k, nstart = nstart),
                 error = function(e) NULL)
  if (is.null(km)) return(NULL)
  relabel(km$cluster)
}

#' Slice labels carrying member counts, e.g. "C1\\n(23)"
#'
#' Turns a cluster-assignment vector into a factor whose levels are the two-line
#' slice titles ComplexHeatmap shows for `row_split`/`column_split`: a `<prefix><id>`
#' line over a `(count)` line. Numeric cluster ids order numerically; levels are
#' ordered so the slices read C1, C2, C3, ...
#'
#' @param clusters A cluster-assignment vector (from [expr_kmeans()]).
#' @param prefix Label prefix before the cluster id (default `"C"` -> "C1").
#' @return A factor aligned to `clusters`, levels = the `"<prefix><id>\n(count)"`
#'   titles.
#' @export
split_with_counts <- function(clusters, prefix = "C") {
  cl <- as.character(clusters)
  cnt <- table(cl)
  lv <- names(cnt)
  num <- suppressWarnings(as.numeric(lv))
  lv <- if (!any(is.na(num))) lv[order(num)] else lv[order(lv)]
  lab_for <- function(x) sprintf("%s%s\n(%d)", prefix, x, as.integer(cnt[x]))
  factor(lab_for(cl), levels = lab_for(lv))
}

#' Cluster membership as a named list (cluster -> ids)
#'
#' Groups the ids (names of a [expr_kmeans()] vector) by cluster, in cluster order
#' -- the basis for "save each row cluster as a gene set".
#'
#' @param clusters A named cluster-assignment vector (names = ids).
#' @return A named list keyed by cluster id, each a character vector of ids.
#' @export
cluster_membership <- function(clusters) {
  ids <- names(clusters)
  if (is.null(ids)) ids <- as.character(seq_along(clusters))
  cl <- as.character(clusters)
  num <- suppressWarnings(as.numeric(unique(cl)))
  keys <- if (!any(is.na(num))) as.character(sort(unique(num))) else sort(unique(cl))
  stats::setNames(lapply(keys, function(k) ids[cl == k]), keys)
}

#' Symmetric colour limits around zero
#'
#' For a z-scored heatmap the divergent ramp should be centred at 0, so the limits
#' are `c(-M, M)` with `M` the largest absolute value. Used when the user leaves the
#' ramp anchors blank in z-score mode.
#'
#' @param values Numeric values (non-finite dropped).
#' @return `c(-M, M)`; `c(-1, 1)` when there is nothing finite / all zero.
#' @export
expr_symmetric_limits <- function(values) {
  v <- values[is.finite(values)]
  if (!length(v)) return(c(-1, 1))
  m <- max(abs(v)); if (m == 0) m <- 1
  c(-m, m)
}
