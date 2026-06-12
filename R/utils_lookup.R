# Feature ID <-> name lookup. Feature IDs (e.g. Ensembl) are not readable, so
# resolve user queries against a feature-name column, defaulting to
# `<feature_type>_name`, falling back to the feature IDs themselves.

#' Resolve queries to feature IDs by name
#'
#' @param query Character vector of names (or IDs) to look up.
#' @param row_data A data.frame-like of feature metadata (one row per feature).
#' @param ids Character vector of feature IDs, one per row of `row_data`
#'   (default `rownames(row_data)`).
#' @param feature_type Feature type (e.g. "gene"); the name column searched is
#'   `paste0(feature_type, "_name")` unless `column` is given.
#' @param column Explicit name column to search; overrides `feature_type`.
#' @return Character vector of matched feature IDs (NA where unmatched), same
#'   length as `query`.
#' @export
lookup_feature <- function(query, row_data, ids = rownames(row_data),
                           feature_type = NULL, column = NULL) {
  row_data <- as.data.frame(row_data, optional = TRUE)
  name_col <- column
  if (is.null(name_col) && !is.null(feature_type)) {
    name_col <- paste0(feature_type, "_name")
  }
  if (!is.null(name_col) && name_col %in% colnames(row_data)) {
    hits <- match(query, row_data[[name_col]])
  } else {
    hits <- match(query, ids)
  }
  ids[hits]
}
