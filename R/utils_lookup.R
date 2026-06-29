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
#' @param case_insensitive Match case-insensitively (default `FALSE`).
#' @return Character vector of matched feature IDs (NA where unmatched), same
#'   length as `query`.
#' @export
lookup_feature <- function(query, row_data, ids = rownames(row_data),
                           feature_type = NULL, column = NULL,
                           case_insensitive = FALSE) {
  row_data <- as.data.frame(row_data, optional = TRUE)
  name_col <- column
  if (is.null(name_col) && !is.null(feature_type)) {
    name_col <- paste0(feature_type, "_name")
  }
  values <- if (!is.null(name_col) && name_col %in% colnames(row_data)) {
    as.character(row_data[[name_col]])
  } else {
    ids
  }
  q <- as.character(query); v <- values
  if (isTRUE(case_insensitive)) { q <- tolower(q); v <- tolower(v) }
  ids[match(q, v)]
}

#' Resolve a single query to a feature id, with a match count
#'
#' Searches `query` against `values` (a per-feature character vector, e.g. a
#' `rowData` column or the feature ids themselves) and returns the *first*
#' matching id together with how many features matched -- so callers can warn
#' "N features matched; showing the first" on a duplicate hit.
#'
#' @param query A single search string.
#' @param values Character vector of searchable values, one per feature.
#' @param ids Feature ids aligned to `values` (default `values`).
#' @param case_insensitive Match case-insensitively (default `FALSE`).
#' @return A list: `id` (first matching id, or `NA_character_`) and `n` (number
#'   of matches).
#' @export
resolve_feature <- function(query, values, ids = values,
                            case_insensitive = FALSE) {
  q <- as.character(query)[1]
  v <- as.character(values)
  if (isTRUE(case_insensitive)) { q <- tolower(q); v <- tolower(v) }
  hit <- which(v == q)
  list(id = if (length(hit)) as.character(ids)[hit[1]] else NA_character_,
       n = length(hit))
}

#' Suggest feature names for a near-miss query ("Did you mean ...?")
#'
#' On an exact-match miss, offers ranked partial matches. Matching is **always
#' case-insensitive** (independent of any case-sensitive lookup toggle), so a
#' case-sensitive search for `"duxf3"` still suggests `"Duxf3"`. Matches are
#' ranked by position -- prefix (`^q`) first, then word-start (`q` after a
#' non-alphanumeric), then any substring -- with ties broken by shorter length
#' then alphabetically. Guards keep it cheap: a query under `min_chars`
#' characters yields nothing, and more than `cap` raw matches returns
#' `over_cap = TRUE` (too broad to list).
#'
#' @param query A single search string.
#' @param values Character vector of searchable values (e.g. feature names).
#' @param n Maximum number of suggestions to return (default 5).
#' @param cap Raw-match ceiling above which the result is "too broad" (default 100).
#' @param min_chars Minimum query length to search (default 2).
#' @return A list: `suggestions` (ranked, de-duplicated, up to `n`), `n_match`
#'   (distinct matches), and `over_cap` (logical, raw matches exceeded `cap`).
#' @export
suggest_features <- function(query, values, n = 5L, cap = 100L, min_chars = 2L) {
  out <- list(suggestions = character(0), n_match = 0L, over_cap = FALSE)
  q <- tolower(trimws(as.character(query)[1]))
  if (length(q) != 1L || is.na(q) || nchar(q) < min_chars) return(out)
  v <- as.character(values)
  v <- v[!is.na(v) & nzchar(v)]
  if (!length(v)) return(out)
  vl <- tolower(v)
  hit <- which(grepl(q, vl, fixed = TRUE))
  if (!length(hit)) return(out)
  if (length(hit) > cap) { out$over_cap <- TRUE; return(out) }
  mv <- v[hit]; ml <- vl[hit]
  keep <- !duplicated(mv); mv <- mv[keep]; ml <- ml[keep]  # de-dup display values
  out$n_match <- length(mv)
  pos <- as.integer(regexpr(q, ml, fixed = TRUE))
  before <- ifelse(pos > 1L, substr(ml, pos - 1L, pos - 1L), "")
  rank1 <- ifelse(pos == 1L, 0L, ifelse(grepl("[^[:alnum:]]", before), 1L, 2L))
  ord <- order(rank1, nchar(mv), mv)
  out$suggestions <- mv[ord][seq_len(min(as.integer(n), length(mv)))]
  out
}

#' Candidate feature-search fields for the gene picker
#'
#' Builds the `selectInput` choices for "Search by" in the PCA gene-colour box:
#' `rowData` columns minus **logical** columns (numeric ids like Entrez are
#' kept), optionally minus columns with duplicate values, plus a always-present
#' "Feature ID" entry (rownames, unique by construction).
#'
#' @param row_data A `rowData`/data.frame of feature metadata.
#' @param include_duplicates Keep columns with duplicate values (default `TRUE`).
#' @return A named character vector (label -> value); the feature-id entry's
#'   value is the sentinel `"__rownames__"`.
#' @export
feature_search_choices <- function(row_data, include_duplicates = TRUE) {
  rd <- as.data.frame(row_data, optional = TRUE)
  cols <- names(rd)
  if (length(cols)) {
    is_logical <- vapply(rd, is.logical, logical(1))
    cols <- cols[!is_logical]
  }
  if (length(cols) && !isTRUE(include_duplicates)) {
    nodup <- vapply(rd[cols], function(x) anyDuplicated(x) == 0L, logical(1))
    cols <- cols[nodup]
  }
  c(stats::setNames(cols, cols), c("Feature ID (row names)" = "__rownames__"))
}
