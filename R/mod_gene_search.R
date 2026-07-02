# Reusable free-text gene-search control -- resolve typed names/IDs to canonical
# feature ids against rowData. Extracted from the inline blocks in the PCA
# (mod_dimreduc) and DE (mod_de) pages (rule of three: Gene Sets is the third
# consumer) so all pages share the explicit "Search by" column picker, the
# case-insensitive toggle, the match-mode selector, and the resolution hints.
#
# Scope is the query -> ids path ONLY: it does NOT own the expression-value
# controls (assay/transform/pseudocount -- that is mod_expr_value) nor any
# page-specific gating (e.g. PCA's "only when colour-by = gene"). Like
# mod_expr_value / mod_plot_subset it runs in the *host* module's namespace: a
# page calls gene_search_ui(ns, suffix, ...) in a sidebar and
# gene_search_server(input, output, session, state, suffix, ...) once, which
# returns a reactive resolving to list(ids, records, unmatched, ...). `suffix`
# distinguishes instances on a page.
#
# Two orthogonal knobs:
#  * `multiple` -- how many results the *consumer* accepts. FALSE (PCA) = a
#    single textInput -> at most one id; TRUE (DE labels, Gene Sets paste) = a
#    multi-line textarea -> one or more ids + an `unmatched` report.
#  * `search_modes` -- the match *strategy* offered: any of "exact" (literal
#    name/id), "contains" (fixed substring), "regex" (pattern). PCA passes
#    "exact" only (no selector shown); DE / Gene Sets pass all three. Exact +
#    contains split the box on commas OR newlines; regex splits on newlines only
#    (a comma is valid regex, e.g. `a{2,3}`) -- one pattern per line.

# The suffixed input/output ids for one instance.
.gene_search_ids <- function(suffix) {
  list(dup         = paste0(suffix, "_dup"),
       searchby    = paste0(suffix, "_searchby"),
       searchby_ui = paste0(suffix, "_searchby_ui"),
       mode        = paste0(suffix, "_mode"),
       ci          = paste0(suffix, "_ci"),
       q           = paste0(suffix, "_q"),
       split_hint  = paste0(suffix, "_split_hint"),
       hint        = paste0(suffix, "_hint"))
}

.gene_search_mode_labels <- c(exact = "Exact", contains = "Contains", regex = "Regex")

# Split the query box into entries for the given mode. Regex splits on newlines
# only (commas are meaningful); exact/contains split on commas or newlines.
.gene_search_split <- function(txt, mode) {
  pat <- if (identical(mode, "regex")) "\n+" else "[,\n]+"
  t <- trimws(strsplit(txt, pat)[[1]])
  unique(t[nzchar(t)])
}

# Row indices matched by a contains/regex pattern (NULL = invalid regex).
.gene_search_pattern_hits <- function(pattern, values, mode, ci) {
  if (identical(mode, "contains")) {
    if (isTRUE(ci)) which(grepl(tolower(pattern), tolower(values), fixed = TRUE))
    else            which(grepl(pattern, values, fixed = TRUE))
  } else {
    # An invalid pattern warns then errors under perl; swallow both -> NULL.
    tryCatch(suppressWarnings(which(grepl(pattern, values, ignore.case = ci, perl = TRUE))),
             error = function(e) NULL)
  }
}

# The controls, in order: (optional) duplicate-columns switch, "Search by"
# selector, (optional) match-mode selector, case-insensitive switch, the query
# box, a mode-aware split hint (multi only), and the resolution hint. `dup_toggle`
# exposes the duplicate-value-columns switch (PCA); `dup_default` is its initial
# value. `multiple` swaps the single textInput for a multi-line textAreaInput.
gene_search_ui <- function(ns, suffix, multiple = FALSE, dup_toggle = FALSE,
                           dup_default = TRUE, search_modes = "exact",
                           label = NULL, placeholder = NULL, rows = 3) {
  ids <- .gene_search_ids(suffix)
  if (is.null(placeholder))
    placeholder <- if (isTRUE(multiple)) "e.g. Actb, Gapdh, ENSG..." else "e.g. Actb"
  query_input <- if (isTRUE(multiple))
    textAreaInput(ns(ids$q), label = label, placeholder = placeholder, rows = rows)
  else
    textInput(ns(ids$q), label = label, placeholder = placeholder)
  modes <- search_modes[search_modes %in% names(.gene_search_mode_labels)]
  mode_sel <- if (length(modes) > 1L)
    selectInput(ns(ids$mode), "Match mode",
                choices  = stats::setNames(modes, unname(.gene_search_mode_labels[modes])),
                selected = modes[1])
  tagList(
    if (isTRUE(dup_toggle)) bslib::input_switch(
      ns(ids$dup),
      bslib::tooltip(
        tags$span("Include columns with duplicate values"),
        "Off hides rowData columns whose values are not unique (e.g. gene names that repeat). When such a column is searched, the first matching feature is used."),
      value = isTRUE(dup_default)),
    uiOutput(ns(ids$searchby_ui)),
    mode_sel,
    bslib::input_switch(ns(ids$ci), "Case-insensitive search", value = FALSE),
    query_input,
    if (isTRUE(multiple)) uiOutput(ns(ids$split_hint)),
    uiOutput(ns(ids$hint))
  )
}

# Wire the controls + return a reactive resolving to
#   list(ids       = <unique matched feature ids, subset of rownames>,
#        records   = <per matched entry; exact: list(query,id,match,n),
#                     contains/regex: list(query, ids, n)>,
#        unmatched = <valid entries with no hit>,
#        n_query   = <number of (deduped) entries>, mode = <current mode>,
#        invalid   = <entries rejected as invalid regex>).
# In single mode `ids`/`records` have length <= 1; records[[1]] carries the
# matched value + hit count for a caption/legend.
gene_search_server <- function(input, output, session, state, suffix,
                               multiple = FALSE, search_modes = "exact") {
  ids <- .gene_search_ids(suffix)
  ns  <- session$ns
  ft  <- function() (state$meta %||% list())$feature_type %||% "feature"

  # "Search by" choices, rebuilt live (robust to conditional display); honours
  # the optional duplicate-columns switch; defaults to <feature_type>_name.
  output[[ids$searchby_ui]] <- renderUI({
    req(state$working)
    inc <- if (is.null(input[[ids$dup]])) TRUE else isTRUE(input[[ids$dup]])
    ch  <- feature_search_choices(SummarizedExperiment::rowData(state$working),
                                  include_duplicates = inc)
    def <- paste0(ft(), "_name")
    sel <- if (def %in% ch) def else "__rownames__"
    cur <- shiny::isolate(input[[ids$searchby]])
    if (!is.null(cur) && cur %in% ch) sel <- cur   # preserve across dup-toggle
    selectInput(ns(ids$searchby), "Search by", choices = ch, selected = sel)
  })

  cur_mode <- reactive({
    m <- input[[ids$mode]] %||% search_modes[1]
    if (!m %in% search_modes) m <- search_modes[1]
    m
  })
  # Mode-aware split guidance under the box (multi only).
  output[[ids$split_hint]] <- renderUI({
    if (!isTRUE(multiple)) return(NULL)
    msg <- if (identical(cur_mode(), "regex")) "One pattern per line."
           else "Separate entries with commas or new lines."
    helpText(class = "small text-muted mt-n2", msg)
  })

  # Debounce the free-text box so the lookup/pattern search don't fire on every
  # keystroke (over all features on a miss).
  gene_q <- debounce(reactive(input[[ids$q]] %||% ""), 300)
  # The character vector searched for the chosen "Search by" field.
  search_values <- reactive({
    req(state$working)
    field <- input[[ids$searchby]] %||% "__rownames__"
    if (identical(field, "__rownames__")) return(rownames(state$working))
    rd <- as.data.frame(SummarizedExperiment::rowData(state$working), optional = TRUE)
    if (field %in% names(rd)) as.character(rd[[field]]) else rownames(state$working)
  })

  # Resolve the (debounced) query into ids/records/unmatched/invalid.
  resolved <- reactive({
    empty <- list(ids = character(0), records = list(), unmatched = character(0),
                  n_query = 0L, mode = cur_mode(), invalid = character(0))
    if (is.null(state$working)) return(empty)
    mode <- cur_mode(); ci <- isTRUE(input[[ids$ci]])
    vals <- search_values(); rn <- rownames(state$working)
    entries <- if (isTRUE(multiple)) .gene_search_split(trimws(gene_q()), mode) else {
      q <- trimws(gene_q()); if (nzchar(q)) q else character(0)
    }
    if (!length(entries)) return(empty)
    recs <- list(); hit_ids <- character(0); miss <- character(0); invalid <- character(0)
    for (t in entries) {
      if (identical(mode, "exact")) {
        rf <- resolve_feature(t, vals, rn, case_insensitive = ci)
        if (!is.na(rf$id)) {
          recs[[length(recs) + 1L]] <- c(list(query = t), rf); hit_ids <- c(hit_ids, rf$id)
        } else miss <- c(miss, t)
      } else {
        hits <- .gene_search_pattern_hits(t, vals, mode, ci)
        if (is.null(hits)) { invalid <- c(invalid, t); next }         # invalid regex
        if (length(hits)) {
          mids <- rn[hits]
          recs[[length(recs) + 1L]] <- list(query = t, ids = mids, n = length(mids))
          hit_ids <- c(hit_ids, mids)
        } else miss <- c(miss, t)
      }
    }
    list(ids = unique(hit_ids), records = recs, unmatched = miss,
         n_query = length(entries), mode = mode, invalid = invalid)
  })

  # Inline feedback under the box.
  output[[ids$hint]] <- renderUI({
    note <- function(cls, ...) tags$div(class = paste("small mb-2", cls), ...)
    q <- trimws(gene_q()); if (!nzchar(q)) return(NULL)
    r <- resolved()

    if (isTRUE(multiple)) {
      nq <- r$n_query; miss <- r$unmatched; nmiss <- length(miss); nid <- length(r$ids)
      if (!nq) return(NULL)
      bits <- list()
      if (length(r$invalid))
        bits <- c(bits, list(note("text-danger",
          sprintf("Invalid regex: %s", paste(utils::head(r$invalid, 5), collapse = ", ")))))
      if (!nmiss) {                                    # everything matched
        bits <- c(bits, list(note("text-muted",
          sprintf("%d feature%s matched.", nid, if (nid == 1L) "" else "s"))))
        return(tagList(bits))
      }
      pct  <- round(100 * nmiss / nq)
      head <- sprintf("%d of %d not found (%d%%).", nmiss, nq, pct)
      # Tier: a large fraction unmatched almost always means the wrong column.
      if (nmiss >= 20L && nmiss >= 0.5 * nq) {
        bits <- c(bits, list(note("text-warning",
          paste0(head, " Are these from a different ID column? Try changing 'Search by'."))))
        return(tagList(bits))
      }
      # Tier: a couple of exact-mode misses -> per-term "did you mean" suggestions.
      if (identical(r$mode, "exact") && nmiss <= 2L) {
        vals <- search_values()
        lines <- vapply(miss, function(t) {
          sg <- suggest_features(t, vals, n = 2L)
          if (length(sg$suggestions))
            sprintf("%s (did you mean: %s?)", t, paste(sg$suggestions, collapse = ", "))
          else t
        }, character(1))
        bits <- c(bits, list(note("text-warning",
          paste0(head, " Not found: ", paste(lines, collapse = "; ")))))
        return(tagList(bits))
      }
      # Default: the capped unmatched list.
      shown <- utils::head(miss, 10L); more <- nmiss - length(shown)
      txt <- paste(shown, collapse = ", ")
      if (more > 0L) txt <- paste0(txt, sprintf(" (+%d more)", more))
      bits <- c(bits, list(note("text-warning", paste0(head, " Not found: ", txt))))
      return(tagList(bits))
    }

    # Single mode: duplicate-hit note, "did you mean ...?", or not found.
    if (length(r$records)) {
      rf <- r$records[[1]]
      if (rf$n > 1L)
        return(note("text-muted", sprintf("%d features matched '%s'; using the first (%s).",
                                           rf$n, q, rf$match)))
      return(NULL)
    }
    sg <- suggest_features(q, search_values())  # always case-insensitive
    if (isTRUE(sg$over_cap))
      return(note("text-primary", "Too many partial matches - type more to narrow it down."))
    if (length(sg$suggestions)) {
      more <- sg$n_match - length(sg$suggestions)
      txt <- paste(sg$suggestions, collapse = ", ")
      if (more > 0L) txt <- paste0(txt, sprintf(" (+%d more)", more))
      return(note("text-primary", sprintf("Feature '%s' not found. Did you mean: %s?", q, txt)))
    }
    note("text-primary", sprintf("Feature '%s' not found.", q))
  })

  resolved
}
