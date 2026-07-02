# Reusable free-text gene-search control -- resolve typed names/IDs to canonical
# feature ids against rowData. Extracted from the inline blocks in the PCA
# (mod_dimreduc) and DE (mod_de) pages (rule of three: Gene Sets is the third
# consumer) so all pages share the explicit "Search by" column picker + the
# case-insensitive toggle + the "did you mean ...?" hint.
#
# Scope is the query -> ids path ONLY: it does NOT own the expression-value
# controls (assay/transform/pseudocount -- that is mod_expr_value) nor any
# page-specific gating (e.g. PCA's "only when colour-by = gene"). Like
# mod_expr_value / mod_plot_subset it runs in the *host* module's namespace: a
# page calls gene_search_ui(ns, suffix, ...) in a sidebar and
# gene_search_server(input, output, session, state, suffix, ...) once, which
# returns a reactive resolving to list(ids, records, unmatched). `suffix`
# distinguishes instances on a page.
#
# `multiple = FALSE` (PCA): a single textInput -> at most one id.
# `multiple = TRUE`  (DE labels, Gene Sets paste): a textAreaInput split on
#   commas/newlines -> one or more ids, with an `unmatched` report.

# The suffixed input/output ids for one instance.
.gene_search_ids <- function(suffix) {
  list(dup         = paste0(suffix, "_dup"),
       searchby    = paste0(suffix, "_searchby"),
       searchby_ui = paste0(suffix, "_searchby_ui"),
       ci          = paste0(suffix, "_ci"),
       q           = paste0(suffix, "_q"),
       hint        = paste0(suffix, "_hint"))
}

# The controls, in the PCA order: (optional) duplicate-columns switch, "Search
# by" selector, case-insensitive switch, the query box, the inline hint. Embed
# inside a sidebar / accordion panel / conditionalPanel as needed. `dup_toggle`
# exposes the duplicate-value-columns switch (PCA); `dup_default` is its initial
# value (the host may compute it from rowData width). `multiple` swaps the single
# textInput for a multi-line textAreaInput.
gene_search_ui <- function(ns, suffix, multiple = FALSE, dup_toggle = FALSE,
                           dup_default = TRUE, label = NULL,
                           placeholder = if (multiple) "e.g. Actb, Gapdh, ENSG..." else "e.g. Actb") {
  ids <- .gene_search_ids(suffix)
  query_input <- if (isTRUE(multiple))
    textAreaInput(ns(ids$q), label = label, placeholder = placeholder, rows = 3)
  else
    textInput(ns(ids$q), label = label, placeholder = placeholder)
  tagList(
    if (isTRUE(dup_toggle)) bslib::input_switch(
      ns(ids$dup),
      bslib::tooltip(
        tags$span("Include columns with duplicate values"),
        "Off hides rowData columns whose values are not unique (e.g. gene names that repeat). When such a column is searched, the first matching feature is used."),
      value = isTRUE(dup_default)),
    uiOutput(ns(ids$searchby_ui)),
    bslib::input_switch(ns(ids$ci), "Case-insensitive search", value = FALSE),
    query_input,
    uiOutput(ns(ids$hint))
  )
}

# Wire the controls + return a reactive resolving to
#   list(ids = <unique matched feature ids, subset of rownames>,
#        records = <per matched query: list(query, id, match, n)>,
#        unmatched = <queries with no hit>).
# In single mode `ids`/`records` have length <= 1; records[[1]] carries the
# matched value + hit count for a caption/legend.
gene_search_server <- function(input, output, session, state, suffix,
                               multiple = FALSE) {
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

  # Debounce the free-text box so the lookup + suggestion search don't fire on
  # every keystroke (the suggestion regex is O(n features) on a miss).
  gene_q <- debounce(reactive(input[[ids$q]] %||% ""), 300)
  # The character vector searched for the chosen "Search by" field.
  search_values <- reactive({
    req(state$working)
    field <- input[[ids$searchby]] %||% "__rownames__"
    if (identical(field, "__rownames__")) return(rownames(state$working))
    rd <- as.data.frame(SummarizedExperiment::rowData(state$working), optional = TRUE)
    if (field %in% names(rd)) as.character(rd[[field]]) else rownames(state$working)
  })

  # Resolve the (debounced) query into ids/records/unmatched.
  resolved <- reactive({
    empty <- list(ids = character(0), records = list(), unmatched = character(0))
    if (is.null(state$working)) return(empty)
    ci   <- isTRUE(input[[ids$ci]])
    vals <- search_values(); rn <- rownames(state$working)
    terms <- if (isTRUE(multiple)) {
      t <- trimws(strsplit(trimws(gene_q()), "[,\n]+")[[1]]); t[nzchar(t)]
    } else {
      q <- trimws(gene_q()); if (nzchar(q)) q else character(0)
    }
    if (!length(terms)) return(empty)
    recs <- list(); hit_ids <- character(0); miss <- character(0)
    for (t in terms) {
      rf <- resolve_feature(t, vals, rn, case_insensitive = ci)
      if (!is.na(rf$id)) {
        recs[[length(recs) + 1L]] <- c(list(query = t), rf)
        hit_ids <- c(hit_ids, rf$id)
      } else miss <- c(miss, t)
    }
    list(ids = unique(hit_ids), records = recs, unmatched = miss)
  })

  # Inline feedback under the box.
  output[[ids$hint]] <- renderUI({
    note <- function(cls, ...) tags$div(class = paste("small mb-2", cls), ...)
    q <- trimws(gene_q()); if (!nzchar(q)) return(NULL)
    r <- resolved()
    if (isTRUE(multiple)) {
      if (!length(r$unmatched)) return(NULL)
      return(note("text-primary", sprintf("Not found: %s", paste(r$unmatched, collapse = ", "))))
    }
    # single mode: duplicate-hit note, "did you mean ...?", or not found.
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
