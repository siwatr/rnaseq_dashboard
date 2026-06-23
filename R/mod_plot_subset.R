# Reusable view-only "Showing:" plot-subset control -- the project standard for
# plot-visualization tabs. It hides samples from a page's plots WITHOUT mutating
# the dds (no data_version bump): purely a display subset. The control is placed
# in each plot sidebar via plot_subset_ui(); all instances on a page are synced
# to one canonical selection, and plot_subset_server() returns a `showing_samples`
# reactive the page filters its plotted data by. These operate in the *host*
# module's namespace (a true Shiny submodule cannot place UI across the host's
# separate sidebars), so a page calls plot_subset_server(input, output, session,
# state, suffixes) once and plot_subset_ui(ns, suffix) in each plot sidebar.

# One sidebar's control: a "show by" column selector + a value multiselect.
# `suffix` distinguishes the instances on a page (e.g. "gen", "rle"). Choices are
# populated by the server on load. Blank value box = show all (no subsetting).
plot_subset_ui <- function(ns, suffix) {
  tagList(
    tags$hr(class = "my-2"),
    selectInput(ns(paste0(suffix, "_show_by")), "Showing (display only)",
                choices = c("All samples" = "__all__"), selected = "__all__"),
    selectizeInput(ns(paste0(suffix, "_show_values")), "Keep (blank = show all)",
                   choices = character(0), multiple = TRUE,
                   options = list(placeholder = "(blank = show all)"))
  )
}

# Wire the synced canonical selection across `suffixes` and return a reactive of
# the samples currently shown (always non-empty; a blank Keep box = all). Call
# once in a page's server. `state` is the shared app-state (see new_app_state()).
plot_subset_server <- function(input, output, session, state, suffixes) {
  show_by_rv     <- reactiveVal("__all__")
  show_values_rv <- reactiveVal(character(0))

  # Value-box choices for the current "show by" column.
  show_val_choices <- function(by) {
    if (is.null(by) || identical(by, "__all__")) return(character(0))
    if (identical(by, "__samples__")) return(colnames(state$working))
    cd <- as.data.frame(SummarizedExperiment::colData(state$working))
    sort(unique(as.character(cd[[by]])))
  }

  # Populate / reset all the per-tab controls when a dataset (re)loads.
  observeEvent(state$working, {
    cols <- colnames(SummarizedExperiment::colData(state$working))
    ch <- c("All samples" = "__all__", stats::setNames(cols, cols),
            "Individual samples" = "__samples__")
    show_by_rv("__all__"); show_values_rv(character(0))
    for (s in suffixes) {
      updateSelectInput(session, paste0(s, "_show_by"), choices = ch, selected = "__all__")
      updateSelectizeInput(session, paste0(s, "_show_values"),
                           choices = character(0), selected = character(0))
    }
  })

  # Per-tab edits -> canonical state (guarded so fan-out echoes are no-ops).
  lapply(suffixes, function(s) {
    observeEvent(input[[paste0(s, "_show_by")]], {
      v <- input[[paste0(s, "_show_by")]]
      if (is.null(v) || identical(v, show_by_rv())) return()
      show_by_rv(v); show_values_rv(character(0))        # reset values on column switch
    }, ignoreInit = TRUE)
    observeEvent(input[[paste0(s, "_show_values")]], {
      v <- input[[paste0(s, "_show_values")]] %||% character(0)
      if (setequal(v, show_values_rv())) return()
      show_values_rv(v)
    }, ignoreNULL = FALSE, ignoreInit = TRUE)
  })

  # Canonical state -> fan out to every tab's controls (keeps them in sync).
  observeEvent(show_by_rv(), {
    ch <- show_val_choices(show_by_rv())
    for (s in suffixes) {
      updateSelectInput(session, paste0(s, "_show_by"), selected = show_by_rv())
      updateSelectizeInput(session, paste0(s, "_show_values"),
                           choices = ch, selected = show_values_rv())
    }
  })
  observeEvent(show_values_rv(), {
    for (s in suffixes) {
      updateSelectizeInput(session, paste0(s, "_show_values"), selected = show_values_rv())
    }
  }, ignoreNULL = FALSE)

  # Samples currently shown (always non-empty: a blank Keep box = show all).
  reactive({
    req(state$working)
    all_s <- colnames(state$working)
    by <- show_by_rv()
    if (identical(by, "__all__")) return(all_s)
    vals <- show_values_rv()
    if (!length(vals)) return(all_s)
    if (identical(by, "__samples__")) return(intersect(all_s, vals))
    cd <- as.data.frame(SummarizedExperiment::colData(state$working))
    all_s[as.character(cd[[by]]) %in% vals]
  })
}
