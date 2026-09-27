#' Drop results that a change has just made obsolete
#'
#' An `HCoCenaExperiment` carries the whole analysis, so changing something
#' near the start leaves everything computed from it standing. The object then
#' looks finished while its network and modules belong to inputs that are gone;
#' worse, the configuration can disagree with the results it is stored next to
#' (a cutoff of 0.95 recorded beside a network filtered at 0.604).
#'
#' The analysis is a chain:
#'
#'   data -> correlation (part1) -> layer network + GFC (part2)
#'        -> integration -> modules -> downstream analyses (satellite)
#'
#' `from` names the first stage that is no longer valid; everything after it is
#' removed, so the next function to need it fails with a message telling the
#' user what to re-run, instead of silently using the previous state.
#'
#' Stages:
#' * `"data"` - counts, annotations or the genes entering the correlation
#'   changed. Nothing computed survives.
#' * `"cutoff"` - only the correlation cutoff changed. The correlations
#'   themselves are still valid, so `part1` is kept and `part2` onwards goes.
#' * `"integration"` - the integrated network changed; modules and downstream
#'   results go.
#' * `"clustering"` - the module assignment changed; downstream results go.
#'
#' @param hc A `HCoCenaExperiment`.
#' @param from One of `"data"`, `"cutoff"`, `"integration"`, `"clustering"`.
#' @param quiet Logical. If `FALSE` (default) a message names what was dropped.
#' @return The object with the obsolete results removed.
#' @noRd
.hc_invalidate_from <- function(hc,
                                from = c("data", "cutoff", "integration", "clustering"),
                                quiet = FALSE) {
  if (!methods::is(hc, "HCoCenaExperiment")) {
    return(hc)
  }
  from <- base::match.arg(from)
  dropped <- base::character(0)

  # Empties `which_part` in every layer result; returns the updated object and
  # whether anything was removed.
  clear_part <- function(hc, which_part) {
    touched <- FALSE
    for (nm in base::names(hc@layer_results)) {
      lr <- hc@layer_results[[nm]]
      if (base::length(methods::slot(lr, which_part)) > 0) {
        methods::slot(lr, which_part) <- S4Vectors::SimpleList()
        hc@layer_results[[nm]] <- lr
        touched <- TRUE
      }
    }
    list(hc = hc, touched = touched)
  }

  if (from %in% c("data")) {
    res <- clear_part(hc, "part1")
    hc <- res$hc
    if (res$touched) dropped <- c(dropped, "correlations")
  }
  if (from %in% c("data", "cutoff")) {
    res <- clear_part(hc, "part2")
    hc <- res$hc
    if (res$touched) dropped <- c(dropped, "layer networks and GFCs")
  }
  if (from %in% c("data", "cutoff", "integration")) {
    if (base::nrow(hc@integration@combined_edgelist) > 0 ||
      !base::is.null(hc@integration@graph) ||
      base::nrow(hc@integration@gfc) > 0) {
      hc@integration@combined_edgelist <- S4Vectors::DataFrame()
      hc@integration@graph <- NULL
      hc@integration@gfc <- S4Vectors::DataFrame()
      dropped <- c(dropped, "integrated network")
    }
  }
  if (base::length(hc@integration@cluster) > 0) {
    hc@integration@cluster <- S4Vectors::SimpleList()
    dropped <- c(dropped, "modules")
  }

  # Everything in `satellite` is derived from the modules, with one exception:
  # `cutoff_selection` records how the cutoff was arrived at and is what
  # `hc_set_cutoff(auto = TRUE)` reads, so a cutoff change must not erase it.
  keep <- if (base::identical(from, "data")) base::character(0) else "cutoff_selection"
  sat <- base::as.list(hc@satellite)
  drop_names <- base::setdiff(base::names(sat), keep)
  if (base::length(drop_names) > 0) {
    hc@satellite <- S4Vectors::SimpleList(sat[base::intersect(base::names(sat), keep)])
    dropped <- c(dropped, "downstream results")
  }

  if (base::length(dropped) > 0 && !isTRUE(quiet)) {
    base::message(
      "Discarding results that no longer match the current state (",
      base::paste(base::unique(dropped), collapse = ", "),
      "). Re-run the analysis from ",
      base::switch(from,
        data = "`hc_run_expression_analysis_1()`",
        cutoff = "`hc_run_expression_analysis_2()`",
        integration = "`hc_build_integrated_network()`",
        clustering = "`hc_cluster_calculation()`"
      ),
      "."
    )
  }
  hc
}

#' A cheap content fingerprint of the assays and annotations
#'
#' Used to decide whether `hc_read_data()` actually changed anything. Gene and
#' sample ids alone are not enough: re-reading the same ids with different
#' values is exactly the case that used to leave a stale network in place.
#' @noRd
.hc_assay_fingerprint <- function(hc) {
  exps <- tryCatch(MultiAssayExperiment::experiments(hc@mae), error = function(e) NULL)
  if (base::is.null(exps) || base::length(exps) == 0) {
    return(NULL)
  }
  base::lapply(base::names(exps), function(nm) {
    se <- exps[[nm]]
    m <- tryCatch(SummarizedExperiment::assay(se), error = function(e) NULL)
    list(
      layer = nm,
      genes = base::rownames(se),
      samples = base::colnames(se),
      values = if (base::is.null(m)) NULL else base::sum(base::as.numeric(m), na.rm = TRUE),
      colData = tryCatch(
        base::as.data.frame(SummarizedExperiment::colData(se)),
        error = function(e) NULL
      )
    )
  })
}
