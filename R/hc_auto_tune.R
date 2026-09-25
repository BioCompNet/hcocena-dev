#' Automatically tune key hCoCena parameters
#'
#' `hc_auto_tune()` recommends and optionally applies:
#' 1. Layer-wise cutoffs based on `run_expression_analysis_1()` output.
#' 2. Clustering settings based on modularity across candidate algorithms.
#'
#' The tuning report is stored in `hc@satellite$auto_tune`.
#'
#' @param hc A `HCoCenaExperiment` object.
#' @param tune_cutoff Logical; tune per-layer correlation cutoffs.
#' @param tune_clustering Logical; tune clustering algorithm / resolution.
#' @param apply Logical; if `TRUE`, apply best recommendations to `hc`.
#' @param cluster_algorithms Character vector of candidate algorithms for
#'   clustering tuning. Supported: `"auto"`, `"cluster_leiden"`,
#'   `"cluster_louvain"`, `"cluster_fast_greedy"`, `"cluster_walktrap"`,
#'   `"cluster_infomap"`, `"cluster_label_prop"`.
#' @param resolution_grid Numeric vector of Leiden resolutions to evaluate for
#'   `"auto"` and `"cluster_leiden"`.
#' @param no_of_iterations Passed to `hc_cluster_calculation()`.
#' @param partition_type Passed to `hc_cluster_calculation()`.
#'   Note: for Leiden with `"ModularityVertexPartition"`, sweeping
#'   `resolution_grid` usually does not change the partition. If you want
#'   resolution-driven module-count tuning, use `"CPMVertexPartition"` or
#'   `"RBConfigurationVertexPartition"`.
#' @param max_cluster_count_per_gene Passed to `hc_cluster_calculation()`.
#' @param module_count_bounds Optional length-2 numeric vector
#'   (`min_modules`, `max_modules`) used as a soft constraint during ranking.
#' @param prefer_resolution Tie-break preference for Leiden resolution.
#' @param verbose Logical; print progress messages.
#'
#' @return Updated `HCoCenaExperiment` with tuning report in
#'   `hc@satellite$auto_tune`.
#' @examples
#' hc <- hc_example_data("after_part1")
#' hc <- hc_auto_tune(hc, apply = FALSE)
#' @export
hc_auto_tune <- function(hc,
                         tune_cutoff = TRUE,
                         tune_clustering = TRUE,
                         apply = TRUE,
                         cluster_algorithms = c("auto", "cluster_leiden", "cluster_louvain"),
                         resolution_grid = base::seq(0.05, 0.5, by = 0.05),
                         no_of_iterations = 2,
                         partition_type = "RBConfigurationVertexPartition",
                         max_cluster_count_per_gene = 1,
                         module_count_bounds = c(3, 60),
                         prefer_resolution = 0.1,
                         verbose = TRUE) {
  if (!inherits(hc, "HCoCenaExperiment")) {
    stop("`hc` must be a `HCoCenaExperiment`.")
  }
  if (!isTRUE(tune_cutoff) && !isTRUE(tune_clustering)) {
    stop("Enable at least one of `tune_cutoff` or `tune_clustering`.")
  }

  report <- list(
    created_at = as.character(base::Sys.time()),
    applied = isTRUE(apply),
    cutoff = NULL,
    clustering = NULL
  )

  if (isTRUE(tune_cutoff)) {
    cutoff_vec <- .hc_auto_collect_cutoffs(hc)
    if (base::length(cutoff_vec) == 0) {
      warning("No optimal cutoffs found in `hc@layer_results`; skipping cutoff tuning.", call. = FALSE)
      report$cutoff <- list(status = "skipped", reason = "no_cutoff_stats")
    } else {
      report$cutoff <- list(
        status = if (isTRUE(apply)) "applied" else "recommended",
        recommended_cutoff_vector = cutoff_vec
      )
      if (isTRUE(apply)) {
        hc <- hc_set_cutoff(hc, cutoff_vector = cutoff_vec)
      }
      if (isTRUE(verbose)) {
        message(
          "hc_auto_tune(): cutoff recommendation = ",
          base::paste0(base::format(cutoff_vec, digits = 4), collapse = ", ")
        )
      }
    }
  }

  if (isTRUE(tune_clustering)) {
    pt_norm <- .hc_auto_normalize_partition_type(partition_type)
    uses_leiden <- any(cluster_algorithms %in% c("auto", "cluster_leiden"))
    resolution_unique <- unique(as.numeric(resolution_grid))
    resolution_unique <- resolution_unique[is.finite(resolution_unique)]
    if (uses_leiden &&
      pt_norm == "modularityvertexpartition" &&
      length(resolution_unique) > 1) {
      resolution_single <- .hc_first_numeric_value(prefer_resolution)
      if (!is.finite(resolution_single)) {
        resolution_single <- resolution_unique[[1]]
      }
      if (isTRUE(verbose)) {
        message(
          "hc_auto_tune(): `partition_type = \"ModularityVertexPartition\"` ",
          "usually yields resolution-invariant Leiden results; using a single ",
          "resolution value (", format(resolution_single, digits = 4), "). ",
          "Use `CPMVertexPartition` or `RBConfigurationVertexPartition` if you ",
          "want resolution-driven module-count tuning."
        )
      }
      resolution_grid <- resolution_single
    }

    graph <- hc@integration@graph
    if (is.null(graph) || !igraph::is_igraph(graph) || igraph::ecount(graph) == 0) {
      warning("No integrated graph found; skipping clustering tuning.", call. = FALSE)
      report$clustering <- list(status = "skipped", reason = "missing_integrated_graph")
    } else {
      bounds <- .hc_auto_parse_module_bounds(module_count_bounds)
      candidates <- .hc_auto_cluster_candidates(cluster_algorithms, resolution_grid)
      run_tbl <- vector("list", base::nrow(candidates))
      best <- list(score = -Inf, n_genes = -Inf, hc = NULL, row = NULL)

      for (i in base::seq_len(base::nrow(candidates))) {
        algo <- as.character(candidates$algorithm[[i]])
        res <- as.numeric(candidates$resolution[[i]])
        res_use <- if (base::is.na(res)) prefer_resolution else res

        candidate <- tryCatch(
          hc_cluster_calculation(
            hc,
            cluster_algo = algo,
            no_of_iterations = no_of_iterations,
            resolution = res_use,
            partition_type = partition_type,
            max_cluster_count_per_gene = max_cluster_count_per_gene,
            return_result = FALSE
          ),
          error = function(e) e
        )

        if (inherits(candidate, "error")) {
          run_tbl[[i]] <- base::data.frame(
            algorithm = algo,
            resolution = if (base::is.na(res)) NA_real_ else res,
            modularity = NA_real_,
            n_modules = NA_real_,
            n_genes_in_modules = NA_real_,
            score = NA_real_,
            status = "error",
            message = conditionMessage(candidate),
            stringsAsFactors = FALSE
          )
          next
        }

        metrics <- .hc_auto_cluster_metrics(candidate)
        score <- .hc_auto_cluster_score(
          modularity = metrics$modularity,
          n_modules = metrics$n_modules,
          resolution = if (base::is.na(res)) prefer_resolution else res,
          lower_modules = bounds$lower,
          upper_modules = bounds$upper,
          prefer_resolution = prefer_resolution
        )

        run_tbl[[i]] <- base::data.frame(
          algorithm = algo,
          resolution = if (base::is.na(res)) NA_real_ else res,
          modularity = metrics$modularity,
          n_modules = metrics$n_modules,
          n_genes_in_modules = metrics$n_genes_in_modules,
          score = score,
          status = "ok",
          message = "",
          stringsAsFactors = FALSE
        )

        is_better <- is.finite(score) && (
          (score > best$score + 1e-12) ||
            (base::abs(score - best$score) <= 1e-12 && metrics$n_genes_in_modules > best$n_genes)
        )
        if (is_better) {
          best$score <- score
          best$n_genes <- metrics$n_genes_in_modules
          best$hc <- candidate
          best$row <- run_tbl[[i]]
        }
      }

      run_tbl <- run_tbl[!base::vapply(run_tbl, is.null, logical(1))]
      run_tbl_df <- if (base::length(run_tbl) == 0) {
        base::data.frame()
      } else {
        base::do.call(base::rbind, run_tbl)
      }
      if (base::nrow(run_tbl_df) > 0) {
        ord <- base::order(
          -ifelse(base::is.finite(run_tbl_df$score), run_tbl_df$score, -Inf),
          -ifelse(base::is.finite(run_tbl_df$n_genes_in_modules), run_tbl_df$n_genes_in_modules, -Inf),
          ifelse(base::is.finite(run_tbl_df$n_modules), run_tbl_df$n_modules, Inf)
        )
        run_tbl_df <- run_tbl_df[ord, , drop = FALSE]
      }

      if (is.null(best$hc)) {
        report$clustering <- list(
          status = "failed",
          candidates = run_tbl_df
        )
      } else {
        if (isTRUE(apply)) {
          hc <- best$hc
        }
        report$clustering <- list(
          status = if (isTRUE(apply)) "applied" else "recommended",
          best = best$row,
          candidates = run_tbl_df
        )
        if (isTRUE(verbose)) {
          msg <- paste0(
            "hc_auto_tune(): clustering recommendation = ", best$row$algorithm,
            if (!is.na(best$row$resolution)) paste0(" (resolution ", best$row$resolution, ")") else "",
            ", modularity=", base::round(best$row$modularity, 4),
            ", modules=", best$row$n_modules
          )
          message(msg)
        }
      }
    }
  }

  sat <- as.list(hc@satellite)
  sat[["auto_tune"]] <- report
  hc@satellite <- S4Vectors::SimpleList(sat)
  methods::validObject(hc)
  hc
}

.hc_auto_normalize_partition_type <- function(x) {
  if (is.null(x) || length(x) == 0) {
    return("")
  }
  pt <- as.character(x[[1]])
  if (!nzchar(pt) || is.na(pt)) {
    return("")
  }
  tolower(gsub("[[:space:]]+", "", pt))
}

.hc_auto_collect_cutoffs <- function(hc) {
  layer_ids <- character(0)
  if (base::nrow(hc@config@layer) > 0 && "layer_id" %in% base::colnames(hc@config@layer)) {
    layer_ids <- as.character(hc@config@layer$layer_id)
  } else if (base::length(hc@layer_results) > 0) {
    layer_ids <- base::names(hc@layer_results)
  }
  if (base::length(layer_ids) == 0) {
    return(numeric(0))
  }

  out <- base::rep(NA_real_, base::length(layer_ids))
  names(out) <- layer_ids
  for (i in base::seq_along(layer_ids)) {
    lid <- layer_ids[[i]]
    lr <- hc@layer_results[[lid]]
    out[[i]] <- .hc_auto_optimal_cutoff_from_layer_result(lr)
  }

  if ("cutoff" %in% base::colnames(hc@config@layer)) {
    existing <- .hc_as_numeric_safely(hc@config@layer$cutoff)
    n <- base::min(base::length(existing), base::length(out))
    repl <- !base::is.finite(out[base::seq_len(n)]) & base::is.finite(existing[base::seq_len(n)])
    out[base::seq_len(n)][repl] <- existing[base::seq_len(n)][repl]
  }

  if (base::any(!base::is.finite(out))) {
    return(numeric(0))
  }
  as.numeric(out)
}

.hc_auto_optimal_cutoff_from_layer_result <- function(layer_result) {
  if (is.null(layer_result)) {
    return(NA_real_)
  }

  part1 <- NULL
  if (inherits(layer_result, "HCoCenaLayerResult")) {
    part1 <- as.list(layer_result@part1)
  } else if (base::is.list(layer_result) && "part1" %in% base::names(layer_result)) {
    part1 <- layer_result[["part1"]]
  }
  if (is.null(part1) || !base::is.list(part1)) {
    return(NA_real_)
  }
  .hc_auto_find_optimal_cutoff(part1)
}

.hc_auto_find_optimal_cutoff <- function(x) {
  if (is.null(x)) {
    return(NA_real_)
  }
  if (base::is.list(x)) {
    if ("optimal_cutoff" %in% base::names(x)) {
      val <- .hc_first_numeric_value(x[["optimal_cutoff"]])
      if (base::is.finite(val)) {
        return(val)
      }
    }
    for (nm in base::names(x)) {
      val <- .hc_auto_find_optimal_cutoff(x[[nm]])
      if (base::is.finite(val)) {
        return(val)
      }
    }
  }
  NA_real_
}

.hc_auto_parse_module_bounds <- function(module_count_bounds) {
  if (is.null(module_count_bounds)) {
    return(list(lower = -Inf, upper = Inf))
  }
  if (!is.numeric(module_count_bounds) || base::length(module_count_bounds) != 2) {
    stop("`module_count_bounds` must be NULL or a numeric vector of length 2.")
  }
  lower <- base::min(module_count_bounds)
  upper <- base::max(module_count_bounds)
  if (!base::is.finite(lower) || !base::is.finite(upper)) {
    stop("`module_count_bounds` values must be finite.")
  }
  list(lower = lower, upper = upper)
}

.hc_auto_cluster_candidates <- function(cluster_algorithms, resolution_grid) {
  supported <- c(
    "auto",
    "cluster_leiden",
    "cluster_louvain",
    "cluster_fast_greedy",
    "cluster_walktrap",
    "cluster_infomap",
    "cluster_label_prop"
  )
  if (!is.character(cluster_algorithms) || base::length(cluster_algorithms) == 0) {
    stop("`cluster_algorithms` must be a non-empty character vector.")
  }
  cluster_algorithms <- unique(cluster_algorithms)
  invalid <- base::setdiff(cluster_algorithms, supported)
  if (base::length(invalid) > 0) {
    stop("Unsupported algorithm(s): ", base::paste0(invalid, collapse = ", "))
  }

  if (!is.numeric(resolution_grid) || base::length(resolution_grid) == 0) {
    stop("`resolution_grid` must be a non-empty numeric vector.")
  }
  resolution_grid <- sort(unique(as.numeric(resolution_grid)))
  resolution_grid <- resolution_grid[base::is.finite(resolution_grid)]
  if (base::length(resolution_grid) == 0) {
    stop("`resolution_grid` contains no finite values.")
  }

  rows <- vector("list", base::length(cluster_algorithms))
  for (i in base::seq_along(cluster_algorithms)) {
    algo <- cluster_algorithms[[i]]
    if (algo %in% c("auto", "cluster_leiden")) {
      rows[[i]] <- base::data.frame(
        algorithm = base::rep(algo, base::length(resolution_grid)),
        resolution = resolution_grid,
        stringsAsFactors = FALSE
      )
    } else {
      rows[[i]] <- base::data.frame(
        algorithm = algo,
        resolution = NA_real_,
        stringsAsFactors = FALSE
      )
    }
  }
  base::do.call(base::rbind, rows)
}

.hc_auto_cluster_metrics <- function(hc) {
  graph <- hc@integration@graph
  if (is.null(graph) || !igraph::is_igraph(graph) || igraph::vcount(graph) == 0) {
    return(list(modularity = NA_real_, n_modules = 0L, n_genes_in_modules = 0L))
  }

  cluster_info <- NULL
  if (base::length(hc@integration@cluster) > 0) {
    if ("cluster_information" %in% base::names(hc@integration@cluster)) {
      cluster_info <- hc@integration@cluster[["cluster_information"]]
    } else {
      cluster_info <- hc@integration@cluster[[1]]
    }
  }

  if (inherits(cluster_info, "DataFrame")) {
    cluster_info <- base::as.data.frame(cluster_info, stringsAsFactors = FALSE)
  }
  if (!base::is.data.frame(cluster_info) || base::nrow(cluster_info) == 0) {
    return(list(modularity = NA_real_, n_modules = 0L, n_genes_in_modules = 0L))
  }

  membership <- .hc_auto_membership_from_cluster_info(cluster_info, graph)
  weights <- tryCatch(igraph::E(graph)$weight, error = function(e) NULL)
  modularity <- tryCatch(
    igraph::modularity(graph, membership = membership, weights = weights),
    error = function(e) NA_real_
  )

  keep <- base::rep(TRUE, base::nrow(cluster_info))
  if ("cluster_included" %in% base::colnames(cluster_info)) {
    keep <- keep & base::tolower(as.character(cluster_info$cluster_included)) == "yes"
  }
  if ("clusters" %in% base::colnames(cluster_info)) {
    keep <- keep & !base::grepl("cluster\\s*0$", base::tolower(as.character(cluster_info$clusters)))
  }

  n_modules <- base::sum(keep, na.rm = TRUE)
  n_genes_in_modules <- 0L
  if ("gene_no" %in% base::colnames(cluster_info)) {
    gene_no <- .hc_as_numeric_safely(cluster_info$gene_no)
    n_genes_in_modules <- base::sum(gene_no[keep], na.rm = TRUE)
  } else if ("gene_n" %in% base::colnames(cluster_info)) {
    genes <- unlist(
      base::strsplit(
        base::paste0(cluster_info$gene_n[keep], collapse = ","),
        ",",
        fixed = TRUE
      ),
      use.names = FALSE
    )
    genes <- base::trimws(genes)
    genes <- genes[genes != ""]
    n_genes_in_modules <- base::length(base::unique(genes))
  }

  list(
    modularity = .hc_first_numeric_value(modularity),
    n_modules = as.integer(n_modules),
    n_genes_in_modules = as.integer(n_genes_in_modules)
  )
}

.hc_auto_membership_from_cluster_info <- function(cluster_info, graph) {
  vertices <- igraph::V(graph)$name
  membership <- base::integer(base::length(vertices))
  names(membership) <- vertices

  if (!"gene_n" %in% base::colnames(cluster_info)) {
    return(membership)
  }

  cluster_id <- base::seq_len(base::nrow(cluster_info))
  if ("clusters" %in% base::colnames(cluster_info)) {
    parsed <- .hc_as_integer_safely(
      base::sub(".*?(\\d+)$", "\\1", as.character(cluster_info$clusters), perl = TRUE)
    )
    ok <- base::is.finite(parsed)
    cluster_id[ok] <- parsed[ok]
  }

  for (i in base::seq_len(base::nrow(cluster_info))) {
    genes_raw <- as.character(cluster_info$gene_n[[i]])
    if (base::is.na(genes_raw) || genes_raw == "") {
      next
    }
    genes <- base::strsplit(genes_raw, ",", fixed = TRUE)[[1]]
    genes <- base::trimws(genes)
    genes <- genes[genes != ""]
    genes <- base::intersect(genes, vertices)
    if (base::length(genes) == 0) {
      next
    }
    membership[genes] <- cluster_id[[i]]
  }
  membership
}

.hc_auto_cluster_score <- function(modularity,
                                   n_modules,
                                   resolution,
                                   lower_modules,
                                   upper_modules,
                                   prefer_resolution) {
  if (!base::is.finite(modularity)) {
    return(-Inf)
  }
  penalty_modules <- 0
  if (base::is.finite(lower_modules) && n_modules < lower_modules) {
    penalty_modules <- penalty_modules + (lower_modules - n_modules)
  }
  if (base::is.finite(upper_modules) && n_modules > upper_modules) {
    penalty_modules <- penalty_modules + (n_modules - upper_modules)
  }
  penalty_resolution <- if (base::is.finite(resolution)) base::abs(resolution - prefer_resolution) else 0
  modularity - 0.03 * penalty_modules - 0.002 * penalty_resolution
}
