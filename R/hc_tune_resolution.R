#' Tune clustering algorithm and resolution
#'
#' Convenience wrapper for clustering-focused auto-tuning.
#' Internally this calls `hc_auto_tune(..., tune_cutoff = FALSE, tune_clustering = TRUE)`.
#' The clustering report is additionally mirrored to `hc@satellite$resolution_tuning`.
#'
#' @param hc A `HCoCenaExperiment` object.
#' @param apply Logical; if `TRUE`, apply the best clustering recommendation.
#' @param cluster_algorithms Character vector of candidate algorithms.
#' @param resolution_grid Numeric vector of Leiden resolutions.
#' @param no_of_iterations Passed to clustering.
#' @param partition_type Passed to clustering.
#'   For Leiden, resolution sweeps are only informative with
#'   `"CPMVertexPartition"` or `"RBConfigurationVertexPartition"`.
#' @param max_cluster_count_per_gene Passed to clustering.
#' @param module_count_bounds Optional length-2 numeric vector
#'   (`min_modules`, `max_modules`) used as a soft constraint.
#' @param prefer_resolution Tie-break preference for resolution.
#' @param verbose Logical; print progress messages.
#'
#' @return Updated `HCoCenaExperiment`.
#' @examples
#' hc <- hc_example_data("clustered")
#' hc <- hc_tune_resolution(hc, apply = FALSE)
#' @export
hc_tune_resolution <- function(hc,
                               apply = FALSE,
                               cluster_algorithms = c("cluster_leiden", "cluster_louvain"),
                               resolution_grid = base::seq(0.05, 0.5, by = 0.05),
                               no_of_iterations = 2,
                               partition_type = "RBConfigurationVertexPartition",
                               max_cluster_count_per_gene = 1,
                               module_count_bounds = c(3, 60),
                               prefer_resolution = 0.1,
                               verbose = TRUE) {
  hc <- hc_auto_tune(
    hc = hc,
    tune_cutoff = FALSE,
    tune_clustering = TRUE,
    apply = apply,
    cluster_algorithms = cluster_algorithms,
    resolution_grid = resolution_grid,
    no_of_iterations = no_of_iterations,
    partition_type = partition_type,
    max_cluster_count_per_gene = max_cluster_count_per_gene,
    module_count_bounds = module_count_bounds,
    prefer_resolution = prefer_resolution,
    verbose = verbose
  )

  sat <- as.list(hc@satellite)
  auto_rep <- sat[["auto_tune"]]
  if (is.list(auto_rep) && "clustering" %in% base::names(auto_rep)) {
    sat[["resolution_tuning"]] <- auto_rep[["clustering"]]
    hc@satellite <- S4Vectors::SimpleList(sat)
    methods::validObject(hc)
  }
  hc
}
