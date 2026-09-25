#' Access `MultiAssayExperiment` from an `HCoCenaExperiment`
#' @param x An `HCoCenaExperiment`.
#' @examples
#' hc <- hc_init()
#' class(hc_mae(hc))
#' @return A `MultiAssayExperiment::MultiAssayExperiment` object.
#' @export
setGeneric("hc_mae", function(x) standardGeneric("hc_mae"))

#' @rdname hc_mae
#' @export
setMethod("hc_mae", "HCoCenaExperiment", function(x) x@mae)

#' Access configuration from an `HCoCenaExperiment`
#' @param x An `HCoCenaExperiment`.
#' @examples
#' hc <- hc_init()
#' class(hc_config(hc))
#' @return A `HCoCenaConfig` object.
#' @export
setGeneric("hc_config", function(x) standardGeneric("hc_config"))

#' @rdname hc_config
#' @export
setMethod("hc_config", "HCoCenaExperiment", function(x) x@config)

#' Access layer results from an `HCoCenaExperiment`
#' @param x An `HCoCenaExperiment`.
#' @examples
#' hc <- hc_init()
#' class(hc_layer_results(hc))
#' @return A `S4Vectors::SimpleList`.
#' @export
setGeneric("hc_layer_results", function(x) standardGeneric("hc_layer_results"))

#' @rdname hc_layer_results
#' @export
setMethod("hc_layer_results", "HCoCenaExperiment", function(x) x@layer_results)

#' Access integration payload from an `HCoCenaExperiment`
#' @param x An `HCoCenaExperiment`.
#' @examples
#' hc <- hc_init()
#' class(hc_integration(hc))
#' @return A `HCoCenaIntegration` object.
#' @export
setGeneric("hc_integration", function(x) standardGeneric("hc_integration"))

#' @rdname hc_integration
#' @export
setMethod("hc_integration", "HCoCenaExperiment", function(x) x@integration)

#' Access cluster information from an `HCoCenaExperiment`
#' @param x An `HCoCenaExperiment`.
#' @examples
#' hc <- hc_init()
#' hc_clusters(hc)
#' @return Cluster information or `NULL` if missing.
#' @export
setGeneric("hc_clusters", function(x) standardGeneric("hc_clusters"))

#' @rdname hc_clusters
#' @export
setMethod("hc_clusters", "HCoCenaExperiment", function(x) {
  if ("cluster_information" %in% base::names(x@integration@cluster)) {
    x@integration@cluster[["cluster_information"]]
  } else {
    NULL
  }
})

#' Access the integrated network of an `HCoCenaExperiment`
#' @param x An `HCoCenaExperiment`.
#' @examples
#' hc <- hc_example_data("clustered")
#' g <- hc_graph(hc)
#' igraph::vcount(g)
#' @return The integrated network as an `igraph` object, or `NULL` if
#'   `hc_build_integrated_network()` has not been run yet.
#' @export
setGeneric("hc_graph", function(x) standardGeneric("hc_graph"))

#' @rdname hc_graph
#' @export
setMethod("hc_graph", "HCoCenaExperiment", function(x) x@integration@graph)

#' Access downstream results stored in an `HCoCenaExperiment`
#'
#' Downstream steps (module gene lists, enrichments, hub genes, module
#' statistics, ...) store their results in a named list inside the object.
#' `hc_satellite()` returns that list, or a single entry of it.
#' @param x An `HCoCenaExperiment`.
#' @param name Optional name of a single entry. If `NULL` (default), the whole
#'   list is returned.
#' @examples
#' hc <- hc_example_data("clustered")
#' names(hc_satellite(hc))
#' @return A `S4Vectors::SimpleList` of stored results, or the entry `name`
#'   (`NULL` if that entry does not exist).
#' @export
setGeneric("hc_satellite", function(x, name = NULL) standardGeneric("hc_satellite"))

#' @rdname hc_satellite
#' @export
setMethod("hc_satellite", "HCoCenaExperiment", function(x, name = NULL) {
  if (base::is.null(name)) {
    return(x@satellite)
  }
  if (!(base::is.character(name) && base::length(name) == 1L)) {
    stop("`name` must be a single character string or NULL.", call. = FALSE)
  }
  if (name %in% base::names(x@satellite)) x@satellite[[name]] else NULL
})
