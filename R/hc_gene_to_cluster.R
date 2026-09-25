#' Gene To Cluster Dictionary
#'
#' The function maps the gene names to their corresponding cluster.
#' @param cluster_information Cluster table, typically
#'   `hcobject[["integrated_output"]][["cluster_calc"]][["cluster_information"]]`.
#' @return A data frame with two columns, the first containing gene names as strings, the second containing cluster colours as strings.
#' @noRd

.hc_gene_to_cluster_impl <- function(cluster_information = hcobject[["integrated_output"]][["cluster_calc"]][["cluster_information"]]) {
  # Nearly every downstream analysis starts here, so this is the right place to
  # notice that there are no modules - after an input or parameter change
  # discarded them, or before clustering has been run at all.
  if (base::is.null(cluster_information) ||
    base::length(base::dim(cluster_information)) != 2L ||
    base::nrow(cluster_information) == 0) {
    stop(
      "No module assignment available. Run `hc_cluster_calculation()` first ",
      "(and `hc_build_integrated_network()` before it). If you changed the ",
      "data, the layer settings or the cutoff, the previous modules were ",
      "discarded because they no longer matched.",
      call. = FALSE
    )
  }
  gtc <- base::do.call(rbind, base::apply(cluster_information, 1, function(x) {
    tmp <- x["gene_n"] %>%
      base::strsplit(., split = ",") %>%
      base::unlist(.)
    base::data.frame(gene = tmp, color = base::rep(x["color"], base::length(tmp)))
  }))
  return(gtc)
}

#' Gene-to-module table (S4 API)
#'
#' Returns the mapping of every network gene to the module it was assigned to.
#'
#' @param hc A `HCoCenaExperiment`.
#' @return A data frame with two columns: `gene` (gene symbol) and `color`
#'   (module colour). Genes that were not assigned to any module carry the
#'   colour `"white"`.
#' @examples
#' hc <- hc_example_data("clustered")
#' head(hc_gene_to_cluster(hc))
#' @export
hc_gene_to_cluster <- function(hc) {
  if (!inherits(hc, "HCoCenaExperiment")) {
    stop("`hc` must be a `HCoCenaExperiment`.")
  }
  cluster_info <- as.list(hc@integration@cluster)[["cluster_information"]]
  if (base::is.null(cluster_info) || base::nrow(cluster_info) == 0) {
    stop("No cluster information found. Run `hc_cluster_calculation()` first.")
  }
  .hc_gene_to_cluster_impl(base::as.data.frame(cluster_info, stringsAsFactors = FALSE))
}
