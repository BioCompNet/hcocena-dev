#' Import Clusters From File
#'
#' Uses an imported clustering model instead of clustering the integrated network.
#' 	The model must be saved as two columns, the first containing genes, the second cluster colors (this is the format .hc_export_clusters_driver() exports to).
#' @param file File path.
#' @param sep The separator of the file. Default is tab-separated.
#' @param header A Boolean. Whether or not the file has headers (column names).
#' @noRd

.hc_import_clusters_driver <- function(file, sep = "\t", header = TRUE) {
  gtc <- readr::read_delim(file = file, delim = sep, col_names = header)
  base::colnames(gtc) <- base::c("gene", "cluster")
  gtc[] <- base::lapply(gtc, base::as.character)

  new_cluster_info <- NULL
  for (c in base::unique(gtc$cluster)) {
    tmp <- gtc[gtc$cluster == c, ]
    gene_no <- base::nrow(tmp)
    gene_n <- base::paste0(tmp$gene, collapse = ",")
    if (c == "white") {
      cluster_included <- "no"
      vertexsize <- 1
    } else {
      cluster_included <- "yes"
      vertexsize <- 3
    }
    color <- c
    conditions <- base::paste0(
      .hc_gfc_condition_names(hcobject[["integrated_output"]][["GFC_all_layers"]]),
      collapse = "#"
    )
    gfc_means <- .hc_gfc_colmeans_for_genes(
      hcobject[["integrated_output"]][["GFC_all_layers"]],
      genes = tmp$gene
    )
    grp_means <- base::paste0(base::round(gfc_means, 3), collapse = ",")

    new_cluster_info <- base::rbind(
      new_cluster_info,
      base::data.frame(
        clusters = c,
        gene_no = gene_no,
        gene_n = gene_n,
        cluster_included = cluster_included,
        color = c,
        conditions = conditions,
        grp_means = grp_means,
        vertexsize = vertexsize,
        stringsAsFactors = FALSE
      )
    )
  }
  .hc_set_bridge_hcobject_slot(c("integrated_output", "cluster_calc", "cluster_information"), new_cluster_info)
}


#' Import an external gene-to-cluster assignment
#'
#' Replaces the current module assignment with one read from a file, e.g. a
#' clustering exported earlier with [hc_export_clusters()].
#'
#' @param hc A `HCoCenaExperiment`.
#' @param file Path to the gene-to-cluster file. Column 1 holds gene symbols,
#'   column 2 the cluster each gene belongs to.
#' @param sep Field separator of the file. Default is `"\t"`.
#' @param header Logical. Whether the file carries a header row. Default `TRUE`.
#' @return Updated `HCoCenaExperiment`.
#' @examples
#' hc <- hc_example_data("clustered")
#' gtc <- hc_gene_to_cluster(hc)
#' f <- tempfile(fileext = ".tsv")
#' utils::write.table(gtc, f, sep = "\t", row.names = FALSE, quote = FALSE)
#' hc <- hc_import_clusters(hc, file = f)
#' @export
hc_import_clusters <- function(hc, file, sep = "\t", header = TRUE) {
  .hc_run_driver(
    hc = hc, fun = .hc_import_clusters_driver,
    file = file, sep = sep, header = header
  )
}
