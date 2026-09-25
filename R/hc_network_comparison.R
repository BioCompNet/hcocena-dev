#' Compare clusters of 2 networks
#'
#' This function calculates and visualizes the Jaccard-Index of all pairs of clusters from two networks.
#' This allows the comparison of the two networks with respect to the clusters they form and how those clusters relate to each other.
#' @param gtc1_path File path to a 'gtc'-file (gene-to-clsuter) of the first network: Essentially just a file with two columns, the first containing gene symbols and the second
#' giving the cluster each gene belongs to. Such a file can be generated during an hCoCena analysis using the function 'hcocena::.hc_export_clusters_driver()', but as long as
#' the described structure is preserved it can also be generated manually elsewhere.
#' @param gtc2_path See 'gtc1_path', only for network 2.
#' @param sep The separator of the 'gtc'-file. Default is tab-separated.
#' @param header A Boolean. Whether or not the file has headers (column names).
#' @param cellsize The size of the cells/tiles in the plotted heatmap. Default is 18, may be adjusted for aestetics reasons.
#' @return The heatmap object for replotting/re-sizing etc. and the result matrix. Output can be found under hcobject$satellite_outputs$network_comparison_1
#' @noRd

.hc_network_comparison_1_driver <- function(gtc1_path, gtc2_path, sep = "\t", header = TRUE, cellsize = 18) {
  gtc1 <- readr::read_delim(file = gtc1_path, delim = sep, col_names = header)
  base::colnames(gtc1) <- c("gene", "color")
  gtc2 <- readr::read_delim(file = gtc2_path, delim = sep, col_names = header)
  base::colnames(gtc2) <- c("gene", "color")

  out <- base::list()

  for (c1 in unique(gtc1$color)) {
    cdf <- NULL

    # get gene set of this cluster:
    set1 <- dplyr::filter(gtc1, color == c1) %>% dplyr::pull(., "gene")

    # for plot indicate which network this cluster came from:
    cname1 <- base::paste0(c1, " network 1 [", length(set1), "]")

    for (c2 in unique(gtc2$color)) {
      # get gene set of this cluster:
      set2 <- dplyr::filter(gtc2, color == c2) %>% dplyr::pull(., "gene")
      # for plot indicate which network this cluster came from:
      cname2 <- base::paste0(c2, " network 2 [", length(set2), "]")
      # calculate Jaccard Index (JI):
      JI <- calc_jaccard(set1, set2)

      rdf <- base::data.frame(V1 = JI)
      base::colnames(rdf) <- c(cname1)
      base::rownames(rdf) <- c(cname2)
      cdf <- base::rbind(cdf, rdf)
    }

    out[[c1]] <- cdf
  }
  out <- rlist::list.cbind(out)
  p <- pheatmap::pheatmap(as.matrix(out),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    cellheight = cellsize,
    cellwidth = cellsize,
    color = RColorBrewer::brewer.pal(name = "Blues", n = 9),
    display_numbers = TRUE,
    fontsize_number = 6,
    number_color = "orange", main = "Jaccard Index of Cluster Pairs"
  )
  .hc_export_single_page_plot(
    file = .hc_output_file("network_comparison_1.pdf"),
    width = 10,
    height = 7,
    draw_fun = function() {
      .hc_display_object(p)
    }
  )

  .hc_set_bridge_hcobject_slot(c("satellite_outputs", "network_comparison_1"), list(heatmap = p, matrix = as.matrix(out)))
}


#' Calculate Jaccard Index of two sets.
#'
#' Only for internal use.
#' @noRd

calc_jaccard <- function(set1, set2) {
  # get size of the intersection of the 2 sets:
  intersection <- base::intersect(set1, set2) %>% length()
  # get size of the union of the 2 sets:
  union <- base::union(set1, set2) %>% length()
  # return the jaccard index of the 2 sets (|intersection|/|union|):
  return(intersection / union)
}


#' Compare degree distribution of gene set
#'
#' This function accepts two networks and a set of genes. It then calculates the degree-distribution of each gene in both networks and also the Jaccard-Index
#' of each gene's neighbourhoods in the two networks. The results are visualized in a 2D dot plot.
#' This allows the comparison of the two networks with respect to the connectivity of selected genes.
#' @param net1 An igraph object of the first network. hCoCena provides its integrated network as igraph objects. You can find it here: hcobject$integrated_output$merged_net.
#' If you want to create a network manually with external data, please refer to the igraph documentation: https://igraph.org/r/.
#' @param net2 See 'net1', only for network 2.
#' @param as In the future we plan to provide either an igraph object or an edge list. For now, only the igraph option is available.
#' @param gene_vec A vector of gene symbols you wish to investigate. Gene symbols must be provided as strings, e.g., c( "YME1L1", "SLC2A5", "SIAH2", "GPI", "IL10RB").
#' @return The ggplot object for replotting/re-sizing/modification etc. and the data used to plot the ggplot. Output can be found under hcobject$satellite_outputs$network_comparison_2
#' @noRd

.hc_network_comparison_2_driver <- function(net1, net2, as = "igraph", gene_vec) {
  nodes1 <- igraph::V(net1)$name
  nodes1 <- nodes1[nodes1 %in% gene_vec]
  nodes2 <- igraph::V(net2)$name
  nodes2 <- nodes2[nodes2 %in% gene_vec]

  intersection <- base::intersect(nodes1, nodes2)
  if (length(intersection) == 0) {
    base::message("None of the given genes are present in both networks. Aborting.")
    return(NULL)
  } else {
    base::message(length(intersection), " of the given genes are present in both networks.")
  }

  out <- NULL
  for (gene in intersection) {
    # neighbours in net1:
    neighbours1 <- igraph::neighbors(net1, gene, mode = c("all"))$name
    # neighbours in net2:
    neighbours2 <- igraph::neighbors(net2, gene, mode = c("all"))$name
    JI_common <- (base::intersect(neighbours1, neighbours2) %>% length()) / (base::union(neighbours1, neighbours2) %>% length())
    rdf <- base::data.frame(
      gene = gene,
      neighbours_in_1 = log(length(neighbours1), base = 2),
      neighbours_in_2 = log(length(neighbours2), base = 2),
      JI_neighbours = JI_common
    )
    out <- base::rbind(out, rdf)
  }

  g <- ggplot2::ggplot(out, ggplot2::aes(x = neighbours_in_1, y = neighbours_in_2, label = gene)) +
    ggplot2::geom_point(ggplot2::aes(x = neighbours_in_1, y = neighbours_in_2, size = JI_neighbours, color = JI_neighbours)) +
    ggplot2::theme_bw() +
    ggplot2::geom_text(hjust = 0.5, vjust = -1, size = 3) +
    ggplot2::ylim(c(min(c(out$neighbours_in_1, out$neighbours_in_2)), max(c(out$neighbours_in_1, out$neighbours_in_2)))) +
    ggplot2::xlim(c(min(c(out$neighbours_in_1, out$neighbours_in_2)), max(c(out$neighbours_in_1, out$neighbours_in_2)))) +
    ggplot2::ylab("node degree in network 2 (log2)") +
    ggplot2::xlab("node degree in network 1 (log2)") +
    ggplot2::geom_abline() +
    ggplot2::guides(
      color = ggplot2::guide_colorbar(order = 1),
      size = ggplot2::guide_legend(order = 2)
    ) +
    ggplot2::labs(color = "Jaccard-Index of Neighbours", size = ggplot2::element_blank())
  graphics::plot(g)
  .hc_export_ggplot_file(
    file = .hc_output_file("network_comparison_2.pdf"),
    plot = g,
    width = 8,
    height = 6
  )
  .hc_set_bridge_hcobject_slot(c("satellite_outputs", "network_comparison_2"), list(plot = g, data = out))
}


#' Compare two gene-to-cluster assignments
#'
#' Compares two externally stored gene-to-cluster tables and visualises how the
#' modules of the first assignment map onto those of the second.
#'
#' @param hc A `HCoCenaExperiment`.
#' @param gtc1_path Path to the first gene-to-cluster file.
#' @param gtc2_path Path to the second gene-to-cluster file.
#' @param sep Field separator of both files. Default is `"\t"`.
#' @param header Logical. Whether the files carry a header row. Default `TRUE`.
#' @param cellsize Cell size of the resulting heatmap. Default is 18.
#' @return Updated `HCoCenaExperiment`.
#' @examples
#' hc <- hc_example_data("clustered")
#' gtc <- hc_gene_to_cluster(hc)
#' f1 <- tempfile(fileext = ".tsv")
#' f2 <- tempfile(fileext = ".tsv")
#' utils::write.table(gtc, f1, sep = "\t", row.names = FALSE, quote = FALSE)
#' utils::write.table(gtc, f2, sep = "\t", row.names = FALSE, quote = FALSE)
#' hc <- hc_network_comparison_1(hc, gtc1_path = f1, gtc2_path = f2)
#' @export
hc_network_comparison_1 <- function(hc, gtc1_path, gtc2_path, sep = "\t",
                                    header = TRUE, cellsize = 18) {
  .hc_run_driver(
    hc = hc, fun = .hc_network_comparison_1_driver,
    gtc1_path = gtc1_path, gtc2_path = gtc2_path, sep = sep,
    header = header, cellsize = cellsize
  )
}

#' Compare the edges of two networks for a set of genes
#'
#' Contrasts the neighbourhoods of `gene_vec` between two networks and reports
#' which edges are shared and which are unique to either network.
#'
#' @param hc A `HCoCenaExperiment`.
#' @param net1 First network as an `igraph` object.
#' @param net2 Second network as an `igraph` object.
#' @param as Object type of `net1`/`net2`. Default is `"igraph"`.
#' @param gene_vec Character vector of genes to restrict the comparison to.
#' @return Updated `HCoCenaExperiment`.
#' @examples
#' hc <- hc_example_data("clustered")
#' g <- hc_graph(hc)
#' hc <- hc_network_comparison_2(
#'   hc,
#'   net1 = g,
#'   net2 = g,
#'   gene_vec = igraph::V(g)$name[1:3]
#' )
#' @export
hc_network_comparison_2 <- function(hc, net1, net2, as = "igraph", gene_vec) {
  .hc_run_driver(
    hc = hc, fun = .hc_network_comparison_2_driver,
    net1 = net1, net2 = net2, as = as, gene_vec = gene_vec
  )
}
