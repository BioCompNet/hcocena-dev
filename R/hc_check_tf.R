#' Check Transcription Factor
#'
#' This function leverages the information collected with hc_tf_overrep_network()
#' 	to allow the user to query specific transcription factors of interest and see how their top targets are spread across modules.
#' 	The goal is to uncover potential co-regulations between clusters.
#' @param TF A string giving the name of the transcription factor to be queried.
#' @noRd

.hc_check_tf_driver <- function(TF) {
  # get edgelist from integrated network:
  edgelist <- hcobject[["integrated_output"]][["combined_edgelist"]][, base::c(1, 2)]
  edgelist[] <- base::lapply(edgelist, as.character)
  edgelist$merged <- base::paste0(edgelist$V1, edgelist$V2)
  edgelist$merged2 <- base::paste0(edgelist$V2, edgelist$V1)

  gtc <- .hc_gene_to_cluster_impl()
  base::colnames(gtc) <- base::c("gene", "cluster")

  # the targets of the transcription factor in question:
  tf_results <- hcobject[["satellite_outputs"]][["tf_network_targets"]] %||%
    hcobject[["integrated_output"]][["enrichall"]]
  if (base::length(tf_results) == 0) {
    stop(
      "No network-wide TF enrichment found. Run `hc_tf_overrep_network()` ",
      "before `hc_check_tf()`.",
      call. = FALSE
    )
  }
  if (!TF %in% base::names(tf_results)) {
    stop(
      "`", TF, "` is not among the TFs returned by `hc_tf_overrep_network()`. ",
      "Available: ", base::paste(base::names(tf_results), collapse = ", "), ".",
      call. = FALSE
    )
  }
  targets <- tf_results[[TF]][["targets"]]

  # create edgelist from TF to it's targets, removing self edges:
  edges <- base::data.frame(from = base::rep(TF, base::length(targets)), to = targets) %>%
    base::unique()
  edges[] <- base::lapply(edges, as.character)
  edges <- edges[!edges$to == TF, ]
  merged <- base::paste0(edges$from, edges$to)
  # if edge exists in network, colour = black, otherwise color = grey:
  edges$color <- base::lapply(merged, function(x) {
    if (x %in% edgelist$merged | x %in% edgelist$merged2) {
      "black"
    } else {
      "grey"
    }
  }) %>% base::unlist()

  # get nodes (TF and targets that are connected in the network) including their cluster colour:
  nodes <- base::data.frame(name = base::unique(base::c(edges$from, edges$to)))
  nodes <- base::merge(nodes, gtc, by.x = "name", by.y = "gene")
  base::colnames(nodes) <- base::c("name", "color")
  nodes <- nodes[base::order(nodes$color), ]

  # create star plot with TF at the center:
  g <- igraph::graph_from_data_frame(d = edges, vertices = nodes$name)
  igraph::V(g)$color <- base::as.character(nodes$color)
  l <- igraph::layout.star(g, center = igraph::V(g)[TF])
  # plot to PDF:
  .hc_export_single_page_plot(
    file = .hc_output_file(base::paste0("TF_", TF, "_starplot.pdf")),
    width = 10,
    height = 10,
    draw_fun = function() {
      igraph::plot.igraph(
        g,
        layout = l,
        edge.arrow.size = 0.5,
        vertex.label.color = "black",
        edge.color = edges$color,
        vertex.label.cex = 0.7,
        vertex.label.font = 2,
        edge.width = 2,
        vertex.frame.color = base::as.character(nodes$color)
      )
    }
  )

  # plot to markdown:
  igraph::plot.igraph(g,
    layout = l, edge.arrow.size = 0.5, vertex.label.color = "black", edge.color = edges$color,
    vertex.label.cex = 0.7, vertex.label.font = 2, edge.width = 2,
    vertex.frame.color = base::as.character(nodes$color)
  )
}


