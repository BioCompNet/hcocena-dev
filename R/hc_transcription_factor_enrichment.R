#' Transcription Factor Over-representation Analysis
#'
#' Runs a ChEA3 transcription factor (TF) enrichment analysis for every selected gene module.
#' 	It filters the ranked enriched TFs for the top highest ranking ones and their top highest ranking targets.
#' 	MeanRank was chosen as a ranking method.
#' Results are visualized as circular plots, one for each module.
#'  The outer circle highlights TFs (grey) and their targets (white). Within the inner circle, the colour of the cell is colored according to the module the respective gene can be found in.
#' 	If a gene is a target of one of the enriched transcription factors, a dashed link is connecting their cells.
#' 	If there is an edge in the constructed co-expression network connecting those two genes, that link is solid.
#' @param topTF Integer. The number of top ranking TFs to return per cluster. Default is 5.
#' @param topTarget Integer. The number of top ranking targets to return per TF. Default is 5.
#' @param clusters Either "all" (default) or a vector of clusters as strings. Defines for which clusters to perform the analysis.
#' @noRd

.hc_TF_overrep_module_driver <- function(clusters = "all", topTF = 5, topTarget = 5) {
  output <- list()
  gtc <- .hc_gene_to_cluster_impl()
  base::colnames(gtc) <- base::c("gene", "cluster")

  all_clusters <- base::unique(hcobject[["integrated_output"]][["cluster_calc"]][["cluster_information"]][["color"]])
  all_clusters <- base::as.character(all_clusters)
  all_clusters <- all_clusters[!base::is.na(all_clusters) & all_clusters != "white" & all_clusters != ""]

  module_label_map <- hcobject[["integrated_output"]][["cluster_calc"]][["module_label_map"]]
  if (!base::is.null(module_label_map) && base::length(module_label_map) > 0) {
    module_label_map <- base::as.character(module_label_map)
    map_names <- base::names(hcobject[["integrated_output"]][["cluster_calc"]][["module_label_map"]])
    if (!base::is.null(map_names) && base::length(map_names) == base::length(module_label_map)) {
      base::names(module_label_map) <- base::as.character(map_names)
    }
    missing_before <- base::setdiff(all_clusters, base::names(module_label_map))
    if (base::length(missing_before) > 0) {
      inverse_map <- stats::setNames(base::names(module_label_map), base::as.character(module_label_map))
      if (base::all(all_clusters %in% base::names(inverse_map))) {
        module_label_map <- inverse_map
      }
    }
  } else {
    module_label_map <- NULL
  }

  display_by_color <- stats::setNames(all_clusters, all_clusters)
  if (!base::is.null(module_label_map) && base::length(module_label_map) > 0) {
    mapped_labels <- base::as.character(module_label_map[all_clusters])
    valid_labels <- !base::is.na(mapped_labels) & base::nzchar(mapped_labels)
    display_by_color[valid_labels] <- mapped_labels[valid_labels]
  }
  label_to_color <- stats::setNames(all_clusters, base::as.character(display_by_color[all_clusters]))

  if (base::length(clusters) >= 1 && clusters[1] == "all") {
    clusters <- all_clusters
  } else {
    requested <- base::as.character(clusters)
    requested <- requested[!base::is.na(requested) & requested != ""]
    resolved_clusters <- base::character(0)
    unresolved <- base::character(0)
    for (cl in requested) {
      if (cl %in% all_clusters) {
        resolved_clusters <- base::c(resolved_clusters, cl)
      } else if (cl %in% base::names(label_to_color)) {
        resolved_clusters <- base::c(resolved_clusters, base::as.character(label_to_color[[cl]]))
      } else {
        unresolved <- base::c(unresolved, cl)
      }
    }
    if (base::length(unresolved) > 0) {
      warning("Ignoring unknown clusters/modules: ", base::paste(base::unique(unresolved), collapse = ", "))
    }
    clusters <- base::unique(resolved_clusters)
    if (base::length(clusters) == 0) {
      stop("No valid clusters/modules selected for TF overrepresentation.")
    }
  }

  tt_list <- list()
  exp_plot_list <- list()
  module_title_by_color <- stats::setNames(base::character(0), base::character(0))
  TFs <- NULL

  for (c in clusters) {
    module_display <- if (c %in% base::names(display_by_color)) {
      base::as.character(display_by_color[[c]])
    } else {
      c
    }

    genes <- dplyr::filter(hcobject[["integrated_output"]][["cluster_calc"]][["cluster_information"]], color == c) %>%
      dplyr::pull(., "gene_n") %>%
      base::strsplit(., split = ",") %>%
      base::unlist(.)

    url <- "https://maayanlab.cloud/chea3/api/enrich/"
    encode <- "json"
    payload <- list(query_name = "myQuery", gene_set = genes)

    # POST to ChEA3 server
    response <- httr::POST(url = url, body = payload, encode = encode)
    json <- httr::content(response, as = "text")

    # results as list of R dataframes
    results <- jsonlite::fromJSON(json)
    results <- results$`Integrated--meanRank`
    gtc_not_white <- gtc[!gtc$cluster == "white", ]
    results <- dplyr::filter(results, TF %in% gtc_not_white$gene)

    # extract those from meanRank since meanRank scored as best method:
    resultlist <- list()
    for (i in base::seq_len(topTF)) {
      if (i > base::length(results$TF)) {
        next
      } else {
        if (results$Overlapping_Genes[i] == "") {
          next
        } else {
          tf <- results$TF[i]
          overlapping_genes <- results$Overlapping_Genes[i] %>%
            base::strsplit(., split = ",") %>%
            base::unlist(.)
          resultlist[[tf]] <- list(
            TF = tf,
            targets = overlapping_genes[base::seq_len(base::min(topTarget, base::length(overlapping_genes)))]
          )
        }
      }
    }
    output_name <- module_display
    if (output_name %in% base::names(output)) {
      output_name <- base::paste0(module_display, " [", c, "]")
    }
    output[[output_name]] <- resultlist

    # prepare TF-target dataframe for visualization
    if (base::length(resultlist) == 0) {
      next
    }

    tt_df <- NULL
    exp_plot_df <- base::names(resultlist)

    for (x in base::names(resultlist)) {
      TFs <- base::c(TFs, x)
      clt <- dplyr::filter(gtc, gene == x) %>% dplyr::pull(., "cluster")
      tmp_df <- base::data.frame(
        TF = base::rep(x, base::length(resultlist[[x]][["targets"]])),
        Target = resultlist[[x]][["targets"]],
        ClusterTF = base::rep(clt, base::length(resultlist[[x]][["targets"]]))
      )
      tmp_df <- tmp_df[stats::complete.cases(tmp_df), ]
      tt_df <- base::rbind(tt_df, tmp_df)
      exp_plot_df <- base::c(exp_plot_df, resultlist[[x]][["targets"]])
    }

    exp_plot_df <- base::as.data.frame(base::unique(exp_plot_df))
    base::colnames(exp_plot_df) <- module_display

    tt_list[[c]] <- tt_df
    exp_plot_list[[c]] <- exp_plot_df
    module_title_by_color[[c]] <- module_display
  }

  # Save TF enrichment results
  .hc_set_bridge_hcobject_slot(c("integrated_output", "TF_overrep_results"), output)


  # Generate plots

  if (length(tt_list) == 0) {
    stop("No transcription factors found to be enriched for any of the modules.")
  }

  TFs <- base::unique(base::as.character(TFs))
  edgelist <- hcobject[["integrated_output"]][["combined_edgelist"]]
  # separator prevents ("MT","CO1") / ("M","TCO1") style key collisions
  edgelist$merged <- base::paste0(base::as.character(edgelist$V1), "\r", base::as.character(edgelist$V2))
  edgelist$merged2 <- base::paste0(base::as.character(edgelist$V2), "\r", base::as.character(edgelist$V1))

  module_keys <- base::names(tt_list)
  if (base::is.null(module_keys)) {
    module_keys <- base::as.character(clusters[base::seq_along(tt_list)])
  }
  module_keys <- base::as.character(module_keys)
  module_titles <- base::vapply(module_keys, function(module_color) {
    if (module_color %in% base::names(module_title_by_color)) {
      base::as.character(module_title_by_color[[module_color]])
    } else {
      module_color
    }
  }, FUN.VALUE = character(1))

  draw_tf_module_panel <- function(idx) {
    module_color <- module_keys[[idx]]
    module_title <- module_titles[[idx]]
    fromto <- tt_list[[idx]]
    fromto <- dplyr::filter(fromto, !Target %in% TFs)
    fromto <- fromto[stats::complete.cases(fromto), ]
    fromto <- base::unique(fromto)

    NodeToColor <- base::rbind(
      base::data.frame(gene = fromto$TF, color = fromto$ClusterTF),
      base::data.frame(gene = fromto$Target, color = base::rep(module_color, base::nrow(fromto)))
    ) %>%
      base::unique()

    factors <- base::unique(base::as.character(NodeToColor$gene))
    circlize::circos.clear()
    circlize::circos.par(points.overflow.warning = FALSE)
    circlize::circos.initialize(factors, xlim = c(0, 1))
    circlize::circos.trackPlotRegion(
      sectors = factors,
      ylim = c(0, 1),
      track.height = 0.05,
      bg.col = base::ifelse(factors %in% TFs, yes = "grey", no = "white")
    )
    circlize::circos.trackPlotRegion(
      sectors = factors,
      ylim = c(0, 1),
      track.height = 0.05,
      bg.col = base::as.character(NodeToColor$color)
    )

    circlize::circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
      xlim <- circlize::get.cell.meta.data("xlim")
      ylim <- circlize::get.cell.meta.data("ylim")
      sector.name <- circlize::get.cell.meta.data("sector.index")
      if (sector.name %in% TFs) {
        circlize::circos.text(base::mean(xlim), base::mean(ylim) + 2.5, sector.name, facing = "inside", niceFacing = TRUE, cex = .9, font = 2)
      } else {
        circlize::circos.text(base::mean(xlim), base::mean(ylim) + 2.5, sector.name, facing = "inside", niceFacing = TRUE, cex = .9)
      }
    })

    for (i in base::seq_len(base::nrow(fromto))) {
      merged <- base::paste0(base::as.character(fromto[i, 1]), "\r", base::as.character(fromto[i, 2]))
      if (merged %in% edgelist$merged | merged %in% edgelist$merged2) {
        circlize::circos.link(
          sector.index1 = base::as.character(fromto[i, 1]), c(0.5),
          sector.index2 = base::as.character(fromto[i, 2]), c(0.5),
          col = base::as.character(fromto[i, 3]),
          lwd = 2, lty = 1,
          directional = 1,
          arr.width = .25,
          arr.length = .25
        )
      } else {
        circlize::circos.link(
          sector.index1 = base::as.character(fromto[i, 1]), c(0.5),
          sector.index2 = base::as.character(fromto[i, 2]), c(0.5),
          col = base::as.character(fromto[i, 3]),
          lwd = 1, lty = 5,
          directional = 1,
          arr.width = .15,
          arr.length = .15
        )
      }
    }
    graphics::title(module_title)
    circlize::circos.clear()
  }

  page_labels <- base::unlist(base::lapply(base::seq_along(module_keys), function(idx) {
    base::c(
      base::paste0(module_titles[[idx]], "_tf"),
      base::paste0(module_titles[[idx]], "_expression")
    )
  }), use.names = FALSE)

  .hc_export_multi_page_plot(
    file = .hc_output_file("TF_overrep_module.pdf"),
    page_labels = page_labels,
    width = 15,
    height = 8,
    display = TRUE,
    draw_page_fun = function(page_idx, page_label) {
      module_idx <- ((page_idx - 1L) %/% 2L) + 1L
      if ((page_idx %% 2L) == 1L) {
        draw_tf_module_panel(module_idx)
      } else {
        .hc_visualize_gene_expression_driver(
          genes = exp_plot_list[[module_idx]] %>% base::unlist(use.names = FALSE),
          name = module_titles[[module_idx]],
          save = FALSE
        )
      }
    }
  )
}


#' Network-Wide Transcription Factor Over-representation Analysis
#'
#' Returns the transcription factors that have the most enriched targets network-wide, including their targets.
#' @param topTF The number of transcription factors with the highest number of enriched targets in the network. Default is 100.
#' @param topTarget Per transcription factor the number of top most enriched targets to return. Default is 30.
#' @noRd


.hc_TF_overrep_network_driver <- function(topTF = 100, topTarget = 30) {
  gtc <- .hc_gene_to_cluster_impl()
  base::colnames(gtc) <- base::c("gene", "cluster")
  genes <- gtc$gene

  # The amp.pharm.mssm.edu host was retired; ChEA3 lives at maayanlab.cloud
  # (the same endpoint .hc_TF_overrep_module_driver() already uses).
  url <- "https://maayanlab.cloud/chea3/api/enrich/"
  encode <- "json"
  payload <- list(query_name = "myQuery", gene_set = genes)

  # POST to ChEA3 server
  response <- httr::POST(url = url, body = payload, encode = encode)
  httr::stop_for_status(response, task = "query the ChEA3 API")
  json <- httr::content(response, as = "text")

  # results as list of R dataframes
  results <- jsonlite::fromJSON(json)
  results <- results$`Integrated--meanRank`
  gtc_not_white <- gtc[!gtc$cluster == "white", ]
  results <- dplyr::filter(results, TF %in% gtc_not_white$gene)

  # extract those from meanRank since meanRank scored as best method:
  resultlist <- list()
  for (i in base::seq_len(topTF)) {
    if (i > length(results$TF)) {
      break
    }
    tf <- results$TF[i]
    overlapping_genes <- results$Overlapping_Genes[i] %>%
      base::strsplit(., split = ",") %>%
      base::unlist(.)

    # genes to which this TF has an edge:
    edgelist <- dplyr::filter(hcobject[["integrated_output"]][["combined_edgelist"]], V1 == tf | V2 == tf)
    edgelist <- base::c(edgelist[, 1] %>% base::as.character(), edgelist[, 2] %>% base::as.character())
    # drop the TF itself from its own neighbour list; this used to compare
    # against the loop counter `i` instead of the TF name.
    edgelist <- edgelist[!edgelist == base::as.character(tf)]

    message("the transcription factor ", tf, " has ", length(overlapping_genes), " targets. It has a co-expression above the cutoff with ", length(overlapping_genes[overlapping_genes %in% edgelist]), " of these targets. The others will be discarded.")
    overlapping_genes <- overlapping_genes[overlapping_genes %in% edgelist]

    if (base::length(overlapping_genes) > topTarget) {
      resultlist[[tf]] <- list(
        TF = tf,
        targets = overlapping_genes[base::seq_len(topTarget)]
      )
    } else {
      resultlist[[tf]] <- list(TF = tf, targets = overlapping_genes)
    }
  }
  .hc_set_bridge_hcobject_slot(c("integrated_output", "enrichall"), resultlist)
  # `integrated_output$enrichall` is not carried back into the S4 object, so
  # keep a copy where it survives; hc_check_tf() reads it from there.
  .hc_set_bridge_hcobject_slot(c("satellite_outputs", "tf_network_targets"), resultlist)
}



