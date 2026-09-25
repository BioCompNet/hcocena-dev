#' Numeric Metadata Correaltion
#'
#' Calculates a matrix where columns are names of the meta categories and rows are modules Cells contain the Pearson correlation value between
#' 		a) the mean expressions of cluster genes in each sample with
#' 		b) the numeric meta value in each sample.
#' @param set An integer. The number of the datset for which to perform the calculation.
#' @param meta A vector of strings. The names of the numeric annotation column(s) which to correlate to the cluster expression patterns.
#' @param p_val The maximum p-value to determine a correlation as significant. Default is 0.05. Non-significant correlations are shown in grey.
#' @param padj Method to use for multiple testing correction. Can be one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".  Default is "BH" (Benjamini-Hochberg).
#' @noRd


.hc_meta_correlation_num_driver <- function(set, meta, p_val = 0.05, padj = "BH") {
  meta_data <- dplyr::select(hcobject[["data"]][[base::paste0("set", set, "_anno")]], dplyr::all_of(meta))
  counts <- sample_wise_cluster_expression(set = set)
  cors <- base::lapply(base::colnames(meta_data), function(x) {
    meta_vals <- dplyr::pull(meta_data, x) %>%
      base::as.character() %>%
      base::as.numeric()
    correlation <- base::lapply(base::rownames(counts), function(y) {
      mean_vals <- counts[y, ] %>%
        base::t() %>%
        base::as.data.frame() %>%
        dplyr::pull(., y)
      cor_and_pval <- stats::cor.test(mean_vals, meta_vals, method = "pearson")
      return(base::data.frame(group_x = x, group_y = y, r = (cor_and_pval$estimate %>% base::round(., digits = 2)), p = cor_and_pval$p.value))
    }) %>%
      rlist::list.rbind() %>%
      base::as.data.frame()
    return(correlation)
  }) %>% rlist::list.rbind()

  cors$p_adj <- stats::p.adjust(cors$p, method = padj)

  cors$p_adj <- .hc_as_numeric_safely(cors$p_adj)
  cors$r <- .hc_as_numeric_safely(cors$r)
  cors$pearson_corr <- base::ifelse(
    base::is.finite(cors$p_adj) & cors$p_adj <= p_val,
    cors$r,
    NA_real_
  )

  heatmap_info <- .hc_heatmap_cache_info(hcobject[["integrated_output"]][["cluster_calc"]])
  row_levels <- heatmap_info$row_order
  if (is.null(row_levels) || length(row_levels) == 0) {
    row_levels <- unique(as.character(cors$group_y))
  }
  row_levels <- c(
    row_levels[row_levels %in% unique(as.character(cors$group_y))],
    setdiff(unique(as.character(cors$group_y)), row_levels)
  )


  g <- ggplot2::ggplot(data = cors, ggplot2::aes(
    x = group_x,
    y = base::factor(group_y, levels = base::rev(row_levels)),
    fill = pearson_corr
  )) +
    ggplot2::geom_tile(color = "black") +
    ggplot2::scale_fill_gradientn(colours = grDevices::colorRampPalette(base::rev(RColorBrewer::brewer.pal(n = 7, name = "BrBG")))(21), limits = c(-1, 1)) +
    ggplot2::geom_text(ggplot2::aes(group_x, group_y, label = r), color = "black", size = 4) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 90)
    ) +
    ggplot2::scale_y_discrete(position = "right") +
    ggplot2::ggtitle(paste0("Correlation of ", meta, " with modules"))

  graphics::plot(g)
  .hc_export_ggplot_file(
    file = .hc_output_file(base::paste0("numerical_meta_correlation_", hcobject[["layers_names"]][set], ".pdf")),
    plot = g,
    width = 10,
    height = 8
  )
}


#' Categorical Metadata Correaltion
#'
#' Calculates a matrix where columns are the different values of the meta category (e.g. if meta = survival and "survival can have values "yes" or "no", then "yes" and "no" will be the columns) and rows are modues.
#' 	Cells contain the Pearson correlation value between
#' 		a) the mean expressions of cluster genes in each voi group with
#' 		b) the counts of the meta value (e.g., "yes" or "no") across voi groups.
#' @param set An integer. The number of the dataset for which to perform the calculation.
#' @param meta A single string. The name of the categorical annotation column which to correlate to the cluster expression patterns.
#' @param p_val The maximum adjusted p-value to determine a correlation as significant. Default is 0.05. Non-significant correlations are shown in grey.
#' @param padj Method to use for multiple testing correction. Can be one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".  Default is "BH" (Benjamini-Hochberg).
#' @noRd

.hc_meta_correlation_cat_driver <- function(meta, set, p_val = 0.05, padj = "BH") {
  set_name <- base::paste0("set", set)
  gfc_all_genes <- hcobject[["layer_specific_outputs"]][[set_name]][["part2"]][["GFC_all_genes"]]

  # Correlation of the pattern the metainfo has across groups with the pattern
  # of the modules across groups. Both vectors have one entry per group, so
  # Pearson needs at least three groups to be defined at all - with two the
  # correlation is degenerate and cor.test() aborts with a cryptic message.
  groups <- base::setdiff(base::colnames(gfc_all_genes), "Gene")
  if (base::length(groups) < 3) {
    stop(
      "At least three groups of `", hcobject[["global_settings"]][["voi"]],
      "` are needed to correlate a categorical variable across groups; found: ",
      base::paste(groups, collapse = ", ")
    )
  }
  # check length of meta
  if (base::length(meta) > 1) stop("Assign only one variable to the parameter 'meta'!")

  #  extract anno of given data set:
  anno <- hcobject[["data"]][[base::paste0("set", set, "_anno")]]

  # for convenience: data frame of 3 columns - one containing sample IDs, one with corresponding voi-value and one with corresponding value of the meta information of interest
  df <- base::data.frame(
    sample = base::rownames(anno),
    voi = anno[hcobject[["global_settings"]][["voi"]]],
    meta = anno[meta]
  )
  base::colnames(df) <- base::c("sample", "voi", "meta")


  # create data frame with moduels as rows and samples as columns. The value of a cell [i,j] is the mean expression of all genes in cluster i in sample j:
  cluster_X_sample_counts <- sample_wise_cluster_expression(set = set)

  # Group the samples into their voi groups to correspond to the module heatmap.
  # Grouping of values happens by calculating their mean:
  cluster_X_sample_mean_counts <- base::lapply(base::unique(df$voi), function(x) {
    samples <- dplyr::filter(df, voi == x) %>% dplyr::pull(., "sample")
    tmp <- dplyr::select(cluster_X_sample_counts, dplyr::all_of(samples)) %>% base::apply(., 1, base::mean)
    tmp <- base::data.frame(v1 = tmp)
    base::colnames(tmp) <- x
    return(tmp)
  }) %>% rlist::list.cbind()

  # For each group of samples (voi) store the values of the meta info of interest for these samples:
  voi_list <- base::lapply(base::unique(df$voi), function(x) {
    tmp <- dplyr::filter(df, voi == x)
    return(tmp["meta"])
  })

  base::names(voi_list) <- base::unique(df$voi)

  # For each value that the meta information can take (e.g., recovered = YES or NO) count its occurences per sample group:
  out <- base::lapply(base::unique(df$meta), function(x) {
    counts <- base::lapply(voi_list, function(y) {
      vec <- y[, 1]
      return(base::length(vec[vec == x]))
    }) %>% base::unlist()
    out <- base::data.frame(counts = counts)
    base::rownames(out) <- base::names(voi_list)
    base::colnames(out) <- base::c(base::as.character(x))

    return(out)
  })
  out <- rlist::list.cbind(out) %>%
    base::t() %>%
    base::as.data.frame()
  out <- dplyr::select(out, dplyr::all_of(base::colnames(cluster_X_sample_mean_counts)))


  # Transform absolute occurrence counts to percentages:
  out <- base::apply(out, 1, function(x) {
    if (base::sum(x) == 0) {
      base::rep(0, base::length(x))
    } else {
      x / base::sum(x)
    }
  }) %>%
    base::t() %>%
    base::as.data.frame()


  # Calculate pearson correlation of those count fractions across groups with cluster mean expressions across groups:
  cors <- base::lapply(base::seq_len(base::nrow(out)), function(x) {
    vec1 <- base::as.numeric(out[x, ])
    tmp <- base::lapply(base::seq_len(base::nrow(cluster_X_sample_mean_counts)), function(y) {
      vec2 <- base::as.numeric(cluster_X_sample_mean_counts[y, ])
      cor_res <- stats::cor.test(x = vec1, y = vec2, method = "pearson")
      return(base::data.frame(r = cor_res$estimate, p = cor_res$p.value))
    }) %>% rlist::list.rbind()
    base::rownames(tmp) <- base::rownames(cluster_X_sample_mean_counts)
    return(tmp)
  })
  base::names(cors) <- base::rownames(out)


  # Prepare for plotting:
  vals <- base::lapply(base::names(cors), function(x) {
    tmp <- base::cbind(base::data.frame(group_x = base::rep(x, base::nrow(cors[[x]])), group_y = base::rownames(cors[[x]])), cors[[x]])
    return(tmp)
  }) %>% rlist::list.rbind()

  vals$p_adj <- stats::p.adjust(vals$p, method = padj)

  # NB: this used to go through apply() over the data frame, which coerces every
  # row to character, so `p_adj > p_val` became a *string* comparison:
  # "1e-08" > "0.05" is TRUE, i.e. the most significant correlations were the
  # ones being discarded. Compare numerically.
  vals$p_adj <- .hc_as_numeric_safely(vals$p_adj)
  vals$r <- .hc_as_numeric_safely(vals$r)
  vals$pearson_corr <- base::ifelse(
    base::is.finite(vals$p_adj) & vals$p_adj <= p_val,
    vals$r,
    NA_real_
  )

  vals["r"] <- base::round(vals["r"], digits = 2)

  heatmap_info <- .hc_heatmap_cache_info(hcobject[["integrated_output"]][["cluster_calc"]])
  row_levels <- heatmap_info$row_order
  if (is.null(row_levels) || length(row_levels) == 0) {
    row_levels <- unique(as.character(vals$group_y))
  }
  row_levels <- c(
    row_levels[row_levels %in% unique(as.character(vals$group_y))],
    setdiff(unique(as.character(vals$group_y)), row_levels)
  )


  # plot
  g <- ggplot2::ggplot(data = vals, ggplot2::aes(
    x = group_x,
    y = base::factor(group_y, levels = base::rev(row_levels)),
    fill = pearson_corr
  )) +
    ggplot2::geom_tile(color = "black") +
    ggplot2::scale_fill_gradientn(colours = grDevices::colorRampPalette(base::rev(RColorBrewer::brewer.pal(n = 7, name = "BrBG")))(21), limits = c(-1, 1)) +
    ggplot2::geom_text(ggplot2::aes(group_x, group_y, label = r), color = "black", size = 4) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank()
    ) +
    ggplot2::scale_y_discrete(position = "right") +
    ggplot2::ggtitle(paste0("Correlation of ", meta, " with modules"))

  graphics::plot(g)
}




#' Correlate numeric metadata with modules (S4 API)
#'
#' Correlates the mean expression of each module's genes per sample with one or
#' more numeric annotation columns, and shows the result as a heatmap.
#' Non-significant correlations are left blank.
#'
#' @param hc A `HCoCenaExperiment`.
#' @param set An integer. Number of the dataset/layer to use.
#' @param meta Character vector of numeric annotation column names to correlate
#'   with the module expression patterns.
#' @param p_val Maximum adjusted p-value for a correlation to count as
#'   significant. Default is 0.05.
#' @param padj Multiple-testing correction method passed to
#'   [stats::p.adjust()]. Default is `"BH"`.
#' @return Updated `HCoCenaExperiment`.
#' @examples
#' hc <- hc_example_data("clustered")
#' hc_meta_correlation_num(hc, set = 1, meta = "age")
#' @export
hc_meta_correlation_num <- function(hc, set, meta, p_val = 0.05,
                                    padj = "BH") {
  .hc_run_driver(
    hc = hc, fun = .hc_meta_correlation_num_driver,
    set = set, meta = meta, p_val = p_val, padj = padj
  )
}
