#' Plot module x cell-type matrix from module cell-type annotation
#'
#' Uses `hc@satellite$celltype_annotation` and draws a heatmap-like matrix
#' (modules x cell types) for easier interpretation than stacked barplots.
#'
#' @param hc A `HCoCenaExperiment`.
#' @param database Optional single database name. If `NULL`, uses combined
#'   `selected_celltypes` across databases.
#' @param value One of `"pct"` or `"count"`.
#' @param top_celltypes Maximum number of cell types to show. Remaining cell
#'   types can be collapsed into `other_label` when `include_other = TRUE`.
#' @param include_other Logical. Collapse non-top cell types into one column.
#' @param other_label Label for the collapsed column.
#' @param cluster_order Optional explicit module order (character vector).
#' @param low_color Low-end fill color.
#' @param high_color High-end fill color.
#' @param return_data Logical. If `TRUE`, return a list with `plot`, `matrix`,
#'   and `long_data`; otherwise return only the plot.
#' @return A `ggplot` object or a list containing plot + data.
#' @examples
#' hc <- hc_example_data("clustered")
#' gmt <- system.file("extdata", "toy_celltype_markers.gmt", package = "hcocena")
#' hc <- hc_celltype_annotation(
#'   hc,
#'   databases = character(0),
#'   custom_gmt_files = c(CellTypes = gmt),
#'   export_excel = FALSE
#' )
#' p <- hc_plot_celltype_annotation_matrix(hc)
#' @export
hc_plot_celltype_annotation_matrix <- function(hc,
                                               database = NULL,
                                               value = c("pct", "count"),
                                               top_celltypes = 15,
                                               include_other = TRUE,
                                               other_label = "Other",
                                               cluster_order = NULL,
                                               low_color = "#f7fbff",
                                               high_color = "#08306b",
                                               return_data = TRUE) {
  if (!inherits(hc, "HCoCenaExperiment")) {
    stop("`hc` must be a `HCoCenaExperiment`.")
  }
  value <- base::match.arg(value)
  if (!base::is.null(database)) {
    if (!base::is.character(database) || base::length(database) != 1 || base::is.na(database) || !base::nzchar(base::trimws(database))) {
      stop("`database` must be NULL or a single non-empty character string.")
    }
    database <- base::trimws(database)
  }
  if (!base::is.numeric(top_celltypes) || base::length(top_celltypes) != 1 || base::is.na(top_celltypes) || top_celltypes < 1) {
    stop("`top_celltypes` must be a single integer >= 1.")
  }
  top_celltypes <- base::as.integer(top_celltypes)
  if (!base::is.logical(include_other) || base::length(include_other) != 1 || base::is.na(include_other)) {
    stop("`include_other` must be TRUE or FALSE.")
  }
  if (!base::is.character(other_label) || base::length(other_label) != 1 || base::is.na(other_label) || !base::nzchar(base::trimws(other_label))) {
    stop("`other_label` must be a non-empty character scalar.")
  }
  other_label <- base::trimws(other_label)

  sat <- hc@satellite
  if (!("celltype_annotation" %in% base::names(sat)) || base::is.null(sat[["celltype_annotation"]])) {
    stop("No `celltype_annotation` results found in `hc@satellite`. Run `hc_celltype_annotation()` first.")
  }
  ct <- sat[["celltype_annotation"]]
  if (!base::is.list(ct)) {
    stop("`hc@satellite$celltype_annotation` is not a valid list.")
  }

  if (base::is.null(database)) {
    df <- ct[["selected_celltypes"]]
  } else {
    by_db <- ct[["selected_celltypes_by_db"]]
    if (base::is.null(by_db) || !(database %in% base::names(by_db))) {
      stop("Database `", database, "` not found in `selected_celltypes_by_db`.")
    }
    df <- by_db[[database]]
  }

  if (base::is.null(df) || !base::is.data.frame(df) || base::nrow(df) == 0) {
    stop("No cell-type rows available for the selected database.")
  }
  if (!all(c("cluster", "cell_type") %in% base::colnames(df))) {
    stop("Selected cell-type table must contain `cluster` and `cell_type` columns.")
  }

  value_col <- if (value == "pct" && "pct" %in% base::colnames(df)) {
    "pct"
  } else if ("count" %in% base::colnames(df)) {
    "count"
  } else {
    stop("Requested value column not available (`pct`/`count`).")
  }

  agg <- stats::aggregate(
    df[[value_col]],
    by = list(cluster = base::as.character(df$cluster), cell_type = base::as.character(df$cell_type)),
    FUN = function(x) base::sum(.hc_as_numeric_safely(x), na.rm = TRUE)
  )
  base::colnames(agg)[3] <- "value"
  agg <- agg[!base::is.na(agg$cluster) & base::nzchar(agg$cluster) &
    !base::is.na(agg$cell_type) & base::nzchar(agg$cell_type), , drop = FALSE]
  if (base::nrow(agg) == 0) {
    stop("No non-empty cluster/cell_type combinations available for plotting.")
  }

  mat <- stats::xtabs(value ~ cluster + cell_type, data = agg)
  mat <- base::as.matrix(mat)

  all_clusters <- base::rownames(mat)
  cluster_calc <- tryCatch(hc@integration@cluster, error = function(e) NULL)
  module_prefix <- "M"
  if (!base::is.null(cluster_calc) && "module_prefix" %in% base::names(cluster_calc)) {
    tmp_prefix <- base::as.character(cluster_calc[["module_prefix"]])
    if (base::length(tmp_prefix) > 0 && !base::is.na(tmp_prefix[[1]]) && base::nzchar(tmp_prefix[[1]])) {
      module_prefix <- tmp_prefix[[1]]
    }
  }
  module_label_map <- if (!base::is.null(cluster_calc) && "module_label_map" %in% base::names(cluster_calc)) {
    cluster_calc[["module_label_map"]]
  } else {
    NULL
  }
  if (!base::is.null(module_label_map) && base::length(module_label_map) > 0) {
    map_names <- base::names(module_label_map)
    module_label_map <- base::as.character(module_label_map)
    if (!base::is.null(map_names) && base::length(map_names) == base::length(module_label_map)) {
      base::names(module_label_map) <- base::as.character(map_names)
    }
  }
  if (base::is.null(module_label_map) || base::length(module_label_map) == 0) {
    module_label_map <- stats::setNames(base::paste0(module_prefix, base::seq_along(all_clusters)), all_clusters)
  }
  missing_map <- base::setdiff(all_clusters, base::names(module_label_map))
  if (base::length(missing_map) > 0) {
    module_label_map <- base::c(
      module_label_map,
      stats::setNames(
        base::paste0(module_prefix, base::seq.int(base::length(module_label_map) + 1, length.out = base::length(missing_map))),
        missing_map
      )
    )
  }
  label_to_cluster <- stats::setNames(base::names(module_label_map), base::as.character(module_label_map))

  contrast_text_color <- function(cols) {
    base::vapply(
      cols,
      function(col) {
        rgb_mat <- tryCatch(grDevices::col2rgb(col) / 255, error = function(e) NULL)
        if (base::is.null(rgb_mat)) {
          return("black")
        }
        luminance <- (0.2126 * rgb_mat[1, 1]) + (0.7152 * rgb_mat[2, 1]) + (0.0722 * rgb_mat[3, 1])
        if (luminance < 0.5) "white" else "black"
      },
      FUN.VALUE = base::character(1)
    )
  }

  if (base::ncol(mat) > top_celltypes) {
    totals <- base::colSums(mat, na.rm = TRUE)
    keep <- base::names(base::sort(totals, decreasing = TRUE))[base::seq_len(top_celltypes)]
    remainder <- base::setdiff(base::colnames(mat), keep)
    mat <- mat[, keep, drop = FALSE]
    if (isTRUE(include_other) && base::length(remainder) > 0) {
      mat <- base::cbind(mat, .other = base::rowSums(stats::xtabs(value ~ cluster + cell_type, data = agg)[, remainder, drop = FALSE], na.rm = TRUE))
      base::colnames(mat)[base::ncol(mat)] <- other_label
    }
  }

  if (base::is.null(cluster_order)) {
    cl <- base::as.character(base::rownames(mat))
    module_labels_default <- base::as.character(module_label_map[cl])
    module_labels_default[base::is.na(module_labels_default) | module_labels_default == ""] <- cl[base::is.na(module_labels_default) | module_labels_default == ""]
    mnum <- .hc_as_numeric_safely(sub("^M([0-9]+).*", "\\1", module_labels_default))
    ord <- base::order(base::is.na(mnum), mnum, cl)
    mat <- mat[ord, , drop = FALSE]
  } else {
    cluster_order <- base::as.character(cluster_order)
    cluster_order_resolved <- ifelse(
      cluster_order %in% base::rownames(mat),
      cluster_order,
      ifelse(cluster_order %in% base::names(label_to_cluster), base::as.character(label_to_cluster[cluster_order]), NA_character_)
    )
    present <- cluster_order_resolved[!base::is.na(cluster_order_resolved) & cluster_order_resolved %in% base::rownames(mat)]
    missing <- base::setdiff(base::rownames(mat), present)
    mat <- mat[base::c(present, missing), , drop = FALSE]
  }

  cluster_colors_current <- base::as.character(base::rownames(mat))
  module_labels <- base::as.character(module_label_map[cluster_colors_current])
  module_labels[base::is.na(module_labels) | module_labels == ""] <- base::rownames(mat)[base::is.na(module_labels) | module_labels == ""]
  base::rownames(mat) <- module_labels

  long_df <- base::as.data.frame(as.table(mat), stringsAsFactors = FALSE)
  base::colnames(long_df) <- c("module", "cell_type", "value")
  long_df$module_color <- stats::setNames(cluster_colors_current, module_labels)[long_df$module]
  long_df$label_color <- contrast_text_color(long_df$module_color)
  x_levels <- base::c(".module", base::colnames(mat))
  long_df$x_key <- factor(base::as.character(long_df$cell_type), levels = x_levels)
  long_df$module <- factor(long_df$module, levels = base::rev(base::rownames(mat)))
  long_df$cell_type <- factor(long_df$cell_type, levels = base::colnames(mat))

  fill_title <- if (value_col == "pct") "Percent" else "Count"
  title_txt <- if (base::is.null(database)) {
    "Module x cell-type matrix (combined DBs)"
  } else {
    base::paste0("Module x cell-type matrix (", database, ")")
  }

  module_strip_df <- base::unique(long_df[, c("module", "module_color", "label_color"), drop = FALSE])
  module_strip_df$x_key <- factor(".module", levels = x_levels)

  p <- ggplot2::ggplot(long_df, ggplot2::aes(x = x_key, y = module, fill = value)) +
    ggplot2::geom_tile(
      data = module_strip_df,
      mapping = ggplot2::aes(x = x_key, y = module),
      inherit.aes = FALSE,
      fill = module_strip_df$module_color,
      color = "white",
      linewidth = 0.25
    ) +
    ggplot2::geom_text(
      data = module_strip_df,
      mapping = ggplot2::aes(x = x_key, y = module, label = module),
      inherit.aes = FALSE,
      color = module_strip_df$label_color,
      fontface = "bold",
      size = 3.4
    ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.25) +
    ggplot2::scale_fill_gradient(low = low_color, high = high_color, name = fill_title) +
    ggplot2::scale_x_discrete(labels = stats::setNames(base::c("", base::colnames(mat)), x_levels), drop = FALSE) +
    ggplot2::labs(title = title_txt, x = "Cell type", y = "Module") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 55, hjust = 1, vjust = 1),
      plot.title = ggplot2::element_text(face = "bold")
    )

  if (isTRUE(return_data)) {
    return(list(plot = p, matrix = mat, long_data = long_df))
  }
  p
}
