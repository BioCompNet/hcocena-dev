#' Regroup samples from module-level fold-change profiles
#'
#' Builds a module-by-sample fold-change matrix for one or more layers, clusters
#' samples hierarchically, and cuts each dendrogram into `k` groups. By default,
#' this is a non-destructive preview: the input object, its variable of interest,
#' module assignments, and annotation files are not changed.
#'
#' Fold changes are calculated for every gene in every sample relative to either
#' that gene's mean across the layer or its mean in the configured control
#' samples. Gene-level values are then averaged within the current modules. This
#' follows the signed-ratio convention used by hCoCena: two-fold down is `-2`, no
#' change is `1`, and two-fold up is `2`.
#'
#' @param hc A `HCoCenaExperiment` with clustered modules.
#' @param layer Layer selection. Use `"all"` (default), a layer index, layer id,
#'   configured layer name, or a vector of these values.
#' @param modules Modules used as clustering features. Use `"all"` (default),
#'   current module labels such as `"M2"` or `"M2.1"`, module colors, or
#'   numeric module indices.
#' @param k Number of sample groups after cutting the dendrogram. Supply one
#'   integer for all selected layers or one integer per selected layer. Named
#'   vectors may use layer ids or configured layer names.
#' @param reference Reference used for sample-level fold changes. `"mean"`
#'   compares each sample with the gene mean across all samples in the layer.
#'   `"control"` compares with the mean of samples matching the configured
#'   control label in the current variable-of-interest column.
#' @param method Hierarchical clustering method passed to [stats::hclust()].
#' @param distance Distance used between samples. Supports methods accepted by
#'   [stats::dist()] and `"correlation"` for one minus Pearson correlation.
#' @param group_col Name of the new annotation column created when
#'   `apply = TRUE`.
#' @param apply Logical. If `FALSE` (default), return a preview without changing
#'   `hc`. If `TRUE`, add the inferred groups to `group_col` in each selected
#'   layer and store the diagnostics in `hc@satellite$sample_regrouping`.
#' @param overwrite Logical. Allow an existing `group_col` to be replaced when
#'   `apply = TRUE`. Defaults to `FALSE`.
#' @param show_values Logical. Display rounded fold-change values in heatmap
#'   cells.
#' @param silent Logical. If `FALSE`, draw each generated heatmap. The default
#'   is `!interactive()` so plots are shown in interactive R sessions.
#'
#' @return A list with `hc`, per-layer `results`, `plots`, `parameters`, and
#'   `applied`. Each layer result contains `sample_gfc`, `module_table`,
#'   `distance`, `tree`, `clusters`, and `plot`.
#'
#' @examples
#' hc <- hc_example_data("clustered")
#' preview <- hc_regroup_samples(hc, layer = 1, k = 2, silent = TRUE)
#' preview$results$set1$clusters
#' preview$results$set1$sample_gfc
#'
#' applied <- hc_regroup_samples(
#'   hc,
#'   layer = "all",
#'   k = 2,
#'   group_col = "module_profile_group",
#'   apply = TRUE,
#'   silent = TRUE
#' )
#' hc <- applied$hc
#'
#' @export
hc_regroup_samples <- function(hc,
                               layer = "all",
                               modules = "all",
                               k = 2,
                               reference = c("mean", "control"),
                               method = "complete",
                               distance = "euclidean",
                               group_col = "hc_regroup",
                               apply = FALSE,
                               overwrite = FALSE,
                               show_values = FALSE,
                               silent = !base::interactive()) {
  if (!inherits(hc, "HCoCenaExperiment")) {
    stop("`hc` must be a `HCoCenaExperiment`.")
  }
  reference <- base::match.arg(reference)
  .hc_regroup_validate_scalar_character(method, "method")
  .hc_regroup_validate_scalar_character(distance, "distance")
  .hc_regroup_validate_scalar_character(group_col, "group_col")
  .hc_regroup_validate_flag(apply, "apply")
  .hc_regroup_validate_flag(overwrite, "overwrite")
  .hc_regroup_validate_flag(show_values, "show_values")
  .hc_regroup_validate_flag(silent, "silent")

  method <- base::trimws(method)
  distance <- base::trimws(distance)
  group_col <- base::trimws(group_col)
  valid_methods <- c(
    "ward.D", "ward.D2", "single", "complete", "average", "mcquitty",
    "median", "centroid"
  )
  if (!(method %in% valid_methods)) {
    stop("Unknown `method`: ", method, ". See `stats::hclust()` for supported methods.")
  }
  valid_distances <- c(
    "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski",
    "correlation"
  )
  if (!(distance %in% valid_distances)) {
    stop(
      "Unknown `distance`: ", distance, ". Use a `stats::dist()` method or `correlation`."
    )
  }

  target_layers <- .hc_longitudinal_step1_target_layers(hc = hc, layer = layer)
  k_by_layer <- .hc_regroup_k_by_layer(hc = hc, k = k, layer_ids = target_layers)
  cluster_calc <- .hc_heatmap_view_cluster_calc(hc)
  module_info <- .hc_regroup_resolve_modules(
    modules = modules,
    cluster_calc = cluster_calc
  )
  gfc_colors <- .hc_heatmap_view_resolve_gfc_colors(
    cluster_calc = cluster_calc,
    gfc_colors = NULL
  )
  gfc_limits <- .hc_heatmap_view_resolve_gfc_scale_limits(
    hc = hc,
    cluster_calc = cluster_calc,
    gfc_scale_limits = NULL
  )

  exps <- MultiAssayExperiment::experiments(hc@mae)
  if (base::isTRUE(apply) && !base::isTRUE(overwrite)) {
    conflicts <- base::vapply(
      target_layers,
      function(layer_id) {
        group_col %in% base::colnames(SummarizedExperiment::colData(exps[[layer_id]]))
      },
      FUN.VALUE = base::logical(1)
    )
    if (base::any(conflicts)) {
      stop(
        "Annotation column `", group_col, "` already exists in layer(s): ",
        base::paste(target_layers[conflicts], collapse = ", "),
        ". Use another `group_col` or set `overwrite = TRUE`."
      )
    }
  }

  results <- base::vector("list", base::length(target_layers))
  base::names(results) <- target_layers
  plots <- base::vector("list", base::length(target_layers))
  base::names(plots) <- target_layers

  for (layer_id in target_layers) {
    layer_result <- .hc_regroup_one_layer(
      hc = hc,
      layer_id = layer_id,
      se = exps[[layer_id]],
      module_info = module_info,
      k = k_by_layer[[layer_id]],
      reference = reference,
      method = method,
      distance = distance,
      group_col = group_col,
      gfc_colors = gfc_colors,
      gfc_limits = gfc_limits,
      show_values = show_values
    )
    results[[layer_id]] <- layer_result
    plots[[layer_id]] <- layer_result$plot
    if (!base::isTRUE(silent)) {
      .hc_draw_regroup_heatmap(layer_result$plot)
    }
  }

  out_hc <- hc
  if (base::isTRUE(apply)) {
    for (layer_id in target_layers) {
      se <- MultiAssayExperiment::experiments(out_hc@mae)[[layer_id]]
      anno <- SummarizedExperiment::colData(se)
      assignments <- results[[layer_id]]$clusters
      idx <- base::match(base::rownames(anno), assignments$sample)
      if (base::any(base::is.na(idx))) {
        stop("Internal sample alignment failed while applying regrouping to layer `", layer_id, "`.")
      }
      anno[[group_col]] <- base::as.character(assignments$group[idx])
      SummarizedExperiment::colData(se) <- anno
      out_hc@mae[[layer_id]] <- se
    }

    stored_layers <- base::lapply(results, function(x) {
      x$plot <- NULL
      x
    })
    stored <- list(
      created_at = base::Sys.time(),
      parameters = list(
        layers = target_layers,
        modules = module_info$labels,
        module_colors = module_info$colors,
        k = k_by_layer,
        reference = reference,
        method = method,
        distance = distance,
        group_col = group_col
      ),
      layers = stored_layers
    )
    satellite <- base::as.list(out_hc@satellite)
    satellite[["sample_regrouping"]] <- stored
    out_hc@satellite <- S4Vectors::SimpleList(satellite)
    methods::validObject(out_hc)
  }

  structure(
    list(
      hc = out_hc,
      results = results,
      plots = plots,
      parameters = list(
        layers = target_layers,
        modules = module_info$labels,
        module_colors = module_info$colors,
        k = k_by_layer,
        reference = reference,
        method = method,
        distance = distance,
        group_col = group_col
      ),
      applied = base::isTRUE(apply)
    ),
    class = c("hc_sample_regrouping", "list")
  )
}

#' Plot sample-regrouping heatmaps
#'
#' Draws one or all heatmaps stored in the result returned by
#' [hc_regroup_samples()]. This uses explicit grid drawing so it also works in
#' RStudio and notebook contexts where `print()` may not display a silent
#' `pheatmap` object.
#'
#' @param x An `hc_sample_regrouping` result.
#' @param layer Layer to draw. Use `"all"`, a result layer id, or a numeric
#'   result-layer index.
#' @param ... Reserved for future plotting options.
#'
#' @return Invisibly returns `x`.
#' @method plot hc_sample_regrouping
#' @examples
#' hc <- hc_example_data("clustered")
#' rg <- hc_regroup_samples(hc, layer = 1, k = 2, apply = FALSE)
#' plot(rg)
#' @export
plot.hc_sample_regrouping <- function(x, layer = "all", ...) {
  if (!inherits(x, "hc_sample_regrouping")) {
    stop("`x` must be an `hc_sample_regrouping` result.")
  }
  plot_names <- base::names(x$plots)
  if (base::is.null(plot_names) || base::length(plot_names) == 0) {
    stop("No regrouping plots are stored in `x`.")
  }

  use_all <- base::is.character(layer) && base::length(layer) == 1 &&
    !base::is.na(layer) && identical(base::tolower(base::trimws(layer)), "all")
  if (use_all) {
    selected <- plot_names
  } else if (base::is.numeric(layer)) {
    if (base::any(!base::is.finite(layer)) ||
      base::any(base::abs(layer - base::round(layer)) > 1e-8) ||
      base::any(layer < 1) || base::any(layer > base::length(plot_names))) {
      stop("Numeric `layer` indices are outside the available regrouping plots.")
    }
    selected <- plot_names[base::as.integer(base::round(layer))]
  } else {
    selected <- base::as.character(layer)
    unknown <- base::setdiff(selected, plot_names)
    if (base::length(unknown) > 0) {
      stop(
        "Unknown regrouping plot layer(s): ", base::paste(unknown, collapse = ", "),
        ". Available layers: ", base::paste(plot_names, collapse = ", "), "."
      )
    }
  }

  for (layer_id in selected) {
    .hc_draw_regroup_heatmap(x$plots[[layer_id]])
  }
  invisible(x)
}

.hc_draw_regroup_heatmap <- function(x) {
  gtable <- tryCatch(x$gtable, error = function(e) NULL)
  if (base::is.null(gtable)) {
    stop("Stored regrouping plot does not contain a drawable `gtable`.")
  }
  grid::grid.newpage()
  grid::grid.draw(gtable)
  invisible(x)
}

.hc_regroup_validate_flag <- function(x, name) {
  if (!base::is.logical(x) || base::length(x) != 1 || base::is.na(x)) {
    stop("`", name, "` must be TRUE or FALSE.")
  }
  invisible(TRUE)
}

.hc_regroup_validate_scalar_character <- function(x, name) {
  if (!base::is.character(x) || base::length(x) != 1 || base::is.na(x) ||
    !base::nzchar(base::trimws(x))) {
    stop("`", name, "` must be a non-empty character scalar.")
  }
  invisible(TRUE)
}

.hc_regroup_k_by_layer <- function(hc, k, layer_ids) {
  if (!base::is.numeric(k) || base::length(k) == 0 || base::any(!base::is.finite(k)) ||
    base::any(base::abs(k - base::round(k)) > 1e-8) || base::any(k < 1)) {
    stop("`k` must contain positive finite integers.")
  }
  k <- base::as.integer(base::round(k))
  original_names <- base::names(k)

  if (base::length(k) == 1) {
    out <- base::rep(k, base::length(layer_ids))
    base::names(out) <- layer_ids
    return(out)
  }
  if (!base::is.null(original_names) && base::all(base::nzchar(original_names))) {
    layer_labels <- base::vapply(
      layer_ids,
      function(layer_id) .hc_longitudinal_layer_label(hc = hc, layer_id = layer_id),
      FUN.VALUE = base::character(1)
    )
    idx <- base::match(layer_ids, original_names)
    missing <- base::is.na(idx)
    idx[missing] <- base::match(layer_labels[missing], original_names)
    if (base::any(base::is.na(idx))) {
      stop("Named `k` values must cover every selected layer id or layer name.")
    }
    out <- k[idx]
    base::names(out) <- layer_ids
    return(out)
  }
  if (base::length(k) != base::length(layer_ids)) {
    stop("Supply one `k` value or one value per selected layer.")
  }
  base::names(k) <- layer_ids
  k
}

.hc_regroup_resolve_modules <- function(modules, cluster_calc) {
  cluster_info <- tryCatch(cluster_calc[["cluster_information"]], error = function(e) NULL)
  if (base::is.null(cluster_info) || !base::is.data.frame(cluster_info) ||
    base::nrow(cluster_info) == 0 ||
    !base::all(c("color", "gene_n") %in% base::colnames(cluster_info))) {
    stop(
      "Current module information is missing. Run `hc_cluster_calculation()` before sample regrouping."
    )
  }
  cluster_info <- base::as.data.frame(cluster_info, stringsAsFactors = FALSE)
  included <- base::rep(TRUE, base::nrow(cluster_info))
  if ("cluster_included" %in% base::colnames(cluster_info)) {
    included <- base::as.character(cluster_info$cluster_included) == "yes"
    included[base::is.na(included)] <- FALSE
  }
  cluster_info <- cluster_info[included, , drop = FALSE]
  available_colors <- base::unique(base::as.character(cluster_info$color))
  available_colors <- available_colors[
    !base::is.na(available_colors) & base::nzchar(available_colors)
  ]
  if (base::length(available_colors) == 0) {
    stop("No included modules are available for sample regrouping.")
  }

  module_prefix <- tryCatch(base::as.character(cluster_calc[["module_prefix"]]), error = function(e) "M")
  if (base::length(module_prefix) != 1 || base::is.na(module_prefix) || !base::nzchar(module_prefix)) {
    module_prefix <- "M"
  }
  label_map <- .hc_normalize_module_label_map_for_split(
    module_label_map = tryCatch(cluster_calc[["module_label_map"]], error = function(e) NULL),
    available_colors = available_colors,
    module_prefix = module_prefix
  )
  module_order <- .hc_split_module_order(
    cluster_calc = cluster_calc,
    available_colors = available_colors
  )

  use_all <- base::is.character(modules) && base::length(modules) == 1 &&
    !base::is.na(modules) && identical(base::tolower(base::trimws(modules)), "all")
  if (use_all) {
    target_colors <- module_order
    target_labels <- base::as.character(label_map[target_colors])
  } else {
    resolved <- .hc_heatmap_view_resolve_modules(
      modules = modules,
      cluster_calc = cluster_calc
    )
    target_colors <- base::as.character(resolved$target_colors)
    target_labels <- base::as.character(resolved$resolved_labels)
  }
  target_labels[base::is.na(target_labels) | !base::nzchar(target_labels)] <-
    target_colors[base::is.na(target_labels) | !base::nzchar(target_labels)]

  gene_map <- .hc_module_gene_map(cluster_info)
  list(
    colors = target_colors,
    labels = target_labels,
    genes = gene_map[target_colors]
  )
}

.hc_regroup_one_layer <- function(hc,
                                  layer_id,
                                  se,
                                  module_info,
                                  k,
                                  reference,
                                  method,
                                  distance,
                                  group_col,
                                  gfc_colors,
                                  gfc_limits,
                                  show_values) {
  assay_names <- SummarizedExperiment::assayNames(se)
  if (base::length(assay_names) == 0) {
    stop("Layer `", layer_id, "` contains no expression assay.")
  }
  assay_name <- if ("counts" %in% assay_names) "counts" else assay_names[[1]]
  expression <- SummarizedExperiment::assay(se, assay_name)
  expression <- base::as.matrix(expression)
  base::storage.mode(expression) <- "double"
  if (base::is.null(base::rownames(expression)) || base::is.null(base::colnames(expression))) {
    stop("Layer `", layer_id, "` expression data must have gene and sample names.")
  }
  if (base::ncol(expression) < 2) {
    stop("Layer `", layer_id, "` needs at least two samples for hierarchical clustering.")
  }
  if (k > base::ncol(expression)) {
    stop("`k` for layer `", layer_id, "` cannot exceed its number of samples.")
  }
  if (base::any(!base::is.finite(expression))) {
    stop("Layer `", layer_id, "` contains non-finite expression values.")
  }

  data_in_log <- .hc_regroup_global_value(hc, "data_in_log", FALSE)
  data_in_log <- base::isTRUE(data_in_log) ||
    base::tolower(base::as.character(data_in_log[[1]])) %in% c("true", "t", "1", "yes")
  if (base::isTRUE(data_in_log)) {
    expression <- 2^expression
  }
  if (base::any(!base::is.finite(expression))) {
    stop("Layer `", layer_id, "` could not be converted from log2 expression safely.")
  }
  if (base::any(expression < 0)) {
    stop(
      "Layer `", layer_id, "` contains negative expression values. ",
      "Set `data_in_log = TRUE` for log2 expression before calculating fold changes."
    )
  }

  reference_values <- .hc_regroup_reference_values(
    hc = hc,
    se = se,
    expression = expression,
    reference = reference,
    layer_id = layer_id
  )
  gene_fc <- .hc_regroup_signed_fold_change(
    expression = expression,
    reference_values = reference_values,
    limits = gfc_limits
  )
  aggregated <- .hc_regroup_aggregate_modules(
    gene_fc = gene_fc,
    module_info = module_info,
    layer_id = layer_id
  )
  sample_gfc <- aggregated$matrix
  if (base::any(!base::is.finite(sample_gfc))) {
    stop("Non-finite module fold changes remained in layer `", layer_id, "`.")
  }

  sample_distance <- .hc_regroup_sample_distance(
    sample_gfc = sample_gfc,
    method = distance,
    layer_id = layer_id
  )
  tree <- stats::hclust(sample_distance, method = method)
  raw_clusters <- stats::cutree(tree, k = k)
  cluster_levels <- base::unique(base::as.character(raw_clusters[tree$order]))
  cluster_map <- stats::setNames(base::seq_along(cluster_levels), cluster_levels)
  cluster_number <- base::as.integer(cluster_map[base::as.character(raw_clusters)])
  base::names(cluster_number) <- base::names(raw_clusters)

  layer_label <- .hc_longitudinal_layer_label(hc = hc, layer_id = layer_id)
  group_prefix <- .hc_longitudinal_safe_suffix(layer_label)
  groups <- base::paste0(group_prefix, "_C", cluster_number)
  clusters <- base::data.frame(
    sample = base::names(cluster_number),
    cluster = cluster_number,
    group = groups,
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  anno <- base::as.data.frame(SummarizedExperiment::colData(se), stringsAsFactors = FALSE)
  anno <- anno[base::colnames(sample_gfc), , drop = FALSE]
  voi <- base::as.character(.hc_regroup_global_value(hc, "voi", NA_character_))
  annotation_col <- base::data.frame(row.names = base::colnames(sample_gfc))
  if (base::length(voi) == 1 && !base::is.na(voi) && base::nzchar(voi) && voi %in% base::colnames(anno)) {
    annotation_col[[voi]] <- base::as.character(anno[[voi]])
  }
  annotation_col[[group_col]] <- groups[base::match(base::rownames(annotation_col), clusters$sample)]

  module_table <- aggregated$module_table
  module_factor <- base::factor(module_table$label, levels = module_table$label)
  annotation_row <- base::data.frame(
    Module = module_factor,
    row.names = base::rownames(sample_gfc),
    check.names = FALSE
  )
  annotation_colors <- list(
    Module = stats::setNames(module_table$color, module_table$label)
  )
  palette <- grDevices::colorRampPalette(gfc_colors)(100)
  breaks <- base::seq(gfc_limits[[1]], gfc_limits[[2]], length.out = base::length(palette) + 1L)
  display_numbers <- if (base::isTRUE(show_values)) base::round(sample_gfc, 2) else FALSE

  heatmap <- pheatmap::pheatmap(
    mat = sample_gfc,
    color = palette,
    breaks = breaks,
    cluster_rows = FALSE,
    cluster_cols = tree,
    cutree_cols = if (k > 1) k else NA_integer_,
    annotation_row = annotation_row,
    annotation_col = annotation_col,
    annotation_colors = annotation_colors,
    show_rownames = TRUE,
    show_colnames = TRUE,
    display_numbers = display_numbers,
    number_format = "%.2f",
    border_color = NA,
    main = base::paste0("Sample module fold changes - ", layer_label),
    silent = TRUE
  )

  list(
    layer_id = layer_id,
    layer_name = layer_label,
    assay = assay_name,
    sample_gfc = sample_gfc,
    module_table = module_table,
    distance = sample_distance,
    tree = tree,
    clusters = clusters,
    plot = heatmap
  )
}

.hc_regroup_global_value <- function(hc, name, default = NULL) {
  global <- hc@config@global
  if (base::nrow(global) == 0 || !(name %in% base::colnames(global))) {
    return(default)
  }
  value <- global[[name]][[1]]
  if (base::is.null(value) || base::length(value) == 0 || base::is.na(value)) {
    return(default)
  }
  value
}

.hc_regroup_reference_values <- function(hc, se, expression, reference, layer_id) {
  if (identical(reference, "mean")) {
    return(base::rowMeans(expression))
  }

  control <- base::as.character(.hc_regroup_global_value(hc, "control", NA_character_))
  voi <- base::as.character(.hc_regroup_global_value(hc, "voi", NA_character_))
  if (base::length(control) != 1 || base::is.na(control) || !base::nzchar(control) ||
    identical(base::tolower(control), "none")) {
    stop("`reference = \"control\"` requires a configured control label.")
  }
  anno <- base::as.data.frame(SummarizedExperiment::colData(se), stringsAsFactors = FALSE)
  if (base::length(voi) != 1 || base::is.na(voi) || !base::nzchar(voi) ||
    !(voi %in% base::colnames(anno))) {
    stop(
      "`reference = \"control\"` requires the configured VOI column in layer `",
      layer_id, "`."
    )
  }
  labels <- base::as.character(anno[[voi]])
  control_samples <- base::grepl(
    base::tolower(control),
    base::tolower(labels),
    fixed = TRUE
  )
  control_samples[base::is.na(control_samples)] <- FALSE
  if (!base::any(control_samples)) {
    stop(
      "No samples matching control label `", control, "` were found in layer `",
      layer_id, "`."
    )
  }
  base::rowMeans(expression[, control_samples, drop = FALSE])
}

.hc_regroup_signed_fold_change <- function(expression, reference_values, limits) {
  ratios <- base::sweep(expression, 1, reference_values, "/")
  out <- ratios
  down <- base::is.finite(ratios) & ratios < 1
  out[down] <- -1 / ratios[down]

  both_zero <- expression == 0 & reference_values[base::row(expression)] == 0
  out[both_zero] <- 1
  out[base::is.infinite(out) & out > 0] <- limits[[2]]
  out[base::is.infinite(out) & out < 0] <- limits[[1]]
  out[out > limits[[2]]] <- limits[[2]]
  out[out < limits[[1]]] <- limits[[1]]
  out
}

.hc_regroup_aggregate_modules <- function(gene_fc, module_info, layer_id) {
  rows <- list()
  table_rows <- list()
  for (i in base::seq_along(module_info$colors)) {
    color <- module_info$colors[[i]]
    label <- module_info$labels[[i]]
    genes <- base::intersect(base::as.character(module_info$genes[[color]]), base::rownames(gene_fc))
    if (base::length(genes) == 0) {
      warning(
        "Skipping module `", label, "` in layer `", layer_id,
        "`: none of its genes are present.",
        call. = FALSE
      )
      next
    }
    rows[[base::length(rows) + 1L]] <- base::colMeans(gene_fc[genes, , drop = FALSE])
    table_rows[[base::length(table_rows) + 1L]] <- base::data.frame(
      label = label,
      color = color,
      n_genes = base::length(genes),
      stringsAsFactors = FALSE
    )
  }
  if (base::length(rows) == 0) {
    stop("No module genes overlap the expression matrix in layer `", layer_id, "`.")
  }
  mat <- base::do.call(base::rbind, rows)
  module_table <- base::do.call(base::rbind, table_rows)
  base::rownames(mat) <- base::make.unique(base::as.character(module_table$label))
  base::colnames(mat) <- base::colnames(gene_fc)
  base::rownames(module_table) <- base::rownames(mat)
  list(matrix = mat, module_table = module_table)
}

.hc_regroup_sample_distance <- function(sample_gfc, method, layer_id) {
  if (!identical(method, "correlation")) {
    return(stats::dist(base::t(sample_gfc), method = method))
  }
  correlations <- stats::cor(sample_gfc, use = "pairwise.complete.obs", method = "pearson")
  if (base::any(!base::is.finite(correlations))) {
    stop(
      "Correlation distance is undefined for at least one sample in layer `",
      layer_id, "`; use `distance = \"euclidean\"` or select more variable modules."
    )
  }
  stats::as.dist(1 - correlations)
}
