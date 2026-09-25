#' Plot a fixed-scale view of selected modules from the module heatmap
#'
#' Creates a module heatmap view for selected modules in the exact order
#' supplied by `modules`. The module labels and module-color annotations are
#' preserved, and the GFC expression color scale is inherited from the full
#' heatmap cache instead of being recomputed from the subset.
#'
#' @param hc A `HCoCenaExperiment`.
#' @param modules Character/numeric vector of modules to show. Accepts current
#'   module labels such as `"M2"` or `"M2.1"`, module colors, or numeric module
#'   indices in the current heatmap order.
#' @param file_name Optional export file name for the view. Use `FALSE` to skip
#'   file export.
#' @param col_order Optional condition order. If `NULL`, the stored full
#'   heatmap column order is reused when available.
#' @param cluster_columns Logical. Whether to cluster columns in the view.
#' @param cluster_rows Logical. Whether to cluster rows in the view. Defaults to
#'   `FALSE` so the `modules` vector controls the displayed order.
#' @param gfc_colors Optional GFC palette override. If `NULL`, the stored full
#'   heatmap palette is reused.
#' @param gfc_scale_limits Optional GFC scale override. If `NULL`, the stored
#'   full heatmap scale is reused; if unavailable, a full-cache/global fallback
#'   is used.
#' @param preserve_full_heatmap Logical. If `TRUE`, restore the pre-existing
#'   full heatmap cache in `hc` after writing the view.
#' @param store_view Logical. If `TRUE`, store lightweight metadata for the
#'   view under `hc@integration@cluster$heatmap_views`.
#' @param view_name Optional name used when `store_view = TRUE`.
#' @param ... Additional plotting arguments forwarded to
#'   [hc_plot_cluster_heatmap()].
#'
#' @examples
#' hc <- hc_example_data("clustered")
#' hc <- hc_plot_cluster_heatmap_view(
#'   hc,
#'   modules = c("M2", "M3"),
#'   file_name = FALSE
#' )
#'
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_plot_cluster_heatmap_view <- function(hc,
                                         modules,
                                         file_name = "Heatmap_modules_view.pdf",
                                         col_order = NULL,
                                         cluster_columns = FALSE,
                                         cluster_rows = FALSE,
                                         gfc_colors = NULL,
                                         gfc_scale_limits = NULL,
                                         preserve_full_heatmap = TRUE,
                                         store_view = TRUE,
                                         view_name = NULL,
                                         ...) {
  if (!inherits(hc, "HCoCenaExperiment")) {
    stop("`hc` must be a `HCoCenaExperiment`.")
  }
  if (missing(modules) || base::is.null(modules) || base::length(modules) == 0) {
    stop("`modules` must contain at least one module label, color, or index.")
  }
  if (!base::is.logical(cluster_columns) || base::length(cluster_columns) != 1 ||
    base::is.na(cluster_columns)) {
    stop("`cluster_columns` must be TRUE or FALSE.")
  }
  if (!base::is.logical(cluster_rows) || base::length(cluster_rows) != 1 ||
    base::is.na(cluster_rows)) {
    stop("`cluster_rows` must be TRUE or FALSE.")
  }
  if (!base::is.logical(preserve_full_heatmap) ||
    base::length(preserve_full_heatmap) != 1 ||
    base::is.na(preserve_full_heatmap)) {
    stop("`preserve_full_heatmap` must be TRUE or FALSE.")
  }
  if (!base::is.logical(store_view) ||
    base::length(store_view) != 1 ||
    base::is.na(store_view)) {
    stop("`store_view` must be TRUE or FALSE.")
  }

  dot_args <- list(...)
  if ("row_order" %in% base::names(dot_args)) {
    stop("Do not pass `row_order`; use `modules` to control the view order.")
  }
  if ("write_module_tables" %in% base::names(dot_args)) {
    stop("`write_module_tables` is managed by `hc_plot_cluster_heatmap_view()`.")
  }

  cluster_calc <- .hc_heatmap_view_cluster_calc(hc)
  resolved_modules <- .hc_heatmap_view_resolve_modules(
    modules = modules,
    cluster_calc = cluster_calc
  )
  if (base::is.null(col_order)) {
    heatmap_info <- tryCatch(.hc_heatmap_cache_info(cluster_calc), error = function(e) NULL)
    if (!base::is.null(heatmap_info) &&
      !base::is.null(heatmap_info$col_order) &&
      base::length(heatmap_info$col_order) > 0) {
      col_order <- base::as.character(heatmap_info$col_order)
    }
  }

  resolved_gfc_colors <- .hc_heatmap_view_resolve_gfc_colors(
    cluster_calc = cluster_calc,
    gfc_colors = gfc_colors
  )
  resolved_gfc_scale_limits <- .hc_heatmap_view_resolve_gfc_scale_limits(
    hc = hc,
    cluster_calc = cluster_calc,
    gfc_scale_limits = gfc_scale_limits
  )

  plot_args <- list(
    hc = hc,
    file_name = file_name,
    row_order = resolved_modules$target_colors,
    col_order = col_order,
    cluster_rows = cluster_rows,
    cluster_columns = cluster_columns,
    gfc_colors = resolved_gfc_colors,
    gfc_scale_limits = resolved_gfc_scale_limits,
    write_module_tables = FALSE
  )
  if (!("module_label_numbering" %in% base::names(dot_args))) {
    plot_args[["module_label_numbering"]] <- "preserve_existing"
  }
  plot_args <- base::c(plot_args, dot_args)

  plotted_hc <- base::do.call(hc_plot_cluster_heatmap, plot_args)
  view_state <- .hc_heatmap_view_capture_state(
    hc = plotted_hc,
    requested_modules = modules,
    resolved_modules = resolved_modules,
    file_name = file_name,
    gfc_colors = resolved_gfc_colors,
    gfc_scale_limits = resolved_gfc_scale_limits
  )

  out <- if (base::isTRUE(preserve_full_heatmap)) {
    .hc_heatmap_view_restore_full_state(
      plotted_hc = plotted_hc,
      original_hc = hc
    )
  } else {
    plotted_hc
  }

  if (base::isTRUE(store_view)) {
    out <- .hc_heatmap_view_store_state(
      hc = out,
      view_name = view_name,
      view_state = view_state
    )
  }

  methods::validObject(out)
  out
}

.hc_heatmap_view_cluster_calc <- function(hc) {
  cluster_calc <- tryCatch(hc@integration@cluster, error = function(e) NULL)
  if (base::is.null(cluster_calc)) {
    return(list())
  }
  base::as.list(cluster_calc)
}

.hc_heatmap_view_available_modules <- function(cluster_calc) {
  cluster_info <- tryCatch(cluster_calc[["cluster_information"]], error = function(e) NULL)
  if (base::is.data.frame(cluster_info) &&
    "color" %in% base::colnames(cluster_info) &&
    base::nrow(cluster_info) > 0) {
    included <- base::rep(TRUE, base::nrow(cluster_info))
    if ("cluster_included" %in% base::colnames(cluster_info)) {
      included <- base::as.character(cluster_info$cluster_included) == "yes"
      included[base::is.na(included)] <- FALSE
    }
    colors <- base::as.character(cluster_info$color[included])
    colors <- colors[!base::is.na(colors) & base::nzchar(colors)]
    if (base::length(colors) > 0) {
      return(base::unique(colors))
    }
  }

  heatmap_info <- tryCatch(.hc_heatmap_cache_info(cluster_calc), error = function(e) NULL)
  if (!base::is.null(heatmap_info) &&
    !base::is.null(heatmap_info$matrix) &&
    !base::is.null(base::rownames(heatmap_info$matrix))) {
    row_ids <- base::as.character(base::rownames(heatmap_info$matrix))
    row_ids <- row_ids[!base::is.na(row_ids) & base::nzchar(row_ids)]
    if (base::length(row_ids) > 0) {
      return(base::unique(row_ids))
    }
  }

  module_label_map <- tryCatch(cluster_calc[["module_label_map"]], error = function(e) NULL)
  module_names <- base::names(module_label_map)
  module_names <- base::as.character(module_names)
  module_names <- module_names[!base::is.na(module_names) & base::nzchar(module_names)]
  base::unique(module_names)
}

.hc_heatmap_view_resolve_modules <- function(modules, cluster_calc) {
  available_colors <- .hc_heatmap_view_available_modules(cluster_calc)
  if (base::length(available_colors) == 0) {
    stop(
      "No module identifiers are available. Run `hc_plot_cluster_heatmap()` ",
      "or check that `hc@integration@cluster` contains module information."
    )
  }

  module_prefix <- tryCatch(cluster_calc[["module_prefix"]], error = function(e) NULL)
  module_prefix <- base::as.character(module_prefix)
  if (base::length(module_prefix) != 1 ||
    base::is.na(module_prefix) ||
    !base::nzchar(module_prefix)) {
    module_prefix <- "M"
  }
  module_label_map <- .hc_normalize_module_label_map_for_split(
    module_label_map = tryCatch(cluster_calc[["module_label_map"]], error = function(e) NULL),
    available_colors = available_colors,
    module_prefix = module_prefix
  )
  module_order <- .hc_split_module_order(
    cluster_calc = cluster_calc,
    available_colors = available_colors
  )
  resolved <- .hc_resolve_modules_for_split(
    modules = modules,
    available_colors = available_colors,
    module_label_map = module_label_map,
    module_order = module_order
  )

  unresolved_tbl <- resolved$resolution_table[
    base::as.character(resolved$resolution_table$status) != "ok", ,
    drop = FALSE
  ]
  if (base::nrow(unresolved_tbl) > 0) {
    unresolved_inputs <- base::unique(base::as.character(unresolved_tbl$input))
    available_labels <- base::unique(base::as.character(module_label_map[module_order]))
    available_labels <- available_labels[!base::is.na(available_labels) & base::nzchar(available_labels)]
    if (base::length(available_labels) > 0) {
      available_preview <- base::paste(utils::head(available_labels, 12L), collapse = ", ")
      if (base::length(available_labels) > 12L) {
        available_preview <- base::paste0(available_preview, ", ...")
      }
      stop(
        "Could not resolve the following module identifier(s): ",
        base::paste(unresolved_inputs, collapse = ", "),
        "\nUse exact current module labels, module colors, or numeric indices.",
        "\nAvailable module labels: ", available_preview
      )
    }
    stop(
      "Could not resolve the following module identifier(s): ",
      base::paste(unresolved_inputs, collapse = ", "),
      "\nUse exact current module labels, module colors, or numeric indices."
    )
  }
  if (base::length(resolved$target_colors) == 0) {
    stop(
      "Could not match any requested module in `modules`.\n",
      "Use module labels, module colors, or module indices."
    )
  }

  resolved
}

.hc_heatmap_view_normalize_gfc_scale_limits <- function(x) {
  if (base::is.null(x)) {
    return(NULL)
  }
  x <- .hc_as_numeric_safely(x)
  if (base::length(x) == 1) {
    if (!base::is.finite(x) || x <= 0) {
      stop("`gfc_scale_limits` as single value must be finite and > 0.")
    }
    return(c(-base::abs(x), base::abs(x)))
  }
  if (base::length(x) != 2 || base::any(!base::is.finite(x))) {
    stop("`gfc_scale_limits` must be NULL, one positive number, or a numeric vector of length 2.")
  }
  x <- base::sort(x)
  if (x[[1]] == x[[2]]) {
    stop("`gfc_scale_limits` must have different min/max values.")
  }
  x
}

.hc_heatmap_view_resolve_gfc_scale_limits <- function(hc,
                                                      cluster_calc,
                                                      gfc_scale_limits = NULL) {
  resolved <- .hc_heatmap_view_normalize_gfc_scale_limits(gfc_scale_limits)
  if (!base::is.null(resolved)) {
    return(resolved)
  }

  stored <- tryCatch(
    .hc_heatmap_view_normalize_gfc_scale_limits(cluster_calc[["gfc_scale_limits"]]),
    error = function(e) NULL
  )
  if (!base::is.null(stored)) {
    return(stored)
  }

  heatmap_info <- tryCatch(.hc_heatmap_cache_info(cluster_calc), error = function(e) NULL)
  mat <- if (!base::is.null(heatmap_info)) heatmap_info$matrix else NULL
  if (!base::is.null(mat)) {
    mat_num <- base::suppressWarnings(base::as.numeric(mat))
    mat_num <- mat_num[base::is.finite(mat_num)]
    if (base::length(mat_num) > 0) {
      lim <- base::max(base::abs(mat_num), na.rm = TRUE)
      if (base::is.finite(lim) && lim > 0) {
        return(c(-base::abs(lim), base::abs(lim)))
      }
    }
  }

  bridge <- tryCatch(.hc_as_bridge_object_for_cluster_plot(hc), error = function(e) NULL)
  fallback_lim <- tryCatch(
    .hc_first_numeric_value(bridge[["global_settings"]][["range_GFC"]]),
    error = function(e) NA_real_
  )
  if (base::length(fallback_lim) != 1 ||
    !base::is.finite(fallback_lim) ||
    fallback_lim <= 0) {
    fallback_lim <- 2
  }
  c(-base::abs(fallback_lim), base::abs(fallback_lim))
}

.hc_heatmap_view_resolve_gfc_colors <- function(cluster_calc, gfc_colors = NULL) {
  colors <- gfc_colors
  if (base::is.null(colors)) {
    colors <- tryCatch(base::as.character(cluster_calc[["gfc_colors"]]), error = function(e) NULL)
  }
  if (base::is.null(colors) ||
    base::length(colors) < 2 ||
    !base::is.character(colors) ||
    base::any(base::is.na(colors)) ||
    base::any(colors == "")) {
    colors <- .hc_default_gfc_colors()
  }
  base::as.character(colors)
}

.hc_heatmap_view_capture_state <- function(hc,
                                           requested_modules,
                                           resolved_modules,
                                           file_name,
                                           gfc_colors,
                                           gfc_scale_limits) {
  cluster_calc <- .hc_heatmap_view_cluster_calc(hc)
  keep_names <- c(
    "heatmap_output_files",
    "heatmap_matrix",
    "heatmap_row_order",
    "heatmap_column_order",
    "heatmap_column_labels_display",
    "heatmap_cluster",
    "heatmap_cluster_raw",
    "module_label_map",
    "module_label_mode",
    "module_label_numbering",
    "gfc_colors",
    "gfc_scale_limits",
    "overall_plot_scale"
  )
  view_state <- cluster_calc[base::intersect(keep_names, base::names(cluster_calc))]
  view_state[["requested_modules"]] <- requested_modules
  view_state[["resolved_modules"]] <- resolved_modules$target_colors
  view_state[["resolved_labels"]] <- resolved_modules$resolved_labels
  view_state[["resolution_table"]] <- resolved_modules$resolution_table
  view_state[["file_name"]] <- file_name
  view_state[["gfc_colors"]] <- gfc_colors
  view_state[["gfc_scale_limits"]] <- gfc_scale_limits
  view_state[["created_at"]] <- base::Sys.time()
  view_state
}

.hc_heatmap_view_restore_full_state <- function(plotted_hc, original_hc) {
  plotted_hc@integration@cluster <- original_hc@integration@cluster
  plotted_hc@satellite <- original_hc@satellite
  plotted_hc
}

.hc_heatmap_view_store_state <- function(hc, view_name = NULL, view_state) {
  if (base::is.null(view_name)) {
    view_name <- base::paste0("view_", base::format(base::Sys.time(), "%Y%m%d_%H%M%S"))
  }
  if (!base::is.character(view_name) ||
    base::length(view_name) != 1 ||
    base::is.na(view_name) ||
    !base::nzchar(base::trimws(view_name))) {
    stop("`view_name` must be NULL or a non-empty character scalar.")
  }
  view_name <- base::trimws(view_name)

  cluster_calc <- .hc_heatmap_view_cluster_calc(hc)
  views <- tryCatch(cluster_calc[["heatmap_views"]], error = function(e) NULL)
  if (base::is.null(views)) {
    views <- list()
  } else {
    views <- base::as.list(views)
  }
  views[[view_name]] <- view_state
  cluster_calc[["heatmap_views"]] <- views

  hc@integration@cluster <- base::do.call(S4Vectors::SimpleList, cluster_calc)
  methods::validObject(hc)
  hc
}
