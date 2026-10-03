#' Upstream regulator and pathway inference
#'
#' Links modules to upstream regulators: TF regulons (CollecTRI, or DoRothEA)
#' and pathway-responsive genes (PROGENy, optional custom GMT files).
#' Significance is an over-representation test of the regulator's targets
#' among the module genes against all network genes, corrected over all
#' module x regulator tests of a resource. The signed score per condition
#' describes how the module's targets move (sign of regulation taken into
#' account). Results are exported as Excel tables and plotted as heatmaps.
#'
#' @param resources Character vector of upstream resources to use.
#'   Allowed values are `"TF"` and `"Pathway"`. Default is both.
#' @param top Integer. Number of top significant regulators/pathways to keep
#'   per module and resource in the selected summary. Default is 5.
#' @param clusters Either `"all"` (default) or a character vector of module
#'   colors to process.
#' @param padj Multiple-testing correction method passed to
#'   [stats::p.adjust()]. Default is `"BH"`.
#' @param qval Adjusted p-value threshold for significance. Default is 0.05.
#' @param tf_resource TF regulon resource: `"auto"` (default) uses CollecTRI
#'   when it can be loaded (needs `OmnipathR` and internet) and DoRothEA
#'   otherwise; `"collectri"` or `"dorothea"` force one of them.
#' @param tf_confidence Character vector of DoRothEA confidence levels to keep
#'   (DoRothEA only). Default is `c("A", "B", "C")`.
#' @param progeny_top Number of most responsive genes per PROGENy pathway.
#'   Default is 100.
#' @param minsize Minimum number of targets (within the network genes) a
#'   regulator needs to be tested, and minimum module size. Default is 5.
#' @param min_overlap Minimum number of a regulator's targets in a module for
#'   the link to count as significant. Default is 3.
#' @param collapse_redundant Logical. If `TRUE` (default), regulators whose
#'   module targets largely repeat a better-ranked regulator of the same module
#'   are flagged in `redundant_with` and left out of the selected (plotted)
#'   regulators.
#' @param redundancy_jaccard Jaccard index of the module targets from which a
#'   regulator counts as redundant. Default is 0.5.
#' @param method Test used to link modules and regulators. Only `"ora"`
#'   (hypergeometric over-representation) is supported.
#' @param activity_input Character scalar selecting the values used for the
#'   signed score per condition:
#'   `"gfc"` (default) uses `integrated_output$GFC_all_layers`,
#'   `"fc"` uses user-defined pairwise log2 fold-changes from
#'   `fc_comparisons`.
#' @param fc_comparisons Character vector of pairwise comparisons used only when
#'   `activity_input = "fc"`. Each entry must be formatted as
#'   `"groupA_vs_groupB"` (numerator vs denominator), e.g.
#'   `c("IFNg_seq_vs_baseline_seq", "IL4_seq_vs_baseline_seq")`.
#' @param custom_pathway_gmt Optional custom pathway GMT file(s) added to the
#'   pathway inference resource. Accepts a character vector or named list of
#'   file paths. Paths can be absolute/relative or file names inside
#'   `dir_reference_files`.
#' @param heatmap_side Position of the hCoCena heatmap in the combined output.
#'   Choose one of `"left"` (default) or `"right"`.
#' @param cluster_columns Logical. If `FALSE` (default), reuse the
#'   column order from the main hCoCena heatmap when available. If `TRUE`,
#'   cluster the columns for this upstream plot instead.
#' @param heatmap_cluster_columns Legacy alias for `cluster_columns`.
#' @param col_order Optional character vector overriding the hCoCena
#'   heatmap column order for this upstream plot only. If `NULL` (default),
#'   the column order from the main module heatmap is reused when available.
#' @param heatmap_col_order Legacy alias for `col_order`.
#' @param gfc_scale_limits Optional numeric vector controlling the module-heatmap
#'   color scale limits used in upstream combined heatmaps (left/right hCoCena
#'   panel). Provide one positive number (`x` -> `c(-x, x)`) or two numbers
#'   (`c(min, max)`). If NULL, uses stored limits from the main heatmap when
#'   available, otherwise falls back to `c(-range_GFC, range_GFC)`.
#' @param plot Logical; if `TRUE` (default), draw plot outputs in the active
#'   graphics device.
#' @param save_pdf Logical; if `TRUE` (default), export plots to
#'   `Upstream_Inference.pdf`.
#' @param pdf_width Optional numeric width (inches) for `Upstream_Inference.pdf`.
#'   If NULL (default), width is auto-estimated from content.
#' @param pdf_height Optional numeric height (inches) for `Upstream_Inference.pdf`.
#'   If NULL (default), height is auto-estimated from content.
#' @param pdf_pointsize Numeric base pointsize used for the upstream PDF device.
#'   Default is 11.
#' @param plot_per_comparison Logical; if `TRUE`, additionally create one
#'   combined upstream heatmap page per activity column (GFC condition or FC
#'   comparison), each with the matching one-column module heatmap.
#' @param consistent_terms Logical; controls term comparability when
#'   `plot_per_comparison = TRUE`.
#'   If `TRUE`, per-comparison pages still show only values from the currently
#'   shown condition, but use a global (all-condition) term axis for
#'   comparability; `*` marks significance for the currently shown condition.
#'   If `FALSE`, each page uses only local selected activities from
#'   the shown condition and no significance marker is drawn.
#' @param overall_plot_scale Numeric scaling factor for plot typography and
#'   marker sizes. Default is 1.
#'
#' @return A named list with selected/significant summaries, per-resource
#'   summaries and plot objects.
#' @noRd
.hc_upstream_inference_driver <- function(resources = c("TF", "Pathway"),
                               top = 5,
                               clusters = c("all"),
                               padj = "BH",
                               qval = 0.05,
                               tf_resource = "auto",
                               tf_confidence = c("A", "B", "C"),
                               progeny_top = 100,
                               minsize = 5,
                               min_overlap = 3,
                               collapse_redundant = TRUE,
                               redundancy_jaccard = 0.5,
                               method = "ora",
                               activity_input = "gfc",
                               fc_comparisons = NULL,
                               custom_pathway_gmt = NULL,
                               heatmap_side = "left",
                               cluster_columns = FALSE,
                               heatmap_cluster_columns = NULL,
                               col_order = NULL,
                               heatmap_col_order = NULL,
                               gfc_scale_limits = NULL,
                               plot = TRUE,
                               save_pdf = TRUE,
                               pdf_width = NULL,
                               pdf_height = NULL,
                               pdf_pointsize = 11,
                               plot_per_comparison = TRUE,
                               consistent_terms = TRUE,
                               overall_plot_scale = 1) {

  if (!requireNamespace("decoupleR", quietly = TRUE)) {
    stop(
      "`.hc_upstream_inference_driver()` requires package `decoupleR`. ",
      "Please install it first (e.g. `BiocManager::install('decoupleR')`)."
    )
  }

  resources <- base::toupper(base::as.character(resources))
  resources <- base::unique(resources)
  allowed_resources <- c("TF", "PATHWAY")
  invalid_resources <- base::setdiff(resources, allowed_resources)
  if (base::length(invalid_resources) > 0) {
    stop(
      "Unknown value(s) in `resources`: ",
      base::paste(invalid_resources, collapse = ", "),
      ". Allowed values are `TF` and `Pathway`."
    )
  }
  if (base::length(resources) == 0) {
    stop("`resources` must contain at least one entry.")
  }
  if (!base::is.numeric(top) || base::length(top) != 1 || base::is.na(top) || top < 1) {
    stop("`top` must be a positive integer.")
  }
  top <- base::as.integer(top)
  if (!(padj %in% stats::p.adjust.methods)) {
    stop(
      "`padj` must be one of: ",
      base::paste(stats::p.adjust.methods, collapse = ", ")
    )
  }
  if (!base::is.numeric(qval) || base::length(qval) != 1 || base::is.na(qval) || qval <= 0 || qval > 1) {
    stop("`qval` must be a numeric value in (0, 1].")
  }
  if (!base::is.numeric(minsize) || base::length(minsize) != 1 || base::is.na(minsize) || minsize < 1) {
    stop("`minsize` must be a positive integer.")
  }
  minsize <- base::as.integer(minsize)
  if (!base::is.numeric(min_overlap) || base::length(min_overlap) != 1 || base::is.na(min_overlap) || min_overlap < 1) {
    stop("`min_overlap` must be a positive integer.")
  }
  min_overlap <- base::as.integer(min_overlap)
  if (!base::is.logical(collapse_redundant) || base::length(collapse_redundant) != 1 || base::is.na(collapse_redundant)) {
    stop("`collapse_redundant` must be TRUE or FALSE.")
  }
  if (!base::is.numeric(redundancy_jaccard) || base::length(redundancy_jaccard) != 1 ||
    base::is.na(redundancy_jaccard) || redundancy_jaccard <= 0 || redundancy_jaccard > 1) {
    stop("`redundancy_jaccard` must be a number in (0, 1].")
  }
  if (!base::is.logical(plot) || base::length(plot) != 1) {
    stop("`plot` must be TRUE or FALSE.")
  }
  if (!base::is.logical(save_pdf) || base::length(save_pdf) != 1) {
    stop("`save_pdf` must be TRUE or FALSE.")
  }
  if (!base::is.null(pdf_width) &&
    (!base::is.numeric(pdf_width) || base::length(pdf_width) != 1 || base::is.na(pdf_width) || pdf_width <= 0)) {
    stop("`pdf_width` must be NULL or a single positive number.")
  }
  if (!base::is.null(pdf_height) &&
    (!base::is.numeric(pdf_height) || base::length(pdf_height) != 1 || base::is.na(pdf_height) || pdf_height <= 0)) {
    stop("`pdf_height` must be NULL or a single positive number.")
  }
  if (!base::is.numeric(pdf_pointsize) ||
    base::length(pdf_pointsize) != 1 ||
    base::is.na(pdf_pointsize) ||
    pdf_pointsize <= 0) {
    stop("`pdf_pointsize` must be a single positive number.")
  }
  if (!base::is.logical(plot_per_comparison) || base::length(plot_per_comparison) != 1 || base::is.na(plot_per_comparison)) {
    stop("`plot_per_comparison` must be TRUE or FALSE.")
  }
  if (!base::is.logical(consistent_terms) || base::length(consistent_terms) != 1 || base::is.na(consistent_terms)) {
    stop("`consistent_terms` must be TRUE or FALSE.")
  }
  if (isTRUE(plot_per_comparison)) {
    if (isTRUE(consistent_terms)) {
      message(
        ".hc_upstream_inference_driver(): per-condition comparison mode enabled ",
        "(`consistent_terms = TRUE`). ",
        "Term axis is fixed across conditions; `*` marks regulators linked to the module and active in the shown condition."
      )
    } else {
      message(
        ".hc_upstream_inference_driver(): per-condition discovery mode enabled ",
        "(`consistent_terms = FALSE`). ",
        "Each condition uses its own term set (not intended for strict cross-condition comparison)."
      )
    }
  }
  if (!base::is.numeric(overall_plot_scale) ||
    base::length(overall_plot_scale) != 1 ||
    base::is.na(overall_plot_scale) ||
    overall_plot_scale <= 0) {
    stop("`overall_plot_scale` must be a positive numeric scalar.")
  }
  overall_plot_scale <- base::max(0.5, base::min(3, overall_plot_scale))
  col_order <- .hc_resolve_col_order_alias(
    col_order = col_order,
    heatmap_col_order = heatmap_col_order,
    col_order_missing = missing(col_order),
    heatmap_col_order_missing = missing(heatmap_col_order),
    context = ".hc_upstream_inference_driver()"
  )
  cluster_columns <- .hc_resolve_cluster_columns_alias(
    cluster_columns = cluster_columns,
    heatmap_cluster_columns = heatmap_cluster_columns,
    cluster_columns_missing = missing(cluster_columns),
    heatmap_cluster_columns_missing = missing(heatmap_cluster_columns),
    context = ".hc_upstream_inference_driver()"
  )
  if (!base::is.logical(cluster_columns) ||
    base::length(cluster_columns) != 1 ||
    base::is.na(cluster_columns)) {
    stop("`cluster_columns` must be TRUE or FALSE.")
  }
  method <- base::tolower(base::as.character(method[[1]]))
  if (identical(method, "ulm")) {
    stop(
      "`method = 'ulm'` was removed: ULM restricted to the genes of one ",
      "co-expressed module compares targets with genes that move the same way. ",
      "Modules are now linked to regulators by over-representation ",
      "(`method = 'ora'`)."
    )
  }
  if (!method %in% "ora") {
    stop("Only `method = 'ora'` is supported.")
  }
  tf_resource <- base::match.arg(
    base::tolower(base::as.character(tf_resource[[1]])),
    choices = c("auto", "collectri", "dorothea")
  )
  if (!base::is.numeric(progeny_top) || base::length(progeny_top) != 1 ||
    base::is.na(progeny_top) || progeny_top < 1) {
    stop("`progeny_top` must be a positive integer.")
  }
  progeny_top <- base::as.integer(progeny_top)
  activity_input <- base::tolower(base::as.character(activity_input[[1]]))
  if (identical(activity_input, "expression")) {
    stop(
      "`activity_input = 'expression'` was removed: absolute expression ",
      "levels describe how abundant the targets are, not how they change. ",
      "Use `'gfc'` or `'fc'`."
    )
  }
  activity_input <- base::match.arg(activity_input, choices = c("gfc", "fc"))
  heatmap_side <- base::match.arg(base::tolower(base::as.character(heatmap_side)), choices = c("left", "right"))

  normalize_scale_limits <- function(x) {
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
    if (base::length(x) != 2 || any(!base::is.finite(x))) {
      stop("`gfc_scale_limits` must be NULL, one positive number, or a numeric vector of length 2.")
    }
    x <- base::sort(x)
    if (x[1] == x[2]) {
      stop("`gfc_scale_limits` must have different min/max values.")
    }
    x
  }

  resolve_gfc_scale_limits <- function(input_limits) {
    lims <- normalize_scale_limits(input_limits)
    if (!base::is.null(lims)) {
      return(lims)
    }
    stored <- tryCatch(
      normalize_scale_limits(hcobject[["integrated_output"]][["cluster_calc"]][["gfc_scale_limits"]]),
      error = function(e) NULL
    )
    if (!base::is.null(stored)) {
      return(stored)
    }
    fallback_lim <- .hc_first_numeric_value(hcobject[["global_settings"]][["range_GFC"]])
    if (!base::is.finite(fallback_lim) || fallback_lim <= 0) {
      fallback_lim <- 2
    }
    c(-base::abs(fallback_lim), base::abs(fallback_lim))
  }

  compute_scale_ticks <- function(lims) {
    br <- pretty(lims, n = 5)
    br <- br[
      br >= lims[1] - .Machine$double.eps^0.5 &
        br <= lims[2] + .Machine$double.eps^0.5
    ]
    if (!any(base::abs(br) < .Machine$double.eps^0.5)) {
      br <- base::sort(base::unique(base::c(br, 0)))
    }
    if (base::length(br) < 3) {
      br <- base::seq(lims[1], lims[2], length.out = 5)
    }
    tick_labels <- base::formatC(br, format = "fg", digits = 3)
    tick_labels <- base::trimws(tick_labels)
    tick_label_width <- base::max(base::nchar(tick_labels), na.rm = TRUE)
    tick_labels <- base::format(tick_labels, width = tick_label_width, justify = "right")
    list(
      breaks = br,
      labels = tick_labels
    )
  }

  gfc_scale_limits <- resolve_gfc_scale_limits(gfc_scale_limits)
  gfc_scale_ticks <- compute_scale_ticks(gfc_scale_limits)

  if (base::is.null(hcobject[["integrated_output"]][["cluster_calc"]][["cluster_information"]])) {
    stop("No cluster information found. Run `hc_cluster_calculation()` first.")
  }
  cluster_info <- hcobject[["integrated_output"]][["cluster_calc"]][["cluster_information"]]
  all_clusters <- base::unique(base::as.character(cluster_info$color))
  all_clusters <- all_clusters[all_clusters != "white" & !base::is.na(all_clusters)]
  if (base::length(all_clusters) == 0) {
    stop("No non-white modules found for upstream inference.")
  }
  if (clusters[1] == "all") {
    clusters <- all_clusters
  }
  clusters <- base::as.character(clusters)
  missing_clusters <- base::setdiff(clusters, all_clusters)
  if (base::length(missing_clusters) > 0) {
    warning(
      "Ignoring unknown cluster(s): ",
      base::paste(missing_clusters, collapse = ", ")
    )
    clusters <- clusters[clusters %in% all_clusters]
  }
  if (base::length(clusters) == 0) {
    stop("No valid clusters selected for upstream inference.")
  }

  cluster_calc <- hcobject[["integrated_output"]][["cluster_calc"]]
  cluster_order <- all_clusters
  heatmap_info <- .hc_heatmap_cache_info(cluster_calc)
  stored_hm <- heatmap_info$heatmap_obj
  main_heatmap_col_order <- heatmap_info$col_order
  if (!base::is.null(heatmap_info$row_order) && base::length(heatmap_info$row_order) > 0) {
    cluster_order <- heatmap_info$row_order
    cluster_order <- cluster_order[cluster_order %in% all_clusters]
    cluster_order <- base::c(cluster_order, base::setdiff(all_clusters, cluster_order))
  }
  cluster_order <- cluster_order[cluster_order %in% clusters]
  if (base::length(cluster_order) == 0) {
    stop("No modules remain after applying ordering and `clusters` filter.")
  }

  module_prefix <- cluster_calc[["module_prefix"]]
  if (base::is.null(module_prefix) || !base::is.character(module_prefix) || base::length(module_prefix) != 1) {
    module_prefix <- "M"
  }
  module_label_map <- cluster_calc[["module_label_map"]]
  if (!base::is.null(module_label_map) && base::length(module_label_map) > 0) {
    map_names <- base::names(cluster_calc[["module_label_map"]])
    module_label_map <- base::as.character(module_label_map)
    if (!base::is.null(map_names) && base::length(map_names) == base::length(module_label_map)) {
      base::names(module_label_map) <- base::as.character(map_names)
    }
    missing_before <- base::setdiff(cluster_order, base::names(module_label_map))
    if (base::length(missing_before) > 0) {
      inverse_map <- stats::setNames(base::names(module_label_map), base::as.character(module_label_map))
      if (base::all(cluster_order %in% base::names(inverse_map))) {
        module_label_map <- inverse_map
      }
    }
  }
  if (base::is.null(module_label_map) || base::length(module_label_map) == 0) {
    module_label_map <- stats::setNames(
      base::paste0(module_prefix, base::seq_along(cluster_order)),
      cluster_order
    )
  }
  missing_map <- base::setdiff(cluster_order, base::names(module_label_map))
  if (base::length(missing_map) > 0) {
    module_label_map <- base::c(
      module_label_map,
      stats::setNames(
        base::paste0(module_prefix, base::seq.int(base::length(module_label_map) + 1, length.out = base::length(missing_map))),
        missing_map
      )
    )
  }
  module_label_map_current <- module_label_map[cluster_order]

  gfc_all <- hcobject[["integrated_output"]][["GFC_all_layers"]]
  if (base::is.null(gfc_all) || base::nrow(gfc_all) == 0 || base::ncol(gfc_all) < 2) {
    stop("Missing or invalid `integrated_output$GFC_all_layers`.")
  }
  .hc_ui_collapse_duplicate_rows <- function(mat) {
    if (!base::anyDuplicated(base::rownames(mat))) {
      return(mat)
    }
    split_idx <- base::split(base::seq_len(base::nrow(mat)), base::rownames(mat))
    out <- base::t(base::vapply(
      split_idx,
      function(idx) {
        vals <- base::colMeans(mat[idx, , drop = FALSE], na.rm = TRUE)
        vals[base::is.nan(vals)] <- NA_real_
        vals
      },
      FUN.VALUE = base::numeric(base::ncol(mat))
    ))
    out <- out %>% base::as.matrix()
    out
  }
  .hc_ui_collapse_duplicate_columns <- function(mat) {
    if (!base::anyDuplicated(base::colnames(mat))) {
      return(mat)
    }
    split_idx <- base::split(base::seq_len(base::ncol(mat)), base::colnames(mat))
    out <- base::vapply(
      split_idx,
      function(idx) {
        if (base::length(idx) == 1) {
          return(mat[, idx])
        }
        vals <- base::rowMeans(mat[, idx, drop = FALSE], na.rm = TRUE)
        vals[base::is.nan(vals)] <- NA_real_
        vals
      },
      FUN.VALUE = base::numeric(base::nrow(mat))
    )
    out <- out %>% base::as.matrix()
    base::rownames(out) <- base::rownames(mat)
    out
  }
  .hc_ui_get_integrated_net_genes <- function() {
    merged_net <- hcobject[["integrated_output"]][["merged_net"]]
    if (base::is.null(merged_net)) {
      stop("Missing integrated network. Run `hc_build_integrated_network()` first.")
    }
    net_genes <- igraph::get.vertex.attribute(merged_net, "name")
    net_genes <- base::unique(base::as.character(net_genes))
    net_genes <- net_genes[!base::is.na(net_genes) & net_genes != ""]
    if (base::length(net_genes) == 0) {
      stop("Could not extract genes from the integrated network.")
    }
    net_genes
  }
  .hc_ui_prepare_group_mean_matrix <- function(set_name) {
    expr_mat <- hcobject[["layer_specific_outputs"]][[set_name]][["part1"]][["topvar"]]
    if (base::is.null(expr_mat)) {
      expr_mat <- hcobject[["data"]][[base::paste0(set_name, "_counts")]]
    }
    anno <- hcobject[["data"]][[base::paste0(set_name, "_anno")]]
    if (base::is.null(expr_mat) || base::is.null(anno)) {
      return(NULL)
    }
    expr_mat <- expr_mat %>% base::as.matrix()
    if (base::nrow(expr_mat) == 0 || base::ncol(expr_mat) == 0) {
      return(NULL)
    }
    mode(expr_mat) <- "numeric"
    samples <- base::intersect(base::colnames(expr_mat), base::rownames(anno))
    if (base::length(samples) == 0) {
      return(NULL)
    }
    expr_mat <- expr_mat[, samples, drop = FALSE]
    anno <- anno[samples, , drop = FALSE]
    grpvar <- .hc_ui_resolve_group_labels(anno)
    grpvar <- base::as.character(grpvar)
    if (base::length(grpvar) != base::length(samples)) {
      return(NULL)
    }
    bad_grp <- base::is.na(grpvar) | grpvar == ""
    if (base::any(bad_grp)) {
      grpvar[bad_grp] <- samples[bad_grp]
    }
    if (isTRUE(hcobject[["global_settings"]][["data_in_log"]])) {
      expr_mat <- .hc_antilog_impl(expr_mat, 2)
    }
    grp_levels <- base::unique(grpvar)
    set_mean_mat <- base::vapply(
      grp_levels,
      function(grp) {
        idx <- base::which(grpvar == grp)
        if (base::length(idx) == 1) {
          return(expr_mat[, idx])
        }
        base::rowMeans(expr_mat[, idx, drop = FALSE], na.rm = TRUE)
      },
      FUN.VALUE = base::numeric(base::nrow(expr_mat))
    )
    set_mean_mat <- set_mean_mat %>% base::as.matrix()
    if (base::ncol(set_mean_mat) == 0) {
      return(NULL)
    }
    base::colnames(set_mean_mat) <- grp_levels
    base::rownames(set_mean_mat) <- base::rownames(expr_mat)
    set_mean_mat
  }
  .hc_ui_parse_fc_comparisons <- function(fc_comparisons) {
    if (base::is.null(fc_comparisons) || base::length(fc_comparisons) == 0) {
      stop(
        "`fc_comparisons` must be provided when `activity_input = 'fc'`. ",
        "Use entries like 'groupA_vs_groupB'."
      )
    }
    parse_one <- function(x) {
      x <- base::trimws(base::as.character(x))
      if (base::is.na(x) || x == "") {
        return(NULL)
      }
      parts <- base::strsplit(x, "\\s*_vs_\\s*", perl = TRUE)[[1]]
      if (base::length(parts) != 2) {
        parts <- base::strsplit(x, "\\s+vs\\s+", perl = TRUE)[[1]]
      }
      if (base::length(parts) != 2) {
        stop(
          "Invalid `fc_comparisons` entry: '", x, "'. ",
          "Use format 'groupA_vs_groupB'."
        )
      }
      num <- base::trimws(parts[[1]])
      den <- base::trimws(parts[[2]])
      if (num == "" || den == "") {
        stop(
          "Invalid `fc_comparisons` entry: '", x, "'. ",
          "Both groups must be non-empty."
        )
      }
      base::data.frame(
        comparison = base::paste0(num, "_vs_", den),
        numerator = num,
        denominator = den,
        stringsAsFactors = FALSE
      )
    }
    out <- base::lapply(fc_comparisons, parse_one)
    out <- out[!base::vapply(out, base::is.null, FUN.VALUE = base::logical(1))]
    if (base::length(out) == 0) {
      stop("No valid entries in `fc_comparisons`.")
    }
    out <- base::do.call(base::rbind, out)
    out <- out[!duplicated(out$comparison), , drop = FALSE]
    base::rownames(out) <- NULL
    out
  }
  .hc_ui_prepare_activity_matrix_from_gfc <- function(gfc_df) {
    gene_col <- if ("Gene" %in% base::colnames(gfc_df)) "Gene" else base::colnames(gfc_df)[base::ncol(gfc_df)]
    # By index, not by name: `GFC_all_layers` repeats its condition columns
    # when the layers share group names, and setdiff() on the names would
    # collapse those duplicates - silently dropping every layer after the
    # first.
    value_idx <- base::setdiff(base::seq_len(base::ncol(gfc_df)),
                               base::which(base::colnames(gfc_df) %in% gene_col))
    value_cols <- base::colnames(gfc_df)[value_idx]
    mat <- gfc_df[, value_idx, drop = FALSE] %>% base::as.matrix()
    base::colnames(mat) <- value_cols
    mode(mat) <- "numeric"
    base::rownames(mat) <- base::as.character(gfc_df[[gene_col]])
    keep_rows <- base::rowSums(!base::is.na(mat)) > 0
    mat <- mat[keep_rows, , drop = FALSE]
    if (base::nrow(mat) == 0) {
      stop("No non-missing rows available in GFC matrix.")
    }
    mat <- .hc_ui_collapse_duplicate_rows(mat)
    mat
  }
  .hc_ui_resolve_group_labels <- function(info_dataset) {
    info_dataset <- info_dataset %>% base::as.data.frame(stringsAsFactors = FALSE)
    if (base::nrow(info_dataset) == 0) {
      return(base::character(0))
    }
    voi <- hcobject[["global_settings"]][["voi"]]
    voi <- base::intersect(voi, base::colnames(info_dataset))
    if (base::length(voi) > 0) {
      return(
        purrr::pmap(info_dataset[, voi, drop = FALSE], paste, sep = "-") %>%
          base::unlist()
      )
    }
    base::as.character(info_dataset[[1]])
  }
  .hc_ui_prepare_activity_matrix_from_fc <- function(fc_comparisons) {
    net_genes <- .hc_ui_get_integrated_net_genes()
    comparison_df <- .hc_ui_parse_fc_comparisons(fc_comparisons)
    comparison_labels <- base::as.character(comparison_df$comparison)
    comparison_found <- stats::setNames(base::rep(FALSE, base::length(comparison_labels)), comparison_labels)
    available_groups <- base::character(0)
    pseudo_count <- 1e-08
    set_indices <- base::seq_len(base::length(hcobject[["layer_specific_outputs"]]))
    set_mats <- list()
    for (z in set_indices) {
      set_name <- base::paste0("set", z)
      set_mean_mat <- .hc_ui_prepare_group_mean_matrix(set_name = set_name)
      if (base::is.null(set_mean_mat)) {
        next
      }
      available_groups <- base::union(available_groups, base::colnames(set_mean_mat))
      overlap <- base::intersect(net_genes, base::rownames(set_mean_mat))
      if (base::length(overlap) == 0) {
        next
      }
      for (i in base::seq_len(base::nrow(comparison_df))) {
        num_grp <- base::as.character(comparison_df$numerator[i])
        den_grp <- base::as.character(comparison_df$denominator[i])
        cmp_label <- base::as.character(comparison_df$comparison[i])
        if (!(num_grp %in% base::colnames(set_mean_mat) && den_grp %in% base::colnames(set_mean_mat))) {
          next
        }
        comparison_found[[cmp_label]] <- TRUE
        num_vals <- .hc_as_numeric_safely(set_mean_mat[overlap, num_grp])
        den_vals <- .hc_as_numeric_safely(set_mean_mat[overlap, den_grp])
        fc_vals <- .hc_log2_safely((num_vals + pseudo_count) / (den_vals + pseudo_count))
        fc_vals[!base::is.finite(fc_vals)] <- NA_real_
        set_full <- base::matrix(
          NA_real_,
          nrow = base::length(net_genes),
          ncol = 1,
          dimnames = list(net_genes, cmp_label)
        )
        set_full[overlap, 1] <- fc_vals
        set_mats[[base::length(set_mats) + 1]] <- set_full
      }
    }
    found_labels <- base::names(comparison_found)[comparison_found]
    missing_labels <- base::names(comparison_found)[!comparison_found]
    if (base::length(found_labels) == 0 || base::length(set_mats) == 0) {
      stop(
        "None of the requested `fc_comparisons` could be computed from current groups. ",
        "Requested: ", base::paste(comparison_labels, collapse = ", "), ". ",
        "Available groups: ", base::paste(base::sort(base::unique(available_groups)), collapse = ", "), "."
      )
    }
    if (base::length(missing_labels) > 0) {
      warning(
        "Some `fc_comparisons` could not be computed and will be ignored: ",
        base::paste(missing_labels, collapse = ", ")
      )
    }
    mat <- base::do.call(base::cbind, set_mats)
    mat <- .hc_ui_collapse_duplicate_rows(mat)
    mat <- .hc_ui_collapse_duplicate_columns(mat)
    keep_rows <- base::rowSums(!base::is.na(mat)) > 0
    mat <- mat[keep_rows, , drop = FALSE]
    if (base::nrow(mat) == 0 || base::ncol(mat) == 0) {
      stop("No non-missing rows/columns available in FC activity matrix.")
    }
    ordered_cols <- comparison_labels[comparison_labels %in% base::colnames(mat)]
    ordered_cols <- base::c(ordered_cols, base::setdiff(base::colnames(mat), ordered_cols))
    mat <- mat[, ordered_cols, drop = FALSE]
    list(
      mat = mat,
      requested = comparison_labels,
      used = ordered_cols,
      missing = missing_labels
    )
  }
  .hc_ui_module_means_from_gene_matrix <- function(gene_mat, cluster_info, cluster_order) {
    out_list <- list()
    for (cl in cluster_order) {
      genes <- dplyr::filter(cluster_info, color == cl) %>%
        dplyr::pull(., "gene_n") %>%
        base::strsplit(split = ",") %>%
        base::unlist()
      genes <- base::unique(base::as.character(genes))
      genes <- base::intersect(genes, base::rownames(gene_mat))
      if (base::length(genes) == 0) {
        next
      }
      vals <- gene_mat[genes, , drop = FALSE]
      out <- base::colMeans(vals, na.rm = TRUE)
      out[base::is.nan(out)] <- NA_real_
      out_list[[cl]] <- out
    }
    if (base::length(out_list) == 0) {
      return(NULL)
    }
    out_mat <- base::do.call(base::rbind, out_list)
    out_mat %>% base::as.matrix()
  }

  fc_summary <- NULL
  activity_mat <- NULL
  activity_label <- "GFC"
  module_heatmap_mat <- NULL
  module_heatmap_col_order <- NULL
  module_heatmap_name <- "GFC"
  activity_module_heatmap_mat <- NULL
  if (identical(activity_input, "fc")) {
    fc_summary <- .hc_ui_prepare_activity_matrix_from_fc(fc_comparisons = fc_comparisons)
    activity_mat <- fc_summary$mat
    activity_label <- "FC"
    module_heatmap_mat <- .hc_ui_module_means_from_gene_matrix(
      gene_mat = activity_mat,
      cluster_info = cluster_info,
      cluster_order = cluster_order
    )
    module_heatmap_col_order <- fc_summary$used
    module_heatmap_name <- "FC"
  } else {
    activity_mat <- .hc_ui_prepare_activity_matrix_from_gfc(gfc_all)
    activity_label <- "GFC"
    # Layers that share group names contribute repeated condition columns.
    # Label those with their layer, as the module heatmap does, so that each
    # layer stays a condition of its own instead of being merged by name.
    gfc_cols <- base::colnames(activity_mat)
    if (base::anyDuplicated(gfc_cols) > 0) {
      base::colnames(activity_mat) <- base::make.unique(
        .hc_gfc_display_col_labels(hcobject, gfc_cols),
        sep = " "
      )
    }
  }
  activity_module_heatmap_mat <- .hc_ui_module_means_from_gene_matrix(
    gene_mat = activity_mat,
    cluster_info = cluster_info,
    cluster_order = cluster_order
  )
  condition_levels <- base::colnames(activity_mat)

  message("...upstream inference activity input: ", activity_label, "...")
  organism <- base::tolower(base::as.character(hcobject[["global_settings"]][["organism"]]))
  if (!organism %in% c("human", "mouse")) {
    warning(
      "Unknown organism setting `", organism, "`. ",
      "Falling back to human priors."
    )
    organism <- "human"
  }

  tf_network <- NULL
  tf_database_label <- "CollecTRI"
  pathway_network <- NULL
  pathway_database_label <- "PROGENy"
  custom_pathway_network <- .hc_load_custom_gmt_pathway_network(
    custom_pathway_gmt = custom_pathway_gmt,
    default_prefix = "CustomPathway"
  )
  custom_pathway_paths <- base::attr(custom_pathway_network, "paths")
  custom_pathway_databases <- base::attr(custom_pathway_network, "databases")
  if (base::is.null(custom_pathway_databases)) {
    custom_pathway_databases <- base::character(0)
  }

  # Test universe: the genes that could have landed in a module at all, i.e.
  # the integrated network genes with a value in the activity matrix.
  universe_genes <- base::rownames(activity_mat)
  net_genes_universe <- tryCatch(.hc_ui_get_integrated_net_genes(), error = function(e) NULL)
  if (!base::is.null(net_genes_universe)) {
    universe_genes <- base::intersect(universe_genes, net_genes_universe)
  }
  if (base::length(universe_genes) == 0) {
    stop("No genes left for the test universe (network genes with ", activity_label, " values).")
  }
  module_genes <- stats::setNames(base::lapply(cluster_order, function(cl) {
    genes <- dplyr::filter(cluster_info, color == cl) %>%
      dplyr::pull(., "gene_n") %>%
      base::strsplit(split = ",") %>%
      base::unlist()
    base::intersect(base::unique(base::as.character(genes)), universe_genes)
  }), cluster_order)

  if ("TF" %in% resources) {
    tf_network <- .hc_ui_load_tf_network(
      organism = organism,
      tf_resource = tf_resource,
      tf_confidence = tf_confidence
    )
    tf_database_label <- base::attr(tf_network, "database")
    tf_network <- tf_network[tf_network$target %in% universe_genes, , drop = FALSE]
    if (base::nrow(tf_network) == 0) {
      warning("TF network contains no targets overlapping with the network genes.")
    }
  }
  if ("PATHWAY" %in% resources) {
    progeny_network <- tryCatch(
      .hc_ui_load_pathway_network(organism = organism, top = progeny_top),
      error = function(e) {
        warning("Could not load PROGENy pathway model: ", conditionMessage(e))
        NULL
      }
    )
    if (!base::is.null(progeny_network) && base::nrow(progeny_network) > 0) {
      progeny_network <- progeny_network[progeny_network$target %in% universe_genes, , drop = FALSE]
      if (base::nrow(progeny_network) == 0) {
        progeny_network <- NULL
      }
    }

    if (!base::is.null(custom_pathway_network) && base::nrow(custom_pathway_network) > 0) {
      custom_pathway_network <- custom_pathway_network[
        custom_pathway_network$target %in% universe_genes, ,
        drop = FALSE
      ]
      if (base::nrow(custom_pathway_network) == 0) {
        custom_pathway_network <- NULL
      }
    }

    pathway_parts <- list()
    if (!base::is.null(progeny_network) && base::nrow(progeny_network) > 0) {
      pathway_parts[[base::length(pathway_parts) + 1]] <- progeny_network
    }
    if (!base::is.null(custom_pathway_network) && base::nrow(custom_pathway_network) > 0) {
      pathway_parts[[base::length(pathway_parts) + 1]] <- custom_pathway_network
    }
    if (base::length(pathway_parts) > 0) {
      pathway_network <- base::do.call(base::rbind, pathway_parts)
      pathway_network <- pathway_network %>% base::as.data.frame(stringsAsFactors = FALSE)
      pathway_network <- base::unique(pathway_network)
      base::rownames(pathway_network) <- NULL
    } else {
      pathway_network <- base::data.frame(
        source = base::character(0),
        target = base::character(0),
        mor = base::numeric(0),
        stringsAsFactors = FALSE
      )
    }

    if (base::nrow(pathway_network) == 0) {
      warning("Pathway network contains no targets overlapping with the network genes.")
    } else {
      has_progeny <- !base::is.null(progeny_network) && base::nrow(progeny_network) > 0
      has_custom <- !base::is.null(custom_pathway_network) && base::nrow(custom_pathway_network) > 0
      pathway_database_label <- if (has_progeny && has_custom) {
        "PROGENy+CustomGMT"
      } else if (has_custom) {
        "CustomGMT"
      } else {
        "PROGENy"
      }
    }
  }

  summary_cols <- c(
    "resource",
    "database",
    "cluster",
    "module_label",
    "rank",
    "term",
    "score",
    "abs_score",
    "pvalue",
    "qvalue",
    "direction",
    "regulation",
    "n_conditions",
    "n_genes",
    "n_targets",
    "n_overlap",
    "fold_enrichment",
    "peak_condition",
    "consistency",
    "n_active_conditions",
    "regulator_module",
    "regulator_in_module",
    "regulator_cor",
    "redundant_with",
    "overlap_genes"
  )
  summary_cols_condition <- base::append(
    base::append(summary_cols, "condition", after = 4),
    c("activity", "activity_pvalue", "activity_qvalue"),
    after = 9
  )
  empty_summary <- function(with_condition = FALSE) {
    cols <- if (with_condition) summary_cols_condition else summary_cols
    num_cols <- c("rank", "score", "abs_score", "pvalue", "qvalue", "n_conditions",
                  "n_genes", "n_targets", "n_overlap", "fold_enrichment", "consistency",
                  "n_active_conditions", "regulator_cor", "activity", "activity_pvalue",
                  "activity_qvalue")
    out <- base::lapply(cols, function(nm) {
      if (nm %in% num_cols) {
        base::numeric(0)
      } else if (nm == "regulator_in_module") {
        base::logical(0)
      } else {
        base::character(0)
      }
    })
    base::names(out) <- cols
    base::as.data.frame(out, stringsAsFactors = FALSE)
  }
  empty_summary_condition <- function() empty_summary(with_condition = TRUE)

  # Module of every network gene (by label), used to check whether a TF is
  # itself co-expressed with its targets.
  gene_module_label <- base::unlist(base::lapply(cluster_order, function(cl) {
    stats::setNames(
      base::rep(module_label_map_current[[cl]], base::length(module_genes[[cl]])),
      module_genes[[cl]]
    )
  }))

  run_one_resource <- function(resource_name, database_name, network_df) {
    empty_out <- list(
      selected = empty_summary(),
      significant = empty_summary(),
      all_by_condition = empty_summary_condition(),
      selected_by_condition = empty_summary_condition(),
      significant_by_condition = empty_summary_condition(),
      raw = base::data.frame(),
      activity = base::data.frame()
    )
    if (base::is.null(network_df) || base::nrow(network_df) == 0) {
      return(empty_out)
    }
    ora <- .hc_ui_module_ora(
      module_genes = module_genes,
      network = network_df,
      universe = universe_genes,
      minsize = minsize,
      padj = padj
    )
    if (base::nrow(ora) == 0) {
      return(empty_out)
    }
    target_stats <- .hc_ui_module_target_scores(ora = ora, network = network_df, values = activity_mat)
    score_mat <- target_stats$score
    consistency_mat <- target_stats$consistency
    peak_idx <- base::apply(base::abs(score_mat), 1, function(v) {
      if (base::all(base::is.na(v))) NA_integer_ else base::which.max(v)
    })
    peak_cells <- base::cbind(base::seq_len(base::nrow(ora)), peak_idx)
    ora$score <- score_mat[peak_cells]
    ora$abs_score <- base::abs(ora$score)
    ora$peak_condition <- base::colnames(score_mat)[peak_idx]
    ora$consistency <- consistency_mat[peak_cells]
    ora$regulation <- .hc_ui_module_regulation(ora = ora, network = network_df)
    ora$n_conditions <- base::rowSums(!base::is.na(score_mat))
    ora$direction <- base::ifelse(ora$score >= 0, "activated", "inhibited")
    ora$resource <- resource_name
    ora$database <- database_name
    ora$module_label <- base::unname(module_label_map_current[ora$cluster])

    # Is the regulator itself a module gene, and does its profile follow the
    # module? Only meaningful for TFs (pathways are not genes).
    if (identical(resource_name, "TF")) {
      ora$regulator_module <- base::unname(gene_module_label[ora$term])
      ora$regulator_module[base::is.na(ora$regulator_module)] <- ""
      ora$regulator_in_module <- ora$regulator_module == ora$module_label
      ora$regulator_cor <- base::vapply(base::seq_len(base::nrow(ora)), function(i) {
        .hc_ui_regulator_cor(
          regulator = ora$term[i],
          cluster = ora$cluster[i],
          values = activity_mat,
          module_means = activity_module_heatmap_mat
        )
      }, base::numeric(1))
    } else {
      ora$regulator_module <- NA_character_
      ora$regulator_in_module <- NA
      ora$regulator_cor <- NA_real_
    }

    # Regulator activity per condition over all universe genes (decoupleR ULM,
    # the use it is designed for): when is the regulator active?
    activity <- .hc_ui_condition_activity(
      values = activity_mat[base::rownames(activity_mat) %in% universe_genes, , drop = FALSE],
      network = network_df,
      minsize = minsize,
      padj = padj
    )
    activity_key <- base::paste(activity$source, activity$condition, sep = "\t")
    active_mat <- base::vapply(condition_levels, function(cond_nm) {
      idx <- base::match(base::paste(ora$term, cond_nm, sep = "\t"), activity_key)
      act <- activity$activity[idx]
      q <- activity$activity_qvalue[idx]
      sc <- score_mat[, cond_nm]
      !base::is.na(q) & q <= qval & !base::is.na(sc) & base::sign(act) == base::sign(sc)
    }, base::logical(base::nrow(ora)))
    active_mat <- base::matrix(active_mat, nrow = base::nrow(ora), dimnames = list(NULL, condition_levels))

    row_order <- base::order(
      base::factor(ora$cluster, levels = cluster_order),
      ora$qvalue,
      -ora$fold_enrichment,
      ora$term
    )
    ora <- ora[row_order, , drop = FALSE]
    score_mat <- score_mat[row_order, , drop = FALSE]
    active_mat <- active_mat[row_order, , drop = FALSE]
    ora$rank <- stats::ave(base::seq_len(base::nrow(ora)), ora$cluster, FUN = base::seq_along)
    base::rownames(ora) <- NULL

    is_sig <- !base::is.na(ora$qvalue) & ora$qvalue <= qval & ora$n_overlap >= min_overlap
    active_mat[!is_sig, ] <- FALSE
    ora$n_active_conditions <- base::rowSums(active_mat)
    ora$redundant_with <- .hc_ui_flag_redundant(ora, keep = is_sig, jaccard = redundancy_jaccard)

    significant <- ora[is_sig, summary_cols, drop = FALSE]
    if (base::nrow(significant) > 0) {
      significant$rank <- stats::ave(base::seq_len(base::nrow(significant)), significant$cluster, FUN = base::seq_along)
    }
    selectable <- if (isTRUE(collapse_redundant)) {
      significant[significant$redundant_with == "", , drop = FALSE]
    } else {
      significant
    }
    if (base::nrow(selectable) > 0) {
      selectable$rank <- stats::ave(base::seq_len(base::nrow(selectable)), selectable$cluster, FUN = base::seq_along)
    }
    selected <- selectable[selectable$rank <= top, , drop = FALSE]

    # One row per condition: same module link (p/q), condition-specific
    # target score and genome-wide regulator activity.
    ora_key <- base::paste(ora$cluster, ora$term, sep = "\t")
    expand_by_condition <- function(df, active_only = FALSE) {
      if (base::nrow(df) == 0) {
        return(empty_summary_condition())
      }
      rows <- base::match(base::paste(df$cluster, df$term, sep = "\t"), ora_key)
      out <- base::lapply(condition_levels, function(cond_nm) {
        x <- df
        x$condition <- cond_nm
        x$score <- score_mat[rows, cond_nm]
        x$abs_score <- base::abs(x$score)
        x$direction <- base::ifelse(x$score >= 0, "activated", "inhibited")
        x$n_conditions <- 1L
        idx <- base::match(base::paste(x$term, cond_nm, sep = "\t"), activity_key)
        x$activity <- activity$activity[idx]
        x$activity_pvalue <- activity$activity_pvalue[idx]
        x$activity_qvalue <- activity$activity_qvalue[idx]
        keep <- !base::is.na(x$score)
        if (isTRUE(active_only)) {
          keep <- keep & active_mat[rows, cond_nm]
        }
        x[keep, summary_cols_condition, drop = FALSE]
      })
      out <- base::do.call(base::rbind, out)
      out <- out[base::order(
        base::factor(out$cluster, levels = cluster_order),
        base::factor(out$condition, levels = condition_levels),
        out$rank
      ), , drop = FALSE]
      base::rownames(out) <- NULL
      out
    }

    raw <- base::cbind(ora[, summary_cols, drop = FALSE], score_mat)
    base::colnames(raw)[base::seq_len(base::ncol(score_mat)) + base::length(summary_cols)] <-
      base::paste0("score_", base::colnames(score_mat))

    list(
      selected = selected,
      significant = significant,
      all_by_condition = expand_by_condition(ora[, summary_cols, drop = FALSE]),
      selected_by_condition = expand_by_condition(selected),
      significant_by_condition = expand_by_condition(significant, active_only = TRUE),
      raw = raw,
      activity = activity
    )
  }

  tf_out <- run_one_resource(
    resource_name = "TF",
    database_name = tf_database_label,
    network_df = tf_network
  )
  pathway_out <- run_one_resource(
    resource_name = "Pathway",
    database_name = pathway_database_label,
    network_df = pathway_network
  )

  combine_rows <- function(...) {
    lst <- list(...)
    lst <- lst[base::vapply(lst, function(x) !base::is.null(x) && base::nrow(x) > 0, FUN.VALUE = base::logical(1))]
    if (base::length(lst) == 0) {
      return(empty_summary())
    }
    out <- base::do.call(base::rbind, lst)
    out <- out[base::order(
      base::factor(out$cluster, levels = cluster_order),
      out$resource,
      out$rank
    ), , drop = FALSE]
    base::rownames(out) <- NULL
    out
  }
  combine_rows_by_condition <- function(...) {
    lst <- list(...)
    lst <- lst[base::vapply(lst, function(x) !base::is.null(x) && base::nrow(x) > 0, FUN.VALUE = base::logical(1))]
    if (base::length(lst) == 0) {
      return(empty_summary_condition())
    }
    out <- base::do.call(base::rbind, lst)
    out <- out[base::order(
      base::factor(out$cluster, levels = cluster_order),
      base::factor(out$condition, levels = condition_levels),
      out$resource,
      out$rank
    ), , drop = FALSE]
    base::rownames(out) <- NULL
    out
  }

  selected_all <- combine_rows(tf_out$selected, pathway_out$selected)
  significant_all <- combine_rows(tf_out$significant, pathway_out$significant)
  all_by_condition_all <- combine_rows_by_condition(tf_out$all_by_condition, pathway_out$all_by_condition)
  selected_by_condition_all <- combine_rows_by_condition(tf_out$selected_by_condition, pathway_out$selected_by_condition)
  significant_by_condition_all <- combine_rows_by_condition(tf_out$significant_by_condition, pathway_out$significant_by_condition)
  if (base::nrow(selected_all) > 0) {
    selected_all$term_with_resource <- base::paste0(selected_all$term, " [", selected_all$resource, "]")
  } else {
    selected_all$term_with_resource <- base::character(0)
  }
  if (base::nrow(significant_all) > 0) {
    significant_all$term_with_resource <- base::paste0(significant_all$term, " [", significant_all$resource, "]")
  } else {
    significant_all$term_with_resource <- base::character(0)
  }
  if (base::nrow(selected_by_condition_all) > 0) {
    selected_by_condition_all$term_with_resource <- base::paste0(selected_by_condition_all$term, " [", selected_by_condition_all$resource, "]")
  } else {
    selected_by_condition_all$term_with_resource <- base::character(0)
  }
  if (base::nrow(significant_by_condition_all) > 0) {
    significant_by_condition_all$term_with_resource <- base::paste0(significant_by_condition_all$term, " [", significant_by_condition_all$resource, "]")
  } else {
    significant_by_condition_all$term_with_resource <- base::character(0)
  }
  if (base::nrow(all_by_condition_all) > 0) {
    all_by_condition_all$term_with_resource <- base::paste0(all_by_condition_all$term, " [", all_by_condition_all$resource, "]")
  } else {
    all_by_condition_all$term_with_resource <- base::character(0)
  }
  term_order_all <- .hc_ui_upstream_term_order(
    selected_all = selected_all,
    cluster_levels = cluster_order
  )
  term_order_by_condition_all <- .hc_ui_upstream_term_order(
    selected_all = selected_by_condition_all,
    cluster_levels = cluster_order
  )
  term_order_tf <- term_order_all[grepl("\\[TF\\]$", term_order_all)]
  term_order_pathway <- term_order_all[grepl("\\[Pathway\\]$", term_order_all)]
  all_scores_for_scale <- .hc_as_numeric_safely(base::c(
    selected_all$score,
    significant_all$score,
    selected_by_condition_all$score,
    significant_by_condition_all$score
  ))
  all_scores_for_scale <- all_scores_for_scale[base::is.finite(all_scores_for_scale)]
  if (base::length(all_scores_for_scale) == 0) {
    activity_score_limit <- 2
  } else {
    activity_score_limit <- .hc_first_numeric_value(
      stats::quantile(
        base::abs(all_scores_for_scale),
        probs = 0.98,
        na.rm = TRUE
      )
    )
    if (!base::is.finite(activity_score_limit) || base::is.na(activity_score_limit)) {
      activity_score_limit <- 2
    }
  }
  activity_score_limit <- base::max(1, base::min(4, activity_score_limit))

  ora_tf <- tf_out$raw
  ora_pathway <- pathway_out$raw

  file_prefix <- base::paste0(
    hcobject[["working_directory"]][["dir_output"]],
    hcobject[["global_settings"]][["save_folder"]]
  )
  excel_path <- base::paste0(file_prefix, "/Upstream_Inference.xlsx")
  pdf_path <- base::paste0(file_prefix, "/Upstream_Inference.pdf")

  .hc_ui_excel_safe_sheet_names <- function(nms) {
    if (base::is.null(nms) || base::length(nms) == 0) {
      return(base::character(0))
    }
    out <- base::as.character(nms)
    out <- gsub("[\\[\\]\\*\\?/\\\\:]", "_", out)
    out[out == "" | base::is.na(out)] <- "sheet"
    out <- base::substr(out, 1, 31)
    used <- base::character(0)
    for (i in base::seq_along(out)) {
      nm <- out[[i]]
      if (!(nm %in% used)) {
        used <- base::c(used, nm)
        next
      }
      base_nm <- nm
      k <- 2L
      repeat {
        suffix <- base::paste0("_", k)
        keep <- 31 - base::nchar(suffix)
        if (keep < 1) {
          keep <- 1
        }
        cand <- base::paste0(base::substr(base_nm, 1, keep), suffix)
        if (!(cand %in% used)) {
          out[[i]] <- cand
          used <- base::c(used, cand)
          break
        }
        k <- k + 1L
      }
    }
    out
  }

  export_tables <- list(
    selected_upstream_all = selected_all,
    significant_upstream_all = significant_all,
    all_upstream_by_cond = all_by_condition_all,
    selected_upstream_by_cond = selected_by_condition_all,
    significant_upstream_by_cond = significant_by_condition_all,
    selected_tf = tf_out$selected,
    significant_tf = tf_out$significant,
    selected_pathway = pathway_out$selected,
    significant_pathway = pathway_out$significant,
    ora_tf = ora_tf,
    ora_pathway = ora_pathway,
    activity_tf = tf_out$activity,
    activity_pathway = pathway_out$activity
  )
  base::names(export_tables) <- .hc_ui_excel_safe_sheet_names(base::names(export_tables))
  .hc_write_xlsx_atomic(
    x = export_tables,
    file = excel_path,
    overwrite = TRUE
  )

  dot_plot <- .hc_ui_build_upstream_dotplot(
    selected_all = selected_all,
    cluster_order = cluster_order,
    module_label_map = module_label_map_current,
    term_levels = term_order_all,
    overall_plot_scale = overall_plot_scale
  )
  heatmap_plot <- .hc_ui_build_upstream_heatmap(
    significant_all = significant_all,
    selected_all = selected_all,
    cluster_order = cluster_order,
    module_label_map = module_label_map_current,
    term_levels = term_order_all,
    activity_score_limit = activity_score_limit,
    overall_plot_scale = overall_plot_scale
  )
  combined_plot <- .hc_ui_build_upstream_combined_heatmap(
    selected_all = selected_all,
    significant_all = significant_all,
    cluster_order = cluster_order,
    module_label_map = module_label_map_current,
    cluster_info = cluster_info,
    gfc_all = gfc_all,
    stored_hm = stored_hm,
    main_heatmap_col_order = main_heatmap_col_order,
    heatmap_cluster_columns = cluster_columns,
    heatmap_col_order = col_order,
    module_heatmap_mat = module_heatmap_mat,
    module_heatmap_col_order = module_heatmap_col_order,
    module_heatmap_name = module_heatmap_name,
    gfc_scale_limits = gfc_scale_limits,
    gfc_scale_ticks = gfc_scale_ticks,
    heatmap_side = heatmap_side,
    term_levels = term_order_all,
    activity_score_limit = activity_score_limit,
    show_horizontal_lines = FALSE,
    overall_plot_scale = overall_plot_scale
  )
  selected_tf_only <- selected_all[selected_all$resource == "TF", , drop = FALSE]
  significant_tf_only <- significant_all[significant_all$resource == "TF", , drop = FALSE]
  selected_pathway_only <- selected_all[selected_all$resource == "Pathway", , drop = FALSE]
  significant_pathway_only <- significant_all[significant_all$resource == "Pathway", , drop = FALSE]
  combined_plot_tf <- .hc_ui_build_upstream_combined_heatmap(
    selected_all = selected_tf_only,
    significant_all = significant_tf_only,
    cluster_order = cluster_order,
    module_label_map = module_label_map_current,
    cluster_info = cluster_info,
    gfc_all = gfc_all,
    stored_hm = stored_hm,
    main_heatmap_col_order = main_heatmap_col_order,
    heatmap_cluster_columns = cluster_columns,
    heatmap_col_order = col_order,
    module_heatmap_mat = module_heatmap_mat,
    module_heatmap_col_order = module_heatmap_col_order,
    module_heatmap_name = module_heatmap_name,
    gfc_scale_limits = gfc_scale_limits,
    gfc_scale_ticks = gfc_scale_ticks,
    heatmap_side = heatmap_side,
    term_levels = term_order_tf,
    activity_score_limit = activity_score_limit,
    show_horizontal_lines = FALSE,
    overall_plot_scale = overall_plot_scale
  )
  combined_plot_pathway <- .hc_ui_build_upstream_combined_heatmap(
    selected_all = selected_pathway_only,
    significant_all = significant_pathway_only,
    cluster_order = cluster_order,
    module_label_map = module_label_map_current,
    cluster_info = cluster_info,
    gfc_all = gfc_all,
    stored_hm = stored_hm,
    main_heatmap_col_order = main_heatmap_col_order,
    heatmap_cluster_columns = cluster_columns,
    heatmap_col_order = col_order,
    module_heatmap_mat = module_heatmap_mat,
    module_heatmap_col_order = module_heatmap_col_order,
    module_heatmap_name = module_heatmap_name,
    gfc_scale_limits = gfc_scale_limits,
    gfc_scale_ticks = gfc_scale_ticks,
    heatmap_side = heatmap_side,
    term_levels = term_order_pathway,
    activity_score_limit = activity_score_limit,
    show_horizontal_lines = FALSE,
    overall_plot_scale = overall_plot_scale
  )

  combined_plots_by_condition <- list()
  if (isTRUE(plot_per_comparison) &&
    base::nrow(selected_by_condition_all) > 0 &&
    base::length(condition_levels) > 0) {
    for (cond_nm in condition_levels) {
      selected_cond <- selected_by_condition_all[selected_by_condition_all$condition == cond_nm, , drop = FALSE]
      if (base::nrow(selected_cond) == 0) {
        next
      }
      significant_cond <- significant_by_condition_all[significant_by_condition_all$condition == cond_nm, , drop = FALSE]
      selected_cond_plot <- selected_cond
      mark_sig_cond <- FALSE
      significant_marks <- NULL
      value_from_significant_cond <- FALSE
      if (isTRUE(consistent_terms)) {
        mark_sig_cond <- TRUE
        significant_marks <- significant_cond
        value_from_significant_cond <- FALSE
      }
      term_levels_cond <- if (isTRUE(consistent_terms) &&
        (base::length(term_order_by_condition_all) > 0 || base::length(term_order_all) > 0)) {
        if (base::length(term_order_by_condition_all) > 0) {
          term_order_by_condition_all
        } else {
          term_order_all
        }
      } else {
        .hc_ui_upstream_term_order(
          selected_all = selected_cond,
          cluster_levels = cluster_order
        )
      }
      condition_module_mat <- NULL
      condition_module_col_order <- NULL
      if (!base::is.null(activity_module_heatmap_mat) &&
        base::nrow(activity_module_heatmap_mat) > 0 &&
        cond_nm %in% base::colnames(activity_module_heatmap_mat)) {
        condition_module_mat <- activity_module_heatmap_mat[, cond_nm, drop = FALSE]
        condition_module_col_order <- cond_nm
      }
      condition_plot <- .hc_ui_build_upstream_combined_heatmap(
        selected_all = selected_cond_plot,
        significant_all = significant_cond,
        cluster_order = cluster_order,
        module_label_map = module_label_map_current,
        cluster_info = cluster_info,
        gfc_all = gfc_all,
        stored_hm = stored_hm,
        main_heatmap_col_order = main_heatmap_col_order,
        heatmap_cluster_columns = cluster_columns,
        heatmap_col_order = col_order,
        module_heatmap_mat = condition_module_mat,
        module_heatmap_col_order = condition_module_col_order,
        module_heatmap_name = module_heatmap_name,
        gfc_scale_limits = gfc_scale_limits,
        gfc_scale_ticks = gfc_scale_ticks,
        heatmap_side = heatmap_side,
        term_levels = term_levels_cond,
        activity_score_limit = activity_score_limit,
        show_horizontal_lines = FALSE,
        mark_significant = mark_sig_cond,
        significant_for_marks = significant_marks,
        value_from_significant = value_from_significant_cond,
        overall_plot_scale = overall_plot_scale
      )
      if (!base::is.null(condition_plot)) {
        combined_plots_by_condition[[cond_nm]] <- condition_plot
      }
    }
  }

  draw_saved_page <- function(page_key) {
    page_key <- base::as.character(page_key[[1]])
    if (isTRUE(plot_per_comparison) && base::length(combined_plots_by_condition) > 0) {
      ComplexHeatmap::draw(
        combined_plots_by_condition[[page_key]],
        merge_legends = TRUE,
        newpage = TRUE
      )
      return(invisible(NULL))
    }
    if (identical(page_key, "combined") && !base::is.null(combined_plot)) {
      ComplexHeatmap::draw(combined_plot, merge_legends = TRUE, newpage = TRUE)
      return(invisible(NULL))
    }
    if (identical(page_key, "combined_tf") && !base::is.null(combined_plot_tf)) {
      ComplexHeatmap::draw(combined_plot_tf, merge_legends = TRUE, newpage = TRUE)
      return(invisible(NULL))
    }
    if (identical(page_key, "combined_pathway") && !base::is.null(combined_plot_pathway)) {
      ComplexHeatmap::draw(combined_plot_pathway, merge_legends = TRUE, newpage = TRUE)
      return(invisible(NULL))
    }
    if (identical(page_key, "dotplot") && !base::is.null(dot_plot)) {
      .hc_display_object(dot_plot)
      return(invisible(NULL))
    }
    if (identical(page_key, "score_heatmap") && !base::is.null(heatmap_plot)) {
      ComplexHeatmap::draw(heatmap_plot, merge_legends = TRUE, newpage = TRUE)
      return(invisible(NULL))
    }
    invisible(NULL)
  }

  saved_page_labels <- if (isTRUE(plot_per_comparison) && base::length(combined_plots_by_condition) > 0) {
    base::names(combined_plots_by_condition)
  } else if (!base::is.null(combined_plot)) {
    base::c(
      "combined",
      if (!base::is.null(combined_plot_tf)) "combined_tf",
      if (!base::is.null(combined_plot_pathway)) "combined_pathway"
    )
  } else if (!base::is.null(dot_plot)) {
    "dotplot"
  } else if (!base::is.null(heatmap_plot)) {
    "score_heatmap"
  } else {
    base::character(0)
  }

  export_files <- list(
    xlsx = excel_path,
    pdf = NULL,
    png = base::character(0)
  )
  if (isTRUE(save_pdf)) {
    n_cols_for_width <- if (base::nrow(selected_all) > 0) {
      base::length(base::unique(selected_all$term_with_resource))
    } else {
      6
    }
    pdf_width_auto <- base::max(11, base::min(26, 9 + 0.15 * n_cols_for_width))
    pdf_height_auto <- base::max(8, base::min(20, 6 + 0.22 * base::length(cluster_order)))
    pdf_width_use <- if (base::is.null(pdf_width)) pdf_width_auto else as.numeric(pdf_width)
    pdf_height_use <- if (base::is.null(pdf_height)) pdf_height_auto else as.numeric(pdf_height)
    export_files <- .hc_export_multi_page_plot(
      file = pdf_path,
      page_labels = saved_page_labels,
      width = pdf_width_use,
      height = pdf_height_use,
      pointsize = pdf_pointsize,
      res = 300,
      draw_page_fun = function(idx, page_key) {
        draw_saved_page(page_key)
      }
    )
    export_files$xlsx <- excel_path
  }

  if (isTRUE(plot)) {
    if (base::length(saved_page_labels) > 0) {
      for (page_key in saved_page_labels) {
        draw_saved_page(page_key)
      }
    }
  }

  output <- list(
    selected_upstream_all = selected_all,
    significant_upstream_all = significant_all,
    all_upstream_by_condition = all_by_condition_all,
    selected_upstream_by_condition = selected_by_condition_all,
    significant_upstream_by_condition = significant_by_condition_all,
    selected_tf = tf_out$selected,
    significant_tf = tf_out$significant,
    selected_pathway = pathway_out$selected,
    significant_pathway = pathway_out$significant,
    raw_tf = tf_out$raw,
    raw_pathway = pathway_out$raw,
    activity_tf = tf_out$activity,
    activity_pathway = pathway_out$activity,
    module_heatmap_matrix = module_heatmap_mat,
    activity_module_heatmap_matrix = activity_module_heatmap_mat,
    module_heatmap_col_order = module_heatmap_col_order,
    condition_levels = condition_levels,
    module_heatmap_name = module_heatmap_name,
    mixed_dotplot = dot_plot,
    combined_heatmap = combined_plot,
    combined_heatmap_tf = combined_plot_tf,
    combined_heatmap_pathway = combined_plot_pathway,
    combined_heatmaps_by_condition = combined_plots_by_condition,
    score_heatmap = heatmap_plot,
    settings = list(
      activity_input = activity_input,
      fc_comparisons = if (!base::is.null(fc_summary)) fc_summary$used else NULL,
      fc_comparisons_missing = if (!base::is.null(fc_summary)) fc_summary$missing else NULL,
      custom_pathway_gmt = custom_pathway_paths,
      pathway_database = pathway_database_label,
      custom_pathway_databases = custom_pathway_databases,
      plot_per_comparison = plot_per_comparison,
      consistent_terms = consistent_terms,
      resources = resources,
      method = method,
      tf_resource = tf_resource,
      tf_database = if ("TF" %in% resources) tf_database_label else NULL,
      tf_confidence = if (identical(tf_database_label, "DoRothEA")) tf_confidence else NULL,
      progeny_top = progeny_top,
      universe_size = base::length(universe_genes),
      resource_date = base::as.character(base::Sys.Date()),
      minsize = minsize,
      min_overlap = min_overlap,
      collapse_redundant = collapse_redundant,
      redundancy_jaccard = redundancy_jaccard,
      activity_score_limit = activity_score_limit,
      gfc_scale_limits = gfc_scale_limits,
      qval = qval,
      padj = padj
    ),
    files = export_files
  )
  .hc_set_bridge_hcobject_slot(c("satellite_outputs", "upstream_inference"), output)
  output
}



# TF regulons. CollecTRI (recommended by the decoupleR authors) comes from
# OmniPath: through decoupleR/OmnipathR when installed, otherwise (or when
# that fails) directly from the OmniPath web service. DoRothEA comes offline
# from the `dorothea` package. The returned network carries the database
# that was actually used in attr(, "database").
.hc_ui_load_tf_network <- function(organism, tf_resource = "auto", tf_confidence = c("A", "B", "C")) {
  org <- if (identical(organism, "mouse")) "mouse" else "human"
  if (tf_resource %in% c("auto", "collectri")) {
    problem <- "OmniPath could not be reached"
    net <- NULL
    if (requireNamespace("OmnipathR", quietly = TRUE)) {
      net <- tryCatch(
        decoupleR::get_collectri(organism = org, split_complexes = FALSE),
        error = function(e) e
      )
    }
    if (base::is.null(net) || inherits(net, "error") || base::nrow(net) == 0) {
      net <- tryCatch(.hc_ui_fetch_collectri(org), error = function(e) e)
    }
    if (inherits(net, "error")) {
      problem <- conditionMessage(net)
      net <- NULL
    }
    if (!base::is.null(net) && base::nrow(net) > 0) {
      out <- .hc_ui_prepare_network(
        net = net,
        source_candidates = c("source", "tf"),
        target_candidates = c("target"),
        mor_candidates = c("mor", "weight")
      )
      base::attr(out, "database") <- "CollecTRI"
      return(out)
    }
    if (identical(tf_resource, "collectri")) {
      stop("Could not load CollecTRI (", problem, "). Check the internet connection, ",
           "or use `tf_resource = \"dorothea\"`.", call. = FALSE)
    }
    message("CollecTRI not available (", problem, "); using DoRothEA (confidence ",
            base::paste(tf_confidence, collapse = ", "), ") instead.")
  }
  net <- .hc_ui_get_data_object("dorothea", if (identical(org, "mouse")) "dorothea_mm" else "dorothea_hs")
  if (base::is.null(net)) {
    stop("Could not load DoRothEA regulons. Install the `dorothea` package ",
         "(or `OmnipathR` for CollecTRI).", call. = FALSE)
  }
  net <- net %>% base::as.data.frame(stringsAsFactors = FALSE)
  if ("confidence" %in% base::colnames(net) && base::length(tf_confidence) > 0) {
    net <- net[base::as.character(net$confidence) %in% tf_confidence, , drop = FALSE]
  }
  out <- .hc_ui_prepare_network(
    net = net,
    source_candidates = c("source", "tf"),
    target_candidates = c("target"),
    mor_candidates = c("mor", "weight")
  )
  if (base::nrow(out) == 0) {
    stop("DoRothEA network is empty after filtering for confidence ",
         base::paste(tf_confidence, collapse = ", "), ".", call. = FALSE)
  }
  base::attr(out, "database") <- "DoRothEA"
  out
}

# CollecTRI straight from the OmniPath web service, prepared like
# decoupleR::get_collectri(split_complexes = FALSE): AP-1 and NF-kB complexes
# become the sources "AP1" and "NFKB"; inhibitory interactions get mor = -1.
.hc_ui_fetch_collectri <- function(organism = "human") {
  taxon <- if (identical(organism, "mouse")) 10090 else 9606
  url <- base::paste0(
    "https://omnipathdb.org/interactions?datasets=collectri&genesymbols=yes",
    "&loops=yes&organisms=", taxon, "&format=tsv"
  )
  old <- base::options(timeout = base::max(120, base::getOption("timeout")))
  on.exit(base::options(old), add = TRUE)
  raw <- utils::read.delim(url, stringsAsFactors = FALSE, quote = "")
  if (base::nrow(raw) == 0 || !"source_genesymbol" %in% base::names(raw)) {
    stop("OmniPath returned no CollecTRI interactions.")
  }
  is_true <- function(x) base::toupper(base::as.character(x)) %in% c("TRUE", "1")
  src <- raw$source_genesymbol
  is_complex <- base::startsWith(base::as.character(raw$source), "COMPLEX:")
  src[is_complex & base::grepl("JUN|FOS", src)] <- "AP1"
  src[is_complex & base::grepl("REL|NFKB", src)] <- "NFKB"
  mor <- base::ifelse(is_true(raw$is_inhibition) & !is_true(raw$is_stimulation), -1, 1)
  out <- base::data.frame(source = src, target = raw$target_genesymbol, mor = mor, stringsAsFactors = FALSE)
  out[!base::duplicated(out[, c("source", "target")]), , drop = FALSE]
}

# PROGENy pathway-responsive genes: the `top` genes per pathway with the
# smallest p-value, signed by their weight. Taken offline from the `progeny`
# package when installed, otherwise from OmniPath via decoupleR.
.hc_ui_load_pathway_network <- function(organism, top = 100) {
  org <- if (identical(organism, "mouse")) "mouse" else "human"
  model <- .hc_ui_get_data_object("progeny", if (identical(org, "mouse")) "model_mouse_full" else "model_human_full")
  if (!base::is.null(model)) {
    model <- model %>% base::as.data.frame(stringsAsFactors = FALSE)
    model <- model[base::order(model$pathway, model$p.value), , drop = FALSE]
    model <- base::do.call(base::rbind, base::lapply(
      base::split(model, model$pathway),
      function(d) utils::head(d, top)
    ))
    net <- base::data.frame(source = model$pathway, target = model$gene, mor = model$weight,
                            stringsAsFactors = FALSE)
  } else {
    net <- tryCatch(decoupleR::get_progeny(organism = org, top = top), error = function(e) NULL)
  }
  if (base::is.null(net)) {
    stop("Could not load the PROGENy model. Install the `progeny` package.", call. = FALSE)
  }
  out <- .hc_ui_prepare_network(
    net = net,
    source_candidates = c("source", "pathway"),
    target_candidates = c("target", "gene"),
    mor_candidates = c("mor", "weight")
  )
  if (base::nrow(out) == 0) {
    stop("PROGENy network is empty after preprocessing.", call. = FALSE)
  }
  out
}

# Hypergeometric over-representation of each regulator's targets among the
# genes of each module, against `universe`. Regulators with fewer than
# `minsize` targets in the universe and modules with fewer than `minsize`
# genes are not tested. The p-values are adjusted over all module x
# regulator tests; rows without any overlap are dropped afterwards.
.hc_ui_module_ora <- function(module_genes, network, universe, minsize = 5, padj = "BH") {
  network <- network[network$target %in% universe, , drop = FALSE]
  targets <- base::lapply(base::split(network$target, network$source), base::unique)
  targets <- targets[base::lengths(targets) >= minsize]
  n_universe <- base::length(universe)
  rows <- base::lapply(base::names(module_genes), function(cl) {
    genes <- base::intersect(module_genes[[cl]], universe)
    n <- base::length(genes)
    if (n < minsize || base::length(targets) == 0) {
      return(NULL)
    }
    overlap <- base::lapply(targets, function(t) t[t %in% genes])
    k <- base::lengths(overlap)
    big_k <- base::lengths(targets)
    base::data.frame(
      cluster = cl,
      term = base::names(targets),
      n_genes = n,
      n_targets = big_k,
      n_overlap = k,
      fold_enrichment = (k / n) / (big_k / n_universe),
      pvalue = stats::phyper(k - 1, big_k, n_universe - big_k, n, lower.tail = FALSE),
      overlap_genes = base::vapply(overlap, base::paste, base::character(1), collapse = ","),
      stringsAsFactors = FALSE
    )
  })
  rows <- rows[!base::vapply(rows, base::is.null, base::logical(1))]
  if (base::length(rows) == 0) {
    return(base::data.frame())
  }
  out <- base::do.call(base::rbind, rows)
  out$qvalue <- stats::p.adjust(out$pvalue, method = padj)
  out <- out[out$n_overlap > 0, , drop = FALSE]
  base::rownames(out) <- NULL
  out
}

# Signed score per condition for each module x regulator row of `ora`:
# weighted mean of mor * value over the regulator's targets in the module.
# Positive = the targets move as expected when the regulator is more active.
# `consistency` is the share of those targets whose mor * value has the same
# sign as the score (1 = all targets agree, 0.5 = no agreement).
.hc_ui_module_target_scores <- function(ora, network, values) {
  network <- network[!base::duplicated(network[, c("source", "target")]), , drop = FALSE]
  mor_by_source <- base::lapply(
    base::split(network, network$source),
    function(d) stats::setNames(d$mor, d$target)
  )
  n_cond <- base::ncol(values)
  per_row <- base::lapply(base::seq_len(base::nrow(ora)), function(i) {
    genes <- base::strsplit(ora$overlap_genes[i], ",", fixed = TRUE)[[1]]
    genes <- genes[genes %in% base::rownames(values)]
    w <- mor_by_source[[ora$term[i]]][genes]
    x <- values[genes, , drop = FALSE]
    present <- !base::is.na(x)
    signed <- x * w
    den <- base::colSums(present * base::abs(w))
    sc <- base::colSums(signed, na.rm = TRUE) / den
    sc[den == 0] <- NA_real_
    agree <- base::sign(signed) == base::matrix(base::sign(sc), nrow = base::nrow(signed), ncol = n_cond, byrow = TRUE)
    cons <- base::colSums(agree & present, na.rm = TRUE) / base::colSums(present)
    cons[den == 0] <- NA_real_
    list(score = sc, consistency = cons)
  })
  to_mat <- function(field) {
    base::matrix(
      base::unlist(base::lapply(per_row, `[[`, field)),
      nrow = base::nrow(ora),
      byrow = TRUE,
      dimnames = list(NULL, base::colnames(values))
    )
  }
  list(score = to_mat("score"), consistency = to_mat("consistency"))
}

# How a regulator acts on its targets in the module, from the sign of the
# mode of regulation (CollecTRI/DoRothEA) or of the pathway response weight
# (PROGENy): "activating" when at least two thirds of the targets are
# positive, "repressing" when at most one third, otherwise "mixed".
.hc_ui_module_regulation <- function(ora, network) {
  network <- network[!base::duplicated(network[, c("source", "target")]), , drop = FALSE]
  mor_by_source <- base::lapply(
    base::split(network, network$source),
    function(d) stats::setNames(d$mor, d$target)
  )
  base::vapply(base::seq_len(base::nrow(ora)), function(i) {
    genes <- base::strsplit(ora$overlap_genes[i], ",", fixed = TRUE)[[1]]
    w <- mor_by_source[[ora$term[i]]][genes]
    w <- w[!base::is.na(w) & w != 0]
    if (base::length(w) == 0) {
      return("mixed")
    }
    share <- base::mean(w > 0)
    if (share >= 2 / 3) "activating" else if (share <= 1 / 3) "repressing" else "mixed"
  }, base::character(1))
}

# Correlation between a regulator's own profile and the mean profile of a
# module across conditions (NA with fewer than 3 conditions or when the
# regulator is not measured).
.hc_ui_regulator_cor <- function(regulator, cluster, values, module_means) {
  if (!regulator %in% base::rownames(values) ||
    base::is.null(module_means) ||
    !cluster %in% base::rownames(module_means) ||
    base::ncol(values) < 3) {
    return(NA_real_)
  }
  cols <- base::intersect(base::colnames(values), base::colnames(module_means))
  a <- values[regulator, cols]
  b <- module_means[cluster, cols]
  ok <- !base::is.na(a) & !base::is.na(b)
  if (base::sum(ok) < 3 || stats::sd(a[ok]) == 0 || stats::sd(b[ok]) == 0) {
    return(NA_real_)
  }
  stats::cor(a[ok], b[ok])
}

# Regulator activity per condition with decoupleR ULM over all genes in
# `values` (one signature per column). P-values are adjusted per condition
# over all regulators.
.hc_ui_condition_activity <- function(values, network, minsize = 5, padj = "BH") {
  empty <- base::data.frame(
    source = base::character(0), condition = base::character(0),
    activity = base::numeric(0), activity_pvalue = base::numeric(0),
    activity_qvalue = base::numeric(0), stringsAsFactors = FALSE
  )
  values <- values[stats::complete.cases(values), , drop = FALSE]
  network <- network[network$target %in% base::rownames(values), , drop = FALSE]
  network <- network[!base::duplicated(network[, c("source", "target")]), , drop = FALSE]
  if (base::nrow(values) == 0 || base::nrow(network) == 0) {
    return(empty)
  }
  raw <- .hc_ui_run_decouple(mat = values, network = network, method = "ulm", minsize = minsize)
  if (base::is.null(raw) || base::nrow(raw) == 0) {
    return(empty)
  }
  out <- base::data.frame(
    source = base::as.character(raw$source),
    condition = base::as.character(raw$condition),
    activity = .hc_as_numeric_safely(raw$score),
    activity_pvalue = .hc_as_numeric_safely(raw$p_value),
    stringsAsFactors = FALSE
  )
  out$activity_qvalue <- stats::ave(
    out$activity_pvalue,
    out$condition,
    FUN = function(p) stats::p.adjust(p, method = padj)
  )
  out
}

# Within each module, mark significant regulators whose overlapping genes
# mostly repeat a better-ranked regulator (Jaccard index >= `jaccard`).
# TFs that are themselves genes of the module are preferred as the
# representative, then the rank. Returns the name of that regulator, or ""
# for representatives. `ora` must be sorted by module and rank.
.hc_ui_flag_redundant <- function(ora, keep, jaccard = 0.5) {
  out <- base::rep("", base::nrow(ora))
  genes <- base::strsplit(ora$overlap_genes, ",", fixed = TRUE)
  in_module <- if ("regulator_in_module" %in% base::names(ora)) {
    ora$regulator_in_module %in% TRUE
  } else {
    base::rep(FALSE, base::nrow(ora))
  }
  for (cl in base::unique(ora$cluster)) {
    idx <- base::which(ora$cluster == cl & keep)
    idx <- idx[base::order(!in_module[idx], base::seq_along(idx))]
    reps <- base::integer(0)
    for (i in idx) {
      hit <- ""
      for (r in reps) {
        j <- base::length(base::intersect(genes[[i]], genes[[r]])) /
          base::length(base::union(genes[[i]], genes[[r]]))
        if (j >= jaccard) {
          hit <- ora$term[r]
          break
        }
      }
      if (nzchar(hit)) {
        out[i] <- hit
      } else {
        reps <- base::c(reps, i)
      }
    }
  }
  out
}

.hc_ui_get_data_object <- function(pkg, obj_name) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    return(NULL)
  }
  ns <- asNamespace(pkg)
  if (exists(obj_name, envir = ns, inherits = FALSE)) {
    return(get(obj_name, envir = ns, inherits = FALSE))
  }
  tmp_env <- new.env(parent = emptyenv())
  loaded <- tryCatch(
    {
      utils::data(list = obj_name, package = pkg, envir = tmp_env)
      TRUE
    },
    error = function(e) FALSE
  )
  if (isTRUE(loaded) && exists(obj_name, envir = tmp_env, inherits = FALSE)) {
    return(get(obj_name, envir = tmp_env, inherits = FALSE))
  }
  NULL
}

.hc_ui_prepare_network <- function(net,
                                   source_candidates,
                                   target_candidates,
                                   mor_candidates = NULL) {
  net <- net %>% base::as.data.frame(stringsAsFactors = FALSE)
  find_col <- function(candidates) {
    hits <- candidates[candidates %in% base::colnames(net)]
    if (base::length(hits) == 0) {
      return(NA_character_)
    }
    hits[[1]]
  }
  source_col <- find_col(source_candidates)
  target_col <- find_col(target_candidates)
  if (base::is.na(source_col) || base::is.na(target_col)) {
    return(base::data.frame(source = base::character(0), target = base::character(0), mor = base::numeric(0)))
  }
  mor_col <- if (!base::is.null(mor_candidates)) find_col(mor_candidates) else NA_character_
  out <- base::data.frame(
    source = base::as.character(net[[source_col]]),
    target = base::as.character(net[[target_col]]),
    stringsAsFactors = FALSE
  )
  if (!base::is.na(mor_col)) {
    out$mor <- .hc_as_numeric_safely(net[[mor_col]])
  } else {
    out$mor <- 1
  }
  out$mor[base::is.na(out$mor)] <- 1
  out <- out[stats::complete.cases(out[, c("source", "target")]), , drop = FALSE]
  out <- out[out$source != "" & out$target != "", , drop = FALSE]
  out <- base::unique(out)
  rownames(out) <- NULL
  out
}

.hc_ui_run_decouple <- function(mat, network, method = "ulm", minsize = 5) {
  run_fn_name <- base::paste0("run_", method)
  run_fn <- get0(run_fn_name, envir = asNamespace("decoupleR"), mode = "function")
  if (base::is.null(run_fn)) {
    stop("`decoupleR::", run_fn_name, "` is not available in this decoupleR version.")
  }
  fn_args <- base::names(base::formals(run_fn))
  args <- list()
  if ("mat" %in% fn_args) {
    args$mat <- mat
  } else if (".mat" %in% fn_args) {
    args$.mat <- mat
  } else {
    args[[fn_args[[1]]]] <- mat
  }
  if ("network" %in% fn_args) {
    args$network <- network
  } else if (".network" %in% fn_args) {
    args$.network <- network
  }
  if (".source" %in% fn_args) {
    args$.source <- "source"
  }
  if (".target" %in% fn_args) {
    args$.target <- "target"
  }
  if (".mor" %in% fn_args) {
    args$.mor <- "mor"
  }
  if ("source" %in% fn_args && !(".source" %in% fn_args)) {
    args$source <- "source"
  }
  if ("target" %in% fn_args && !(".target" %in% fn_args)) {
    args$target <- "target"
  }
  if ("mor" %in% fn_args && !(".mor" %in% fn_args)) {
    args$mor <- "mor"
  }
  if ("minsize" %in% fn_args) {
    args$minsize <- minsize
  }
  out <- tryCatch(base::do.call(run_fn, args), error = function(e) {
    warning("decoupleR inference failed: ", conditionMessage(e))
    NULL
  })
  if (base::is.null(out)) {
    return(NULL)
  }
  out %>% base::as.data.frame(stringsAsFactors = FALSE)
}

.hc_ui_summarize_decouple_result <- function(df, padj = "BH") {
  if (base::is.null(df) || base::nrow(df) == 0) {
    return(base::data.frame())
  }
  col_pick <- function(cands) {
    hits <- cands[cands %in% base::colnames(df)]
    if (base::length(hits) == 0) {
      return(NA_character_)
    }
    hits[[1]]
  }
  source_col <- col_pick(c("source", "tf", "pathway", "regulator"))
  score_col <- col_pick(c("score", "statistic", "ulm", "norm_ulm", "estimate", "activity"))
  condition_col <- col_pick(c("condition", "sample", "group", "contrast"))
  pvalue_col <- col_pick(c("p_value", "pvalue", "p.val", "p_val", "p"))
  if (base::is.na(source_col) || base::is.na(score_col)) {
    return(base::data.frame())
  }
  if (base::is.na(condition_col)) {
    df$condition_hc <- "all"
    condition_col <- "condition_hc"
  }
  source_vals <- base::as.character(df[[source_col]])
  split_idx <- base::split(base::seq_len(base::nrow(df)), source_vals)

  out <- base::lapply(base::names(split_idx), function(src) {
    idx <- split_idx[[src]]
    scores <- .hc_as_numeric_safely(df[[score_col]][idx])
    scores <- scores[!base::is.na(scores)]
    if (base::length(scores) == 0) {
      return(NULL)
    }
    if (!base::is.na(pvalue_col)) {
      pv <- .hc_as_numeric_safely(df[[pvalue_col]][idx])
      pv <- pv[!base::is.na(pv) & is.finite(pv)]
    } else {
      pv <- base::numeric(0)
    }
    p_min <- if (base::length(pv) == 0) {
      2 * stats::pnorm(base::abs(base::mean(scores)), lower.tail = FALSE)
    } else {
      base::min(pv)
    }
    cond_vals <- base::as.character(df[[condition_col]][idx])
    base::data.frame(
      term = src,
      score = base::mean(scores),
      abs_score = base::mean(base::abs(scores)),
      pvalue = p_min,
      n_conditions = base::length(base::unique(cond_vals)),
      stringsAsFactors = FALSE
    )
  })
  out <- out[!base::vapply(out, base::is.null, FUN.VALUE = base::logical(1))]
  if (base::length(out) == 0) {
    return(base::data.frame())
  }
  out <- base::do.call(base::rbind, out)
  out$qvalue <- stats::p.adjust(out$pvalue, method = padj)
  out$direction <- base::ifelse(out$score >= 0, "activated", "inhibited")
  out <- out[base::order(out$qvalue, -out$abs_score, out$term), , drop = FALSE]
  rownames(out) <- NULL
  out
}

.hc_ui_summarize_decouple_by_condition <- function(df, padj = "BH") {
  if (base::is.null(df) || base::nrow(df) == 0) {
    return(base::data.frame())
  }
  col_pick <- function(cands) {
    hits <- cands[cands %in% base::colnames(df)]
    if (base::length(hits) == 0) {
      return(NA_character_)
    }
    hits[[1]]
  }
  source_col <- col_pick(c("source", "tf", "pathway", "regulator"))
  score_col <- col_pick(c("score", "statistic", "ulm", "norm_ulm", "estimate", "activity"))
  condition_col <- col_pick(c("condition", "sample", "group", "contrast"))
  pvalue_col <- col_pick(c("p_value", "pvalue", "p.val", "p_val", "p"))
  if (base::is.na(source_col) || base::is.na(score_col)) {
    return(base::data.frame())
  }
  if (base::is.na(condition_col)) {
    df$condition_hc <- "all"
    condition_col <- "condition_hc"
  }

  cond_split <- base::split(base::seq_len(base::nrow(df)), base::as.character(df[[condition_col]]))
  by_condition <- base::lapply(base::names(cond_split), function(cond_nm) {
    idx_cond <- cond_split[[cond_nm]]
    sub <- df[idx_cond, , drop = FALSE]
    src_vals <- base::as.character(sub[[source_col]])
    src_split <- base::split(base::seq_len(base::nrow(sub)), src_vals)
    out_src <- base::lapply(base::names(src_split), function(src) {
      idx_src <- src_split[[src]]
      scores <- .hc_as_numeric_safely(sub[[score_col]][idx_src])
      scores <- scores[!base::is.na(scores)]
      if (base::length(scores) == 0) {
        return(NULL)
      }
      if (!base::is.na(pvalue_col)) {
        pv <- .hc_as_numeric_safely(sub[[pvalue_col]][idx_src])
        pv <- pv[!base::is.na(pv) & is.finite(pv)]
      } else {
        pv <- base::numeric(0)
      }
      p_min <- if (base::length(pv) == 0) {
        2 * stats::pnorm(base::abs(base::mean(scores)), lower.tail = FALSE)
      } else {
        base::min(pv)
      }
      base::data.frame(
        condition = cond_nm,
        term = src,
        score = base::mean(scores),
        abs_score = base::mean(base::abs(scores)),
        pvalue = p_min,
        n_conditions = 1L,
        stringsAsFactors = FALSE
      )
    })
    out_src <- out_src[!base::vapply(out_src, base::is.null, FUN.VALUE = base::logical(1))]
    if (base::length(out_src) == 0) {
      return(NULL)
    }
    cond_df <- base::do.call(base::rbind, out_src)
    cond_df$qvalue <- stats::p.adjust(cond_df$pvalue, method = padj)
    cond_df$direction <- base::ifelse(cond_df$score >= 0, "activated", "inhibited")
    cond_df <- cond_df[base::order(cond_df$qvalue, -cond_df$abs_score, cond_df$term), , drop = FALSE]
    cond_df
  })
  by_condition <- by_condition[!base::vapply(by_condition, base::is.null, FUN.VALUE = base::logical(1))]
  if (base::length(by_condition) == 0) {
    return(base::data.frame())
  }
  out <- base::do.call(base::rbind, by_condition)
  base::rownames(out) <- NULL
  out
}

.hc_ui_upstream_term_order <- function(selected_all, cluster_levels) {
  if (base::is.null(selected_all) || base::nrow(selected_all) == 0) {
    return(base::character(0))
  }
  needed_cols <- c("term_with_resource", "resource", "cluster", "qvalue", "rank")
  if (!all(needed_cols %in% base::colnames(selected_all))) {
    return(base::unique(base::as.character(selected_all$term_with_resource)))
  }
  df <- selected_all
  df$term_with_resource <- base::as.character(df$term_with_resource)
  df$resource <- base::as.character(df$resource)
  df$cluster <- base::as.character(df$cluster)
  df$qvalue <- .hc_as_numeric_safely(df$qvalue)
  df$rank <- .hc_as_numeric_safely(df$rank)
  cluster_idx_map <- stats::setNames(base::seq_along(cluster_levels), cluster_levels)
  df$cluster_idx <- cluster_idx_map[df$cluster]
  df$cluster_idx[base::is.na(df$cluster_idx)] <- base::length(cluster_levels) + 1

  idx_by_term <- base::split(base::seq_len(base::nrow(df)), df$term_with_resource)
  term_meta <- base::do.call(base::rbind, base::lapply(base::names(idx_by_term), function(term_nm) {
    idx <- idx_by_term[[term_nm]]
    sub <- df[idx, , drop = FALSE]
    qv <- sub$qvalue[!base::is.na(sub$qvalue)]
    rk <- sub$rank[!base::is.na(sub$rank)]
    base::data.frame(
      term_with_resource = term_nm,
      resource = sub$resource[1],
      first_cluster_idx = base::min(sub$cluster_idx, na.rm = TRUE),
      mean_cluster_idx = base::mean(sub$cluster_idx, na.rm = TRUE),
      hit_count = base::length(base::unique(sub$cluster)),
      best_q = if (base::length(qv) == 0) Inf else base::min(qv),
      best_rank = if (base::length(rk) == 0) Inf else base::min(rk),
      stringsAsFactors = FALSE
    )
  }))
  if (base::is.null(term_meta) || base::nrow(term_meta) == 0) {
    return(base::character(0))
  }
  term_meta$first_cluster_idx[!is.finite(term_meta$first_cluster_idx)] <- base::length(cluster_levels) + 1
  term_meta$mean_cluster_idx[!is.finite(term_meta$mean_cluster_idx)] <- base::length(cluster_levels) + 1
  term_meta$best_q[!is.finite(term_meta$best_q)] <- Inf
  term_meta$best_rank[!is.finite(term_meta$best_rank)] <- Inf

  group_order <- base::unique(df$resource)
  cluster_buckets <- base::sort(base::unique(term_meta$first_cluster_idx))
  ordered_terms <- base::character(0)
  for (bucket_idx in cluster_buckets) {
    bucket_meta <- term_meta[term_meta$first_cluster_idx == bucket_idx, , drop = FALSE]
    if (base::nrow(bucket_meta) == 0) {
      next
    }
    group_lists <- base::lapply(group_order, function(gr) {
      x <- bucket_meta[bucket_meta$resource == gr, , drop = FALSE]
      if (base::nrow(x) == 0) {
        return(base::character(0))
      }
      x <- x[base::order(x$best_rank, x$best_q, -x$hit_count, x$mean_cluster_idx, x$term_with_resource), , drop = FALSE]
      x$term_with_resource
    })
    base::names(group_lists) <- group_order
    max_len <- base::max(base::lengths(group_lists))
    for (k in base::seq_len(max_len)) {
      for (gr in group_order) {
        cur <- group_lists[[gr]]
        if (base::length(cur) >= k) {
          ordered_terms <- base::c(ordered_terms, cur[k])
        }
      }
    }
  }
  ordered_terms <- base::unique(ordered_terms)
  missing_terms <- base::setdiff(term_meta$term_with_resource, ordered_terms)
  if (base::length(missing_terms) > 0) {
    missing_meta <- term_meta[base::match(missing_terms, term_meta$term_with_resource), , drop = FALSE]
    missing_meta <- missing_meta[base::order(missing_meta$best_rank, missing_meta$best_q, -missing_meta$hit_count, missing_meta$mean_cluster_idx, missing_meta$term_with_resource), , drop = FALSE]
    ordered_terms <- base::c(ordered_terms, missing_meta$term_with_resource)
  }
  ordered_terms
}

.hc_ui_build_upstream_combined_heatmap <- function(selected_all,
                                                   significant_all,
                                                   cluster_order,
                                                   module_label_map,
                                                   cluster_info,
                                                   gfc_all,
                                                   stored_hm = NULL,
                                                   main_heatmap_col_order = NULL,
                                                   heatmap_cluster_columns = FALSE,
                                                   heatmap_col_order = NULL,
                                                   module_heatmap_mat = NULL,
                                                   module_heatmap_col_order = NULL,
                                                   module_heatmap_name = "GFC",
                                                   gfc_scale_limits = c(-2, 2),
                                                   gfc_scale_ticks = NULL,
                                                   heatmap_side = "left",
                                                   term_levels = NULL,
                                                   activity_score_limit = 2,
                                                   show_horizontal_lines = FALSE,
                                                   mark_significant = FALSE,
                                                   significant_for_marks = NULL,
                                                   value_from_significant = TRUE,
                                                   overall_plot_scale = 1) {
  if (base::is.null(selected_all) || base::nrow(selected_all) == 0) {
    return(NULL)
  }
  if (!base::is.numeric(activity_score_limit) ||
    base::length(activity_score_limit) != 1 ||
    base::is.na(activity_score_limit)) {
    activity_score_limit <- 2
  }
  activity_score_limit <- base::abs(base::as.numeric(activity_score_limit))
  if (!base::is.finite(activity_score_limit) || activity_score_limit <= 0) {
    activity_score_limit <- 2
  }
  activity_score_limit <- base::max(1, base::min(4, activity_score_limit))

  if (!base::is.numeric(gfc_scale_limits) || base::length(gfc_scale_limits) != 2 || any(!base::is.finite(gfc_scale_limits))) {
    gfc_scale_limits <- c(-2, 2)
  }
  gfc_scale_limits <- base::sort(base::as.numeric(gfc_scale_limits))
  if (gfc_scale_limits[1] == gfc_scale_limits[2]) {
    lim_abs <- base::abs(gfc_scale_limits[1])
    if (!base::is.finite(lim_abs) || lim_abs <= 0) {
      lim_abs <- 2
    }
    gfc_scale_limits <- c(-lim_abs, lim_abs)
  }
  if (base::is.null(gfc_scale_ticks) ||
    !base::is.list(gfc_scale_ticks) ||
    !all(c("breaks", "labels") %in% base::names(gfc_scale_ticks))) {
    br <- pretty(gfc_scale_limits, n = 5)
    br <- br[
      br >= gfc_scale_limits[1] - .Machine$double.eps^0.5 &
        br <= gfc_scale_limits[2] + .Machine$double.eps^0.5
    ]
    if (!any(base::abs(br) < .Machine$double.eps^0.5)) {
      br <- base::sort(base::unique(base::c(br, 0)))
    }
    if (base::length(br) < 3) {
      br <- base::seq(gfc_scale_limits[1], gfc_scale_limits[2], length.out = 5)
    }
    tick_labels <- base::formatC(br, format = "fg", digits = 3)
    tick_labels <- base::trimws(tick_labels)
    tick_label_width <- base::max(base::nchar(tick_labels), na.rm = TRUE)
    tick_labels <- base::format(tick_labels, width = tick_label_width, justify = "right")
    gfc_scale_ticks <- list(
      breaks = br,
      labels = tick_labels
    )
  }

  extract_module_heatmap <- function(stored_hm, gfc_all, cluster_info, cluster_order, module_label_map, module_heatmap_mat) {
    if (!base::is.null(module_heatmap_mat) && base::nrow(module_heatmap_mat) > 0 && base::ncol(module_heatmap_mat) > 0) {
      hm_mat <- module_heatmap_mat %>% base::as.matrix()
      mode(hm_mat) <- "numeric"
      return(hm_mat)
    }
    hm_mat <- tryCatch(stored_hm@ht_list[[1]]@matrix, error = function(e) NULL)
    if (!base::is.null(hm_mat) && base::nrow(hm_mat) > 0) {
      hm_mat <- hm_mat %>% base::as.matrix()
      inv_map <- stats::setNames(base::names(module_label_map), base::as.character(module_label_map))
      mapped_rows <- ifelse(base::rownames(hm_mat) %in% base::names(inv_map), inv_map[base::rownames(hm_mat)], base::rownames(hm_mat))
      base::rownames(hm_mat) <- mapped_rows
      return(hm_mat)
    }

    if (base::is.null(gfc_all) || base::nrow(gfc_all) == 0 || base::ncol(gfc_all) < 2) {
      return(NULL)
    }
    gene_col <- if ("Gene" %in% base::colnames(gfc_all)) "Gene" else base::colnames(gfc_all)[base::ncol(gfc_all)]
    # By index, not by name: `GFC_all_layers` repeats its condition columns
    # when the layers share group names, and setdiff() on the names would
    # collapse those duplicates - silently dropping every layer after the
    # first.
    value_idx <- base::setdiff(base::seq_len(base::ncol(gfc_all)),
                               base::which(base::colnames(gfc_all) %in% gene_col))
    value_cols <- base::colnames(gfc_all)[value_idx]
    out_list <- list()
    for (cl in cluster_order) {
      genes <- dplyr::filter(cluster_info, color == cl) %>%
        dplyr::pull(., "gene_n") %>%
        base::strsplit(split = ",") %>%
        base::unlist()
      genes <- base::unique(base::as.character(genes))
      if (base::length(genes) == 0) {
        next
      }
      tmp <- gfc_all[gfc_all[[gene_col]] %in% genes, value_idx, drop = FALSE]
      if (base::nrow(tmp) == 0) {
        next
      }
      vals <- tmp %>% base::as.matrix()
      mode(vals) <- "numeric"
      out_list[[cl]] <- base::colMeans(vals, na.rm = TRUE)
    }
    if (base::length(out_list) == 0) {
      return(NULL)
    }
    hm_mat <- base::do.call(base::rbind, out_list)
    hm_mat %>% base::as.matrix()
  }

  heatmap_mat <- extract_module_heatmap(
    stored_hm = stored_hm,
    gfc_all = gfc_all,
    cluster_info = cluster_info,
    cluster_order = cluster_order,
    module_label_map = module_label_map,
    module_heatmap_mat = module_heatmap_mat
  )
  if (base::is.null(heatmap_mat) || base::nrow(heatmap_mat) == 0) {
    return(NULL)
  }

  prepared_cols <- .hc_prepare_plot_heatmap_columns(
    mat = heatmap_mat,
    cluster_columns = heatmap_cluster_columns,
    plot_order = heatmap_col_order,
    main_order = main_heatmap_col_order,
    fallback_order = module_heatmap_col_order,
    context = "upstream inference heatmap"
  )
  heatmap_mat <- prepared_cols$mat
  hc_col_dend <- prepared_cols$col_dend

  keep_clusters <- cluster_order[cluster_order %in% base::rownames(heatmap_mat)]
  if (base::length(keep_clusters) == 0) {
    return(NULL)
  }
  heatmap_mat <- heatmap_mat[keep_clusters, , drop = FALSE]
  module_labels <- base::as.character(module_label_map[keep_clusters])
  module_labels[is.na(module_labels) | module_labels == ""] <- keep_clusters[is.na(module_labels) | module_labels == ""]
  base::rownames(heatmap_mat) <- module_labels
  row_name_by_cluster <- stats::setNames(base::rownames(heatmap_mat), keep_clusters)
  row_color_by_name <- stats::setNames(keep_clusters, base::rownames(heatmap_mat))

  cluster_counts_df <- dplyr::filter(cluster_info, color %in% keep_clusters) %>%
    dplyr::distinct(color, .keep_all = TRUE)
  if ("gene_no" %in% base::colnames(cluster_counts_df)) {
    gene_count_vals <- .hc_as_numeric_safely(cluster_counts_df$gene_no)
  } else if ("gene_n" %in% base::colnames(cluster_counts_df)) {
    gene_count_vals <- base::vapply(
      base::as.character(cluster_counts_df$gene_n),
      function(x) {
        if (base::is.na(x) || x == "") {
          return(0)
        }
        parts <- base::trimws(base::unlist(base::strsplit(x, ",", fixed = TRUE)))
        base::sum(parts != "")
      },
      FUN.VALUE = base::numeric(1)
    )
  } else {
    gene_count_vals <- base::rep(NA_real_, base::nrow(cluster_counts_df))
  }
  gene_count_map <- stats::setNames(gene_count_vals, cluster_counts_df$color)
  gene_counts <- gene_count_map[keep_clusters]
  gene_counts[base::is.na(gene_counts)] <- 0

  font_axis <- 10 * overall_plot_scale
  font_annotation <- 11 * overall_plot_scale
  font_module <- 8.2 * overall_plot_scale
  cluster_calc <- hcobject[["integrated_output"]][["cluster_calc"]]
  stored_module_label_fontsize <- cluster_calc[["module_label_fontsize"]]
  if (!base::is.numeric(stored_module_label_fontsize) ||
    base::length(stored_module_label_fontsize) != 1 ||
    base::is.na(stored_module_label_fontsize) ||
    stored_module_label_fontsize <= 0) {
    stored_module_label_fontsize <- NULL
  }
  stored_module_label_pt_size <- cluster_calc[["module_label_pt_size"]]
  if (!base::is.numeric(stored_module_label_pt_size) ||
    base::length(stored_module_label_pt_size) != 1 ||
    base::is.na(stored_module_label_pt_size) ||
    stored_module_label_pt_size <= 0) {
    stored_module_label_pt_size <- NULL
  }
  stored_module_box_width_cm <- cluster_calc[["module_box_width_cm"]]
  if (!base::is.numeric(stored_module_box_width_cm) ||
    base::length(stored_module_box_width_cm) != 1 ||
    base::is.na(stored_module_box_width_cm) ||
    stored_module_box_width_cm <= 0) {
    stored_module_box_width_cm <- NULL
  }
  stored_module_box_to_cell_ratio <- cluster_calc[["module_box_to_cell_ratio"]]
  if (!base::is.numeric(stored_module_box_to_cell_ratio) ||
    base::length(stored_module_box_to_cell_ratio) != 1 ||
    base::is.na(stored_module_box_to_cell_ratio) ||
    stored_module_box_to_cell_ratio <= 0) {
    stored_module_box_to_cell_ratio <- NULL
  }
  stored_heatmap_cell_size_mm <- cluster_calc[["heatmap_cell_size_mm"]]
  if (!base::is.numeric(stored_heatmap_cell_size_mm) ||
    base::length(stored_heatmap_cell_size_mm) != 1 ||
    base::is.na(stored_heatmap_cell_size_mm) ||
    stored_heatmap_cell_size_mm <= 0) {
    stored_heatmap_cell_size_mm <- NULL
  }
  if (base::is.null(stored_module_box_to_cell_ratio) &&
    !base::is.null(stored_module_box_width_cm) &&
    !base::is.null(stored_heatmap_cell_size_mm)) {
    stored_module_box_to_cell_ratio <- (stored_module_box_width_cm * 10) / stored_heatmap_cell_size_mm
    if (!base::is.finite(stored_module_box_to_cell_ratio) ||
      stored_module_box_to_cell_ratio <= 0) {
      stored_module_box_to_cell_ratio <- NULL
    }
  }
  stored_gene_count_fontsize <- cluster_calc[["gene_count_fontsize"]]
  if (!base::is.numeric(stored_gene_count_fontsize) ||
    base::length(stored_gene_count_fontsize) != 1 ||
    base::is.na(stored_gene_count_fontsize) ||
    stored_gene_count_fontsize <= 0) {
    stored_gene_count_fontsize <- NULL
  }
  stored_gene_count_renderer <- cluster_calc[["gene_count_renderer"]]
  if (base::is.null(stored_gene_count_renderer) ||
    !base::is.character(stored_gene_count_renderer) ||
    base::length(stored_gene_count_renderer) != 1 ||
    !(stored_gene_count_renderer %in% c("pch", "text"))) {
    stored_gene_count_renderer <- "pch"
  }
  stored_gene_count_pt_size <- cluster_calc[["gene_count_pt_size"]]
  if (!base::is.numeric(stored_gene_count_pt_size) ||
    base::length(stored_gene_count_pt_size) != 1 ||
    base::is.na(stored_gene_count_pt_size) ||
    stored_gene_count_pt_size <= 0) {
    stored_gene_count_pt_size <- NULL
  }
  stored_gfc_colors <- cluster_calc[["gfc_colors"]]
  if (!base::is.character(stored_gfc_colors) ||
    base::length(stored_gfc_colors) < 2 ||
    any(base::is.na(stored_gfc_colors)) ||
    any(stored_gfc_colors == "")) {
    stored_gfc_colors <- NULL
  } else {
    stored_gfc_colors <- base::as.character(stored_gfc_colors)
  }
  gfc_colors <- if (!base::is.null(stored_gfc_colors)) {
    stored_gfc_colors
  } else {
    .hc_default_gfc_colors()
  }

  n_hc_rows <- base::nrow(heatmap_mat)
  n_hc_cols <- base::ncol(heatmap_mat)
  if (!base::is.null(stored_heatmap_cell_size_mm)) {
    hc_cell_mm <- stored_heatmap_cell_size_mm * overall_plot_scale
  } else {
    hc_cell_mm_base <- if (n_hc_rows <= 10) {
      7
    } else if (n_hc_rows <= 18) {
      6
    } else if (n_hc_rows <= 30) {
      5
    } else {
      4.2
    }
    hc_cell_mm <- base::max(3.6, base::min(10.5, hc_cell_mm_base)) * overall_plot_scale
  }

  module_color_map <- stats::setNames(keep_clusters, keep_clusters)
  max_label_chars <- base::max(base::nchar(module_labels), na.rm = TRUE)
  label_fontsize <- if (!base::is.null(stored_module_label_fontsize)) {
    stored_module_label_fontsize
  } else {
    base::max(6.3, base::min(12.0, font_module * 1.15))
  }
  module_box_width_cm <- if (!base::is.null(stored_module_box_width_cm)) {
    stored_module_box_width_cm
  } else {
    base::max(0.72, base::min(3.2, 0.34 + (0.14 * max_label_chars) + (0.035 * label_fontsize)))
  }
  module_label_pt_size <- if (!base::is.null(stored_module_label_pt_size)) {
    stored_module_label_pt_size
  } else {
    base::max(0.22, base::min(0.90, 0.075 * label_fontsize))
  }
  module_label_fit <- .hc_module_label_fit_pt(
    module_label_pt_size = module_label_pt_size,
    module_box_width_cm = module_box_width_cm,
    module_labels_display = module_labels,
    n_heat_rows = n_hc_rows,
    cell_size_mm = hc_cell_mm,
    module_label_fontsize = label_fontsize,
    use_fontsize_request = FALSE,
    fontface = "bold"
  )
  module_label_pt_size_unit <- if (base::length(module_labels) > 0 &&
    base::length(module_label_fit$pt_size) == base::length(module_labels)) {
    grid::unit(module_label_fit$pt_size, "pt")
  } else {
    grid::unit(module_label_pt_size, "snpc")
  }
  module_label_fontsize_draw <- .hc_module_label_effective_fontsize(
    module_label_fit = module_label_fit,
    fallback_fontsize = label_fontsize,
    fallback_pt_size = module_label_pt_size
  )
  module_box_anno <- .hc_module_label_box_annotation(
    values = keep_clusters,
    colors = module_color_map,
    labels = module_labels,
    label_color = "white",
    label_fontsize_pt = module_label_fit$base_pt_size,
    fontface = "bold",
    width_cm = module_box_width_cm,
    border_gp = grid::gpar(col = "black"),
    which = "row"
  )
  right_anno <- ComplexHeatmap::HeatmapAnnotation(
    modules = module_box_anno,
    which = "row",
    show_legend = FALSE,
    show_annotation_name = FALSE,
    gap = grid::unit(1.6 * overall_plot_scale, "mm")
  )

  column_gap_spec <- .hc_heatmap_column_gap_spec(
    hcobject = hcobject,
    cols = base::colnames(heatmap_mat),
    cluster_columns = heatmap_cluster_columns,
    gap_mm = 0.6 * overall_plot_scale
  )
  hc_body_w_mm <- if (n_hc_cols <= 1) {
    hc_cell_mm
  } else {
    base::max(18, (n_hc_cols * hc_cell_mm) + column_gap_spec$total_gap_mm)
  }
  hc_body_h_mm <- base::max(20, n_hc_rows * hc_cell_mm)
  gfc_palette <- grDevices::colorRampPalette(gfc_colors)(51)
  gfc_col_fun <- circlize::colorRamp2(
    seq(gfc_scale_limits[1], gfc_scale_limits[2], length.out = base::length(gfc_palette)),
    gfc_palette
  )
  heatmap_column_labels_display <- .hc_gfc_display_col_labels(hcobject, base::colnames(heatmap_mat))
  hc_ht_args <- list(
    matrix = heatmap_mat,
    name = module_heatmap_name,
    right_annotation = right_anno,
    col = gfc_col_fun,
    cluster_rows = FALSE,
    cluster_columns = if (base::is.null(hc_col_dend)) FALSE else hc_col_dend,
    show_row_names = FALSE,
    show_column_dend = !base::is.null(hc_col_dend),
    show_row_dend = FALSE,
    column_names_rot = 90,
    column_labels = heatmap_column_labels_display,
    column_names_gp = grid::gpar(fontsize = font_axis),
    width = grid::unit(hc_body_w_mm, "mm"),
    height = grid::unit(hc_body_h_mm, "mm"),
    rect_gp = grid::gpar(col = "black"),
    heatmap_legend_param = list(
      title = module_heatmap_name,
      at = gfc_scale_ticks$breaks,
      labels = gfc_scale_ticks$labels,
      title_gp = grid::gpar(fontsize = 7.6 * overall_plot_scale, fontface = "bold"),
      labels_gp = grid::gpar(fontsize = 6.6 * overall_plot_scale)
    )
  )
  hc_ht_args <- .hc_heatmap_add_column_gap_args(
    hc_ht_args,
    column_gap_spec = column_gap_spec,
    title_gp = grid::gpar(fontsize = font_axis, fontface = "bold")
  )
  hc_ht <- do.call(ComplexHeatmap::Heatmap, hc_ht_args)

  df <- selected_all
  df <- df[df$cluster %in% keep_clusters, , drop = FALSE]
  if (base::nrow(df) == 0) {
    return(hc_ht)
  }
  if (base::is.null(term_levels) || base::length(term_levels) == 0) {
    term_levels <- .hc_ui_upstream_term_order(
      selected_all = df,
      cluster_levels = keep_clusters
    )
  } else {
    term_levels <- base::as.character(term_levels)
    term_levels <- term_levels[term_levels %in% base::unique(base::as.character(df$term_with_resource))]
  }
  if (base::length(term_levels) == 0) {
    return(hc_ht)
  }

  upstream_mat <- base::matrix(
    NA_real_,
    nrow = base::nrow(heatmap_mat),
    ncol = base::length(term_levels),
    dimnames = list(base::rownames(heatmap_mat), term_levels)
  )
  q_mat <- base::matrix(
    NA_real_,
    nrow = base::nrow(heatmap_mat),
    ncol = base::length(term_levels),
    dimnames = list(base::rownames(heatmap_mat), term_levels)
  )
  sig_mat <- base::matrix(
    FALSE,
    nrow = base::nrow(heatmap_mat),
    ncol = base::length(term_levels),
    dimnames = list(base::rownames(heatmap_mat), term_levels)
  )
  # TF that is itself a gene of the module (co-expressed with its targets).
  member_mat <- sig_mat
  src <- if (isTRUE(value_from_significant) && !base::is.null(significant_all) && base::nrow(significant_all) > 0) {
    significant_all
  } else {
    selected_all
  }
  src <- src[src$cluster %in% keep_clusters & src$term_with_resource %in% term_levels, , drop = FALSE]
  if (base::nrow(src) > 0) {
    src <- src[base::order(src$qvalue, -src$abs_score), , drop = FALSE]
    for (k in base::seq_len(base::nrow(src))) {
      rn <- row_name_by_cluster[[base::as.character(src$cluster[k])]]
      cn <- base::as.character(src$term_with_resource[k])
      if (!base::is.null(rn) && rn %in% base::rownames(upstream_mat) && cn %in% base::colnames(upstream_mat)) {
        old_q <- q_mat[rn, cn]
        new_q <- .hc_first_numeric_value(src$qvalue[k])
        if (base::is.na(old_q) || (!base::is.na(new_q) && new_q < old_q)) {
          upstream_mat[rn, cn] <- .hc_first_numeric_value(src$score[k])
          q_mat[rn, cn] <- new_q
          if ("regulator_in_module" %in% base::colnames(src)) {
            member_mat[rn, cn] <- isTRUE(base::as.logical(src$regulator_in_module[k]))
          }
        }
      }
    }
  }
  sig_src <- if (!base::is.null(significant_for_marks)) {
    significant_for_marks
  } else if (!base::is.null(significant_all) && base::nrow(significant_all) > 0) {
    significant_all
  } else {
    base::data.frame()
  }
  if (base::nrow(sig_src) > 0) {
    sig_src <- sig_src[sig_src$cluster %in% keep_clusters & sig_src$term_with_resource %in% term_levels, , drop = FALSE]
    if (base::nrow(sig_src) > 0) {
      for (k in base::seq_len(base::nrow(sig_src))) {
        rn <- row_name_by_cluster[[base::as.character(sig_src$cluster[k])]]
        cn <- base::as.character(sig_src$term_with_resource[k])
        if (!base::is.null(rn) && rn %in% base::rownames(sig_mat) && cn %in% base::colnames(sig_mat)) {
          sig_mat[rn, cn] <- TRUE
        }
      }
    }
  }

  resource_by_term <- selected_all[base::match(term_levels, selected_all$term_with_resource), c("term_with_resource", "resource"), drop = FALSE]
  resource_by_term <- stats::setNames(base::as.character(resource_by_term$resource), base::as.character(resource_by_term$term_with_resource))
  resource_colors <- c(TF = "#3B7EA1", Pathway = "#B26B2C")
  top_anno <- ComplexHeatmap::HeatmapAnnotation(
    resource = resource_by_term[term_levels],
    col = list(resource = resource_colors),
    show_annotation_name = FALSE,
    annotation_legend_param = list(title = "Resource")
  )

  score_lim <- activity_score_limit
  activity_col_fun <- circlize::colorRamp2(c(-score_lim, 0, score_lim), c("#2c7bb6", "#f7f7f7", "#d7191c"))
  n_cols_up <- base::ncol(upstream_mat)
  n_rows_up <- base::nrow(upstream_mat)
  q_to_pt_size <- function(q) {
    if (base::is.na(q) || q <= 0) {
      return(grid::unit(0.6, "mm"))
    }
    sig <- -base::log10(base::pmax(q, 1e-300))
    val <- base::min(1, base::max(0, (sig - 1) / 6))
    base_max_mm <- if (n_cols_up <= 10) {
      3.2
    } else if (n_cols_up <= 20) {
      2.6
    } else if (n_cols_up <= 35) {
      2.2
    } else {
      1.8
    }
    if (n_rows_up <= 10) {
      base_max_mm <- base_max_mm + 0.2
    }
    grid::unit(0.6 + (base_max_mm - 0.6) * val, "mm")
  }
  upstream_cell_w_mm_base <- if (n_cols_up <= 8) {
    8
  } else if (n_cols_up <= 20) {
    6
  } else {
    4.8
  }
  upstream_body_w_mm <- base::max(24, n_cols_up * upstream_cell_w_mm_base * overall_plot_scale)
  upstream_marks_title <- base::paste(c(
    if (isTRUE(mark_significant) && base::any(sig_mat)) "* active in this condition",
    if (base::any(member_mat)) "framed: TF is itself a gene of the module"
  ), collapse = "    ")
  if (!base::nzchar(upstream_marks_title)) {
    upstream_marks_title <- NULL
  }
  upstream_ht <- ComplexHeatmap::Heatmap(
    upstream_mat,
    name = "Target score",
    col = activity_col_fun,
    na_col = "grey96",
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    column_names_rot = 90,
    column_names_gp = grid::gpar(fontsize = font_axis),
    width = grid::unit(upstream_body_w_mm, "mm"),
    top_annotation = top_anno,
    rect_gp = grid::gpar(col = "grey85"),
    column_title = upstream_marks_title,
    column_title_side = "bottom",
    column_title_gp = grid::gpar(fontsize = font_axis * 0.9, col = "grey30"),
    cell_fun = function(j, i, x, y, width, height, fill) {
      val <- upstream_mat[i, j]
      if (!base::is.na(val)) {
        row_nm <- base::rownames(upstream_mat)[i]
        module_col <- row_color_by_name[[row_nm]]
        if (base::is.null(module_col) || base::is.na(module_col) || module_col == "") {
          module_col <- "grey70"
        }
        if (isTRUE(show_horizontal_lines)) {
          grid::grid.segments(
            x0 = grid::unit(0, "npc"),
            y0 = y,
            x1 = x,
            y1 = y,
            gp = grid::gpar(
              col = grDevices::adjustcolor(module_col, alpha.f = 0.55),
              lwd = 0.8 * overall_plot_scale
            )
          )
        }
        pt_col <- grDevices::adjustcolor(activity_col_fun(val), alpha.f = 0.85)
        grid::grid.points(
          x = x,
          y = y,
          pch = 16,
          size = q_to_pt_size(q_mat[i, j]),
          gp = grid::gpar(col = pt_col, fill = pt_col)
        )
        if (isTRUE(member_mat[i, j])) {
          grid::grid.rect(
            x = x,
            y = y,
            width = width,
            height = height,
            gp = grid::gpar(col = "black", fill = NA, lwd = 1.6 * overall_plot_scale)
          )
        }
        if (isTRUE(mark_significant) && isTRUE(sig_mat[i, j])) {
          grid::grid.text(
            label = "*",
            x = x,
            y = y,
            gp = grid::gpar(col = "black", fontsize = base::max(7, 8 * overall_plot_scale), fontface = "bold")
          )
        }
      }
    }
  )

  if (identical(heatmap_side, "left")) {
    hc_ht + upstream_ht
  } else {
    upstream_ht + hc_ht
  }
}

.hc_ui_build_upstream_dotplot <- function(selected_all,
                                          cluster_order,
                                          module_label_map,
                                          term_levels = NULL,
                                          overall_plot_scale = 1) {
  if (base::is.null(selected_all) || base::nrow(selected_all) == 0) {
    warning("No significant upstream hits found for selected resources/modules.")
    return(NULL)
  }
  df <- selected_all
  if (base::is.null(term_levels) || base::length(term_levels) == 0) {
    term_levels <- .hc_ui_upstream_term_order(
      selected_all = df,
      cluster_levels = cluster_order
    )
  } else {
    term_levels <- base::as.character(term_levels)
  }

  module_levels <- base::rev(module_label_map[cluster_order])
  df$module_label <- base::factor(df$module_label, levels = module_levels)
  df$term_with_resource <- base::factor(df$term_with_resource, levels = term_levels)
  df$neglog10_q <- -base::log10(base::pmax(df$qvalue, 1e-300))
  df$neglog10_q[df$neglog10_q > 20] <- 20

  ggplot2::ggplot(df, ggplot2::aes(x = term_with_resource, y = module_label)) +
    ggplot2::geom_point(
      ggplot2::aes(size = neglog10_q, color = score, shape = resource),
      alpha = 0.9
    ) +
    ggplot2::scale_color_gradient2(
      low = "#2c7bb6",
      mid = "#f7f7f7",
      high = "#d7191c",
      midpoint = 0,
      name = "Target score"
    ) +
    ggplot2::scale_size_continuous(name = "-log10(FDR)") +
    ggplot2::scale_shape_manual(values = c(TF = 16, Pathway = 17)) +
    ggplot2::theme_bw(base_size = 10 * overall_plot_scale) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
      panel.grid.major = ggplot2::element_line(color = "grey90"),
      panel.grid.minor = ggplot2::element_blank()
    ) +
    ggplot2::ggtitle("Upstream inference (TF + Pathway, mixed)")
}

.hc_ui_build_upstream_heatmap <- function(significant_all,
                                          selected_all,
                                          cluster_order,
                                          module_label_map,
                                          term_levels = NULL,
                                          activity_score_limit = 2,
                                          overall_plot_scale = 1) {
  if (base::is.null(selected_all) || base::nrow(selected_all) == 0) {
    return(NULL)
  }
  if (!base::is.numeric(activity_score_limit) ||
    base::length(activity_score_limit) != 1 ||
    base::is.na(activity_score_limit)) {
    activity_score_limit <- 2
  }
  activity_score_limit <- base::abs(base::as.numeric(activity_score_limit))
  if (!base::is.finite(activity_score_limit) || activity_score_limit <= 0) {
    activity_score_limit <- 2
  }
  activity_score_limit <- base::max(1, base::min(4, activity_score_limit))
  if (base::is.null(term_levels) || base::length(term_levels) == 0) {
    term_levels <- base::unique(base::as.character(selected_all$term_with_resource))
  } else {
    term_levels <- base::as.character(term_levels)
    available_terms <- base::unique(base::as.character(selected_all$term_with_resource))
    term_levels <- term_levels[term_levels %in% available_terms]
  }
  if (base::length(term_levels) == 0) {
    return(NULL)
  }
  module_levels <- module_label_map[cluster_order]
  hm <- base::matrix(
    NA_real_,
    nrow = base::length(cluster_order),
    ncol = base::length(term_levels),
    dimnames = list(module_levels, term_levels)
  )

  src <- if (!base::is.null(significant_all) && base::nrow(significant_all) > 0) {
    significant_all
  } else {
    selected_all
  }
  src <- src[src$term_with_resource %in% term_levels, , drop = FALSE]
  if (base::nrow(src) > 0) {
    src <- src[base::order(src$qvalue, -src$abs_score), , drop = FALSE]
    key <- base::paste(src$module_label, src$term_with_resource, sep = "\t")
    src <- src[!duplicated(key), , drop = FALSE]
    for (i in base::seq_len(base::nrow(src))) {
      if (src$module_label[i] %in% base::rownames(hm) &&
        src$term_with_resource[i] %in% base::colnames(hm)) {
        hm[src$module_label[i], src$term_with_resource[i]] <- src$score[i]
      }
    }
  }

  resource_by_term <- selected_all[base::match(term_levels, selected_all$term_with_resource), c("term_with_resource", "resource"), drop = FALSE]
  resource_by_term <- stats::setNames(base::as.character(resource_by_term$resource), base::as.character(resource_by_term$term_with_resource))
  resource_colors <- c(TF = "#3B7EA1", Pathway = "#B26B2C")
  top_anno <- ComplexHeatmap::HeatmapAnnotation(
    resource = resource_by_term[base::colnames(hm)],
    col = list(resource = resource_colors),
    show_annotation_name = FALSE,
    annotation_legend_param = list(title = "Resource")
  )
  col_fun <- circlize::colorRamp2(
    c(-activity_score_limit, 0, activity_score_limit),
    c("#2c7bb6", "#f7f7f7", "#d7191c")
  )
  ComplexHeatmap::Heatmap(
    hm,
    name = "Target score",
    col = col_fun,
    na_col = "grey95",
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    top_annotation = top_anno,
    row_names_gp = grid::gpar(fontsize = 10 * overall_plot_scale),
    column_names_gp = grid::gpar(fontsize = 9 * overall_plot_scale),
    column_names_rot = 90,
    rect_gp = grid::gpar(col = "grey80")
  )
}
