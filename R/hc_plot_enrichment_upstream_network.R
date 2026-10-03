#' Plot a module knowledge network from enrichment + upstream inference
#'
#' Builds a three-column network: modules -> enrichment terms -> upstream
#' regulators (TF/Pathway). A term is linked to a regulator when a relevant
#' share of the regulator's targets in the module belongs to the term
#' (`link_min_share`); regulators without such a term are linked to the module
#' directly. Arrows show how a regulator acts on its targets in the module
#' (activating ->, repressing -|, mixed without end). Regulators that are
#' themselves genes of a module are filled with that module's colour.
#' The output is shown together with the hCoCena module heatmap to preserve the
#' direct connection to module-level expression patterns. The heatmap panel
#' automatically follows the last upstream inference settings (e.g. GFC vs FC).
#'
#' @param enrichment_mode Character scalar. One of `"selected"` (default) or
#'   `"significant"`.
#' @param upstream_mode Character scalar. One of `"selected"` (default) or
#'   `"significant"`.
#' @param clusters Either `"all"` (default) or a character vector of module
#'   colors to include.
#' @param max_enrichment_per_module Optional positive integer. If set, keeps at
#'   most this many enrichment edges per module (best q-values first).
#' @param max_upstream_per_module Optional positive integer. If set, keeps at
#'   most this many upstream edges per module (best q-values first).
#' @param collapse_redundant_terms Logical. If `TRUE` (default), enrichment
#'   terms whose module genes lie mostly (>= 80 %) in a better term of the same
#'   module are left out.
#' @param link_min_share Minimum share of a regulator's targets in the module
#'   that must belong to a term to draw a term -> regulator line. Default 0.25.
#' @param label_mode Character scalar controlling term label density:
#'   `"both"` (default), `"upstream_only"`, or `"focus_only"`.
#' @param show_plot Logical; if `TRUE` (default), prints the combined overview
#'   (heatmap + network) in the active graphics device.
#' @param save_pdf Logical; if `TRUE` (default), writes a multi-page PDF to the
#'   current hCoCena save folder (overview + per-module focus pages).
#' @param pdf_name Output PDF filename.
#' @param gfc_scale_limits Optional numeric vector controlling the left module
#'   heatmap color scale limits. Provide one positive number (`x` -> `c(-x, x)`)
#'   or two numbers (`c(min, max)`). If NULL, uses upstream inference settings
#'   first (if available), then main heatmap settings, then `c(-range_GFC, range_GFC)`.
#' @param col_order Optional character vector overriding the hCoCena
#'   heatmap column order for this knowledge-network plot only. If `NULL`
#'   (default), the column order from the main module heatmap is reused when
#'   available.
#' @param heatmap_col_order Legacy alias for `col_order`.
#' @param cluster_columns Logical. If `FALSE` (default), reuse the
#'   column order from the main hCoCena heatmap when available. If `TRUE`,
#'   cluster the columns for this knowledge-network plot instead.
#' @param heatmap_cluster_columns Legacy alias for `cluster_columns`.
#' @param pdf_width Optional numeric width (inches) for network PDFs.
#'   If NULL (default), width is auto-estimated from content.
#' @param pdf_height Optional numeric height (inches) for network PDFs.
#'   If NULL (default), height is auto-estimated from content.
#' @param pdf_pointsize Numeric base pointsize used for PDF export devices.
#'   Default is 11.
#' @param overall_plot_scale Numeric scaling factor for text and line sizes.
#'
#' @return A list with overview/focus plots, nodes, edges, and output file.
#' @noRd
.hc_plot_enrichment_upstream_network_driver <- function(enrichment_mode = "selected",
                                             upstream_mode = "selected",
                                             clusters = c("all"),
                                             max_enrichment_per_module = NULL,
                                             max_upstream_per_module = NULL,
                                             collapse_redundant_terms = TRUE,
                                             link_min_share = 0.25,
                                             label_mode = "both",
                                             show_plot = TRUE,
                                             save_pdf = TRUE,
                                             pdf_name = "Module_Knowledge_Network.pdf",
                                             gfc_scale_limits = NULL,
                                             col_order = NULL,
                                             heatmap_col_order = NULL,
                                             cluster_columns = FALSE,
                                             heatmap_cluster_columns = NULL,
                                             pdf_width = NULL,
                                             pdf_height = NULL,
                                             pdf_pointsize = 11,
                                             overall_plot_scale = 1) {

  enrichment_mode <- base::match.arg(
    base::tolower(base::as.character(enrichment_mode)),
    choices = c("selected", "significant")
  )
  upstream_mode <- base::match.arg(
    base::tolower(base::as.character(upstream_mode)),
    choices = c("selected", "significant")
  )
  label_mode <- base::match.arg(
    base::tolower(base::as.character(label_mode)),
    choices = c("both", "upstream_only", "focus_only")
  )
  if (!base::is.logical(show_plot) || base::length(show_plot) != 1) {
    stop("`show_plot` must be TRUE or FALSE.")
  }
  if (!base::is.logical(save_pdf) || base::length(save_pdf) != 1) {
    stop("`save_pdf` must be TRUE or FALSE.")
  }
  if (!base::is.character(pdf_name) || base::length(pdf_name) != 1 || base::is.na(pdf_name) || pdf_name == "") {
    stop("`pdf_name` must be a non-empty filename.")
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
    context = ".hc_plot_enrichment_upstream_network_driver()"
  )
  cluster_columns <- .hc_resolve_cluster_columns_alias(
    cluster_columns = cluster_columns,
    heatmap_cluster_columns = heatmap_cluster_columns,
    cluster_columns_missing = missing(cluster_columns),
    heatmap_cluster_columns_missing = missing(heatmap_cluster_columns),
    context = ".hc_plot_enrichment_upstream_network_driver()"
  )
  if (!base::is.logical(cluster_columns) ||
    base::length(cluster_columns) != 1 ||
    base::is.na(cluster_columns)) {
    stop("`cluster_columns` must be TRUE or FALSE.")
  }

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

  check_optional_limit <- function(x, arg) {
    if (base::is.null(x)) {
      return(NULL)
    }
    if (!base::is.numeric(x) || base::length(x) != 1 || base::is.na(x) || x < 1) {
      stop("`", arg, "` must be NULL or a positive integer.")
    }
    base::as.integer(x)
  }
  max_enrichment_per_module <- check_optional_limit(max_enrichment_per_module, "max_enrichment_per_module")
  max_upstream_per_module <- check_optional_limit(max_upstream_per_module, "max_upstream_per_module")

  file_prefix <- base::paste0(
    hcobject[["working_directory"]][["dir_output"]],
    hcobject[["global_settings"]][["save_folder"]]
  )

  cluster_info <- hcobject[["integrated_output"]][["cluster_calc"]][["cluster_information"]]
  if (base::is.null(cluster_info) || base::nrow(cluster_info) == 0) {
    stop("No cluster information found. Run `hc_cluster_calculation()` first.")
  }
  all_clusters <- base::unique(base::as.character(cluster_info$color))
  all_clusters <- all_clusters[all_clusters != "white" & !base::is.na(all_clusters)]
  if (base::length(all_clusters) == 0) {
    stop("No non-white modules available for plotting.")
  }
  if (clusters[1] == "all") {
    clusters <- all_clusters
  }
  clusters <- base::as.character(clusters)
  clusters <- clusters[clusters %in% all_clusters]
  if (base::length(clusters) == 0) {
    stop("No valid clusters selected.")
  }

  cluster_calc <- hcobject[["integrated_output"]][["cluster_calc"]]
  heatmap_info <- .hc_heatmap_cache_info(cluster_calc)
  stored_hm <- heatmap_info$heatmap_obj
  main_heatmap_col_order <- heatmap_info$col_order

  cluster_order <- all_clusters
  if (!base::is.null(heatmap_info$row_order) && base::length(heatmap_info$row_order) > 0) {
    cluster_order <- heatmap_info$row_order
    cluster_order <- cluster_order[cluster_order %in% all_clusters]
    cluster_order <- base::c(cluster_order, base::setdiff(all_clusters, cluster_order))
  }
  cluster_order <- cluster_order[cluster_order %in% clusters]
  if (base::length(cluster_order) == 0) {
    stop("No modules remain after applying order and `clusters` filter.")
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
  module_prefix <- cluster_calc[["module_prefix"]]
  if (base::is.null(module_prefix) || !base::is.character(module_prefix) || base::length(module_prefix) != 1) {
    module_prefix <- "M"
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
  module_label_map <- module_label_map[cluster_order]

  read_xlsx_sheet <- function(file, sheet) {
    if (!requireNamespace("openxlsx", quietly = TRUE)) {
      return(NULL)
    }
    if (!base::file.exists(file)) {
      return(NULL)
    }
    out <- tryCatch(
      openxlsx::read.xlsx(xlsxFile = file, sheet = sheet),
      error = function(e) NULL
    )
    if (base::is.null(out) || !base::is.data.frame(out) || base::nrow(out) == 0) {
      return(NULL)
    }
    out
  }
  read_enrichment_xlsx_by_mode <- function(file, mode = c("selected", "significant")) {
    mode <- base::match.arg(mode)
    if (!requireNamespace("openxlsx", quietly = TRUE)) {
      return(NULL)
    }
    if (!base::file.exists(file)) {
      return(NULL)
    }
    sheet_names <- tryCatch(openxlsx::getSheetNames(file), error = function(e) NULL)
    if (base::is.null(sheet_names) || base::length(sheet_names) == 0) {
      return(NULL)
    }
    if (identical(mode, "selected")) {
      top_like <- base::grep("^top_[0-9]+_enrichments_all_dbs$", sheet_names, value = TRUE)
      if (base::length(top_like) > 1) {
        ord <- base::order(
          .hc_as_integer_safely(base::sub("^top_([0-9]+)_.*$", "\\1", top_like)),
          decreasing = TRUE,
          na.last = TRUE
        )
        top_like <- top_like[ord]
      }
      candidate_sheets <- base::unique(base::c("selected_enrichments_all_dbs", top_like))
    } else {
      candidate_sheets <- "significant_enrichments_all_dbs"
    }
    for (nm in candidate_sheets) {
      if (!(nm %in% sheet_names)) {
        next
      }
      tmp <- read_xlsx_sheet(file = file, sheet = nm)
      if (!base::is.null(tmp) && base::nrow(tmp) > 0) {
        return(tmp)
      }
    }
    NULL
  }

  normalize_enrichment_df <- function(df) {
    if (base::is.null(df) || !base::is.data.frame(df) || base::nrow(df) == 0) {
      return(NULL)
    }
    out <- base::as.data.frame(df, stringsAsFactors = FALSE)
    if (!("database" %in% base::colnames(out))) {
      if ("DB" %in% base::colnames(out)) {
        out$database <- base::as.character(out$DB)
      } else {
        out$database <- NA_character_
      }
    }
    if (!("cluster" %in% base::colnames(out))) {
      return(NULL)
    }
    if (!("term" %in% base::colnames(out))) {
      if ("Description" %in% base::colnames(out)) {
        out$term <- base::as.character(out$Description)
      } else {
        return(NULL)
      }
    }
    if (!("qvalue" %in% base::colnames(out))) {
      if ("p.adjust" %in% base::colnames(out)) {
        out$qvalue <- .hc_as_numeric_safely(out$p.adjust)
      } else if ("pvalue" %in% base::colnames(out)) {
        out$qvalue <- .hc_as_numeric_safely(out$pvalue)
      } else {
        out$qvalue <- NA_real_
      }
    }
    out$database <- base::as.character(out$database)
    out$cluster <- base::as.character(out$cluster)
    out$term <- base::as.character(out$term)
    out$qvalue <- .hc_as_numeric_safely(out$qvalue)
    out <- out[!base::is.na(out$cluster) & !base::is.na(out$term) & out$term != "", , drop = FALSE]
    if (base::nrow(out) == 0) {
      return(NULL)
    }
    out
  }

  normalize_upstream_df <- function(df) {
    if (base::is.null(df) || !base::is.data.frame(df) || base::nrow(df) == 0) {
      return(NULL)
    }
    out <- base::as.data.frame(df, stringsAsFactors = FALSE)
    if (!("resource" %in% base::colnames(out))) {
      if ("Resource" %in% base::colnames(out)) {
        out$resource <- base::as.character(out$Resource)
      } else {
        out$resource <- NA_character_
      }
    }
    if (!("cluster" %in% base::colnames(out))) {
      return(NULL)
    }
    if (!("term" %in% base::colnames(out))) {
      if ("source" %in% base::colnames(out)) {
        out$term <- base::as.character(out$source)
      } else {
        return(NULL)
      }
    }
    if (!("qvalue" %in% base::colnames(out))) {
      if ("p.adjust" %in% base::colnames(out)) {
        out$qvalue <- .hc_as_numeric_safely(out$p.adjust)
      } else if ("pvalue" %in% base::colnames(out)) {
        out$qvalue <- .hc_as_numeric_safely(out$pvalue)
      } else {
        out$qvalue <- NA_real_
      }
    }
    if (!("direction" %in% base::colnames(out))) {
      out$direction <- NA_character_
    }
    out$resource <- base::as.character(out$resource)
    out$resource[base::tolower(out$resource) == "pathway"] <- "Pathway"
    out$resource[base::tolower(out$resource) == "tf"] <- "TF"
    out$cluster <- base::as.character(out$cluster)
    out$term <- base::as.character(out$term)
    out$qvalue <- .hc_as_numeric_safely(out$qvalue)
    out$direction <- base::as.character(out$direction)
    out <- out[!base::is.na(out$cluster) & !base::is.na(out$term) & out$term != "", , drop = FALSE]
    if (base::nrow(out) == 0) {
      return(NULL)
    }
    out
  }

  resolve_enrichment_df <- function(store, mode) {
    key <- if (identical(mode, "selected")) "selected_enrichments_all_dbs" else "significant_enrichments_all_dbs"
    direct <- normalize_enrichment_df(store[[key]])
    if (!base::is.null(direct) && base::nrow(direct) > 0) {
      return(direct)
    }
    if (identical(mode, "selected") && base::is.list(store[["top_all_dbs"]])) {
      cand <- normalize_enrichment_df(store[["top_all_dbs"]][["result"]])
      if (!base::is.null(cand) && base::nrow(cand) > 0) {
        return(cand)
      }
    }
    per_db_keys <- base::grep("^top_", base::names(store), value = TRUE)
    per_db_keys <- per_db_keys[!per_db_keys %in% c("top_all_dbs", "top_all_dbs_mixed")]
    if (base::length(per_db_keys) > 0) {
      collected <- base::lapply(per_db_keys, function(k) {
        entry <- store[[k]]
        if (!base::is.list(entry)) {
          return(NULL)
        }
        subkey <- if (identical(mode, "selected")) "selected_enrichments" else "significant_enrichments"
        normalize_enrichment_df(entry[[subkey]])
      })
      collected <- collected[!base::vapply(collected, base::is.null, FUN.VALUE = base::logical(1))]
      if (base::length(collected) > 0) {
        out <- base::do.call(base::rbind, collected)
        base::rownames(out) <- NULL
        return(out)
      }
    }
    NULL
  }

  resolve_upstream_df <- function(store, mode) {
    key <- if (identical(mode, "selected")) "selected_upstream_all" else "significant_upstream_all"
    normalize_upstream_df(store[[key]])
  }

  enrich_store <- hcobject[["integrated_output"]][["enrichments"]]
  if (base::is.null(enrich_store)) {
    enrich_store <- hcobject[["satellite_outputs"]][["enrichments"]]
  }
  if (base::is.null(enrich_store)) {
    enrich_store <- list()
  }
  enrich_df <- resolve_enrichment_df(enrich_store, enrichment_mode)
  if (base::is.null(enrich_df) || base::nrow(enrich_df) == 0) {
    enrich_file <- base::paste0(file_prefix, "/Enrichment_Selected_All_DBs.xlsx")
    enrich_df <- normalize_enrichment_df(read_enrichment_xlsx_by_mode(enrich_file, enrichment_mode))
  }
  if (base::is.null(enrich_df) || base::nrow(enrich_df) == 0) {
    stop("Missing enrichment results. Run `.hc_functional_enrichment_driver()` first.")
  }

  upstream_store <- hcobject[["integrated_output"]][["upstream_inference"]]
  if (base::is.null(upstream_store)) {
    upstream_store <- hcobject[["satellite_outputs"]][["upstream_inference"]]
  }
  if (base::is.null(upstream_store)) {
    upstream_store <- list()
  }
  upstream_key <- if (identical(upstream_mode, "selected")) "selected_upstream_all" else "significant_upstream_all"
  upstream_df <- resolve_upstream_df(upstream_store, upstream_mode)
  if (base::is.null(upstream_df) || base::nrow(upstream_df) == 0) {
    upstream_file <- base::paste0(file_prefix, "/Upstream_Inference.xlsx")
    upstream_df <- normalize_upstream_df(read_xlsx_sheet(upstream_file, upstream_key))
  }
  if (base::is.null(upstream_df) || base::nrow(upstream_df) == 0) {
    stop("Missing upstream results. Run `.hc_upstream_inference_driver()` first.")
  }
  upstream_settings <- upstream_store[["settings"]]
  upstream_activity_input <- "gfc"
  if (base::is.list(upstream_settings) &&
    "activity_input" %in% base::names(upstream_settings) &&
    base::length(upstream_settings$activity_input) > 0 &&
    !base::is.na(upstream_settings$activity_input[[1]])) {
    upstream_activity_input <- base::tolower(base::as.character(upstream_settings$activity_input[[1]]))
  }
  if (!upstream_activity_input %in% c("gfc", "fc", "expression")) {
    upstream_activity_input <- "gfc"
  }
  upstream_module_heatmap_mat <- upstream_store[["module_heatmap_matrix"]]
  if (!base::is.null(upstream_module_heatmap_mat)) {
    upstream_module_heatmap_mat <- upstream_module_heatmap_mat %>% base::as.matrix()
    if (base::nrow(upstream_module_heatmap_mat) == 0 || base::ncol(upstream_module_heatmap_mat) == 0) {
      upstream_module_heatmap_mat <- NULL
    }
  }
  upstream_module_heatmap_col_order <- upstream_store[["module_heatmap_col_order"]]
  if (base::is.null(upstream_module_heatmap_col_order)) {
    upstream_module_heatmap_col_order <- base::character(0)
  } else {
    upstream_module_heatmap_col_order <- base::as.character(upstream_module_heatmap_col_order)
  }
  upstream_module_heatmap_name <- upstream_store[["module_heatmap_name"]]
  if (base::is.null(upstream_module_heatmap_name) || base::length(upstream_module_heatmap_name) == 0 || base::is.na(upstream_module_heatmap_name[[1]])) {
    upstream_module_heatmap_name <- if (identical(upstream_activity_input, "fc")) "FC" else "GFC"
  }
  upstream_module_heatmap_name <- base::as.character(upstream_module_heatmap_name[[1]])
  resolved_gfc_scale_limits <- normalize_scale_limits(gfc_scale_limits)
  if (base::is.null(resolved_gfc_scale_limits) &&
    base::is.list(upstream_settings) &&
    "gfc_scale_limits" %in% base::names(upstream_settings)) {
    resolved_gfc_scale_limits <- tryCatch(
      normalize_scale_limits(upstream_settings[["gfc_scale_limits"]]),
      error = function(e) NULL
    )
  }
  if (base::is.null(resolved_gfc_scale_limits)) {
    resolved_gfc_scale_limits <- tryCatch(
      normalize_scale_limits(cluster_calc[["gfc_scale_limits"]]),
      error = function(e) NULL
    )
  }
  if (base::is.null(resolved_gfc_scale_limits)) {
    fallback_lim <- .hc_first_numeric_value(hcobject[["global_settings"]][["range_GFC"]])
    if (!base::is.finite(fallback_lim) || fallback_lim <= 0) {
      fallback_lim <- 2
    }
    resolved_gfc_scale_limits <- c(-base::abs(fallback_lim), base::abs(fallback_lim))
  }
  resolved_gfc_scale_ticks <- compute_scale_ticks(resolved_gfc_scale_limits)

  build_fc_module_heatmap_from_settings <- function(fc_comparisons, cluster_order, cluster_info) {
    if (base::is.null(fc_comparisons) || base::length(fc_comparisons) == 0) {
      stop(
        "Upstream was run with `activity_input = 'fc'`, but no `fc_comparisons` are stored. ",
        "Please re-run `hc_upstream_inference(..., activity_input = 'fc', fc_comparisons = ...)`."
      )
    }
    merged_net <- hcobject[["integrated_output"]][["merged_net"]]
    if (base::is.null(merged_net)) {
      stop("Missing integrated network for FC heatmap reconstruction.")
    }
    net_genes <- igraph::get.vertex.attribute(merged_net, "name")
    net_genes <- base::unique(base::as.character(net_genes))
    net_genes <- net_genes[!base::is.na(net_genes) & net_genes != ""]
    if (base::length(net_genes) == 0) {
      stop("Could not extract genes from integrated network for FC heatmap reconstruction.")
    }
    parse_fc <- function(x) {
      x <- base::trimws(base::as.character(x))
      parts <- base::strsplit(x, "\\s*_vs_\\s*", perl = TRUE)[[1]]
      if (base::length(parts) != 2) {
        parts <- base::strsplit(x, "\\s+vs\\s+", perl = TRUE)[[1]]
      }
      if (base::length(parts) != 2) {
        stop("Invalid FC comparison format: '", x, "'. Use 'groupA_vs_groupB'.")
      }
      base::data.frame(
        comparison = base::paste0(base::trimws(parts[[1]]), "_vs_", base::trimws(parts[[2]])),
        numerator = base::trimws(parts[[1]]),
        denominator = base::trimws(parts[[2]]),
        stringsAsFactors = FALSE
      )
    }
    comparisons <- base::lapply(fc_comparisons, parse_fc)
    comparisons <- base::do.call(base::rbind, comparisons)
    comparisons <- comparisons[!duplicated(comparisons$comparison), , drop = FALSE]
    rownames(comparisons) <- NULL
    collapse_rows <- function(mat) {
      if (!base::anyDuplicated(base::rownames(mat))) {
        return(mat)
      }
      idx <- base::split(base::seq_len(base::nrow(mat)), base::rownames(mat))
      out <- base::t(base::vapply(
        idx,
        function(i) {
          vals <- base::colMeans(mat[i, , drop = FALSE], na.rm = TRUE)
          vals[base::is.nan(vals)] <- NA_real_
          vals
        },
        FUN.VALUE = base::numeric(base::ncol(mat))
      ))
      out <- out %>% base::as.matrix()
      out
    }
    collapse_cols <- function(mat) {
      if (!base::anyDuplicated(base::colnames(mat))) {
        return(mat)
      }
      idx <- base::split(base::seq_len(base::ncol(mat)), base::colnames(mat))
      out <- base::vapply(
        idx,
        function(i) {
          if (base::length(i) == 1) {
            return(mat[, i])
          }
          vals <- base::rowMeans(mat[, i, drop = FALSE], na.rm = TRUE)
          vals[base::is.nan(vals)] <- NA_real_
          vals
        },
        FUN.VALUE = base::numeric(base::nrow(mat))
      )
      out <- out %>% base::as.matrix()
      base::rownames(out) <- base::rownames(mat)
      out
    }
    resolve_group_labels <- function(info_dataset) {
      info_dataset <- info_dataset %>% base::as.data.frame(stringsAsFactors = FALSE)
      if (base::nrow(info_dataset) == 0) {
        return(base::character(0))
      }
      voi <- hcobject[["global_settings"]][["voi"]]
      voi <- base::intersect(voi, base::colnames(info_dataset))
      if (base::length(voi) > 0) {
        return(purrr::pmap(info_dataset[, voi, drop = FALSE], paste, sep = "-") %>% base::unlist())
      }
      base::as.character(info_dataset[[1]])
    }
    fc_limit <- .hc_first_numeric_value(hcobject[["global_settings"]][["range_GFC"]])
    if (!base::is.finite(fc_limit) || base::is.na(fc_limit) || fc_limit <= 0) {
      fc_limit <- 2
    }
    pseudo_count <- 1e-08
    set_indices <- base::seq_len(base::length(hcobject[["layer_specific_outputs"]]))
    set_mats <- list()
    found <- stats::setNames(base::rep(FALSE, base::nrow(comparisons)), comparisons$comparison)
    for (z in set_indices) {
      set_name <- base::paste0("set", z)
      expr_mat <- hcobject[["layer_specific_outputs"]][[set_name]][["part1"]][["topvar"]]
      if (base::is.null(expr_mat)) {
        expr_mat <- hcobject[["data"]][[base::paste0(set_name, "_counts")]]
      }
      anno <- hcobject[["data"]][[base::paste0(set_name, "_anno")]]
      if (base::is.null(expr_mat) || base::is.null(anno)) {
        next
      }
      expr_mat <- expr_mat %>% base::as.matrix()
      if (base::nrow(expr_mat) == 0 || base::ncol(expr_mat) == 0) {
        next
      }
      mode(expr_mat) <- "numeric"
      samples <- base::intersect(base::colnames(expr_mat), base::rownames(anno))
      if (base::length(samples) == 0) {
        next
      }
      expr_mat <- expr_mat[, samples, drop = FALSE]
      anno <- anno[samples, , drop = FALSE]
      grpvar <- resolve_group_labels(anno)
      grpvar <- base::as.character(grpvar)
      if (base::length(grpvar) != base::length(samples)) {
        next
      }
      bad <- base::is.na(grpvar) | grpvar == ""
      if (base::any(bad)) {
        grpvar[bad] <- samples[bad]
      }
      if (isTRUE(hcobject[["global_settings"]][["data_in_log"]])) {
        expr_mat <- .hc_antilog_impl(expr_mat, 2)
      }
      grp_levels <- base::unique(grpvar)
      set_mean <- base::vapply(
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
      set_mean <- set_mean %>% base::as.matrix()
      if (base::ncol(set_mean) == 0) {
        next
      }
      base::colnames(set_mean) <- grp_levels
      base::rownames(set_mean) <- base::rownames(expr_mat)
      overlap <- base::intersect(net_genes, base::rownames(set_mean))
      if (base::length(overlap) == 0) {
        next
      }
      for (i in base::seq_len(base::nrow(comparisons))) {
        num_grp <- comparisons$numerator[[i]]
        den_grp <- comparisons$denominator[[i]]
        cmp <- comparisons$comparison[[i]]
        if (!(num_grp %in% base::colnames(set_mean) && den_grp %in% base::colnames(set_mean))) {
          next
        }
        found[[cmp]] <- TRUE
        num_vals <- .hc_as_numeric_safely(set_mean[overlap, num_grp])
        den_vals <- .hc_as_numeric_safely(set_mean[overlap, den_grp])
        fc_vals <- .hc_log2_safely((num_vals + pseudo_count) / (den_vals + pseudo_count))
        fc_vals[!base::is.finite(fc_vals)] <- NA_real_
        fc_vals[fc_vals > fc_limit] <- fc_limit
        fc_vals[fc_vals < (-fc_limit)] <- -fc_limit
        set_full <- base::matrix(NA_real_, nrow = base::length(net_genes), ncol = 1, dimnames = list(net_genes, cmp))
        set_full[overlap, 1] <- fc_vals
        set_mats[[base::length(set_mats) + 1]] <- set_full
      }
    }
    if (base::length(set_mats) == 0 || !base::any(found)) {
      stop(
        "Could not reconstruct FC heatmap for knowledge network. ",
        "Please re-run `hc_upstream_inference(..., activity_input = 'fc', fc_comparisons = ...)` on this object."
      )
    }
    gene_fc <- base::do.call(base::cbind, set_mats)
    gene_fc <- collapse_rows(gene_fc)
    gene_fc <- collapse_cols(gene_fc)
    keep_rows <- base::rowSums(!base::is.na(gene_fc)) > 0
    gene_fc <- gene_fc[keep_rows, , drop = FALSE]
    requested_order <- base::as.character(comparisons$comparison)
    col_order <- requested_order[requested_order %in% base::colnames(gene_fc)]
    col_order <- base::c(col_order, base::setdiff(base::colnames(gene_fc), col_order))
    gene_fc <- gene_fc[, col_order, drop = FALSE]

    module_out <- list()
    for (cl in cluster_order) {
      genes <- dplyr::filter(cluster_info, color == cl) %>%
        dplyr::pull(., "gene_n") %>%
        base::strsplit(split = ",") %>%
        base::unlist()
      genes <- base::unique(base::as.character(genes))
      genes <- base::intersect(genes, base::rownames(gene_fc))
      if (base::length(genes) == 0) {
        next
      }
      vals <- gene_fc[genes, , drop = FALSE]
      m <- base::colMeans(vals, na.rm = TRUE)
      m[base::is.nan(m)] <- NA_real_
      module_out[[cl]] <- m
    }
    if (base::length(module_out) == 0) {
      stop("No module-level FC means could be computed for the knowledge network heatmap.")
    }
    module_mat <- base::do.call(base::rbind, module_out)
    module_mat <- module_mat %>% base::as.matrix()
    list(mat = module_mat, col_order = col_order)
  }

  if (identical(upstream_activity_input, "fc")) {
    if (base::is.null(upstream_module_heatmap_name) || base::length(upstream_module_heatmap_name) == 0) {
      upstream_module_heatmap_name <- "FC"
    }
    upstream_module_heatmap_name <- "FC"
    if (base::is.null(upstream_module_heatmap_mat) || base::nrow(upstream_module_heatmap_mat) == 0 || base::ncol(upstream_module_heatmap_mat) == 0) {
      rebuilt_fc <- build_fc_module_heatmap_from_settings(
        fc_comparisons = if (base::is.list(upstream_settings)) upstream_settings[["fc_comparisons"]] else NULL,
        cluster_order = cluster_order,
        cluster_info = cluster_info
      )
      upstream_module_heatmap_mat <- rebuilt_fc$mat
      upstream_module_heatmap_col_order <- rebuilt_fc$col_order
    }
  }

  enrich_df <- enrich_df[base::as.character(enrich_df$cluster) %in% cluster_order, , drop = FALSE]
  upstream_df <- upstream_df[base::as.character(upstream_df$cluster) %in% cluster_order, , drop = FALSE]
  if (base::nrow(enrich_df) == 0 && base::nrow(upstream_df) == 0) {
    stop("No enrichment/upstream rows remain after applying cluster filter.")
  }

  clamp_q <- function(x) {
    x <- .hc_as_numeric_safely(x)
    x[base::is.na(x) | x <= 0] <- 1e-300
    x
  }

  limit_per_module <- function(df, cluster_col, q_col, n_max) {
    if (base::is.null(n_max) || base::nrow(df) == 0) {
      return(df)
    }
    split_idx <- base::split(base::seq_len(base::nrow(df)), base::as.character(df[[cluster_col]]))
    kept <- base::lapply(split_idx, function(idx) {
      sub <- df[idx, , drop = FALSE]
      qv <- clamp_q(sub[[q_col]])
      ord <- base::order(qv)
      sub[ord[base::seq_len(base::min(base::length(ord), n_max))], , drop = FALSE]
    })
    out <- base::do.call(base::rbind, kept)
    base::rownames(out) <- NULL
    out
  }

  if (isTRUE(collapse_redundant_terms) && base::nrow(enrich_df) > 0) {
    enrich_df <- .hc_kn_collapse_terms(enrich_df, contained = 0.8)
  }
  enrich_df <- limit_per_module(enrich_df, "cluster", "qvalue", max_enrichment_per_module)
  upstream_df <- limit_per_module(upstream_df, "cluster", "qvalue", max_upstream_per_module)
  if (base::nrow(enrich_df) == 0 && base::nrow(upstream_df) == 0) {
    stop("No rows remain after applying per-module limits.")
  }

  cluster_idx_map <- stats::setNames(base::seq_along(cluster_order), cluster_order)
  x_module <- 0
  x_enrichment <- 1.05
  x_upstream <- 2.75
  x_limits <- c(-0.25, 3.75)
  module_nodes <- base::data.frame(
    id = base::paste0("M::", cluster_order),
    key = cluster_order,
    label = base::as.character(module_label_map[cluster_order]),
    node_type = "module",
    x = x_module,
    y = base::rev(base::seq_len(base::length(cluster_order))),
    node_color = cluster_order,
    stringsAsFactors = FALSE
  )
  module_y_map <- stats::setNames(module_nodes$y, module_nodes$key)

  interleaved_order <- function(meta_df, group_col) {
    if (base::is.null(meta_df) || base::nrow(meta_df) == 0) {
      return(base::character(0))
    }
    group_levels <- base::unique(base::as.character(meta_df[[group_col]]))
    bucket_levels <- base::sort(base::unique(meta_df$first_cluster_idx))
    ordered <- base::character(0)
    for (b in bucket_levels) {
      bucket <- meta_df[meta_df$first_cluster_idx == b, , drop = FALSE]
      group_lists <- base::lapply(group_levels, function(g) {
        x <- bucket[base::as.character(bucket[[group_col]]) == g, , drop = FALSE]
        if (base::nrow(x) == 0) {
          return(base::character(0))
        }
        x <- x[base::order(x$best_q, -x$hit_count, x$mean_cluster_idx, x$term), , drop = FALSE]
        x$node_key
      })
      base::names(group_lists) <- group_levels
      max_len <- base::max(base::lengths(group_lists))
      for (k in base::seq_len(max_len)) {
        for (g in group_levels) {
          if (base::length(group_lists[[g]]) >= k) {
            ordered <- base::c(ordered, group_lists[[g]][k])
          }
        }
      }
    }
    ordered <- base::unique(ordered)
    base::c(ordered, base::setdiff(meta_df$node_key, ordered))
  }

  enrich_nodes <- base::data.frame(
    id = base::character(0),
    key = base::character(0),
    label = base::character(0),
    node_type = base::character(0),
    node_subtype = base::character(0),
    x = base::numeric(0),
    y = base::numeric(0),
    stringsAsFactors = FALSE
  )
  enrich_edges <- base::data.frame(
    from_id = base::character(0),
    to_id = base::character(0),
    from_key = base::character(0),
    to_key = base::character(0),
    edge_type = base::character(0),
    edge_subtype = base::character(0),
    qvalue = base::numeric(0),
    direction = base::character(0),
    edge_from = base::character(0),
    x = base::numeric(0),
    y = base::numeric(0),
    xend = base::numeric(0),
    yend = base::numeric(0),
    stringsAsFactors = FALSE
  )
  if (base::nrow(enrich_df) > 0) {
    enrich_df$database <- base::as.character(enrich_df$database)
    enrich_df$cluster <- base::as.character(enrich_df$cluster)
    enrich_df$term <- base::as.character(enrich_df$term)
    enrich_df$qvalue <- clamp_q(enrich_df$qvalue)
    enrich_df$cluster_idx <- cluster_idx_map[enrich_df$cluster]
    enrich_df$cluster_idx[base::is.na(enrich_df$cluster_idx)] <- base::length(cluster_order) + 1
    enrich_df$node_key <- base::paste(enrich_df$database, enrich_df$term, sep = "||")

    split_en <- base::split(base::seq_len(base::nrow(enrich_df)), enrich_df$node_key)
    enrich_meta <- base::do.call(base::rbind, base::lapply(base::names(split_en), function(k) {
      idx <- split_en[[k]]
      sub <- enrich_df[idx, , drop = FALSE]
      base::data.frame(
        node_key = k,
        database = sub$database[[1]],
        term = sub$term[[1]],
        first_cluster_idx = base::min(sub$cluster_idx, na.rm = TRUE),
        mean_cluster_idx = base::mean(sub$cluster_idx, na.rm = TRUE),
        best_q = base::min(sub$qvalue, na.rm = TRUE),
        hit_count = base::length(base::unique(sub$cluster)),
        stringsAsFactors = FALSE
      )
    }))
    enrich_meta <- enrich_meta[base::order(enrich_meta$first_cluster_idx, enrich_meta$best_q, -enrich_meta$hit_count, enrich_meta$term), , drop = FALSE]
    en_order <- interleaved_order(enrich_meta, "database")
    enrich_meta <- enrich_meta[base::match(en_order, enrich_meta$node_key), , drop = FALSE]
    n_en <- base::nrow(enrich_meta)
    y_en <- if (n_en <= 1) {
      base::mean(module_nodes$y)
    } else {
      base::seq(base::max(module_nodes$y), base::min(module_nodes$y), length.out = n_en)
    }
    enrich_nodes <- base::data.frame(
      id = base::paste0("E::", enrich_meta$node_key),
      key = enrich_meta$node_key,
      label = base::paste0("[", enrich_meta$database, "] ", enrich_meta$term),
      node_type = "enrichment",
      node_subtype = enrich_meta$database,
      x = x_enrichment,
      y = y_en,
      stringsAsFactors = FALSE
    )
    enrich_y_map <- stats::setNames(enrich_nodes$y, enrich_nodes$key)
    enrich_edges <- base::data.frame(
      from_id = base::paste0("M::", enrich_df$cluster),
      to_id = base::paste0("E::", enrich_df$node_key),
      from_key = enrich_df$cluster,
      to_key = enrich_df$node_key,
      edge_type = "enrichment",
      edge_subtype = enrich_df$database,
      qvalue = enrich_df$qvalue,
      direction = NA_character_,
      edge_from = "module",
      stringsAsFactors = FALSE
    )
    enrich_edges$x <- x_module
    enrich_edges$y <- module_y_map[enrich_edges$from_key]
    enrich_edges$xend <- x_enrichment
    enrich_edges$yend <- enrich_y_map[enrich_edges$to_key]
  }

  links <- base::data.frame()
  upstream_nodes <- base::data.frame(
    id = base::character(0),
    key = base::character(0),
    label = base::character(0),
    node_type = base::character(0),
    node_subtype = base::character(0),
    x = base::numeric(0),
    y = base::numeric(0),
    stringsAsFactors = FALSE
  )
  upstream_edges <- base::data.frame(
    from_id = base::character(0),
    to_id = base::character(0),
    from_key = base::character(0),
    to_key = base::character(0),
    edge_type = base::character(0),
    edge_subtype = base::character(0),
    qvalue = base::numeric(0),
    direction = base::character(0),
    edge_from = base::character(0),
    x = base::numeric(0),
    y = base::numeric(0),
    xend = base::numeric(0),
    yend = base::numeric(0),
    stringsAsFactors = FALSE
  )
  if (base::nrow(upstream_df) > 0) {
    upstream_df$resource <- base::as.character(upstream_df$resource)
    upstream_df$cluster <- base::as.character(upstream_df$cluster)
    upstream_df$term <- base::as.character(upstream_df$term)
    upstream_df$direction <- base::as.character(upstream_df$direction)
    if ("regulation" %in% base::colnames(upstream_df)) {
      # How the regulator acts on its targets in the module; the direction of
      # change belongs to the per-condition heatmaps.
      reg <- base::as.character(upstream_df$regulation)
      upstream_df$direction <- base::ifelse(
        reg == "activating", "activated",
        base::ifelse(reg == "repressing", "inhibited", "mixed")
      )
    }
    if (!"regulator_module" %in% base::colnames(upstream_df)) {
      upstream_df$regulator_module <- NA_character_
    }
    upstream_df$qvalue <- clamp_q(upstream_df$qvalue)
    upstream_df$cluster_idx <- cluster_idx_map[upstream_df$cluster]
    upstream_df$cluster_idx[base::is.na(upstream_df$cluster_idx)] <- base::length(cluster_order) + 1
    upstream_df$node_key <- base::paste(upstream_df$resource, upstream_df$term, sep = "||")

    split_up <- base::split(base::seq_len(base::nrow(upstream_df)), upstream_df$node_key)
    upstream_meta <- base::do.call(base::rbind, base::lapply(base::names(split_up), function(k) {
      idx <- split_up[[k]]
      sub <- upstream_df[idx, , drop = FALSE]
      base::data.frame(
        node_key = k,
        resource = sub$resource[[1]],
        term = sub$term[[1]],
        first_cluster_idx = base::min(sub$cluster_idx, na.rm = TRUE),
        mean_cluster_idx = base::mean(sub$cluster_idx, na.rm = TRUE),
        best_q = base::min(sub$qvalue, na.rm = TRUE),
        hit_count = base::length(base::unique(sub$cluster)),
        regulator_module = {
          rm <- base::as.character(sub$regulator_module)
          rm <- rm[!base::is.na(rm) & base::nzchar(rm)]
          if (base::length(rm) > 0) rm[[1]] else NA_character_
        },
        stringsAsFactors = FALSE
      )
    }))
    upstream_meta <- upstream_meta[base::order(upstream_meta$first_cluster_idx, upstream_meta$best_q, -upstream_meta$hit_count, upstream_meta$term), , drop = FALSE]
    up_order <- interleaved_order(upstream_meta, "resource")
    upstream_meta <- upstream_meta[base::match(up_order, upstream_meta$node_key), , drop = FALSE]
    n_up <- base::nrow(upstream_meta)
    y_up <- if (n_up <= 1) {
      base::mean(module_nodes$y)
    } else {
      base::seq(base::max(module_nodes$y), base::min(module_nodes$y), length.out = n_up)
    }
    upstream_nodes <- base::data.frame(
      id = base::paste0("U::", upstream_meta$node_key),
      key = upstream_meta$node_key,
      label = base::paste0(upstream_meta$term, " [", upstream_meta$resource, "]"),
      node_type = "upstream",
      node_subtype = upstream_meta$resource,
      x = x_upstream,
      y = y_up,
      regulator_module = upstream_meta$regulator_module,
      stringsAsFactors = FALSE
    )
    upstream_y_map <- stats::setNames(upstream_nodes$y, upstream_nodes$key)
    upstream_edges <- base::data.frame(
      from_id = base::paste0("M::", upstream_df$cluster),
      to_id = base::paste0("U::", upstream_df$node_key),
      from_key = upstream_df$cluster,
      to_key = upstream_df$node_key,
      edge_type = "upstream",
      edge_subtype = upstream_df$resource,
      qvalue = upstream_df$qvalue,
      direction = upstream_df$direction,
      edge_from = "module",
      stringsAsFactors = FALSE
    )
    upstream_edges$x <- x_module
    upstream_edges$y <- module_y_map[upstream_edges$from_key]
    upstream_edges$xend <- x_upstream
    upstream_edges$yend <- upstream_y_map[upstream_edges$to_key]

    # Term -> regulator lines where the regulator's targets in the module
    # overlap the term; such regulators are no longer linked to the module
    # directly.
    links <- if (base::nrow(enrich_df) > 0) {
      .hc_kn_term_regulator_links(enrich_df, upstream_df, min_share = link_min_share, min_genes = 3)
    } else {
      base::data.frame()
    }
    if (base::nrow(links) > 0) {
      up_row <- base::match(base::paste(links$cluster, links$up_key), base::paste(upstream_df$cluster, upstream_df$node_key))
      link_edges <- base::data.frame(
        from_id = base::paste0("E::", links$term_key),
        to_id = base::paste0("U::", links$up_key),
        from_key = links$cluster,
        to_key = links$up_key,
        edge_type = "upstream",
        edge_subtype = upstream_df$resource[up_row],
        qvalue = upstream_df$qvalue[up_row],
        direction = upstream_df$direction[up_row],
        edge_from = "term",
        stringsAsFactors = FALSE
      )
      link_edges$x <- x_enrichment
      link_edges$y <- enrich_y_map[links$term_key]
      link_edges$xend <- x_upstream
      link_edges$yend <- upstream_y_map[link_edges$to_key]
      linked <- base::unique(base::paste(links$cluster, links$up_key))
      upstream_edges <- upstream_edges[!base::paste(upstream_edges$from_key, upstream_edges$to_key) %in% linked, , drop = FALSE]
      upstream_edges <- base::rbind(upstream_edges, link_edges)
    }
  }

  edges <- base::rbind(enrich_edges, upstream_edges)
  if (base::nrow(edges) == 0) {
    stop("No edges to plot after filtering.")
  }

  edges$qvalue <- clamp_q(edges$qvalue)
  edges$neglog10_q <- -base::log10(base::pmax(edges$qvalue, 1e-300))
  edges$neglog10_q[!base::is.finite(edges$neglog10_q)] <- 0
  edges$neglog10_q <- base::pmin(12, edges$neglog10_q)
  edges$edge_alpha <- 0.24 + (0.56 * base::pmin(1, edges$neglog10_q / 6))

  normalize_db <- function(x) {
    x <- base::as.character(x)
    low <- base::tolower(x)
    out <- x
    out[low %in% c("go", "gobp", "go_bp")] <- "Go"
    out[low %in% c("kegg")] <- "Kegg"
    out[low %in% c("hallmark")] <- "Hallmark"
    out[low %in% c("reactome")] <- "Reactome"
    out
  }
  edges$edge_subtype <- normalize_db(edges$edge_subtype)

  to_dir <- function(x) {
    x <- base::tolower(base::as.character(x))
    if (x %in% c("activated", "up", "positive")) {
      return("activated")
    }
    if (x %in% c("inhibited", "down", "negative")) {
      return("inhibited")
    }
    ""
  }

  edges$direction_class <- ifelse(
    edges$edge_type == "upstream",
    base::vapply(edges$direction, to_dir, FUN.VALUE = base::character(1)),
    ""
  )

  edges$line_group <- base::vapply(base::seq_len(base::nrow(edges)), function(i) {
    if (identical(edges$edge_type[[i]], "enrichment")) {
      return(base::paste0(edges$edge_subtype[[i]], " enrichment"))
    }
    d <- edges$direction_class[[i]]
    if (d == "activated") {
      return(base::paste0(edges$edge_subtype[[i]], " activated"))
    }
    if (d == "inhibited") {
      return(base::paste0(edges$edge_subtype[[i]], " inhibited"))
    }
    base::paste0(edges$edge_subtype[[i]], " (no direction)")
  }, FUN.VALUE = base::character(1))

  line_color_map <- c(
    "Go enrichment" = "#4E79A7",
    "Kegg enrichment" = "#F28E2B",
    "Hallmark enrichment" = "#59A14F",
    "Reactome enrichment" = "#E15759",
    "TF activated" = "#D62728",
    "TF inhibited" = "#1F77B4",
    "Pathway activated" = "#C03D3E",
    "Pathway inhibited" = "#2C73B8",
    "TF (no direction)" = "#3B7EA1",
    "Pathway (no direction)" = "#B26B2C"
  )
  missing_groups <- base::setdiff(base::unique(edges$line_group), base::names(line_color_map))
  if (base::length(missing_groups) > 0) {
    line_color_map <- base::c(
      line_color_map,
      stats::setNames(base::rep("grey65", base::length(missing_groups)), missing_groups)
    )
  }
  line_group_order <- c(
    "Go enrichment",
    "Kegg enrichment",
    "Hallmark enrichment",
    "Reactome enrichment",
    "TF activated",
    "TF inhibited",
    "Pathway activated",
    "Pathway inhibited",
    "TF (no direction)",
    "Pathway (no direction)"
  )
  line_group_order <- base::c(
    line_group_order[line_group_order %in% base::names(line_color_map)],
    base::setdiff(base::names(line_color_map), line_group_order)
  )

  build_hc_heatmap_data <- function(cluster_order,
                                    module_label_map,
                                    cluster_info,
                                    stored_hm,
                                    main_heatmap_col_order = NULL,
                                    heatmap_col_order = NULL,
                                    heatmap_cluster_columns = FALSE,
                                    override_mat = NULL,
                                    override_col_order = base::character(0),
                                    value_name = "GFC",
                                    gfc_scale_limits = c(-2, 2),
                                    gfc_scale_ticks = NULL) {
    if (base::is.null(gfc_scale_ticks) ||
      !base::is.list(gfc_scale_ticks) ||
      !all(c("breaks", "labels") %in% base::names(gfc_scale_ticks))) {
      gfc_scale_ticks <- compute_scale_ticks(gfc_scale_limits)
    }

    extract_heatmap <- function(stored_hm, cluster_info, cluster_order, module_label_map, override_mat) {
      if (!base::is.null(override_mat) && base::nrow(override_mat) > 0 && base::ncol(override_mat) > 0) {
        hm_mat <- override_mat %>% base::as.matrix()
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
      gfc_all <- hcobject[["integrated_output"]][["GFC_all_layers"]]
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
      out <- base::lapply(cluster_order, function(cl) {
        genes <- dplyr::filter(cluster_info, color == cl) %>%
          dplyr::pull(., "gene_n") %>%
          base::strsplit(split = ",") %>%
          base::unlist()
        genes <- base::unique(base::as.character(genes))
        if (base::length(genes) == 0) {
          return(NULL)
        }
        tmp <- gfc_all[gfc_all[[gene_col]] %in% genes, value_idx, drop = FALSE]
        if (base::nrow(tmp) == 0) {
          return(NULL)
        }
        vals <- tmp %>% base::as.matrix()
        mode(vals) <- "numeric"
        base::colMeans(vals, na.rm = TRUE)
      })
      names(out) <- cluster_order
      out <- out[!base::vapply(out, base::is.null, FUN.VALUE = base::logical(1))]
      if (base::length(out) == 0) {
        return(NULL)
      }
      hm <- base::do.call(base::rbind, out)
      hm %>% base::as.matrix()
    }

    heatmap_mat <- extract_heatmap(stored_hm, cluster_info, cluster_order, module_label_map, override_mat)
    if (base::is.null(heatmap_mat) || base::nrow(heatmap_mat) == 0) {
      return(NULL)
    }

    prepared_cols <- .hc_prepare_plot_heatmap_columns(
      mat = heatmap_mat,
      cluster_columns = heatmap_cluster_columns,
      plot_order = heatmap_col_order,
      main_order = main_heatmap_col_order,
      fallback_order = override_col_order,
      context = "knowledge network heatmap"
    )
    heatmap_mat <- prepared_cols$mat

    keep_clusters <- cluster_order[cluster_order %in% base::rownames(heatmap_mat)]
    if (base::length(keep_clusters) == 0) {
      return(NULL)
    }
    heatmap_mat <- heatmap_mat[keep_clusters, , drop = FALSE]
    module_labels <- base::as.character(module_label_map[keep_clusters])
    module_labels[is.na(module_labels) | module_labels == ""] <- keep_clusters[is.na(module_labels) | module_labels == ""]
    base::rownames(heatmap_mat) <- module_labels

    cdf <- dplyr::filter(cluster_info, color %in% keep_clusters) %>%
      dplyr::distinct(color, .keep_all = TRUE)
    if ("gene_no" %in% base::colnames(cdf)) {
      gene_counts <- .hc_as_numeric_safely(cdf$gene_no)
    } else if ("gene_n" %in% base::colnames(cdf)) {
      gene_counts <- base::vapply(
        base::as.character(cdf$gene_n),
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
      gene_counts <- base::rep(0, base::nrow(cdf))
    }
    gene_counts <- stats::setNames(gene_counts, base::as.character(cdf$color))[keep_clusters]
    gene_counts[base::is.na(gene_counts)] <- 0
    value_name <- base::as.character(value_name[[1]])
    if (base::is.na(value_name) || value_name == "") {
      value_name <- "GFC"
    }
    value_name_upper <- base::toupper(value_name)
    if (!identical(value_name_upper, "FC")) {
      value_name <- "GFC"
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
    if (base::is.null(stored_gfc_colors)) {
      stored_gfc_colors <- .hc_default_gfc_colors()
    }

    list(
      mat = heatmap_mat,
      keep_clusters = keep_clusters,
      module_labels = module_labels,
      module_colors = keep_clusters,
      gene_counts = as.numeric(gene_counts),
      value_name = value_name,
      gfc_colors = stored_gfc_colors,
      scale_limits = gfc_scale_limits,
      scale_breaks = gfc_scale_ticks$breaks,
      scale_labels = gfc_scale_ticks$labels
    )
  }

  hc_heatmap_data <- build_hc_heatmap_data(
    cluster_order = cluster_order,
    module_label_map = module_label_map,
    cluster_info = cluster_info,
    stored_hm = stored_hm,
    main_heatmap_col_order = main_heatmap_col_order,
    heatmap_col_order = col_order,
    heatmap_cluster_columns = cluster_columns,
    override_mat = upstream_module_heatmap_mat,
    override_col_order = upstream_module_heatmap_col_order,
    value_name = upstream_module_heatmap_name,
    gfc_scale_limits = resolved_gfc_scale_limits,
    gfc_scale_ticks = resolved_gfc_scale_ticks
  )
  if (base::is.null(hc_heatmap_data)) {
    stop("Unable to build hCoCena heatmap panel for the network plot.")
  }
  # Force exact y-alignment between heatmap rows and module rows in the network.
  cluster_order <- base::as.character(hc_heatmap_data$keep_clusters)
  module_label_map <- stats::setNames(
    base::as.character(hc_heatmap_data$module_labels),
    cluster_order
  )
  module_nodes <- module_nodes[module_nodes$key %in% cluster_order, , drop = FALSE]
  module_nodes <- module_nodes[base::match(cluster_order, module_nodes$key), , drop = FALSE]
  module_nodes$label <- base::as.character(module_label_map[module_nodes$key])
  module_nodes$y <- base::rev(base::seq_len(base::length(cluster_order)))
  module_y_map <- stats::setNames(module_nodes$y, module_nodes$key)
  if (base::nrow(enrich_edges) > 0) {
    enrich_edges$y <- module_y_map[base::as.character(enrich_edges$from_key)]
  }
  from_module <- edges$edge_from == "module"
  edges$y[from_module] <- module_y_map[base::as.character(edges$from_key[from_module])]
  missing_edge_y <- base::is.na(edges$y) | !(base::as.character(edges$from_key) %in% cluster_order)
  if (base::any(missing_edge_y)) {
    edges <- edges[!missing_edge_y, , drop = FALSE]
    if (base::nrow(edges) == 0) {
      stop("No edges remain after aligning module rows with heatmap rows.")
    }
  }

  hm_cols <- base::colnames(hc_heatmap_data$mat)
  page_column_gap_spec <- .hc_heatmap_column_gap_spec(
    hcobject = hcobject,
    cols = hm_cols,
    cluster_columns = heatmap_cluster_columns,
    gap_mm = 0.6 * overall_plot_scale
  )
  page_column_layout <- .hc_heatmap_ggplot_column_layout(
    cols = hm_cols,
    column_gap_spec = page_column_gap_spec,
    default_cell_mm = 5
  )
  page_column_labels <- .hc_gfc_display_col_labels(hcobject, hm_cols)
  enrich_page <- if (base::nrow(enrich_df) > 0) enrich_df else
    base::data.frame(cluster = base::character(0), node_key = base::character(0), database = base::character(0),
                     term = base::character(0), qvalue = base::numeric(0))
  upstream_page <- if (base::nrow(upstream_df) > 0) upstream_df else
    base::data.frame(cluster = base::character(0), node_key = base::character(0), resource = base::character(0),
                     term = base::character(0), qvalue = base::numeric(0), direction = base::character(0),
                     regulator_module = base::character(0))
  overview_title <- "Module knowledge network"
  build_page <- function(focus_cluster = NULL) {
    .hc_kn_page_plot(
      hm = hc_heatmap_data,
      column_labels = page_column_labels,
      column_layout = page_column_layout,
      enrich_df = enrich_page[enrich_page$cluster %in% cluster_order, , drop = FALSE],
      upstream_df = upstream_page[upstream_page$cluster %in% cluster_order, , drop = FALSE],
      links = links,
      focus = focus_cluster,
      title = if (base::is.null(focus_cluster)) overview_title else
        base::paste0(overview_title, " - ", module_label_map[[focus_cluster]]),
      scale = overall_plot_scale
    )
  }
  overview_page <- build_page()
  overview_plot <- overview_page$plot
  focus_plots <- stats::setNames(
    lapply(cluster_order, function(cl) build_page(cl)$plot),
    base::as.character(module_label_map[cluster_order])
  )

  pdf_width_use <- if (base::is.null(pdf_width)) overview_page$width else as.numeric(pdf_width)
  pdf_height_use <- if (base::is.null(pdf_height)) overview_page$height else as.numeric(pdf_height)

  out_file <- NULL
  overview_png_files <- stats::setNames(base::character(0), base::character(0))
  focus_files <- stats::setNames(base::character(0), base::character(0))
  focus_png_files <- stats::setNames(base::character(0), base::character(0))
  if (isTRUE(save_pdf)) {
    out_file <- base::paste0(file_prefix, "/", pdf_name)
    overview_page_labels <- base::c("overview", base::names(focus_plots))
    overview_export_files <- .hc_export_multi_page_plot(
      file = out_file,
      page_labels = overview_page_labels,
      width = pdf_width_use,
      height = pdf_height_use,
      pointsize = pdf_pointsize,
      res = 300,
      draw_page_fun = function(idx, page_key) {
        if (identical(page_key, "overview")) {
          print(overview_plot)
        } else {
          print(focus_plots[[page_key]])
        }
      }
    )
    overview_png_files <- overview_export_files$png

    focus_dir <- base::paste0(file_prefix, "/Module_Knowledge_Network_by_module")
    if (!base::dir.exists(focus_dir)) {
      base::dir.create(focus_dir, recursive = TRUE, showWarnings = FALSE)
    }
    focus_names <- base::names(focus_plots)
    focus_files <- stats::setNames(base::character(base::length(focus_names)), focus_names)
    focus_png_files <- stats::setNames(base::character(base::length(focus_names)), focus_names)
    for (i in base::seq_along(focus_plots)) {
      mod_nm <- focus_names[[i]]
      mod_file <- base::paste0(
        focus_dir,
        "/Module_Knowledge_Network_",
        .hc_export_sanitize_stem(mod_nm, default = "module"),
        ".pdf"
      )
      focus_export_files <- .hc_export_single_page_plot(
        file = mod_file,
        width = pdf_width_use,
        height = pdf_height_use,
        pointsize = pdf_pointsize,
        res = 300,
        draw_fun = function() {
          print(focus_plots[[i]])
        }
      )
      focus_files[[i]] <- focus_export_files$pdf
      focus_png_files[[i]] <- focus_export_files$png
    }
  }

  if (isTRUE(show_plot)) {
    print(overview_plot)
    for (i in base::seq_along(focus_plots)) {
      print(focus_plots[[i]])
    }
  }

  nodes <- base::rbind(
    module_nodes[, c("id", "key", "label", "node_type", "x", "y"), drop = FALSE],
    if (base::nrow(enrich_nodes) > 0) {
      enrich_nodes[, c("id", "key", "label", "node_type", "x", "y"), drop = FALSE]
    } else {
      base::data.frame()
    },
    if (base::nrow(upstream_nodes) > 0) {
      upstream_nodes[, c("id", "key", "label", "node_type", "x", "y"), drop = FALSE]
    } else {
      base::data.frame()
    }
  )

  output <- list(
    plot = overview_plot,
    focus_plots = focus_plots,
    focus_files = focus_files,
    focus_png_files = focus_png_files,
    nodes = nodes,
    edges = edges,
    file = out_file,
    png_files = overview_png_files,
    settings = list(
      enrichment_mode = enrichment_mode,
      upstream_mode = upstream_mode,
      upstream_activity_input = upstream_activity_input,
      heatmap_value_name = hc_heatmap_data$value_name,
      heatmap_scale_limits = hc_heatmap_data$scale_limits,
      label_mode = label_mode,
      clusters = cluster_order,
      max_enrichment_per_module = max_enrichment_per_module,
      max_upstream_per_module = max_upstream_per_module
    )
  )
  .hc_set_bridge_hcobject_slot(c("satellite_outputs", "knowledge_network"), output)
  output
}




# Within each module, drop enrichment terms whose genes lie mostly (share
# >= `contained`) in a better-ranked term that is kept, e.g. KEGG cell cycle
# inside GO cell cycle. Terms without gene lists stay.
.hc_kn_collapse_terms <- function(enrich_df, contained = 0.8) {
  if (!"geneID" %in% base::colnames(enrich_df)) {
    return(enrich_df)
  }
  genes <- base::strsplit(base::as.character(enrich_df$geneID), "/", fixed = TRUE)
  keep <- base::rep(TRUE, base::nrow(enrich_df))
  for (cl in base::unique(base::as.character(enrich_df$cluster))) {
    idx <- base::which(base::as.character(enrich_df$cluster) == cl)
    idx <- idx[base::order(.hc_as_numeric_safely(enrich_df$qvalue[idx]))]
    kept <- base::integer(0)
    for (i in idx) {
      g <- genes[[i]]
      if (base::length(g) == 0 || base::all(base::is.na(g))) {
        next
      }
      inside <- base::any(base::vapply(kept, function(k) {
        base::length(base::intersect(g, genes[[k]])) / base::length(g) >= contained
      }, base::logical(1)))
      if (inside) keep[i] <- FALSE else kept <- c(kept, i)
    }
  }
  enrich_df[keep, , drop = FALSE]
}

# Term -> regulator links within a module. `share` is the fraction of the
# regulator's targets in the module that belong to the term; `ratio` compares
# it with the fraction of all module genes in the term, so that large generic
# terms are not linked to every regulator. Links need `min_genes` shared genes,
# `share >= min_share` and `ratio >= min_ratio`; per regulator and module the
# `top` links with the highest ratio are kept.
.hc_kn_term_regulator_links <- function(enrich_df, upstream_df, min_share = 0.25,
                                        min_genes = 3, min_ratio = 1.5, top = 3) {
  empty <- base::data.frame(cluster = base::character(0), term_key = base::character(0),
                            up_key = base::character(0), n_shared = base::integer(0),
                            share = base::numeric(0), ratio = base::numeric(0),
                            stringsAsFactors = FALSE)
  if (!"geneID" %in% base::colnames(enrich_df) || !"overlap_genes" %in% base::colnames(upstream_df)) {
    return(empty)
  }
  term_genes <- base::strsplit(base::as.character(enrich_df$geneID), "/", fixed = TRUE)
  reg_genes <- base::strsplit(base::as.character(upstream_df$overlap_genes), ",", fixed = TRUE)
  module_size <- if ("n_genes" %in% base::colnames(upstream_df)) {
    .hc_as_numeric_safely(upstream_df$n_genes)
  } else {
    base::rep(NA_real_, base::nrow(upstream_df))
  }
  out <- list()
  for (i in base::seq_len(base::nrow(upstream_df))) {
    targets <- reg_genes[[i]]
    targets <- targets[!base::is.na(targets) & base::nzchar(targets)]
    if (base::length(targets) == 0 || !base::is.finite(module_size[[i]])) next
    same <- base::which(base::as.character(enrich_df$cluster) == base::as.character(upstream_df$cluster[[i]]))
    rows <- base::lapply(same, function(j) {
      n <- base::length(base::intersect(targets, term_genes[[j]]))
      share <- n / base::length(targets)
      ratio <- share / (base::length(term_genes[[j]]) / module_size[[i]])
      if (n < min_genes || share < min_share || ratio < min_ratio) return(NULL)
      base::data.frame(
        cluster = base::as.character(upstream_df$cluster[[i]]),
        term_key = base::as.character(enrich_df$node_key[[j]]),
        up_key = base::as.character(upstream_df$node_key[[i]]),
        n_shared = n, share = share, ratio = ratio,
        stringsAsFactors = FALSE
      )
    })
    rows <- base::do.call(base::rbind, rows[!base::vapply(rows, base::is.null, base::logical(1))])
    if (!base::is.null(rows) && base::nrow(rows) > 0) {
      out[[base::length(out) + 1]] <- utils::head(rows[base::order(-rows$ratio), , drop = FALSE], top)
    }
  }
  if (base::length(out) == 0) empty else base::do.call(base::rbind, out)
}

# ---- page layout ----------------------------------------------------------
# One page of the knowledge network, drawn as a single ggplot with a fixed
# aspect ratio so that heatmap cells and module boxes stay square:
#   module heatmap | upstream regulators -> modules -> enriched terms
# Regulators point at the module they act on (arrow = activates its targets
# in the module, T-end = represses, plain = mixed). With `focus`, the other
# modules are greyed out and the regulator -> term links (shared genes) of the
# focus module are drawn as dashed curves.

.hc_kn_pretty_term <- function(term) {
  term <- base::as.character(term)
  term <- base::sub("^(HALLMARK|GOBP|GOCC|GOMF|GO|KEGG|REACTOME|WP|BIOCARTA|PID)_", "", term)
  base::gsub("_", " ", term, fixed = TRUE)
}

.hc_kn_trim <- function(x, n) {
  x <- base::as.character(x)
  long <- base::nchar(x) > n
  x[long] <- base::paste0(base::substr(x[long], 1, n - 1), "\u2026")
  x
}

# Smooth S-shaped path between two points: horizontal at both ends, so arrow
# heads and T-ends sit straight on the target.
.hc_kn_curves <- function(x0, y0, x1, y1, id, n = 30) {
  if (base::length(x0) == 0) {
    return(base::data.frame(x = base::numeric(0), y = base::numeric(0), id = base::character(0)))
  }
  t <- base::seq(0, 1, length.out = n)
  s <- 3 * t^2 - 2 * t^3
  base::do.call(base::rbind, base::lapply(base::seq_along(x0), function(i) {
    base::data.frame(
      x = x0[[i]] + (x1[[i]] - x0[[i]]) * t,
      y = y0[[i]] + (y1[[i]] - y0[[i]]) * s,
      id = id[[i]],
      stringsAsFactors = FALSE
    )
  }))
}

.hc_kn_spread <- function(n, height) {
  if (n <= 1) {
    return(base::rep(height / 2, n))
  }
  base::seq(height - 0.5, 0.5, length.out = n)
}

.hc_kn_page_plot <- function(hm, column_labels, column_layout, enrich_df, upstream_df,
                             links, focus = NULL, title = NULL, scale = 1) {
  mat <- hm$mat
  n_r <- base::nrow(mat)
  n_c <- base::ncol(mat)
  clusters <- base::as.character(hm$keep_clusters)
  labels <- stats::setNames(base::as.character(hm$module_labels), clusters)
  if (!base::is.null(focus) && !focus %in% clusters) focus <- NULL
  db_colors <- c(Go = "#4E79A7", Kegg = "#F28E2B", Hallmark = "#59A14F", Reactome = "#B07AA1")
  reg_colors <- c(activates = "#C0392B", represses = "#2471A3", mixed = "#8C8C8C")

  # -- nodes -----------------------------------------------------------------
  terms <- if (base::nrow(enrich_df) > 0) {
    t <- enrich_df[!base::duplicated(enrich_df$node_key), c("node_key", "database", "term"), drop = FALSE]
    first <- base::tapply(base::match(enrich_df$cluster, clusters), enrich_df$node_key, base::min)
    bestq <- base::tapply(enrich_df$qvalue, enrich_df$node_key, base::min)
    t$first <- first[t$node_key]
    t$bestq <- bestq[t$node_key]
    t[base::order(t$first, t$bestq), , drop = FALSE]
  } else {
    base::data.frame(node_key = base::character(0), database = base::character(0), term = base::character(0))
  }
  regs <- if (base::nrow(upstream_df) > 0) {
    r <- upstream_df[!base::duplicated(upstream_df$node_key), c("node_key", "resource", "term"), drop = FALSE]
    first <- base::tapply(base::match(upstream_df$cluster, clusters), upstream_df$node_key, base::min)
    bestq <- base::tapply(upstream_df$qvalue, upstream_df$node_key, base::min)
    in_mod <- base::tapply(base::as.character(upstream_df$regulator_module), upstream_df$node_key, function(v) {
      v <- v[!base::is.na(v) & base::nzchar(v)]
      if (base::length(v) == 0) NA_character_ else v[[1]]
    })
    r$first <- first[r$node_key]
    r$bestq <- bestq[r$node_key]
    r$home <- base::names(labels)[base::match(in_mod[r$node_key], labels)]
    r[base::order(r$first, r$resource != "TF", r$bestq), , drop = FALSE]
  } else {
    base::data.frame(node_key = base::character(0), resource = base::character(0), term = base::character(0), home = base::character(0))
  }

  height <- base::max(n_r, 0.5 * base::nrow(terms), 0.55 * base::nrow(regs))
  offset <- (height - n_r) / 2
  row_y <- stats::setNames(offset + n_r - base::seq_len(n_r) + 0.5, clusters)
  char_u <- 0.215
  reg_lab <- base::paste0(regs$term, base::ifelse(regs$resource == "TF", "", " (pathway)"))
  term_lab <- .hc_kn_trim(.hc_kn_pretty_term(terms$term), 40)
  hm_right <- base::max(column_layout$x)
  x_reg <- hm_right + 2.6 + char_u * base::max(c(8, base::nchar(reg_lab))) + 0.5
  x_mod <- x_reg + 4.5
  x_term <- x_mod + 4.5
  x_max <- x_term + 0.9 + 0.245 * base::max(c(10, base::nchar(term_lab)))
  regs$y <- .hc_kn_spread(base::nrow(regs), height)
  terms$y <- .hc_kn_spread(base::nrow(terms), height)
  reg_y <- stats::setNames(regs$y, regs$node_key)
  term_y <- stats::setNames(terms$y, terms$node_key)

  is_focus_cluster <- function(cl) base::is.null(focus) | cl %in% focus
  q_width <- function(q) {
    v <- base::pmin(10, -base::log10(base::pmax(.hc_as_numeric_safely(q), 1e-300)))
    0.35 + 1.6 * v / 10
  }

  # -- edges -------------------------------------------------------------------
  up <- upstream_df
  up$group <- base::ifelse(up$direction == "activated", "activates",
                           base::ifelse(up$direction == "inhibited", "represses", "mixed"))
  up$focus <- is_focus_cluster(up$cluster)
  up$id <- base::paste0("U", base::seq_len(base::nrow(up)))
  up_x1 <- x_mod - 0.62
  reg_paths <- .hc_kn_curves(base::rep(x_reg + 0.25, base::nrow(up)), reg_y[up$node_key],
                             base::rep(up_x1, base::nrow(up)), row_y[up$cluster], up$id)
  reg_paths <- base::merge(reg_paths, up[, c("id", "group", "focus", "qvalue")], by = "id", sort = FALSE)
  reg_paths$lw <- q_width(reg_paths$qvalue)

  en <- enrich_df
  en$group <- base::ifelse(en$database %in% base::names(db_colors), en$database, "Other")
  en$focus <- is_focus_cluster(en$cluster)
  en$id <- base::paste0("E", base::seq_len(base::nrow(en)))
  term_paths <- .hc_kn_curves(base::rep(x_mod + 0.5, base::nrow(en)), row_y[en$cluster],
                              base::rep(x_term - 0.25, base::nrow(en)), term_y[en$node_key], en$id)
  term_paths <- base::merge(term_paths, en[, c("id", "group", "focus", "qvalue")], by = "id", sort = FALSE)
  term_paths$lw <- q_width(term_paths$qvalue)

  link_paths <- base::data.frame()
  if (!base::is.null(focus) && base::nrow(links) > 0) {
    lk <- links[links$cluster == focus & links$term_key %in% terms$node_key & links$up_key %in% regs$node_key, , drop = FALSE]
    if (base::nrow(lk) > 0) {
      lk$id <- base::paste0("L", base::seq_len(base::nrow(lk)))
      link_paths <- .hc_kn_curves(base::rep(x_reg + 0.25, base::nrow(lk)), reg_y[lk$up_key],
                                  base::rep(x_term - 0.25, base::nrow(lk)), term_y[lk$term_key], lk$id, n = 40)
    }
  }

  grey_out <- function(df) {
    if (base::nrow(df) == 0) return(df)
    df$alpha <- base::ifelse(df$focus, 0.85, 0.08)
    df
  }
  reg_paths <- grey_out(reg_paths)
  term_paths <- grey_out(term_paths)
  # draw focused edges on top
  reg_paths <- reg_paths[base::order(reg_paths$focus), , drop = FALSE]
  term_paths <- term_paths[base::order(term_paths$focus), , drop = FALSE]

  tbars <- up[up$group == "represses", , drop = FALSE]
  tbars <- base::data.frame(
    x = base::rep(up_x1, base::nrow(tbars)),
    y = row_y[tbars$cluster] - 0.22,
    yend = row_y[tbars$cluster] + 0.22,
    group = tbars$group,
    alpha = base::ifelse(tbars$focus, 0.95, 0.08),
    stringsAsFactors = FALSE
  )

  # -- node attributes -----------------------------------------------------
  focus_regs <- base::unique(up$node_key[up$focus])
  focus_terms <- base::unique(en$node_key[en$focus])
  regs$active <- regs$node_key %in% focus_regs
  terms$active <- terms$node_key %in% focus_terms
  # Node colour = how the regulator acts on its targets (as the links);
  # "mixed" when it differs between modules.
  effect <- base::tapply(up$group, up$node_key, function(g) {
    g <- base::unique(g)
    if (base::length(g) == 1) g else "mixed"
  })
  regs$effect <- base::as.character(effect[regs$node_key])
  regs$fill <- base::ifelse(regs$active, reg_colors[regs$effect], "grey92")
  # Regulator genes that belong to a module get a tag in that module's colour.
  tags <- regs[!base::is.na(regs$home), , drop = FALSE]
  tags$label <- labels[tags$home]
  tags$fill <- base::ifelse(tags$active, tags$home, "grey88")
  tags$text <- base::ifelse(tags$active, "white", "grey60")
  reg_label_x <- base::ifelse(base::is.na(regs$home), x_reg - 0.4, x_reg - 1.15)
  regs$shape <- base::ifelse(regs$resource == "TF", 21, 23)
  terms$fill <- base::ifelse(terms$active,
                             base::ifelse(terms$database %in% base::names(db_colors), db_colors[terms$database], "grey50"),
                             "grey90")
  mod <- base::data.frame(cluster = clusters, label = labels[clusters], y = row_y[clusters], stringsAsFactors = FALSE)
  mod$fill <- base::ifelse(is_focus_cluster(mod$cluster), mod$cluster, "grey88")
  mod$text <- base::ifelse(is_focus_cluster(mod$cluster), "white", "grey60")
  counts <- .hc_as_numeric_safely(hm$gene_counts)
  mod$genes <- if (base::length(counts) == base::length(clusters)) {
    base::formatC(counts, format = "d", big.mark = ",")
  } else {
    ""
  }
  lab_col <- function(active) base::ifelse(active, "grey10", "grey75")

  hm_long <- base::data.frame(
    x = base::rep(column_layout$x - 0.5, each = n_r),
    y = base::rep(row_y[clusters], times = n_c),
    value = base::as.vector(mat),
    stringsAsFactors = FALSE
  )
  gfc_pal <- grDevices::colorRampPalette(hm$gfc_colors)(51)

  # -- legends via invisible keys --------------------------------------------
  edge_values <- c(reg_colors, db_colors, Other = "grey50")
  edge_labels <- c(
    activates = "Regulator activates its targets  ->",
    represses = "Regulator represses its targets  -|",
    mixed = "Regulator, mixed effect",
    Go = "GO term", Kegg = "KEGG pathway", Hallmark = "Hallmark gene set",
    Reactome = "Reactome pathway", Other = "Other term"
  )
  used_groups <- base::intersect(base::names(edge_values), base::unique(c(up$group, en$group)))
  key_df <- base::data.frame(
    x = x_mod, y = height / 2,
    node = c("Module", "TF", "Signalling pathway (PROGENy)", "TF gene lies in this module"),
    stringsAsFactors = FALSE
  )

  fs <- 7.2 * scale
  p <- ggplot2::ggplot() +
    ggplot2::geom_tile(data = hm_long, ggplot2::aes(x = x, y = y, fill = value),
                       width = 0.94, height = 0.94, color = NA) +
    ggplot2::scale_fill_gradientn(
      colors = gfc_pal, limits = hm$scale_limits, breaks = hm$scale_breaks,
      labels = hm$scale_labels, oob = scales::squish, name = hm$value_name,
      guide = ggplot2::guide_colorbar(order = 1, barheight = grid::unit(28 * scale, "mm"),
                                      barwidth = grid::unit(3.2 * scale, "mm"))
    ) +
    ggplot2::geom_tile(data = mod, ggplot2::aes(x = hm_right + 0.35, y = y), width = 0.3, height = 0.94,
                       fill = mod$fill, color = NA) +
    ggplot2::geom_text(data = mod, ggplot2::aes(x = hm_right + 0.62, y = y, label = genes),
                       hjust = 0, size = fs * 0.9 / ggplot2::.pt,
                       color = base::ifelse(is_focus_cluster(mod$cluster), "grey25", "grey70")) +
    ggplot2::annotate("text", x = hm_right + 0.62, y = offset + n_r + 0.3, label = "genes",
                      hjust = 0, vjust = 0, size = fs * 0.85 / ggplot2::.pt, color = "grey40",
                      fontface = "italic") +
    ggplot2::geom_path(data = term_paths, ggplot2::aes(x = x, y = y, group = id, color = group,
                                                       linewidth = lw, alpha = alpha),
                       lineend = "round") +
    {
      if (base::nrow(link_paths) > 0) {
        ggplot2::geom_path(data = link_paths, ggplot2::aes(x = x, y = y, group = id),
                           color = "grey45", linewidth = 0.35, linetype = "22", alpha = 0.8)
      }
    } +
    ggplot2::geom_path(data = reg_paths[reg_paths$group != "activates", , drop = FALSE],
                       ggplot2::aes(x = x, y = y, group = id, color = group, linewidth = lw, alpha = alpha),
                       lineend = "round") +
    ggplot2::geom_path(data = reg_paths[reg_paths$group == "activates", , drop = FALSE],
                       ggplot2::aes(x = x, y = y, group = id, color = group, linewidth = lw, alpha = alpha),
                       lineend = "round",
                       arrow = grid::arrow(type = "closed", length = grid::unit(1.9 * scale, "mm"))) +
    ggplot2::geom_segment(data = tbars, ggplot2::aes(x = x, xend = x, y = y, yend = yend, color = group, alpha = alpha),
                          linewidth = 1.1, show.legend = FALSE) +
    ggplot2::scale_color_manual(values = edge_values, breaks = used_groups, labels = edge_labels[used_groups],
                                name = "Regulator effect / links", guide = ggplot2::guide_legend(order = 2, override.aes = list(linewidth = 1.4, alpha = 1))) +
    ggplot2::scale_linewidth_identity() +
    ggplot2::scale_alpha_identity() +
    ggplot2::geom_tile(data = mod, ggplot2::aes(x = x_mod, y = y), width = 0.94, height = 0.94,
                       fill = mod$fill, color = "grey20", linewidth = 0.25) +
    ggplot2::geom_text(data = mod, ggplot2::aes(x = x_mod, y = y, label = label),
                       color = mod$text, fontface = "bold", size = fs * 0.95 / ggplot2::.pt) +
    ggplot2::geom_point(data = regs, ggplot2::aes(x = x_reg, y = y), shape = regs$shape,
                        fill = regs$fill, color = base::ifelse(regs$active, "grey15", "grey80"),
                        size = 2.8 * scale, stroke = 0.4) +
    ggplot2::geom_tile(data = tags, ggplot2::aes(x = x_reg - 0.72, y = y), width = 0.62, height = 0.34,
                       fill = tags$fill, color = NA) +
    ggplot2::geom_text(data = tags, ggplot2::aes(x = x_reg - 0.72, y = y, label = label),
                       color = tags$text, fontface = "bold", size = fs * 0.72 / ggplot2::.pt) +
    ggplot2::geom_text(data = regs, ggplot2::aes(x = reg_label_x, y = y, label = reg_lab),
                       hjust = 1, size = fs / ggplot2::.pt, color = lab_col(regs$active)) +
    ggplot2::geom_point(data = terms, ggplot2::aes(x = x_term, y = y), shape = 21,
                        fill = terms$fill, color = "white", size = 2.3 * scale, stroke = 0.3) +
    ggplot2::geom_text(data = terms, ggplot2::aes(x = x_term + 0.4, y = y, label = term_lab),
                       hjust = 0, size = fs / ggplot2::.pt, color = lab_col(terms$active)) +
    ggplot2::annotate("text", x = column_layout$x - 0.5, y = offset - 0.3, label = column_labels,
                      angle = 90, hjust = 1, size = fs / ggplot2::.pt, color = "grey15") +
    ggplot2::annotate("text", x = c(base::mean(c(0, hm_right)), x_reg, x_mod, x_term),
                      y = height + 0.9,
                      label = c(hm$value_name, "Upstream regulators", "Modules", "Enriched terms"),
                      hjust = c(0.5, 1, 0.5, 0), fontface = "bold", size = fs * 1.15 / ggplot2::.pt) +
    ggplot2::geom_point(data = key_df, ggplot2::aes(x = x, y = y, shape = node), alpha = 0) +
    ggplot2::scale_shape_manual(
      name = "Nodes", values = c(Module = 22, TF = 21, `Signalling pathway (PROGENy)` = 23, `TF gene lies in this module` = 22),
      breaks = key_df$node,
      guide = ggplot2::guide_legend(order = 3, override.aes = list(
        alpha = 1, size = c(3.2, 3.2, 3.2, 2.4), fill = c("grey55", "grey85", "grey85", "grey55"),
        color = c("grey15", "grey15", "grey15", NA)
      ))
    ) +
    ggplot2::coord_fixed(ratio = 1, xlim = c(-0.3, x_max), ylim = c(-0.3, height + 1.4), clip = "off", expand = FALSE) +
    ggplot2::theme_void(base_size = 9 * scale) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 11.5 * scale, hjust = 0,
                                         margin = ggplot2::margin(b = 6)),
      plot.subtitle = ggplot2::element_text(size = 7.8 * scale, color = "grey35", hjust = 0,
                                            margin = ggplot2::margin(b = 4)),
      legend.title = ggplot2::element_text(face = "bold", size = 8.5 * scale),
      legend.text = ggplot2::element_text(size = 7.8 * scale),
      legend.key.height = grid::unit(4.2 * scale, "mm"),
      legend.spacing.y = grid::unit(2 * scale, "mm"),
      plot.margin = ggplot2::margin(10, 10, 6 + 4.3 * base::max(base::nchar(column_labels)) * scale, 10),
      plot.background = ggplot2::element_rect(fill = "white", color = NA)
    ) +
    ggplot2::labs(
      title = title,
      subtitle = if (base::nrow(link_paths) > 0) {
        "Dashed lines: the regulator's targets in this module are over-represented among the term's genes."
      } else {
        NULL
      }
    )

  unit_in <- 0.27 * scale
  list(
    plot = p,
    width = (x_max + 0.6) * unit_in + 3.4 * scale,
    height = (height + 1.7) * unit_in + 0.55 + 0.062 * base::max(base::nchar(column_labels)) * scale
  )
}
