.hc_normalize_longitudinal_k <- function(k, arg_name) {
  if (is.null(k)) {
    return(NULL)
  }
  k <- base::as.integer(k)
  k <- k[base::is.finite(k) & !base::is.na(k)]
  k <- base::sort(base::unique(k))
  if (base::length(k) == 0) {
    stop("`", arg_name, "` must contain at least one finite integer.")
  }
  if (base::any(k < 2)) {
    stop("`", arg_name, "` must contain integers >= 2.")
  }
  k
}

.hc_normalize_time_levels <- function(time_levels) {
  if (missing(time_levels)) {
    stop("`time_levels` is required and must be an explicit ordered timepoint vector.")
  }
  if (!base::is.atomic(time_levels) || base::length(time_levels) == 0) {
    stop("`time_levels` must be a non-empty vector.")
  }
  time_levels <- base::as.character(time_levels)
  time_levels <- time_levels[!base::is.na(time_levels) & base::nzchar(time_levels)]
  if (base::length(time_levels) == 0) {
    stop("`time_levels` must contain at least one non-empty value.")
  }
  base::unique(time_levels)
}

.hc_longitudinal_unwrap_hc <- function(hc, arg_name = "hc") {
  if (inherits(hc, "HCoCenaExperiment")) {
    return(hc)
  }
  if (base::is.list(hc) && !base::is.null(hc$hc) && inherits(hc$hc, "HCoCenaExperiment")) {
    return(hc$hc)
  }
  stop(
    "`", arg_name, "` must be a `HCoCenaExperiment` or a longitudinal-step result list containing `$hc`."
  )
}

.hc_longitudinal_quick_preset_defaults <- function(preset = c("standard", "consensus")) {
  preset <- base::match.arg(preset)
  if (identical(preset, "standard")) {
    return(list(
      module_k = 2:4,
      clustering_nstart = 100,
      cap_runs = 50,
      cap_nstart = 30,
      meta_k = 2:8,
      meta_consensus_runs = 250,
      meta_consensus_sample_fraction = 0.8,
      meta_consensus_feature_fraction = 0.85,
      include_consensus = FALSE,
      cap_donor_order = "none",
      meta_save_tables = TRUE
    ))
  }

  list(
    module_k = 2:4,
    clustering_nstart = 100,
    cap_runs = 50,
    cap_nstart = 30,
    meta_k = 3:8,
    meta_consensus_runs = 500,
    meta_consensus_sample_fraction = 0.9,
    meta_consensus_feature_fraction = 0.9,
    include_consensus = TRUE,
    cap_donor_order = "module_cluster_heatmap",
    meta_save_tables = FALSE
  )
}

.hc_longitudinal_step1_target_layers <- function(hc, layer = NULL) {
  exp_names <- base::names(MultiAssayExperiment::experiments(hc@mae))
  if (base::length(exp_names) == 0) {
    stop("No layers found in `hc@mae`. Run `hc_read_data()` first.")
  }

  if (base::is.null(layer)) {
    return(.hc_resolve_layer_id(hc = hc, layer = NULL))
  }

  if (base::is.character(layer) && base::length(layer) == 1 && identical(layer, "all")) {
    return(exp_names)
  }

  if (base::length(layer) > 1) {
    out <- base::vapply(
      base::seq_along(layer),
      function(i) .hc_resolve_layer_id(hc = hc, layer = layer[[i]]),
      FUN.VALUE = base::character(1)
    )
    return(base::unique(out))
  }

  .hc_resolve_layer_id(hc = hc, layer = layer)
}

.hc_longitudinal_layer_label <- function(hc, layer_id) {
  cfg <- hc@config@layer
  if (base::nrow(cfg) > 0 &&
    base::all(c("layer_id", "layer_name") %in% base::colnames(cfg))) {
    idx <- base::which(base::as.character(cfg$layer_id) == base::as.character(layer_id))
    if (base::length(idx) > 0) {
      lbl <- base::as.character(cfg$layer_name[[idx[[1]]]])
      if (!base::is.na(lbl) && base::nzchar(lbl)) {
        return(lbl)
      }
    }
  }
  base::as.character(layer_id[[1]])
}

.hc_longitudinal_safe_suffix <- function(x) {
  x <- base::trimws(base::as.character(x[[1]]))
  x <- gsub("[^A-Za-z0-9]+", "_", x, perl = TRUE)
  x <- gsub("^_+|_+$", "", x, perl = TRUE)
  if (!base::nzchar(x)) {
    return("layer")
  }
  x
}

#' Longitudinal step 1: module/donor clustering
#'
#' Computes longitudinal module means, then clusters donors separately within
#' each module directly on the donor-by-time trajectory matrix. The downstream
#' quick steps remain unchanged and can be run via
#' [hc_longitudinal_step2_meta_clustering()],
#' [hc_longitudinal_step3_meta_module_trajectories()], and
#' [hc_plot_longitudinal_enrichment_meta_waves()].
#'
#' @param hc A `HCoCenaExperiment`.
#' @param donor_col Annotation column containing donor IDs.
#' @param time_col Annotation column containing timepoint labels.
#' @param layer Optional layer index, layer id, layer name, or `"all"`. `NULL`
#'   uses the first available layer.
#' @param time_levels Required explicit timepoint order vector.
#' @param k Candidate donor-trajectory cluster numbers per module.
#' @param method One of `"kmeans"` or `"hclust"` for direct clustering on the
#'   donor-by-time matrix.
#' @param nstart Number of random starts for direct `kmeans`.
#' @param cap_runs Number of repeated clustering runs used to build CAP.
#' @param impute Logical. If `TRUE`, impute missing donor-time values before
#'   direct clustering.
#' @param impute_method Missing-value method passed to `mice`.
#' @param ntree Number of trees for `rfcont` imputation.
#' @param min_cluster_fraction Minimum allowed fraction of donors in the
#'   smallest module cluster.
#' @param score_method One of `"calinski_harabasz"` or `"tot_withinss"` for
#'   direct selection.
#' @param scale_features Logical. If `TRUE`, z-score timepoint columns before
#'   clustering.
#' @param na_impute One of `"median"` or `"zero"` for final NA handling after
#'   optional imputation.
#' @param cap_scale_features Logical. If `TRUE`, z-score features for CAP
#'   reruns.
#' @param cap_na_impute One of `"median"` or `"zero"` for CAP feature handling.
#' @param cap_jitter_sd Optional Gaussian jitter SD added per CAP rerun when
#'   `method = "kmeans"`.
#' @param seed Random seed.
#' @param means_slot Satellite slot for module means.
#' @param output_slot Satellite slot for endotype outputs.
#'
#' @return A list with updated `hc`, step plots, and diagnostics. When
#'   `layer = "all"` (or multiple layers are supplied), returns per-layer plot
#'   and diagnostic lists keyed by layer id.
#' @examples
#' hc <- hc_example_data("clustered")
#' res <- hc_longitudinal_step1_module_donor(
#'   hc,
#'   donor_col = "donor",
#'   time_col = "timepoint",
#'   time_levels = c("T1", "T2"),
#'   k = 2,
#'   nstart = 1,
#'   cap_runs = 1,
#'   impute = FALSE,
#'   seed = 1
#' )
#' hc <- res$hc
#' @export
hc_longitudinal_step1_module_donor <- function(hc,
                                               donor_col = "Subject",
                                               time_col = "Time_token",
                                               layer = NULL,
                                               time_levels,
                                               k = 2:6,
                                               method = c("kmeans", "hclust"),
                                               nstart = 50,
                                               cap_runs = 50,
                                               impute = TRUE,
                                               impute_method = "rfcont",
                                               ntree = 10,
                                               min_cluster_fraction = 0.1,
                                               score_method = c("calinski_harabasz", "tot_withinss"),
                                               scale_features = FALSE,
                                               na_impute = c("median", "zero"),
                                               cap_scale_features = scale_features,
                                               cap_na_impute = na_impute,
                                               cap_jitter_sd = 0,
                                               seed = 42,
                                               means_slot = "longitudinal_module_means",
                                               output_slot = "longitudinal_endotypes") {
  hc <- .hc_longitudinal_unwrap_hc(hc)
  time_levels <- .hc_normalize_time_levels(time_levels)
  k <- .hc_normalize_longitudinal_k(k, "k")
  method <- base::match.arg(method)
  score_method <- base::match.arg(score_method)
  na_impute <- base::match.arg(na_impute)
  # Explicit choices: the formal default is `na_impute`, which has already
  # been collapsed to a single value by the line above. Without them an
  # explicitly supplied vector - as passed down by
  # hc_longitudinal_workflow_direct() - is matched against one choice only
  # and match.arg() aborts with "'arg' must be of length 1".
  cap_na_impute <- base::match.arg(cap_na_impute, choices = c("median", "zero"))

  if (!base::is.numeric(nstart) || base::length(nstart) != 1 || !base::is.finite(nstart) || nstart < 1) {
    stop("`nstart` must be a single integer >= 1.")
  }
  if (!base::is.numeric(cap_runs) || base::length(cap_runs) != 1 || !base::is.finite(cap_runs) || cap_runs < 1) {
    stop("`cap_runs` must be a single integer >= 1.")
  }
  if (!base::is.logical(impute) || base::length(impute) != 1 || base::is.na(impute)) {
    stop("`impute` must be TRUE or FALSE.")
  }
  if (!base::is.character(impute_method) || base::length(impute_method) != 1 || !base::nzchar(impute_method)) {
    stop("`impute_method` must be a non-empty character scalar.")
  }
  if (isTRUE(impute) && identical(impute_method, "rfcont")) {
    .hc_require_namespace("CALIBERrfimpute", "`rfcont` longitudinal imputation")
  }
  if (!base::is.numeric(ntree) || base::length(ntree) != 1 || !base::is.finite(ntree) || ntree < 1) {
    stop("`ntree` must be a single integer >= 1.")
  }
  if (!base::is.numeric(min_cluster_fraction) || base::length(min_cluster_fraction) != 1 ||
    !base::is.finite(min_cluster_fraction) || min_cluster_fraction < 0 || min_cluster_fraction > 1) {
    stop("`min_cluster_fraction` must be between 0 and 1.")
  }
  if (!base::is.logical(scale_features) || base::length(scale_features) != 1 || base::is.na(scale_features)) {
    stop("`scale_features` must be TRUE or FALSE.")
  }
  if (!base::is.logical(cap_scale_features) || base::length(cap_scale_features) != 1 || base::is.na(cap_scale_features)) {
    stop("`cap_scale_features` must be TRUE or FALSE.")
  }
  if (!base::is.numeric(cap_jitter_sd) || base::length(cap_jitter_sd) != 1 || !base::is.finite(cap_jitter_sd) || cap_jitter_sd < 0) {
    stop("`cap_jitter_sd` must be a single numeric value >= 0.")
  }

  requested_all_layers <- base::is.character(layer) &&
    base::length(layer) == 1 &&
    identical(layer, "all")
  target_layers <- .hc_longitudinal_step1_target_layers(hc = hc, layer = layer)
  multi_layer <- base::length(target_layers) > 1

  if (!isTRUE(multi_layer) && !isTRUE(requested_all_layers)) {
    return(.hc_longitudinal_step1_run_single_direct(
      hc = hc,
      donor_col = donor_col,
      time_col = time_col,
      layer_id = target_layers[[1]],
      time_levels = time_levels,
      k = k,
      method = method,
      nstart = base::as.integer(nstart),
      cap_runs = base::as.integer(cap_runs),
      impute = impute,
      impute_method = impute_method,
      ntree = base::as.integer(ntree),
      min_cluster_fraction = min_cluster_fraction,
      score_method = score_method,
      scale_features = scale_features,
      na_impute = na_impute,
      cap_scale_features = cap_scale_features,
      cap_na_impute = cap_na_impute,
      cap_jitter_sd = cap_jitter_sd,
      seed = seed,
      means_slot = means_slot,
      output_slot = output_slot,
      module_means_prefix = "Longitudinal_ModuleMeans",
      module_clusters_prefix = "Longitudinal_ModuleClusters",
      cap_prefix = "Longitudinal_CAP"
    ))
  }

  plots_by_layer <- list()
  diagnostics_by_layer <- list()
  layer_info <- base::vector("list", base::length(target_layers))

  for (i in base::seq_along(target_layers)) {
    lid <- target_layers[[i]]
    layer_label <- .hc_longitudinal_layer_label(hc = hc, layer_id = lid)
    layer_suffix <- .hc_longitudinal_safe_suffix(lid)
    means_slot_i <- base::paste0(means_slot, "_", layer_suffix)
    output_slot_i <- base::paste0(output_slot, "_", layer_suffix)
    file_suffix <- .hc_longitudinal_safe_suffix(layer_label)

    res_i <- .hc_longitudinal_step1_run_single_direct(
      hc = hc,
      donor_col = donor_col,
      time_col = time_col,
      layer_id = lid,
      time_levels = time_levels,
      k = k,
      method = method,
      nstart = base::as.integer(nstart),
      cap_runs = base::as.integer(cap_runs),
      impute = impute,
      impute_method = impute_method,
      ntree = base::as.integer(ntree),
      min_cluster_fraction = min_cluster_fraction,
      score_method = score_method,
      scale_features = scale_features,
      na_impute = na_impute,
      cap_scale_features = cap_scale_features,
      cap_na_impute = cap_na_impute,
      cap_jitter_sd = cap_jitter_sd,
      seed = seed,
      means_slot = means_slot_i,
      output_slot = output_slot_i,
      module_means_prefix = base::paste0("Longitudinal_ModuleMeans_", file_suffix),
      module_clusters_prefix = base::paste0("Longitudinal_ModuleClusters_", file_suffix),
      cap_prefix = base::paste0("Longitudinal_CAP_", file_suffix)
    )
    hc <- res_i$hc

    plots_by_layer[[lid]] <- res_i$plots
    diagnostics_by_layer[[lid]] <- res_i$diagnostics
    layer_info[[i]] <- base::data.frame(
      layer_id = lid,
      layer_name = layer_label,
      means_slot = means_slot_i,
      output_slot = output_slot_i,
      stringsAsFactors = FALSE
    )
  }

  list(
    hc = hc,
    plots = plots_by_layer,
    diagnostics = diagnostics_by_layer,
    layer_info = base::do.call(base::rbind, layer_info)
  )
}

.hc_longitudinal_has_cap_matrix_slot <- function(obj) {
  base::is.list(obj) && !base::is.null(obj$cap_matrix)
}

.hc_longitudinal_step2_target_slots <- function(hc, slot_name = "longitudinal_endotypes") {
  sat <- as.list(hc@satellite)
  sat_names <- base::names(sat)
  if (base::is.null(sat_names) || base::length(sat_names) == 0) {
    stop("No satellite outputs found in `hc@satellite`. Run longitudinal step 1 first.")
  }

  cap_slots <- sat_names[base::vapply(sat, .hc_longitudinal_has_cap_matrix_slot, FUN.VALUE = base::logical(1))]
  if (base::length(cap_slots) == 0) {
    stop("No longitudinal step 1 outputs with `cap_matrix` found in `hc@satellite`.", call. = FALSE)
  }

  if (base::length(slot_name) > 1) {
    out <- base::unlist(
      lapply(
        base::as.list(slot_name),
        function(x) .hc_longitudinal_step2_target_slots(hc = hc, slot_name = x)
      ),
      use.names = FALSE
    )
    return(base::unique(out))
  }

  slot_name <- base::as.character(slot_name[[1]])
  if (!base::nzchar(slot_name)) {
    stop("`slot_name` must be a non-empty character value.")
  }

  if (identical(slot_name, "all")) {
    return(cap_slots)
  }

  family_prefix <- base::paste0(slot_name, "_")
  family_matches <- cap_slots[base::startsWith(cap_slots, family_prefix)]
  if (base::length(family_matches) > 0) {
    return(family_matches)
  }

  if (slot_name %in% cap_slots) {
    return(slot_name)
  }

  stop(
    "Could not resolve `slot_name` to a longitudinal step 1 output. ",
    "Checked exact slot `", slot_name, "` and family matches like `", slot_name, "_*`.",
    call. = FALSE
  )
}

.hc_longitudinal_requested_slot_family <- function(slot_name) {
  if (base::is.null(slot_name) || base::length(slot_name) != 1) {
    return(FALSE)
  }
  slot_name <- base::as.character(slot_name[[1]])
  if (!base::nzchar(slot_name)) {
    return(FALSE)
  }
  identical(slot_name, "all") ||
    identical(slot_name, "longitudinal_endotypes") ||
    base::startsWith(slot_name, "longitudinal_endotypes_")
}

.hc_longitudinal_step2_slot_context <- function(hc, slot_name, requested_slot_name = NULL) {
  sat <- as.list(hc@satellite)
  obj <- sat[[slot_name]]
  source_slot <- if (!base::is.null(obj$source_slot)) {
    base::as.character(obj$source_slot[[1]])
  } else {
    NA_character_
  }

  suffix_candidates <- base::character(0)
  if (!base::is.null(requested_slot_name) &&
    !identical(requested_slot_name, "all") &&
    base::length(requested_slot_name) == 1) {
    requested_slot_name <- base::as.character(requested_slot_name[[1]])
    req_prefix <- base::paste0(requested_slot_name, "_")
    if (base::startsWith(slot_name, req_prefix)) {
      suffix_candidates <- c(suffix_candidates, base::substring(slot_name, base::nchar(req_prefix) + 1L))
    }
  }
  if (base::startsWith(slot_name, "longitudinal_endotypes_")) {
    suffix_candidates <- c(
      suffix_candidates,
      base::substring(slot_name, base::nchar("longitudinal_endotypes_") + 1L)
    )
  }
  if (!base::is.na(source_slot) && base::startsWith(source_slot, "longitudinal_module_means_")) {
    suffix_candidates <- c(
      suffix_candidates,
      base::substring(source_slot, base::nchar("longitudinal_module_means_") + 1L)
    )
  }
  suffix_candidates <- base::unique(suffix_candidates[!base::is.na(suffix_candidates) & base::nzchar(suffix_candidates)])

  layer_id <- NA_character_
  layer_name <- slot_name
  output_name <- slot_name
  if (base::length(suffix_candidates) > 0) {
    layer_id <- base::as.character(suffix_candidates[[1]])
    layer_name <- .hc_longitudinal_layer_label(hc = hc, layer_id = layer_id)
    output_name <- layer_id
  }

  list(
    output_name = output_name,
    layer_id = layer_id,
    layer_name = layer_name,
    source_slot = source_slot,
    file_suffix = .hc_longitudinal_safe_suffix(layer_name)
  )
}

.hc_longitudinal_step2_run_single <- function(hc,
                                              slot_name,
                                              dimensions,
                                              graph_method,
                                              knn_method,
                                              graph_k,
                                              resolution,
                                              leiden_method,
                                              cluster_prefix,
                                              compute_umap,
                                              umap_neighbors,
                                              umap_min_dist,
                                              seed,
                                              file_prefix = "Longitudinal_Meta") {
  hc <- .hc_run_longitudinal_step2_graph(
    hc = hc,
    slot_name = slot_name,
    dimensions = as.integer(dimensions),
    graph_method = graph_method,
    knn_method = knn_method,
    graph_k = as.integer(graph_k),
    resolution = as.numeric(resolution),
    leiden_method = leiden_method,
    cluster_prefix = cluster_prefix,
    seed = seed,
    compute_umap = compute_umap,
    umap_neighbors = umap_neighbors,
    umap_min_dist = umap_min_dist
  )

  meta_embeddings <- hc_plot_longitudinal_meta_embeddings(
    hc,
    slot_name = slot_name,
    save_pdf = TRUE,
    file_prefix = file_prefix,
    show_endotype_crosstab = FALSE,
    show_cluster_labels = TRUE,
    save_tables = TRUE,
    table_format = "xlsx",
    table_detail = "simple"
  )
  meta_obj <- hc@satellite[[slot_name]]

  list(
    hc = hc,
    plots = list(
      pca = meta_embeddings$pca,
      umap = meta_embeddings$umap
    ),
    diagnostics = list(
      k_criterion = NULL,
      meta_consensus = NULL,
      cross_tab = meta_embeddings$cross_tab,
      tables = meta_embeddings$tables,
      graph_method = graph_method,
      knn_method = knn_method,
      graph_k = as.integer(graph_k),
      resolution = as.numeric(resolution),
      dimensions = as.integer(dimensions),
      leiden_method = leiden_method,
      meta_cluster = meta_obj$meta_cluster,
      meta_method_comparison = meta_obj$meta_method_comparison,
      meta_score_table = meta_obj$meta_score_table
    )
  )
}

.hc_longitudinal_has_meta_cluster_slot <- function(obj) {
  base::is.list(obj) && !base::is.null(obj$meta_cluster)
}

.hc_longitudinal_step3_target_slots <- function(hc, slot_name = "longitudinal_endotypes") {
  sat <- as.list(hc@satellite)
  sat_names <- base::names(sat)
  if (base::is.null(sat_names) || base::length(sat_names) == 0) {
    stop("No satellite outputs found in `hc@satellite`. Run longitudinal step 2 first.")
  }

  meta_slots <- sat_names[base::vapply(sat, .hc_longitudinal_has_meta_cluster_slot, FUN.VALUE = base::logical(1))]
  if (base::length(meta_slots) == 0) {
    stop("No longitudinal step 2 outputs with `meta_cluster` found in `hc@satellite`.", call. = FALSE)
  }

  if (base::length(slot_name) > 1) {
    out <- base::unlist(
      lapply(
        base::as.list(slot_name),
        function(x) .hc_longitudinal_step3_target_slots(hc = hc, slot_name = x)
      ),
      use.names = FALSE
    )
    return(base::unique(out))
  }

  slot_name <- base::as.character(slot_name[[1]])
  if (!base::nzchar(slot_name)) {
    stop("`slot_name` must be a non-empty character value.")
  }

  if (identical(slot_name, "all")) {
    return(meta_slots)
  }

  family_prefix <- base::paste0(slot_name, "_")
  family_matches <- meta_slots[base::startsWith(meta_slots, family_prefix)]
  if (base::length(family_matches) > 0) {
    return(family_matches)
  }

  if (slot_name %in% meta_slots) {
    return(slot_name)
  }

  stop(
    "Could not resolve `slot_name` to a longitudinal step 2 output. ",
    "Checked exact slot `", slot_name, "` and family matches like `", slot_name, "_*`.",
    call. = FALSE
  )
}

.hc_longitudinal_step3_run_single <- function(hc,
                                              slot_name,
                                              facet_ncol,
                                              free_y,
                                              square_panels,
                                              value_mode,
                                              value_range,
                                              file_prefix = "Longitudinal_Meta_ModuleWaves") {
  meta_module_waves <- hc_plot_longitudinal_meta_module_waves(
    hc,
    slot_name = slot_name,
    save_pdf = TRUE,
    file_prefix = file_prefix,
    facet_ncol = facet_ncol,
    free_y = free_y,
    square_panels = square_panels,
    value_mode = value_mode,
    value_range = value_range,
    save_width = 10,
    save_height = 10
  )

  list(
    hc = hc,
    plots = list(
      meta_module_waves = meta_module_waves$meta_module_waves
    )
  )
}

#' Longitudinal step 2: meta-clustering
#'
#' Runs the exact CAP -> PCA -> graph -> Leiden donor meta-clustering and
#' returns the PCA/UMAP embeddings.
#'
#' @param hc A `HCoCenaExperiment`.
#' @param slot_name Satellite slot containing step 1 outputs. Use a concrete
#'   slot name, a shared prefix such as `"longitudinal_endotypes"` to process
#'   matching suffixed slots, or `"all"` to run over every compatible step 1
#'   slot in `hc@satellite`.
#' @param dimensions Number of PCA dimensions used to build the donor graph.
#' @param graph_method Graph type used before Leiden clustering.
#' @param knn_method Nearest-neighbor backend used for graph construction.
#' @param graph_k Number of neighbors used for the donor graph.
#' @param resolution Leiden resolution parameter.
#' @param leiden_method Leiden partition method.
#' @param cluster_prefix Prefix for meta-cluster labels.
#' @param compute_umap Logical. If `TRUE`, compute an additional UMAP from the
#'   PCA subspace used for graph clustering.
#' @param umap_neighbors UMAP neighborhood size.
#' @param umap_min_dist UMAP minimum distance.
#' @param seed Random seed.
#'
#' @return A list with updated `hc`, step plots, and diagnostics. When multiple
#'   slots are processed, `plots` and `diagnostics` are returned as named lists
#'   keyed by layer/slot, plus a `slot_info` table describing the resolved
#'   slots.
#' @examples
#' hc <- hc_example_data("clustered")
#' hc <- hc_longitudinal_step1_module_donor(
#'   hc,
#'   donor_col = "donor",
#'   time_col = "timepoint",
#'   time_levels = c("T1", "T2"),
#'   k = 2,
#'   nstart = 1,
#'   cap_runs = 1,
#'   impute = FALSE,
#'   seed = 1
#' )$hc
#' res <- hc_longitudinal_step2_meta_clustering(hc, graph_k = 5, dimensions = 2)
#' hc <- res$hc
#' @export
hc_longitudinal_step2_meta_clustering <- function(hc,
                                                  slot_name = "longitudinal_endotypes",
                                                  dimensions = 4,
                                                  graph_method = c("knn", "snn"),
                                                  knn_method = c("annoy", "rann"),
                                                  graph_k = 7,
                                                  resolution = 0.4,
                                                  leiden_method = "RBConfigurationVertexPartition",
                                                  cluster_prefix = "MC",
                                                  compute_umap = TRUE,
                                                  umap_neighbors = 15,
                                                  umap_min_dist = 0.3,
                                                  seed = 42) {
  hc <- .hc_longitudinal_unwrap_hc(hc)
  graph_method <- base::match.arg(graph_method)
  knn_method <- base::match.arg(knn_method)
  if (!base::is.numeric(dimensions) || base::length(dimensions) != 1 || !base::is.finite(dimensions) || dimensions < 1) {
    stop("`dimensions` must be a single integer >= 1.")
  }
  if (!base::is.numeric(graph_k) || base::length(graph_k) != 1 || !base::is.finite(graph_k) || graph_k < 1) {
    stop("`graph_k` must be a single integer >= 1.")
  }
  if (!base::is.numeric(resolution) || base::length(resolution) != 1 || !base::is.finite(resolution) || resolution <= 0) {
    stop("`resolution` must be a single positive number.")
  }

  target_slots <- .hc_longitudinal_step2_target_slots(hc = hc, slot_name = slot_name)
  requested_slot_family <- .hc_longitudinal_requested_slot_family(slot_name)
  if (base::length(target_slots) == 1 && !isTRUE(requested_slot_family)) {
    return(.hc_longitudinal_step2_run_single(
      hc = hc,
      slot_name = target_slots[[1]],
      dimensions = dimensions,
      graph_method = graph_method,
      knn_method = knn_method,
      graph_k = graph_k,
      resolution = resolution,
      leiden_method = leiden_method,
      cluster_prefix = cluster_prefix,
      compute_umap = compute_umap,
      umap_neighbors = umap_neighbors,
      umap_min_dist = umap_min_dist,
      seed = seed,
      file_prefix = "Longitudinal_Meta"
    ))
  }

  plots_by_slot <- list()
  diagnostics_by_slot <- list()
  slot_info <- base::vector("list", base::length(target_slots))

  for (i in base::seq_along(target_slots)) {
    slot_i <- target_slots[[i]]
    slot_ctx <- .hc_longitudinal_step2_slot_context(
      hc = hc,
      slot_name = slot_i,
      requested_slot_name = slot_name
    )
    output_name <- slot_ctx$output_name
    if (!base::nzchar(output_name) || output_name %in% base::names(plots_by_slot)) {
      output_name <- slot_i
    }

    res_i <- .hc_longitudinal_step2_run_single(
      hc = hc,
      slot_name = slot_i,
      dimensions = dimensions,
      graph_method = graph_method,
      knn_method = knn_method,
      graph_k = graph_k,
      resolution = resolution,
      leiden_method = leiden_method,
      cluster_prefix = cluster_prefix,
      compute_umap = compute_umap,
      umap_neighbors = umap_neighbors,
      umap_min_dist = umap_min_dist,
      seed = seed,
      file_prefix = base::paste0("Longitudinal_Meta_", slot_ctx$file_suffix)
    )
    hc <- res_i$hc

    plots_by_slot[[output_name]] <- res_i$plots
    diagnostics_by_slot[[output_name]] <- res_i$diagnostics
    slot_info[[i]] <- base::data.frame(
      output_name = output_name,
      slot_name = slot_i,
      layer_id = slot_ctx$layer_id,
      layer_name = slot_ctx$layer_name,
      source_slot = slot_ctx$source_slot,
      stringsAsFactors = FALSE
    )
  }

  list(
    hc = hc,
    plots = plots_by_slot,
    diagnostics = diagnostics_by_slot,
    slot_info = base::do.call(base::rbind, slot_info)
  )
}

#' Longitudinal step 3: module trajectories by meta-cluster
#'
#' Plots module trajectories after meta-clustering using donor and meta-cluster
#' mean waves.
#'
#' @param hc A `HCoCenaExperiment`.
#' @param slot_name Satellite slot containing meta-clustering outputs. Use a
#'   concrete slot name, a shared prefix such as `"longitudinal_endotypes"` to
#'   process matching suffixed slots, or `"all"` to run over every compatible
#'   step 2 slot in `hc@satellite`.
#' @param facet_ncol Number of facet columns.
#' @param free_y Logical. If `TRUE`, use free y-scales per module.
#' @param square_panels Logical. If `TRUE`, draw square panels.
#' @param value_mode One of `"scaled_mean_vst"` or `"expression"`.
#' @param value_range Optional numeric vector of length 2 for the y-axis range
#'   when `value_mode = "scaled_mean_vst"`. Default is `c(-2, 2)`. Use `NULL`
#'   for automatic scaling.
#'
#' @return A list with unchanged `hc` and the step 3 plot. When multiple slots
#'   are processed, `plots` is returned as a named list keyed by layer/slot,
#'   plus a `slot_info` table describing the resolved slots.
#' @examples
#' hc <- hc_example_data("clustered")
#' hc <- hc_longitudinal_step1_module_donor(
#'   hc,
#'   donor_col = "donor",
#'   time_col = "timepoint",
#'   time_levels = c("T1", "T2"),
#'   k = 2,
#'   nstart = 1,
#'   cap_runs = 1,
#'   impute = FALSE,
#'   seed = 1
#' )$hc
#' hc <- hc_longitudinal_endotype_clustering(hc, k = 2, nstart = 1, seed = 1)
#' hc <- hc_longitudinal_meta_clustering(
#'   hc,
#'   feature_source = "feature_matrix_used",
#'   k = 2:3,
#'   method = "kmeans",
#'   consensus = TRUE,
#'   seed = 1
#' )
#' res <- hc_longitudinal_step3_meta_module_trajectories(hc)
#' @export
hc_longitudinal_step3_meta_module_trajectories <- function(hc,
                                                           slot_name = "longitudinal_endotypes",
                                                           facet_ncol = 4,
                                                           free_y = TRUE,
                                                           square_panels = TRUE,
                                                           value_mode = c("scaled_mean_vst", "expression"),
                                                           value_range = c(-2, 2)) {
  hc <- .hc_longitudinal_unwrap_hc(hc)
  value_mode <- base::match.arg(value_mode)
  target_slots <- .hc_longitudinal_step3_target_slots(hc = hc, slot_name = slot_name)
  requested_slot_family <- .hc_longitudinal_requested_slot_family(slot_name)
  if (base::length(target_slots) == 1 && !isTRUE(requested_slot_family)) {
    return(.hc_longitudinal_step3_run_single(
      hc = hc,
      slot_name = target_slots[[1]],
      facet_ncol = facet_ncol,
      free_y = free_y,
      square_panels = square_panels,
      value_mode = value_mode,
      value_range = value_range,
      file_prefix = "Longitudinal_Meta_ModuleWaves"
    ))
  }

  plots_by_slot <- list()
  slot_info <- base::vector("list", base::length(target_slots))

  for (i in base::seq_along(target_slots)) {
    slot_i <- target_slots[[i]]
    slot_ctx <- .hc_longitudinal_step2_slot_context(
      hc = hc,
      slot_name = slot_i,
      requested_slot_name = slot_name
    )
    output_name <- slot_ctx$output_name
    if (!base::nzchar(output_name) || output_name %in% base::names(plots_by_slot)) {
      output_name <- slot_i
    }

    res_i <- .hc_longitudinal_step3_run_single(
      hc = hc,
      slot_name = slot_i,
      facet_ncol = facet_ncol,
      free_y = free_y,
      square_panels = square_panels,
      value_mode = value_mode,
      value_range = value_range,
      file_prefix = base::paste0("Longitudinal_Meta_ModuleWaves_", slot_ctx$file_suffix)
    )

    plots_by_slot[[output_name]] <- res_i$plots
    slot_info[[i]] <- base::data.frame(
      output_name = output_name,
      slot_name = slot_i,
      layer_id = slot_ctx$layer_id,
      layer_name = slot_ctx$layer_name,
      source_slot = slot_ctx$source_slot,
      stringsAsFactors = FALSE
    )
  }

  list(
    hc = hc,
    plots = plots_by_slot,
    slot_info = base::do.call(base::rbind, slot_info)
  )
}

#' Print longitudinal endotype outputs
#'
#' Convenience printer for outputs from the stepwise longitudinal functions.
#'
#' @param x Output list from one of the longitudinal step functions, or a list
#'   containing nested `step1_module_donor`, `step2_meta_clustering`, and
#'   `step3_meta_module_trajectories` plot blocks.
#' @param show_tables Logical. If `TRUE`, print optional meta-clustering tables
#'   when available.
#'
#' @return Invisibly returns `x`.
#' @examples
#' hc <- hc_example_data("clustered")
#' long <- hc_longitudinal_step1_module_donor(
#'   hc,
#'   donor_col = "donor",
#'   time_col = "timepoint",
#'   time_levels = c("T1", "T2"),
#'   k = 2,
#'   nstart = 1,
#'   cap_runs = 1,
#'   impute = FALSE,
#'   seed = 1
#' )
#' hc <- long$hc
#' hc_print_longitudinal_endotypes(long)
#' @export
hc_print_longitudinal_endotypes <- function(x, show_tables = TRUE) {
  if (base::is.null(x$plots)) {
    return(base::invisible(x))
  }
  p <- x$plots

  is_step1_plot_block <- function(obj) {
    base::is.list(obj) &&
      base::any(base::c("module_means_waves", "module_cluster_waves", "module_cluster_heatmap", "cap_heatmap") %in% base::names(obj))
  }
  is_step2_plot_block <- function(obj) {
    base::is.list(obj) &&
      base::any(base::c("pca", "umap") %in% base::names(obj))
  }
  is_step3_plot_block <- function(obj) {
    base::is.list(obj) &&
      "meta_module_waves" %in% base::names(obj)
  }

  if (!base::is.null(base::names(p)) &&
    base::length(p) > 0 &&
    base::all(base::vapply(p, is_step1_plot_block, FUN.VALUE = base::logical(1)))) {
    for (nm in base::names(p)) {
      layer_block <- p[[nm]]
      if (!base::is.null(layer_block$module_means_waves)) .hc_display_object(layer_block$module_means_waves)
      if (!base::is.null(layer_block$module_cluster_waves)) .hc_display_object(layer_block$module_cluster_waves)
      if (!base::is.null(layer_block$module_cluster_heatmap)) .hc_display_object(layer_block$module_cluster_heatmap)
      if (!base::is.null(layer_block$cap_heatmap)) .hc_display_object(layer_block$cap_heatmap)
    }
    return(base::invisible(x))
  }
  if (!base::is.null(base::names(p)) &&
    base::length(p) > 0 &&
    base::all(base::vapply(p, is_step2_plot_block, FUN.VALUE = base::logical(1)))) {
    for (nm in base::names(p)) {
      slot_block <- p[[nm]]
      diag_block <- if (!base::is.null(x$diagnostics[[nm]])) x$diagnostics[[nm]] else NULL
      if (!base::is.null(slot_block$pca)) .hc_display_object(slot_block$pca)
      if (!base::is.null(slot_block$umap)) .hc_display_object(slot_block$umap)
      if (!base::is.null(diag_block$cross_tab)) .hc_display_object(diag_block$cross_tab)
      if (isTRUE(show_tables) && !base::is.null(diag_block$tables$Method_Clusters)) {
        .hc_display_object(diag_block$tables$Method_Clusters)
      }
      if (isTRUE(show_tables) && !base::is.null(diag_block$tables$Consensus_Decision)) {
        .hc_display_object(diag_block$tables$Consensus_Decision)
      }
    }
    return(base::invisible(x))
  }
  if (!base::is.null(base::names(p)) &&
    base::length(p) > 0 &&
    base::all(base::vapply(p, is_step3_plot_block, FUN.VALUE = base::logical(1)))) {
    for (nm in base::names(p)) {
      slot_block <- p[[nm]]
      if (!base::is.null(slot_block$meta_module_waves)) .hc_display_object(slot_block$meta_module_waves)
    }
    return(base::invisible(x))
  }

  step1 <- p$step1_module_donor
  step2 <- p$step2_meta_clustering
  step3 <- p$step3_meta_module_trajectories

  if (base::is.null(step1) &&
    base::any(base::c("module_means_waves", "module_cluster_waves", "module_cluster_heatmap", "cap_heatmap") %in% base::names(p))) {
    step1 <- p
  }
  if (base::is.null(step2) &&
    base::any(base::c("pca", "umap") %in% base::names(p))) {
    step2 <- p
  }
  if (base::is.null(step3) &&
    "meta_module_waves" %in% base::names(p)) {
    step3 <- p
  }

  if (!base::is.null(step1$module_means_waves)) .hc_display_object(step1$module_means_waves)
  if (!base::is.null(step1$module_cluster_waves)) .hc_display_object(step1$module_cluster_waves)
  if (!base::is.null(step1$module_cluster_heatmap)) .hc_display_object(step1$module_cluster_heatmap)
  if (!base::is.null(step1$cap_heatmap)) .hc_display_object(step1$cap_heatmap)
  if (!base::is.null(step2$pca)) .hc_display_object(step2$pca)
  if (!base::is.null(step2$umap)) .hc_display_object(step2$umap)
  if (!base::is.null(step3$meta_module_waves)) .hc_display_object(step3$meta_module_waves)

  meta_embeddings <- p$meta_embeddings
  if (base::is.null(meta_embeddings) && !base::is.null(x$diagnostics)) {
    meta_embeddings <- list(
      cross_tab = x$diagnostics$cross_tab,
      tables = x$diagnostics$tables
    )
  }
  k_criterion <- if (!base::is.null(p$k_criterion)) p$k_criterion else x$diagnostics$k_criterion
  meta_consensus <- if (!base::is.null(p$meta_consensus)) p$meta_consensus else x$diagnostics$meta_consensus

  if (!base::is.null(meta_embeddings$cross_tab)) .hc_display_object(meta_embeddings$cross_tab)
  if (isTRUE(show_tables) && !base::is.null(k_criterion$global_k)) .hc_display_object(k_criterion$global_k)
  if (isTRUE(show_tables) && !base::is.null(k_criterion$module_k)) .hc_display_object(k_criterion$module_k)
  if (isTRUE(show_tables) && !base::is.null(meta_consensus$consensus_k)) .hc_display_object(meta_consensus$consensus_k)
  if (isTRUE(show_tables) && !base::is.null(meta_consensus$stability)) .hc_display_object(meta_consensus$stability)
  if (isTRUE(show_tables) && !base::is.null(meta_consensus$consensus_heatmap)) .hc_display_object(meta_consensus$consensus_heatmap)
  if (isTRUE(show_tables) && !base::is.null(meta_embeddings$tables$Method_Clusters)) {
    .hc_display_object(meta_embeddings$tables$Method_Clusters)
  }
  if (isTRUE(show_tables) && !base::is.null(meta_embeddings$tables$Consensus_Decision)) {
    .hc_display_object(meta_embeddings$tables$Consensus_Decision)
  }

  base::invisible(x)
}

#' Add meta-cluster x time grouping column for heatmap regrouping
#'
#' Creates an annotation column that combines donor meta-cluster assignment and
#' timepoint labels (e.g. `MC1__24h`). This column can be used as `grouping_v`
#' in `hc_run_expression_analysis_2()` to build a cluster heatmap by
#' meta-cluster x timepoint.
#'
#' @param hc A `HCoCenaExperiment`.
#' @param donor_col Annotation column containing donor IDs.
#' @param time_col Annotation column containing timepoint labels.
#' @param grouping_col Name of the new annotation column to create.
#' @param slot_name Satellite slot holding longitudinal outputs. Must contain
#'   `meta_cluster` table. Use a concrete slot name, a shared prefix such as
#'   `"longitudinal_endotypes"` to process matching suffixed slots, or
#'   `"all"` to use every compatible step 2 slot in `hc@satellite`.
#' @param layer Layer index, layer id/name, or `"all"` (default) to add the
#'   grouping column to all layers.
#' @param time_levels Optional explicit timepoint order vector.
#' @param meta_cluster_levels Optional explicit meta-cluster order vector.
#' @param separator String between meta-cluster and timepoint labels.
#' @param append_layer_suffix Logical. If `TRUE`, append layer id to each group
#'   label (recommended when multiple layers share group names).
#'
#' @return Updated `HCoCenaExperiment`. The computed order is stored in
#'   `hc@satellite[[slot_name]]$meta_time_grouping` and/or the resolved
#'   per-layer step 2 slots.
#' @examples
#' hc <- hc_example_data("clustered")
#' hc <- hc_longitudinal_step1_module_donor(
#'   hc,
#'   donor_col = "donor",
#'   time_col = "timepoint",
#'   time_levels = c("T1", "T2"),
#'   k = 2,
#'   nstart = 1,
#'   cap_runs = 1,
#'   impute = FALSE,
#'   seed = 1
#' )$hc
#' hc <- hc_longitudinal_endotype_clustering(hc, k = 2, nstart = 1, seed = 1)
#' hc <- hc_longitudinal_meta_clustering(
#'   hc,
#'   feature_source = "feature_matrix_used",
#'   k = 2:3,
#'   method = "kmeans",
#'   consensus = TRUE,
#'   seed = 1
#' )
#' hc <- hc_add_meta_time_grouping(hc, donor_col = "donor", time_col = "timepoint")
#' @export
hc_add_meta_time_grouping <- function(hc,
                                      donor_col = "Subject",
                                      time_col = "Time_token",
                                      grouping_col = "meta_cluster_time",
                                      slot_name = "longitudinal_endotypes",
                                      layer = "all",
                                      time_levels = NULL,
                                      meta_cluster_levels = NULL,
                                      separator = "__",
                                      append_layer_suffix = FALSE) {
  if (!inherits(hc, "HCoCenaExperiment")) {
    stop("`hc` must be a `HCoCenaExperiment`.")
  }
  if (!base::is.character(grouping_col) || base::length(grouping_col) != 1 || !base::nzchar(grouping_col)) {
    stop("`grouping_col` must be a non-empty string.")
  }
  if (!base::is.character(slot_name) || base::length(slot_name) != 1 || !base::nzchar(slot_name)) {
    stop("`slot_name` must be a non-empty string.")
  }
  if (!base::is.character(separator) || base::length(separator) != 1) {
    stop("`separator` must be a single string.")
  }

  sat <- as.list(hc@satellite)
  target_slots <- .hc_longitudinal_step3_target_slots(hc = hc, slot_name = slot_name)

  exps <- MultiAssayExperiment::experiments(hc@mae)
  exp_names <- base::names(exps)
  if (base::length(exp_names) == 0) {
    stop("No layers found in `hc@mae`.")
  }

  target_layers <- NULL
  if (base::is.character(layer) && base::length(layer) == 1 && identical(layer, "all")) {
    target_layers <- exp_names
  } else if (base::is.null(layer)) {
    target_layers <- exp_names
  } else {
    target_layers <- .hc_resolve_layer_id(hc = hc, layer = layer)
  }
  target_layers <- base::as.character(target_layers)
  target_layers <- target_layers[!base::is.na(target_layers) & base::nzchar(target_layers)]
  target_layers <- base::unique(target_layers)

  slot_contexts <- base::lapply(
    target_slots,
    function(sn) {
      .hc_longitudinal_step2_slot_context(
        hc = hc,
        slot_name = sn,
        requested_slot_name = slot_name
      )
    }
  )
  base::names(slot_contexts) <- target_slots

  slot_by_layer <- stats::setNames(base::rep(NA_character_, base::length(target_layers)), target_layers)
  if (base::length(target_slots) == 1) {
    only_ctx <- slot_contexts[[1]]
    if (base::is.na(only_ctx$layer_id) || !base::nzchar(only_ctx$layer_id)) {
      slot_by_layer[] <- target_slots[[1]]
    }
  }
  for (lid in target_layers) {
    if (!base::is.na(slot_by_layer[[lid]]) && base::nzchar(slot_by_layer[[lid]])) {
      next
    }
    matched_slots <- base::names(slot_contexts)[
      base::vapply(
        slot_contexts,
        function(ctx) {
          !base::is.na(ctx$layer_id) && identical(base::as.character(ctx$layer_id), base::as.character(lid))
        },
        FUN.VALUE = base::logical(1)
      )
    ]
    if (base::length(matched_slots) == 0) {
      stop(
        "Could not find a longitudinal meta-clustering slot for layer `", lid, "`.\n",
        "Run `hc_longitudinal_step2_meta_clustering()` for this layer first, or pass a matching `slot_name`."
      )
    }
    slot_by_layer[[lid]] <- matched_slots[[1]]
  }

  col_order_by_layer <- list()
  time_levels_by_layer <- list()
  group_levels_by_layer <- list()
  meta_cluster_levels_by_layer <- list()
  slot_used_by_layer <- list()

  for (lid in target_layers) {
    slot_i <- slot_by_layer[[lid]]
    obj <- sat[[slot_i]]
    if (base::is.null(obj) || !base::is.list(obj) || base::is.null(obj$meta_cluster)) {
      stop("Slot `", slot_i, "` must contain `meta_cluster`. Run `hc_longitudinal_step2_meta_clustering()` first.")
    }

    meta_tbl <- base::as.data.frame(obj$meta_cluster, stringsAsFactors = FALSE)
    if (!base::all(c("donor", "meta_cluster") %in% base::colnames(meta_tbl))) {
      stop("`meta_cluster` table in slot `", slot_i, "` must contain columns `donor` and `meta_cluster`.")
    }
    meta_tbl$donor <- base::trimws(base::as.character(meta_tbl$donor))
    meta_tbl$meta_cluster <- base::trimws(base::as.character(meta_tbl$meta_cluster))
    keep_meta <- !base::is.na(meta_tbl$donor) & base::nzchar(meta_tbl$donor) &
      !base::is.na(meta_tbl$meta_cluster) & base::nzchar(meta_tbl$meta_cluster)
    meta_tbl <- meta_tbl[keep_meta, , drop = FALSE]
    meta_tbl <- meta_tbl[!base::duplicated(meta_tbl$donor), , drop = FALSE]
    if (base::nrow(meta_tbl) == 0) {
      stop("`meta_cluster` table in slot `", slot_i, "` does not contain valid donor assignments.")
    }

    if (base::is.null(meta_cluster_levels)) {
      meta_cluster_levels_use <- base::sort(base::unique(meta_tbl$meta_cluster))
    } else {
      meta_cluster_levels_use <- base::trimws(base::as.character(meta_cluster_levels))
      meta_cluster_levels_use <- meta_cluster_levels_use[!base::is.na(meta_cluster_levels_use) & base::nzchar(meta_cluster_levels_use)]
      obs_meta <- base::sort(base::unique(meta_tbl$meta_cluster))
      miss_meta <- obs_meta[!obs_meta %in% meta_cluster_levels_use]
      if (base::length(miss_meta) > 0) {
        warning(
          "`meta_cluster_levels` missed observed clusters for layer `", lid, "`: ",
          base::paste(miss_meta, collapse = ", "), ". Appending at the end."
        )
        meta_cluster_levels_use <- base::c(meta_cluster_levels_use, miss_meta)
      }
    }
    meta_cluster_levels_use <- base::unique(meta_cluster_levels_use)

    se <- exps[[lid]]
    anno <- base::as.data.frame(SummarizedExperiment::colData(se), stringsAsFactors = FALSE)
    if (!base::all(c(donor_col, time_col) %in% base::colnames(anno))) {
      stop("Layer `", lid, "` is missing required annotation columns: `", donor_col, "` and/or `", time_col, "`.")
    }

    donor_vec <- base::trimws(base::as.character(anno[[donor_col]]))
    time_vec <- base::trimws(base::as.character(anno[[time_col]]))
    donor_vec[!base::nzchar(donor_vec)] <- NA_character_
    time_vec[!base::nzchar(time_vec)] <- NA_character_

    if (base::is.null(time_levels)) {
      raw_time <- SummarizedExperiment::colData(se)[[time_col]]
      if (base::is.factor(raw_time)) {
        tl <- base::trimws(base::levels(raw_time))
        tl <- tl[!base::is.na(tl) & base::nzchar(tl)]
        obs_t <- base::unique(time_vec[!base::is.na(time_vec) & base::nzchar(time_vec)])
        time_levels_use <- tl[tl %in% obs_t]
        if (base::length(time_levels_use) == 0) {
          time_levels_use <- base::sort(obs_t)
        }
      } else {
        time_levels_use <- base::sort(base::unique(time_vec[!base::is.na(time_vec) & base::nzchar(time_vec)]))
      }
    } else {
      time_levels_use <- base::trimws(base::as.character(time_levels))
      time_levels_use <- time_levels_use[!base::is.na(time_levels_use) & base::nzchar(time_levels_use)]
      obs_t <- base::unique(time_vec[!base::is.na(time_vec) & base::nzchar(time_vec)])
      miss_t <- obs_t[!obs_t %in% time_levels_use]
      if (base::length(miss_t) > 0) {
        warning(
          "Layer `", lid, "` has time values missing in `time_levels`: ",
          base::paste(miss_t, collapse = ", "), ". Appending at the end."
        )
        time_levels_use <- base::c(time_levels_use, miss_t)
      }
    }
    time_levels_use <- base::unique(time_levels_use)
    if (base::length(time_levels_use) == 0) {
      stop("No valid time levels found for layer `", lid, "`.")
    }

    meta_vec <- meta_tbl$meta_cluster[base::match(donor_vec, meta_tbl$donor)]
    time_idx <- base::match(time_vec, time_levels_use)
    time_use <- base::ifelse(base::is.na(time_idx), NA_character_, time_levels_use[time_idx])

    valid_pair <- !base::is.na(meta_vec) & !base::is.na(time_use) & base::nzchar(time_use)
    combined <- base::rep(NA_character_, base::length(time_use))
    combined[valid_pair] <- base::paste0(meta_vec[valid_pair], separator, time_use[valid_pair])

    combo_levels <- base::as.vector(base::outer(
      meta_cluster_levels_use,
      time_levels_use,
      FUN = function(mc, tm) base::paste0(mc, separator, tm)
    ))
    layer_label <- .hc_longitudinal_layer_label(hc = hc, layer_id = lid)

    if (isTRUE(append_layer_suffix)) {
      combined[!base::is.na(combined)] <- base::paste0(combined[!base::is.na(combined)], "_", layer_label)
      combo_levels <- base::paste0(combo_levels, "_", layer_label)
    }

    anno[[grouping_col]] <- base::factor(combined, levels = combo_levels, ordered = TRUE)
    n_assigned <- base::sum(!base::is.na(anno[[grouping_col]]))
    if (n_assigned == 0) {
      donor_examples <- base::unique(donor_vec[!base::is.na(donor_vec)])
      donor_examples <- utils::head(donor_examples, 8)
      stop(
        "No sample in layer `", lid, "` could be assigned to `", grouping_col, "`.\n",
        "Likely mismatch between annotation donors/times and longitudinal meta-cluster/time levels.\n",
        "Check donor values in `", donor_col, "` and time values in `", time_col, "`.\n",
        "Example donors in this layer: ", base::paste(donor_examples, collapse = ", ")
      )
    }
    n_unassigned <- base::length(anno[[grouping_col]]) - n_assigned
    if (n_unassigned > 0) {
      warning(
        "Layer `", lid, "`: ", n_unassigned, " sample(s) could not be assigned to `",
        grouping_col, "` and remain NA."
      )
    }

    SummarizedExperiment::colData(se) <- S4Vectors::DataFrame(anno, check.names = FALSE)
    exps[[lid]] <- se

    present_levels <- combo_levels[combo_levels %in% base::unique(base::as.character(stats::na.omit(anno[[grouping_col]])))]
    group_levels_by_layer[[lid]] <- present_levels
    col_order_by_layer[[lid]] <- base::paste0(present_levels, "_", layer_label)
    time_levels_by_layer[[lid]] <- time_levels_use
    meta_cluster_levels_by_layer[[lid]] <- meta_cluster_levels_use
    slot_used_by_layer[[lid]] <- slot_i
  }

  hc@mae <- MultiAssayExperiment::MultiAssayExperiment(experiments = exps)

  used_slots <- base::unique(base::as.character(base::unlist(slot_used_by_layer, use.names = FALSE)))
  store_info_for_slot <- function(sat, slot_key, layer_ids) {
    if (base::is.null(sat[[slot_key]]) || !base::is.list(sat[[slot_key]])) {
      sat[[slot_key]] <- list()
    }
    sat[[slot_key]][["meta_time_grouping"]] <- list(
      grouping_col = grouping_col,
      donor_col = donor_col,
      time_col = time_col,
      separator = separator,
      append_layer_suffix = isTRUE(append_layer_suffix),
      layer_ids = layer_ids,
      layer_labels = stats::setNames(
        base::vapply(layer_ids, function(x) .hc_longitudinal_layer_label(hc = hc, layer_id = x), FUN.VALUE = base::character(1)),
        layer_ids
      ),
      slot_name = slot_key,
      requested_slot_name = slot_name,
      slot_by_layer = slot_used_by_layer[layer_ids],
      meta_cluster_levels_by_layer = meta_cluster_levels_by_layer[layer_ids],
      time_levels_by_layer = time_levels_by_layer[layer_ids],
      group_levels_by_layer = group_levels_by_layer[layer_ids],
      col_order_by_layer = col_order_by_layer[layer_ids],
      col_order = base::unlist(col_order_by_layer[layer_ids], use.names = FALSE)
    )
    sat
  }

  for (slot_i in used_slots) {
    slot_layers <- base::names(slot_used_by_layer)[base::unlist(slot_used_by_layer, use.names = FALSE) == slot_i]
    sat <- store_info_for_slot(sat, slot_i, slot_layers)
  }
  if (!identical(slot_name, "all") && base::length(used_slots) > 1) {
    sat <- store_info_for_slot(sat, slot_name, target_layers)
  }
  hc@satellite <- S4Vectors::SimpleList(sat)

  methods::validObject(hc)
  hc
}

#' Get column order for meta-cluster x timepoint heatmaps
#'
#' Returns the column-order vector stored by `hc_add_meta_time_grouping()`.
#'
#' @param hc A `HCoCenaExperiment`.
#' @param slot_name Satellite slot containing `meta_time_grouping`, or a base
#'   prefix such as `"longitudinal_endotypes"` that resolves matching suffixed
#'   slots.
#' @param layer Layer index, layer id/name, or `NULL` for combined order across
#'   all resolved layers.
#'
#' @return Character vector with heatmap column order as expected by
#'   `hc_change_grouping_parameter()`.
#' @examples
#' hc <- hc_example_data("clustered")
#' hc <- hc_longitudinal_step1_module_donor(
#'   hc,
#'   donor_col = "donor",
#'   time_col = "timepoint",
#'   time_levels = c("T1", "T2"),
#'   k = 2,
#'   nstart = 1,
#'   cap_runs = 1,
#'   impute = FALSE,
#'   seed = 1
#' )$hc
#' hc <- hc_longitudinal_endotype_clustering(hc, k = 2, nstart = 1, seed = 1)
#' hc <- hc_longitudinal_meta_clustering(
#'   hc,
#'   feature_source = "feature_matrix_used",
#'   k = 2:3,
#'   method = "kmeans",
#'   consensus = TRUE,
#'   seed = 1
#' )
#' hc <- hc_add_meta_time_grouping(hc, donor_col = "donor", time_col = "timepoint")
#' hc_get_meta_time_col_order(hc)
#' @export
hc_get_meta_time_col_order <- function(hc,
                                       slot_name = "longitudinal_endotypes",
                                       layer = NULL) {
  if (!inherits(hc, "HCoCenaExperiment")) {
    stop("`hc` must be a `HCoCenaExperiment`.")
  }
  sat <- as.list(hc@satellite)

  info_from_obj <- function(obj) {
    if (base::is.null(obj) || !base::is.list(obj) || base::is.null(obj$meta_time_grouping)) {
      return(NULL)
    }
    obj$meta_time_grouping
  }

  target_slots <- .hc_longitudinal_step3_target_slots(hc = hc, slot_name = slot_name)
  info <- info_from_obj(sat[[slot_name]])
  if (base::is.null(info) && base::length(target_slots) == 1) {
    info <- info_from_obj(sat[[target_slots[[1]]]])
  }
  if (base::is.null(info) && base::length(target_slots) > 1) {
    info <- list(
      col_order_by_layer = list(),
      col_order = character(0)
    )
    for (slot_i in target_slots) {
      slot_info <- info_from_obj(sat[[slot_i]])
      if (base::is.null(slot_info)) {
        next
      }
      if (!base::is.null(slot_info$col_order_by_layer)) {
        info$col_order_by_layer <- c(info$col_order_by_layer, slot_info$col_order_by_layer)
      }
      slot_col_order <- base::as.character(slot_info$col_order)
      slot_col_order <- slot_col_order[!base::is.na(slot_col_order) & base::nzchar(slot_col_order)]
      info$col_order <- base::c(info$col_order, slot_col_order)
    }
  }
  if (base::is.null(info)) {
    stop("No `meta_time_grouping` found for `", slot_name, "`. Run `hc_add_meta_time_grouping()` first.")
  }

  if (base::is.null(layer)) {
    out <- base::as.character(info$col_order)
    out <- out[!base::is.na(out) & base::nzchar(out)]
    return(base::unique(out))
  }

  lid <- .hc_resolve_layer_id(hc = hc, layer = layer)
  if (base::is.null(info$col_order_by_layer[[lid]])) {
    stop("No stored col-order for layer `", lid, "`.")
  }
  out <- base::as.character(info$col_order_by_layer[[lid]])
  out <- out[!base::is.na(out) & base::nzchar(out)]
  base::unique(out)
}
