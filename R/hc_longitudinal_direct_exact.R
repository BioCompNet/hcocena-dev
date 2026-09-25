#' Prepare one module matrix for direct trajectory clustering
#' @noRd
.hc_direct_prepare_module_input <- function(time_data,
                                            donor_col = "donor",
                                            impute = TRUE,
                                            impute_method = "rfcont",
                                            ntree = 10,
                                            na_impute = c("median", "zero"),
                                            scale_features = TRUE,
                                            seed = 42) {
  na_impute <- base::match.arg(na_impute)

  donor_ids_all <- base::as.character(time_data[[donor_col]])
  time_data_use <- base::as.data.frame(time_data, stringsAsFactors = FALSE)

  if (isTRUE(impute)) {
    time_data_use <- .hc_longitudinal_impute_time_data(
      time_data = time_data_use,
      donor_col = donor_col,
      method = impute_method,
      ntree = ntree,
      seed = seed
    )
  }

  time_data_use <- time_data_use[base::match(donor_ids_all, base::as.character(time_data_use[[donor_col]])), , drop = FALSE]
  mat_raw <- base::as.matrix(time_data_use[, base::setdiff(base::colnames(time_data_use), donor_col), drop = FALSE])
  storage.mode(mat_raw) <- "numeric"
  base::rownames(mat_raw) <- donor_ids_all

  if (identical(na_impute, "median")) {
    for (cc in base::seq_len(base::ncol(mat_raw))) {
      v <- mat_raw[, cc]
      if (base::any(base::is.na(v))) {
        med <- stats::median(v, na.rm = TRUE)
        if (!base::is.finite(med)) {
          med <- 0
        }
        v[base::is.na(v)] <- med
        mat_raw[, cc] <- v
      }
    }
  } else {
    mat_raw[base::is.na(mat_raw)] <- 0
  }

  keep_cols <- base::apply(mat_raw, 2, stats::sd, na.rm = TRUE) > 0
  mat_use <- mat_raw[, keep_cols, drop = FALSE]
  if (base::ncol(mat_use) == 0) {
    mat_use <- mat_raw
  }

  if (isTRUE(scale_features)) {
    mat_use <- base::as.matrix(base::scale(mat_use))
    mat_use[!base::is.finite(mat_use)] <- 0
  }

  list(
    donors_all = donor_ids_all,
    donors_used = donor_ids_all,
    donors_dropped = character(0),
    time_data_complete = time_data_use,
    feature_matrix_raw = mat_raw,
    feature_matrix = mat_use
  )
}

#' Run one fixed-k direct clustering
#' @noRd
.hc_direct_cluster_fixed_k <- function(x,
                                       k,
                                       method = c("kmeans", "hclust"),
                                       nstart = 50,
                                       seed = 42) {
  method <- base::match.arg(method)
  x <- base::as.matrix(x)

  if (k <= 1L || base::nrow(x) <= 1L) {
    return(list(
      labels = base::rep(1L, base::nrow(x)),
      tot_withinss = NA_real_,
      fit = NULL
    ))
  }

  if (identical(method, "kmeans")) {
    fit <- .hc_with_seed(seed, tryCatch(
      stats::kmeans(x, centers = base::as.integer(k), nstart = base::as.integer(nstart), iter.max = 100),
      error = function(e) NULL
    ))
    if (!base::is.null(fit) && !base::is.null(fit$cluster)) {
      return(list(
        labels = base::as.integer(fit$cluster),
        tot_withinss = .hc_first_numeric_value(fit$tot.withinss),
        fit = fit
      ))
    }
  }

  h_tree <- stats::hclust(stats::dist(x), method = "ward.D2")
  list(
    labels = base::as.integer(stats::cutree(h_tree, k = base::as.integer(k))),
    tot_withinss = NA_real_,
    fit = h_tree
  )
}

#' Select best k for one module using direct trajectory clustering
#' @noRd
.hc_direct_select_k <- function(feature_matrix,
                                k = 2:6,
                                method = c("kmeans", "hclust"),
                                nstart = 50,
                                min_cluster_fraction = 0.1,
                                score_method = c("calinski_harabasz", "tot_withinss"),
                                seed = 42) {
  method <- base::match.arg(method)
  score_method <- base::match.arg(score_method)
  feature_matrix <- base::as.matrix(feature_matrix)

  donors_n <- base::nrow(feature_matrix)
  valid_k <- base::sort(base::unique(base::as.integer(k)))
  valid_k <- valid_k[base::is.finite(valid_k) & !base::is.na(valid_k)]
  valid_k <- valid_k[valid_k >= 2L & valid_k < donors_n]

  zero_var <- if (base::ncol(feature_matrix) > 0) {
    base::all(base::apply(feature_matrix, 2, function(v) stats::sd(v, na.rm = TRUE) == 0))
  } else {
    TRUE
  }

  if (donors_n <= 2L || base::ncol(feature_matrix) == 0L || zero_var || base::length(valid_k) == 0L) {
    return(list(
      best_k = 1L,
      best_labels = if (donors_n > 0) base::rep(1L, donors_n) else integer(0),
      score_table = base::data.frame(
        k = 1L,
        ch_index = NA_real_,
        tot_withinss = NA_real_,
        cluster_min_fraction = if (donors_n > 0) 1 else NA_real_,
        eligible = TRUE,
        selection_pool = "singleton",
        stringsAsFactors = FALSE
      )
    ))
  }

  score_rows <- list()
  labels_by_k <- list()

  for (kk in valid_k) {
    fit_kk <- .hc_direct_cluster_fixed_k(
      x = feature_matrix,
      k = kk,
      method = method,
      nstart = nstart,
      seed = seed + kk
    )
    lbl_kk <- base::as.integer(fit_kk$labels)
    labels_by_k[[base::as.character(kk)]] <- lbl_kk

    min_fraction_kk <- if (base::length(lbl_kk) > 0) {
      base::min(base::table(lbl_kk)) / base::length(lbl_kk)
    } else {
      NA_real_
    }
    eligible_kk <- base::is.null(min_cluster_fraction) ||
      !base::is.finite(min_cluster_fraction) ||
      min_fraction_kk >= min_cluster_fraction

    ch_kk <- if (kk < base::length(lbl_kk)) {
      .hc_ch_index(feature_matrix, lbl_kk)
    } else {
      NA_real_
    }

    score_rows[[base::length(score_rows) + 1L]] <- base::data.frame(
      k = base::as.integer(kk),
      ch_index = ch_kk,
      tot_withinss = fit_kk$tot_withinss,
      cluster_min_fraction = min_fraction_kk,
      eligible = isTRUE(eligible_kk),
      stringsAsFactors = FALSE
    )
  }

  score_tbl <- do.call(base::rbind, score_rows)
  if (identical(score_method, "tot_withinss") &&
    identical(method, "kmeans") &&
    base::any(base::is.finite(score_tbl$tot_withinss))) {
    score_vec <- -score_tbl$tot_withinss
  } else {
    score_vec <- score_tbl$ch_index
  }

  eligible_finite_rows <- base::which(score_tbl$eligible & base::is.finite(score_vec))
  any_finite_rows <- base::which(base::is.finite(score_vec))

  if (base::length(eligible_finite_rows) > 0L) {
    score_tbl$selection_pool <- "eligible"
    choose_rows <- eligible_finite_rows
  } else if (base::length(any_finite_rows) > 0L) {
    score_tbl$selection_pool <- "fallback_all"
    choose_rows <- any_finite_rows
  } else {
    score_tbl$selection_pool <- "fallback_smallest_k"
    choose_rows <- integer(0)
  }

  if (base::length(choose_rows) > 0L) {
    sub_tbl <- score_tbl[choose_rows, , drop = FALSE]
    sub_score <- score_vec[choose_rows]
    sub_tbl$.score <- base::replace(sub_score, !base::is.finite(sub_score), -Inf)
    sub_tbl <- sub_tbl[base::order(-sub_tbl$.score, sub_tbl$k), , drop = FALSE]
    best_k <- base::as.integer(sub_tbl$k[[1]])
  } else {
    best_k <- base::min(valid_k)
  }

  list(
    best_k = base::as.integer(best_k),
    best_labels = base::as.integer(labels_by_k[[base::as.character(best_k)]]),
    score_table = score_tbl
  )
}

#' Run longitudinal step 1 using direct trajectory clustering
#' @noRd
.hc_run_direct_step1_exact <- function(hc,
                                       means_slot = "longitudinal_module_means",
                                       output_slot = "longitudinal_endotypes",
                                       k = 2:6,
                                       method = c("kmeans", "hclust"),
                                       nstart = 50,
                                       cap_runs = 50,
                                       impute = TRUE,
                                       impute_method = "rfcont",
                                       ntree = 10,
                                       min_cluster_fraction = 0.1,
                                       score_method = c("calinski_harabasz", "tot_withinss"),
                                       scale_features = TRUE,
                                       na_impute = c("median", "zero"),
                                       cap_scale_features = scale_features,
                                       cap_na_impute = na_impute,
                                       cap_jitter_sd = 0,
                                       seed = 42) {
  method <- base::match.arg(method)
  score_method <- base::match.arg(score_method)
  na_impute <- base::match.arg(na_impute)
  # Explicit choices: the formal default is `na_impute`, which has already
  # been collapsed to a single value by the line above. Without them an
  # explicitly supplied vector - as passed down by
  # hc_longitudinal_workflow_direct() - is matched against one choice only
  # and match.arg() aborts with "'arg' must be of length 1".
  cap_na_impute <- base::match.arg(cap_na_impute, choices = c("median", "zero"))

  sat <- as.list(hc@satellite)
  inp <- sat[[means_slot]]
  req <- c("module_lookup", "donor_time_module", "time_levels")
  if (is.null(inp) || !base::all(req %in% base::names(inp))) {
    stop("Step 1 requires outputs from `hc_longitudinal_module_means()`.", call. = FALSE)
  }

  donor_time_module <- base::as.data.frame(inp$donor_time_module, stringsAsFactors = FALSE)
  if (base::nrow(donor_time_module) == 0) {
    stop("No `donor_time_module` data found in slot `", means_slot, "`.", call. = FALSE)
  }

  module_lookup <- base::as.data.frame(inp$module_lookup, stringsAsFactors = FALSE)
  module_lookup <- .hc_reorder_module_lookup_natural(module_lookup)
  if (!"module_color" %in% base::colnames(module_lookup)) {
    module_lookup$module_color <- "grey60"
  }

  module_ids <- base::as.character(module_lookup$module)
  donors <- if (!base::is.null(inp$donor_feature_matrix)) {
    base::rownames(base::as.matrix(inp$donor_feature_matrix))
  } else {
    base::sort(base::unique(base::as.character(donor_time_module$donor)))
  }
  donors <- base::as.character(donors)
  time_levels <- base::as.character(inp$time_levels)

  module_score_rows <- list()
  module_best_rows <- list()
  module_assignment_rows <- list()
  filter_rows <- list()

  for (m in module_ids) {
    time_data <- reshape2::dcast(
      donor_time_module[donor_time_module$module == m, , drop = FALSE],
      donor ~ time,
      value.var = "value"
    )
    time_data <- base::merge(
      base::data.frame(donor = donors, stringsAsFactors = FALSE),
      time_data,
      by = "donor",
      all.x = TRUE,
      sort = FALSE
    )
    time_data <- time_data[base::match(donors, base::as.character(time_data$donor)), , drop = FALSE]
    base::colnames(time_data) <- base::as.character(base::colnames(time_data))

    keep_cols <- c("donor", time_levels)
    missing_tp <- base::setdiff(keep_cols, base::colnames(time_data))
    for (tp in missing_tp) {
      time_data[[tp]] <- NA_real_
    }
    time_data <- time_data[, keep_cols, drop = FALSE]

    prep <- .hc_direct_prepare_module_input(
      time_data = time_data,
      donor_col = "donor",
      impute = impute,
      impute_method = impute_method,
      ntree = ntree,
      na_impute = na_impute,
      scale_features = scale_features,
      seed = seed
    )

    select_kk <- .hc_direct_select_k(
      feature_matrix = prep$feature_matrix,
      k = k,
      method = method,
      nstart = nstart,
      min_cluster_fraction = min_cluster_fraction,
      score_method = score_method,
      seed = seed
    )

    donors_used <- base::as.character(prep$donors_used)
    best_labels <- base::as.integer(select_kk$best_labels)
    if (base::length(donors_used) > 0 && base::length(best_labels) > 0) {
      module_assignment_rows[[base::length(module_assignment_rows) + 1L]] <- base::data.frame(
        donor = donors_used,
        module = m,
        cluster = best_labels,
        stringsAsFactors = FALSE
      )
    }

    score_tbl_module <- base::as.data.frame(select_kk$score_table, stringsAsFactors = FALSE)
    if (base::nrow(score_tbl_module) > 0) {
      score_tbl_module$module <- m
    }
    module_score_rows[[base::length(module_score_rows) + 1L]] <- score_tbl_module
    module_best_rows[[base::length(module_best_rows) + 1L]] <- base::data.frame(
      module = m,
      best_k = base::as.integer(select_kk$best_k),
      selected_run = 1L,
      stringsAsFactors = FALSE
    )
    filter_rows[[base::length(filter_rows) + 1L]] <- base::data.frame(
      module = m,
      donors_total = base::length(prep$donors_all),
      donors_used = base::length(prep$donors_used),
      donors_dropped = base::length(prep$donors_dropped),
      stringsAsFactors = FALSE
    )
  }

  module_cluster_assignments <- if (base::length(module_assignment_rows) > 0) {
    base::do.call(base::rbind, module_assignment_rows)
  } else {
    base::data.frame(
      donor = character(0),
      module = character(0),
      cluster = integer(0),
      stringsAsFactors = FALSE
    )
  }
  if (base::nrow(module_cluster_assignments) > 0) {
    module_cluster_assignments$module_cluster_key <- base::paste0(
      module_cluster_assignments$module,
      "__",
      module_cluster_assignments$cluster
    )
  } else {
    module_cluster_assignments$module_cluster_key <- character(0)
  }

  module_cluster_matrix <- base::matrix(
    NA_integer_,
    nrow = base::length(donors),
    ncol = base::length(module_ids),
    dimnames = list(donors, module_ids)
  )
  if (base::nrow(module_cluster_assignments) > 0) {
    for (i in base::seq_len(base::nrow(module_cluster_assignments))) {
      d <- module_cluster_assignments$donor[[i]]
      m <- module_cluster_assignments$module[[i]]
      module_cluster_matrix[d, m] <- base::as.integer(module_cluster_assignments$cluster[[i]])
    }
  }

  donor_time_module_with_clusters <- base::merge(
    donor_time_module,
    module_cluster_assignments[, c("donor", "module", "cluster", "module_cluster_key"), drop = FALSE],
    by = c("donor", "module"),
    all.x = TRUE,
    sort = FALSE
  )

  module_cluster_trajectory_mean <- stats::aggregate(
    value ~ module + time + cluster + module_cluster_key,
    data = donor_time_module_with_clusters,
    FUN = base::mean
  )

  module_cluster_best_k <- if (base::length(module_best_rows) > 0) {
    base::do.call(base::rbind, module_best_rows)
  } else {
    base::data.frame(module = character(0), best_k = integer(0), selected_run = integer(0), stringsAsFactors = FALSE)
  }
  module_cluster_score_table <- if (base::length(module_score_rows) > 0) {
    base::do.call(base::rbind, module_score_rows)
  } else {
    base::data.frame(
      module = character(0),
      k = integer(0),
      ch_index = numeric(0),
      tot_withinss = numeric(0),
      cluster_min_fraction = numeric(0),
      eligible = logical(0),
      selection_pool = character(0),
      stringsAsFactors = FALSE
    )
  }
  module_cluster_palette <- .hc_module_cluster_palette(
    module_lookup = module_lookup,
    module_best_k = module_cluster_best_k[, c("module", "best_k"), drop = FALSE]
  )

  sat[[output_slot]] <- list(
    created_at = as.character(base::Sys.time()),
    source_slot = means_slot,
    donor_col = inp$donor_col,
    time_col = inp$time_col,
    group_col = inp$group_col,
    time_levels = time_levels,
    module_lookup = module_lookup,
    donor_time_module = donor_time_module,
    value_label = inp$value_label,
    module_cluster_matrix = base::as.matrix(module_cluster_matrix),
    module_cluster_assignments = module_cluster_assignments,
    donor_time_module_with_module_clusters = donor_time_module_with_clusters,
    module_cluster_trajectory_mean = module_cluster_trajectory_mean,
    module_cluster_palette = module_cluster_palette,
    module_cluster_best_k = module_cluster_best_k,
    module_cluster_score = module_cluster_score_table,
    module_cluster_score_table = module_cluster_score_table,
    module_cluster_score_method = score_method,
    module_cluster_min_fraction = min_cluster_fraction,
    direct_method = method,
    direct_candidate_k = base::as.integer(k),
    direct_nstart = base::as.integer(nstart),
    direct_cap_runs = base::as.integer(cap_runs),
    direct_scale_features = isTRUE(scale_features),
    direct_na_impute = na_impute,
    direct_module_filtering = base::do.call(base::rbind, filter_rows),
    master_donors = donors
  )
  hc@satellite <- S4Vectors::SimpleList(sat)
  methods::validObject(hc)

  hc <- hc_longitudinal_module_cap(
    hc = hc,
    slot_name = output_slot,
    output_slot = output_slot,
    runs = cap_runs,
    method = method,
    nstart = nstart,
    seed = seed,
    scale_features = cap_scale_features,
    na_impute = cap_na_impute,
    jitter_sd = cap_jitter_sd
  )

  sat <- as.list(hc@satellite)
  sat[[output_slot]]$direct_cap_scale_features <- isTRUE(cap_scale_features)
  sat[[output_slot]]$direct_cap_na_impute <- cap_na_impute
  sat[[output_slot]]$direct_cap_jitter_sd <- as.numeric(cap_jitter_sd)
  hc@satellite <- S4Vectors::SimpleList(sat)
  methods::validObject(hc)
  hc
}

#' Run one single-layer direct longitudinal step 1
#' @noRd
.hc_longitudinal_step1_run_single_direct <- function(hc,
                                                     donor_col,
                                                     time_col,
                                                     layer_id,
                                                     time_levels,
                                                     k,
                                                     method,
                                                     nstart,
                                                     cap_runs,
                                                     impute,
                                                     impute_method,
                                                     ntree,
                                                     min_cluster_fraction,
                                                     score_method,
                                                     scale_features,
                                                     na_impute,
                                                     cap_scale_features,
                                                     cap_na_impute,
                                                     cap_jitter_sd,
                                                     seed,
                                                     means_slot,
                                                     output_slot,
                                                     module_means_prefix = "Longitudinal_Direct_ModuleMeans",
                                                     module_clusters_prefix = "Longitudinal_Direct_ModuleClusters",
                                                     cap_prefix = "Longitudinal_Direct_CAP") {
  lmm_args <- list(
    hc = hc,
    donor_col = donor_col,
    time_col = time_col,
    layer = layer_id,
    group_col = NULL,
    use_module_labels = TRUE,
    time_levels = time_levels,
    impute_missing = "none",
    slot_name = means_slot
  )
  if ("value_label" %in% base::names(base::formals(hc_longitudinal_module_means))) {
    lmm_args$value_label <- NULL
  }
  hc <- do.call(hc_longitudinal_module_means, lmm_args)

  hc <- .hc_run_direct_step1_exact(
    hc = hc,
    means_slot = means_slot,
    output_slot = output_slot,
    k = k,
    method = method,
    nstart = nstart,
    cap_runs = cap_runs,
    impute = impute,
    impute_method = impute_method,
    ntree = ntree,
    min_cluster_fraction = min_cluster_fraction,
    score_method = score_method,
    scale_features = scale_features,
    na_impute = na_impute,
    cap_scale_features = cap_scale_features,
    cap_na_impute = cap_na_impute,
    cap_jitter_sd = cap_jitter_sd,
    seed = seed
  )

  module_means <- hc_plot_longitudinal_module_means(
    hc,
    slot_name = means_slot,
    save_pdf = TRUE,
    file_prefix = module_means_prefix,
    square_panels = TRUE,
    save_width = 10,
    save_height = 10
  )

  module_clusters <- hc_plot_longitudinal_module_clusters(
    hc,
    slot_name = output_slot,
    save_pdf = TRUE,
    file_prefix = module_clusters_prefix,
    show_heatmap_numbers = TRUE,
    square_panels = TRUE,
    save_waves_width = 10,
    save_waves_height = 10
  )

  donor_order <- NULL
  heatmap_obj <- module_clusters$module_cluster_heatmap
  if (!base::is.null(heatmap_obj) &&
    !base::is.null(heatmap_obj$data) &&
    "donor" %in% base::colnames(heatmap_obj$data)) {
    donor_vec <- heatmap_obj$data$donor
    donor_levels <- base::levels(donor_vec)
    donor_order <- if (!base::is.null(donor_levels) && base::length(donor_levels) > 0) {
      donor_levels
    } else {
      base::unique(base::as.character(donor_vec))
    }
  }

  cap_args <- list(
    hc = hc,
    slot_name = output_slot,
    save_pdf = TRUE,
    file_prefix = cap_prefix,
    show_values = FALSE
  )
  if (!base::is.null(donor_order)) {
    cap_args$donor_order <- donor_order
  }
  cap <- do.call(hc_plot_longitudinal_cap, cap_args)

  list(
    hc = hc,
    plots = list(
      module_means_waves = module_means$module_means,
      module_cluster_waves = module_clusters$module_cluster_waves,
      module_cluster_heatmap = module_clusters$module_cluster_heatmap,
      cap_heatmap = cap$cap_heatmap
    ),
    diagnostics = list(
      k = k,
      method = method,
      nstart = base::as.integer(nstart),
      cap_runs = base::as.integer(cap_runs),
      impute = isTRUE(impute),
      impute_method = impute_method,
      ntree = base::as.integer(ntree),
      min_cluster_fraction = as.numeric(min_cluster_fraction),
      score_method = score_method,
      scale_features = isTRUE(scale_features),
      na_impute = na_impute,
      module_cluster_score = hc@satellite[[output_slot]]$module_cluster_score,
      module_cluster_best_k = hc@satellite[[output_slot]]$module_cluster_best_k,
      direct_module_filtering = hc@satellite[[output_slot]]$direct_module_filtering
    )
  )
}

#' Longitudinal step 1: module/donor clustering with direct trajectories
#'
#' Computes longitudinal module means, then clusters donors separately within
#' each module directly on the donor-by-time trajectory matrix, especially for
#' short trajectories with only a few timepoints.
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
#' @param impute_method Missing-value method passed to `mice`. Use `"rfcont"`
#'   to match the previous workflow most closely. For `"rfcont"`,
#'   `CALIBERrfimpute` must be installed; it is loaded automatically.
#' @param ntree Number of trees for `rfcont` imputation.
#' @param min_cluster_fraction Minimum allowed fraction of donors in the
#'   smallest module cluster.
#' @param score_method One of `"calinski_harabasz"` or `"tot_withinss"` for
#'   selecting the final per-module `k`.
#' @param scale_features Logical. If `TRUE`, z-score timepoint columns before
#'   clustering.
#' @param na_impute One of `"median"` or `"zero"` for final NA handling after
#'   optional imputation.
#' @param cap_scale_features Logical. If `TRUE`, z-score features for the CAP
#'   reruns.
#' @param cap_na_impute One of `"median"` or `"zero"` for CAP feature handling.
#' @param cap_jitter_sd Optional Gaussian jitter SD added per CAP rerun when
#'   `method = "kmeans"`.
#' @param seed Random seed.
#' @param means_slot Satellite slot for module means.
#' @param output_slot Satellite slot for direct endotype outputs.
#'
#' @return A list with updated `hc`, step plots, and diagnostics. When
#'   `layer = "all"` (or multiple layers are supplied), returns per-layer plot
#'   and diagnostic lists keyed by layer id.
#' @examples
#' hc <- hc_example_data("clustered")
#' res <- hc_longitudinal_step1_module_donor_direct(
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
hc_longitudinal_step1_module_donor_direct <- function(hc,
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
                                                      scale_features = TRUE,
                                                      na_impute = c("median", "zero"),
                                                      cap_scale_features = scale_features,
                                                      cap_na_impute = na_impute,
                                                      cap_jitter_sd = 0,
                                                      seed = 42,
                                                      means_slot = "longitudinal_module_means_direct",
                                                      output_slot = "longitudinal_endotypes_direct") {
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

  target_layers <- .hc_longitudinal_step1_target_layers(hc = hc, layer = layer)
  multi_layer <- base::length(target_layers) > 1

  if (!isTRUE(multi_layer)) {
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
      output_slot = output_slot
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
      module_means_prefix = base::paste0("Longitudinal_Direct_ModuleMeans_", file_suffix),
      module_clusters_prefix = base::paste0("Longitudinal_Direct_ModuleClusters_", file_suffix),
      cap_prefix = base::paste0("Longitudinal_Direct_CAP_", file_suffix)
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

#' Convenience longitudinal workflow using the direct step-1 backend
#'
#' Runs the direct donor-by-time step 1 backend and then reuses the existing
#' step 2 and step 3 quick wrappers on the resulting output slot family.
#'
#' @inheritParams hc_longitudinal_step1_module_donor_direct
#'
#' @return A list with updated `hc`, nested `plots`, and per-step diagnostics.
#' @examples
#' hc <- hc_example_data("clustered")
#' res <- hc_longitudinal_workflow_direct(
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
hc_longitudinal_workflow_direct <- function(hc,
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
                                            scale_features = TRUE,
                                            na_impute = c("median", "zero"),
                                            cap_scale_features = scale_features,
                                            cap_na_impute = na_impute,
                                            cap_jitter_sd = 0,
                                            seed = 42,
                                            means_slot = "longitudinal_module_means_direct",
                                            output_slot = "longitudinal_endotypes_direct") {
  step1 <- hc_longitudinal_step1_module_donor_direct(
    hc = hc,
    donor_col = donor_col,
    time_col = time_col,
    layer = layer,
    time_levels = time_levels,
    k = k,
    method = method,
    nstart = nstart,
    cap_runs = cap_runs,
    impute = impute,
    impute_method = impute_method,
    ntree = ntree,
    min_cluster_fraction = min_cluster_fraction,
    score_method = score_method,
    scale_features = scale_features,
    na_impute = na_impute,
    cap_scale_features = cap_scale_features,
    cap_na_impute = cap_na_impute,
    cap_jitter_sd = cap_jitter_sd,
    seed = seed,
    means_slot = means_slot,
    output_slot = output_slot
  )
  hc <- step1$hc

  step2 <- hc_longitudinal_step2_meta_clustering(
    hc = hc,
    slot_name = output_slot,
    seed = seed
  )
  hc <- step2$hc

  step3 <- hc_longitudinal_step3_meta_module_trajectories(
    hc = hc,
    slot_name = output_slot
  )
  hc <- step3$hc

  list(
    hc = hc,
    plots = list(
      step1_module_donor = step1$plots,
      step2_meta_clustering = step2$plots,
      step3_meta_module_trajectories = step3$plots
    ),
    diagnostics = list(
      step1_module_donor = step1$diagnostics,
      step2_meta_clustering = step2$diagnostics,
      step3_meta_module_trajectories = step3$slot_info
    ),
    slot_info = list(
      step1 = step1$layer_info,
      step2 = step2$slot_info,
      step3 = step3$slot_info
    )
  )
}
