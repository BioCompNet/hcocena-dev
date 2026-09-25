#' Tiered and robust cutoff recommendation
#'
#' Alternative cutoff tuning method that can be compared against the current
#' simple cutoff recommendation.
#'
#' For each layer, the function:
#' 1. Summarizes cutoff statistics.
#' 2. Computes additional diagnostics (`slope`, `mean_k`).
#' 3. Applies a staged rule set (`strict`, `relaxed`, `best_available`).
#'
#' Results are stored in `hc@satellite$cutoff_tuning`
#' (and mirrored to `hc@satellite$auto_tune_cutoff_tiered` for compatibility).
#'
#' @param hc A `HCoCenaExperiment` object.
#' @param apply Logical; if `TRUE`, apply the resulting cutoff vector with
#'   `hc_set_cutoff()`. Missing tiered recommendations fall back to existing
#'   cutoffs from `hc@config@layer$cutoff` where possible.
#' @param tier1_sft_min Numeric minimum SFT fit (`R.squared`) for the `strict`
#'   rule.
#' @param tier1_mean_k_min Numeric minimum mean connectivity (`mean_k`) for the
#'   `strict` rule.
#' @param tier2_sft_min Numeric minimum SFT fit (`R.squared`) for the
#'   `relaxed` rule.
#' @param tier2_mean_k_min Numeric minimum mean connectivity (`mean_k`) for the
#'   `relaxed` rule.
#' @param fallback_sft_min Numeric lower bound for `best_available` selection
#'   using the maximum available SFT fit.
#' @param require_negative_slope Logical; require a negative slope (`slope < 0`)
#'   in Tier 1/Tier 2 (and optionally fallback, see
#'   `fallback_respect_base_filters`).
#' @param max_no_networks Optional numeric scalar. If set, only cutoffs with
#'   `no_of_networks <= max_no_networks` are accepted in Tier 1/Tier 2.
#' @param max_components Deprecated alias for `max_no_networks` (kept for
#'   backward compatibility).
#' @param min_node_fraction Optional numeric in `(0, 1]`. If set, require
#'   `no_nodes >= min_node_fraction * max(no_nodes)` in Tier 1/Tier 2.
#' @param fallback_respect_base_filters Logical; if `TRUE`, `best_available`
#'   selection is restricted by base filters (`slope`, `max_no_networks`,
#'   `min_node_fraction`).
#' @param verbose Logical; print per-layer recommendation messages.
#'
#' @return Updated `HCoCenaExperiment` with a report in
#'   `hc@satellite$cutoff_tuning`.
#' @examples
#' hc <- hc_example_data("after_part1")
#' hc <- hc_auto_tune_cutoff_tiered(hc, apply = FALSE, verbose = FALSE)
#' @export
hc_auto_tune_cutoff_tiered <- function(hc,
                                       apply = FALSE,
                                       tier1_sft_min = 0.85,
                                       tier1_mean_k_min = 10,
                                       tier2_sft_min = 0.80,
                                       tier2_mean_k_min = 5,
                                       fallback_sft_min = 0.70,
                                       require_negative_slope = TRUE,
                                       max_no_networks = NULL,
                                       max_components = NULL,
                                       min_node_fraction = NULL,
                                       fallback_respect_base_filters = TRUE,
                                       verbose = TRUE) {
  if (!inherits(hc, "HCoCenaExperiment")) {
    stop("`hc` must be a `HCoCenaExperiment`.")
  }
  if (!is.logical(apply) || length(apply) != 1 || is.na(apply)) {
    stop("`apply` must be TRUE or FALSE.")
  }
  if (!is.numeric(tier1_sft_min) || length(tier1_sft_min) != 1 || !is.finite(tier1_sft_min)) {
    stop("`tier1_sft_min` must be a finite numeric scalar.")
  }
  if (!is.numeric(tier1_mean_k_min) || length(tier1_mean_k_min) != 1 || !is.finite(tier1_mean_k_min)) {
    stop("`tier1_mean_k_min` must be a finite numeric scalar.")
  }
  if (!is.numeric(tier2_sft_min) || length(tier2_sft_min) != 1 || !is.finite(tier2_sft_min)) {
    stop("`tier2_sft_min` must be a finite numeric scalar.")
  }
  if (!is.numeric(tier2_mean_k_min) || length(tier2_mean_k_min) != 1 || !is.finite(tier2_mean_k_min)) {
    stop("`tier2_mean_k_min` must be a finite numeric scalar.")
  }
  if (!is.numeric(fallback_sft_min) || length(fallback_sft_min) != 1 || !is.finite(fallback_sft_min)) {
    stop("`fallback_sft_min` must be a finite numeric scalar.")
  }
  if (!is.logical(require_negative_slope) || length(require_negative_slope) != 1 || is.na(require_negative_slope)) {
    stop("`require_negative_slope` must be TRUE or FALSE.")
  }
  if (!is.null(max_no_networks)) {
    if (!is.numeric(max_no_networks) || length(max_no_networks) != 1 || !is.finite(max_no_networks) || max_no_networks < 1) {
      stop("`max_no_networks` must be NULL or a positive numeric scalar.")
    }
  }
  if (!is.null(max_components)) {
    if (!is.numeric(max_components) || length(max_components) != 1 || !is.finite(max_components) || max_components < 1) {
      stop("`max_components` (deprecated) must be NULL or a positive numeric scalar.")
    }
  }
  if (!is.null(max_no_networks) && !is.null(max_components) &&
    abs(as.numeric(max_no_networks) - as.numeric(max_components)) > .Machine$double.eps^0.5) {
    stop("Please set only one of `max_no_networks` or `max_components` (deprecated alias), or set both to the same value.")
  }
  if (is.null(max_no_networks) && !is.null(max_components)) {
    max_no_networks <- as.numeric(max_components)
    warning("`max_components` is deprecated; please use `max_no_networks`.", call. = FALSE)
  }
  if (!is.null(min_node_fraction)) {
    if (!is.numeric(min_node_fraction) || length(min_node_fraction) != 1 || !is.finite(min_node_fraction) ||
      min_node_fraction <= 0 || min_node_fraction > 1) {
      stop("`min_node_fraction` must be NULL or a numeric scalar in (0, 1].")
    }
  }
  if (!is.logical(fallback_respect_base_filters) || length(fallback_respect_base_filters) != 1 ||
    is.na(fallback_respect_base_filters)) {
    stop("`fallback_respect_base_filters` must be TRUE or FALSE.")
  }
  if (!is.logical(verbose) || length(verbose) != 1 || is.na(verbose)) {
    stop("`verbose` must be TRUE or FALSE.")
  }

  extract_part1 <- function(layer_result) {
    if (is.null(layer_result)) {
      return(NULL)
    }
    if (inherits(layer_result, "HCoCenaLayerResult")) {
      return(as.list(layer_result@part1))
    }
    if (is.list(layer_result) && "part1" %in% names(layer_result)) {
      return(layer_result[["part1"]])
    }
    NULL
  }

  find_cutoff_stats <- function(x) {
    if (is.null(x)) {
      return(NULL)
    }
    if (is.data.frame(x)) {
      has_cols <- all(c("cutoff", "R.squared", "no_edges", "no_nodes", "no_of_networks") %in% colnames(x))
      if (isTRUE(has_cols)) {
        return(x)
      }
      return(NULL)
    }
    if (is.list(x)) {
      if ("cutoff_stats" %in% names(x) && is.data.frame(x[["cutoff_stats"]])) {
        return(x[["cutoff_stats"]])
      }
      for (nm in names(x)) {
        out <- find_cutoff_stats(x[[nm]])
        if (!is.null(out)) {
          return(out)
        }
      }
    }
    NULL
  }

  first_finite <- function(x) {
    x <- .hc_as_numeric_safely(x)
    x <- x[is.finite(x)]
    if (length(x) == 0) {
      return(NA_real_)
    }
    x[[1]]
  }

  compute_slope <- function(degree, probs) {
    d <- .hc_as_numeric_safely(degree)
    p <- .hc_as_numeric_safely(probs)
    keep <- is.finite(d) & is.finite(p) & d > 0 & p > 0
    d <- d[keep]
    p <- p[keep]
    if (length(d) < 2) {
      return(NA_real_)
    }
    fit <- tryCatch(stats::lm(log(p) ~ log(d)), error = function(e) NULL)
    if (is.null(fit)) {
      return(NA_real_)
    }
    cf <- tryCatch(stats::coef(fit), error = function(e) NULL)
    if (is.null(cf) || length(cf) < 2) {
      return(NA_real_)
    }
    slope <- .hc_first_numeric_value(cf[[2]])
    if (!is.finite(slope)) {
      return(NA_real_)
    }
    slope
  }

  summarize_cutoffs <- function(cutoff_stats) {
    if (is.null(cutoff_stats) || !is.data.frame(cutoff_stats) || nrow(cutoff_stats) == 0) {
      return(NULL)
    }
    df <- as.data.frame(cutoff_stats, stringsAsFactors = FALSE)
    for (nm in c("cutoff", "R.squared", "no_edges", "no_nodes", "no_of_networks", "degree", "Probs")) {
      if (nm %in% colnames(df)) {
        df[[nm]] <- .hc_as_numeric_safely(df[[nm]])
      }
    }
    if (!("cutoff" %in% colnames(df))) {
      return(NULL)
    }
    df <- df[is.finite(df$cutoff), , drop = FALSE]
    if (nrow(df) == 0) {
      return(NULL)
    }

    idx_by_cutoff <- split(seq_len(nrow(df)), as.character(df$cutoff))
    rows <- lapply(names(idx_by_cutoff), function(k) {
      sub <- df[idx_by_cutoff[[k]], , drop = FALSE]
      no_edges <- first_finite(sub$no_edges)
      no_nodes <- first_finite(sub$no_nodes)
      mean_k <- if (is.finite(no_edges) && is.finite(no_nodes) && no_nodes > 0) {
        (2 * no_edges) / no_nodes
      } else {
        NA_real_
      }
      data.frame(
        cutoff = .hc_first_numeric_value(k),
        sft_rsq = first_finite(sub$R.squared),
        slope = if (all(c("degree", "Probs") %in% colnames(sub))) compute_slope(sub$degree, sub$Probs) else NA_real_,
        mean_k = mean_k,
        no_edges = no_edges,
        no_nodes = no_nodes,
        no_of_networks = first_finite(sub$no_of_networks),
        stringsAsFactors = FALSE
      )
    })
    out <- do.call(rbind, rows)
    out <- out[is.finite(out$cutoff), , drop = FALSE]
    if (nrow(out) == 0) {
      return(NULL)
    }
    out <- out[order(out$cutoff), , drop = FALSE]
    rownames(out) <- NULL
    out
  }

  select_tiered <- function(summary_df) {
    if (is.null(summary_df) || !is.data.frame(summary_df) || nrow(summary_df) == 0) {
      return(list(
        cutoff = NA_real_,
        tier = "none",
        reason = "no_cutoff_stats",
        strict_cutoff = NA_real_,
        relaxed_cutoff = NA_real_,
        best_available_cutoff = NA_real_,
        table = summary_df
      ))
    }

    df <- summary_df
    neg_ok <- if (isTRUE(require_negative_slope)) is.finite(df$slope) & df$slope < 0 else rep(TRUE, nrow(df))
    comp_ok <- if (is.null(max_no_networks)) rep(TRUE, nrow(df)) else is.finite(df$no_of_networks) & df$no_of_networks <= max_no_networks
    if (is.null(min_node_fraction)) {
      nodes_ok <- rep(TRUE, nrow(df))
      min_nodes_abs <- NA_real_
    } else {
      max_nodes <- .hc_max_finite(df$no_nodes)
      if (!is.finite(max_nodes) || max_nodes <= 0) {
        min_nodes_abs <- NA_real_
        nodes_ok <- rep(FALSE, nrow(df))
      } else {
        min_nodes_abs <- max_nodes * min_node_fraction
        nodes_ok <- is.finite(df$no_nodes) & df$no_nodes >= min_nodes_abs
      }
    }
    base_ok <- neg_ok & comp_ok & nodes_ok

    strict_ok <- base_ok &
      is.finite(df$sft_rsq) & df$sft_rsq >= tier1_sft_min &
      is.finite(df$mean_k) & df$mean_k > tier1_mean_k_min
    relaxed_ok <- base_ok &
      is.finite(df$sft_rsq) & df$sft_rsq >= tier2_sft_min &
      is.finite(df$mean_k) & df$mean_k > tier2_mean_k_min

    df$base_ok <- base_ok
    df$strict_ok <- strict_ok
    df$relaxed_ok <- relaxed_ok
    # keep legacy flags for compatibility with older code paths
    df$tier1_ok <- strict_ok
    df$tier2_ok <- relaxed_ok
    df$min_nodes_required <- min_nodes_abs

    strict_idx <- which(strict_ok)
    strict_cutoff <- if (length(strict_idx) > 0) df$cutoff[[strict_idx[[1]]]] else NA_real_

    relaxed_idx <- which(relaxed_ok)
    relaxed_cutoff <- if (length(relaxed_idx) > 0) df$cutoff[[relaxed_idx[[1]]]] else NA_real_

    best_available_cutoff <- NA_real_
    if (isTRUE(fallback_respect_base_filters)) {
      pool <- which(base_ok & is.finite(df$sft_rsq))
    } else {
      pool <- which(is.finite(df$sft_rsq))
    }
    if (length(pool) > 0) {
      best_pool_idx <- pool[[which.max(df$sft_rsq[pool])]]
      best_sft <- df$sft_rsq[[best_pool_idx]]
      if (is.finite(best_sft) && best_sft > fallback_sft_min) {
        best_available_cutoff <- df$cutoff[[best_pool_idx]]
      }
    }

    if (is.finite(strict_cutoff)) {
      return(list(
        cutoff = strict_cutoff,
        tier = "strict",
        reason = "first_strict_in_ascending_cutoff_order",
        strict_cutoff = strict_cutoff,
        relaxed_cutoff = relaxed_cutoff,
        best_available_cutoff = best_available_cutoff,
        table = df
      ))
    }
    if (is.finite(relaxed_cutoff)) {
      return(list(
        cutoff = relaxed_cutoff,
        tier = "relaxed",
        reason = "first_relaxed_in_ascending_cutoff_order",
        strict_cutoff = strict_cutoff,
        relaxed_cutoff = relaxed_cutoff,
        best_available_cutoff = best_available_cutoff,
        table = df
      ))
    }
    if (is.finite(best_available_cutoff)) {
      return(list(
        cutoff = best_available_cutoff,
        tier = "best_available",
        reason = "max_sft_above_fallback_threshold",
        strict_cutoff = strict_cutoff,
        relaxed_cutoff = relaxed_cutoff,
        best_available_cutoff = best_available_cutoff,
        table = df
      ))
    }
    if (length(pool) == 0) {
      return(list(
        cutoff = NA_real_,
        tier = "none",
        reason = "no_best_available_candidates",
        strict_cutoff = strict_cutoff,
        relaxed_cutoff = relaxed_cutoff,
        best_available_cutoff = best_available_cutoff,
        table = df
      ))
    }

    list(
      cutoff = NA_real_,
      tier = "none",
      reason = "best_available_sft_below_threshold",
      strict_cutoff = strict_cutoff,
      relaxed_cutoff = relaxed_cutoff,
      best_available_cutoff = best_available_cutoff,
      table = df
    )
  }

  layer_ids <- character(0)
  if (nrow(hc@config@layer) > 0 && "layer_id" %in% colnames(hc@config@layer)) {
    layer_ids <- as.character(hc@config@layer$layer_id)
  } else if (length(hc@layer_results) > 0) {
    layer_ids <- names(hc@layer_results)
  }
  if (length(layer_ids) == 0) {
    stop("No layer IDs found in `hc`.")
  }

  recommended <- rep(NA_real_, length(layer_ids))
  names(recommended) <- layer_ids
  strict_cutoffs <- rep(NA_real_, length(layer_ids))
  relaxed_cutoffs <- rep(NA_real_, length(layer_ids))
  best_available_cutoffs <- rep(NA_real_, length(layer_ids))
  names(strict_cutoffs) <- layer_ids
  names(relaxed_cutoffs) <- layer_ids
  names(best_available_cutoffs) <- layer_ids
  tiers <- rep(NA_character_, length(layer_ids))
  reasons <- rep(NA_character_, length(layer_ids))
  names(tiers) <- layer_ids
  names(reasons) <- layer_ids
  details <- vector("list", length(layer_ids))
  names(details) <- layer_ids

  for (i in seq_along(layer_ids)) {
    lid <- layer_ids[[i]]
    layer_result <- hc@layer_results[[lid]]
    part1 <- extract_part1(layer_result)
    cutoff_stats <- find_cutoff_stats(part1)
    summary_df <- summarize_cutoffs(cutoff_stats)
    sel <- select_tiered(summary_df)

    recommended[[i]] <- sel$cutoff
    strict_cutoffs[[i]] <- sel$strict_cutoff
    relaxed_cutoffs[[i]] <- sel$relaxed_cutoff
    best_available_cutoffs[[i]] <- sel$best_available_cutoff
    tiers[[i]] <- sel$tier
    reasons[[i]] <- sel$reason
    details[[i]] <- list(
      selected_cutoff = sel$cutoff,
      selected_policy = sel$tier,
      selected_tier = sel$tier,
      selected_reason = sel$reason,
      strict_cutoff = sel$strict_cutoff,
      relaxed_cutoff = sel$relaxed_cutoff,
      best_available_cutoff = sel$best_available_cutoff,
      cutoff_summary = sel$table
    )

    if (isTRUE(verbose)) {
      message(
        "hc_tune_cutoff(): ", lid,
        " -> cutoff=", if (is.finite(sel$cutoff)) format(sel$cutoff, digits = 4) else "NA",
        " (", sel$tier, "; ", sel$reason, ")"
      )
    }
  }

  simple_vec <- tryCatch(.hc_auto_collect_cutoffs(hc), error = function(e) numeric(0))
  simple_aligned <- rep(NA_real_, length(layer_ids))
  names(simple_aligned) <- layer_ids
  if (length(simple_vec) > 0) {
    n_copy <- min(length(simple_vec), length(simple_aligned))
    simple_aligned[seq_len(n_copy)] <- simple_vec[seq_len(n_copy)]
  }

  comparison <- data.frame(
    layer_id = layer_ids,
    cutoff_tiered = as.numeric(recommended),
    selected_policy = as.character(tiers),
    tier = as.character(tiers),
    reason = as.character(reasons),
    cutoff_strict = as.numeric(strict_cutoffs),
    cutoff_relaxed = as.numeric(relaxed_cutoffs),
    cutoff_best_available = as.numeric(best_available_cutoffs),
    cutoff_simple = as.numeric(simple_aligned),
    delta_tiered_minus_simple = as.numeric(recommended) - as.numeric(simple_aligned),
    stringsAsFactors = FALSE
  )

  applied_cutoffs <- NULL
  apply_status <- "not_applied"
  if (isTRUE(apply)) {
    cut_to_apply <- as.numeric(recommended)
    existing <- rep(NA_real_, length(cut_to_apply))
    if (nrow(hc@config@layer) > 0 && "cutoff" %in% colnames(hc@config@layer)) {
      existing <- .hc_as_numeric_safely(hc@config@layer$cutoff)
      if (length(existing) < length(cut_to_apply)) {
        existing <- c(existing, rep(NA_real_, length(cut_to_apply) - length(existing)))
      }
      existing <- existing[seq_len(length(cut_to_apply))]
    }

    miss <- !is.finite(cut_to_apply)
    if (any(miss)) {
      cut_to_apply[miss] <- existing[miss]
    }

    if (all(is.finite(cut_to_apply))) {
      hc <- hc_set_cutoff(hc, cutoff_vector = cut_to_apply)
      applied_cutoffs <- cut_to_apply
      apply_status <- "applied"
      if (isTRUE(verbose)) {
        message(
          "hc_tune_cutoff(): applied cutoff vector = ",
          paste0(format(cut_to_apply, digits = 4), collapse = ", ")
        )
      }
    } else {
      apply_status <- "apply_failed_missing_cutoffs"
      warning(
        "Could not apply tiered cutoff vector because at least one layer has no finite cutoff ",
        "and no finite fallback from existing `hc@config@layer$cutoff`.",
        call. = FALSE
      )
    }
  }

  report <- list(
    created_at = as.character(Sys.time()),
    status = apply_status,
    parameters = list(
      tier1_sft_min = tier1_sft_min,
      tier1_mean_k_min = tier1_mean_k_min,
      tier2_sft_min = tier2_sft_min,
      tier2_mean_k_min = tier2_mean_k_min,
      fallback_sft_min = fallback_sft_min,
      require_negative_slope = require_negative_slope,
      max_no_networks = max_no_networks,
      max_components_legacy = max_components,
      min_node_fraction = min_node_fraction,
      fallback_respect_base_filters = fallback_respect_base_filters
    ),
    recommended_cutoff_vector = as.numeric(recommended),
    applied_cutoff_vector = applied_cutoffs,
    comparison_with_simple = comparison,
    layer_details = details
  )
  names(report$recommended_cutoff_vector) <- layer_ids
  if (!is.null(report$applied_cutoff_vector)) {
    names(report$applied_cutoff_vector) <- layer_ids
  }

  sat <- as.list(hc@satellite)
  sat[["cutoff_tuning"]] <- report
  sat[["auto_tune_cutoff_tiered"]] <- report
  hc@satellite <- S4Vectors::SimpleList(sat)
  methods::validObject(hc)
  hc
}

#' Robust tiered correlation-cutoff tuning
#'
#' Short-name wrapper around `hc_auto_tune_cutoff_tiered()` for daily use.
#' The report is stored in `hc@satellite$cutoff_tuning`.
#'
#' Key tuning knobs:
#' - `tier1_sft_min`, `tier1_mean_k_min`: strict-tier thresholds.
#' - `tier2_sft_min`, `tier2_mean_k_min`: relaxed-tier thresholds.
#' - `fallback_sft_min`: minimum fit for `best_available` selection.
#' - `require_negative_slope`: enforce negative slope in tier filters.
#' - `max_no_networks`: optional upper bound for connected-component count.
#' - `min_node_fraction`: optional retained-node requirement relative to best
#'   available node count across tested cutoffs.
#' - `fallback_respect_base_filters`: whether `best_available` must respect the same
#'   base filters (`slope/components/nodes`).
#'
#' @inheritParams hc_auto_tune_cutoff_tiered
#' @return Updated `HCoCenaExperiment`.
#' @examples
#' hc <- hc_example_data("after_part1")
#' hc <- hc_tune_cutoff(hc, apply = FALSE, verbose = FALSE)
#' @export
hc_tune_cutoff <- function(hc,
                           apply = FALSE,
                           tier1_sft_min = 0.85,
                           tier1_mean_k_min = 10,
                           tier2_sft_min = 0.80,
                           tier2_mean_k_min = 5,
                           fallback_sft_min = 0.70,
                           require_negative_slope = TRUE,
                           max_no_networks = NULL,
                           max_components = NULL,
                           min_node_fraction = NULL,
                           fallback_respect_base_filters = TRUE,
                           verbose = TRUE) {
  hc_auto_tune_cutoff_tiered(
    hc = hc,
    apply = apply,
    tier1_sft_min = tier1_sft_min,
    tier1_mean_k_min = tier1_mean_k_min,
    tier2_sft_min = tier2_sft_min,
    tier2_mean_k_min = tier2_mean_k_min,
    fallback_sft_min = fallback_sft_min,
    require_negative_slope = require_negative_slope,
    max_no_networks = max_no_networks,
    max_components = max_components,
    min_node_fraction = min_node_fraction,
    fallback_respect_base_filters = fallback_respect_base_filters,
    verbose = verbose
  )
}
