#' Resolve a layer identifier for `HCoCenaExperiment`
#' @noRd
.hc_resolve_layer_id <- function(hc, layer = NULL) {
  exp_names <- base::names(MultiAssayExperiment::experiments(hc@mae))
  if (base::length(exp_names) == 0) {
    stop("No layers found in `hc@mae`. Run `hc_read_data()` first.")
  }

  if (base::is.null(layer)) {
    return(exp_names[[1]])
  }

  if (base::is.numeric(layer) && base::length(layer) == 1) {
    idx <- as.integer(layer)
    if (!is.finite(idx) || idx < 1 || idx > base::length(exp_names)) {
      stop("`layer` index out of range. Available layers: 1..", base::length(exp_names), ".")
    }
    return(exp_names[[idx]])
  }

  layer_chr <- as.character(layer[[1]])
  if (layer_chr %in% exp_names) {
    return(layer_chr)
  }

  cfg <- hc@config@layer
  if (base::nrow(cfg) > 0 &&
    "layer_name" %in% base::colnames(cfg) &&
    "layer_id" %in% base::colnames(cfg)) {
    idx <- base::which(as.character(cfg$layer_name) == layer_chr)
    if (base::length(idx) > 0) {
      cand <- as.character(cfg$layer_id[[idx[[1]]]])
      if (cand %in% exp_names) {
        return(cand)
      }
    }
  }

  stop(
    "Unknown `layer` value: ", layer_chr, ". ",
    "Use a layer index, layer id, or configured layer name."
  )
}

#' Parse module genes from cluster information
#' @noRd
.hc_module_gene_map <- function(cluster_info) {
  req <- c("color", "gene_n")
  if (!base::all(req %in% base::colnames(cluster_info))) {
    stop("Cluster information must contain columns: ", base::paste(req, collapse = ", "), ".")
  }

  split_rows <- base::split(cluster_info, cluster_info$color)
  out <- list()
  for (nm in base::names(split_rows)) {
    block <- split_rows[[nm]]
    gene_chunks <- base::as.character(block$gene_n)
    gene_chunks <- gene_chunks[!base::is.na(gene_chunks)]
    genes <- base::unlist(base::strsplit(gene_chunks, ",", fixed = TRUE), use.names = FALSE)
    genes <- base::trimws(base::as.character(genes))
    genes <- genes[base::nzchar(genes)]
    out[[nm]] <- base::unique(genes)
  }
  out
}

#' Impute missing values along one ordered time series
#' @noRd
.hc_impute_series <- function(x, method = c("none", "linear", "locf")) {
  method <- base::match.arg(method)
  x <- .hc_as_numeric_safely(x)
  if (base::all(base::is.na(x)) || identical(method, "none")) {
    return(x)
  }

  obs <- base::which(!base::is.na(x))
  if (base::length(obs) == 1) {
    x[base::is.na(x)] <- x[[obs[[1]]]]
    return(x)
  }

  if (identical(method, "linear")) {
    idx <- base::seq_along(x)
    return(stats::approx(x = obs, y = x[obs], xout = idx, method = "linear", rule = 2)$y)
  }

  # LOCF + backward fill
  xf <- x
  for (i in base::seq_along(xf)) {
    if (base::is.na(xf[[i]]) && i > 1) {
      xf[[i]] <- xf[[i - 1]]
    }
  }
  if (base::is.na(xf[[1]])) {
    non_na <- base::which(!base::is.na(xf))
    if (base::length(non_na) > 0) {
      xf[[1]] <- xf[[non_na[[1]]]]
    }
  }
  for (i in base::length(xf):2) {
    if (base::is.na(xf[[i - 1]]) && !base::is.na(xf[[i]])) {
      xf[[i - 1]] <- xf[[i]]
    }
  }
  xf
}

#' Calinski-Harabasz index
#' @noRd
.hc_ch_index <- function(x, cluster) {
  x <- base::as.matrix(x)
  cluster <- as.integer(cluster)

  n <- base::nrow(x)
  k <- base::length(base::unique(cluster))
  if (n <= 2 || k <= 1 || k >= n) {
    return(NA_real_)
  }

  grand <- base::colMeans(x)
  wss <- 0
  bss <- 0
  for (cl in base::sort(base::unique(cluster))) {
    idx <- base::which(cluster == cl)
    if (base::length(idx) == 0) {
      next
    }
    xi <- x[idx, , drop = FALSE]
    ci <- base::colMeans(xi)
    centered <- base::sweep(xi, 2, ci, "-")
    wss <- wss + base::sum(centered * centered)
    d <- ci - grand
    bss <- bss + base::length(idx) * base::sum(d * d)
  }
  if (!is.finite(wss) || wss <= 0) {
    return(NA_real_)
  }
  (bss / (k - 1)) / (wss / (n - k))
}

#' Publication-oriented ggplot theme for longitudinal plots
#' @noRd
.hc_theme_pub <- function(base_size = 11) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_line(color = "grey92", linewidth = 0.25),
      panel.grid.major.y = ggplot2::element_line(color = "grey90", linewidth = 0.3),
      strip.background = ggplot2::element_rect(fill = "grey96", color = "grey70"),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold"),
      legend.title = ggplot2::element_text(face = "bold"),
      legend.key = ggplot2::element_rect(fill = "white", color = NA)
    )
}

#' Colorblind-friendly distinct palette
#' @noRd
.hc_distinct_palette <- function(n) {
  base_pal <- c(
    "#1F77B4", # blue
    "#FF7F0E", # orange
    "#2CA02C", # green
    "#D62728", # red
    "#9467BD", # purple
    "#8C564B", # brown
    "#E377C2", # pink
    "#17BECF", # cyan
    "#BCBD22", # olive
    "#7F7F7F" # gray
  )
  n <- base::as.integer(n[[1]])
  if (!base::is.finite(n) || n <= 0L) {
    return(character(0))
  }
  if (n <= base::length(base_pal)) {
    return(base_pal[base::seq_len(n)])
  }
  extra_n <- n - base::length(base_pal)
  extra <- grDevices::hcl.colors(extra_n, palette = "Dark 3")
  c(base_pal, extra)
}

#' Build module-to-color map
#' @noRd
.hc_module_color_map <- function(module_lookup,
                                 module_levels = NULL,
                                 default_col = "grey60") {
  if (base::is.null(module_levels)) {
    module_levels <- character(0)
  } else {
    module_levels <- base::as.character(module_levels)
  }

  if (is.null(module_lookup) || !base::is.data.frame(module_lookup) || base::nrow(module_lookup) == 0) {
    if (base::length(module_levels) == 0) {
      return(stats::setNames(character(0), character(0)))
    }
    cols <- base::rep(default_col, base::length(module_levels))
    return(stats::setNames(cols, module_levels))
  }

  module_lookup <- base::as.data.frame(module_lookup, stringsAsFactors = FALSE)
  if (!"module" %in% base::colnames(module_lookup)) {
    if (base::length(module_levels) == 0) {
      return(stats::setNames(character(0), character(0)))
    }
    cols <- base::rep(default_col, base::length(module_levels))
    return(stats::setNames(cols, module_levels))
  }
  if (!"module_color" %in% base::colnames(module_lookup)) {
    module_lookup$module_color <- default_col
  }

  mods <- base::as.character(module_lookup$module)
  cols <- base::as.character(module_lookup$module_color)
  cols[base::is.na(cols) | !base::nzchar(cols)] <- default_col
  map <- stats::setNames(cols, mods)

  if (base::length(module_levels) > 0) {
    map <- map[module_levels]
    map[base::is.na(map) | !base::nzchar(base::as.character(map))] <- default_col
    map <- stats::setNames(base::as.character(map), module_levels)
  }
  map
}

#' Natural ordering for module labels including split suffixes
#' @noRd
.hc_natural_module_order <- function(modules) {
  modules <- base::unique(base::as.character(modules))
  modules <- modules[!base::is.na(modules) & base::nzchar(base::trimws(modules))]
  if (base::length(modules) <= 1) {
    return(modules)
  }

  parsed <- lapply(modules, function(x) {
    x <- base::as.character(x[[1]])
    prefix <- if (base::grepl("[0-9]", x)) {
      base::sub("^([^0-9]*).*", "\\1", x, perl = TRUE)
    } else {
      x
    }
    number_hits <- regmatches(x, gregexpr("[0-9]+", x, perl = TRUE))[[1]]
    numbers <- .hc_as_integer_safely(number_hits)
    numbers <- numbers[base::is.finite(numbers)]
    suffix <- base::tolower(gsub("[0-9]+", "", x, perl = TRUE))
    list(
      prefix = base::tolower(prefix),
      numbers = numbers,
      suffix = suffix
    )
  })

  max_numbers <- base::max(vapply(parsed, function(x) base::length(x$numbers), integer(1)))
  num_mat <- if (max_numbers > 0) {
    base::matrix(-1, nrow = base::length(parsed), ncol = max_numbers)
  } else {
    NULL
  }
  if (!base::is.null(num_mat)) {
    for (i in base::seq_along(parsed)) {
      nums <- parsed[[i]]$numbers
      if (base::length(nums) > 0) {
        num_mat[i, base::seq_along(nums)] <- nums
      }
    }
  }

  order_args <- list(
    base::vapply(parsed, function(x) x$prefix, FUN.VALUE = base::character(1))
  )
  if (!base::is.null(num_mat)) {
    for (j in base::seq_len(base::ncol(num_mat))) {
      order_args[[base::length(order_args) + 1L]] <- num_mat[, j]
    }
  }
  order_args[[base::length(order_args) + 1L]] <- base::vapply(parsed, function(x) x$suffix, FUN.VALUE = base::character(1))
  order_args[[base::length(order_args) + 1L]] <- base::tolower(modules)
  order_args[[base::length(order_args) + 1L]] <- modules

  ord <- base::do.call(base::order, order_args)
  modules[ord]
}

#' Reorder module lookup by natural module order
#' @noRd
.hc_reorder_module_lookup_natural <- function(module_lookup) {
  if (base::is.null(module_lookup) || !base::is.data.frame(module_lookup) || base::nrow(module_lookup) == 0) {
    return(module_lookup)
  }
  module_lookup <- base::as.data.frame(module_lookup, stringsAsFactors = FALSE)
  if (!"module" %in% base::colnames(module_lookup)) {
    return(module_lookup)
  }
  module_lookup$module <- base::as.character(module_lookup$module)
  ord_levels <- .hc_natural_module_order(module_lookup$module)
  ord <- base::match(ord_levels, module_lookup$module)
  ord <- ord[!base::is.na(ord)]
  module_lookup[ord, , drop = FALSE]
}

#' Check whether labels look like structured module ids
#' @noRd
.hc_is_structured_module_label <- function(x) {
  x <- base::as.character(x)
  !base::is.na(x) &
    base::nzchar(base::trimws(x)) &
    !.hc_is_color_token(x) &
    base::grepl("^[A-Za-z]+[0-9]+(?:\\.[0-9]+)*$", x)
}

#' Normalize longitudinal module labels while preserving split suffixes
#' @noRd
.hc_normalize_longitudinal_module_labels <- function(module_labels,
                                                     module_colors) {
  module_labels <- base::as.character(module_labels)
  module_colors <- base::as.character(module_colors)
  label_ok <- .hc_is_structured_module_label(module_labels)
  if (base::all(label_ok)) {
    return(module_labels)
  }

  prefix <- "M"
  if (base::any(label_ok)) {
    pref_guess <- base::sub(
      "^([A-Za-z]+)[0-9]+(?:\\.[0-9]+)*$",
      "\\1",
      module_labels[label_ok][[1]]
    )
    if (base::length(pref_guess) == 1 && base::nzchar(pref_guess)) {
      prefix <- pref_guess
    }
  }

  base::paste0(prefix, base::seq_along(module_colors))
}

#' Resolve module label map orientation for a given set of module colors
#' @noRd
.hc_resolve_module_label_map_for_colors <- function(label_map,
                                                    module_colors = NULL) {
  if (base::is.null(label_map) || base::length(label_map) == 0) {
    return(NULL)
  }

  out <- base::as.character(label_map)
  map_names <- base::names(label_map)
  if (!base::is.null(map_names) && base::length(map_names) == base::length(out)) {
    base::names(out) <- base::as.character(map_names)
  }

  if (base::is.null(module_colors) || base::length(module_colors) == 0 || base::is.null(base::names(out))) {
    return(out)
  }

  module_colors <- base::as.character(module_colors)
  missing_before <- base::setdiff(module_colors, base::names(out))
  if (base::length(missing_before) > 0) {
    inverse_map <- stats::setNames(base::names(out), base::as.character(out))
    if (base::all(module_colors %in% base::names(inverse_map))) {
      out <- inverse_map
    }
  }

  out
}

#' Preserve module color order from cluster information
#' @noRd
.hc_module_colors_in_cluster_order <- function(cluster_info,
                                               module_genes = NULL) {
  if (!"color" %in% base::colnames(cluster_info)) {
    stop("Cluster information must contain a `color` column.")
  }

  out <- base::unique(base::as.character(cluster_info$color))
  out <- out[!base::is.na(out) & base::nzchar(base::trimws(out))]
  if (!base::is.null(module_genes)) {
    out <- out[out %in% base::names(module_genes)]
  }
  out
}

#' Align module labels/colors with observed enrichment tables
#' @noRd
.hc_harmonize_module_lookup_with_enrichment <- function(module_lookup,
                                                        enrich_df) {
  base_lookup <- if (base::is.data.frame(module_lookup) && base::nrow(module_lookup) > 0) {
    base::as.data.frame(module_lookup, stringsAsFactors = FALSE)
  } else {
    base::data.frame(
      module = base::character(0),
      module_color = base::character(0),
      stringsAsFactors = FALSE
    )
  }
  if (!"module" %in% base::colnames(base_lookup)) {
    base_lookup$module <- base::character(base::nrow(base_lookup))
  }
  if (!"module_color" %in% base::colnames(base_lookup)) {
    base_lookup$module_color <- base::character(base::nrow(base_lookup))
  }
  base_lookup$module <- base::as.character(base_lookup$module)
  base_lookup$module_color <- base::as.character(base_lookup$module_color)
  base_lookup <- base_lookup[
    !base::is.na(base_lookup$module) & base::nzchar(base::trimws(base_lookup$module)) &
      !base::is.na(base_lookup$module_color) & base::nzchar(base::trimws(base_lookup$module_color)),
    c("module", "module_color"),
    drop = FALSE
  ]
  base_lookup <- base_lookup[!base::duplicated(base::paste(base_lookup$module, base_lookup$module_color, sep = "\r")), , drop = FALSE]

  if (!base::is.data.frame(enrich_df) || base::nrow(enrich_df) == 0) {
    return(base_lookup)
  }

  obs <- base::as.data.frame(enrich_df, stringsAsFactors = FALSE)
  if (!base::all(c("cluster", "module_label") %in% base::colnames(obs))) {
    return(base_lookup)
  }
  obs <- obs[, c("cluster", "module_label"), drop = FALSE]
  base::colnames(obs) <- c("module_color", "module")
  obs$module_color <- base::as.character(obs$module_color)
  obs$module <- base::trimws(base::as.character(obs$module))
  obs$module[base::is.na(obs$module)] <- ""
  obs$module[base::vapply(obs$module, .hc_is_color_token, FUN.VALUE = base::logical(1))] <- ""

  if (base::nrow(base_lookup) > 0) {
    fill_by_color <- !base::nzchar(obs$module)
    if (base::any(fill_by_color)) {
      color_map <- stats::setNames(base_lookup$module, base_lookup$module_color)
      obs$module[fill_by_color] <- base::as.character(color_map[obs$module_color[fill_by_color]])
    }
    fill_by_module <- base::is.na(obs$module_color) | !base::nzchar(base::trimws(obs$module_color))
    if (base::any(fill_by_module)) {
      module_map <- stats::setNames(base_lookup$module_color, base_lookup$module)
      obs$module_color[fill_by_module] <- base::as.character(module_map[obs$module[fill_by_module]])
    }
  }

  obs <- obs[
    !base::is.na(obs$module) & base::nzchar(base::trimws(obs$module)) &
      !base::is.na(obs$module_color) & base::nzchar(base::trimws(obs$module_color)), ,
    drop = FALSE
  ]
  if (base::nrow(obs) == 0) {
    return(base_lookup)
  }
  obs <- obs[!base::duplicated(obs$module), , drop = FALSE]

  if (base::nrow(base_lookup) == 0) {
    base::rownames(obs) <- NULL
    return(obs)
  }

  keep_base <- !(base_lookup$module %in% obs$module | base_lookup$module_color %in% obs$module_color)
  out <- base::rbind(obs, base_lookup[keep_base, c("module", "module_color"), drop = FALSE])
  out <- out[!base::duplicated(out$module), , drop = FALSE]
  base::rownames(out) <- NULL
  out
}

#' Detect whether labels look like literal color tokens
#' @noRd
.hc_is_color_token <- function(x) {
  x <- base::as.character(x)
  x <- base::trimws(x)
  x_low <- base::tolower(x)
  is_hex <- grepl("^#(?:[0-9a-f]{3}|[0-9a-f]{6})$", x_low)
  is_named <- x_low %in% base::tolower(grDevices::colors())
  is_hex | is_named
}

#' Remove common enrichment DB prefixes from term labels
#' @noRd
.hc_clean_enrichment_term_label <- function(x) {
  x_in <- base::as.character(x)
  x_out <- base::trimws(x_in)
  pats <- c(
    "^HALLMARK_",
    "^KEGG_",
    "^REACTOME_",
    "^WIKIPATHWAYS_",
    "^GOBP_",
    "^GO_BP_",
    "^GOCC_",
    "^GOMF_",
    "^GO_"
  )
  for (pat in pats) {
    x_out <- gsub(pat, "", x_out, perl = TRUE)
  }
  x_out <- gsub("_+", " ", x_out, perl = TRUE)
  x_out <- base::trimws(x_out)
  empty_idx <- base::is.na(x_out) | x_out == ""
  x_out[empty_idx] <- x_in[empty_idx]
  x_out
}

#' Match enrichment terms against user-provided custom term filters
#' @noRd
.hc_match_custom_enrichment_terms <- function(term,
                                              module,
                                              custom_terms = NULL,
                                              match_mode = c("exact", "contains", "regex")) {
  match_mode <- base::match.arg(match_mode)
  term <- base::as.character(term)
  module <- base::as.character(module)
  term[base::is.na(term)] <- ""
  module[base::is.na(module)] <- ""
  term_clean <- .hc_clean_enrichment_term_label(term)

  if (base::is.null(custom_terms)) {
    return(base::rep(TRUE, base::length(term)))
  }

  .norm <- function(x) {
    base::toupper(base::trimws(base::as.character(x)))
  }
  .prepare_queries <- function(q) {
    q <- base::as.character(q)
    q <- q[!base::is.na(q)]
    q <- base::trimws(q)
    q[q != ""]
  }
  .match_vec <- function(raw_vec, clean_vec, queries) {
    q <- .prepare_queries(queries)
    if (base::length(q) == 0) {
      return(base::rep(FALSE, base::length(raw_vec)))
    }
    if (identical(match_mode, "exact")) {
      q_raw <- .norm(q)
      q_clean <- .norm(.hc_clean_enrichment_term_label(q))
      raw_u <- .norm(raw_vec)
      clean_u <- .norm(clean_vec)
      return(raw_u %in% q_raw | clean_u %in% q_raw | raw_u %in% q_clean | clean_u %in% q_clean)
    }
    out <- base::rep(FALSE, base::length(raw_vec))
    for (pat in q) {
      cur <- tryCatch(
        {
          if (identical(match_mode, "contains")) {
            grepl(pat, raw_vec, ignore.case = TRUE, fixed = TRUE) |
              grepl(pat, clean_vec, ignore.case = TRUE, fixed = TRUE)
          } else {
            grepl(pat, raw_vec, ignore.case = TRUE, perl = TRUE) |
              grepl(pat, clean_vec, ignore.case = TRUE, perl = TRUE)
          }
        },
        error = function(e) {
          base::rep(FALSE, base::length(raw_vec))
        }
      )
      out <- out | cur
    }
    out
  }

  if (base::is.character(custom_terms)) {
    return(.match_vec(term, term_clean, custom_terms))
  }

  if (!base::is.list(custom_terms)) {
    stop("`custom_terms` must be NULL, a character vector, or a named list.")
  }

  ct_names <- base::names(custom_terms)
  if (base::is.null(ct_names) || !base::any(base::nzchar(base::trimws(ct_names)))) {
    flat <- base::unlist(custom_terms, use.names = FALSE)
    return(.match_vec(term, term_clean, flat))
  }

  keep <- base::rep(FALSE, base::length(term))
  nm <- base::trimws(base::as.character(ct_names))
  nm_low <- base::tolower(nm)
  global_idx <- base::which(!base::nzchar(nm) | nm_low %in% c("all", "*", "default", ".default"))
  if (base::length(global_idx) > 0) {
    q_global <- base::unlist(custom_terms[global_idx], use.names = FALSE)
    keep <- keep | .match_vec(term, term_clean, q_global)
  }

  mod_low <- base::tolower(module)
  local_idx <- base::setdiff(base::seq_along(custom_terms), global_idx)
  for (ii in local_idx) {
    key <- nm[[ii]]
    if (!base::nzchar(key)) {
      next
    }
    idx_mod <- mod_low == base::tolower(key)
    if (!base::any(idx_mod)) {
      next
    }
    keep[idx_mod] <- keep[idx_mod] | .match_vec(term[idx_mod], term_clean[idx_mod], custom_terms[[ii]])
  }
  keep
}

#' Build module faceting with module-colored strip titles
#' @noRd
.hc_module_facet_wrap <- function(module_levels,
                                  mod_col,
                                  ncol = 4,
                                  scales = "free_y",
                                  facets = ~module) {
  module_levels <- base::as.character(module_levels)
  mod_col <- base::as.character(mod_col[module_levels])
  mod_col[base::is.na(mod_col) | !base::nzchar(mod_col)] <- "grey60"

  txt_col <- base::rep("black", base::length(mod_col))
  if (base::length(mod_col) > 0) {
    rgb <- grDevices::col2rgb(mod_col) / 255
    lum <- 0.2126 * rgb[1, ] + 0.7152 * rgb[2, ] + 0.0722 * rgb[3, ]
    txt_col <- ifelse(lum < 0.55, "white", "black")
  }

  if (requireNamespace("ggh4x", quietly = TRUE)) {
    return(
      ggh4x::facet_wrap2(
        facets = facets,
        ncol = ncol,
        scales = scales,
        strip = ggh4x::strip_themed(
          background_x = ggh4x::elem_list_rect(fill = mod_col, color = "grey35"),
          text_x = ggh4x::elem_list_text(colour = txt_col, face = "bold")
        )
      )
    )
  }

  ggplot2::facet_wrap(facets = facets, ncol = ncol, scales = scales)
}

#' Apply module-colored strip backgrounds to faceted plots
#' @noRd
.hc_colorize_module_strips <- function(p, module_levels, mod_col) {
  if (!inherits(p, "ggplot")) {
    return(p)
  }
  if (requireNamespace("ggh4x", quietly = TRUE)) {
    return(p)
  }

  module_levels <- base::as.character(module_levels)
  if (base::length(module_levels) == 0) {
    return(p)
  }

  col_map <- stats::setNames(base::as.character(mod_col[module_levels]), module_levels)
  col_map[base::is.na(col_map) | !base::nzchar(col_map)] <- "grey60"
  rgb <- grDevices::col2rgb(col_map) / 255
  lum <- 0.2126 * rgb[1, ] + 0.7152 * rgb[2, ] + 0.0722 * rgb[3, ]
  txt_map <- stats::setNames(ifelse(lum < 0.55, "white", "black"), base::names(col_map))

  g <- ggplot2::ggplotGrob(p)
  lay <- g[["layout"]]
  if (is.null(lay) || !"name" %in% base::colnames(lay)) {
    return(p)
  }

  strip_idx <- base::which(grepl("^strip-", lay[["name"]]))
  if (base::length(strip_idx) == 0) {
    return(p)
  }

  for (ii in strip_idx) {
    sg <- g[["grobs"]][[ii]]
    if (!inherits(sg, "gtable") || base::length(sg[["grobs"]]) == 0) {
      next
    }
    gt <- sg[["grobs"]][[1]]
    if (is.null(gt[["children"]]) || base::length(gt[["children"]]) == 0) {
      next
    }

    child_names <- base::names(gt[["children"]])
    i_bg <- base::which(grepl("strip\\.background", child_names))[1]
    i_tx <- base::which(grepl("strip\\.text", child_names))[1]
    if (base::is.na(i_bg) || base::is.na(i_tx)) {
      next
    }

    label <- NA_character_
    text_g <- gt[["children"]][[i_tx]]
    if (!is.null(text_g[["children"]]) && base::length(text_g[["children"]]) > 0) {
      tx_child <- text_g[["children"]][[1]]
      if (!is.null(tx_child[["label"]])) {
        label <- base::as.character(tx_child[["label"]][[1]])
      }
    }
    if (is.na(label) || !base::nzchar(label)) {
      next
    }
    label_key <- label
    if (!(label_key %in% base::names(col_map))) {
      # Support strip labels like "M1: TERM_NAME" by mapping via the prefix.
      prefix <- base::trimws(base::strsplit(label, ":", fixed = TRUE)[[1]][1])
      if (base::nzchar(prefix) && prefix %in% base::names(col_map)) {
        label_key <- prefix
      }
    }
    if (!(label_key %in% base::names(col_map))) {
      next
    }

    fill_col <- col_map[[label_key]]
    txt_col <- txt_map[[label_key]]
    if (is.na(fill_col) || !base::nzchar(fill_col)) {
      fill_col <- "grey60"
    }
    if (is.na(txt_col) || !base::nzchar(txt_col)) {
      txt_col <- "black"
    }

    bg_g <- gt[["children"]][[i_bg]]
    if (is.null(bg_g[["gp"]])) {
      bg_g[["gp"]] <- grid::gpar()
    }
    bg_g[["gp"]][["fill"]] <- fill_col
    bg_g[["gp"]][["col"]] <- "grey35"
    gt[["children"]][[i_bg]] <- bg_g

    if (!is.null(text_g[["children"]]) && base::length(text_g[["children"]]) > 0) {
      for (kk in base::seq_along(text_g[["children"]])) {
        tx <- text_g[["children"]][[kk]]
        if (is.null(tx[["gp"]])) {
          tx[["gp"]] <- grid::gpar()
        }
        tx[["gp"]][["col"]] <- txt_col
        tx[["gp"]][["font"]] <- 2
        text_g[["children"]][[kk]] <- tx
      }
    }
    gt[["children"]][[i_tx]] <- text_g
    sg[["grobs"]][[1]] <- gt
    g[["grobs"]][[ii]] <- sg
  }

  if (requireNamespace("patchwork", quietly = TRUE)) {
    return(patchwork::wrap_elements(full = g))
  }
  g
}

#' Keep only groups with enough points for ellipse overlays
#' @noRd
.hc_filter_groups_for_ellipse <- function(df,
                                          group_col,
                                          min_points = 3L,
                                          x_col = NULL,
                                          y_col = NULL,
                                          min_unique = 3L) {
  if (!base::is.data.frame(df) || base::nrow(df) == 0 || !(group_col %in% base::colnames(df))) {
    return(df[0, , drop = FALSE])
  }

  grp <- base::as.character(df[[group_col]])
  cnt <- stats::ave(base::rep_len(1L, base::length(grp)), grp, FUN = base::length)
  keep <- !base::is.na(grp) & !base::is.na(cnt) & cnt >= base::as.integer(min_points[[1]])

  if (!base::is.null(x_col) && !base::is.null(y_col) &&
    x_col %in% base::colnames(df) && y_col %in% base::colnames(df)) {
    x <- .hc_as_numeric_safely(df[[x_col]])
    y <- .hc_as_numeric_safely(df[[y_col]])
    finite_xy <- base::is.finite(x) & base::is.finite(y)
    keep <- keep & finite_xy

    if (base::any(keep)) {
      keep_groups <- base::unique(grp[keep])
      ok_groups <- keep_groups[base::vapply(
        keep_groups,
        function(g) {
          idx <- keep & grp == g
          xx <- x[idx]
          yy <- y[idx]
          if (base::length(xx) < base::as.integer(min_points[[1]])) {
            return(FALSE)
          }
          if (base::length(base::unique(xx)) < base::as.integer(min_unique[[1]]) ||
            base::length(base::unique(yy)) < base::as.integer(min_unique[[1]])) {
            return(FALSE)
          }
          vx <- stats::var(xx, na.rm = TRUE)
          vy <- stats::var(yy, na.rm = TRUE)
          if (!(base::is.finite(vx) && base::is.finite(vy) && vx > 0 && vy > 0)) {
            return(FALSE)
          }
          cov_xy <- stats::cov(base::cbind(xx, yy), use = "complete.obs")
          det_xy <- base::det(cov_xy)
          base::is.finite(det_xy) && det_xy > 1e-10
        },
        logical(1)
      )]
      keep <- keep & (grp %in% ok_groups)
    }
  }

  out <- df[keep, , drop = FALSE]
  if (group_col %in% base::colnames(out) && base::is.factor(out[[group_col]])) {
    out[[group_col]] <- base::droplevels(out[[group_col]])
  }
  out
}

#' Compute 2D cluster centers for label placement
#' @noRd
.hc_cluster_centers_2d <- function(df, x_col, y_col, group_col) {
  req <- c(x_col, y_col, group_col)
  if (!base::is.data.frame(df) || base::nrow(df) == 0 || !base::all(req %in% base::colnames(df))) {
    return(base::data.frame(
      x = numeric(0),
      y = numeric(0),
      group = character(0),
      stringsAsFactors = FALSE
    ))
  }
  out <- stats::aggregate(
    df[, c(x_col, y_col), drop = FALSE],
    by = list(group = base::as.character(df[[group_col]])),
    FUN = base::mean,
    na.rm = TRUE
  )
  base::colnames(out)[base::colnames(out) == x_col] <- "x"
  base::colnames(out)[base::colnames(out) == y_col] <- "y"
  out
}

.hc_square_limits_2d <- function(df, x_col, y_col, pad_frac = 0.06) {
  if (!base::is.data.frame(df) || base::nrow(df) == 0 ||
    !(x_col %in% base::colnames(df)) || !(y_col %in% base::colnames(df))) {
    return(list(x = NULL, y = NULL))
  }
  x <- .hc_as_numeric_safely(df[[x_col]])
  y <- .hc_as_numeric_safely(df[[y_col]])
  keep <- base::is.finite(x) & base::is.finite(y)
  if (!base::any(keep)) {
    return(list(x = NULL, y = NULL))
  }
  x <- x[keep]
  y <- y[keep]
  xr <- base::range(x, na.rm = TRUE)
  yr <- base::range(y, na.rm = TRUE)
  xmid <- base::mean(xr)
  ymid <- base::mean(yr)
  xspan <- diff(xr)
  yspan <- diff(yr)
  span <- base::max(xspan, yspan, na.rm = TRUE)
  if (!base::is.finite(span) || span <= 0) {
    span <- 1
  }
  span <- span * (1 + pad_frac)
  half <- span / 2
  list(
    x = c(xmid - half, xmid + half),
    y = c(ymid - half, ymid + half)
  )
}

.hc_attach_right_legend <- function(plot_obj, legend_width = 0.32) {
  if (!inherits(plot_obj, "ggplot")) {
    return(plot_obj)
  }
  if (!requireNamespace("cowplot", quietly = TRUE)) {
    return(plot_obj)
  }
  legend_plot <- plot_obj + ggplot2::theme(legend.position = "right")
  legend_grob <- cowplot::get_legend(legend_plot)
  if (is.null(legend_grob)) {
    return(plot_obj)
  }
  plot_noleg <- plot_obj + ggplot2::theme(legend.position = "none")
  cowplot::plot_grid(
    plot_noleg,
    cowplot::ggdraw(legend_grob),
    nrow = 1,
    rel_widths = c(1, legend_width),
    align = "h"
  )
}

.hc_order_module_levels <- function(x) {
  x <- base::as.character(x)
  x <- x[!base::is.na(x) & base::nzchar(x)]
  x <- base::unique(x)
  if (base::length(x) <= 1L) {
    return(x)
  }
  num_part <- .hc_as_integer_safely(gsub("^[^0-9]*([0-9]+).*$", "\\1", x))
  has_num <- base::is.finite(num_part)
  ord <- base::order(
    !has_num,
    ifelse(has_num, gsub("^([A-Za-z]+).*", "\\1", x), x),
    ifelse(has_num, num_part, Inf),
    x
  )
  x[ord]
}

.hc_panel_mean_limits <- function(mean_df,
                                  panel_col = "panel_label",
                                  value_col = "score",
                                  lower_col = NULL,
                                  upper_col = NULL,
                                  pad_frac = 0.12,
                                  min_span = 0.04) {
  if (!base::is.data.frame(mean_df) || base::nrow(mean_df) == 0 ||
    !(panel_col %in% base::colnames(mean_df)) ||
    !(value_col %in% base::colnames(mean_df))) {
    return(base::data.frame(
      panel_label = character(0),
      y_lower = numeric(0),
      y_upper = numeric(0),
      stringsAsFactors = FALSE
    ))
  }
  panels <- base::as.character(mean_df[[panel_col]])
  val <- .hc_as_numeric_safely(mean_df[[value_col]])
  low <- if (!base::is.null(lower_col) && lower_col %in% base::colnames(mean_df)) {
    .hc_as_numeric_safely(mean_df[[lower_col]])
  } else {
    val
  }
  high <- if (!base::is.null(upper_col) && upper_col %in% base::colnames(mean_df)) {
    .hc_as_numeric_safely(mean_df[[upper_col]])
  } else {
    val
  }
  keep <- !base::is.na(panels) & base::nzchar(panels) &
    (base::is.finite(val) | base::is.finite(low) | base::is.finite(high))
  if (!base::any(keep)) {
    return(base::data.frame(
      panel_label = character(0),
      y_lower = numeric(0),
      y_upper = numeric(0),
      stringsAsFactors = FALSE
    ))
  }
  tmp <- base::data.frame(
    panel_label = panels[keep],
    value = val[keep],
    lower = low[keep],
    upper = high[keep],
    stringsAsFactors = FALSE
  )
  split_df <- base::split(tmp, tmp$panel_label)
  out <- lapply(split_df, function(df) {
    lower_now <- min(df$lower[base::is.finite(df$lower)], df$value[base::is.finite(df$value)], na.rm = TRUE)
    upper_now <- max(df$upper[base::is.finite(df$upper)], df$value[base::is.finite(df$value)], na.rm = TRUE)
    span <- upper_now - lower_now
    if (!base::is.finite(span) || span < min_span) {
      center <- base::mean(c(lower_now, upper_now), na.rm = TRUE)
      half <- min_span / 2
      lower_now <- center - half
      upper_now <- center + half
      span <- upper_now - lower_now
    }
    pad <- span * pad_frac
    base::data.frame(
      panel_label = base::as.character(df$panel_label[[1]]),
      y_lower = lower_now - pad,
      y_upper = upper_now + pad,
      stringsAsFactors = FALSE
    )
  })
  base::do.call(base::rbind, out)
}

.hc_normalize_emmeans_summary <- function(df,
                                          estimate_candidates = c("emmean", "response"),
                                          lower_candidates = c("lower.CL", "asymp.LCL", "lower.HPD", "lower"),
                                          upper_candidates = c("upper.CL", "asymp.UCL", "upper.HPD", "upper")) {
  if (!base::is.data.frame(df) || base::nrow(df) == 0) {
    return(df)
  }

  pick_first_col <- function(candidates) {
    hits <- candidates[candidates %in% base::colnames(df)]
    if (base::length(hits) == 0) {
      return(NULL)
    }
    hits[[1]]
  }

  est_col <- pick_first_col(estimate_candidates)
  if (base::is.null(est_col)) {
    stop("Could not find an estimate column in `emmeans` summary output.")
  }
  low_col <- pick_first_col(lower_candidates)
  up_col <- pick_first_col(upper_candidates)

  df$emmean <- .hc_as_numeric_safely(df[[est_col]])
  if (!base::is.null(low_col)) {
    df$lower.CL <- .hc_as_numeric_safely(df[[low_col]])
  } else {
    df$lower.CL <- df$emmean
  }
  if (!base::is.null(up_col)) {
    df$upper.CL <- .hc_as_numeric_safely(df[[up_col]])
  } else {
    df$upper.CL <- df$emmean
  }

  df
}

.hc_align_longitudinal_panel_factors <- function(df,
                                                 panel_levels,
                                                 term_display_levels = NULL) {
  if (!base::is.data.frame(df) || base::nrow(df) == 0) {
    return(df)
  }

  is_valid <- function(x) {
    x <- base::as.character(x)
    !base::is.na(x) & base::nzchar(base::trimws(x)) & x != "NA"
  }

  if ("panel_label" %in% base::colnames(df)) {
    panel_chr <- base::as.character(df$panel_label)
    keep <- is_valid(panel_chr)
    if (!base::is.null(panel_levels) && base::length(panel_levels) > 0) {
      keep <- keep & panel_chr %in% base::as.character(panel_levels)
    }
    df <- df[keep, , drop = FALSE]
    panel_chr <- base::as.character(df$panel_label)
    df$panel_label <- base::factor(panel_chr, levels = base::as.character(panel_levels))
    df <- df[!base::is.na(df$panel_label), , drop = FALSE]
  }

  if ("term_display" %in% base::colnames(df) &&
    !base::is.null(term_display_levels) &&
    base::length(term_display_levels) > 0) {
    term_chr <- base::as.character(df$term_display)
    term_chr[!is_valid(term_chr)] <- "Term unavailable"
    df$term_display <- base::factor(term_chr, levels = base::as.character(term_display_levels))
  }

  df
}

#' Build longitudinal meta-clustering report tables
#' @noRd
.hc_meta_report_tables <- function(obj, detail = c("simple", "full")) {
  detail <- base::match.arg(detail)
  if (is.null(obj) || is.null(obj$meta_cluster)) {
    return(list())
  }

  meta_cluster <- base::as.data.frame(obj$meta_cluster, stringsAsFactors = FALSE)
  n_donors <- base::nrow(meta_cluster)
  n_meta <- base::length(base::unique(meta_cluster$meta_cluster))
  selected_k <- if (!is.null(obj$meta_best_k)) as.integer(obj$meta_best_k[[1]]) else NA_integer_
  selected_method <- if (!is.null(obj$meta_method)) as.character(obj$meta_method[[1]]) else NA_character_
  use_consensus <- isTRUE(obj$meta_consensus)

  score_table <- if (!is.null(obj$meta_score_table)) {
    base::as.data.frame(obj$meta_score_table, stringsAsFactors = FALSE)
  } else {
    NULL
  }

  selected_score <- NA_real_
  score_name <- "ch_index"
  selection_rule <- "Maximum CH index (tie-break: smaller k)."
  consensus_explanation <- "Consensus not used."

  if (!is.null(score_table) && base::nrow(score_table) > 0) {
    if (use_consensus && "consensus_delta" %in% base::colnames(score_table)) {
      score_name <- "consensus_delta"
      selection_rule <- "Maximum consensus_delta (within-between); tie-break: lower PAC, then smaller k."
      consensus_explanation <- base::paste0(
        "Consensus uses repeated donor/feature subsampling per k (runs=",
        as.integer(obj$meta_consensus_runs[[1]]),
        "). For each k, higher within-cluster and lower between-cluster co-assignment are preferred."
      )
      idx <- base::match(selected_k, score_table$k)
      if (!base::is.na(idx)) {
        selected_score <- as.numeric(score_table$consensus_delta[[idx]])
      }
    } else {
      idx <- base::match(selected_k, score_table$k)
      if (!base::is.na(idx) && "ch_index" %in% base::colnames(score_table)) {
        selected_score <- as.numeric(score_table$ch_index[[idx]])
      }
    }
  }

  overview <- base::data.frame(
    selected_method = selected_method,
    selected_k = selected_k,
    n_meta_clusters = n_meta,
    n_donors = n_donors,
    feature_source = if (!is.null(obj$meta_feature_source)) as.character(obj$meta_feature_source[[1]]) else NA_character_,
    consensus_used = use_consensus,
    selected_score_name = score_name,
    selected_score = selected_score,
    selection_rule = selection_rule,
    consensus_runs = if (use_consensus && !is.null(obj$meta_consensus_runs)) as.integer(obj$meta_consensus_runs[[1]]) else NA_integer_,
    consensus_sample_fraction = if (use_consensus && !is.null(obj$meta_consensus_sample_fraction)) as.numeric(obj$meta_consensus_sample_fraction[[1]]) else NA_real_,
    consensus_feature_fraction = if (use_consensus && !is.null(obj$meta_consensus_feature_fraction)) as.numeric(obj$meta_consensus_feature_fraction[[1]]) else NA_real_,
    consensus_linkage = if (use_consensus && !is.null(obj$meta_consensus_linkage)) as.character(obj$meta_consensus_linkage[[1]]) else NA_character_,
    consensus_explanation = consensus_explanation,
    stringsAsFactors = FALSE
  )

  method_comparison <- NULL
  if (!is.null(obj$meta_method_comparison)) {
    method_comparison <- base::as.data.frame(obj$meta_method_comparison, stringsAsFactors = FALSE)
  } else if (!is.null(obj$meta_feature_matrix_used) && !is.null(obj$meta_candidate_k)) {
    x_use <- base::as.matrix(obj$meta_feature_matrix_used)
    k_vec <- base::as.integer(obj$meta_candidate_k)
    nstart <- if (!is.null(obj$meta_nstart)) as.integer(obj$meta_nstart[[1]]) else 50L
    seed <- if (!is.null(obj$meta_seed)) as.integer(obj$meta_seed[[1]]) else 42L
    mr <- list()
    for (mm in c("kmeans", "hclust")) {
      cmp <- tryCatch(
        .hc_cluster_with_scoring(
          x = x_use,
          k = k_vec,
          method = mm,
          nstart = nstart,
          seed = seed
        ),
        error = function(e) e
      )
      if (inherits(cmp, "error")) {
        mr[[base::length(mr) + 1]] <- base::data.frame(
          method = mm,
          best_k_ch = NA_integer_,
          best_ch_index = NA_real_,
          n_clusters_at_best_k = NA_integer_,
          status = base::conditionMessage(cmp),
          stringsAsFactors = FALSE
        )
      } else {
        st <- base::as.data.frame(cmp$score_table, stringsAsFactors = FALSE)
        idx <- base::match(as.integer(cmp$best_k), st$k)
        mr[[base::length(mr) + 1]] <- base::data.frame(
          method = mm,
          best_k_ch = as.integer(cmp$best_k),
          best_ch_index = if (!base::is.na(idx)) as.numeric(st$ch_index[[idx]]) else NA_real_,
          n_clusters_at_best_k = base::length(base::unique(cmp$cluster)),
          status = "ok",
          stringsAsFactors = FALSE
        )
      }
    }
    method_comparison <- base::do.call(base::rbind, mr)
    method_comparison$selected_method <- method_comparison$method == selected_method
    method_comparison$selected_k_pipeline <- ifelse(method_comparison$selected_method, selected_k, NA_integer_)
    method_comparison$n_meta_clusters_pipeline <- ifelse(method_comparison$selected_method, n_meta, NA_integer_)
    method_comparison$selection_mode <- if (use_consensus) "consensus" else "ch_index"
  }

  k_candidates <- NULL
  if (!is.null(score_table)) {
    k_candidates <- score_table
    if (use_consensus && "consensus_delta" %in% base::colnames(k_candidates)) {
      sel <- as.numeric(k_candidates$consensus_delta)
      sel[!base::is.finite(sel)] <- -Inf
      pac <- if ("pac" %in% base::colnames(k_candidates)) as.numeric(k_candidates$pac) else base::rep(Inf, base::nrow(k_candidates))
      pac[!base::is.finite(pac)] <- Inf
      ord <- base::order(-sel, pac, as.numeric(k_candidates$k))
      k_candidates <- k_candidates[ord, , drop = FALSE]
      k_candidates$selection_score <- sel[ord]
      k_candidates$selection_score_name <- "consensus_delta"
    } else if ("ch_index" %in% base::colnames(k_candidates)) {
      sel <- as.numeric(k_candidates$ch_index)
      sel[!base::is.finite(sel)] <- -Inf
      ord <- base::order(-sel, as.numeric(k_candidates$k))
      k_candidates <- k_candidates[ord, , drop = FALSE]
      k_candidates$selection_score <- sel[ord]
      k_candidates$selection_score_name <- "ch_index"
    }
    k_candidates$selection_rank <- base::seq_len(base::nrow(k_candidates))
    k_candidates$selected <- as.integer(k_candidates$k) == selected_k
  }

  cluster_sizes <- NULL
  if (base::nrow(meta_cluster) > 0) {
    size_df <- stats::aggregate(
      donor ~ meta_cluster,
      data = meta_cluster,
      FUN = base::length
    )
    base::colnames(size_df)[base::colnames(size_df) == "donor"] <- "n_donors"
    if ("stability" %in% base::colnames(meta_cluster)) {
      stab_df <- stats::aggregate(
        stability ~ meta_cluster,
        data = meta_cluster,
        FUN = function(x) base::mean(as.numeric(x), na.rm = TRUE)
      )
      base::colnames(stab_df)[base::colnames(stab_df) == "stability"] <- "mean_stability"
      size_df <- base::merge(size_df, stab_df, by = "meta_cluster", all.x = TRUE, sort = FALSE)
    }
    cluster_sizes <- size_df[base::order(size_df$meta_cluster), , drop = FALSE]
  }

  donor_clusterings <- NULL
  if (!is.null(obj$meta_donor_clusterings)) {
    donor_clusterings <- base::as.data.frame(obj$meta_donor_clusterings, stringsAsFactors = FALSE)
  } else if (base::nrow(meta_cluster) > 0) {
    donor_clusterings <- base::data.frame(
      donor = base::as.character(meta_cluster$donor),
      selected_method = selected_method,
      selected_k = selected_k,
      selected_cluster_id = base::as.integer(meta_cluster$cluster),
      selected_cluster = base::as.character(meta_cluster$meta_cluster),
      stringsAsFactors = FALSE
    )
    if ("stability" %in% base::colnames(meta_cluster)) {
      donor_clusterings$selected_stability <- .hc_as_numeric_safely(meta_cluster$stability)
    }
  }

  method_clusters_tbl <- NULL
  if (!base::is.null(method_comparison) && base::nrow(method_comparison) > 0) {
    method_clusters_tbl <- method_comparison[, base::c(
      "method",
      "best_k_ch",
      "n_clusters_at_best_k",
      "selected_method"
    ), drop = FALSE]
    base::colnames(method_clusters_tbl) <- c(
      "clustering",
      "best_k",
      "meta_clusters_found",
      "used_in_pipeline"
    )
  }

  consensus_decision_tbl <- base::data.frame(
    selected_clustering = selected_method,
    selected_k = selected_k,
    selected_meta_clusters = n_meta,
    consensus_used = use_consensus,
    selection_mode = if (use_consensus) "consensus" else "ch_index",
    decision_rule = if (use_consensus) {
      "max(consensus_delta), tie-break: min(PAC), then smaller k"
    } else {
      "max(CH index), tie-break: smaller k"
    },
    stringsAsFactors = FALSE
  )

  if (identical(detail, "simple")) {
    return(list(
      Method_Clusters = method_clusters_tbl,
      Consensus_Decision = consensus_decision_tbl
    ))
  }

  list(
    Meta_Overview = overview,
    Method_Comparison = method_comparison,
    K_Candidates = k_candidates,
    Cluster_Sizes = cluster_sizes,
    Donor_Clusterings = donor_clusterings
  )
}

#' Resolve output directory from hc config
#' @noRd
.hc_resolve_output_dir <- function(hc) {
  out_dir <- "."
  if (base::nrow(hc@config@paths) > 0 && "dir_output" %in% base::colnames(hc@config@paths)) {
    cand <- as.character(hc@config@paths$dir_output[[1]])
    if (base::length(cand) == 1 && !base::is.na(cand) && base::nzchar(cand) && !identical(cand, "FALSE")) {
      out_dir <- cand
    }
  }
  if (base::nrow(hc@config@global) > 0 && "save_folder" %in% base::colnames(hc@config@global)) {
    sf <- as.character(hc@config@global$save_folder[[1]])
    if (base::length(sf) == 1 && !base::is.na(sf) && base::nzchar(sf)) {
      out_dir <- base::file.path(out_dir, sf)
    }
  }
  if (!base::dir.exists(out_dir)) {
    base::dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }
  out_dir
}

#' Cluster matrix rows with CH-based k selection
#' @noRd
.hc_cluster_with_scoring <- function(x,
                                     k,
                                     method = c("kmeans", "hclust"),
                                     nstart = 50,
                                     seed = 42) {
  method <- base::match.arg(method)
  x <- base::as.matrix(x)
  n <- base::nrow(x)

  k <- base::sort(base::unique(base::as.integer(k)))
  k <- k[is.finite(k) & k >= 2 & k < n]
  if (base::length(k) == 0) {
    stop("No valid `k` left. Need values in [2, n_rows-1].")
  }

  score_rows <- list()
  fit_list <- list()
  h_tree <- NULL
  if (identical(method, "hclust")) {
    h_tree <- stats::hclust(stats::dist(x), method = "ward.D2")
  }

  for (kk in k) {
    if (identical(method, "kmeans")) {
      fit <- .hc_with_seed(seed + kk, tryCatch(
        stats::kmeans(x, centers = kk, nstart = nstart, iter.max = 100),
        error = function(e) NULL
      ))

      # Robust fallback when kmeans cannot place centers (e.g. duplicate rows).
      if (is.null(fit)) {
        if (is.null(h_tree)) {
          h_tree <- stats::hclust(stats::dist(x), method = "ward.D2")
        }
        cl <- stats::cutree(h_tree, k = kk)
        tw <- NA_real_
        fit_list[[as.character(kk)]] <- cl
      } else {
        cl <- fit$cluster
        tw <- fit$tot.withinss
        fit_list[[as.character(kk)]] <- fit
      }
    } else {
      cl <- stats::cutree(h_tree, k = kk)
      tw <- NA_real_
      fit_list[[as.character(kk)]] <- cl
    }

    score_rows[[base::length(score_rows) + 1]] <- base::data.frame(
      k = kk,
      ch_index = .hc_ch_index(x, cl),
      tot_withinss = tw,
      stringsAsFactors = FALSE
    )
  }

  score_df <- base::do.call(base::rbind, score_rows)
  score_df$score_for_rank <- score_df$ch_index
  score_df$score_for_rank[!is.finite(score_df$score_for_rank)] <- -Inf
  score_df <- score_df[base::order(-score_df$score_for_rank, score_df$k), , drop = FALSE]
  score_df$score_for_rank <- NULL

  best_k <- score_df$k[[1]]
  if (!is.finite(score_df$ch_index[[1]])) {
    warning(
      "CH index is non-finite for all candidate k. Using smallest candidate k = ",
      best_k, "."
    )
  }

  if (identical(method, "kmeans")) {
    best_fit <- fit_list[[as.character(best_k)]]
    if (base::is.list(best_fit) && !is.null(best_fit$cluster)) {
      best_cluster <- best_fit$cluster
    } else {
      best_cluster <- base::as.integer(best_fit)
    }
  } else {
    best_fit <- fit_list[[as.character(best_k)]]
    best_cluster <- base::as.integer(best_fit)
  }

  list(
    score_table = score_df,
    best_k = best_k,
    cluster = base::as.integer(best_cluster),
    fit = best_fit
  )
}

#' Build module-cluster color palette (shades per module)
#' @noRd
.hc_module_cluster_palette <- function(module_lookup, module_best_k) {
  pal <- c()
  for (i in base::seq_len(base::nrow(module_best_k))) {
    mod <- as.character(module_best_k$module[[i]])
    kk <- as.integer(module_best_k$best_k[[i]])
    if (!is.finite(kk) || kk < 1) {
      next
    }
    base_col <- module_lookup$module_color[base::match(mod, module_lookup$module)]
    if (base::length(base_col) == 0 || is.na(base_col)) {
      base_col <- "#4E79A7"
    }
    shades <- grDevices::colorRampPalette(c("#F5F5F5", base_col))(kk + 1)
    shades <- shades[2:(kk + 1)]
    base::names(shades) <- base::paste0(mod, "__", base::seq_len(kk))
    pal <- base::c(pal, shades)
  }
  pal
}

#' Build donor x time matrix for one module
#' @noRd
.hc_module_time_matrix <- function(donor_time_module,
                                   donors,
                                   module_id,
                                   time_levels,
                                   traj_impute_mode = c("none", "linear", "locf"),
                                   na_impute = c("median", "zero")) {
  traj_impute_mode <- base::match.arg(traj_impute_mode)
  na_impute <- base::match.arg(na_impute)

  m_df <- donor_time_module[donor_time_module$module == module_id, , drop = FALSE]
  mat <- base::matrix(
    NA_real_,
    nrow = base::length(donors),
    ncol = base::length(time_levels),
    dimnames = list(donors, time_levels)
  )
  if (base::nrow(m_df) > 0) {
    row_idx <- base::match(base::as.character(m_df$donor), donors)
    col_idx <- base::match(base::as.character(m_df$time), time_levels)
    ok <- !base::is.na(row_idx) & !base::is.na(col_idx)
    mat[base::cbind(row_idx[ok], col_idx[ok])] <- as.numeric(m_df$value[ok])
  }

  for (rr in base::seq_len(base::nrow(mat))) {
    mat[rr, ] <- .hc_impute_series(mat[rr, ], method = traj_impute_mode)
  }

  if (identical(na_impute, "median")) {
    for (cc in base::seq_len(base::ncol(mat))) {
      v <- mat[, cc]
      if (base::any(base::is.na(v))) {
        med <- stats::median(v, na.rm = TRUE)
        if (!base::is.finite(med)) {
          med <- 0
        }
        v[base::is.na(v)] <- med
        mat[, cc] <- v
      }
    }
  } else {
    mat[base::is.na(mat)] <- 0
  }

  keep_cols <- base::apply(mat, 2, stats::sd, na.rm = TRUE) > 0
  mat_use <- mat[, keep_cols, drop = FALSE]
  if (base::ncol(mat_use) == 0) {
    mat_use <- mat
  }

  list(raw = mat, use = mat_use)
}

#' Cluster one matrix into k labels
#' @noRd
.hc_cluster_labels <- function(x,
                               k,
                               method = c("kmeans", "hclust"),
                               seed = 42,
                               nstart = 30) {
  method <- base::match.arg(method)
  if (k <= 1 || base::nrow(x) <= 1) {
    return(base::rep(1L, base::nrow(x)))
  }
  if (identical(method, "kmeans")) {
    fit <- .hc_with_seed(seed, tryCatch(
      stats::kmeans(x, centers = k, nstart = nstart, iter.max = 100),
      error = function(e) NULL
    ))
    if (!base::is.null(fit) && !base::is.null(fit$cluster)) {
      return(base::as.integer(fit$cluster))
    }
  }
  h_tree <- stats::hclust(stats::dist(x), method = "ward.D2")
  base::as.integer(stats::cutree(h_tree, k = k))
}

#' Compute cluster centroids
#' @noRd
.hc_cluster_centroids <- function(x, labels, k) {
  x <- base::as.matrix(x)
  labels <- base::as.integer(labels)
  out <- base::matrix(
    NA_real_,
    nrow = k,
    ncol = base::ncol(x),
    dimnames = list(base::as.character(base::seq_len(k)), base::colnames(x))
  )
  for (cl in base::seq_len(k)) {
    idx <- base::which(labels == cl)
    if (base::length(idx) > 0) {
      out[cl, ] <- base::colMeans(x[idx, , drop = FALSE], na.rm = TRUE)
    }
  }
  out
}

#' Build PCA axis labels with explained variance
#' @noRd
.hc_pca_axis_labels <- function(var_explained = NULL) {
  x_lab <- "PC1"
  y_lab <- "PC2"
  if (!is.null(var_explained)) {
    vv <- .hc_as_numeric_safely(var_explained)
    if (base::length(vv) >= 2 && base::all(base::is.finite(vv[base::seq_len(2L)]))) {
      x_lab <- base::sprintf("PC1 (%.1f%%)", vv[[1]])
      y_lab <- base::sprintf("PC2 (%.1f%%)", vv[[2]])
    }
  }
  list(x = x_lab, y = y_lab)
}

#' Align cluster labels to reference centroids
#' @noRd
.hc_align_labels_to_reference <- function(labels, x, ref_centroids) {
  labels <- base::as.integer(labels)
  k <- base::nrow(ref_centroids)
  if (k <= 1) {
    return(base::rep(1L, base::length(labels)))
  }

  cur_centroids <- .hc_cluster_centroids(x = x, labels = labels, k = k)
  if (base::any(!base::is.finite(ref_centroids)) || base::any(!base::is.finite(cur_centroids))) {
    return(labels)
  }

  perms <- gtools::permutations(n = k, r = k, v = base::seq_len(k))
  best_cost <- Inf
  best_perm <- base::seq_len(k)
  for (i in base::seq_len(base::nrow(perms))) {
    p <- perms[i, ]
    cst <- 0
    for (cl in base::seq_len(k)) {
      d <- cur_centroids[cl, ] - ref_centroids[p[[cl]], ]
      cst <- cst + base::sum(d * d)
    }
    if (base::is.finite(cst) && cst < best_cost) {
      best_cost <- cst
      best_perm <- p
    }
  }

  out <- labels
  for (cl in base::seq_len(k)) {
    out[labels == cl] <- base::as.integer(best_perm[[cl]])
  }
  base::as.integer(out)
}

#' Interpolate a module color by numeric proportion
#' @noRd
.hc_mix_with_white <- function(base_color, p) {
  p <- base::pmin(base::pmax(as.numeric(p), 0), 1)
  rgb_base <- grDevices::col2rgb(base_color) / 255
  rgb_out <- (1 - p) * 1 + p * rgb_base
  grDevices::rgb(rgb_out[1], rgb_out[2], rgb_out[3])
}

#' Build consensus clustering for a single k
#' @noRd
.hc_consensus_one_k <- function(x,
                                k,
                                method = c("kmeans", "hclust"),
                                runs = 200,
                                sample_fraction = 0.8,
                                feature_fraction = 0.8,
                                nstart = 50,
                                seed = 42,
                                linkage = "average") {
  method <- base::match.arg(method)
  x <- base::as.matrix(x)
  n <- base::nrow(x)
  p <- base::ncol(x)
  if (n < 3 || p < 1 || k < 2 || k >= n) {
    stop("Invalid dimensions/k for consensus clustering.")
  }

  runs <- base::max(1L, as.integer(runs))
  sample_fraction <- base::min(1, base::max(0.3, as.numeric(sample_fraction)))
  feature_fraction <- base::min(1, base::max(0.2, as.numeric(feature_fraction)))

  row_n <- base::min(n, base::max(k + 1L, as.integer(base::floor(n * sample_fraction))))
  col_n <- base::min(p, base::max(1L, as.integer(base::floor(p * feature_fraction))))

  co_count <- base::matrix(0, nrow = n, ncol = n)
  co_total <- base::matrix(0, nrow = n, ncol = n)

  for (rr in base::seq_len(runs)) {
    sample_idx <- .hc_with_seed(seed + rr, {
      list(
        rows = base::sort(base::sample.int(n = n, size = row_n, replace = FALSE)),
        cols = base::sort(base::sample.int(n = p, size = col_n, replace = FALSE))
      )
    })
    idx_rows <- sample_idx$rows
    idx_cols <- sample_idx$cols
    x_sub <- x[idx_rows, idx_cols, drop = FALSE]

    lbl <- .hc_cluster_labels(
      x = x_sub,
      k = k,
      method = method,
      seed = seed + rr * 31L,
      nstart = nstart
    )

    co_total[idx_rows, idx_rows] <- co_total[idx_rows, idx_rows] + 1

    for (cl in base::sort(base::unique(lbl))) {
      ids <- idx_rows[lbl == cl]
      if (base::length(ids) > 0) {
        co_count[ids, ids] <- co_count[ids, ids] + 1
      }
    }
  }

  cons <- base::matrix(0, nrow = n, ncol = n, dimnames = list(base::rownames(x), base::rownames(x)))
  obs <- co_total > 0
  cons[obs] <- co_count[obs] / co_total[obs]
  cons <- (cons + base::t(cons)) / 2
  base::diag(cons) <- 1

  d <- stats::as.dist(1 - cons)
  h <- stats::hclust(d, method = linkage)
  lbl_final <- base::as.integer(stats::cutree(h, k = k))

  up <- upper.tri(cons, diag = FALSE)
  same <- up & outer(lbl_final, lbl_final, FUN = "==")
  diff <- up & !outer(lbl_final, lbl_final, FUN = "==")

  within <- if (base::any(same)) base::mean(cons[same], na.rm = TRUE) else NA_real_
  between <- if (base::any(diff)) base::mean(cons[diff], na.rm = TRUE) else NA_real_
  delta <- within - between
  pac <- if (base::any(up)) base::mean(cons[up] > 0.1 & cons[up] < 0.9, na.rm = TRUE) else NA_real_
  observed_pair_fraction <- if (base::any(up)) base::mean(obs[up], na.rm = TRUE) else NA_real_

  stability <- base::vapply(base::seq_len(n), function(i) {
    idx <- base::which(lbl_final == lbl_final[[i]] & base::seq_len(n) != i)
    if (base::length(idx) == 0) {
      return(1)
    }
    base::mean(cons[i, idx], na.rm = TRUE)
  }, numeric(1))

  list(
    labels = lbl_final,
    consensus_matrix = cons,
    within = within,
    between = between,
    delta = delta,
    pac = pac,
    observed_pair_fraction = observed_pair_fraction,
    stability = stability
  )
}

#' Build longitudinal module means per donor and timepoint
#'
#' Computes module-level mean expression trajectories from the clustered modules
#' in an `HCoCenaExperiment`. This creates a donor-by-(module x timepoint)
#' feature matrix that can be used for longitudinal endotype discovery.
#'
#' @param hc A `HCoCenaExperiment`.
#' @param donor_col Annotation column holding donor/patient IDs.
#' @param time_col Annotation column holding ordered timepoint labels.
#' @param layer Layer to use (index, layer id, or layer name). Defaults to the
#'   first available layer.
#' @param group_col Optional annotation column (e.g. outcome group) to carry.
#' @param modules Optional subset of modules (module labels or color names).
#' @param use_module_labels Logical. If `TRUE`, use `M`-style labels from
#'   `hc_plot_cluster_heatmap()` when available; otherwise use color names.
#' @param time_levels Optional explicit timepoint order.
#' @param aggregate_fun Function used to aggregate technical replicates per
#'   donor/timepoint/module (default `mean`).
#' @param value_label Optional y-axis label base (e.g. `"Module mean GFC"`).
#'   If `NULL`, an automatic label is inferred from the selected layer id.
#' @param impute_missing One of `"none"`, `"linear"`, `"locf"` to impute missing
#'   donor-timepoint values within each donor/module trajectory.
#' @param min_genes_per_module Minimum number of overlapping genes required to
#'   keep a module.
#' @param slot_name Name of the satellite slot to store outputs.
#'
#' @return Updated `HCoCenaExperiment` with longitudinal data in
#'   `hc@satellite[[slot_name]]`.
#' @examples
#' hc <- hc_example_data("clustered")
#' hc <- hc_longitudinal_module_means(hc, donor_col = "donor", time_col = "timepoint")
#' @export
hc_longitudinal_module_means <- function(hc,
                                         donor_col = "donor",
                                         time_col = "time_point",
                                         layer = NULL,
                                         group_col = NULL,
                                         modules = NULL,
                                         use_module_labels = TRUE,
                                         time_levels = NULL,
                                         aggregate_fun = function(x) base::mean(x, na.rm = TRUE),
                                         value_label = NULL,
                                         impute_missing = c("linear", "locf", "none"),
                                         min_genes_per_module = 5,
                                         slot_name = "longitudinal_module_means") {
  if (!inherits(hc, "HCoCenaExperiment")) {
    stop("`hc` must be a `HCoCenaExperiment`.")
  }
  impute_missing <- base::match.arg(impute_missing)

  cluster_obj <- as.list(hc@integration@cluster)
  cluster_info <- cluster_obj[["cluster_information"]]
  if (is.null(cluster_info)) {
    stop("No cluster information found. Run `hc_cluster_calculation()` first.")
  }
  cluster_info <- base::as.data.frame(cluster_info, stringsAsFactors = FALSE)
  if ("cluster_included" %in% base::colnames(cluster_info)) {
    cluster_info <- dplyr::filter(cluster_info, cluster_included == "yes")
  }
  if (base::nrow(cluster_info) == 0) {
    stop("No included modules found in cluster information.")
  }

  layer_id <- .hc_resolve_layer_id(hc = hc, layer = layer)
  se <- MultiAssayExperiment::experiments(hc@mae)[[layer_id]]
  counts <- SummarizedExperiment::assay(se, "counts")
  counts <- base::as.matrix(counts)
  if (!is.numeric(counts)) {
    counts <- base::matrix(
      data = .hc_as_numeric_safely(counts),
      nrow = base::nrow(counts),
      ncol = base::ncol(counts),
      dimnames = base::dimnames(counts)
    )
  }
  if (is.null(base::rownames(counts)) || is.null(base::colnames(counts))) {
    stop("Counts matrix must have gene row names and sample column names.")
  }

  anno <- base::as.data.frame(SummarizedExperiment::colData(se), stringsAsFactors = FALSE)
  if (is.null(base::rownames(anno))) {
    base::rownames(anno) <- base::colnames(counts)
  }
  sample_ids <- base::colnames(counts)
  anno_idx <- base::match(sample_ids, base::rownames(anno))
  if (base::any(base::is.na(anno_idx))) {
    miss <- sample_ids[base::is.na(anno_idx)]
    stop(
      "Annotation rows missing for sample IDs: ",
      base::paste(utils::head(miss, 10), collapse = ", "),
      if (base::length(miss) > 10) " ..."
    )
  }
  anno <- anno[anno_idx, , drop = FALSE]

  req_cols <- c(donor_col, time_col)
  if (!base::is.null(group_col)) {
    req_cols <- base::c(req_cols, group_col)
  }
  missing_cols <- req_cols[!req_cols %in% base::colnames(anno)]
  if (base::length(missing_cols) > 0) {
    stop("Missing annotation columns in selected layer: ", base::paste(missing_cols, collapse = ", "), ".")
  }

  module_genes <- .hc_module_gene_map(cluster_info = cluster_info)
  module_colors <- .hc_module_colors_in_cluster_order(
    cluster_info = cluster_info,
    module_genes = module_genes
  )

  label_map <- cluster_obj[["module_label_map"]]
  label_map <- if (is.null(label_map)) NULL else base::unlist(label_map)
  label_map <- .hc_resolve_module_label_map_for_colors(
    label_map = label_map,
    module_colors = module_colors
  )

  module_labels <- module_colors
  if (isTRUE(use_module_labels) && !is.null(label_map) && base::length(label_map) > 0) {
    idx <- base::match(module_colors, base::names(label_map))
    ok <- !base::is.na(idx)
    module_labels[ok] <- label_map[idx[ok]]
  }
  if (isTRUE(use_module_labels)) {
    module_labels <- .hc_normalize_longitudinal_module_labels(
      module_labels = module_labels,
      module_colors = module_colors
    )
    natural_labels <- .hc_natural_module_order(module_labels)
    ord <- base::match(natural_labels, module_labels)
    ord <- ord[!base::is.na(ord)]
    module_labels <- module_labels[ord]
    module_colors <- module_colors[ord]
  }

  module_lookup <- base::data.frame(
    module_color = module_colors,
    module = module_labels,
    stringsAsFactors = FALSE
  )

  if (!base::is.null(modules)) {
    modules <- base::as.character(modules)
    keep <- module_lookup$module_color %in% modules | module_lookup$module %in% modules
    module_lookup <- module_lookup[keep, , drop = FALSE]
    module_colors <- module_lookup$module_color
    module_labels <- module_lookup$module
    if (base::length(module_colors) == 0) {
      stop("None of the requested `modules` were found.")
    }
  }

  # Compute module means per sample.
  module_sample <- base::matrix(
    NA_real_,
    nrow = base::length(module_colors),
    ncol = base::length(sample_ids),
    dimnames = list(module_labels, sample_ids)
  )
  gene_overlap_n <- base::integer(base::length(module_colors))

  for (i in base::seq_along(module_colors)) {
    mcol <- module_colors[[i]]
    genes <- module_genes[[mcol]]
    genes <- base::intersect(genes, base::rownames(counts))
    gene_overlap_n[[i]] <- base::length(genes)
    if (base::length(genes) < min_genes_per_module) {
      next
    }
    vals <- counts[genes, , drop = FALSE]
    module_sample[i, ] <- base::colMeans(vals, na.rm = TRUE)
  }

  module_lookup$genes_in_counts <- gene_overlap_n
  keep_rows <- gene_overlap_n >= min_genes_per_module
  if (!base::any(keep_rows)) {
    stop("No module passed `min_genes_per_module = ", min_genes_per_module, "`.")
  }
  module_lookup <- module_lookup[keep_rows, , drop = FALSE]
  module_sample <- module_sample[keep_rows, , drop = FALSE]

  donor_vec <- as.character(anno[[donor_col]])
  time_vec <- as.character(anno[[time_col]])
  group_vec <- if (!base::is.null(group_col)) as.character(anno[[group_col]]) else NULL
  time_vec <- base::trimws(time_vec)
  observed_time_levels <- base::sort(base::unique(time_vec[!base::is.na(time_vec) & base::nzchar(time_vec)]))

  if (base::is.null(time_levels)) {
    if (base::is.factor(anno[[time_col]])) {
      time_levels <- base::levels(anno[[time_col]])
      time_levels <- time_levels[time_levels %in% base::unique(time_vec)]
    } else {
      time_levels <- base::sort(base::unique(time_vec))
    }
  } else {
    time_levels <- base::trimws(base::as.character(time_levels))
    time_levels <- time_levels[!base::is.na(time_levels) & base::nzchar(time_levels)]
  }

  if (base::length(time_levels) == 0) {
    stop("`time_levels` is empty after normalization.")
  }
  overlap_time_levels <- base::intersect(time_levels, observed_time_levels)
  if (base::length(overlap_time_levels) == 0) {
    stop(
      "`time_levels` does not match the observed values in `", time_col, "`.\n",
      "Provided: ", base::paste(time_levels, collapse = ", "), "\n",
      "Observed: ", base::paste(observed_time_levels, collapse = ", "),
      call. = FALSE
    )
  }

  if (base::is.null(value_label) || !base::nzchar(as.character(value_label[[1]]))) {
    if (grepl("gfc", layer_id, ignore.case = TRUE)) {
      value_label <- "Module mean GFC"
    } else if (grepl("fc", layer_id, ignore.case = TRUE)) {
      value_label <- "Module mean FC"
    } else {
      value_label <- "Module mean expression"
    }
  } else {
    value_label <- as.character(value_label[[1]])
  }

  long_list <- base::lapply(base::seq_len(base::nrow(module_sample)), function(i) {
    base::data.frame(
      sample_id = sample_ids,
      donor = donor_vec,
      time = time_vec,
      module = base::rownames(module_sample)[[i]],
      value = as.numeric(module_sample[i, ]),
      stringsAsFactors = FALSE
    )
  })
  long_df <- base::do.call(base::rbind, long_list)
  if (!base::is.null(group_col)) {
    long_df$group <- base::rep(group_vec, times = base::nrow(module_sample))
  }

  # Aggregate per donor/time/module.
  if (!base::is.null(group_col)) {
    agg <- stats::aggregate(
      value ~ donor + time + module + group,
      data = long_df,
      FUN = aggregate_fun
    )
  } else {
    agg <- stats::aggregate(
      value ~ donor + time + module,
      data = long_df,
      FUN = aggregate_fun
    )
  }

  donors <- base::sort(base::unique(agg$donor))
  modules_out <- base::unique(module_lookup$module)
  full_grid <- base::expand.grid(
    donor = donors,
    module = modules_out,
    time = time_levels,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  merge_cols <- c("donor", "module", "time")
  if (!base::is.null(group_col) && "group" %in% base::colnames(agg)) {
    donor_group <- stats::aggregate(group ~ donor, data = agg, FUN = function(z) z[[1]])
    full_grid <- base::merge(full_grid, donor_group, by = "donor", all.x = TRUE, sort = FALSE)
  }

  full_df <- base::merge(full_grid, agg, by = merge_cols, all.x = TRUE, sort = FALSE)
  full_df$time <- base::factor(full_df$time, levels = time_levels, ordered = TRUE)
  full_df <- full_df[base::order(full_df$donor, full_df$module, full_df$time), , drop = FALSE]

  if (!identical(impute_missing, "none")) {
    split_idx <- base::split(
      base::seq_len(base::nrow(full_df)),
      base::paste(full_df$donor, full_df$module, sep = "||")
    )
    for (idx in split_idx) {
      full_df$value[idx] <- .hc_impute_series(full_df$value[idx], method = impute_missing)
    }
  }

  feature_order <- base::as.vector(
    base::unlist(
      base::lapply(modules_out, function(m) base::paste0(m, "__", time_levels)),
      use.names = FALSE
    )
  )
  donor_feature <- base::matrix(
    NA_real_,
    nrow = base::length(donors),
    ncol = base::length(feature_order),
    dimnames = list(donors, feature_order)
  )

  full_df$feature <- base::paste0(full_df$module, "__", base::as.character(full_df$time))
  for (i in base::seq_len(base::nrow(full_df))) {
    d <- full_df$donor[[i]]
    f <- full_df$feature[[i]]
    if (d %in% base::rownames(donor_feature) && f %in% base::colnames(donor_feature)) {
      donor_feature[d, f] <- full_df$value[[i]]
    }
  }

  sat <- as.list(hc@satellite)
  sat[[slot_name]] <- list(
    created_at = as.character(base::Sys.time()),
    layer_id = layer_id,
    donor_col = donor_col,
    time_col = time_col,
    group_col = group_col,
    time_levels = time_levels,
    module_lookup = module_lookup,
    module_sample_means = module_sample,
    donor_time_module = full_df,
    donor_feature_matrix = donor_feature,
    value_label = value_label,
    impute_missing = impute_missing,
    min_genes_per_module = as.integer(min_genes_per_module)
  )
  hc@satellite <- S4Vectors::SimpleList(sat)
  methods::validObject(hc)
  hc
}

#' Cluster donors into longitudinal endotypes
#'
#' Clusters donors based on module-time features produced by
#' [hc_longitudinal_module_means()]. Candidate cluster numbers are scored via the
#' Calinski-Harabasz index.
#'
#' @param hc A `HCoCenaExperiment`.
#' @param input_slot Satellite slot produced by `hc_longitudinal_module_means()`.
#' @param output_slot Satellite slot to store clustering results.
#' @param k Integer vector of candidate cluster counts.
#' @param method One of `"kmeans"` or `"hclust"`.
#' @param scale_features Logical. If `TRUE`, z-score features before clustering.
#' @param nstart Number of random starts for `kmeans`.
#' @param seed Random seed.
#' @param na_impute One of `"median"` or `"zero"` for feature NA imputation.
#' @param endotype_prefix Prefix for final endotype labels.
#' @param module_k Candidate cluster counts for module-wise trajectory clustering.
#' @param module_method One of `"kmeans"` or `"hclust"` for module-wise
#'   clustering.
#' @param module_scale_features Logical. If `TRUE`, z-score module trajectory
#'   features before module-wise clustering.
#' @param module_nstart Number of random starts for module-wise `kmeans`.
#' @param module_seed Random seed for module-wise clustering.
#' @param module_na_impute One of `"median"` or `"zero"` for module-wise NA
#'   imputation.
#'
#' @return Updated `HCoCenaExperiment` with endotype clustering in
#'   `hc@satellite[[output_slot]]`.
#' @examples
#' hc <- hc_example_data("clustered")
#' hc <- hc_longitudinal_module_means(hc, donor_col = "donor", time_col = "timepoint")
#' hc <- hc_longitudinal_endotype_clustering(hc, k = 2)
#' @export
hc_longitudinal_endotype_clustering <- function(hc,
                                                input_slot = "longitudinal_module_means",
                                                output_slot = "longitudinal_endotypes",
                                                k = 2:6,
                                                method = c("kmeans", "hclust"),
                                                scale_features = TRUE,
                                                nstart = 50,
                                                seed = 42,
                                                na_impute = c("median", "zero"),
                                                endotype_prefix = "E",
                                                module_k = NULL,
                                                module_method = c("kmeans", "hclust"),
                                                module_scale_features = TRUE,
                                                module_nstart = 30,
                                                module_seed = 42,
                                                module_na_impute = c("median", "zero")) {
  if (!inherits(hc, "HCoCenaExperiment")) {
    stop("`hc` must be a `HCoCenaExperiment`.")
  }
  method <- base::match.arg(method)
  na_impute <- base::match.arg(na_impute)
  module_method <- base::match.arg(module_method)
  module_na_impute <- base::match.arg(module_na_impute)
  if (base::is.null(module_k)) {
    module_k <- k
  }

  sat <- as.list(hc@satellite)
  inp <- sat[[input_slot]]
  if (is.null(inp) || is.null(inp$donor_feature_matrix)) {
    stop("Input slot `", input_slot, "` is missing or does not contain `donor_feature_matrix`.")
  }

  x <- base::as.matrix(inp$donor_feature_matrix)
  if (base::nrow(x) < 3) {
    stop("Need at least 3 donors for endotype clustering.")
  }

  # Drop unusable features.
  keep_cols <- base::colSums(!base::is.na(x)) > 0
  x <- x[, keep_cols, drop = FALSE]
  if (base::ncol(x) == 0) {
    stop("No usable features left after removing all-NA columns.")
  }

  if (identical(na_impute, "median")) {
    for (j in base::seq_len(base::ncol(x))) {
      v <- x[, j]
      if (base::any(base::is.na(v))) {
        med <- stats::median(v, na.rm = TRUE)
        if (!is.finite(med)) {
          med <- 0
        }
        v[base::is.na(v)] <- med
        x[, j] <- v
      }
    }
  } else {
    x[base::is.na(x)] <- 0
  }

  # Drop zero-variance features.
  keep_var <- base::apply(x, 2, stats::sd, na.rm = TRUE) > 0
  x <- x[, keep_var, drop = FALSE]
  if (base::ncol(x) == 0) {
    stop("No variable features available for clustering.")
  }

  x_use <- base::as.matrix(x)
  if (isTRUE(scale_features)) {
    x_use <- base::as.matrix(base::scale(x_use))
    x_use[!is.finite(x_use)] <- 0
  }

  global_fit <- .hc_cluster_with_scoring(
    x = x_use,
    k = k,
    method = method,
    nstart = nstart,
    seed = seed
  )

  score_df <- global_fit$score_table
  best_k <- global_fit$best_k
  final_cluster <- global_fit$cluster
  endotype <- base::paste0(endotype_prefix, final_cluster)
  donor_ids <- base::rownames(x_use)

  donor_cluster <- base::data.frame(
    donor = donor_ids,
    cluster = as.integer(final_cluster),
    endotype = endotype,
    stringsAsFactors = FALSE
  )

  pca_fit <- stats::prcomp(x_use, center = TRUE, scale. = FALSE)
  pca_var <- 100 * (pca_fit$sdev^2 / sum(pca_fit$sdev^2))
  pc1 <- pca_fit$x[, 1]
  pc2 <- if (base::ncol(pca_fit$x) >= 2) pca_fit$x[, 2] else base::rep(0, base::nrow(pca_fit$x))
  pca_df <- base::data.frame(
    donor = donor_ids,
    PC1 = pc1,
    PC2 = pc2,
    endotype = endotype,
    stringsAsFactors = FALSE
  )

  donor_time_module <- inp$donor_time_module
  donor_time_module$endotype <- donor_cluster$endotype[base::match(donor_time_module$donor, donor_cluster$donor)]
  donor_time_module <- donor_time_module[!base::is.na(donor_time_module$endotype), , drop = FALSE]

  traj_mean <- stats::aggregate(
    value ~ endotype + module + time,
    data = donor_time_module,
    FUN = base::mean
  )

  # Module-wise trajectory clustering (each module may have a different best k).
  module_lookup <- inp$module_lookup
  if (is.null(module_lookup) || base::nrow(module_lookup) == 0) {
    module_lookup <- base::data.frame(
      module = base::sort(base::unique(donor_time_module$module)),
      module_color = "grey60",
      stringsAsFactors = FALSE
    )
  }
  module_lookup <- .hc_reorder_module_lookup_natural(module_lookup)
  if (!"module_color" %in% base::colnames(module_lookup)) {
    module_lookup$module_color <- "grey60"
  }

  time_levels <- inp$time_levels
  if (is.null(time_levels) || base::length(time_levels) == 0) {
    time_levels <- base::sort(base::unique(base::as.character(donor_time_module$time)))
  } else {
    time_levels <- base::as.character(time_levels)
  }

  module_ids <- base::unique(base::as.character(module_lookup$module))
  module_score_rows <- list()
  module_best_rows <- list()
  module_assignment_rows <- list()

  for (m in module_ids) {
    m_df <- donor_time_module[donor_time_module$module == m, , drop = FALSE]
    if (base::nrow(m_df) == 0) {
      next
    }

    mat <- base::matrix(
      NA_real_,
      nrow = base::length(donor_ids),
      ncol = base::length(time_levels),
      dimnames = list(donor_ids, time_levels)
    )

    row_idx <- base::match(base::as.character(m_df$donor), donor_ids)
    col_idx <- base::match(base::as.character(m_df$time), time_levels)
    ok <- !base::is.na(row_idx) & !base::is.na(col_idx)
    mat[base::cbind(row_idx[ok], col_idx[ok])] <- as.numeric(m_df$value[ok])

    traj_impute_mode <- inp$impute_missing
    if (is.null(traj_impute_mode) || !traj_impute_mode %in% c("none", "linear", "locf")) {
      traj_impute_mode <- "none"
    }
    for (rr in base::seq_len(base::nrow(mat))) {
      mat[rr, ] <- .hc_impute_series(mat[rr, ], method = traj_impute_mode)
    }

    if (identical(module_na_impute, "median")) {
      for (cc in base::seq_len(base::ncol(mat))) {
        v <- mat[, cc]
        if (base::any(base::is.na(v))) {
          med <- stats::median(v, na.rm = TRUE)
          if (!is.finite(med)) {
            med <- 0
          }
          v[base::is.na(v)] <- med
          mat[, cc] <- v
        }
      }
    } else {
      mat[base::is.na(mat)] <- 0
    }

    keep_cols_mod <- base::apply(mat, 2, stats::sd, na.rm = TRUE) > 0
    mat_use <- mat[, keep_cols_mod, drop = FALSE]
    if (base::ncol(mat_use) == 0) {
      mat_use <- mat
    }

    if (isTRUE(module_scale_features)) {
      mat_use <- base::as.matrix(base::scale(mat_use))
      mat_use[!is.finite(mat_use)] <- 0
    }

    # if rows are completely identical, set one cluster.
    valid_module_k <- base::sort(base::unique(base::as.integer(module_k)))
    valid_module_k <- valid_module_k[base::is.finite(valid_module_k) & valid_module_k >= 2 & valid_module_k < base::nrow(mat_use)]

    if (base::nrow(mat_use) <= 2 ||
      base::all(base::apply(mat_use, 2, function(v) stats::sd(v, na.rm = TRUE) == 0)) ||
      base::length(valid_module_k) == 0) {
      mod_cluster <- base::rep(1L, base::nrow(mat_use))
      mod_score <- base::data.frame(
        k = 1L,
        ch_index = NA_real_,
        tot_withinss = NA_real_,
        stringsAsFactors = FALSE
      )
      mod_best_k <- 1L
    } else {
      mod_fit <- .hc_cluster_with_scoring(
        x = mat_use,
        k = valid_module_k,
        method = module_method,
        nstart = module_nstart,
        seed = module_seed
      )
      mod_cluster <- base::as.integer(mod_fit$cluster)
      mod_score <- mod_fit$score_table
      mod_best_k <- as.integer(mod_fit$best_k)
    }

    mod_score$module <- m
    module_score_rows[[base::length(module_score_rows) + 1]] <- mod_score
    module_best_rows[[base::length(module_best_rows) + 1]] <- base::data.frame(
      module = m,
      best_k = mod_best_k,
      stringsAsFactors = FALSE
    )
    module_assignment_rows[[base::length(module_assignment_rows) + 1]] <- base::data.frame(
      donor = donor_ids,
      module = m,
      cluster = mod_cluster,
      stringsAsFactors = FALSE
    )
  }

  if (base::length(module_assignment_rows) > 0) {
    module_cluster_assign <- base::do.call(base::rbind, module_assignment_rows)
  } else {
    module_cluster_assign <- base::data.frame(
      donor = character(0),
      module = character(0),
      cluster = integer(0),
      stringsAsFactors = FALSE
    )
  }
  if (base::length(module_score_rows) > 0) {
    module_score_table <- base::do.call(base::rbind, module_score_rows)
  } else {
    module_score_table <- base::data.frame(
      k = integer(0),
      ch_index = numeric(0),
      tot_withinss = numeric(0),
      module = character(0),
      stringsAsFactors = FALSE
    )
  }
  if (base::length(module_best_rows) > 0) {
    module_best_k <- base::do.call(base::rbind, module_best_rows)
  } else {
    module_best_k <- base::data.frame(module = character(0), best_k = integer(0), stringsAsFactors = FALSE)
  }

  module_cluster_assign$module_cluster_key <- base::paste0(
    module_cluster_assign$module, "__", module_cluster_assign$cluster
  )

  module_cluster_wide <- base::matrix(
    NA_integer_,
    nrow = base::length(donor_ids),
    ncol = base::length(module_ids),
    dimnames = list(donor_ids, module_ids)
  )
  if (base::nrow(module_cluster_assign) > 0) {
    for (i in base::seq_len(base::nrow(module_cluster_assign))) {
      d <- module_cluster_assign$donor[[i]]
      m <- module_cluster_assign$module[[i]]
      module_cluster_wide[d, m] <- as.integer(module_cluster_assign$cluster[[i]])
    }
  }

  donor_time_module_w_clusters <- base::merge(
    donor_time_module,
    module_cluster_assign[, c("donor", "module", "cluster", "module_cluster_key"), drop = FALSE],
    by = c("donor", "module"),
    all.x = TRUE,
    sort = FALSE
  )

  module_cluster_traj_mean <- stats::aggregate(
    value ~ module + time + cluster + module_cluster_key,
    data = donor_time_module_w_clusters,
    FUN = base::mean
  )

  module_cluster_palette <- .hc_module_cluster_palette(
    module_lookup = module_lookup,
    module_best_k = module_best_k
  )

  sat[[output_slot]] <- list(
    created_at = as.character(base::Sys.time()),
    source_slot = input_slot,
    method = method,
    candidate_k = k,
    score_table = score_df,
    best_k = best_k,
    donor_cluster = donor_cluster,
    pca = pca_df,
    pca_variance = pca_var,
    feature_matrix_used = x_use,
    donor_time_module = donor_time_module,
    trajectory_mean = traj_mean,
    time_levels = time_levels,
    module_method = module_method,
    module_k_candidates = module_k,
    module_cluster_score_table = module_score_table,
    module_cluster_best_k = module_best_k,
    module_cluster_assignments = module_cluster_assign,
    module_cluster_matrix = module_cluster_wide,
    donor_time_module_with_module_clusters = donor_time_module_w_clusters,
    module_cluster_trajectory_mean = module_cluster_traj_mean,
    module_cluster_palette = module_cluster_palette,
    module_lookup = module_lookup,
    value_label = inp$value_label
  )

  hc@satellite <- S4Vectors::SimpleList(sat)
  methods::validObject(hc)
  hc
}

#' Plot donor-level module trajectories (faint donors + mean line)
#'
#' Uses outputs from [hc_longitudinal_module_means()] to visualize donor
#' trajectories per module. Individual donors are shown in a light stroke, with
#' a highlighted module mean trajectory on top.
#'
#' @param hc A `HCoCenaExperiment`.
#' @param slot_name Satellite slot containing results from
#'   `hc_longitudinal_module_means()`.
#' @param save_pdf Logical. If `TRUE`, save PDF output.
#' @param file_prefix Prefix used for saved PDF file names.
#' @param donor_alpha Alpha for individual donor trajectories.
#' @param donor_linewidth Line width for individual donor trajectories.
#' @param mean_linewidth Line width for mean trajectories.
#' @param mean_point_size Point size for mean trajectory markers.
#' @param facet_ncol Number of facet columns.
#' @param free_y Logical. If `TRUE`, use free y-axis scales per module.
#' @param square_panels Logical. If `TRUE`, force square panel aspect ratio.
#' @param save_width,save_height Optional PDF width/height in inches. If
#'   `NULL`, defaults are chosen automatically.
#'
#' @return A list with one ggplot object (`module_means`).
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
#' p <- hc_plot_longitudinal_module_means(hc, save_pdf = FALSE)
#' @export
hc_plot_longitudinal_module_means <- function(hc,
                                              slot_name = "longitudinal_module_means",
                                              save_pdf = TRUE,
                                              file_prefix = "Longitudinal_ModuleMeans",
                                              donor_alpha = 0.16,
                                              donor_linewidth = 0.35,
                                              mean_linewidth = 1.1,
                                              mean_point_size = 1.8,
                                              facet_ncol = 4,
                                              free_y = TRUE,
                                              square_panels = FALSE,
                                              save_width = NULL,
                                              save_height = NULL) {
  if (!inherits(hc, "HCoCenaExperiment")) {
    stop("`hc` must be a `HCoCenaExperiment`.")
  }

  sat <- as.list(hc@satellite)
  obj <- sat[[slot_name]]
  if (is.null(obj) || is.null(obj$donor_time_module) || is.null(obj$module_lookup)) {
    stop("Slot `", slot_name, "` does not contain module trajectory data.")
  }

  df <- base::as.data.frame(obj$donor_time_module, stringsAsFactors = FALSE)
  module_lookup <- base::as.data.frame(obj$module_lookup, stringsAsFactors = FALSE)
  module_lookup <- .hc_reorder_module_lookup_natural(module_lookup)
  module_levels <- base::as.character(module_lookup$module)
  mod_col <- .hc_module_color_map(
    module_lookup = module_lookup,
    module_levels = module_levels
  )

  if (!is.null(obj$time_levels)) {
    df$time <- base::factor(base::as.character(df$time), levels = base::as.character(obj$time_levels), ordered = TRUE)
  }
  df$module <- base::factor(base::as.character(df$module), levels = module_levels)
  df$donor <- base::as.character(df$donor)

  mean_df <- stats::aggregate(value ~ module + time, data = df, FUN = base::mean)
  y_lab <- if (!is.null(obj$value_label) && base::nzchar(as.character(obj$value_label[[1]]))) {
    as.character(obj$value_label[[1]])
  } else {
    "Module mean expression"
  }

  p_mod <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = time, y = value, group = donor, color = module)
  ) +
    ggplot2::geom_line(alpha = donor_alpha, linewidth = donor_linewidth, show.legend = FALSE) +
    ggplot2::geom_line(
      data = mean_df,
      mapping = ggplot2::aes(group = module),
      linewidth = mean_linewidth,
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      data = mean_df,
      mapping = ggplot2::aes(group = module),
      size = mean_point_size,
      show.legend = FALSE
    ) +
    .hc_module_facet_wrap(
      module_levels = module_levels,
      mod_col = mod_col,
      ncol = facet_ncol,
      scales = if (isTRUE(free_y)) "free_y" else "fixed",
      facets = ~module
    ) +
    ggplot2::scale_color_manual(values = mod_col, drop = FALSE) +
    .hc_theme_pub(base_size = 11) +
    ggplot2::theme(
      aspect.ratio = if (isTRUE(square_panels)) 1 else NULL,
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    ) +
    ggplot2::labs(
      title = "Module trajectories across time (donors + module mean)",
      x = "Timepoint",
      y = y_lab
    )
  p_mod <- .hc_colorize_module_strips(
    p = p_mod,
    module_levels = module_levels,
    mod_col = mod_col
  )

  if (isTRUE(save_pdf)) {
    out_dir <- .hc_resolve_output_dir(hc)
    w <- if (is.null(save_width)) {
      if (isTRUE(square_panels)) 10 else 12
    } else {
      as.numeric(save_width[[1]])
    }
    h <- if (is.null(save_height)) {
      if (isTRUE(square_panels)) w else 8
    } else {
      as.numeric(save_height[[1]])
    }
    .hc_ggsave_pdf_png(
      filename = base::file.path(out_dir, base::paste0(file_prefix, ".pdf")),
      plot = p_mod,
      width = w,
      height = h,
      units = "in",
      device = grDevices::cairo_pdf
    )
  }

  list(module_means = p_mod)
}

#' Plot module-wise trajectory clusters and cluster heatmap
#'
#' Uses module-wise clustering results from [hc_longitudinal_endotype_clustering()]
#' and produces:
#' \enumerate{
#'   \item donor trajectories with module-cluster means per module
#'   \item a module-by-donor cluster heatmap
#' }
#'
#' @param hc A `HCoCenaExperiment`.
#' @param slot_name Satellite slot containing results from
#'   `hc_longitudinal_endotype_clustering()`.
#' @param save_pdf Logical. If `TRUE`, save PDF outputs.
#' @param file_prefix Prefix used for saved PDF file names.
#' @param donor_alpha Alpha for donor trajectories.
#' @param donor_linewidth Line width for donor trajectories.
#' @param mean_linewidth Line width for cluster mean trajectories.
#' @param mean_point_size Point size for cluster mean markers.
#' @param facet_ncol Number of facet columns.
#' @param free_y Logical. If `TRUE`, use free y-axis scales per module.
#' @param show_heatmap_numbers Logical. If `TRUE`, overlay cluster IDs as text.
#' @param square_panels Logical. If `TRUE`, force square panel aspect ratio for
#'   trajectory facets.
#' @param save_waves_width,save_waves_height Optional PDF size for trajectory
#'   plot in inches.
#' @param save_heatmap_width,save_heatmap_height Optional PDF size for CAP
#'   heatmap in inches.
#'
#' @return A list with ggplot objects: `module_cluster_waves` and
#'   `module_cluster_heatmap`.
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
#' p <- hc_plot_longitudinal_module_clusters(hc, save_pdf = FALSE)
#' @export
hc_plot_longitudinal_module_clusters <- function(hc,
                                                 slot_name = "longitudinal_endotypes",
                                                 save_pdf = TRUE,
                                                 file_prefix = "Longitudinal_ModuleClusters",
                                                 donor_alpha = 0.14,
                                                 donor_linewidth = 0.28,
                                                 mean_linewidth = 1.15,
                                                 mean_point_size = 1.9,
                                                 facet_ncol = 4,
                                                 free_y = TRUE,
                                                 show_heatmap_numbers = TRUE,
                                                 square_panels = FALSE,
                                                 save_waves_width = NULL,
                                                 save_waves_height = NULL,
                                                 save_heatmap_width = NULL,
                                                 save_heatmap_height = NULL) {
  if (!inherits(hc, "HCoCenaExperiment")) {
    stop("`hc` must be a `HCoCenaExperiment`.")
  }

  sat <- as.list(hc@satellite)
  obj <- sat[[slot_name]]
  req <- c(
    "module_cluster_assignments",
    "donor_time_module_with_module_clusters",
    "module_cluster_trajectory_mean",
    "module_cluster_palette",
    "module_cluster_matrix",
    "module_lookup"
  )
  if (is.null(obj) || !base::all(req %in% base::names(obj))) {
    stop("Slot `", slot_name, "` does not contain module-wise clustering outputs.")
  }

  pal <- obj$module_cluster_palette
  if (is.null(pal) || base::length(pal) == 0) {
    stop("No module-cluster palette found in slot `", slot_name, "`.")
  }

  module_lookup <- base::as.data.frame(obj$module_lookup, stringsAsFactors = FALSE)
  module_lookup <- .hc_reorder_module_lookup_natural(module_lookup)
  mod_levels <- base::as.character(module_lookup$module)
  mod_col <- .hc_module_color_map(
    module_lookup = module_lookup,
    module_levels = mod_levels
  )

  t_levels <- if (!is.null(obj$time_levels)) {
    base::as.character(obj$time_levels)
  } else {
    NULL
  }

  traj_df <- base::as.data.frame(obj$donor_time_module_with_module_clusters, stringsAsFactors = FALSE)
  mean_df <- base::as.data.frame(obj$module_cluster_trajectory_mean, stringsAsFactors = FALSE)
  y_lab <- if (!is.null(obj$value_label) && base::nzchar(as.character(obj$value_label[[1]]))) {
    as.character(obj$value_label[[1]])
  } else {
    "Module mean expression"
  }
  if (!is.null(t_levels)) {
    traj_df$time <- base::factor(base::as.character(traj_df$time), levels = t_levels, ordered = TRUE)
    mean_df$time <- base::factor(base::as.character(mean_df$time), levels = t_levels, ordered = TRUE)
  }

  traj_df <- traj_df[!base::is.na(traj_df$module_cluster_key), , drop = FALSE]
  traj_df$module <- base::factor(base::as.character(traj_df$module), levels = mod_levels)
  mean_df$module <- base::factor(base::as.character(mean_df$module), levels = mod_levels)

  p_waves <- ggplot2::ggplot(
    traj_df,
    ggplot2::aes(
      x = time,
      y = value,
      group = donor,
      color = module_cluster_key
    )
  ) +
    ggplot2::geom_line(alpha = donor_alpha, linewidth = donor_linewidth, show.legend = FALSE) +
    ggplot2::geom_line(
      data = mean_df,
      mapping = ggplot2::aes(group = module_cluster_key),
      linewidth = mean_linewidth,
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      data = mean_df,
      mapping = ggplot2::aes(group = module_cluster_key),
      size = mean_point_size,
      show.legend = FALSE
    ) +
    .hc_module_facet_wrap(
      module_levels = mod_levels,
      mod_col = mod_col,
      ncol = facet_ncol,
      scales = if (isTRUE(free_y)) "free_y" else "fixed",
      facets = ~module
    ) +
    ggplot2::scale_color_manual(values = pal, drop = FALSE) +
    .hc_theme_pub(base_size = 11) +
    ggplot2::theme(
      aspect.ratio = if (isTRUE(square_panels)) 1 else NULL,
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    ) +
    ggplot2::labs(
      title = "Module-wise trajectory clusters",
      x = "Timepoint",
      y = y_lab
    )
  p_waves <- .hc_colorize_module_strips(
    p = p_waves,
    module_levels = mod_levels,
    mod_col = mod_col
  )

  cl_mat <- base::as.matrix(obj$module_cluster_matrix)
  cl_mat_num <- cl_mat
  for (cc in base::seq_len(base::ncol(cl_mat_num))) {
    v <- as.numeric(cl_mat_num[, cc])
    if (base::any(base::is.na(v))) {
      med <- stats::median(v, na.rm = TRUE)
      if (!is.finite(med)) {
        med <- 0
      }
      v[base::is.na(v)] <- med
      cl_mat_num[, cc] <- v
    }
  }
  donor_order <- base::rownames(cl_mat_num)
  if (!base::is.null(obj$master_donors)) {
    stored_order <- base::as.character(obj$master_donors)
    stored_order <- stored_order[stored_order %in% base::rownames(cl_mat_num)]
    donor_order <- base::c(stored_order, base::setdiff(base::rownames(cl_mat_num), stored_order))
  } else if (base::nrow(cl_mat_num) > 1) {
    donor_order <- base::rownames(cl_mat_num)[stats::hclust(stats::dist(cl_mat_num), method = "ward.D2")$order]
  }
  module_order <- base::colnames(cl_mat)
  if (base::length(mod_levels) > 0) {
    module_order <- mod_levels[mod_levels %in% module_order]
  }

  hm_df <- base::as.data.frame(base::as.table(cl_mat), stringsAsFactors = FALSE)
  base::colnames(hm_df) <- c("donor", "module", "cluster")
  hm_df$cluster <- .hc_as_integer_safely(hm_df$cluster)
  hm_df$module_cluster_key <- ifelse(
    base::is.na(hm_df$cluster),
    NA_character_,
    base::paste0(hm_df$module, "__", hm_df$cluster)
  )
  hm_fill_cols <- pal[base::as.character(hm_df$module_cluster_key)]
  hm_fill_cols[base::is.na(hm_fill_cols)] <- "#F2F2F2"
  hm_rgb <- grDevices::col2rgb(hm_fill_cols) / 255
  hm_lum <- 0.2126 * hm_rgb[1, ] + 0.7152 * hm_rgb[2, ] + 0.0722 * hm_rgb[3, ]
  hm_df$text_color <- ifelse(hm_lum < 0.55, "white", "black")
  hm_df$donor <- base::factor(base::as.character(hm_df$donor), levels = donor_order)
  hm_df$module <- base::factor(base::as.character(hm_df$module), levels = module_order)

  p_heat <- ggplot2::ggplot(
    hm_df,
    ggplot2::aes(x = donor, y = module, fill = module_cluster_key)
  ) +
    ggplot2::geom_tile(color = "grey95", linewidth = 0.25) +
    ggplot2::scale_fill_manual(values = pal, na.value = "grey90", guide = "none") +
    .hc_theme_pub(base_size = 10) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),
      axis.text.y = ggplot2::element_text(face = "bold"),
      axis.title = ggplot2::element_text(face = "bold")
    ) +
    ggplot2::labs(
      title = "Module cluster membership heatmap",
      x = "Donor",
      y = "Module"
    )

  if (isTRUE(show_heatmap_numbers)) {
    p_heat <- p_heat +
      ggplot2::geom_text(
        ggplot2::aes(
          label = ifelse(base::is.na(cluster), "", cluster),
          color = text_color
        ),
        size = 2.4,
        show.legend = FALSE
      ) +
      ggplot2::scale_color_identity()
  }

  if (isTRUE(save_pdf)) {
    out_dir <- .hc_resolve_output_dir(hc)
    w_waves <- if (is.null(save_waves_width)) {
      if (isTRUE(square_panels)) 10 else 12
    } else {
      as.numeric(save_waves_width[[1]])
    }
    h_waves <- if (is.null(save_waves_height)) {
      if (isTRUE(square_panels)) w_waves else 8
    } else {
      as.numeric(save_waves_height[[1]])
    }
    w_heat <- if (is.null(save_heatmap_width)) 13 else as.numeric(save_heatmap_width[[1]])
    h_heat <- if (is.null(save_heatmap_height)) 6.5 else as.numeric(save_heatmap_height[[1]])
    .hc_ggsave_pdf_png(
      filename = base::file.path(out_dir, base::paste0(file_prefix, "_Waves.pdf")),
      plot = p_waves,
      width = w_waves,
      height = h_waves,
      units = "in",
      device = grDevices::cairo_pdf
    )
    .hc_ggsave_pdf_png(
      filename = base::file.path(out_dir, base::paste0(file_prefix, "_Heatmap.pdf")),
      plot = p_heat,
      width = w_heat,
      height = h_heat,
      units = "in",
      device = grDevices::cairo_pdf
    )
  }

  list(
    module_cluster_waves = p_waves,
    module_cluster_heatmap = p_heat
  )
}

#' Plot longitudinal endotype results
#'
#' Creates a PCA donor plot and a module-wave plot (mean module value over time
#' per endotype) from `hc_longitudinal_endotype_clustering()` outputs.
#'
#' @param hc A `HCoCenaExperiment`.
#' @param slot_name Satellite slot containing endotype clustering output.
#' @param save_pdf Logical. If `TRUE`, save PDFs into the configured output
#'   folder.
#' @param file_prefix Prefix for generated PDF names.
#' @param pca_width,pca_height PDF size for PCA plot.
#' @param waves_width,waves_height PDF size for trajectory plot.
#'
#' @return A named list with `pca` and `waves` ggplot objects.
#' @examples
#' hc <- hc_example_data("clustered")
#' hc <- hc_longitudinal_module_means(hc, donor_col = "donor", time_col = "timepoint")
#' hc <- hc_longitudinal_endotype_clustering(hc, k = 2)
#' p <- hc_plot_longitudinal_endotypes(hc, save_pdf = FALSE)
#' @export
hc_plot_longitudinal_endotypes <- function(hc,
                                           slot_name = "longitudinal_endotypes",
                                           save_pdf = TRUE,
                                           file_prefix = "Longitudinal_Endotypes",
                                           pca_width = 7,
                                           pca_height = 5,
                                           waves_width = 12,
                                           waves_height = 8) {
  if (!inherits(hc, "HCoCenaExperiment")) {
    stop("`hc` must be a `HCoCenaExperiment`.")
  }

  sat <- as.list(hc@satellite)
  obj <- sat[[slot_name]]
  if (is.null(obj) || is.null(obj$pca) || is.null(obj$trajectory_mean)) {
    stop("Slot `", slot_name, "` does not contain endotype plotting data.")
  }

  pca_df <- base::as.data.frame(obj$pca, stringsAsFactors = FALSE)
  pca_df$endotype <- base::factor(pca_df$endotype, levels = base::sort(base::unique(pca_df$endotype)))
  endo_lv <- base::levels(pca_df$endotype)
  endo_pal <- .hc_distinct_palette(base::max(3L, base::length(endo_lv)))[base::seq_along(endo_lv)]
  endo_pal <- stats::setNames(endo_pal, endo_lv)
  pca_ell <- .hc_filter_groups_for_ellipse(
    pca_df,
    group_col = "endotype",
    min_points = 4L,
    x_col = "PC1",
    y_col = "PC2",
    min_unique = 3L
  )
  pca_groups_all <- base::sort(base::unique(base::as.character(stats::na.omit(pca_df$endotype))))
  pca_groups_ell <- if (base::nrow(pca_ell) > 0) {
    base::sort(base::unique(base::as.character(pca_ell$endotype)))
  } else {
    character(0)
  }
  show_pca_ell <- base::length(pca_groups_all) > 0 && base::identical(pca_groups_all, pca_groups_ell)
  if (!isTRUE(show_pca_ell)) {
    pca_ell <- pca_ell[0, , drop = FALSE]
  }
  pca_ctr <- .hc_cluster_centers_2d(pca_df, x_col = "PC1", y_col = "PC2", group_col = "endotype")
  pca_axis_lab <- .hc_pca_axis_labels(obj$pca_variance)

  p_pca <- ggplot2::ggplot(
    pca_df,
    ggplot2::aes(x = PC1, y = PC2)
  ) +
    ggplot2::geom_point(
      ggplot2::aes(color = endotype),
      shape = 16,
      size = 3.1,
      alpha = 0.95
    ) +
    ggplot2::scale_fill_manual(values = endo_pal, drop = FALSE) +
    ggplot2::scale_color_manual(values = endo_pal, drop = FALSE) +
    ggplot2::coord_equal() +
    .hc_theme_pub(base_size = 11.5) +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(color = "grey35", linewidth = 0.55, fill = NA)
    ) +
    ggplot2::guides(
      fill = "none",
      color = ggplot2::guide_legend(
        override.aes = list(shape = 16, size = 4, alpha = 1)
      )
    ) +
    ggplot2::labs(
      title = "Longitudinal endotypes (PCA)",
      x = pca_axis_lab$x,
      y = pca_axis_lab$y,
      color = "Endotype"
    )
  if (base::nrow(pca_ell) > 0) {
    p_pca <- p_pca +
      ggplot2::stat_ellipse(
        data = pca_ell,
        mapping = ggplot2::aes(x = PC1, y = PC2, fill = endotype, group = endotype),
        geom = "polygon",
        type = "norm",
        level = 0.72,
        alpha = 0.10,
        color = NA,
        inherit.aes = FALSE,
        show.legend = FALSE
      ) +
      ggplot2::stat_ellipse(
        data = pca_ell,
        mapping = ggplot2::aes(x = PC1, y = PC2, color = endotype, group = endotype),
        type = "norm",
        level = 0.72,
        alpha = 0.85,
        linewidth = 0.8,
        inherit.aes = FALSE,
        show.legend = FALSE
      )
  }
  if (base::nrow(pca_ctr) > 0) {
    p_pca <- p_pca +
      ggplot2::geom_label(
        data = pca_ctr,
        mapping = ggplot2::aes(x = x, y = y, label = group, fill = group),
        color = "black",
        fontface = "bold",
        size = 3.1,
        label.padding = grid::unit(0.12, "lines"),
        label.r = grid::unit(0.08, "lines"),
        label.size = 0.25,
        alpha = 0.9,
        inherit.aes = FALSE,
        show.legend = FALSE
      )
  }

  tr <- base::as.data.frame(obj$trajectory_mean, stringsAsFactors = FALSE)
  y_lab <- if (!is.null(obj$value_label) && base::nzchar(as.character(obj$value_label[[1]]))) {
    as.character(obj$value_label[[1]])
  } else {
    "Module mean"
  }
  module_levels <- base::sort(base::unique(base::as.character(tr$module)))
  module_lookup <- NULL
  if (!is.null(obj$module_lookup) && "module" %in% base::colnames(obj$module_lookup)) {
    module_lookup <- base::as.data.frame(obj$module_lookup, stringsAsFactors = FALSE)
    module_lookup <- .hc_reorder_module_lookup_natural(module_lookup)
    ml <- base::as.character(module_lookup$module)
    module_levels <- ml[ml %in% module_levels]
    if (base::length(module_levels) == 0) {
      module_levels <- .hc_natural_module_order(base::unique(base::as.character(tr$module)))
    }
  }
  if (is.null(module_lookup) || base::nrow(module_lookup) == 0) {
    module_lookup <- base::data.frame(
      module = module_levels,
      module_color = "grey60",
      stringsAsFactors = FALSE
    )
  }
  mod_col <- .hc_module_color_map(
    module_lookup = module_lookup,
    module_levels = module_levels
  )
  time_levels <- obj$time_levels
  if (!is.null(time_levels)) {
    tr$time <- base::factor(tr$time, levels = base::as.character(time_levels), ordered = TRUE)
  }
  tr$endotype <- base::factor(tr$endotype, levels = base::sort(base::unique(tr$endotype)))
  tr$module <- base::factor(base::as.character(tr$module), levels = module_levels)

  p_waves <- ggplot2::ggplot(
    tr,
    ggplot2::aes(
      x = time,
      y = value,
      color = endotype,
      group = endotype
    )
  ) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_point(size = 2) +
    .hc_module_facet_wrap(
      module_levels = module_levels,
      mod_col = mod_col,
      ncol = 4,
      scales = "free_y",
      facets = ~module
    ) +
    .hc_theme_pub(base_size = 11) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    ) +
    ggplot2::labs(
      title = "Module trajectories per endotype",
      x = "Timepoint",
      y = y_lab,
      color = "Endotype"
    )
  p_waves <- .hc_colorize_module_strips(
    p = p_waves,
    module_levels = module_levels,
    mod_col = mod_col
  )

  if (isTRUE(save_pdf)) {
    out_dir <- .hc_resolve_output_dir(hc)

    .hc_ggsave_pdf_png(
      filename = base::file.path(out_dir, base::paste0(file_prefix, "_PCA.pdf")),
      plot = p_pca,
      width = pca_width,
      height = pca_height,
      units = "in",
      device = grDevices::cairo_pdf
    )
    .hc_ggsave_pdf_png(
      filename = base::file.path(out_dir, base::paste0(file_prefix, "_ModuleWaves.pdf")),
      plot = p_waves,
      width = waves_width,
      height = waves_height,
      units = "in",
      device = grDevices::cairo_pdf
    )
  }

  list(pca = p_pca, waves = p_waves)
}

#' Compute module-level CAP matrix from repeated trajectory clustering
#'
#' CAP (cluster assignment proportions) quantifies for each donor the
#' proportion of repeated clustering runs assigned to each module-specific
#' cluster. This is the generalized, data-agnostic analogue of the CAP concept
#' used in the legacy longitudinal scripts.
#'
#' @param hc A `HCoCenaExperiment`.
#' @param slot_name Satellite slot containing outputs from
#'   `hc_longitudinal_endotype_clustering()`.
#' @param output_slot Satellite slot to write results to. Defaults to
#'   `slot_name` (in-place update).
#' @param runs Number of repeated clustering runs per module (for `kmeans`).
#' @param method One of `"kmeans"` or `"hclust"` for repeated module
#'   clustering.
#' @param nstart `kmeans` random starts per run.
#' @param seed Random seed.
#' @param scale_features Logical. If `TRUE`, z-score module trajectories before
#'   repeated clustering.
#' @param na_impute One of `"median"` or `"zero"` for NA imputation in module
#'   trajectory matrices.
#' @param jitter_sd Optional Gaussian jitter SD added to features per run
#'   (`kmeans` only) to probe local assignment stability.
#'
#' @return Updated `HCoCenaExperiment` with CAP outputs in
#'   `hc@satellite[[output_slot]]`.
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
#' hc <- hc_longitudinal_module_cap(hc, runs = 1)
#' @export
hc_longitudinal_module_cap <- function(hc,
                                       slot_name = "longitudinal_endotypes",
                                       output_slot = slot_name,
                                       runs = 50,
                                       method = c("kmeans", "hclust"),
                                       nstart = 30,
                                       seed = 42,
                                       scale_features = TRUE,
                                       na_impute = c("median", "zero"),
                                       jitter_sd = 0) {
  if (!inherits(hc, "HCoCenaExperiment")) {
    stop("`hc` must be a `HCoCenaExperiment`.")
  }
  method <- base::match.arg(method)
  na_impute <- base::match.arg(na_impute)

  runs <- as.integer(runs)
  if (!base::is.finite(runs) || runs < 1) {
    runs <- 1L
  }
  if (identical(method, "hclust")) {
    runs <- 1L
  }

  sat <- as.list(hc@satellite)
  obj <- sat[[slot_name]]
  req <- c("module_cluster_best_k", "donor_time_module", "module_lookup", "time_levels")
  if (is.null(obj) || !base::all(req %in% base::names(obj))) {
    stop("Slot `", slot_name, "` does not contain required longitudinal clustering outputs.")
  }

  donor_time_module <- base::as.data.frame(obj$donor_time_module, stringsAsFactors = FALSE)
  if (base::nrow(donor_time_module) == 0) {
    stop("No `donor_time_module` data found in slot `", slot_name, "`.")
  }
  donors <- if (!is.null(obj$module_cluster_matrix)) {
    base::rownames(base::as.matrix(obj$module_cluster_matrix))
  } else {
    base::sort(base::unique(base::as.character(donor_time_module$donor)))
  }
  donors <- base::as.character(donors)
  time_levels <- base::as.character(obj$time_levels)
  if (base::length(time_levels) == 0) {
    time_levels <- base::sort(base::unique(base::as.character(donor_time_module$time)))
  }

  module_lookup <- base::as.data.frame(obj$module_lookup, stringsAsFactors = FALSE)
  module_ids <- base::as.character(module_lookup$module)
  module_best_k <- base::as.data.frame(obj$module_cluster_best_k, stringsAsFactors = FALSE)
  if (base::nrow(module_best_k) == 0) {
    stop("No `module_cluster_best_k` values found in slot `", slot_name, "`.")
  }

  traj_impute_mode <- "none"
  if (!is.null(obj$source_slot) && obj$source_slot %in% base::names(sat)) {
    src <- sat[[obj$source_slot]]
    if (!is.null(src$impute_missing) && src$impute_missing %in% c("none", "linear", "locf")) {
      traj_impute_mode <- src$impute_missing
    }
  }

  cap_parts <- list()
  cap_assign_rows <- list()

  for (i in base::seq_along(module_ids)) {
    m <- module_ids[[i]]
    kk <- module_best_k$best_k[base::match(m, module_best_k$module)]
    if (base::length(kk) == 0 || base::is.na(kk) || !base::is.finite(kk)) {
      kk <- 1L
    }
    kk <- base::as.integer(kk[[1]])
    if (kk < 1) {
      kk <- 1L
    }

    mm <- .hc_module_time_matrix(
      donor_time_module = donor_time_module,
      donors = donors,
      module_id = m,
      time_levels = time_levels,
      traj_impute_mode = traj_impute_mode,
      na_impute = na_impute
    )
    x <- mm$use
    if (isTRUE(scale_features)) {
      x <- base::as.matrix(base::scale(x))
      x[!base::is.finite(x)] <- 0
    }

    counts <- base::matrix(
      0,
      nrow = base::length(donors),
      ncol = kk,
      dimnames = list(donors, base::seq_len(kk))
    )

    if (kk == 1L || base::nrow(x) <= 1) {
      counts[, 1] <- runs
    } else {
      ref_labels <- NULL
      if (!is.null(obj$module_cluster_assignments)) {
        ma <- obj$module_cluster_assignments
        ma <- ma[ma$module == m, , drop = FALSE]
        if (base::nrow(ma) > 0) {
          ref_labels <- ma$cluster[base::match(donors, ma$donor)]
        }
      }
      if (is.null(ref_labels) || base::any(base::is.na(ref_labels))) {
        ref_labels <- .hc_cluster_labels(
          x = x,
          k = kk,
          method = method,
          seed = seed + i,
          nstart = nstart
        )
      }
      ref_labels <- base::as.integer(ref_labels)
      ref_centroids <- .hc_cluster_centroids(x = x, labels = ref_labels, k = kk)

      for (rr in base::seq_len(runs)) {
        x_run <- x
        if (identical(method, "kmeans") && base::is.finite(jitter_sd) && jitter_sd > 0) {
          x_run <- x_run + stats::rnorm(n = base::length(x_run), mean = 0, sd = jitter_sd)
        }
        lbl <- .hc_cluster_labels(
          x = x_run,
          k = kk,
          method = method,
          seed = seed + rr + i * 1000,
          nstart = nstart
        )
        lbl <- .hc_align_labels_to_reference(
          labels = lbl,
          x = x_run,
          ref_centroids = ref_centroids
        )
        for (cl in base::seq_len(kk)) {
          counts[lbl == cl, cl] <- counts[lbl == cl, cl] + 1
        }
      }
    }

    prop <- counts / runs
    keys <- base::paste0(m, "__", base::seq_len(kk))
    base::colnames(prop) <- keys
    cap_parts[[m]] <- prop

    winner <- base::max.col(prop, ties.method = "first")
    cap_assign_rows[[base::length(cap_assign_rows) + 1]] <- base::data.frame(
      donor = donors,
      module = m,
      cluster = base::as.integer(winner),
      module_cluster_key = base::paste0(m, "__", winner),
      max_prop = prop[base::cbind(base::seq_len(base::nrow(prop)), winner)],
      stringsAsFactors = FALSE
    )
  }

  non_empty <- cap_parts[base::vapply(cap_parts, function(x) !base::is.null(x), logical(1))]
  if (base::length(non_empty) == 0) {
    stop("No CAP components were generated.")
  }
  cap_matrix <- base::do.call(base::cbind, non_empty)
  cap_matrix <- base::as.matrix(cap_matrix)
  cap_matrix <- cap_matrix[donors, , drop = FALSE]

  cap_long <- base::as.data.frame(base::as.table(cap_matrix), stringsAsFactors = FALSE)
  base::colnames(cap_long) <- c("donor", "module_cluster_key", "proportion")
  cap_long$module <- sub("__.*$", "", base::as.character(cap_long$module_cluster_key))
  cap_long$cluster <- .hc_as_integer_safely(sub("^.*__", "", base::as.character(cap_long$module_cluster_key)))
  cap_long$proportion <- as.numeric(cap_long$proportion)

  cap_assign <- if (base::length(cap_assign_rows) > 0) {
    base::do.call(base::rbind, cap_assign_rows)
  } else {
    base::data.frame(
      donor = character(0),
      module = character(0),
      cluster = integer(0),
      module_cluster_key = character(0),
      max_prop = numeric(0),
      stringsAsFactors = FALSE
    )
  }

  obj$cap_matrix <- cap_matrix
  obj$cap_long <- cap_long
  obj$cap_consensus_assignments <- cap_assign
  obj$cap_method <- method
  obj$cap_runs <- as.integer(runs)
  obj$cap_nstart <- as.integer(nstart)
  obj$cap_seed <- as.integer(seed)
  obj$cap_scale_features <- isTRUE(scale_features)
  obj$cap_na_impute <- na_impute
  obj$cap_jitter_sd <- as.numeric(jitter_sd)

  sat[[output_slot]] <- obj
  hc@satellite <- S4Vectors::SimpleList(sat)
  methods::validObject(hc)
  hc
}

#' Meta-cluster donors based on longitudinal CAP features
#'
#' Clusters donors on CAP features (or alternative longitudinal feature
#' matrices), providing CH-based model selection and PCA/UMAP embeddings.
#'
#' @param hc A `HCoCenaExperiment`.
#' @param slot_name Satellite slot containing longitudinal outputs.
#' @param output_slot Satellite slot to write results to (default in-place).
#' @param feature_source One of `"cap_matrix"`, `"module_cluster_matrix"`, or
#'   `"feature_matrix_used"`.
#' @param k Candidate meta-cluster numbers.
#' @param method One of `"kmeans"` or `"hclust"`.
#' @param scale_features Logical. If `TRUE`, z-score features before clustering.
#' @param nstart Number of random starts for `kmeans`.
#' @param seed Random seed.
#' @param na_impute One of `"median"` or `"zero"`.
#' @param cluster_prefix Prefix for meta-cluster labels.
#' @param consensus Logical. If `TRUE`, run consensus clustering with repeated
#'   subsampling for each candidate `k` and select the best `k` via consensus
#'   separation.
#' @param consensus_runs Number of consensus resampling runs per candidate `k`.
#' @param consensus_sample_fraction Fraction of donors sampled per consensus run.
#' @param consensus_feature_fraction Fraction of features sampled per consensus
#'   run.
#' @param consensus_linkage Linkage method for final hierarchical clustering on
#'   the consensus distance (`1 - consensus`).
#' @param compute_umap Logical. If `TRUE`, compute UMAP coordinates when
#'   `uwot` or `umap` is available.
#' @param umap_neighbors UMAP neighborhood size.
#' @param umap_min_dist UMAP minimum distance.
#'
#' @return Updated `HCoCenaExperiment`.
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
#'   seed = 1
#' )
#' @export
hc_longitudinal_meta_clustering <- function(hc,
                                            slot_name = "longitudinal_endotypes",
                                            output_slot = slot_name,
                                            feature_source = c("cap_matrix", "module_cluster_matrix", "feature_matrix_used"),
                                            k = 2:8,
                                            method = c("kmeans", "hclust"),
                                            scale_features = TRUE,
                                            nstart = 50,
                                            seed = 42,
                                            na_impute = c("median", "zero"),
                                            cluster_prefix = "MC",
                                            consensus = TRUE,
                                            consensus_runs = 200,
                                            consensus_sample_fraction = 0.8,
                                            consensus_feature_fraction = 0.8,
                                            consensus_linkage = "average",
                                            compute_umap = TRUE,
                                            umap_neighbors = 15,
                                            umap_min_dist = 0.3) {
  if (!inherits(hc, "HCoCenaExperiment")) {
    stop("`hc` must be a `HCoCenaExperiment`.")
  }
  feature_source <- base::match.arg(feature_source)
  method <- base::match.arg(method)
  na_impute <- base::match.arg(na_impute)

  sat <- as.list(hc@satellite)
  obj <- sat[[slot_name]]
  if (is.null(obj)) {
    stop("Slot `", slot_name, "` not found.")
  }

  if (!feature_source %in% base::names(obj)) {
    stop("Feature source `", feature_source, "` not found in slot `", slot_name, "`.")
  }

  x <- base::as.matrix(obj[[feature_source]])
  if (base::nrow(x) < 3) {
    stop("Need at least 3 donors for meta-clustering.")
  }

  if (identical(na_impute, "median")) {
    for (j in base::seq_len(base::ncol(x))) {
      v <- x[, j]
      if (base::any(base::is.na(v))) {
        med <- stats::median(v, na.rm = TRUE)
        if (!base::is.finite(med)) {
          med <- 0
        }
        v[base::is.na(v)] <- med
        x[, j] <- v
      }
    }
  } else {
    x[base::is.na(x)] <- 0
  }

  keep_var <- base::apply(x, 2, stats::sd, na.rm = TRUE) > 0
  x <- x[, keep_var, drop = FALSE]
  if (base::ncol(x) == 0) {
    stop("No variable features available for meta-clustering.")
  }

  x_use <- base::as.matrix(x)
  if (isTRUE(scale_features)) {
    x_use <- base::as.matrix(base::scale(x_use))
    x_use[!base::is.finite(x_use)] <- 0
  }

  donor_ids <- base::rownames(x_use)
  fit <- NULL
  cl <- NULL
  score_table <- NULL
  best_k <- NULL
  consensus_matrix <- NULL
  donor_stability <- NULL
  method_comparison <- NULL
  method_donor_assign <- list()

  if (!isTRUE(consensus)) {
    fit <- .hc_cluster_with_scoring(
      x = x_use,
      k = k,
      method = method,
      nstart = nstart,
      seed = seed
    )
    score_table <- fit$score_table
    best_k <- fit$best_k
    cl <- base::as.integer(fit$cluster)
    donor_stability <- base::rep(NA_real_, base::length(cl))
  } else {
    k_vec <- base::sort(base::unique(base::as.integer(k)))
    k_vec <- k_vec[base::is.finite(k_vec) & k_vec >= 2 & k_vec < base::nrow(x_use)]
    if (base::length(k_vec) == 0) {
      stop("No valid `k` left for consensus meta-clustering. Need values in [2, n_donors-1].")
    }

    cons_results <- list()
    score_rows <- list()
    for (kk in k_vec) {
      cr <- .hc_consensus_one_k(
        x = x_use,
        k = kk,
        method = method,
        runs = consensus_runs,
        sample_fraction = consensus_sample_fraction,
        feature_fraction = consensus_feature_fraction,
        nstart = nstart,
        seed = seed + kk * 1000L,
        linkage = consensus_linkage
      )
      cons_results[[as.character(kk)]] <- cr
      score_rows[[base::length(score_rows) + 1]] <- base::data.frame(
        k = kk,
        ch_index = cr$delta,
        tot_withinss = NA_real_,
        consensus_within = cr$within,
        consensus_between = cr$between,
        consensus_delta = cr$delta,
        pac = cr$pac,
        observed_pair_fraction = cr$observed_pair_fraction,
        stringsAsFactors = FALSE
      )
    }
    score_table <- base::do.call(base::rbind, score_rows)
    rank_score <- score_table$consensus_delta
    rank_score[!base::is.finite(rank_score)] <- -Inf
    rank_tie <- score_table$pac
    rank_tie[!base::is.finite(rank_tie)] <- Inf
    ord <- base::order(-rank_score, rank_tie, score_table$k)
    score_table <- score_table[ord, , drop = FALSE]
    best_k <- as.integer(score_table$k[[1]])

    best_res <- cons_results[[as.character(best_k)]]
    cl <- base::as.integer(best_res$labels)
    consensus_matrix <- best_res$consensus_matrix
    donor_stability <- best_res$stability
  }

  method_rows <- list()
  for (mm in c("kmeans", "hclust")) {
    cmp <- tryCatch(
      .hc_cluster_with_scoring(
        x = x_use,
        k = k,
        method = mm,
        nstart = nstart,
        seed = seed
      ),
      error = function(e) e
    )
    if (inherits(cmp, "error")) {
      method_rows[[base::length(method_rows) + 1]] <- base::data.frame(
        method = mm,
        best_k_ch = NA_integer_,
        best_ch_index = NA_real_,
        n_clusters_at_best_k = NA_integer_,
        status = base::conditionMessage(cmp),
        stringsAsFactors = FALSE
      )
    } else {
      st <- base::as.data.frame(cmp$score_table, stringsAsFactors = FALSE)
      idx <- base::match(as.integer(cmp$best_k), st$k)
      best_ch <- if (!base::is.na(idx)) as.numeric(st$ch_index[[idx]]) else NA_real_
      cmp_best_k <- as.integer(cmp$best_k)
      cmp_cluster <- base::as.integer(cmp$cluster)
      method_rows[[base::length(method_rows) + 1]] <- base::data.frame(
        method = mm,
        best_k_ch = cmp_best_k,
        best_ch_index = best_ch,
        n_clusters_at_best_k = base::length(base::unique(cmp_cluster)),
        status = "ok",
        stringsAsFactors = FALSE
      )
      method_donor_assign[[mm]] <- base::data.frame(
        donor = donor_ids,
        method = mm,
        best_k_ch = cmp_best_k,
        cluster_id = cmp_cluster,
        cluster_label = base::paste0(cluster_prefix, cmp_cluster),
        stringsAsFactors = FALSE
      )
    }
  }
  method_comparison <- base::do.call(base::rbind, method_rows)
  method_comparison$selected_method <- method_comparison$method == method
  method_comparison$selected_k_pipeline <- ifelse(
    method_comparison$selected_method,
    as.integer(best_k),
    NA_integer_
  )
  n_meta_pipeline <- base::length(base::unique(cl))
  method_comparison$n_meta_clusters_pipeline <- ifelse(
    method_comparison$selected_method,
    n_meta_pipeline,
    NA_integer_
  )
  method_comparison$selection_mode <- if (isTRUE(consensus)) "consensus" else "ch_index"

  donor_clusterings <- base::data.frame(
    donor = donor_ids,
    selected_method = method,
    selected_k = as.integer(best_k),
    selected_cluster_id = as.integer(cl),
    selected_cluster = base::paste0(cluster_prefix, cl),
    stringsAsFactors = FALSE
  )
  if (!base::is.null(donor_stability)) {
    donor_clusterings$selected_stability <- as.numeric(donor_stability)
  }
  if (base::length(method_donor_assign) > 0) {
    for (mm in base::names(method_donor_assign)) {
      tmp <- method_donor_assign[[mm]]
      if (base::is.null(tmp) || base::nrow(tmp) == 0) {
        next
      }
      k_mm <- base::unique(base::as.integer(tmp$best_k_ch))
      k_mm <- k_mm[base::is.finite(k_mm)]
      k_lbl <- if (base::length(k_mm) > 0) k_mm[[1]] else NA_integer_
      col_id <- base::paste0("cluster_id_", mm, "_k", k_lbl)
      col_lb <- base::paste0("cluster_", mm, "_k", k_lbl)
      donor_clusterings[[col_id]] <- tmp$cluster_id[base::match(donor_clusterings$donor, tmp$donor)]
      donor_clusterings[[col_lb]] <- tmp$cluster_label[base::match(donor_clusterings$donor, tmp$donor)]
    }
  }

  meta_cluster <- base::data.frame(
    donor = donor_ids,
    cluster = cl,
    meta_cluster = base::paste0(cluster_prefix, cl),
    stability = donor_stability,
    stringsAsFactors = FALSE
  )

  pca_fit <- stats::prcomp(x_use, center = TRUE, scale. = FALSE)
  meta_pca_var <- 100 * (pca_fit$sdev^2 / sum(pca_fit$sdev^2))
  meta_pca <- base::data.frame(
    donor = donor_ids,
    PC1 = pca_fit$x[, 1],
    PC2 = if (base::ncol(pca_fit$x) >= 2) pca_fit$x[, 2] else base::rep(0, base::nrow(pca_fit$x)),
    meta_cluster = meta_cluster$meta_cluster,
    stringsAsFactors = FALSE
  )

  meta_umap <- NULL
  if (isTRUE(compute_umap)) {
    if (requireNamespace("uwot", quietly = TRUE)) {
      um <- .hc_with_seed(seed, uwot::umap(
        x_use,
        n_neighbors = as.integer(umap_neighbors),
        min_dist = as.numeric(umap_min_dist),
        metric = "euclidean",
        verbose = FALSE
      ))
      meta_umap <- base::data.frame(
        donor = donor_ids,
        UMAP1 = um[, 1],
        UMAP2 = um[, 2],
        meta_cluster = meta_cluster$meta_cluster,
        stringsAsFactors = FALSE
      )
    } else if (requireNamespace("umap", quietly = TRUE)) {
      um <- .hc_with_seed(seed, umap::umap(x_use))
      meta_umap <- base::data.frame(
        donor = donor_ids,
        UMAP1 = um$layout[, 1],
        UMAP2 = um$layout[, 2],
        meta_cluster = meta_cluster$meta_cluster,
        stringsAsFactors = FALSE
      )
    } else {
      warning("UMAP not computed: neither `uwot` nor `umap` is installed.")
    }
  }

  obj$meta_feature_source <- feature_source
  obj$meta_feature_matrix_used <- x_use
  obj$meta_method <- method
  obj$meta_candidate_k <- k
  obj$meta_nstart <- as.integer(nstart)
  obj$meta_seed <- as.integer(seed)
  obj$meta_score_table <- score_table
  obj$meta_best_k <- best_k
  obj$meta_cluster <- meta_cluster
  obj$meta_donor_clusterings <- donor_clusterings
  obj$meta_pca <- meta_pca
  obj$meta_pca_variance <- meta_pca_var
  obj$meta_umap <- meta_umap
  obj$meta_method_comparison <- method_comparison
  obj$meta_consensus <- isTRUE(consensus)
  obj$meta_consensus_matrix <- consensus_matrix
  obj$meta_consensus_runs <- as.integer(consensus_runs)
  obj$meta_consensus_sample_fraction <- as.numeric(consensus_sample_fraction)
  obj$meta_consensus_feature_fraction <- as.numeric(consensus_feature_fraction)
  obj$meta_consensus_linkage <- as.character(consensus_linkage)
  obj$meta_selection_tables <- .hc_meta_report_tables(obj, detail = "full")

  if (!is.null(obj$donor_cluster)) {
    obj$donor_cluster$meta_cluster <- meta_cluster$meta_cluster[
      base::match(obj$donor_cluster$donor, meta_cluster$donor)
    ]
  }

  sat[[output_slot]] <- obj
  hc@satellite <- S4Vectors::SimpleList(sat)
  methods::validObject(hc)
  hc
}

#' Plot k-criterion diagnostics for longitudinal clustering
#'
#' @param hc A `HCoCenaExperiment`.
#' @param slot_name Satellite slot with longitudinal clustering outputs.
#' @param save_pdf Logical. If `TRUE`, save PDFs.
#' @param file_prefix Prefix for output files.
#' @param facet_ncol Number of facet columns for module-wise diagnostics.
#' @param global_scope One of `"auto"`, `"meta"`, `"endotype"`.
#'   Controls which global k-selection table is shown.
#'
#' @return List with ggplots `global_k` and `module_k`.
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
#' p <- hc_plot_longitudinal_k_criterion(hc, save_pdf = FALSE)
#' @export
hc_plot_longitudinal_k_criterion <- function(hc,
                                             slot_name = "longitudinal_endotypes",
                                             save_pdf = TRUE,
                                             file_prefix = "Longitudinal_KCriterion",
                                             facet_ncol = 4,
                                             global_scope = c("auto", "meta", "endotype")) {
  if (!inherits(hc, "HCoCenaExperiment")) {
    stop("`hc` must be a `HCoCenaExperiment`.")
  }
  global_scope <- base::match.arg(global_scope)
  sat <- as.list(hc@satellite)
  obj <- sat[[slot_name]]
  if (is.null(obj) || is.null(obj$score_table) || is.null(obj$module_cluster_score_table)) {
    stop("Slot `", slot_name, "` does not contain k-criterion tables.")
  }

  use_meta <- FALSE
  if (identical(global_scope, "meta")) {
    use_meta <- TRUE
  } else if (identical(global_scope, "auto")) {
    use_meta <- !is.null(obj$meta_score_table) && !is.null(obj$meta_best_k)
  }

  if (use_meta) {
    if (is.null(obj$meta_score_table) || is.null(obj$meta_best_k)) {
      stop("Meta score table missing. Run `hc_longitudinal_meta_clustering()` first or use `global_scope = \"endotype\"`.")
    }
    gl <- base::as.data.frame(obj$meta_score_table, stringsAsFactors = FALSE)
    best_k_global <- as.integer(obj$meta_best_k)
  } else {
    gl <- base::as.data.frame(obj$score_table, stringsAsFactors = FALSE)
    best_k_global <- as.integer(obj$best_k)
  }
  gl$selected <- gl$k == best_k_global

  is_consensus_global <- use_meta && isTRUE(obj$meta_consensus) && "consensus_delta" %in% base::colnames(gl)
  y_lab_global <- if (is_consensus_global) {
    "Consensus separation (within - between)"
  } else {
    "Calinski-Harabasz"
  }
  title_global <- if (is_consensus_global) {
    "Global meta-cluster k selection (consensus)"
  } else {
    "Global meta-cluster k selection (CH index)"
  }

  p_global <- ggplot2::ggplot(gl, ggplot2::aes(x = k, y = ch_index)) +
    ggplot2::geom_line(color = "#2F4F4F", linewidth = 0.8) +
    ggplot2::geom_point(ggplot2::aes(size = selected, color = selected), show.legend = FALSE) +
    ggplot2::geom_vline(xintercept = as.numeric(best_k_global), linetype = "dashed", color = "#B22222", linewidth = 0.7) +
    ggplot2::scale_color_manual(values = c("FALSE" = "#2F4F4F", "TRUE" = "#B22222")) +
    ggplot2::scale_size_manual(values = c("FALSE" = 2.3, "TRUE" = 3.3)) +
    .hc_theme_pub(base_size = 11) +
    ggplot2::labs(
      title = title_global,
      x = "k",
      y = y_lab_global
    )

  ms <- base::as.data.frame(obj$module_cluster_score_table, stringsAsFactors = FALSE)
  mb <- base::as.data.frame(obj$module_cluster_best_k, stringsAsFactors = FALSE)
  lookup <- base::as.data.frame(obj$module_lookup, stringsAsFactors = FALSE)
  module_levels <- base::as.character(lookup$module)
  mod_col <- .hc_module_color_map(
    module_lookup = lookup,
    module_levels = module_levels
  )
  ms$module <- base::factor(base::as.character(ms$module), levels = module_levels)
  ms$best_k <- mb$best_k[base::match(base::as.character(ms$module), mb$module)]
  ms$selected <- ms$k == ms$best_k

  p_module <- ggplot2::ggplot(ms, ggplot2::aes(x = k, y = ch_index, group = module, color = module)) +
    ggplot2::geom_line(alpha = 0.4, linewidth = 0.5, show.legend = FALSE) +
    ggplot2::geom_point(ggplot2::aes(size = selected), show.legend = FALSE) +
    .hc_module_facet_wrap(
      module_levels = module_levels,
      mod_col = mod_col,
      ncol = facet_ncol,
      scales = "free_y",
      facets = ~module
    ) +
    ggplot2::scale_color_manual(values = mod_col, drop = FALSE) +
    ggplot2::scale_size_manual(values = c("FALSE" = 1.5, "TRUE" = 2.6)) +
    .hc_theme_pub(base_size = 10) +
    ggplot2::labs(
      title = "Module-wise trajectory k selection (CH index)",
      x = "k",
      y = "Calinski-Harabasz"
    )
  p_module <- .hc_colorize_module_strips(
    p = p_module,
    module_levels = module_levels,
    mod_col = mod_col
  )

  if (isTRUE(save_pdf)) {
    out_dir <- .hc_resolve_output_dir(hc)
    .hc_ggsave_pdf_png(
      filename = base::file.path(out_dir, base::paste0(file_prefix, "_Global.pdf")),
      plot = p_global,
      width = 7,
      height = 5,
      units = "in",
      device = grDevices::cairo_pdf
    )
    .hc_ggsave_pdf_png(
      filename = base::file.path(out_dir, base::paste0(file_prefix, "_Modules.pdf")),
      plot = p_module,
      width = 12,
      height = 8,
      units = "in",
      device = grDevices::cairo_pdf
    )
  }

  list(global_k = p_global, module_k = p_module)
}

#' Plot CAP heatmap (donor x module-cluster proportions)
#'
#' @param hc A `HCoCenaExperiment`.
#' @param slot_name Satellite slot with CAP outputs.
#' @param save_pdf Logical. If `TRUE`, save PDF.
#' @param file_prefix Prefix for output filename.
#' @param show_values Logical. If `TRUE`, show proportions as text.
#' @param min_label_prop Minimum proportion threshold for text labels.
#' @param donor_order Optional donor order for x-axis. If `NULL`, donors are
#'   ordered by hierarchical clustering on the CAP matrix.
#'
#' @return List with `cap_heatmap` ggplot.
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
#' p <- hc_plot_longitudinal_cap(hc, save_pdf = FALSE)
#' @export
hc_plot_longitudinal_cap <- function(hc,
                                     slot_name = "longitudinal_endotypes",
                                     save_pdf = TRUE,
                                     file_prefix = "Longitudinal_CAP",
                                     show_values = FALSE,
                                     min_label_prop = 0.5,
                                     donor_order = NULL) {
  if (!inherits(hc, "HCoCenaExperiment")) {
    stop("`hc` must be a `HCoCenaExperiment`.")
  }
  sat <- as.list(hc@satellite)
  obj <- sat[[slot_name]]
  if (is.null(obj) || is.null(obj$cap_matrix) || is.null(obj$module_lookup)) {
    stop("Slot `", slot_name, "` does not contain CAP outputs. Run `hc_longitudinal_module_cap()` first.")
  }

  cap_mat <- base::as.matrix(obj$cap_matrix)
  if (base::nrow(cap_mat) == 0 || base::ncol(cap_mat) == 0) {
    stop("CAP matrix is empty.")
  }

  donor_order_auto <- base::rownames(cap_mat)
  if (base::nrow(cap_mat) > 1) {
    donor_order_auto <- base::rownames(cap_mat)[stats::hclust(stats::dist(cap_mat), method = "ward.D2")$order]
  }
  if (base::is.null(donor_order)) {
    donor_order <- donor_order_auto
  } else {
    donor_order <- base::as.character(donor_order)
    donor_order <- donor_order[!base::is.na(donor_order) & base::nzchar(donor_order)]
    donor_order <- donor_order[donor_order %in% base::rownames(cap_mat)]
    donor_order <- base::unique(donor_order)
    if (base::length(donor_order) == 0) {
      donor_order <- donor_order_auto
    } else {
      donor_order <- base::c(donor_order, base::setdiff(base::rownames(cap_mat), donor_order))
    }
  }

  module_lookup <- base::as.data.frame(obj$module_lookup, stringsAsFactors = FALSE)
  module_lookup <- .hc_reorder_module_lookup_natural(module_lookup)
  module_lookup$module_color <- ifelse(
    base::is.na(module_lookup$module_color) | !base::nzchar(module_lookup$module_color),
    "grey60",
    module_lookup$module_color
  )
  module_levels <- base::as.character(module_lookup$module)
  mod_col <- stats::setNames(base::as.character(module_lookup$module_color), module_levels)

  key_df <- base::data.frame(
    module_cluster_key = base::colnames(cap_mat),
    stringsAsFactors = FALSE
  )
  key_df$module <- sub("__.*$", "", key_df$module_cluster_key)
  key_df$cluster <- .hc_as_integer_safely(sub("^.*__", "", key_df$module_cluster_key))
  key_df$cluster[base::is.na(key_df$cluster)] <- 1L
  key_df <- key_df[base::order(base::match(key_df$module, module_levels), key_df$cluster), , drop = FALSE]
  key_order <- key_df$module_cluster_key

  hm_df <- base::as.data.frame(base::as.table(cap_mat), stringsAsFactors = FALSE)
  base::colnames(hm_df) <- c("donor", "module_cluster_key", "proportion")
  hm_df$donor <- base::factor(base::as.character(hm_df$donor), levels = donor_order)
  hm_df$module_cluster_key <- base::factor(base::as.character(hm_df$module_cluster_key), levels = key_order)
  hm_df$module <- sub("__.*$", "", base::as.character(hm_df$module_cluster_key))
  hm_df$base_color <- mod_col[base::match(hm_df$module, base::names(mod_col))]
  hm_df$base_color[base::is.na(hm_df$base_color)] <- "grey60"
  hm_df$proportion <- as.numeric(hm_df$proportion)
  hm_df$fill_color <- base::mapply(
    FUN = .hc_mix_with_white,
    base_color = hm_df$base_color,
    p = hm_df$proportion,
    USE.NAMES = FALSE
  )
  hm_rgb <- grDevices::col2rgb(hm_df$fill_color) / 255
  hm_lum <- 0.2126 * hm_rgb[1, ] + 0.7152 * hm_rgb[2, ] + 0.0722 * hm_rgb[3, ]
  hm_df$text_color <- ifelse(hm_lum < 0.55, "white", "black")

  p_cap <- ggplot2::ggplot(
    hm_df,
    ggplot2::aes(x = donor, y = module_cluster_key, fill = fill_color)
  ) +
    ggplot2::geom_tile(color = "grey95", linewidth = 0.2) +
    ggplot2::scale_fill_identity() +
    .hc_theme_pub(base_size = 10) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),
      axis.text.y = ggplot2::element_text(size = 7),
      axis.title = ggplot2::element_text(face = "bold")
    ) +
    ggplot2::labs(
      title = "Cluster assignment proportions (CAP) per donor and module",
      x = "Donor",
      y = "Module cluster"
    )

  if (isTRUE(show_values)) {
    p_cap <- p_cap +
      ggplot2::geom_text(
        data = hm_df[hm_df$proportion >= min_label_prop, , drop = FALSE],
        ggplot2::aes(label = base::sprintf("%.2f", proportion), color = text_color),
        size = 2.1,
        show.legend = FALSE
      ) +
      ggplot2::scale_color_identity()
  }

  if (isTRUE(save_pdf)) {
    out_dir <- .hc_resolve_output_dir(hc)
    .hc_ggsave_pdf_png(
      filename = base::file.path(out_dir, base::paste0(file_prefix, "_Heatmap.pdf")),
      plot = p_cap,
      width = 13,
      height = 7,
      units = "in",
      device = grDevices::cairo_pdf
    )
  }

  list(cap_heatmap = p_cap)
}

#' Plot meta-clustering embeddings and endotype cross-tab
#'
#' @param hc A `HCoCenaExperiment`.
#' @param slot_name Satellite slot with meta-clustering outputs.
#' @param save_pdf Logical. If `TRUE`, save plots as PDFs.
#' @param file_prefix Output filename prefix.
#' @param show_endotype_crosstab Logical. If `TRUE`, include an
#'   endotype-vs-meta-cluster overlap heatmap when endotype labels are present.
#' @param show_cluster_labels Logical. If `TRUE`, label cluster centers in PCA
#'   and UMAP.
#' @param save_tables Logical. If `TRUE`, export meta-clustering summary tables.
#' @param table_file_prefix Optional prefix for summary tables. Defaults to
#'   `file_prefix`.
#' @param table_format One of `"xlsx"` or `"csv"` for summary tables.
#' @param table_detail One of `"simple"` or `"full"`. `"simple"` exports only
#'   method-wise meta-cluster counts plus consensus decision summary.
#'
#' @return List with `pca`, optional `umap`, optional `cross_tab`, and `tables`.
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
#' p <- hc_plot_longitudinal_meta_embeddings(
#'   hc,
#'   save_pdf = FALSE,
#'   save_tables = FALSE
#' )
#' @export
hc_plot_longitudinal_meta_embeddings <- function(hc,
                                                 slot_name = "longitudinal_endotypes",
                                                 save_pdf = TRUE,
                                                 file_prefix = "Longitudinal_Meta",
                                                 show_endotype_crosstab = FALSE,
                                                 show_cluster_labels = TRUE,
                                                 save_tables = TRUE,
                                                 table_file_prefix = NULL,
                                                 table_format = c("xlsx", "csv"),
                                                 table_detail = c("simple", "full")) {
  if (!inherits(hc, "HCoCenaExperiment")) {
    stop("`hc` must be a `HCoCenaExperiment`.")
  }
  table_format <- base::match.arg(table_format)
  table_detail <- base::match.arg(table_detail)
  sat <- as.list(hc@satellite)
  obj <- sat[[slot_name]]
  if (is.null(obj) || is.null(obj$meta_cluster) || is.null(obj$meta_pca)) {
    stop("Slot `", slot_name, "` does not contain meta-clustering outputs. Run `hc_longitudinal_meta_clustering()` first.")
  }

  mc <- base::as.data.frame(obj$meta_cluster, stringsAsFactors = FALSE)
  lv <- base::sort(base::unique(base::as.character(mc$meta_cluster)))
  lv <- lv[!base::is.na(lv) & base::nzchar(lv)]
  if (base::length(lv) == 0L) {
    stop("No non-missing meta-cluster labels found in slot `", slot_name, "`.")
  }
  pal <- .hc_distinct_palette(base::max(3L, base::length(lv)))[base::seq_along(lv)]
  pal <- stats::setNames(pal, lv)

  pca_df <- base::as.data.frame(obj$meta_pca, stringsAsFactors = FALSE)
  pca_df$meta_cluster <- base::factor(pca_df$meta_cluster, levels = lv)
  pca_axis_lab <- .hc_pca_axis_labels(obj$meta_pca_variance)
  pca_ell <- .hc_filter_groups_for_ellipse(
    pca_df,
    group_col = "meta_cluster",
    min_points = 4L,
    x_col = "PC1",
    y_col = "PC2",
    min_unique = 3L
  )
  pca_groups_all <- base::sort(base::unique(base::as.character(stats::na.omit(pca_df$meta_cluster))))
  pca_groups_ell <- if (base::nrow(pca_ell) > 0) {
    base::sort(base::unique(base::as.character(pca_ell$meta_cluster)))
  } else {
    character(0)
  }
  show_pca_ell <- base::length(pca_groups_all) > 0 && base::identical(pca_groups_all, pca_groups_ell)
  if (!isTRUE(show_pca_ell)) {
    pca_ell <- pca_ell[0, , drop = FALSE]
  }
  pca_ctr <- .hc_cluster_centers_2d(pca_df, x_col = "PC1", y_col = "PC2", group_col = "meta_cluster")
  pca_lim <- .hc_square_limits_2d(pca_df, x_col = "PC1", y_col = "PC2")
  p_pca <- ggplot2::ggplot(pca_df, ggplot2::aes(x = PC1, y = PC2)) +
    ggplot2::geom_point(
      ggplot2::aes(color = meta_cluster),
      shape = 16,
      size = 3.1,
      alpha = 0.95
    ) +
    ggplot2::scale_fill_manual(values = pal, drop = FALSE) +
    ggplot2::scale_color_manual(values = pal, drop = FALSE) +
    ggplot2::coord_equal(xlim = pca_lim$x, ylim = pca_lim$y) +
    .hc_theme_pub(base_size = 11.5) +
    ggplot2::theme(
      legend.position = "right",
      aspect.ratio = 1,
      panel.border = ggplot2::element_rect(color = "grey35", linewidth = 0.55, fill = NA)
    ) +
    ggplot2::guides(
      fill = "none",
      color = ggplot2::guide_legend(
        override.aes = list(shape = 16, size = 4, alpha = 1)
      )
    ) +
    ggplot2::labs(
      title = "Meta-clusters on CAP space (PCA)",
      x = pca_axis_lab$x,
      y = pca_axis_lab$y,
      color = "Meta cluster"
    )
  if (base::nrow(pca_ell) > 0) {
    p_pca <- p_pca +
      ggplot2::stat_ellipse(
        data = pca_ell,
        mapping = ggplot2::aes(x = PC1, y = PC2, fill = meta_cluster, group = meta_cluster),
        geom = "polygon",
        type = "norm",
        level = 0.72,
        alpha = 0.08,
        color = NA,
        inherit.aes = FALSE,
        show.legend = FALSE
      ) +
      ggplot2::stat_ellipse(
        data = pca_ell,
        mapping = ggplot2::aes(x = PC1, y = PC2, color = meta_cluster, group = meta_cluster),
        type = "norm",
        level = 0.72,
        alpha = 0.95,
        linewidth = 0.9,
        inherit.aes = FALSE,
        show.legend = FALSE
      )
  }
  if (isTRUE(show_cluster_labels) && base::nrow(pca_ctr) > 0) {
    p_pca <- p_pca +
      ggplot2::geom_label(
        data = pca_ctr,
        mapping = ggplot2::aes(x = x, y = y, label = group, fill = group),
        color = "black",
        fontface = "bold",
        size = 3.3,
        label.padding = grid::unit(0.12, "lines"),
        label.r = grid::unit(0.08, "lines"),
        label.size = 0.25,
        alpha = 0.9,
        inherit.aes = FALSE,
        show.legend = FALSE
      )
  }
  p_pca <- .hc_attach_right_legend(p_pca)

  p_umap <- NULL
  if (!is.null(obj$meta_umap)) {
    um_df <- base::as.data.frame(obj$meta_umap, stringsAsFactors = FALSE)
    um_df$meta_cluster <- base::factor(um_df$meta_cluster, levels = lv)
    um_ell <- .hc_filter_groups_for_ellipse(
      um_df,
      group_col = "meta_cluster",
      min_points = 4L,
      x_col = "UMAP1",
      y_col = "UMAP2",
      min_unique = 3L
    )
    um_groups_all <- base::sort(base::unique(base::as.character(stats::na.omit(um_df$meta_cluster))))
    um_groups_ell <- if (base::nrow(um_ell) > 0) {
      base::sort(base::unique(base::as.character(um_ell$meta_cluster)))
    } else {
      character(0)
    }
    show_um_ell <- base::length(um_groups_all) > 0 && base::identical(um_groups_all, um_groups_ell)
    if (!isTRUE(show_um_ell)) {
      um_ell <- um_ell[0, , drop = FALSE]
    }
    um_ctr <- .hc_cluster_centers_2d(um_df, x_col = "UMAP1", y_col = "UMAP2", group_col = "meta_cluster")
    um_lim <- .hc_square_limits_2d(um_df, x_col = "UMAP1", y_col = "UMAP2")
    p_umap <- ggplot2::ggplot(um_df, ggplot2::aes(x = UMAP1, y = UMAP2)) +
      ggplot2::geom_point(
        ggplot2::aes(color = meta_cluster),
        shape = 16,
        size = 3.1,
        alpha = 0.95
      ) +
      ggplot2::scale_fill_manual(values = pal, drop = FALSE) +
      ggplot2::scale_color_manual(values = pal, drop = FALSE) +
      ggplot2::coord_equal(xlim = um_lim$x, ylim = um_lim$y) +
      .hc_theme_pub(base_size = 11.5) +
      ggplot2::theme(
        legend.position = "right",
        aspect.ratio = 1,
        panel.border = ggplot2::element_rect(color = "grey35", linewidth = 0.55, fill = NA)
      ) +
      ggplot2::guides(
        fill = "none",
        color = ggplot2::guide_legend(
          override.aes = list(shape = 16, size = 4, alpha = 1)
        )
      ) +
      ggplot2::labs(
        title = "Meta-clusters on CAP space (UMAP)",
        x = "UMAP1",
        y = "UMAP2",
        color = "Meta cluster"
      )
    if (base::nrow(um_ell) > 0) {
      p_umap <- p_umap +
        ggplot2::stat_ellipse(
          data = um_ell,
          mapping = ggplot2::aes(x = UMAP1, y = UMAP2, fill = meta_cluster, group = meta_cluster),
          geom = "polygon",
          type = "norm",
          level = 0.72,
          alpha = 0.08,
          color = NA,
          inherit.aes = FALSE,
          show.legend = FALSE
        ) +
        ggplot2::stat_ellipse(
          data = um_ell,
          mapping = ggplot2::aes(x = UMAP1, y = UMAP2, color = meta_cluster, group = meta_cluster),
          type = "norm",
          level = 0.72,
          alpha = 0.95,
          linewidth = 0.9,
          inherit.aes = FALSE,
          show.legend = FALSE
        )
    }
    if (isTRUE(show_cluster_labels) && base::nrow(um_ctr) > 0) {
      p_umap <- p_umap +
        ggplot2::geom_label(
          data = um_ctr,
          mapping = ggplot2::aes(x = x, y = y, label = group, fill = group),
          color = "black",
          fontface = "bold",
          size = 3.3,
          label.padding = grid::unit(0.12, "lines"),
          label.r = grid::unit(0.08, "lines"),
          label.size = 0.25,
          alpha = 0.9,
          inherit.aes = FALSE,
          show.legend = FALSE
        )
    }
    p_umap <- .hc_attach_right_legend(p_umap)
  }

  p_cross <- NULL
  if (isTRUE(show_endotype_crosstab) &&
    !is.null(obj$donor_cluster) &&
    "endotype" %in% base::colnames(obj$donor_cluster)) {
    tmp <- base::merge(
      obj$donor_cluster[, c("donor", "endotype"), drop = FALSE],
      mc[, c("donor", "meta_cluster"), drop = FALSE],
      by = "donor",
      all = FALSE,
      sort = FALSE
    )
    if (base::nrow(tmp) > 0) {
      tab <- base::as.data.frame(base::table(tmp$endotype, tmp$meta_cluster), stringsAsFactors = FALSE)
      base::colnames(tab) <- c("endotype", "meta_cluster", "n")
      p_cross <- ggplot2::ggplot(tab, ggplot2::aes(x = endotype, y = meta_cluster, fill = n)) +
        ggplot2::geom_tile(color = "white", linewidth = 0.3) +
        ggplot2::geom_text(ggplot2::aes(label = n), size = 3) +
        ggplot2::scale_fill_gradient(low = "grey95", high = "#1B9E77") +
        .hc_theme_pub(base_size = 11) +
        ggplot2::labs(
          title = "Endotype vs meta-cluster overlap",
          x = "Endotype",
          y = "Meta cluster",
          fill = "Donors"
        )
    }
  }

  meta_tables <- .hc_meta_report_tables(obj, detail = table_detail)
  meta_tables <- meta_tables[
    base::vapply(
      meta_tables,
      function(x) base::is.data.frame(x) && !base::is.null(x),
      logical(1)
    )
  ]

  out_dir <- NULL
  if (isTRUE(save_pdf) || isTRUE(save_tables)) {
    out_dir <- .hc_resolve_output_dir(hc)
  }

  if (isTRUE(save_pdf)) {
    .hc_ggsave_pdf_png(
      filename = base::file.path(out_dir, base::paste0(file_prefix, "_PCA.pdf")),
      plot = p_pca,
      width = 7,
      height = 5,
      units = "in",
      device = grDevices::cairo_pdf
    )
    if (!is.null(p_umap)) {
      .hc_ggsave_pdf_png(
        filename = base::file.path(out_dir, base::paste0(file_prefix, "_UMAP.pdf")),
        plot = p_umap,
        width = 7,
        height = 5,
        units = "in",
        device = grDevices::cairo_pdf
      )
    }
    if (!is.null(p_cross)) {
      .hc_ggsave_pdf_png(
        filename = base::file.path(out_dir, base::paste0(file_prefix, "_CrossTab.pdf")),
        plot = p_cross,
        width = 6.5,
        height = 4.8,
        units = "in",
        device = grDevices::cairo_pdf
      )
    }
  }

  if (isTRUE(save_tables) && base::length(meta_tables) > 0) {
    if (is.null(out_dir)) {
      out_dir <- .hc_resolve_output_dir(hc)
    }
    tbl_prefix <- if (base::is.null(table_file_prefix) || !base::nzchar(base::as.character(table_file_prefix[[1]]))) {
      file_prefix
    } else {
      base::as.character(table_file_prefix[[1]])
    }
    if (identical(table_format, "xlsx") && requireNamespace("openxlsx", quietly = TRUE)) {
      .hc_write_xlsx_atomic(
        x = meta_tables,
        file = base::file.path(out_dir, base::paste0(tbl_prefix, "_Meta_Clustering_Summary.xlsx")),
        overwrite = TRUE
      )
    } else {
      if (identical(table_format, "xlsx") && !requireNamespace("openxlsx", quietly = TRUE)) {
        warning("Package `openxlsx` not installed. Falling back to CSV for meta summary tables.")
      }
      for (nm in base::names(meta_tables)) {
        utils::write.csv(
          x = meta_tables[[nm]],
          file = base::file.path(out_dir, base::paste0(tbl_prefix, "_", nm, ".csv")),
          row.names = FALSE
        )
      }
    }
  }

  list(pca = p_pca, umap = p_umap, cross_tab = p_cross, tables = meta_tables)
}

#' Plot consensus diagnostics for longitudinal meta-clustering
#'
#' Creates three consensus diagnostics:
#' 1) donor-by-donor consensus matrix heatmap (selected `k`)
#' 2) consensus metric trajectories across tested `k` (delta/PAC/coverage)
#' 3) donor stability distribution per meta-cluster.
#'
#' @param hc A `HCoCenaExperiment`.
#' @param slot_name Satellite slot with meta-clustering outputs.
#' @param save_pdf Logical. If `TRUE`, save plots as PDFs.
#' @param file_prefix Output filename prefix.
#' @param show_donor_labels Logical. If `TRUE`, print donor IDs on consensus heatmap axes.
#' @param save_table Logical. If `TRUE`, export consensus score table.
#' @param table_format One of `"xlsx"` or `"csv"`.
#'
#' @return A list with `consensus_heatmap`, optional `consensus_k`,
#'   optional `stability`, and optional `summary_table`.
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
#' p <- hc_plot_longitudinal_meta_consensus(
#'   hc,
#'   save_pdf = FALSE,
#'   save_table = FALSE
#' )
#' @export
hc_plot_longitudinal_meta_consensus <- function(hc,
                                                slot_name = "longitudinal_endotypes",
                                                save_pdf = TRUE,
                                                file_prefix = "Longitudinal_Meta_Consensus",
                                                show_donor_labels = FALSE,
                                                save_table = TRUE,
                                                table_format = c("xlsx", "csv")) {
  if (!inherits(hc, "HCoCenaExperiment")) {
    stop("`hc` must be a `HCoCenaExperiment`.")
  }
  table_format <- base::match.arg(table_format)
  sat <- as.list(hc@satellite)
  obj <- sat[[slot_name]]
  if (is.null(obj) || !isTRUE(obj$meta_consensus) || is.null(obj$meta_consensus_matrix)) {
    stop("Slot `", slot_name, "` has no consensus matrix. Run `hc_longitudinal_meta_clustering(..., consensus = TRUE)` first.")
  }
  if (is.null(obj$meta_cluster)) {
    stop("Slot `", slot_name, "` has no `meta_cluster` table.")
  }

  cons <- base::as.matrix(obj$meta_consensus_matrix)
  if (base::nrow(cons) == 0 || base::ncol(cons) == 0) {
    stop("Consensus matrix is empty.")
  }

  meta_cluster <- base::as.data.frame(obj$meta_cluster, stringsAsFactors = FALSE)
  if (!"donor" %in% base::colnames(meta_cluster)) {
    stop("`meta_cluster` table must contain a `donor` column.")
  }
  meta_cluster$donor <- base::as.character(meta_cluster$donor)
  meta_cluster <- meta_cluster[!base::is.na(meta_cluster$donor) & base::nzchar(meta_cluster$donor), , drop = FALSE]
  meta_cluster <- meta_cluster[!base::duplicated(meta_cluster$donor), , drop = FALSE]
  donors_meta <- base::as.character(meta_cluster$donor)

  if (base::is.null(base::rownames(cons)) || base::is.null(base::colnames(cons))) {
    if (base::length(donors_meta) == base::nrow(cons)) {
      base::rownames(cons) <- donors_meta
      base::colnames(cons) <- donors_meta
    } else {
      auto_ids <- base::paste0("D", base::seq_len(base::nrow(cons)))
      base::rownames(cons) <- auto_ids
      base::colnames(cons) <- auto_ids
    }
  }

  # Defensive clean-up for duplicate or mismatched row/column donor IDs.
  if (base::is.null(base::rownames(cons)) || base::is.null(base::colnames(cons))) {
    stop("Consensus matrix must have row and column names.")
  }
  cons <- cons[!base::duplicated(base::rownames(cons)), , drop = FALSE]
  cons <- cons[, !base::duplicated(base::colnames(cons)), drop = FALSE]
  donors_cons <- base::intersect(base::rownames(cons), base::colnames(cons))
  donors <- base::intersect(donors_cons, donors_meta)
  if (base::length(donors) < 3) {
    stop("Too few overlapping donors between consensus matrix and meta-cluster table.")
  }

  cons <- cons[donors, donors, drop = FALSE]
  cons[!base::is.finite(cons)] <- 0
  cons <- (cons + base::t(cons)) / 2
  diag(cons) <- 1

  meta_cluster <- meta_cluster[base::match(donors, meta_cluster$donor), , drop = FALSE]
  meta_cluster$meta_cluster <- base::as.character(meta_cluster$meta_cluster)
  lv <- base::sort(base::unique(meta_cluster$meta_cluster))
  pal <- .hc_distinct_palette(base::max(3L, base::length(lv)))[base::seq_along(lv)]
  pal <- stats::setNames(pal, lv)

  ord_idx <- integer(0)
  for (cl in lv) {
    idx <- base::which(meta_cluster$meta_cluster == cl)
    if (base::length(idx) <= 2) {
      ord_idx <- base::c(ord_idx, idx)
    } else {
      subc <- cons[idx, idx, drop = FALSE]
      subc[!base::is.finite(subc)] <- 0
      subc <- (subc + base::t(subc)) / 2
      dmat <- base::pmax(1 - subc, 0)
      dmat <- base::as.matrix(dmat)
      diag(dmat) <- 0
      local_order <- tryCatch(
        {
          d <- stats::as.dist(dmat)
          expected_len <- (base::nrow(dmat) * (base::nrow(dmat) - 1)) / 2
          if (base::length(d) != expected_len) {
            stop("Invalid dist length for within-cluster ordering.")
          }
          stats::hclust(d, method = "average")$order
        },
        error = function(e) base::seq_len(base::length(idx))
      )
      ord_idx <- base::c(ord_idx, idx[local_order])
    }
  }
  if (base::length(ord_idx) == 0) {
    ord_idx <- base::seq_along(donors)
  }
  donor_order <- donors[ord_idx]
  cons_ord <- cons[donor_order, donor_order, drop = FALSE]
  cluster_order <- meta_cluster$meta_cluster[base::match(donor_order, meta_cluster$donor)]

  hm_df <- base::as.data.frame(base::as.table(cons_ord), stringsAsFactors = FALSE)
  base::colnames(hm_df) <- c("donor_row", "donor_col", "consensus")
  hm_df$donor_row <- base::factor(base::as.character(hm_df$donor_row), levels = donor_order)
  hm_df$donor_col <- base::factor(base::as.character(hm_df$donor_col), levels = donor_order)
  hm_df$consensus <- as.numeric(hm_df$consensus)

  p_heat <- ggplot2::ggplot(hm_df, ggplot2::aes(x = donor_col, y = donor_row, fill = consensus)) +
    ggplot2::geom_tile(color = NA) +
    ggplot2::coord_fixed() +
    ggplot2::scale_fill_gradientn(
      colors = c("#F7FBFF", "#C6DBEF", "#6BAED6", "#2171B5", "#08306B"),
      limits = c(0, 1),
      oob = function(x, ...) base::pmin(base::pmax(x, 0), 1)
    ) +
    .hc_theme_pub(base_size = 10.5) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      axis.text.x = if (isTRUE(show_donor_labels)) ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5) else ggplot2::element_blank(),
      axis.text.y = if (isTRUE(show_donor_labels)) ggplot2::element_text(size = 5) else ggplot2::element_blank(),
      axis.ticks = if (isTRUE(show_donor_labels)) ggplot2::element_line() else ggplot2::element_blank(),
      legend.position = "right"
    ) +
    ggplot2::labs(
      title = "Consensus co-assignment matrix (selected meta-k)",
      subtitle = base::paste0(
        "k = ", as.integer(obj$meta_best_k[[1]]),
        " | runs = ", as.integer(obj$meta_consensus_runs[[1]]),
        " | sample frac = ", base::formatC(as.numeric(obj$meta_consensus_sample_fraction[[1]]), format = "f", digits = 2),
        " | feature frac = ", base::formatC(as.numeric(obj$meta_consensus_feature_fraction[[1]]), format = "f", digits = 2)
      ),
      fill = "Consensus"
    )

  rl <- base::rle(cluster_order)
  if (base::length(rl$lengths) > 1) {
    b <- base::cumsum(rl$lengths)
    b <- b[-base::length(b)]
    if (base::length(b) > 0) {
      p_heat <- p_heat +
        ggplot2::geom_vline(xintercept = b + 0.5, color = "white", linewidth = 0.45) +
        ggplot2::geom_hline(yintercept = b + 0.5, color = "white", linewidth = 0.45)
    }
  }

  p_k <- NULL
  summary_tbl <- NULL
  if (!is.null(obj$meta_score_table)) {
    sc <- base::as.data.frame(obj$meta_score_table, stringsAsFactors = FALSE)
    if ("k" %in% base::colnames(sc)) {
      sc$selected <- as.integer(sc$k) == as.integer(obj$meta_best_k[[1]])
      summary_tbl <- sc

      metric_cols <- c("consensus_delta", "pac", "observed_pair_fraction")
      metric_cols <- metric_cols[metric_cols %in% base::colnames(sc)]
      if (base::length(metric_cols) > 0) {
        long_list <- base::lapply(
          metric_cols,
          function(m) {
            base::data.frame(
              k = as.numeric(sc$k),
              metric = m,
              value = as.numeric(sc[[m]]),
              selected = sc$selected,
              stringsAsFactors = FALSE
            )
          }
        )
        k_long <- base::do.call(base::rbind, long_list)
        metric_labels <- c(
          consensus_delta = "Delta (within - between, higher better)",
          pac = "PAC (lower better)",
          observed_pair_fraction = "Observed pair fraction (higher better)"
        )
        k_long$metric <- base::factor(
          k_long$metric,
          levels = metric_cols,
          labels = metric_labels[metric_cols]
        )

        p_k <- ggplot2::ggplot(k_long, ggplot2::aes(x = k, y = value)) +
          ggplot2::geom_line(color = "#2F4F4F", linewidth = 0.85) +
          ggplot2::geom_point(
            ggplot2::aes(color = selected, size = selected),
            show.legend = FALSE
          ) +
          ggplot2::geom_vline(
            xintercept = as.numeric(obj$meta_best_k[[1]]),
            linetype = "dashed",
            color = "#B22222",
            linewidth = 0.7
          ) +
          ggplot2::scale_color_manual(values = c("FALSE" = "#2F4F4F", "TRUE" = "#B22222")) +
          ggplot2::scale_size_manual(values = c("FALSE" = 2.1, "TRUE" = 3.1)) +
          ggplot2::facet_wrap(~metric, ncol = 1, scales = "free_y") +
          .hc_theme_pub(base_size = 10.5) +
          ggplot2::labs(
            title = "Consensus diagnostics across candidate k",
            x = "k",
            y = NULL
          )
      }
    }
  }

  p_stability <- NULL
  if ("stability" %in% base::colnames(meta_cluster)) {
    stb <- base::data.frame(
      donor = base::as.character(meta_cluster$donor),
      meta_cluster = base::factor(base::as.character(meta_cluster$meta_cluster), levels = lv),
      stability = .hc_as_numeric_safely(meta_cluster$stability),
      stringsAsFactors = FALSE
    )
    stb <- stb[base::is.finite(stb$stability), , drop = FALSE]
    if (base::nrow(stb) > 0) {
      p_stability <- ggplot2::ggplot(stb, ggplot2::aes(x = meta_cluster, y = stability, fill = meta_cluster)) +
        ggplot2::geom_boxplot(width = 0.62, alpha = 0.40, outlier.shape = NA, color = "black", linewidth = 0.45) +
        ggplot2::geom_jitter(shape = 21, color = "black", size = 2.2, alpha = 0.90, width = 0.12, height = 0) +
        ggplot2::scale_fill_manual(values = pal, drop = FALSE) +
        .hc_theme_pub(base_size = 11) +
        ggplot2::theme(legend.position = "none") +
        ggplot2::labs(
          title = "Donor consensus stability by meta-cluster",
          x = "Meta cluster",
          y = "Stability"
        )
    }
  }

  out_dir <- NULL
  if (isTRUE(save_pdf) || isTRUE(save_table)) {
    out_dir <- .hc_resolve_output_dir(hc)
  }

  if (isTRUE(save_pdf)) {
    .hc_ggsave_pdf_png(
      filename = base::file.path(out_dir, base::paste0(file_prefix, "_Matrix.pdf")),
      plot = p_heat,
      width = 7.5,
      height = 6.5,
      units = "in",
      device = grDevices::cairo_pdf
    )
    if (!is.null(p_k)) {
      .hc_ggsave_pdf_png(
        filename = base::file.path(out_dir, base::paste0(file_prefix, "_KDiagnostics.pdf")),
        plot = p_k,
        width = 7.2,
        height = 7.2,
        units = "in",
        device = grDevices::cairo_pdf
      )
    }
    if (!is.null(p_stability)) {
      .hc_ggsave_pdf_png(
        filename = base::file.path(out_dir, base::paste0(file_prefix, "_Stability.pdf")),
        plot = p_stability,
        width = 6.4,
        height = 4.8,
        units = "in",
        device = grDevices::cairo_pdf
      )
    }
  }

  if (isTRUE(save_table) && !is.null(summary_tbl)) {
    if (is.null(out_dir)) {
      out_dir <- .hc_resolve_output_dir(hc)
    }
    if (identical(table_format, "xlsx") && requireNamespace("openxlsx", quietly = TRUE)) {
      .hc_write_xlsx_atomic(
        x = list(Consensus_K_Table = summary_tbl),
        file = base::file.path(out_dir, base::paste0(file_prefix, "_Summary.xlsx")),
        overwrite = TRUE
      )
    } else {
      if (identical(table_format, "xlsx") && !requireNamespace("openxlsx", quietly = TRUE)) {
        warning("Package `openxlsx` not installed. Falling back to CSV for consensus summary table.")
      }
      utils::write.csv(
        x = summary_tbl,
        file = base::file.path(out_dir, base::paste0(file_prefix, "_Summary.csv")),
        row.names = FALSE
      )
    }
  }

  list(
    consensus_heatmap = p_heat,
    consensus_k = p_k,
    stability = p_stability,
    summary_table = summary_tbl
  )
}

#' Plot module trajectories by donor meta-cluster
#'
#' Creates module-wise wave plots after meta-clustering: faint donor trajectories
#' plus bold meta-cluster mean trajectories.
#'
#' @param hc A `HCoCenaExperiment`.
#' @param slot_name Satellite slot with meta-clustering outputs.
#' @param save_pdf Logical. If `TRUE`, save PDF.
#' @param file_prefix Output filename prefix.
#' @param donor_alpha Alpha for donor-level trajectories.
#' @param donor_linewidth Line width for donor-level trajectories.
#' @param mean_linewidth Line width for meta-cluster mean lines.
#' @param mean_point_size Point size for meta-cluster mean points.
#' @param facet_ncol Number of facet columns.
#' @param free_y Logical. If `TRUE`, use free y-axis scales per module.
#' @param square_panels Logical. If `TRUE`, use square facet panels.
#' @param value_mode One of `"scaled_mean_vst"` or `"expression"`. The default
#'   scales module values within each module before plotting.
#' @param value_range Optional numeric vector of length 2 giving the y-axis
#'   range for `value_mode = "scaled_mean_vst"`. Defaults to `c(-2, 2)`. Use
#'   `NULL` to keep the automatic range.
#' @param save_width,save_height Optional PDF width/height in inches.
#'
#' @return A list with one ggplot object (`meta_module_waves`).
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
#' p <- hc_plot_longitudinal_meta_module_waves(hc, save_pdf = FALSE)
#' @export
hc_plot_longitudinal_meta_module_waves <- function(hc,
                                                   slot_name = "longitudinal_endotypes",
                                                   save_pdf = TRUE,
                                                   file_prefix = "Longitudinal_Meta_ModuleWaves",
                                                   donor_alpha = 0.13,
                                                   donor_linewidth = 0.28,
                                                   mean_linewidth = 1.2,
                                                   mean_point_size = 1.9,
                                                   facet_ncol = 4,
                                                   free_y = TRUE,
                                                   square_panels = TRUE,
                                                   value_mode = c("scaled_mean_vst", "expression"),
                                                   value_range = c(-2, 2),
                                                   save_width = NULL,
                                                   save_height = NULL) {
  if (!inherits(hc, "HCoCenaExperiment")) {
    stop("`hc` must be a `HCoCenaExperiment`.")
  }
  value_mode <- base::match.arg(value_mode)
  if (!is.null(value_range)) {
    value_range <- .hc_as_numeric_safely(value_range)
    if (base::length(value_range) != 2 || !base::all(base::is.finite(value_range))) {
      stop("`value_range` must be NULL or a numeric vector of length 2.")
    }
    value_range <- sort(value_range)
  }
  sat <- as.list(hc@satellite)
  obj <- sat[[slot_name]]
  if (is.null(obj) || is.null(obj$meta_cluster) || is.null(obj$donor_time_module)) {
    stop("Slot `", slot_name, "` does not contain required meta-clustering trajectory data.")
  }

  tr <- base::as.data.frame(obj$donor_time_module, stringsAsFactors = FALSE)
  mc <- base::as.data.frame(obj$meta_cluster, stringsAsFactors = FALSE)
  tr <- base::merge(
    tr,
    mc[, c("donor", "meta_cluster"), drop = FALSE],
    by = "donor",
    all = FALSE,
    sort = FALSE
  )
  if (base::nrow(tr) == 0) {
    stop("No overlap between trajectory donors and meta-cluster donors.")
  }

  if (!is.null(obj$time_levels)) {
    tr$time <- base::factor(base::as.character(tr$time), levels = base::as.character(obj$time_levels), ordered = TRUE)
  }

  module_levels <- base::sort(base::unique(base::as.character(tr$module)))
  module_lookup <- NULL
  if (!is.null(obj$module_lookup) && "module" %in% base::colnames(obj$module_lookup)) {
    module_lookup <- base::as.data.frame(obj$module_lookup, stringsAsFactors = FALSE)
    module_lookup <- .hc_reorder_module_lookup_natural(module_lookup)
    ml <- base::as.character(module_lookup$module)
    module_levels <- ml[ml %in% module_levels]
  }
  if (is.null(module_lookup) || base::nrow(module_lookup) == 0) {
    module_lookup <- base::data.frame(
      module = module_levels,
      module_color = "grey60",
      stringsAsFactors = FALSE
    )
  }
  mod_col <- .hc_module_color_map(
    module_lookup = module_lookup,
    module_levels = module_levels
  )

  tr$module <- base::factor(base::as.character(tr$module), levels = module_levels)
  tr$meta_cluster <- base::as.character(tr$meta_cluster)
  tr$donor <- base::as.character(tr$donor)

  tr$value_plot <- as.numeric(tr$value)
  if (identical(value_mode, "scaled_mean_vst")) {
    split_idx <- base::split(base::seq_len(base::nrow(tr)), base::as.character(tr$module))
    for (idx in split_idx) {
      vv <- .hc_as_numeric_safely(tr$value_plot[idx])
      if (base::length(vv) == 0) {
        next
      }
      sd_v <- stats::sd(vv, na.rm = TRUE)
      if (is.finite(sd_v) && sd_v > 0) {
        tr$value_plot[idx] <- as.numeric(base::scale(vv))
      } else {
        tr$value_plot[idx] <- 0
      }
    }
  }

  mean_df <- stats::aggregate(value_plot ~ module + time + meta_cluster, data = tr, FUN = base::mean)
  mc_levels <- base::sort(base::unique(base::as.character(tr$meta_cluster)))
  tr$meta_cluster <- base::factor(tr$meta_cluster, levels = mc_levels)
  mean_df$meta_cluster <- base::factor(base::as.character(mean_df$meta_cluster), levels = mc_levels)

  pal <- grDevices::hcl.colors(base::max(3, base::length(mc_levels)), palette = "Dark 3")[base::seq_along(mc_levels)]
  pal <- stats::setNames(pal, mc_levels)

  y_lab <- if (identical(value_mode, "scaled_mean_vst")) {
    "Scaled mean variance-stabilized expression"
  } else {
    if (!is.null(obj$value_label) && base::nzchar(as.character(obj$value_label[[1]]))) {
      as.character(obj$value_label[[1]])
    } else {
      "Module mean expression"
    }
  }

  p <- ggplot2::ggplot(
    tr,
    ggplot2::aes(x = time, y = value_plot, group = donor, color = meta_cluster)
  ) +
    ggplot2::geom_line(alpha = donor_alpha, linewidth = donor_linewidth, show.legend = FALSE) +
    ggplot2::geom_line(
      data = mean_df,
      mapping = ggplot2::aes(y = value_plot, group = meta_cluster),
      linewidth = mean_linewidth
    ) +
    ggplot2::geom_point(
      data = mean_df,
      mapping = ggplot2::aes(y = value_plot, group = meta_cluster),
      size = mean_point_size
    ) +
    .hc_module_facet_wrap(
      module_levels = module_levels,
      mod_col = mod_col,
      ncol = facet_ncol,
      scales = if (isTRUE(free_y)) "free_y" else "fixed",
      facets = ~module
    ) +
    ggplot2::scale_color_manual(values = pal, drop = FALSE) +
    .hc_theme_pub(base_size = 11) +
    ggplot2::theme(
      aspect.ratio = if (isTRUE(square_panels)) 1 else NULL,
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    ) +
    ggplot2::labs(
      title = "Module trajectories by meta-cluster",
      x = "Timepoint",
      y = y_lab,
      color = "Meta cluster"
    )
  if (identical(value_mode, "scaled_mean_vst") && !is.null(value_range)) {
    p <- p + ggplot2::coord_cartesian(ylim = value_range)
  }
  p <- .hc_colorize_module_strips(
    p = p,
    module_levels = module_levels,
    mod_col = mod_col
  )

  if (isTRUE(save_pdf)) {
    out_dir <- .hc_resolve_output_dir(hc)
    w <- if (is.null(save_width)) {
      if (isTRUE(square_panels)) 10 else 12
    } else {
      as.numeric(save_width[[1]])
    }
    h <- if (is.null(save_height)) {
      if (isTRUE(square_panels)) w else 8
    } else {
      as.numeric(save_height[[1]])
    }
    .hc_ggsave_pdf_png(
      filename = base::file.path(out_dir, base::paste0(file_prefix, ".pdf")),
      plot = p,
      width = w,
      height = h,
      units = "in",
      device = grDevices::cairo_pdf
    )
  }

  list(meta_module_waves = p)
}

#' Normalize enrichment summary table shape
#' @noRd
.hc_normalize_enrichment_summary <- function(df) {
  if (base::is.null(df) || !base::is.data.frame(df) || base::nrow(df) == 0) {
    return(NULL)
  }
  out <- base::as.data.frame(df, stringsAsFactors = FALSE)

  if (!"database" %in% base::colnames(out) && "DB" %in% base::colnames(out)) {
    out$database <- base::as.character(out$DB)
  }
  if (!"cluster" %in% base::colnames(out) && "color" %in% base::colnames(out)) {
    out$cluster <- base::as.character(out$color)
  }
  if (!"term" %in% base::colnames(out) && "Description" %in% base::colnames(out)) {
    out$term <- base::as.character(out$Description)
  }
  if (!"qvalue" %in% base::colnames(out)) {
    if ("p.adjust" %in% base::colnames(out)) {
      out$qvalue <- .hc_as_numeric_safely(out$p.adjust)
    } else if ("pvalue" %in% base::colnames(out)) {
      out$qvalue <- .hc_as_numeric_safely(out$pvalue)
    } else {
      out$qvalue <- NA_real_
    }
  }
  if (!"rank" %in% base::colnames(out)) {
    out$rank <- NA_integer_
  }
  if (!"module_label" %in% base::colnames(out) && "module" %in% base::colnames(out)) {
    out$module_label <- base::as.character(out$module)
  }
  if (!"geneID" %in% base::colnames(out)) {
    out$geneID <- NA_character_
  }

  need <- c("database", "cluster", "term")
  if (!base::all(need %in% base::colnames(out))) {
    return(NULL)
  }

  out$database <- base::as.character(out$database)
  out$cluster <- base::as.character(out$cluster)
  out$term <- base::as.character(out$term)
  out$module_label <- base::as.character(out$module_label)
  out$rank <- .hc_as_numeric_safely(out$rank)
  out$qvalue <- .hc_as_numeric_safely(out$qvalue)
  out$geneID <- base::as.character(out$geneID)

  keep <- !base::is.na(out$database) &
    !base::is.na(out$cluster) &
    !base::is.na(out$term) &
    out$database != "" &
    out$cluster != "" &
    out$term != ""
  out <- out[keep, , drop = FALSE]
  if (base::nrow(out) == 0) {
    return(NULL)
  }

  cols <- c("database", "cluster", "module_label", "term", "rank", "qvalue", "geneID")
  miss <- cols[!cols %in% base::colnames(out)]
  if (base::length(miss) > 0) {
    for (nm in miss) {
      out[[nm]] <- NA
    }
  }
  out <- out[, cols, drop = FALSE]
  base::rownames(out) <- NULL
  out
}

#' Resolve enrichment summary table from satellite storage
#' @noRd
.hc_collect_enrichment_summary <- function(hc,
                                           enrichment_table = c("selected", "all", "significant")) {
  enrichment_table <- base::match.arg(enrichment_table)
  sat <- as.list(hc@satellite)
  store <- sat[["enrichments"]]
  if (base::is.null(store) || !base::is.list(store)) {
    store <- list()
  }

  table_key <- switch(enrichment_table,
    selected = "selected_enrichments_all_dbs",
    all = "all_enrichments_all_dbs",
    significant = "significant_enrichments_all_dbs"
  )

  collected <- list()

  direct <- .hc_normalize_enrichment_summary(store[[table_key]])
  if (base::is.null(direct)) {
    direct <- .hc_normalize_enrichment_summary(sat[[table_key]])
  }
  if (!base::is.null(direct) && base::nrow(direct) > 0) {
    collected[[base::length(collected) + 1L]] <- direct
  }

  if (identical(enrichment_table, "selected") && base::is.list(store[["top_all_dbs"]])) {
    top_all <- .hc_normalize_enrichment_summary(store[["top_all_dbs"]][["result"]])
    if (!base::is.null(top_all) && base::nrow(top_all) > 0) {
      collected[[base::length(collected) + 1L]] <- top_all
    }
  }

  subkey <- switch(enrichment_table,
    selected = "selected_enrichments",
    all = "all_enrichments",
    significant = "significant_enrichments"
  )
  per_db_keys <- base::grep("^top_", base::names(store), value = TRUE)
  per_db_keys <- per_db_keys[!per_db_keys %in% c("top_all_dbs", "top_all_dbs_mixed")]
  if (base::length(per_db_keys) == 0) {
    return(NULL)
  }

  collected <- base::lapply(per_db_keys, function(k) {
    ent <- store[[k]]
    if (base::is.null(ent) || !base::is.list(ent)) {
      return(NULL)
    }
    .hc_normalize_enrichment_summary(ent[[subkey]])
  })
  collected <- collected[!base::vapply(collected, base::is.null, FUN.VALUE = base::logical(1))]
  if (base::length(collected) == 0) {
    return(NULL)
  }
  out <- base::do.call(base::rbind, collected)
  out <- out[!base::duplicated(base::do.call(base::paste, c(out, sep = "\r"))), , drop = FALSE]
  base::rownames(out) <- NULL
  out
}

#' Build safe wrapped term labels for longitudinal enrichment plots
#' @noRd
.hc_safe_longitudinal_term_label <- function(term, wrap_term) {
  term <- base::as.character(term)
  base::vapply(term, function(x) {
    term_clean <- .hc_clean_enrichment_term_label(x)
    label <- tryCatch(
      wrap_term(term_clean),
      error = function(e) NA_character_
    )
    label <- base::as.character(label)[1]
    if (base::is.na(label) || !base::nzchar(base::trimws(label))) {
      label <- base::as.character(term_clean)[1]
    }
    if (base::is.na(label) || !base::nzchar(base::trimws(label))) {
      label <- base::as.character(x)[1]
    }
    if (base::is.na(label) || !base::nzchar(base::trimws(label))) {
      label <- "Term unavailable"
    }
    label
  }, FUN.VALUE = base::character(1))
}

#' Build safe panel labels for longitudinal enrichment plots
#' @noRd
.hc_longitudinal_panel_label <- function(module, rank_label, term_display) {
  module <- base::as.character(module)
  rank_label <- base::as.character(rank_label)
  term_display <- base::as.character(term_display)

  module[base::is.na(module) | !base::nzchar(base::trimws(module))] <- "Module"
  rank_label[base::is.na(rank_label) | !base::nzchar(base::trimws(rank_label))] <- "Top"
  term_display[base::is.na(term_display) | !base::nzchar(base::trimws(term_display))] <- "Term unavailable"

  base::paste0(module, ": ", rank_label, "\n", term_display)
}

#' Rebuild contiguous per-module top-term ranks after downstream filtering
#' @noRd
.hc_longitudinal_rank_lookup <- function(top_df,
                                         wrap_term,
                                         valid_terms = NULL) {
  out_cols <- c("module", "term", "term_rank", "term_rank_label", "term_label", "module_color", "rank", "qvalue")
  empty <- stats::setNames(
    base::vector("list", base::length(out_cols)),
    out_cols
  )
  empty$module <- base::character(0)
  empty$term <- base::character(0)
  empty$term_rank <- base::integer(0)
  empty$term_rank_label <- base::character(0)
  empty$term_label <- base::character(0)
  empty$module_color <- base::character(0)
  empty$rank <- base::numeric(0)
  empty$qvalue <- base::numeric(0)
  empty <- base::as.data.frame(empty, stringsAsFactors = FALSE)

  top_df <- base::as.data.frame(top_df, stringsAsFactors = FALSE)
  if (base::nrow(top_df) == 0) {
    return(empty)
  }

  miss <- base::setdiff(c("module", "term", "module_color", "rank", "qvalue"), base::colnames(top_df))
  if (base::length(miss) > 0) {
    for (nm in miss) {
      top_df[[nm]] <- NA
    }
  }

  top_df$module <- base::as.character(top_df$module)
  top_df$term <- base::as.character(top_df$term)
  top_df$module_color <- base::as.character(top_df$module_color)
  top_df <- top_df[
    !base::is.na(top_df$module) & base::nzchar(base::trimws(top_df$module)) &
      !base::is.na(top_df$term) & base::nzchar(base::trimws(top_df$term)), ,
    drop = FALSE
  ]
  if (base::nrow(top_df) == 0) {
    return(empty)
  }

  if (!base::is.null(valid_terms)) {
    valid_terms <- base::as.data.frame(valid_terms, stringsAsFactors = FALSE)
    if (base::nrow(valid_terms) == 0 || !base::all(c("module", "term") %in% base::colnames(valid_terms))) {
      return(empty)
    }
    valid_key <- base::paste(base::as.character(valid_terms$module), base::as.character(valid_terms$term), sep = "\r")
    top_key <- base::paste(top_df$module, top_df$term, sep = "\r")
    top_df <- top_df[top_key %in% valid_key, , drop = FALSE]
    if (base::nrow(top_df) == 0) {
      return(empty)
    }
  }

  top_df$rank_num <- .hc_as_numeric_safely(top_df$rank)
  top_df$q_num <- .hc_as_numeric_safely(top_df$qvalue)
  out <- base::lapply(base::split(top_df, top_df$module), function(x) {
    x <- x[!base::duplicated(x$term), , drop = FALSE]
    x <- x[
      base::order(
        ifelse(base::is.na(x$rank_num), Inf, x$rank_num),
        ifelse(base::is.na(x$q_num), Inf, x$q_num),
        x$term
      ), ,
      drop = FALSE
    ]
    x$term_rank <- base::seq_len(base::nrow(x))
    x$term_rank_label <- base::paste0("Top ", x$term_rank)
    x$term_label <- .hc_safe_longitudinal_term_label(x$term, wrap_term = wrap_term)
    x
  })
  out <- out[!base::vapply(out, function(x) base::is.null(x) || base::nrow(x) == 0, FUN.VALUE = base::logical(1))]
  if (base::length(out) == 0) {
    return(empty)
  }
  out <- base::do.call(base::rbind, out)
  out <- out[, c("module", "term", "term_rank", "term_rank_label", "term_label", "module_color", "rank_num", "q_num"), drop = FALSE]
  base::colnames(out)[base::colnames(out) == "rank_num"] <- "rank"
  base::colnames(out)[base::colnames(out) == "q_num"] <- "qvalue"
  out <- out[!base::duplicated(base::paste(out$module, out$term, sep = "\r")), , drop = FALSE]
  base::rownames(out) <- NULL
  out
}

#' Try GSVA-based ssGSEA scoring across API variants
#' @noRd
.hc_try_gsva_ssgsea <- function(expr_mat,
                                gene_sets,
                                alpha = 0.25,
                                normalize = TRUE) {
  if (!requireNamespace("GSVA", quietly = TRUE)) {
    return(NULL)
  }
  err_state <- new.env(parent = emptyenv())
  err_state$message <- NULL
  out <- NULL

  if (base::exists("ssgseaParam", where = asNamespace("GSVA"), inherits = FALSE)) {
    out <- tryCatch(
      {
        param <- GSVA::ssgseaParam(
          exprData = expr_mat,
          geneSets = gene_sets,
          alpha = alpha,
          normalize = normalize
        )
        tryCatch(
          GSVA::gsva(param, verbose = FALSE),
          error = function(e) GSVA::gsva(param)
        )
      },
      error = function(e) {
        err_state$message <- conditionMessage(e)
        NULL
      }
    )
  }

  if (base::is.null(out)) {
    out <- tryCatch(
      {
        GSVA::gsva(
          expr = expr_mat,
          gset.idx.list = gene_sets,
          method = "ssgsea",
          alpha = alpha,
          ssgsea.norm = normalize,
          verbose = FALSE
        )
      },
      error = function(e1) {
        tryCatch(
          {
            GSVA::gsva(
              expr_mat,
              gene_sets,
              method = "ssgsea",
              alpha = alpha,
              ssgsea.norm = normalize,
              verbose = FALSE
            )
          },
          error = function(e2) {
            err_state$message <- base::paste(conditionMessage(e1), conditionMessage(e2), sep = " | ")
            NULL
          }
        )
      }
    )
  }

  if (base::is.null(out)) {
    return(NULL)
  }
  out <- base::as.matrix(out)
  if (base::nrow(out) == 0 || base::ncol(out) == 0) {
    return(NULL)
  }
  base::attr(out, "error_message") <- err_state$message
  out
}

#' Try GSVA-based z-score scoring across API variants
#' @noRd
.hc_try_gsva_zscore <- function(expr_mat, gene_sets) {
  if (!requireNamespace("GSVA", quietly = TRUE)) {
    return(NULL)
  }
  err_state <- new.env(parent = emptyenv())
  err_state$message <- NULL
  out <- NULL

  if (base::exists("zscoreParam", where = asNamespace("GSVA"), inherits = FALSE)) {
    out <- tryCatch(
      {
        param <- GSVA::zscoreParam(
          exprData = expr_mat,
          geneSets = gene_sets
        )
        tryCatch(
          GSVA::gsva(param, verbose = FALSE),
          error = function(e) GSVA::gsva(param)
        )
      },
      error = function(e) {
        err_state$message <- conditionMessage(e)
        NULL
      }
    )
  }

  if (base::is.null(out)) {
    out <- tryCatch(
      {
        GSVA::gsva(
          expr = expr_mat,
          gset.idx.list = gene_sets,
          method = "zscore",
          verbose = FALSE
        )
      },
      error = function(e1) {
        tryCatch(
          {
            GSVA::gsva(
              expr_mat,
              gene_sets,
              method = "zscore",
              verbose = FALSE
            )
          },
          error = function(e2) {
            err_state$message <- base::paste(conditionMessage(e1), conditionMessage(e2), sep = " | ")
            NULL
          }
        )
      }
    )
  }

  if (base::is.null(out)) {
    return(NULL)
  }
  out <- base::as.matrix(out)
  if (base::nrow(out) == 0 || base::ncol(out) == 0) {
    return(NULL)
  }
  base::attr(out, "error_message") <- err_state$message
  out
}

#' Z-scale scores within groups for plotting
#' @noRd
.hc_add_groupwise_zscores <- function(df,
                                      value_col = "score",
                                      group_cols = c("module", "term"),
                                      out_col = "score_plot") {
  if (!base::is.data.frame(df) || base::nrow(df) == 0) {
    df[[out_col]] <- numeric(0)
    return(df)
  }
  if (!base::all(c(value_col, group_cols) %in% base::colnames(df))) {
    stop("Missing columns needed for grouped z-scaling: ", base::paste(base::setdiff(c(value_col, group_cols), base::colnames(df)), collapse = ", "))
  }

  df[[out_col]] <- .hc_as_numeric_safely(df[[value_col]])
  split_idx <- base::split(
    base::seq_len(base::nrow(df)),
    base::do.call(base::paste, c(df[group_cols], sep = "\r"))
  )
  for (idx in split_idx) {
    vals <- .hc_as_numeric_safely(df[[value_col]][idx])
    keep <- base::is.finite(vals)
    if (!base::any(keep)) {
      df[[out_col]][idx] <- NA_real_
      next
    }
    if (base::sum(keep) <= 1 || !base::is.finite(stats::sd(vals[keep])) || stats::sd(vals[keep]) <= 0) {
      z <- base::rep(NA_real_, base::length(vals))
      z[keep] <- 0
      df[[out_col]][idx] <- z
      next
    }
    z <- base::rep(NA_real_, base::length(vals))
    z[keep] <- as.numeric(base::scale(vals[keep]))
    df[[out_col]][idx] <- z
  }
  df
}

#' Fast rank-mean enrichment score (fallback when GSVA is unavailable)
#' @noRd
.hc_rank_mean_scores <- function(expr_mat, gene_sets) {
  gene_sets <- gene_sets[base::lengths(gene_sets) > 0]
  out <- base::matrix(
    NA_real_,
    nrow = base::length(gene_sets),
    ncol = base::ncol(expr_mat),
    dimnames = list(base::names(gene_sets), base::colnames(expr_mat))
  )
  if (base::length(gene_sets) == 0) {
    return(out)
  }

  ranks <- base::apply(
    expr_mat,
    2,
    function(v) {
      base::rank(v, ties.method = "average", na.last = "keep")
    }
  )
  if (base::is.null(base::dim(ranks))) {
    ranks <- base::matrix(
      ranks,
      nrow = base::nrow(expr_mat),
      ncol = 1,
      dimnames = list(base::rownames(expr_mat), base::colnames(expr_mat))
    )
  }

  center <- (base::nrow(expr_mat) + 1) / 2
  denom <- base::max(1, base::nrow(expr_mat) / 2)

  for (i in base::seq_along(gene_sets)) {
    genes <- gene_sets[[i]]
    idx <- base::match(genes, base::rownames(expr_mat))
    idx <- base::unique(idx[!base::is.na(idx)])
    if (base::length(idx) == 0) {
      next
    }
    out[i, ] <- (base::colMeans(ranks[idx, , drop = FALSE], na.rm = TRUE) - center) / denom
  }
  out
}

#' Plot longitudinal enrichment waves for top terms per module
#'
#' Uses enrichment outputs from [hc_functional_enrichment()] to select top terms
#' per module for one or multiple databases (for example `Kegg`), then computes
#' sample-level term activity scores and visualizes trajectories across
#' donor/timepoint in module-faceted wave plots.
#'
#' Scoring uses GSVA ssGSEA when available (`score_method = "ssgsea"`), with a
#' rank-based fallback if GSVA is unavailable.
#'
#' @param hc A `HCoCenaExperiment`.
#' @param slot_name Satellite slot with longitudinal settings/module lookup.
#'   Usually `"longitudinal_module_means"`.
#' @param databases Character vector of databases to plot.
#' @param top Number of top terms per module.
#' @param custom_terms Optional custom term filter. Can be a character vector
#'   (applied to all modules) or a named list with module names (e.g. `M1`) as
#'   keys and character vectors of terms as values. If set, matched terms are
#'   plotted instead of automatic top-term selection.
#' @param term_match Matching mode for `custom_terms`: `"exact"`,
#'   `"contains"`, or `"regex"`.
#' @param enrichment_table Which enrichment table to use:
#'   `"selected"` (default), `"all"`, or `"significant"`.
#' @param score_method One of `"ssgsea"` (default), `"zscore"`, or `"rank_mean"`.
#' @param score_scale Plot-only scaling: `"none"` (default) or `"z"` to
#'   z-standardize scores within each module-term before plotting.
#' @param layer Optional layer selector (index/id/name). If `NULL`, uses
#'   `layer_id` from the longitudinal slot.
#' @param donor_col,time_col Optional donor/time annotation columns. If `NULL`,
#'   values from the longitudinal slot are used.
#' @param time_levels Optional explicit time ordering.
#' @param min_term_genes Minimum number of matched genes per term.
#' @param impute_missing One of `"auto"`, `"linear"`, `"locf"`, `"none"`.
#' @param qvalue_max Optional q-value filter before top-term selection.
#' @param show_donor_lines Logical; add faint donor lines behind term means.
#' @param facet_ncol Number of facet columns.
#' @param free_y Logical; free y-axis per module facet.
#' @param save_pdf Logical; save one PDF per database.
#' @param file_prefix Prefix used for output files.
#' @param save_width,save_height PDF size in inches.
#' @param export_excel Logical; export top-term + trajectory tables to Excel.
#'
#' @return A list with `plots`, `top_terms`, `mean_trajectories`,
#'   `donor_trajectories`, and `score_method_used`.
#' @examples
#' hc <- hc_example_data("clustered")
#' gmt <- system.file("extdata", "toy_celltype_markers.gmt", package = "hcocena")
#' hc <- hc_functional_enrichment(
#'   hc,
#'   gene_sets = character(0),
#'   custom_gmt_files = c(CellTypes = gmt),
#'   universe = "network"
#' )
#' hc <- hc_longitudinal_module_means(hc, donor_col = "donor",
#'                                    time_col = "timepoint")
#' # Database names are title-cased by hc_functional_enrichment().
#' waves <- hc_plot_longitudinal_enrichment_waves(
#'   hc,
#'   databases = "Celltypes",
#'   top = 1,
#'   min_term_genes = 3,
#'   score_method = "rank_mean",
#'   save_pdf = FALSE,
#'   export_excel = FALSE
#' )
#' waves$top_terms
#' @export
hc_plot_longitudinal_enrichment_waves <- function(hc,
                                                  slot_name = "longitudinal_module_means",
                                                  databases = "Kegg",
                                                  top = 5,
                                                  custom_terms = NULL,
                                                  term_match = c("exact", "contains", "regex"),
                                                  enrichment_table = c("selected", "all", "significant"),
                                                  score_method = c("ssgsea", "zscore", "rank_mean"),
                                                  score_scale = c("none", "z"),
                                                  layer = NULL,
                                                  donor_col = NULL,
                                                  time_col = NULL,
                                                  time_levels = NULL,
                                                  min_term_genes = 5,
                                                  impute_missing = c("auto", "linear", "locf", "none"),
                                                  qvalue_max = NULL,
                                                  show_donor_lines = TRUE,
                                                  facet_ncol = 4,
                                                  free_y = TRUE,
                                                  save_pdf = TRUE,
                                                  file_prefix = "Longitudinal_Enrichment_Waves",
                                                  save_width = 13,
                                                  save_height = 9,
                                                  export_excel = TRUE) {
  if (!inherits(hc, "HCoCenaExperiment")) {
    stop("`hc` must be a `HCoCenaExperiment`.")
  }
  enrichment_table <- base::match.arg(enrichment_table)
  score_method <- base::match.arg(score_method)
  score_scale <- base::match.arg(score_scale)
  impute_missing <- base::match.arg(impute_missing)
  term_match <- base::match.arg(term_match)

  if (!base::is.numeric(top) || base::length(top) != 1 || !base::is.finite(top) || top < 1) {
    stop("`top` must be a single numeric value >= 1.")
  }
  top <- base::as.integer(base::round(top))
  min_term_genes <- .hc_as_integer_safely(min_term_genes[[1]])
  if (!base::is.finite(min_term_genes) || min_term_genes < 1) {
    stop("`min_term_genes` must be >= 1.")
  }

  sat <- as.list(hc@satellite)
  obj <- sat[[slot_name]]
  if (base::is.null(obj) || !base::is.list(obj)) {
    stop("Slot `", slot_name, "` not found in `hc@satellite`.")
  }
  base_obj <- obj
  if ((base::is.null(base_obj$layer_id) || base::is.null(base_obj$donor_col) || base::is.null(base_obj$time_col)) &&
    !base::is.null(obj$source_slot) &&
    obj$source_slot %in% base::names(sat)) {
    src <- sat[[obj$source_slot]]
    if (base::is.list(src)) {
      base_obj <- src
    }
  }

  module_lookup <- if (!base::is.null(obj$module_lookup)) obj$module_lookup else base_obj$module_lookup
  if (base::is.null(module_lookup)) {
    stop("No `module_lookup` available in slot `", slot_name, "`.")
  }
  module_lookup <- base::as.data.frame(module_lookup, stringsAsFactors = FALSE)
  if (!base::all(c("module", "module_color") %in% base::colnames(module_lookup))) {
    stop("`module_lookup` must contain `module` and `module_color`.")
  }
  module_lookup$module <- base::as.character(module_lookup$module)
  module_lookup$module_color <- base::as.character(module_lookup$module_color)
  module_levels <- module_lookup$module
  module_col <- .hc_module_color_map(module_lookup = module_lookup, module_levels = module_levels)
  color_to_module <- stats::setNames(module_lookup$module, module_lookup$module_color)

  donor_col <- if (base::is.null(donor_col)) base_obj$donor_col else donor_col
  time_col <- if (base::is.null(time_col)) base_obj$time_col else time_col
  if (base::is.null(donor_col) || base::is.null(time_col)) {
    stop("`donor_col`/`time_col` missing. Provide them or run `hc_longitudinal_module_means()` first.")
  }

  layer_id <- if (base::is.null(layer)) base_obj$layer_id else .hc_resolve_layer_id(hc = hc, layer = layer)
  if (base::is.null(layer_id) || !base::nzchar(base::as.character(layer_id[[1]]))) {
    layer_id <- .hc_resolve_layer_id(hc = hc, layer = NULL)
  }
  layer_id <- base::as.character(layer_id[[1]])

  se <- MultiAssayExperiment::experiments(hc@mae)[[layer_id]]
  expr_mat <- SummarizedExperiment::assay(se, "counts")
  expr_mat <- base::as.matrix(expr_mat)
  if (!base::is.numeric(expr_mat)) {
    expr_mat <- base::matrix(
      data = .hc_as_numeric_safely(expr_mat),
      nrow = base::nrow(expr_mat),
      ncol = base::ncol(expr_mat),
      dimnames = base::dimnames(expr_mat)
    )
  }
  if (base::is.null(base::rownames(expr_mat)) || base::is.null(base::colnames(expr_mat))) {
    stop("Selected layer must have gene row names and sample column names.")
  }

  anno <- base::as.data.frame(SummarizedExperiment::colData(se), stringsAsFactors = FALSE)
  if (base::is.null(base::rownames(anno))) {
    base::rownames(anno) <- base::colnames(expr_mat)
  }
  sample_ids <- base::colnames(expr_mat)
  idx <- base::match(sample_ids, base::rownames(anno))
  if (base::any(base::is.na(idx))) {
    stop("Annotation rows missing for some expression samples.")
  }
  anno <- anno[idx, , drop = FALSE]
  if (!base::all(c(donor_col, time_col) %in% base::colnames(anno))) {
    stop("Missing `donor_col` or `time_col` in selected layer annotation.")
  }

  anno_df <- base::data.frame(
    sample_id = sample_ids,
    donor = base::as.character(anno[[donor_col]]),
    time = base::as.character(anno[[time_col]]),
    stringsAsFactors = FALSE
  )
  keep_samples <- !base::is.na(anno_df$donor) &
    !base::is.na(anno_df$time) &
    anno_df$donor != "" &
    anno_df$time != ""
  anno_df <- anno_df[keep_samples, , drop = FALSE]
  if (base::nrow(anno_df) == 0) {
    stop("No samples with valid donor/time labels.")
  }
  expr_mat <- expr_mat[, anno_df$sample_id, drop = FALSE]

  if (base::is.null(time_levels) || base::length(time_levels) == 0) {
    if (!base::is.null(obj$time_levels)) {
      time_levels <- base::as.character(obj$time_levels)
    } else if (!base::is.null(base_obj$time_levels)) {
      time_levels <- base::as.character(base_obj$time_levels)
    } else if (base::is.factor(anno[[time_col]])) {
      time_levels <- base::levels(anno[[time_col]])
    } else {
      time_levels <- base::sort(base::unique(anno_df$time))
    }
  } else {
    time_levels <- base::as.character(time_levels)
  }

  if (identical(impute_missing, "auto")) {
    impute_mode <- base_obj$impute_missing
    if (base::is.null(impute_mode) || !impute_mode %in% c("none", "linear", "locf")) {
      impute_mode <- "none"
    }
  } else {
    impute_mode <- impute_missing
  }

  enrich_df <- .hc_collect_enrichment_summary(hc = hc, enrichment_table = enrichment_table)
  if (base::is.null(enrich_df) || base::nrow(enrich_df) == 0) {
    stop("No enrichment summary found. Run `hc_functional_enrichment()` first.")
  }
  module_lookup <- .hc_harmonize_module_lookup_with_enrichment(
    module_lookup = module_lookup,
    enrich_df = enrich_df
  )
  module_levels <- module_lookup$module
  module_col <- .hc_module_color_map(module_lookup = module_lookup, module_levels = module_levels)
  color_to_module <- stats::setNames(module_lookup$module, module_lookup$module_color)

  databases <- base::unique(base::trimws(base::as.character(databases)))
  databases <- databases[!base::is.na(databases) & databases != ""]
  available_db <- base::unique(base::as.character(enrich_df$database))
  db_use <- available_db[base::match(base::tolower(databases), base::tolower(available_db))]
  db_use <- base::unique(db_use[!base::is.na(db_use)])
  if (base::length(db_use) == 0) {
    stop("None of the requested databases are available in enrichment outputs.")
  }

  if (!base::is.null(qvalue_max)) {
    qvalue_max <- .hc_first_numeric_value(qvalue_max[[1]])
    if (!base::is.finite(qvalue_max) || qvalue_max <= 0) {
      stop("`qvalue_max` must be NULL or a positive numeric value.")
    }
  }

  refs <- as.list(hc@references@data)
  ref_names <- base::names(refs)
  expr_genes <- base::rownames(expr_mat)
  expr_genes_upper <- base::toupper(expr_genes)

  wrap_term <- function(x) {
    x <- .hc_clean_enrichment_term_label(x)
    y <- base::strwrap(x, width = 34, simplify = FALSE)
    if (base::length(y) == 0) x else base::paste(y[[1]], collapse = "\n")
  }

  out_plots <- list()
  out_top <- list()
  out_donor <- list()
  out_mean <- list()
  out_method <- list()

  for (db_nm in db_use) {
    db_tbl <- enrich_df[base::tolower(enrich_df$database) == base::tolower(db_nm), , drop = FALSE]
    if (base::nrow(db_tbl) == 0) {
      next
    }

    module_from_cluster <- color_to_module[db_tbl$cluster]
    module_from_label <- base::trimws(base::as.character(db_tbl$module_label))
    use_label <- !base::is.na(module_from_label) &
      base::nzchar(module_from_label) &
      !base::vapply(module_from_label, .hc_is_color_token, FUN.VALUE = base::logical(1))
    module_id <- module_from_label
    use_cluster <- base::is.na(module_id) | !base::nzchar(module_id)
    module_id[use_cluster] <- base::as.character(module_from_cluster[use_cluster])
    db_tbl$module <- base::as.character(module_id)
    db_tbl <- db_tbl[!base::is.na(db_tbl$module) & db_tbl$module != "", , drop = FALSE]
    if (base::nrow(db_tbl) == 0) {
      next
    }
    db_tbl$module_color_resolved <- base::as.character(db_tbl$cluster)
    fill_mod_col <- base::is.na(db_tbl$module_color_resolved) | !base::nzchar(base::trimws(db_tbl$module_color_resolved))
    if (base::any(fill_mod_col)) {
      db_tbl$module_color_resolved[fill_mod_col] <- base::as.character(module_col[db_tbl$module[fill_mod_col]])
    }
    db_tbl$module_color_resolved[base::is.na(db_tbl$module_color_resolved) | !base::nzchar(base::trimws(db_tbl$module_color_resolved))] <- "grey60"

    if (!base::is.null(qvalue_max)) {
      db_tbl <- db_tbl[base::is.na(db_tbl$qvalue) | db_tbl$qvalue <= qvalue_max, , drop = FALSE]
      if (base::nrow(db_tbl) == 0) {
        next
      }
    }

    db_tbl$rank_num <- .hc_as_numeric_safely(db_tbl$rank)
    db_tbl$q_num <- .hc_as_numeric_safely(db_tbl$qvalue)
    use_custom_terms <- !base::is.null(custom_terms)
    if (use_custom_terms) {
      keep_custom <- .hc_match_custom_enrichment_terms(
        term = db_tbl$term,
        module = db_tbl$module,
        custom_terms = custom_terms,
        match_mode = term_match
      )
      db_tbl <- db_tbl[keep_custom, , drop = FALSE]
      if (base::nrow(db_tbl) == 0) {
        next
      }
    }

    top_rows <- base::lapply(base::split(db_tbl, db_tbl$module), function(x) {
      x <- x[!base::is.na(x$term) & x$term != "", , drop = FALSE]
      if (base::nrow(x) == 0) {
        return(NULL)
      }
      x <- x[!base::duplicated(x$term), , drop = FALSE]
      ord <- base::order(
        ifelse(base::is.na(x$rank_num), Inf, x$rank_num),
        ifelse(base::is.na(x$q_num), Inf, x$q_num),
        x$term
      )
      x <- x[ord, , drop = FALSE]
      if (use_custom_terms) {
        x
      } else {
        x[base::seq_len(base::min(top, base::nrow(x))), , drop = FALSE]
      }
    })
    top_rows <- top_rows[!base::vapply(top_rows, base::is.null, FUN.VALUE = base::logical(1))]
    if (base::length(top_rows) == 0) {
      next
    }
    top_df <- base::do.call(base::rbind, top_rows)
    top_df$module <- base::as.character(top_df$module)
    module_col_db <- db_tbl[!base::duplicated(db_tbl$module), c("module", "module_color_resolved"), drop = FALSE]
    module_col_db <- stats::setNames(base::as.character(module_col_db$module_color_resolved), base::as.character(module_col_db$module))
    top_df$module_color <- module_col_db[top_df$module]
    fill_top_col <- base::is.na(top_df$module_color) | !base::nzchar(base::trimws(top_df$module_color))
    if (base::any(fill_top_col)) {
      top_df$module_color[fill_top_col] <- base::as.character(module_col[top_df$module[fill_top_col]])
    }
    top_df$module_color[base::is.na(top_df$module_color) | !base::nzchar(base::trimws(top_df$module_color))] <- "grey60"
    base::rownames(top_df) <- NULL

    term_ids <- base::unique(base::as.character(top_df$term))
    ref_idx <- base::match(base::tolower(db_nm), base::tolower(ref_names))
    ref_db <- if (base::is.na(ref_idx)) NA_character_ else ref_names[[ref_idx]]
    term_genes <- list()
    if (!base::is.na(ref_db)) {
      ref_df <- refs[[ref_db]]
      if (base::is.data.frame(ref_df) && base::all(c("term", "gene") %in% base::colnames(ref_df))) {
        ref_df <- base::as.data.frame(ref_df, stringsAsFactors = FALSE)
        ref_df <- ref_df[ref_df$term %in% term_ids, c("term", "gene"), drop = FALSE]
        if (base::nrow(ref_df) > 0) {
          term_genes <- base::split(base::as.character(ref_df$gene), base::as.character(ref_df$term))
        }
      }
    }
    missing_terms <- base::setdiff(term_ids, base::names(term_genes))
    if (base::length(missing_terms) > 0) {
      for (tm in missing_terms) {
        ids <- top_df$geneID[top_df$term == tm]
        ids <- ids[!base::is.na(ids) & ids != ""]
        if (base::length(ids) == 0) next
        g <- base::unlist(base::strsplit(base::paste(ids, collapse = "/"), "/", fixed = TRUE), use.names = FALSE)
        g <- base::trimws(base::as.character(g))
        g <- g[g != ""]
        if (base::length(g) > 0) term_genes[[tm]] <- base::unique(g)
      }
    }

    term_gene_sets <- list()
    term_n <- stats::setNames(base::integer(0), base::character(0))
    for (tm in base::names(term_genes)) {
      gu <- base::toupper(base::as.character(term_genes[[tm]]))
      idx_g <- base::match(gu, expr_genes_upper)
      idx_g <- base::unique(idx_g[!base::is.na(idx_g)])
      if (base::length(idx_g) < min_term_genes) next
      term_gene_sets[[tm]] <- expr_genes[idx_g]
      term_n[[tm]] <- base::length(idx_g)
    }
    if (base::length(term_gene_sets) == 0) {
      next
    }
    top_df <- top_df[top_df$term %in% base::names(term_gene_sets), , drop = FALSE]
    if (base::nrow(top_df) == 0) {
      next
    }

    method_now <- score_method
    if (identical(score_method, "ssgsea")) {
      score_mat <- .hc_try_gsva_ssgsea(expr_mat = expr_mat, gene_sets = term_gene_sets, alpha = 0.25, normalize = TRUE)
      if (base::is.null(score_mat)) {
        score_mat <- .hc_rank_mean_scores(expr_mat = expr_mat, gene_sets = term_gene_sets)
        method_now <- "rank_mean"
      }
    } else if (identical(score_method, "zscore")) {
      score_mat <- .hc_try_gsva_zscore(expr_mat = expr_mat, gene_sets = term_gene_sets)
      if (base::is.null(score_mat)) {
        score_mat <- .hc_rank_mean_scores(expr_mat = expr_mat, gene_sets = term_gene_sets)
        method_now <- "rank_mean"
      }
    } else {
      score_mat <- .hc_rank_mean_scores(expr_mat = expr_mat, gene_sets = term_gene_sets)
      method_now <- "rank_mean"
    }

    score_mat <- score_mat[base::rownames(score_mat) %in% base::unique(top_df$term), , drop = FALSE]
    if (base::nrow(score_mat) == 0) {
      next
    }

    score_long <- base::as.data.frame(base::as.table(score_mat), stringsAsFactors = FALSE)
    base::colnames(score_long) <- c("term", "sample_id", "score")
    score_long <- base::merge(score_long, anno_df, by = "sample_id", all.x = TRUE, sort = FALSE)

    top_annot <- top_df[, c("database", "module", "module_color", "term", "rank_num", "q_num"), drop = FALSE]
    base::colnames(top_annot) <- c("database", "module", "module_color", "term", "rank", "qvalue")
    plot_df <- base::merge(top_annot, score_long, by = "term", all.x = TRUE, sort = FALSE)
    if (base::nrow(plot_df) == 0) {
      next
    }

    agg <- stats::aggregate(
      score ~ donor + time + database + module + module_color + term + rank + qvalue,
      data = plot_df,
      FUN = function(x) base::mean(base::as.numeric(x), na.rm = TRUE)
    )
    agg$score[!base::is.finite(agg$score)] <- NA_real_

    donors <- base::sort(base::unique(anno_df$donor))
    comb <- base::unique(agg[, c("database", "module", "module_color", "term", "rank", "qvalue"), drop = FALSE])
    grid_rows <- base::lapply(base::seq_len(base::nrow(comb)), function(i) {
      g <- base::expand.grid(donor = donors, time = time_levels, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
      cbind(g, comb[base::rep(i, base::nrow(g)), , drop = FALSE], stringsAsFactors = FALSE)
    })
    full_df <- base::do.call(base::rbind, grid_rows)
    full_df <- base::merge(
      full_df,
      agg,
      by = c("donor", "time", "database", "module", "module_color", "term", "rank", "qvalue"),
      all.x = TRUE,
      sort = FALSE
    )
    full_df$time <- base::factor(base::as.character(full_df$time), levels = time_levels, ordered = TRUE)
    full_df$module <- base::factor(base::as.character(full_df$module), levels = module_levels)
    full_df <- full_df[!base::is.na(full_df$module), , drop = FALSE]
    full_df <- full_df[base::order(full_df$module, full_df$term, full_df$donor, full_df$time), , drop = FALSE]
    if (!identical(impute_mode, "none")) {
      split_idx <- base::split(
        base::seq_len(base::nrow(full_df)),
        base::paste(full_df$donor, full_df$module, full_df$term, sep = "||")
      )
      for (idv in split_idx) {
        full_df$score[idv] <- .hc_impute_series(full_df$score[idv], method = impute_mode)
      }
    }
    full_df$score_plot <- full_df$score
    if (identical(score_scale, "z")) {
      full_df <- .hc_add_groupwise_zscores(
        df = full_df,
        value_col = "score",
        group_cols = c("module", "term"),
        out_col = "score_plot"
      )
    }

    mean_df <- stats::aggregate(
      score ~ database + module + module_color + term + rank + qvalue + time,
      data = full_df,
      FUN = function(x) base::mean(base::as.numeric(x), na.rm = TRUE)
    )
    mean_df$score[!base::is.finite(mean_df$score)] <- NA_real_
    mean_plot_df <- stats::aggregate(
      score_plot ~ database + module + module_color + term + rank + qvalue + time,
      data = full_df,
      FUN = function(x) base::mean(base::as.numeric(x), na.rm = TRUE)
    )
    mean_plot_df$score_plot[!base::is.finite(mean_plot_df$score_plot)] <- NA_real_

    ord_df <- top_annot
    ord_df$module <- base::factor(base::as.character(ord_df$module), levels = module_levels)
    ord_df <- ord_df[base::order(ord_df$module, ord_df$rank, ord_df$term), , drop = FALSE]
    term_levels <- base::unique(base::as.character(ord_df$term))
    term_labels <- base::make.unique(base::vapply(term_levels, wrap_term, FUN.VALUE = base::character(1)), sep = " ")
    term_map <- stats::setNames(term_labels, term_levels)
    full_df$term_display <- base::factor(term_map[base::as.character(full_df$term)], levels = term_labels)
    mean_df$term_display <- base::factor(term_map[base::as.character(mean_df$term)], levels = term_labels)
    mean_plot_df$term_display <- base::factor(term_map[base::as.character(mean_plot_df$term)], levels = term_labels)
    term_pal <- .hc_distinct_palette(base::length(term_labels))
    base::names(term_pal) <- term_labels

    module_levels_db <- .hc_order_module_levels(
      base::c(
        module_levels[module_levels %in% base::unique(base::as.character(full_df$module))],
        base::unique(base::as.character(full_df$module))
      )
    )
    module_cols_db <- .hc_module_color_map(
      module_lookup = rbind(
        module_lookup[, c("module", "module_color"), drop = FALSE],
        top_df[, c("module", "module_color"), drop = FALSE]
      ),
      module_levels = module_levels_db
    )

    p <- ggplot2::ggplot()
    if (isTRUE(show_donor_lines)) {
      p <- p + ggplot2::geom_line(
        data = full_df,
        mapping = ggplot2::aes(
          x = time,
          y = score_plot,
          group = interaction(donor, term_display),
          color = term_display
        ),
        alpha = 0.06,
        linewidth = 0.22,
        show.legend = FALSE,
        na.rm = TRUE
      )
    }
    p <- p +
      ggplot2::geom_line(
        data = mean_plot_df,
        mapping = ggplot2::aes(
          x = time,
          y = score_plot,
          group = term_display,
          color = term_display
        ),
        linewidth = 1.15,
        na.rm = TRUE
      ) +
      ggplot2::geom_point(
        data = mean_plot_df,
        mapping = ggplot2::aes(
          x = time,
          y = score_plot,
          group = term_display,
          color = term_display
        ),
        size = 1.7,
        na.rm = TRUE
      ) +
      .hc_module_facet_wrap(
        module_levels = module_levels_db,
        mod_col = module_cols_db,
        ncol = facet_ncol,
        scales = if (isTRUE(free_y)) "free_y" else "fixed",
        facets = ~module
      ) +
      ggplot2::scale_color_manual(values = term_pal, drop = FALSE) +
      .hc_theme_pub(base_size = 11) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        legend.position = "right"
      ) +
      ggplot2::guides(
        color = ggplot2::guide_legend(
          title = "Enriched term",
          override.aes = list(alpha = 1, linewidth = 1.1, size = 2.3)
        )
      ) +
      ggplot2::labs(
        title = if (base::is.null(custom_terms)) {
          base::paste0(
            db_nm,
            " top ",
            top,
            " terms per module (",
            if (identical(method_now, "ssgsea")) "ssGSEA" else "rank score",
            ")"
          )
        } else {
          base::paste0(
            db_nm,
            " selected terms per module (",
            if (identical(method_now, "ssgsea")) "ssGSEA" else "rank score",
            ")"
          )
        },
        x = "Timepoint",
        y = if (identical(score_scale, "z")) "Enrichment score (z-scaled)" else if (identical(method_now, "ssgsea")) "ssGSEA score" else if (identical(method_now, "zscore")) "GSVA z-score" else "Rank enrichment score"
      )
    p <- .hc_colorize_module_strips(
      p = p,
      module_levels = module_levels_db,
      mod_col = module_cols_db
    )

    db_token <- gsub("[^A-Za-z0-9]+", "_", db_nm)
    db_token <- gsub("^_+|_+$", "", db_token)
    if (!base::nzchar(db_token)) db_token <- "DB"

    if (isTRUE(save_pdf)) {
      out_dir <- .hc_resolve_output_dir(hc)
      .hc_ggsave_pdf_png(
        filename = base::file.path(out_dir, base::paste0(file_prefix, "_", db_token, ".pdf")),
        plot = p,
        width = save_width,
        height = save_height,
        units = "in",
        device = grDevices::cairo_pdf
      )
    }

    top_export <- top_annot
    top_export$n_genes_used <- as.integer(term_n[base::as.character(top_export$term)])
    top_export$score_method <- method_now
    top_export$score_scale <- score_scale
    top_export$enrichment_table <- enrichment_table
    top_export <- top_export[
      base::order(base::match(base::as.character(top_export$module), module_levels), top_export$rank, top_export$term), ,
      drop = FALSE
    ]
    base::rownames(top_export) <- NULL

    if (isTRUE(export_excel) && requireNamespace("openxlsx", quietly = TRUE)) {
      out_dir <- .hc_resolve_output_dir(hc)
      .hc_write_xlsx_atomic(
        x = list(
          top_terms = top_export,
          mean_trajectories = mean_df,
          mean_trajectories_plot = mean_plot_df,
          donor_trajectories = full_df
        ),
        file = base::file.path(out_dir, base::paste0(file_prefix, "_", db_token, ".xlsx")),
        overwrite = TRUE
      )
    }

    out_plots[[db_nm]] <- p
    out_top[[db_nm]] <- top_export
    out_donor[[db_nm]] <- full_df
    mean_df <- base::merge(
      mean_df,
      mean_plot_df,
      by = c("database", "module", "module_color", "term", "rank", "qvalue", "time", "term_display"),
      all.x = TRUE,
      sort = FALSE
    )
    mean_df$mean_enrichment_score <- mean_df$score
    mean_df$mean_enrichment_score_plot <- mean_df$score_plot
    mean_df$score_scale <- score_scale
    out_mean[[db_nm]] <- mean_df
    out_method[[db_nm]] <- method_now
  }

  if (base::length(out_plots) == 0) {
    stop("No database produced plottable longitudinal enrichment waves.")
  }

  list(
    plots = out_plots,
    top_terms = out_top,
    donor_trajectories = out_donor,
    mean_trajectories = out_mean,
    score_method_used = out_method
  )
}

#' Plot longitudinal enrichment waves by meta-cluster
#'
#' Creates one plot per enrichment database. For each module, the top `x` terms
#' are placed into dedicated side-by-side panels (`Top 1`, `Top 2`, ...), and
#' meta-cluster trajectories are drawn as longitudinal waves.
#'
#' @param hc A `HCoCenaExperiment`.
#' @param slot_name Satellite slot with `meta_cluster` results. Default:
#'   `"longitudinal_endotypes"`. Use a concrete slot name, a shared prefix
#'   such as `"longitudinal_endotypes"` to process matching suffixed slots, or
#'   `"all"` to run over every compatible step 2 slot in `hc@satellite`.
#' @param databases Character vector of databases to visualize (e.g. `"Kegg"`).
#' @param top Integer number of top terms per module.
#' @param custom_terms Optional custom term filter (character vector or named
#'   module list). If set, matched terms are used instead of automatic top-term
#'   selection.
#' @param term_match Matching mode for `custom_terms`: `"exact"`,
#'   `"contains"`, or `"regex"`.
#' @param enrichment_table Source enrichment table: `"all"` (default),
#'   `"selected"`, or `"significant"`.
#' @param score_method `"ssgsea"` (default), `"zscore"`, or `"rank_mean"`.
#' @param score_scale Plot-only scaling: `"none"` (default) or `"z"` to
#'   z-standardize scores within each module-term before plotting.
#' @param layer Optional layer selector forwarded to
#'   [hc_plot_longitudinal_enrichment_waves()].
#' @param donor_col,time_col,time_levels Optional overrides forwarded to
#'   [hc_plot_longitudinal_enrichment_waves()].
#' @param min_term_genes Minimum matched genes per term.
#' @param impute_missing Missing-value handling.
#' @param qvalue_max Optional q-value cutoff before top-term selection.
#' @param show_donor_lines Logical; include faint donor trajectories.
#' @param donor_alpha,donor_linewidth Styling for donor lines.
#' @param mean_linewidth,mean_point_size Styling for meta-cluster means.
#' @param show_ci Logical; if `TRUE`, draw subtle confidence ribbons around
#'   meta-cluster mean trajectories.
#' @param ci_level Confidence level in `(0,1)` for ribbons (default `0.95`).
#' @param ci_alpha Ribbon alpha for confidence bands.
#' @param free_y Logical; use free y scales per module row.
#' @param save_pdf Logical; save one PDF per selected database.
#' @param file_prefix Prefix for saved files.
#' @param save_width,save_height PDF size.
#' @param export_excel Logical; export tables per database.
#'
#' @return A list with `plots`, `top_terms`, `donor_trajectories`,
#'   `mean_trajectories`, and `score_method_used`. When multiple slots are
#'   processed, these top-level components are flattened with slot-prefixed
#'   names for convenient iteration, and nested per-slot outputs are also
#'   returned in `results_by_slot`, together with `slot_info`.
#' @examples
#' hc <- hc_example_data("clustered")
#' gmt <- system.file("extdata", "toy_celltype_markers.gmt", package = "hcocena")
#' hc <- hc_functional_enrichment(
#'   hc,
#'   gene_sets = character(0),
#'   custom_gmt_files = c(CellTypes = gmt),
#'   universe = "network"
#' )
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
#' waves <- hc_plot_longitudinal_enrichment_meta_waves(
#'   hc,
#'   databases = "Celltypes",
#'   top = 1,
#'   min_term_genes = 3,
#'   score_method = "rank_mean",
#'   save_pdf = FALSE,
#'   export_excel = FALSE
#' )
#' @export
hc_plot_longitudinal_enrichment_meta_waves <- function(hc,
                                                       slot_name = "longitudinal_endotypes",
                                                       databases = "Kegg",
                                                       top = 3,
                                                       custom_terms = NULL,
                                                       term_match = c("exact", "contains", "regex"),
                                                       enrichment_table = c("all", "selected", "significant"),
                                                       score_method = c("ssgsea", "zscore", "rank_mean"),
                                                       score_scale = c("none", "z"),
                                                       layer = NULL,
                                                       donor_col = NULL,
                                                       time_col = NULL,
                                                       time_levels = NULL,
                                                       min_term_genes = 5,
                                                       impute_missing = c("auto", "linear", "locf", "none"),
                                                       qvalue_max = NULL,
                                                       show_donor_lines = TRUE,
                                                       donor_alpha = 0.04,
                                                       donor_linewidth = 0.18,
                                                       mean_linewidth = 2.0,
                                                       mean_point_size = 1.8,
                                                       show_ci = FALSE,
                                                       ci_level = 0.95,
                                                       ci_alpha = 0.12,
                                                       free_y = TRUE,
                                                       save_pdf = TRUE,
                                                       file_prefix = "Longitudinal_Enrichment_MetaWaves",
                                                       save_width = 12,
                                                       save_height = 12,
                                                       export_excel = TRUE) {
  hc <- .hc_longitudinal_unwrap_hc(hc)
  enrichment_table <- base::match.arg(enrichment_table)
  score_method <- base::match.arg(score_method)
  score_scale <- base::match.arg(score_scale)
  impute_missing <- base::match.arg(impute_missing)
  term_match <- base::match.arg(term_match)
  free_y <- TRUE
  ci_level <- .hc_first_numeric_value(ci_level[[1]])
  ci_alpha <- .hc_first_numeric_value(ci_alpha[[1]])
  if (!base::is.finite(ci_level) || ci_level <= 0 || ci_level >= 1) {
    stop("`ci_level` must be a numeric value in (0, 1).")
  }
  if (!base::is.finite(ci_alpha) || ci_alpha < 0 || ci_alpha > 1) {
    stop("`ci_alpha` must be a numeric value between 0 and 1.")
  }

  target_slots <- .hc_longitudinal_step3_target_slots(hc = hc, slot_name = slot_name)
  requested_slot_family <- .hc_longitudinal_requested_slot_family(slot_name)
  if (base::length(target_slots) > 1) {
    slot_results <- list()
    slot_info <- base::vector("list", base::length(target_slots))

    flatten_component <- function(results, field) {
      out <- list()
      for (nm in base::names(results)) {
        comp <- results[[nm]][[field]]
        if (base::is.null(comp)) {
          next
        }
        comp_names <- base::names(comp)
        if (base::is.null(comp_names) || base::any(!base::nzchar(comp_names))) {
          comp_names <- base::paste0(field, "_", base::seq_along(comp))
        }
        base::names(comp) <- base::paste0(nm, "__", comp_names)
        out <- base::c(out, comp)
      }
      out
    }

    for (i in base::seq_along(target_slots)) {
      slot_i <- target_slots[[i]]
      slot_ctx <- .hc_longitudinal_step2_slot_context(
        hc = hc,
        slot_name = slot_i,
        requested_slot_name = slot_name
      )
      output_name <- slot_ctx$output_name
      if (!base::nzchar(output_name) || output_name %in% base::names(slot_results)) {
        output_name <- slot_i
      }

      res_i <- hc_plot_longitudinal_enrichment_meta_waves(
        hc = hc,
        slot_name = slot_i,
        databases = databases,
        top = top,
        custom_terms = custom_terms,
        term_match = term_match,
        enrichment_table = enrichment_table,
        score_method = score_method,
        score_scale = score_scale,
        layer = layer,
        donor_col = donor_col,
        time_col = time_col,
        time_levels = time_levels,
        min_term_genes = min_term_genes,
        impute_missing = impute_missing,
        qvalue_max = qvalue_max,
        show_donor_lines = show_donor_lines,
        donor_alpha = donor_alpha,
        donor_linewidth = donor_linewidth,
        mean_linewidth = mean_linewidth,
        mean_point_size = mean_point_size,
        show_ci = show_ci,
        ci_level = ci_level,
        ci_alpha = ci_alpha,
        free_y = free_y,
        save_pdf = save_pdf,
        file_prefix = base::paste0(file_prefix, "_", slot_ctx$file_suffix),
        save_width = save_width,
        save_height = save_height,
        export_excel = export_excel
      )

      slot_results[[output_name]] <- res_i
      slot_info[[i]] <- base::data.frame(
        output_name = output_name,
        slot_name = slot_i,
        layer_id = slot_ctx$layer_id,
        layer_name = slot_ctx$layer_name,
        source_slot = slot_ctx$source_slot,
        stringsAsFactors = FALSE
      )
    }

    return(list(
      plots = flatten_component(slot_results, "plots"),
      top_terms = flatten_component(slot_results, "top_terms"),
      donor_trajectories = flatten_component(slot_results, "donor_trajectories"),
      mean_trajectories = flatten_component(slot_results, "mean_trajectories"),
      score_method_used = flatten_component(slot_results, "score_method_used"),
      results_by_slot = slot_results,
      slot_info = base::do.call(base::rbind, slot_info)
    ))
  }

  if (base::length(target_slots) == 1 && isTRUE(requested_slot_family)) {
    slot_ctx <- .hc_longitudinal_step2_slot_context(
      hc = hc,
      slot_name = target_slots[[1]],
      requested_slot_name = slot_name
    )
    slot_name <- target_slots[[1]]
    if (!base::is.null(slot_ctx$file_suffix) &&
      base::nzchar(slot_ctx$file_suffix) &&
      !base::grepl(base::paste0("_", slot_ctx$file_suffix, "$"), file_prefix)) {
      file_prefix <- base::paste0(file_prefix, "_", slot_ctx$file_suffix)
    }
  }

  sq_size <- .hc_max_finite(base::c(save_width[[1]], save_height[[1]]))
  if (!base::is.finite(sq_size) || sq_size <= 0) {
    sq_size <- 12
  }
  save_width <- sq_size
  save_height <- sq_size

  sat <- as.list(hc@satellite)
  obj <- sat[[slot_name]]
  if (base::is.null(obj) || !base::is.list(obj) || base::is.null(obj$meta_cluster)) {
    stop("Slot `", slot_name, "` must contain `meta_cluster`. Run `hc_longitudinal_step2_meta_clustering()` first.")
  }
  meta_tbl <- base::as.data.frame(obj$meta_cluster, stringsAsFactors = FALSE)
  if (!base::all(c("donor", "meta_cluster") %in% base::colnames(meta_tbl))) {
    stop("`meta_cluster` table must contain columns `donor` and `meta_cluster`.")
  }
  meta_tbl$donor <- base::as.character(meta_tbl$donor)
  meta_tbl$meta_cluster <- base::as.character(meta_tbl$meta_cluster)
  meta_tbl <- meta_tbl[!base::is.na(meta_tbl$donor) & !base::is.na(meta_tbl$meta_cluster), , drop = FALSE]
  if (base::nrow(meta_tbl) == 0) {
    stop("No valid donor/meta-cluster assignments found.")
  }

  # Reuse base enrichment scorer (term selection + ssGSEA/rank scoring).
  base_out <- hc_plot_longitudinal_enrichment_waves(
    hc = hc,
    slot_name = slot_name,
    databases = databases,
    top = top,
    custom_terms = custom_terms,
    term_match = term_match,
    enrichment_table = enrichment_table,
    score_method = score_method,
    score_scale = "none",
    layer = layer,
    donor_col = donor_col,
    time_col = time_col,
    time_levels = time_levels,
    min_term_genes = min_term_genes,
    impute_missing = impute_missing,
    qvalue_max = qvalue_max,
    show_donor_lines = FALSE,
    facet_ncol = 4,
    free_y = free_y,
    save_pdf = FALSE,
    file_prefix = file_prefix,
    save_width = save_width,
    save_height = save_height,
    export_excel = FALSE
  )

  out_plots <- list()
  out_top <- list()
  out_donor <- list()
  out_mean <- list()
  out_method <- list()

  wrap_term <- function(x) {
    x <- .hc_clean_enrichment_term_label(x)
    y <- base::strwrap(base::as.character(x), width = 26, simplify = FALSE)
    if (base::length(y) == 0) {
      base::as.character(x)
    } else {
      base::paste(y[[1]], collapse = "\n")
    }
  }

  for (db_nm in base::names(base_out$donor_trajectories)) {
    donor_df <- base::as.data.frame(base_out$donor_trajectories[[db_nm]], stringsAsFactors = FALSE)
    top_df <- base::as.data.frame(base_out$top_terms[[db_nm]], stringsAsFactors = FALSE)
    if (base::nrow(donor_df) == 0 || base::nrow(top_df) == 0) {
      next
    }

    donor_df$donor <- base::as.character(donor_df$donor)
    donor_df$module <- base::as.character(donor_df$module)
    donor_df$term <- base::as.character(donor_df$term)
    donor_df <- base::merge(
      donor_df,
      meta_tbl[, c("donor", "meta_cluster"), drop = FALSE],
      by = "donor",
      all.x = TRUE,
      sort = FALSE
    )
    donor_df <- donor_df[!base::is.na(donor_df$meta_cluster), , drop = FALSE]
    if (base::nrow(donor_df) == 0) {
      next
    }

    pair_key <- function(module, term) {
      base::paste(base::as.character(module), base::as.character(term), sep = "\r")
    }

    valid_terms <- if ("score" %in% base::colnames(donor_df)) {
      keep_terms <- stats::aggregate(
        score ~ module + term,
        data = donor_df,
        FUN = function(x) {
          base::any(base::is.finite(base::as.numeric(x)))
        }
      )
      keep_terms <- keep_terms[base::isTRUE(keep_terms$score) | keep_terms$score %in% TRUE, c("module", "term"), drop = FALSE]
      keep_terms
    } else {
      base::unique(donor_df[, c("module", "term"), drop = FALSE])
    }

    rank_lookup <- .hc_longitudinal_rank_lookup(
      top_df = top_df,
      wrap_term = wrap_term,
      valid_terms = valid_terms
    )
    if (base::is.null(custom_terms)) {
      rank_lookup <- rank_lookup[rank_lookup$term_rank <= as.integer(top), , drop = FALSE]
    }
    if (base::nrow(rank_lookup) == 0) {
      next
    }

    drop_cols <- base::intersect(
      c("term_rank", "term_rank_label", "term_label", "term_display", "panel_label"),
      base::colnames(donor_df)
    )
    if (base::length(drop_cols) > 0) {
      donor_df <- donor_df[, base::setdiff(base::colnames(donor_df), drop_cols), drop = FALSE]
    }
    donor_df <- base::merge(
      donor_df,
      rank_lookup[, c("module", "term", "term_rank", "term_rank_label", "term_label"), drop = FALSE],
      by = c("module", "term"),
      all.x = TRUE,
      sort = FALSE
    )
    donor_df <- donor_df[!base::is.na(donor_df$term_rank), , drop = FALSE]
    if (base::nrow(donor_df) == 0) {
      next
    }

    # Keep one canonical module color column even if upstream tables changed.
    module_col_src <- rank_lookup[!base::duplicated(base::as.character(rank_lookup$module)), c("module", "module_color"), drop = FALSE]
    module_col_lookup <- stats::setNames(
      base::as.character(module_col_src$module_color),
      base::as.character(module_col_src$module)
    )
    if ("module_color" %in% base::colnames(donor_df)) {
      donor_df$module_color <- base::as.character(donor_df$module_color)
    } else {
      donor_df$module_color <- module_col_lookup[base::as.character(donor_df$module)]
    }
    fill_idx <- base::is.na(donor_df$module_color) | donor_df$module_color == ""
    if (base::any(fill_idx)) {
      donor_df$module_color[fill_idx] <- module_col_lookup[base::as.character(donor_df$module[fill_idx])]
    }
    donor_df$module_color[base::is.na(donor_df$module_color) | donor_df$module_color == ""] <- "grey60"

    # Ordering
    module_levels <- .hc_order_module_levels(rank_lookup$module)
    module_col <- stats::setNames(base::as.character(rank_lookup$module_color), base::as.character(rank_lookup$module))
    module_col <- module_col[module_levels]
    module_col[base::is.na(module_col) | !base::nzchar(module_col)] <- "grey60"

    ord_df <- rank_lookup
    ord_df$module <- base::factor(base::as.character(ord_df$module), levels = module_levels)
    ord_df <- ord_df[base::order(ord_df$module, ord_df$term_rank, ord_df$term), , drop = FALSE]
    ord_keys <- pair_key(ord_df$module, ord_df$term)
    term_labels <- base::make.unique(base::as.character(ord_df$term_label), sep = " ")
    term_map <- stats::setNames(term_labels, ord_keys)
    rank_label_map <- stats::setNames(base::as.character(ord_df$term_rank_label), ord_keys)
    safe_term_display <- function(module, term) {
      keys <- pair_key(module, term)
      out <- term_map[keys]
      fill_idx <- base::is.na(out) | !base::nzchar(base::trimws(out))
      if (base::any(fill_idx)) {
        out[fill_idx] <- .hc_safe_longitudinal_term_label(term[fill_idx], wrap_term = wrap_term)
      }
      out
    }
    is_valid_panel_field <- function(x) {
      x <- base::as.character(x)
      !base::is.na(x) & base::nzchar(base::trimws(x)) & x != "NA"
    }

    donor_keys <- pair_key(donor_df$module, donor_df$term)
    donor_df$term_display <- safe_term_display(donor_df$module, donor_df$term)
    donor_df$term_display <- base::factor(base::as.character(donor_df$term_display), levels = base::unique(term_labels))
    donor_df$term_rank_label <- rank_label_map[donor_keys]
    donor_df$panel_label <- .hc_longitudinal_panel_label(
      module = donor_df$module,
      rank_label = donor_df$term_rank_label,
      term_display = donor_df$term_display
    )

    if (base::is.factor(donor_df$time)) {
      t_levels <- base::levels(donor_df$time)
    } else if (!base::is.null(time_levels) && base::length(time_levels) > 0) {
      t_levels <- base::as.character(time_levels)
    } else {
      t_levels <- base::sort(base::unique(base::as.character(donor_df$time)))
    }
    donor_df$time <- base::factor(base::as.character(donor_df$time), levels = t_levels, ordered = TRUE)
    donor_df$module <- base::factor(base::as.character(donor_df$module), levels = module_levels)
    donor_df$score_plot <- donor_df$score
    if (identical(score_scale, "z")) {
      donor_df <- .hc_add_groupwise_zscores(
        df = donor_df,
        value_col = "score",
        group_cols = c("module", "term"),
        out_col = "score_plot"
      )
    }

    meta_levels <- base::sort(base::unique(base::as.character(donor_df$meta_cluster)))
    donor_df$meta_cluster <- base::factor(base::as.character(donor_df$meta_cluster), levels = meta_levels)
    meta_pal <- .hc_distinct_palette(base::length(meta_levels))
    base::names(meta_pal) <- meta_levels
    mean_agg <- stats::aggregate(
      score ~ module + module_color + term + term_rank + term_rank_label + term_label + time + meta_cluster,
      data = donor_df,
      FUN = function(x) {
        xv <- base::as.numeric(x)
        xv <- xv[base::is.finite(xv)]
        n <- base::length(xv)
        if (n == 0) {
          return(c(mean = NA_real_, lower = NA_real_, upper = NA_real_, n = 0))
        }
        m <- base::mean(xv)
        if (n == 1) {
          return(c(mean = m, lower = NA_real_, upper = NA_real_, n = 1))
        }
        se <- stats::sd(xv) / base::sqrt(n)
        err <- stats::qt((1 + ci_level) / 2, df = n - 1) * se
        c(mean = m, lower = m - err, upper = m + err, n = n)
      }
    )
    mean_df <- mean_agg
    score_stats <- mean_agg$score
    mean_df$score <- .hc_as_numeric_safely(score_stats[, "mean"])
    mean_df$ci_lower <- .hc_as_numeric_safely(score_stats[, "lower"])
    mean_df$ci_upper <- .hc_as_numeric_safely(score_stats[, "upper"])
    mean_df$n_donor <- .hc_as_integer_safely(score_stats[, "n"])
    mean_df$score[!base::is.finite(mean_df$score)] <- NA_real_
    mean_df$ci_lower[!base::is.finite(mean_df$ci_lower)] <- NA_real_
    mean_df$ci_upper[!base::is.finite(mean_df$ci_upper)] <- NA_real_
    mean_df$mean_enrichment_score <- mean_df$score
    mean_df$lower_enrichment_score <- mean_df$ci_lower
    mean_df$upper_enrichment_score <- mean_df$ci_upper
    mean_df$database <- db_nm
    mean_plot_agg <- stats::aggregate(
      score_plot ~ module + module_color + term + term_rank + term_rank_label + term_label + time + meta_cluster,
      data = donor_df,
      FUN = function(x) {
        xv <- base::as.numeric(x)
        xv <- xv[base::is.finite(xv)]
        n <- base::length(xv)
        if (n == 0) {
          return(c(mean = NA_real_, lower = NA_real_, upper = NA_real_, n = 0))
        }
        m <- base::mean(xv)
        if (n == 1) {
          return(c(mean = m, lower = NA_real_, upper = NA_real_, n = 1))
        }
        se <- stats::sd(xv) / base::sqrt(n)
        err <- stats::qt((1 + ci_level) / 2, df = n - 1) * se
        c(mean = m, lower = m - err, upper = m + err, n = n)
      }
    )
    mean_plot_df <- mean_plot_agg
    score_plot_stats <- mean_plot_agg$score_plot
    mean_plot_df$score_plot <- .hc_as_numeric_safely(score_plot_stats[, "mean"])
    mean_plot_df$ci_lower_plot <- .hc_as_numeric_safely(score_plot_stats[, "lower"])
    mean_plot_df$ci_upper_plot <- .hc_as_numeric_safely(score_plot_stats[, "upper"])
    mean_plot_df$n_donor_plot <- .hc_as_integer_safely(score_plot_stats[, "n"])
    mean_plot_df$score_plot[!base::is.finite(mean_plot_df$score_plot)] <- NA_real_
    mean_plot_df$ci_lower_plot[!base::is.finite(mean_plot_df$ci_lower_plot)] <- NA_real_
    mean_plot_df$ci_upper_plot[!base::is.finite(mean_plot_df$ci_upper_plot)] <- NA_real_
    mean_keys <- pair_key(mean_df$module, mean_df$term)
    mean_plot_keys <- pair_key(mean_plot_df$module, mean_plot_df$term)
    mean_df$term_display <- base::factor(
      safe_term_display(mean_df$module, mean_df$term),
      levels = base::levels(donor_df$term_display)
    )
    mean_plot_df$term_display <- base::factor(
      safe_term_display(mean_plot_df$module, mean_plot_df$term),
      levels = base::levels(donor_df$term_display)
    )
    mean_df$term_rank_label <- rank_label_map[mean_keys]
    mean_df$panel_label <- .hc_longitudinal_panel_label(
      module = mean_df$module,
      rank_label = mean_df$term_rank_label,
      term_display = mean_df$term_display
    )
    mean_plot_df$term_rank_label <- rank_label_map[mean_plot_keys]
    mean_plot_df$panel_label <- .hc_longitudinal_panel_label(
      module = mean_plot_df$module,
      rank_label = mean_plot_df$term_rank_label,
      term_display = mean_plot_df$term_display
    )

    panel_lookup <- ord_df[, c("module", "module_color", "term", "term_rank", "term_rank_label", "term_label"), drop = FALSE]
    panel_lookup$term_display <- safe_term_display(panel_lookup$module, panel_lookup$term)
    panel_lookup <- panel_lookup[
      is_valid_panel_field(panel_lookup$module) &
        is_valid_panel_field(panel_lookup$term) &
        is_valid_panel_field(panel_lookup$term_display), ,
      drop = FALSE
    ]
    present_keys <- base::unique(base::c(
      pair_key(donor_df$module, donor_df$term),
      pair_key(mean_df$module, mean_df$term),
      pair_key(mean_plot_df$module, mean_plot_df$term)
    ))
    present_keys <- present_keys[is_valid_panel_field(present_keys)]
    panel_lookup$key <- pair_key(panel_lookup$module, panel_lookup$term)
    panel_lookup <- panel_lookup[panel_lookup$key %in% present_keys, , drop = FALSE]
    if (base::nrow(panel_lookup) > 0) {
      empty_panel_lookup <- panel_lookup[0, , drop = FALSE]
      panel_lookup <- base::lapply(base::split(panel_lookup, panel_lookup$module), function(x) {
        x <- x[base::order(x$term_rank, x$term), , drop = FALSE]
        x <- x[!base::duplicated(x$term), , drop = FALSE]
        x$term_rank <- base::seq_len(base::nrow(x))
        x$term_rank_label <- base::paste0("Top ", x$term_rank)
        x$panel_label <- .hc_longitudinal_panel_label(
          module = x$module,
          rank_label = x$term_rank_label,
          term_display = x$term_display
        )
        x
      })
      panel_lookup <- panel_lookup[!base::vapply(panel_lookup, function(x) base::is.null(x) || base::nrow(x) == 0, FUN.VALUE = base::logical(1))]
      if (base::length(panel_lookup) > 0) {
        panel_lookup <- base::do.call(base::rbind, panel_lookup)
      } else {
        panel_lookup <- empty_panel_lookup
      }
    }
    if (base::nrow(panel_lookup) == 0) {
      next
    }
    panel_lookup <- panel_lookup[base::order(panel_lookup$module, panel_lookup$term_rank, panel_lookup$term), , drop = FALSE]
    panel_lookup <- panel_lookup[is_valid_panel_field(panel_lookup$panel_label), , drop = FALSE]
    if (base::nrow(panel_lookup) == 0) {
      next
    }
    lookup_cols <- c("term_rank", "term_rank_label", "term_display", "panel_label")
    donor_df <- donor_df[, base::setdiff(base::colnames(donor_df), lookup_cols), drop = FALSE]
    mean_df <- mean_df[, base::setdiff(base::colnames(mean_df), lookup_cols), drop = FALSE]
    mean_plot_df <- mean_plot_df[, base::setdiff(base::colnames(mean_plot_df), lookup_cols), drop = FALSE]
    panel_merge <- panel_lookup[, c("module", "term", "term_rank", "term_rank_label", "term_display", "panel_label"), drop = FALSE]
    donor_df <- base::merge(donor_df, panel_merge, by = c("module", "term"), all.x = FALSE, sort = FALSE)
    mean_df <- base::merge(mean_df, panel_merge, by = c("module", "term"), all.x = FALSE, sort = FALSE)
    mean_plot_df <- base::merge(mean_plot_df, panel_merge, by = c("module", "term"), all.x = FALSE, sort = FALSE)
    if (base::nrow(mean_df) == 0 || base::nrow(mean_plot_df) == 0) {
      next
    }
    panel_levels <- base::unique(base::as.character(panel_lookup$panel_label))
    panel_levels <- panel_levels[is_valid_panel_field(panel_levels)]
    if (base::length(panel_levels) == 0) {
      next
    }
    term_display_levels <- base::unique(base::as.character(panel_lookup$term_display))
    term_display_levels <- term_display_levels[is_valid_panel_field(term_display_levels)]
    donor_df <- .hc_align_longitudinal_panel_factors(
      df = donor_df,
      panel_levels = panel_levels,
      term_display_levels = term_display_levels
    )
    mean_df <- .hc_align_longitudinal_panel_factors(
      df = mean_df,
      panel_levels = panel_levels,
      term_display_levels = term_display_levels
    )
    mean_plot_df <- .hc_align_longitudinal_panel_factors(
      df = mean_plot_df,
      panel_levels = panel_levels,
      term_display_levels = term_display_levels
    )
    if (base::nrow(donor_df) == 0 || base::nrow(mean_df) == 0 || base::nrow(mean_plot_df) == 0) {
      next
    }
    panels_per_module <- base::max(panel_lookup$term_rank, na.rm = TRUE)
    if (!base::is.finite(panels_per_module) || panels_per_module < 1) {
      panels_per_module <- 1
    }
    modules_per_row <- base::max(1L, as.integer(panels_per_module))
    facet_ncol <- base::max(1L, as.integer(panels_per_module * modules_per_row))
    panel_col <- stats::setNames(
      base::as.character(panel_lookup$module_color),
      base::as.character(panel_lookup$panel_label)
    )
    panel_col[base::is.na(panel_col) | !base::nzchar(panel_col)] <- "grey60"

    panel_limits <- .hc_panel_mean_limits(
      mean_df = if (identical(score_scale, "z")) mean_plot_df else mean_df,
      panel_col = "panel_label",
      value_col = if (identical(score_scale, "z")) "score_plot" else "score",
      lower_col = if (isTRUE(show_ci) && identical(score_scale, "z")) "ci_lower_plot" else if (isTRUE(show_ci)) "ci_lower" else NULL,
      upper_col = if (isTRUE(show_ci) && identical(score_scale, "z")) "ci_upper_plot" else if (isTRUE(show_ci)) "ci_upper" else NULL,
      pad_frac = 0.12,
      min_span = 0.04
    )
    if (base::nrow(panel_limits) > 0) {
      panel_limits <- panel_limits[
        base::as.character(panel_limits$panel_label) %in% panel_levels, ,
        drop = FALSE
      ]
    }
    if (base::nrow(panel_limits) > 0) {
      donor_df <- base::merge(
        donor_df,
        panel_limits,
        by = "panel_label",
        all.x = TRUE,
        sort = FALSE
      )
      mean_df <- base::merge(
        mean_df,
        panel_limits,
        by = "panel_label",
        all.x = TRUE,
        sort = FALSE
      )
      lim_ok <- base::is.finite(donor_df$y_lower) & base::is.finite(donor_df$y_upper)
      donor_df$score_plot[lim_ok] <- base::pmin(
        base::pmax(donor_df$score_plot[lim_ok], donor_df$y_lower[lim_ok]),
        donor_df$y_upper[lim_ok]
      )
    } else {
      mean_df$y_lower <- NA_real_
      mean_df$y_upper <- NA_real_
    }
    ribbon_df <- if (identical(score_scale, "z")) mean_plot_df else mean_df
    if (identical(score_scale, "z")) {
      ribbon_df$ci_lower_use <- ribbon_df$ci_lower_plot
      ribbon_df$ci_upper_use <- ribbon_df$ci_upper_plot
    } else {
      ribbon_df$ci_lower_use <- ribbon_df$ci_lower
      ribbon_df$ci_upper_use <- ribbon_df$ci_upper
    }

    p <- ggplot2::ggplot()
    if (isTRUE(show_ci)) {
      p <- p + ggplot2::geom_ribbon(
        data = ribbon_df,
        mapping = ggplot2::aes(
          x = time,
          ymin = ci_lower_use,
          ymax = ci_upper_use,
          group = meta_cluster,
          fill = meta_cluster
        ),
        alpha = ci_alpha,
        color = NA,
        show.legend = FALSE,
        na.rm = TRUE
      )
    }
    if (isTRUE(show_donor_lines)) {
      p <- p + ggplot2::geom_line(
        data = donor_df,
        mapping = ggplot2::aes(
          x = time,
          y = score_plot,
          group = interaction(donor, meta_cluster),
          color = meta_cluster
        ),
        alpha = donor_alpha,
        linewidth = donor_linewidth,
        show.legend = FALSE,
        na.rm = TRUE
      )
    }
    p <- p +
      ggplot2::geom_line(
        data = mean_plot_df,
        mapping = ggplot2::aes(
          x = time,
          y = score_plot,
          group = meta_cluster,
          color = meta_cluster
        ),
        linewidth = mean_linewidth,
        na.rm = TRUE
      ) +
      ggplot2::geom_point(
        data = mean_plot_df,
        mapping = ggplot2::aes(
          x = time,
          y = score_plot,
          group = meta_cluster,
          color = meta_cluster
        ),
        size = mean_point_size,
        na.rm = TRUE
      ) +
      .hc_module_facet_wrap(
        module_levels = panel_levels,
        mod_col = panel_col,
        ncol = facet_ncol,
        scales = if (isTRUE(free_y)) "free_y" else "fixed",
        facets = ~panel_label
      ) +
      ggplot2::scale_fill_manual(values = meta_pal, drop = FALSE, guide = "none") +
      ggplot2::scale_color_manual(values = meta_pal, drop = FALSE) +
      .hc_theme_pub(base_size = 11) +
      ggplot2::theme(
        aspect.ratio = 1,
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        legend.position = "right"
      ) +
      ggplot2::guides(
        color = ggplot2::guide_legend(
          title = "Meta cluster",
          override.aes = list(alpha = 1, linewidth = 1.2, size = 2.2)
        )
      ) +
      ggplot2::labs(
        title = base::paste0(db_nm, " enrichment trajectories by meta-cluster"),
        x = "Timepoint",
        y = if (identical(score_scale, "z")) "Mean enrichment score (z-scaled)" else "Mean enrichment score"
      )
    if (!base::is.null(custom_terms)) {
      p <- p + ggplot2::labs(title = base::paste0(db_nm, " selected-term enrichment trajectories by meta-cluster"))
    }
    p <- .hc_colorize_module_strips(
      p = p,
      module_levels = panel_levels,
      mod_col = panel_col
    )

    db_token <- gsub("[^A-Za-z0-9]+", "_", db_nm)
    db_token <- gsub("^_+|_+$", "", db_token)
    if (!base::nzchar(db_token)) {
      db_token <- "DB"
    }

    if (isTRUE(save_pdf)) {
      out_dir <- .hc_resolve_output_dir(hc)
      .hc_ggsave_pdf_png(
        filename = base::file.path(out_dir, base::paste0(file_prefix, "_", db_token, ".pdf")),
        plot = p,
        width = save_width,
        height = save_height,
        units = "in",
        device = grDevices::cairo_pdf
      )
    }

    top_export <- rank_lookup[, c("module", "module_color", "term_rank", "term_rank_label", "term_label", "term", "rank", "qvalue"), drop = FALSE]
    top_export$score_method <- base_out$score_method_used[[db_nm]]
    top_export$score_scale <- score_scale
    top_export$enrichment_table <- enrichment_table
    top_export <- top_export[
      base::order(
        base::match(base::as.character(top_export$module), module_levels),
        top_export$term_rank
      ), ,
      drop = FALSE
    ]
    base::rownames(top_export) <- NULL

    if (isTRUE(export_excel) && requireNamespace("openxlsx", quietly = TRUE)) {
      out_dir <- .hc_resolve_output_dir(hc)
      .hc_write_xlsx_atomic(
        x = list(
          top_terms = top_export,
          meta_mean = mean_df,
          meta_mean_plot = mean_plot_df,
          donor_traj = donor_df
        ),
        file = base::file.path(out_dir, base::paste0(file_prefix, "_", db_token, ".xlsx")),
        overwrite = TRUE
      )
    }

    out_plots[[db_nm]] <- p
    out_top[[db_nm]] <- top_export
    out_donor[[db_nm]] <- donor_df
    mean_df <- base::merge(
      mean_df,
      mean_plot_df,
      by = c("module", "module_color", "term", "term_rank", "term_rank_label", "term_label", "time", "meta_cluster", "term_display", "panel_label"),
      all.x = TRUE,
      sort = FALSE
    )
    mean_df$mean_enrichment_score_plot <- mean_df$score_plot
    mean_df$score_scale <- score_scale
    out_mean[[db_nm]] <- mean_df
    out_method[[db_nm]] <- base_out$score_method_used[[db_nm]]
  }

  if (base::length(out_plots) == 0) {
    stop("No database produced plottable longitudinal enrichment meta-wave plots.")
  }

  list(
    plots = out_plots,
    top_terms = out_top,
    donor_trajectories = out_donor,
    mean_trajectories = out_mean,
    score_method_used = out_method
  )
}

#' Plot endotype-to-meta-cluster flow as Sankey diagram
#'
#' @param hc A `HCoCenaExperiment`.
#' @param slot_name Satellite slot with `donor_cluster` and `meta_cluster`.
#' @param save_html Logical. If `TRUE`, save widget to HTML.
#' @param file_prefix Output filename prefix.
#'
#' @return A list with `sankey` htmlwidget (or `NULL` if unavailable) and
#'   `links` table.
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
#' p <- hc_plot_longitudinal_endotype_meta_flow(hc, save_html = FALSE)
#' @export
hc_plot_longitudinal_endotype_meta_flow <- function(hc,
                                                    slot_name = "longitudinal_endotypes",
                                                    save_html = TRUE,
                                                    file_prefix = "Longitudinal_Endotype_Meta_Flow") {
  if (!inherits(hc, "HCoCenaExperiment")) {
    stop("`hc` must be a `HCoCenaExperiment`.")
  }
  sat <- as.list(hc@satellite)
  obj <- sat[[slot_name]]
  if (is.null(obj) || is.null(obj$donor_cluster) || is.null(obj$meta_cluster)) {
    stop("Slot `", slot_name, "` must contain `donor_cluster` and `meta_cluster`.")
  }

  tmp <- base::merge(
    obj$donor_cluster[, c("donor", "endotype"), drop = FALSE],
    obj$meta_cluster[, c("donor", "meta_cluster"), drop = FALSE],
    by = "donor",
    all = FALSE
  )
  if (base::nrow(tmp) == 0) {
    stop("No overlapping donors between endotypes and meta-clusters.")
  }
  links <- stats::aggregate(
    donor ~ endotype + meta_cluster,
    data = tmp,
    FUN = base::length
  )
  base::colnames(links)[base::colnames(links) == "donor"] <- "value"

  sankey <- NULL
  if (requireNamespace("networkD3", quietly = TRUE)) {
    lhs <- base::sort(base::unique(base::as.character(links$endotype)))
    rhs <- base::sort(base::unique(base::as.character(links$meta_cluster)))
    node_names <- base::c(base::paste0("Endotype: ", lhs), base::paste0("Meta: ", rhs))
    nodes <- base::data.frame(name = node_names, stringsAsFactors = FALSE)
    src <- base::match(base::paste0("Endotype: ", links$endotype), nodes$name) - 1L
    trg <- base::match(base::paste0("Meta: ", links$meta_cluster), nodes$name) - 1L
    link_df <- base::data.frame(
      source = as.integer(src),
      target = as.integer(trg),
      value = as.numeric(links$value),
      stringsAsFactors = FALSE
    )
    sankey <- networkD3::sankeyNetwork(
      Links = link_df,
      Nodes = nodes,
      Source = "source",
      Target = "target",
      Value = "value",
      NodeID = "name",
      fontSize = 12,
      nodeWidth = 30
    )

    if (isTRUE(save_html) && requireNamespace("htmlwidgets", quietly = TRUE)) {
      out_dir <- .hc_resolve_output_dir(hc)
      htmlwidgets::saveWidget(
        widget = sankey,
        file = base::file.path(out_dir, base::paste0(file_prefix, ".html")),
        selfcontained = TRUE
      )
    }
  } else {
    warning("Package `networkD3` not installed; Sankey widget was not created.")
  }

  list(sankey = sankey, links = links)
}
