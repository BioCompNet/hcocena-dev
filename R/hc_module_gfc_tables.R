#' Resolve the save folder of an `HCoCenaExperiment`
#'
#' Combines `dir_output` and `save_folder` the same way the legacy drivers do.
#' @param hc A `HCoCenaExperiment`.
#' @return A single directory path.
#' @noRd
.hc_module_gfc_output_dir <- function(hc) {
  paths <- base::as.data.frame(hc@config@paths)
  gl <- base::as.data.frame(hc@config@global)
  dir_output <- base::as.character(paths$dir_output[[1]])
  save_folder <- if ("save_folder" %in% base::colnames(gl)) {
    base::as.character(gl$save_folder[[1]])
  } else {
    ""
  }
  if (base::is.na(dir_output) || !base::nzchar(dir_output)) {
    stop("No `dir_output` set. Run `hc_set_paths()` first.", call. = FALSE)
  }
  if (base::is.na(save_folder) || !base::nzchar(save_folder)) {
    return(dir_output)
  }
  base::file.path(dir_output, save_folder)
}

#' Signed fold change against a reference, clamped to the GFC range
#'
#' Mirrors `gtools::foldchange()` (`x / ref` when `x >= ref`, otherwise
#' `-ref / x`) followed by the same clamping `gfc_calc()` applies. Kept local so
#' the helper does not depend on `gtools` being attached.
#' @param x Numeric matrix or vector of values.
#' @param ref Numeric reference, recycled along the rows of `x`.
#' @param range_gfc Positive clamp limit.
#' @return Numeric object of the same shape as `x`.
#' @noRd
.hc_module_gfc_foldchange <- function(x, ref, range_gfc) {
  out <- base::ifelse(x >= ref, x / ref, -ref / x)
  base::pmin(base::pmax(out, -range_gfc), range_gfc)
}

#' Module GFC per group and per sample
#'
#' Returns the Group-Fold-Change values behind the module heatmap as tables:
#' once aggregated per group (the numbers shown in the heatmap) and once
#' resolved per individual sample.
#'
#' hCoCena itself only stores the group level (`module_gfc_means`, modules by
#' groups). This function computes the per-sample counterpart using hCoCena's
#' own GFC definition, so both tables live on the same scale:
#'
#' 1. when `data_in_log` is `TRUE`, values are anti-logged (base 2) and all
#'    averaging happens on the linear scale,
#' 2. group means (or the raw sample values) are formed per gene,
#' 3. the reference follows the `control_keyword` global setting, exactly as
#'    [hc_run_expression_analysis_2()] does. With `"none"` it is the unweighted
#'    mean of the **group** means -- not the mean across all samples, which
#'    differs whenever the groups are unequally sized. With a control keyword
#'    it is the mean of the matching control group, and that group is dropped
#'    from the group table just like in the heatmap,
#' 4. `foldchange(x, ref)` is `x / ref` when `x >= ref` and `-ref / x`
#'    otherwise, a signed ratio rather than a difference of logs,
#' 5. values are clamped to `+/- range_GFC`,
#' 6. the module value is the mean across the genes of that module.
#'
#' Because step 4 is non-linear and step 5 truncates, **the group table is not
#' the average of the sample table** (`fc(mean(x)) != mean(fc(x))`). Both are
#' correct but answer different questions: `by_group` is the value plotted in
#' the heatmap, `by_sample` shows how uniformly the samples carry it. Use
#' `by_sample` for sample-level statistics and `by_group` for figure labels.
#'
#' As a self-check the function compares its own `by_group` against the
#' `grp_means` stored by hCoCena and reports the largest absolute difference,
#' which is expected to be on the order of `1e-3` -- the rounding hCoCena
#' applies per gene.
#'
#' @param hc A `HCoCenaExperiment`.
#' @param layer Layer to use, either its index or its name. Defaults to the
#'   first layer.
#' @param group_col Annotation column used for grouping. Defaults to the
#'   `variable_of_interest` from the global settings.
#' @param modules Either `"all"` (default) or a vector of module labels.
#' @param annotate Logical; if `TRUE` (default), the layer annotation is
#'   appended to the long-format sample table so it can be grouped directly.
#' @param export_xlsx Logical; if `TRUE` (default), all tables are written to
#'   an Excel workbook.
#' @param out_dir Target folder for the workbook. If `NULL` (default), the save
#'   folder of `hc` is used.
#' @param file_name File name of the workbook.
#' @param verbose Logical; print the self-check and the export path.
#' @examples
#' hc <- hc_example_data("clustered")
#' hc <- hc_module_gfc_tables(hc, export_xlsx = FALSE)
#' hc_satellite(hc, "module_gfc_tables")$by_group
#' @return Updated `HCoCenaExperiment`. The tables are stored in
#'   `hc@satellite$module_gfc_tables` as `by_group`, `by_group_long`,
#'   `by_sample`, `by_sample_long`, plus `reference` and `settings`.
#' @seealso [hc_plot_cluster_heatmap()], which produces the group-level
#'   `module_gfc_means` table.
#' @export
hc_module_gfc_tables <- function(hc,
                                 layer = 1,
                                 group_col = NULL,
                                 modules = "all",
                                 annotate = TRUE,
                                 export_xlsx = TRUE,
                                 out_dir = NULL,
                                 file_name = "Module_GFC_group_and_sample.xlsx",
                                 verbose = TRUE) {
  if (!methods::is(hc, "HCoCenaExperiment")) {
    stop("`hc` must be a `HCoCenaExperiment`.", call. = FALSE)
  }
  if (!base::is.logical(annotate) || base::length(annotate) != 1 ||
    base::is.na(annotate)) {
    stop("`annotate` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!base::is.logical(export_xlsx) || base::length(export_xlsx) != 1 ||
    base::is.na(export_xlsx)) {
    stop("`export_xlsx` must be TRUE or FALSE.", call. = FALSE)
  }

  gl <- base::as.data.frame(hc@config@global)
  range_gfc <- base::as.numeric(gl$range_GFC[[1]])
  if (!base::is.finite(range_gfc) || range_gfc <= 0) {
    stop("`range_GFC` in the global settings must be a positive number.",
      call. = FALSE
    )
  }
  data_in_log <- base::isTRUE(base::as.logical(gl$data_in_log[[1]]))
  control_keyword <- if ("control" %in% base::colnames(gl)) {
    base::as.character(gl$control[[1]])
  } else {
    "none"
  }
  if (base::is.na(control_keyword) || !base::nzchar(control_keyword)) {
    control_keyword <- "none"
  }
  if (base::is.null(group_col)) {
    group_col <- base::as.character(gl$voi[[1]])
  }
  group_col <- base::as.character(group_col[[1]])

  layer_results <- hc@layer_results
  if (base::is.null(layer_results) || base::length(layer_results) == 0) {
    stop(
      "No layer results found. Run `hc_run_expression_analysis_1()` first.",
      call. = FALSE
    )
  }
  if (base::is.character(layer)) {
    if (!layer %in% base::names(layer_results)) {
      stop("Unknown `layer`: '", layer, "'.", call. = FALSE)
    }
  } else if (layer < 1 || layer > base::length(layer_results)) {
    stop("`layer` index out of range.", call. = FALSE)
  }

  mat <- base::as.matrix(layer_results[[layer]]@part1[["topvar"]])
  if (base::is.null(mat) || base::nrow(mat) == 0) {
    stop(
      "No `topvar` matrix in this layer. Run `hc_run_expression_analysis_1()` first.",
      call. = FALSE
    )
  }

  anno <- base::as.data.frame(SummarizedExperiment::colData(
    MultiAssayExperiment::experiments(hc@mae)[[layer]]
  ))
  anno <- anno[base::colnames(mat), , drop = FALSE]
  if (!group_col %in% base::colnames(anno)) {
    stop("`group_col` '", group_col, "' not found in the layer annotation.",
      call. = FALSE
    )
  }

  cluster_calc <- hc@integration@cluster
  cluster_info <- base::as.data.frame(cluster_calc[["cluster_information"]])
  if (base::nrow(cluster_info) == 0) {
    stop("No cluster information found. Run `hc_cluster_calculation()` first.",
      call. = FALSE
    )
  }
  label_map <- .hc_resolve_module_label_map_for_colors(
    label_map = cluster_calc[["module_label_map"]],
    module_colors = base::as.character(cluster_info$color)
  )
  module_label <- function(colour) {
    lbl <- if (!base::is.null(label_map)) label_map[[colour]] else NULL
    if (base::is.null(lbl) || base::length(lbl) == 0 || base::is.na(lbl[[1]]) ||
      !base::nzchar(base::as.character(lbl[[1]]))) {
      colour
    } else {
      base::as.character(lbl[[1]])
    }
  }
  module_genes <- stats::setNames(
    base::lapply(base::seq_len(base::nrow(cluster_info)), function(i) {
      base::unique(base::trimws(base::unlist(
        base::strsplit(base::as.character(cluster_info$gene_n[i]), ",")
      )))
    }),
    base::vapply(base::as.character(cluster_info$color), module_label,
      FUN.VALUE = base::character(1)
    )
  )
  if (!base::identical(modules, "all")) {
    module_genes <- module_genes[
      base::intersect(base::names(module_genes), base::as.character(modules))
    ]
  }
  if (base::length(module_genes) == 0) {
    stop("No modules selected.", call. = FALSE)
  }

  use_control <- !base::identical(base::tolower(control_keyword), "none")

  # The control branch of `GFC_calculation()` substitutes zeros before the
  # antilog; mirror that so both branches reproduce hCoCena exactly.
  log_mat <- mat
  if (use_control && data_in_log) {
    log_mat[log_mat == 0] <- 1
  }
  linear <- if (data_in_log) 2^log_mat else log_mat

  groups <- base::factor(base::as.character(anno[[group_col]]))
  group_means <- base::vapply(
    base::levels(groups),
    function(g) base::rowMeans(linear[, groups == g, drop = FALSE]),
    FUN.VALUE = base::numeric(base::nrow(linear))
  )

  control_group <- NA_character_
  if (use_control) {
    is_ctrl <- base::grepl(control_keyword, base::colnames(group_means),
      ignore.case = TRUE, fixed = FALSE
    )
    if (base::sum(is_ctrl) == 0) {
      stop(
        "The control keyword '", control_keyword, "' does not match any sample group. ",
        "Available groups: ",
        base::paste(base::colnames(group_means), collapse = ", "), ".",
        call. = FALSE
      )
    }
    if (base::sum(is_ctrl) > 1) {
      stop(
        "The control keyword '", control_keyword, "' matches more than one group: ",
        base::paste(base::colnames(group_means)[is_ctrl], collapse = ", "), ".",
        call. = FALSE
      )
    }
    if (base::ncol(group_means) < 2) {
      stop("Only one sample group present; no fold changes can be computed.",
        call. = FALSE
      )
    }
    control_group <- base::colnames(group_means)[is_ctrl]
    reference <- group_means[, is_ctrl]
    group_means <- group_means[, !is_ctrl, drop = FALSE]
  } else {
    reference <- base::rowMeans(group_means)
  }

  module_means <- function(gfc_mat) {
    vals <- base::vapply(module_genes, function(genes) {
      genes <- base::intersect(genes, base::rownames(gfc_mat))
      if (base::length(genes) == 0) {
        base::rep(NA_real_, base::ncol(gfc_mat))
      } else {
        base::colMeans(gfc_mat[genes, , drop = FALSE], na.rm = TRUE)
      }
    }, FUN.VALUE = base::numeric(base::ncol(gfc_mat)))
    # `vapply()` drops to a plain vector when a single column is left (e.g. two
    # groups with a control), which would transpose the table. Force the shape.
    out <- base::matrix(
      vals,
      nrow = base::ncol(gfc_mat),
      ncol = base::length(module_genes),
      dimnames = base::list(base::colnames(gfc_mat), base::names(module_genes))
    )
    base::t(out)
  }
  as_wide <- function(m) {
    base::data.frame(
      module = base::rownames(m),
      base::as.data.frame(m),
      check.names = FALSE,
      row.names = NULL,
      stringsAsFactors = FALSE
    )
  }
  as_long <- function(m, id_name) {
    out <- base::data.frame(
      module = base::rep(base::rownames(m), times = base::ncol(m)),
      id = base::rep(base::colnames(m), each = base::nrow(m)),
      gfc = base::as.vector(m),
      check.names = FALSE,
      row.names = NULL,
      stringsAsFactors = FALSE
    )
    base::names(out)[2] <- id_name
    out
  }

  by_group <- module_means(
    .hc_module_gfc_foldchange(group_means, reference, range_gfc)
  )
  by_sample <- module_means(
    .hc_module_gfc_foldchange(linear, reference, range_gfc)
  )

  by_sample_long <- as_long(by_sample, "sample")
  if (base::isTRUE(annotate)) {
    keep <- base::setdiff(base::colnames(anno), c("module", "sample", "gfc"))
    by_sample_long <- base::cbind(
      by_sample_long,
      anno[base::match(by_sample_long$sample, base::rownames(anno)), keep,
        drop = FALSE
      ]
    )
    base::rownames(by_sample_long) <- NULL
  }

  # `grp_means` concatenates the conditions of ALL layers, so pick the block
  # belonging to the requested layer before comparing.
  layer_index <- if (base::is.character(layer)) {
    base::match(layer, base::names(layer_results))
  } else {
    base::as.integer(layer)
  }
  layer_widths <- base::vapply(base::seq_along(layer_results), function(j) {
    gfc_j <- tryCatch(layer_results[[j]]@part2[["GFC_all_genes"]],
      error = function(e) NULL
    )
    if (base::is.null(gfc_j)) {
      NA_integer_
    } else {
      base::sum(base::colnames(base::as.data.frame(gfc_j)) != "Gene")
    }
  }, FUN.VALUE = base::integer(1))
  offset <- if (layer_index > 1 && !base::anyNA(layer_widths)) {
    base::sum(layer_widths[base::seq_len(layer_index - 1)])
  } else {
    0L
  }
  take <- offset + base::seq_len(base::ncol(by_group))

  max_abs_diff <- NA_real_
  if ("grp_means" %in% base::colnames(cluster_info)) {
    stored <- tryCatch(
      {
        m <- base::t(base::vapply(
          base::seq_len(base::nrow(cluster_info)),
          function(i) {
            vals <- base::as.numeric(base::strsplit(
              base::as.character(cluster_info$grp_means[i]), ","
            )[[1]])
            vals[take]
          },
          FUN.VALUE = base::numeric(base::ncol(by_group))
        ))
        if (base::ncol(by_group) == 1) {
          m <- base::matrix(base::as.vector(m), ncol = 1)
        }
        base::rownames(m) <- base::vapply(
          base::as.character(cluster_info$color), module_label,
          FUN.VALUE = base::character(1)
        )
        m
      },
      error = function(e) NULL
    )
    if (!base::is.null(stored)) {
      common <- base::intersect(base::rownames(stored), base::rownames(by_group))
      if (base::length(common) > 0) {
        max_abs_diff <- base::max(base::abs(
          stored[common, , drop = FALSE] - by_group[common, , drop = FALSE]
        ), na.rm = TRUE)
        if (base::isTRUE(verbose)) {
          base::message(
            "Module GFC tables: max. deviation from hCoCena's stored ",
            "`grp_means` = ", base::signif(max_abs_diff, 3)
          )
        }
      }
    }
  }

  result <- base::list(
    by_group = as_wide(by_group),
    by_group_long = as_long(by_group, group_col),
    by_sample = as_wide(by_sample),
    by_sample_long = by_sample_long,
    reference = reference,
    settings = base::list(
      layer = layer,
      group_col = group_col,
      control_keyword = control_keyword,
      control_group = control_group,
      range_GFC = range_gfc,
      data_in_log = data_in_log,
      n_modules = base::length(module_genes),
      n_samples = base::ncol(mat),
      max_abs_diff_vs_hcocena = max_abs_diff
    )
  )

  if (base::isTRUE(export_xlsx)) {
    if (base::is.null(out_dir)) {
      out_dir <- .hc_module_gfc_output_dir(hc)
    }
    out_dir <- base::as.character(out_dir[[1]])
    if (!base::dir.exists(out_dir)) {
      stop("Output folder does not exist: ", out_dir, call. = FALSE)
    }
    target <- base::file.path(out_dir, base::as.character(file_name[[1]]))
    tryCatch(
      .hc_write_xlsx_atomic(
        x = result[c(
          "by_group", "by_sample", "by_group_long", "by_sample_long"
        )],
        file = target,
        overwrite = TRUE
      ),
      error = function(e) {
        base::warning("Could not write ", base::basename(target), ": ",
          base::conditionMessage(e),
          call. = FALSE
        )
      }
    )
    if (base::isTRUE(verbose)) {
      base::message("Module GFC tables written to: ", target)
    }
  }

  sat <- tryCatch(base::as.list(hc@satellite), error = function(e) base::list())
  if (base::is.null(sat) || !base::is.list(sat)) {
    sat <- base::list()
  }
  sat[["module_gfc_tables"]] <- result
  hc@satellite <- S4Vectors::SimpleList(sat)
  methods::validObject(hc)
  hc
}
