#' Initialize an empty `HCoCenaExperiment`
#'
#' @examples
#' hc <- hc_init()
#' methods::is(hc, "HCoCenaExperiment")
#' @return A valid, empty `HCoCenaExperiment`.
#' @export
hc_init <- function() {
  methods::new("HCoCenaExperiment")
}

.hc_normalize_path_value <- function(x, arg_name) {
  if (base::isFALSE(x) || base::identical(x, FALSE)) {
    return(FALSE)
  }

  if (base::length(x) != 1) {
    stop("`", arg_name, "` must be a scalar path string or FALSE.")
  }

  x_chr <- base::as.character(x[[1]])
  if (base::is.na(x_chr) || !base::nzchar(x_chr)) {
    stop("`", arg_name, "` must be a non-empty path string or FALSE.")
  }

  x_chr
}

#' Set input/output paths in an `HCoCenaExperiment`
#'
#' @param hc A `HCoCenaExperiment`.
#' @param dir_count_data Path to count data.
#' @param dir_annotation Path to annotation data.
#' @param dir_reference_files Path to reference files.
#' @param dir_output Output directory.
#' @examples
#' hc <- hc_init()
#' hc <- hc_set_paths(
#'   hc,
#'   dir_count_data = FALSE,
#'   dir_annotation = FALSE,
#'   dir_reference_files = tempdir(),
#'   dir_output = tempdir()
#' )
#' as.data.frame(methods::slot(hc_config(hc), "paths"))
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_set_paths <- function(hc, dir_count_data, dir_annotation, dir_reference_files, dir_output) {
  if (!inherits(hc, "HCoCenaExperiment")) {
    stop("`hc` must be a `HCoCenaExperiment`.")
  }

  dir_count_data <- .hc_normalize_path_value(dir_count_data, "dir_count_data")
  dir_annotation <- .hc_normalize_path_value(dir_annotation, "dir_annotation")
  dir_reference_files <- .hc_normalize_path_value(dir_reference_files, "dir_reference_files")
  dir_output <- .hc_normalize_path_value(dir_output, "dir_output")

  hc@config@paths <- S4Vectors::DataFrame(
    dir_count_data = dir_count_data,
    dir_annotation = dir_annotation,
    dir_reference_files = dir_reference_files,
    dir_output = dir_output
  )
  methods::validObject(hc)
  hc
}

.hc_path_with_trailing_slash <- function(x) {
  x <- base::as.character(x[[1]])
  x <- base::gsub("\\\\", "/", x)
  if (base::grepl("/$", x)) x else base::paste0(x, "/")
}

.hc_scalar_or_null <- function(x, arg_name) {
  if (is.null(x)) {
    return(NULL)
  }
  if (base::length(x) != 1) {
    stop("`", arg_name, "` must be NULL or a scalar string.")
  }
  out <- base::as.character(x[[1]])
  if (base::is.na(out) || !base::nzchar(out)) {
    return(NULL)
  }
  out
}

.hc_find_local_dir <- function(project_root, subdir) {
  if (is.null(subdir) || !base::nzchar(subdir)) {
    return(NULL)
  }
  candidate <- base::file.path(project_root, subdir)
  if (!base::dir.exists(candidate)) {
    return(NULL)
  }
  .hc_path_with_trailing_slash(base::normalizePath(candidate, winslash = "/", mustWork = TRUE))
}

#' Resolve hCoCena input/output directories for Docker or local projects
#'
#' @param project_root_hint Optional local project root to prioritize if it exists.
#'   If missing or not found, the current working directory is used.
#' @param local_count_subdir Local count-data subfolder below `project_root_hint`.
#' @param local_annotation_subdir Local annotation-data subfolder below `project_root_hint`.
#' @param local_reference_subdir Local reference-files subfolder below `project_root_hint`.
#' @param local_output_subdir Local output subfolder below `project_root_hint`.
#' @param docker_reference_dir Docker reference-files directory used to detect Docker mode.
#' @param docker_count_dir Docker count-data directory.
#' @param docker_annotation_dir Docker annotation-data directory.
#' @param docker_output_dir Docker output directory.
#' @param create_docker_dirs Boolean. If TRUE and Docker mode is detected, create the
#'   Docker directories when missing.
#' @param fallback_count_data Fallback count-data directory string used outside Docker
#'   when no local folder was detected.
#' @param fallback_annotation Fallback annotation directory string used outside Docker
#'   when no local folder was detected.
#' @param fallback_reference_files Fallback reference-files directory string used outside
#'   Docker when no local folder was detected.
#' @param fallback_output Fallback output directory string used outside Docker when no
#'   local folder was detected.
#' @examples
#' root <- file.path(tempdir(), "hcocena-paths")
#' dir.create(file.path(root, "reference_files"), recursive = TRUE, showWarnings = FALSE)
#' dir.create(file.path(root, "output"), recursive = TRUE, showWarnings = FALSE)
#' paths <- hc_resolve_paths(
#'   project_root_hint = root,
#'   local_count_subdir = NULL,
#'   local_annotation_subdir = NULL
#' )
#' names(paths)
#' @return Named list with `docker_mode`, `project_root`, `dir_count_data`,
#'   `dir_annotation`, `dir_reference_files`, and `dir_output`.
#' @export
hc_resolve_paths <- function(
  project_root_hint = NULL,
  local_count_subdir = "count_data",
  local_annotation_subdir = "annotation_data",
  local_reference_subdir = "reference_files",
  local_output_subdir = "output",
  docker_reference_dir = "/home/rstudio/reference_files",
  docker_count_dir = "/home/rstudio/count_data",
  docker_annotation_dir = "/home/rstudio/annotation_data",
  docker_output_dir = "/home/rstudio/output",
  create_docker_dirs = TRUE,
  fallback_count_data = "PATH_TO_FOLDER_THAT_HOLDS_YOUR_GENE_EXPRESSION_TABLES/",
  fallback_annotation = "PATH_TO_FOLDER_THAT_HOLDS_YOUR_ANNOTATION_TABLES/",
  fallback_reference_files = "PATH_TO_REFERENCE_FILES_FOLDER/",
  fallback_output = "PATH_TO_FOLDER_WHERE_HCOCENA_SHOULD_SAVE_ALL_ANALYSIS_OUTPUTS/"
) {
  project_root_hint <- .hc_scalar_or_null(project_root_hint, "project_root_hint")
  local_count_subdir <- .hc_scalar_or_null(local_count_subdir, "local_count_subdir")
  local_annotation_subdir <- .hc_scalar_or_null(local_annotation_subdir, "local_annotation_subdir")
  local_reference_subdir <- .hc_scalar_or_null(local_reference_subdir, "local_reference_subdir")
  local_output_subdir <- .hc_scalar_or_null(local_output_subdir, "local_output_subdir")
  docker_reference_dir <- .hc_scalar_or_null(docker_reference_dir, "docker_reference_dir")
  docker_count_dir <- .hc_scalar_or_null(docker_count_dir, "docker_count_dir")
  docker_annotation_dir <- .hc_scalar_or_null(docker_annotation_dir, "docker_annotation_dir")
  docker_output_dir <- .hc_scalar_or_null(docker_output_dir, "docker_output_dir")
  fallback_count_data <- .hc_scalar_or_null(fallback_count_data, "fallback_count_data")
  fallback_annotation <- .hc_scalar_or_null(fallback_annotation, "fallback_annotation")
  fallback_reference_files <- .hc_scalar_or_null(fallback_reference_files, "fallback_reference_files")
  fallback_output <- .hc_scalar_or_null(fallback_output, "fallback_output")

  docker_mode <- !is.null(docker_reference_dir) && base::dir.exists(docker_reference_dir)
  if (isTRUE(docker_mode) && isTRUE(create_docker_dirs)) {
    for (x in base::c(docker_count_dir, docker_annotation_dir, docker_reference_dir, docker_output_dir)) {
      if (!is.null(x) && base::nzchar(x)) {
        base::dir.create(x, recursive = TRUE, showWarnings = FALSE)
      }
    }
  }

  project_root <- if (!is.null(project_root_hint) && base::dir.exists(project_root_hint)) {
    base::normalizePath(project_root_hint, winslash = "/", mustWork = TRUE)
  } else {
    base::normalizePath(".", winslash = "/", mustWork = TRUE)
  }

  local_count_dir <- .hc_find_local_dir(project_root, local_count_subdir)
  local_annotation_dir <- .hc_find_local_dir(project_root, local_annotation_subdir)
  local_reference_dir <- .hc_find_local_dir(project_root, local_reference_subdir)
  local_output_dir <- .hc_find_local_dir(project_root, local_output_subdir)

  dir_count_data <- if (!is.null(local_count_dir)) {
    local_count_dir
  } else if (isTRUE(docker_mode) && !is.null(docker_count_dir)) {
    .hc_path_with_trailing_slash(docker_count_dir)
  } else {
    .hc_path_with_trailing_slash(fallback_count_data)
  }

  dir_annotation <- if (!is.null(local_annotation_dir)) {
    local_annotation_dir
  } else if (isTRUE(docker_mode) && !is.null(docker_annotation_dir)) {
    .hc_path_with_trailing_slash(docker_annotation_dir)
  } else {
    .hc_path_with_trailing_slash(fallback_annotation)
  }

  dir_reference_files <- if (!is.null(local_reference_dir)) {
    local_reference_dir
  } else if (isTRUE(docker_mode) && !is.null(docker_reference_dir)) {
    .hc_path_with_trailing_slash(docker_reference_dir)
  } else {
    .hc_path_with_trailing_slash(fallback_reference_files)
  }

  dir_output <- if (!is.null(local_output_dir)) {
    local_output_dir
  } else if (isTRUE(docker_mode) && !is.null(docker_output_dir)) {
    .hc_path_with_trailing_slash(docker_output_dir)
  } else {
    .hc_path_with_trailing_slash(fallback_output)
  }

  list(
    docker_mode = isTRUE(docker_mode),
    project_root = project_root,
    dir_count_data = dir_count_data,
    dir_annotation = dir_annotation,
    dir_reference_files = dir_reference_files,
    dir_output = dir_output
  )
}

#' Resolve and apply standard paths in one call
#'
#' @param hc A `HCoCenaExperiment`.
#' @param dir_count_data Optional explicit override for count-data directory.
#' @param dir_annotation Optional explicit override for annotation directory.
#' @param dir_reference_files Optional explicit override for reference-files directory.
#' @param dir_output Optional explicit override for output directory.
#' @inheritParams hc_resolve_paths
#' @examples
#' root <- file.path(tempdir(), "hcocena-auto-paths")
#' dir.create(file.path(root, "reference_files"), recursive = TRUE, showWarnings = FALSE)
#' dir.create(file.path(root, "output"), recursive = TRUE, showWarnings = FALSE)
#' hc <- hc_init()
#' hc <- hc_auto_set_paths(
#'   hc,
#'   project_root_hint = root,
#'   local_count_subdir = NULL,
#'   local_annotation_subdir = NULL
#' )
#' as.data.frame(methods::slot(hc_config(hc), "paths"))
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_auto_set_paths <- function(
  hc,
  project_root_hint = NULL,
  local_count_subdir = "count_data",
  local_annotation_subdir = "annotation_data",
  local_reference_subdir = "reference_files",
  local_output_subdir = "output",
  docker_reference_dir = "/home/rstudio/reference_files",
  docker_count_dir = "/home/rstudio/count_data",
  docker_annotation_dir = "/home/rstudio/annotation_data",
  docker_output_dir = "/home/rstudio/output",
  create_docker_dirs = TRUE,
  fallback_count_data = "PATH_TO_FOLDER_THAT_HOLDS_YOUR_GENE_EXPRESSION_TABLES/",
  fallback_annotation = "PATH_TO_FOLDER_THAT_HOLDS_YOUR_ANNOTATION_TABLES/",
  fallback_reference_files = "PATH_TO_REFERENCE_FILES_FOLDER/",
  fallback_output = "PATH_TO_FOLDER_WHERE_HCOCENA_SHOULD_SAVE_ALL_ANALYSIS_OUTPUTS/",
  dir_count_data = NULL,
  dir_annotation = NULL,
  dir_reference_files = NULL,
  dir_output = NULL
) {
  resolved <- hc_resolve_paths(
    project_root_hint = project_root_hint,
    local_count_subdir = local_count_subdir,
    local_annotation_subdir = local_annotation_subdir,
    local_reference_subdir = local_reference_subdir,
    local_output_subdir = local_output_subdir,
    docker_reference_dir = docker_reference_dir,
    docker_count_dir = docker_count_dir,
    docker_annotation_dir = docker_annotation_dir,
    docker_output_dir = docker_output_dir,
    create_docker_dirs = create_docker_dirs,
    fallback_count_data = fallback_count_data,
    fallback_annotation = fallback_annotation,
    fallback_reference_files = fallback_reference_files,
    fallback_output = fallback_output
  )

  final_count <- if (is.null(dir_count_data)) {
    resolved$dir_count_data
  } else {
    .hc_normalize_path_value(dir_count_data, "dir_count_data")
  }
  final_anno <- if (is.null(dir_annotation)) {
    resolved$dir_annotation
  } else {
    .hc_normalize_path_value(dir_annotation, "dir_annotation")
  }
  final_ref <- if (is.null(dir_reference_files)) {
    resolved$dir_reference_files
  } else {
    .hc_normalize_path_value(dir_reference_files, "dir_reference_files")
  }
  final_out <- if (is.null(dir_output)) {
    resolved$dir_output
  } else {
    .hc_normalize_path_value(dir_output, "dir_output")
  }

  hc_set_paths(
    hc,
    dir_count_data = final_count,
    dir_annotation = final_anno,
    dir_reference_files = final_ref,
    dir_output = final_out
  )
}

#' Define layers in an `HCoCenaExperiment`
#'
#' @param hc A `HCoCenaExperiment`.
#' @param data_sets Named list of layers, each value a pair (count, annotation).
#' @examples
#' hc <- hc_init()
#' hc <- hc_define_layers(
#'   hc,
#'   data_sets = list(
#'     Layer1 = c("counts.tsv", "anno.tsv")
#'   )
#' )
#' as.data.frame(methods::slot(hc_config(hc), "layer"))
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_define_layers <- function(hc, data_sets = list()) {
  if (!inherits(hc, "HCoCenaExperiment")) {
    stop("`hc` must be a `HCoCenaExperiment`.")
  }
  if (length(data_sets) == 0) {
    hc@config@layer <- S4Vectors::DataFrame()
    methods::validObject(hc)
    return(hc)
  }

  layer_names <- names(data_sets)
  if (is.null(layer_names)) {
    stop("`data_sets` must be a named list.")
  }

  rows <- vector("list", length(data_sets))
  for (i in seq_along(data_sets)) {
    pair <- data_sets[[i]]
    if (length(pair) < 2) {
      stop(
        "Layer `", layer_names[[i]], "` must provide exactly two entries: ",
        "count source and annotation source."
      )
    }

    count_src <- pair[[1]]
    anno_src <- pair[[2]]

    if (!is.character(count_src) || length(count_src) != 1 || is.na(count_src) || !nzchar(count_src)) {
      stop(
        "Layer `", layer_names[[i]], "`: count source must be a single string.\n",
        "If you use objects from the environment, pass quoted object names, e.g. ",
        "`c(\"counts_df\", \"anno_df\")` (not `c(counts_df, anno_df)`)."
      )
    }
    if (!is.character(anno_src) || length(anno_src) != 1 || is.na(anno_src) || !nzchar(anno_src)) {
      stop(
        "Layer `", layer_names[[i]], "`: annotation source must be a single string.\n",
        "If you use objects from the environment, pass quoted object names, e.g. ",
        "`c(\"counts_df\", \"anno_df\")` (not `c(counts_df, anno_df)`)."
      )
    }

    rows[[i]] <- list(
      layer_id = paste0("set", i),
      layer_name = as.character(layer_names[[i]]),
      count_source = as.character(count_src),
      annotation_source = as.character(anno_src)
    )
  }
  hc@config@layer <- .hc_rows_to_data_frame(rows)
  methods::validObject(hc)
  hc
}

#' Read expression and annotation data (S4 API)
#'
#' @noRd
# Internal helpers shared by modern and bridged legacy import code.
.hc_has_output_dir <- function(hc) {
  paths_cfg <- hc@config@paths
  if (!(base::nrow(paths_cfg) > 0 && "dir_output" %in% base::colnames(paths_cfg))) {
    return(FALSE)
  }

  out_dir <- as.character(paths_cfg$dir_output[[1]])
  base::length(out_dir) == 1 &&
    !base::is.na(out_dir) &&
    base::nzchar(out_dir) &&
    !base::identical(out_dir, "FALSE")
}

.hc_auto_prepare_output <- function(hc, create_output_dir = TRUE, project_folder = NULL) {
  if (!inherits(hc, "HCoCenaExperiment")) {
    stop("`hc` must be a `HCoCenaExperiment`.")
  }
  if (!isTRUE(.hc_has_output_dir(hc))) {
    return(hc)
  }

  hc <- hc_check_dirs(hc, create_output_dir = create_output_dir)

  save_folder <- project_folder
  if (is.null(save_folder)) {
    global_cfg <- hc@config@global
    if (base::nrow(global_cfg) > 0 && "save_folder" %in% base::colnames(global_cfg)) {
      current <- as.character(global_cfg$save_folder[[1]])
      if (base::length(current) == 1 && !base::is.na(current)) {
        save_folder <- current
      }
    }
  }

  if (is.null(save_folder) || base::length(save_folder) == 0 || base::is.na(save_folder)) {
    save_folder <- ""
  }
  save_folder <- as.character(save_folder[[1]])

  hc_init_save_folder(
    hc,
    name = save_folder,
    use_output_dir = base::identical(save_folder, "")
  )
}

.hc_normalize_layer_source <- function(src, kind, layer_idx) {
  if (!base::is.character(src) || base::length(src) != 1 || base::is.na(src) || !base::nzchar(src)) {
    stop(
      "Layer ", layer_idx, ": ", kind, " source must be a single non-empty string.\n",
      "If you provide objects, pass quoted names in `hc_define_layers()`, e.g. ",
      "`c(\"counts_df\", \"anno_df\")`."
    )
  }
  base::as.character(src)
}

.hc_legacy_dir_is_false <- function(x) {
  base::isFALSE(x) || (base::is.character(x) && base::identical(x, "FALSE"))
}

.hc_resolve_source_path <- function(dir_path, source_name) {
  if (base::file.exists(source_name)) {
    return(source_name)
  }
  if (.hc_legacy_dir_is_false(dir_path)) {
    return(source_name)
  }
  base::file.path(base::as.character(dir_path), source_name)
}

.hc_normalize_count_object <- function(obj, source_name, gene_symbol_col = NULL) {
  if (base::is.matrix(obj)) {
    obj <- base::as.data.frame(obj, stringsAsFactors = FALSE, check.names = FALSE)
  }
  if (!base::is.data.frame(obj)) {
    stop("Count source `", source_name, "` exists but is not a data.frame or matrix.")
  }
  obj <- base::as.data.frame(obj, stringsAsFactors = FALSE, check.names = FALSE)

  if (!base::is.null(gene_symbol_col) && gene_symbol_col %in% base::colnames(obj)) {
    obj <- make_rownames_unique(counts = obj, gene_symbol_col = gene_symbol_col)
  } else {
    if (base::is.null(base::rownames(obj)) || base::any(!base::nzchar(base::rownames(obj)))) {
      stop(
        "Count object `", source_name, "` has no usable rownames.\n",
        "Provide `gene_symbol_col` (matching a column in the count object) ",
        "or set rownames to gene symbols before calling `hc_read_data()`."
      )
    }
    numeric_cols <- base::vapply(obj, base::is.numeric, logical(1))
    if (!base::all(numeric_cols)) {
      obj <- obj[, numeric_cols, drop = FALSE]
    }
    if (base::ncol(obj) == 0) {
      stop("Count object `", source_name, "` has no numeric sample columns.")
    }
  }

  obj
}

.hc_normalize_anno_object <- function(obj, source_name, sample_col = NULL) {
  if (base::is.matrix(obj)) {
    obj <- base::as.data.frame(obj, stringsAsFactors = FALSE, check.names = FALSE)
  }
  if (!base::is.data.frame(obj)) {
    stop("Annotation source `", source_name, "` exists but is not a data.frame or matrix.")
  }
  obj <- base::as.data.frame(obj, stringsAsFactors = FALSE, check.names = FALSE)

  if (!base::is.null(sample_col) && sample_col %in% base::colnames(obj)) {
    base::rownames(obj) <- base::as.character(obj[[sample_col]])
  }
  if (base::is.null(base::rownames(obj)) || base::any(!base::nzchar(base::rownames(obj)))) {
    stop(
      "Annotation object `", source_name, "` has no usable rownames.\n",
      "Provide `sample_col` (matching a column in the annotation object) ",
      "or set rownames to sample IDs before calling `hc_read_data()`."
    )
  }

  obj[] <- base::lapply(obj, base::factor)
  obj
}

.hc_load_count_source <- function(source_name,
                                  dir_count_data,
                                  gene_symbol_col,
                                  count_has_rn,
                                  sep_counts,
                                  source_env) {
  if (base::exists(source_name, envir = source_env, inherits = TRUE)) {
    return(.hc_normalize_count_object(
      base::get(source_name, envir = source_env, inherits = TRUE),
      source_name = source_name,
      gene_symbol_col = gene_symbol_col
    ))
  }

  if (base::is.null(gene_symbol_col)) {
    stop("You must provide the 'gene_symbol_col' parameter.")
  }

  count_file <- .hc_resolve_source_path(dir_count_data, source_name)
  if (!base::file.exists(count_file)) {
    stop(
      "Count source `", source_name, "` was not found as an object, and file `",
      count_file, "` does not exist.\n",
      "If you use object input, pass quoted object names in `hc_define_layers()`."
    )
  }

  read_expression_data(
    file = count_file,
    rown = count_has_rn,
    sep = sep_counts,
    gene_symbol_col = gene_symbol_col
  )
}

.hc_load_anno_source <- function(source_name,
                                 dir_annotation,
                                 sample_col,
                                 anno_has_rn,
                                 sep_anno,
                                 source_env) {
  if (base::exists(source_name, envir = source_env, inherits = TRUE)) {
    return(.hc_normalize_anno_object(
      base::get(source_name, envir = source_env, inherits = TRUE),
      source_name = source_name,
      sample_col = sample_col
    ))
  }

  if (base::is.null(sample_col)) {
    stop("You must provide the 'sample_col' parameter.")
  }

  anno_file <- .hc_resolve_source_path(dir_annotation, source_name)
  if (!base::file.exists(anno_file)) {
    stop(
      "Annotation source `", source_name, "` was not found as an object, and file `",
      anno_file, "` does not exist.\n",
      "If you use object input, pass quoted object names in `hc_define_layers()`."
    )
  }

  read_anno(
    file = anno_file,
    rown = anno_has_rn,
    sep = sep_anno,
    sample_col = sample_col
  )
}

.hc_read_data_impl <- function(hc,
                               sep_counts = "\t",
                               sep_anno = "\t",
                               gene_symbol_col = NULL,
                               sample_col = NULL,
                               count_has_rn = TRUE,
                               anno_has_rn = TRUE,
                               auto_setup_output = TRUE,
                               create_output_dir = TRUE,
                               project_folder = NULL,
                               source_env = parent.frame()) {
  if (!inherits(hc, "HCoCenaExperiment")) {
    stop("`hc` must be a `HCoCenaExperiment`.")
  }
  if (isTRUE(auto_setup_output)) {
    hc <- .hc_auto_prepare_output(
      hc,
      create_output_dir = create_output_dir,
      project_folder = project_folder
    )
  }

  layer_cfg <- hc@config@layer
  required_cols <- c("layer_id", "count_source", "annotation_source")
  if (!(base::nrow(layer_cfg) > 0 && base::all(required_cols %in% base::colnames(layer_cfg)))) {
    stop("No valid layers found. Run `hc_define_layers()` before `hc_read_data()`.")
  }

  paths <- .hc_row_to_list(hc@config@paths)
  data <- list()
  for (i in base::seq_len(base::nrow(layer_cfg))) {
    lid <- as.character(layer_cfg$layer_id[[i]])
    count_source <- .hc_normalize_layer_source(layer_cfg$count_source[[i]], "count", i)
    anno_source <- .hc_normalize_layer_source(layer_cfg$annotation_source[[i]], "annotation", i)

    counts <- .hc_load_count_source(
      source_name = count_source,
      dir_count_data = paths[["dir_count_data"]],
      gene_symbol_col = gene_symbol_col,
      count_has_rn = count_has_rn,
      sep_counts = sep_counts,
      source_env = source_env
    )

    var.df <- rank_variance(counts)
    if (base::all(base::is.na(var.df$variance))) {
      stop(
        "Count data for layer ", i, " contains no valid numeric expression values after preprocessing."
      )
    }
    zero_var <- var.df$gene[!base::is.na(var.df$variance) & var.df$variance == 0]
    if (base::length(zero_var) > 0) {
      message("Detected genes with 0 variance in dataset ", i, ".")
      counts <- counts[!(base::rownames(counts) %in% zero_var), , drop = FALSE]
      message(base::length(zero_var), " gene(s) were removed from dataset ", i, ".")
    }

    anno <- .hc_load_anno_source(
      source_name = anno_source,
      dir_annotation = paths[["dir_annotation"]],
      sample_col = sample_col,
      anno_has_rn = anno_has_rn,
      sep_anno = sep_anno,
      source_env = source_env
    )

    if (!base::ncol(counts) == base::nrow(anno)) {
      stop(
        "The count table has ", base::ncol(counts), " columns but the annotation has ",
        base::nrow(anno), " rows. These values are required to be the same since they\n",
        "                 should correspond to the number of samples. THE LOADING OF THE DATA WILL BE TERMINATED."
      )
    }
    if (!base::all(base::as.character(base::colnames(counts)) %in% base::as.character(base::rownames(anno)))) {
      stop("The column names of the count file do not all match the rownames of the annotation. Please make sure they contain the same samples.")
    }
    anno <- anno[base::colnames(counts), , drop = FALSE]

    data[[base::paste0(lid, "_counts")]] <- counts
    data[[base::paste0(lid, "_anno")]] <- anno
  }

  previous <- tryCatch(.hc_assay_fingerprint(hc), error = function(e) NULL)
  hc@mae <- .hc_build_mae(list(data = data), layer_cfg)
  methods::validObject(hc)
  # Same gene ids with different values still passes validation, so compare the
  # contents rather than the shape.
  if (!identical(previous, tryCatch(.hc_assay_fingerprint(hc), error = function(e) NULL))) {
    hc <- .hc_invalidate_from(hc, from = "data")
  }
  hc
}
#
#' Function To Read All Count And Annotation Files
#'
#' The function loads the count and annotation data for each layer and saves it in the hCoCena-Object's "data" slot.
#' Genes with a variance of 0 will be automatically be excluded from the analysis.
#' @param hc A `HCoCenaExperiment`.
#' @param auto_setup_output Boolean. If `TRUE`, run `hc_check_dirs()` and
#'   initialize the save folder before reading data (if `dir_output` is set).
#' @param create_output_dir Boolean passed to `hc_check_dirs()` when
#'   `auto_setup_output = TRUE`.
#' @param project_folder Optional save subfolder name. Use `""` to write
#'   directly into `dir_output`. If `NULL`, keep an already configured
#'   `save_folder` or default to `""`.
#' @details Count inputs are normalized to numeric expression matrices during
#'   import. For object-based count inputs, non-numeric helper columns are
#'   dropped. Genes with zero variance are removed during loading because they
#'   can destabilize downstream variance-based heuristics such as
#'   `hc_suggest_topvar()`. As a result, suggested inflection points can differ
#'   slightly from older releases on the same raw input.
#' @examples
#' extdir <- paste0(
#'   normalizePath(system.file("extdata", package = "hcocena"), winslash = "/"),
#'   "/"
#' )
#' outdir <- file.path(tempdir(), "hcocena-read-data")
#' dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
#' hc <- hc_init()
#' hc <- hc_set_paths(
#'   hc,
#'   dir_count_data = extdir,
#'   dir_annotation = extdir,
#'   dir_reference_files = extdir,
#'   dir_output = outdir
#' )
#' hc <- hc_define_layers(
#'   hc,
#'   data_sets = list(
#'     Layer1 = c("toy_layer1_counts.tsv", "toy_layer1_anno.tsv"),
#'     Layer2 = c("toy_layer2_counts.tsv", "toy_layer2_anno.tsv")
#'   )
#' )
#' hc <- hc_read_data(
#'   hc,
#'   sep_counts = "\t",
#'   sep_anno = "\t",
#'   gene_symbol_col = "SYMBOL",
#'   sample_col = "SampleID",
#'   count_has_rn = FALSE,
#'   anno_has_rn = FALSE
#' )
#' names(MultiAssayExperiment::experiments(hc_mae(hc)))
#' @return Updated `HCoCenaExperiment`.
#' @param sep_counts The separator of the count files. Default is tab separated files. Ignore when loading data from objects instead of files.
#' @param sep_anno The separator of the annotation files. Default is tab separated files. Ignore when loading data from objects instead of files.
#' @param gene_symbol_col A String. Name of the column that contains the gene symbols. Ignore when loading data from objects instead of files.
#' @param sample_col A String. Name of the column that contains the sample IDs. Ignore when loading data from objects instead of files.
#' @param count_has_rn A Boolean. Whether or not the count file has rownames. Default is TRUE. Ignore when loading data from objects instead of files.
#' @param anno_has_rn A Boolean. Whether or not the annotation file has rownames. Default is TRUE. Ignore when loading data from objects instead of files.
#' @export
hc_read_data <- function(hc,
                         sep_counts = "\t",
                         sep_anno = "\t",
                         gene_symbol_col = NULL,
                         sample_col = NULL,
                         count_has_rn = TRUE,
                         anno_has_rn = TRUE,
                         auto_setup_output = TRUE,
                         create_output_dir = TRUE,
                         project_folder = NULL) {
  .hc_read_data_impl(
    hc = hc,
    sep_counts = sep_counts,
    sep_anno = sep_anno,
    gene_symbol_col = gene_symbol_col,
    sample_col = sample_col,
    count_has_rn = count_has_rn,
    anno_has_rn = anno_has_rn,
    auto_setup_output = auto_setup_output,
    create_output_dir = create_output_dir,
    project_folder = project_folder,
    source_env = parent.frame()
  )
}

#' Set global settings (S4 API)
#'
#' @noRd
# Internal implementation shared by S4 and legacy entry points.
.hc_set_global_settings_impl <- function(hc,
                                         organism = "human",
                                         control_keyword = "none",
                                         variable_of_interest = "merged",
                                         min_nodes_number_for_network = 50,
                                         min_nodes_number_for_cluster = 50,
                                         range_GFC = 2.0,
                                         layout_algorithm = "layout_with_stress",
                                         data_in_log = TRUE) {
  if (!inherits(hc, "HCoCenaExperiment")) {
    stop("`hc` must be a `HCoCenaExperiment`.")
  }

  layout_algorithm <- .hc_normalize_layout_algorithm(layout_algorithm)
  if (layout_algorithm %in% c("layout_with_stress", "layout_with_sparse_stress") &&
    !requireNamespace("graphlayouts", quietly = TRUE)) {
    warning(
      "Selected `layout_algorithm = '", layout_algorithm, "'`, but package `graphlayouts` is not installed.\n",
      "Network plotting will fall back to `layout_with_fr` until `graphlayouts` is available.",
      call. = FALSE
    )
  }

  global_cfg <- .hc_row_to_list(hc@config@global)
  global_cfg[["organism"]] <- organism
  global_cfg[["control"]] <- control_keyword
  global_cfg[["voi"]] <- variable_of_interest
  global_cfg[["min_nodes_number_for_network"]] <- min_nodes_number_for_network
  global_cfg[["min_nodes_number_for_cluster"]] <- min_nodes_number_for_cluster
  global_cfg[["range_GFC"]] <- range_GFC
  global_cfg[["layout_algorithm"]] <- layout_algorithm
  global_cfg[["data_in_log"]] <- data_in_log

  hc@config@global <- .hc_to_data_frame(global_cfg)
  methods::validObject(hc)
  hc
}
#
#' Define Global Settings
#'
#' Receives all settings that are globally valid, i.e., that are not dataset specific.
#' @param hc A `HCoCenaExperiment`.
#' @examples
#' extdir <- paste0(
#'   normalizePath(system.file("extdata", package = "hcocena"), winslash = "/"),
#'   "/"
#' )
#' outdir <- file.path(tempdir(), "hcocena-global-settings")
#' dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
#' hc <- hc_init()
#' hc <- hc_set_paths(
#'   hc,
#'   dir_count_data = extdir,
#'   dir_annotation = extdir,
#'   dir_reference_files = extdir,
#'   dir_output = outdir
#' )
#' hc <- hc_define_layers(
#'   hc,
#'   data_sets = list(
#'     Layer1 = c("toy_layer1_counts.tsv", "toy_layer1_anno.tsv"),
#'     Layer2 = c("toy_layer2_counts.tsv", "toy_layer2_anno.tsv")
#'   )
#' )
#' hc <- hc_read_data(
#'   hc,
#'   sep_counts = "\t",
#'   sep_anno = "\t",
#'   gene_symbol_col = "SYMBOL",
#'   sample_col = "SampleID",
#'   count_has_rn = FALSE,
#'   anno_has_rn = FALSE
#' )
#' hc <- hc_set_global_settings(
#'   hc,
#'   organism = "human",
#'   control_keyword = "control",
#'   variable_of_interest = "group",
#'   data_in_log = TRUE
#' )
#' as.data.frame(methods::slot(hc_config(hc), "global"))
#' @return Updated `HCoCenaExperiment`.
#' @param organism Specification of the organism needed for all TF-related analyses. Set to "human" for human data and "mouse" for mouse data. Currently, hCoCena only supports human and mouse data, however, all analysis steps not related to TFs can be applied to other organisms.
#' @param control_keyword Either 'none' if no controls are present (only possible when analysing one single dataset) or a string contained in the control sample descriptor of all annotation files, e.g. "healthy".
#' 	The string must only be contained in the descriptor, e.g. "healthy" would work for "rhinovirusSetHealthy" and "influenzaSetHealthy", it does not have to match it perfectly.
#' @param variable_of_interest The name of the column that mus be rpesent in all annotation files and which will be used for grouping samples, e.g., "condition"".
#' @param min_nodes_number_for_network An integer. The minimum number of nodes in the subsequently created network that can define a graph component. Graph components with less nodes will be discarded. Default is 50.
#' @param min_nodes_number_for_cluster An integer. The minimum number of nodes that constitute a module/cluster when detecting community structures in the network. Default is 50.
#' @param range_GFC A float. Defines the maximum value the group fold changes (GFCs) can acquire, all values above this value or beneath its negative will be truncated. Default is 2.0.
#' @param layout_algorithm Layout algorithm used for the network. Supported values are:
#'  `"layout_with_stress"` (default), `"layout_with_sparse_stress"`, `"layout_with_fr"`,
#'  `"layout_with_drl"`, `"layout_with_kk"` and `"cytoscape"`.
#'  The stress-based layouts require the optional `graphlayouts` package and are typically
#'  more stable/readable than Fruchterman-Reingold for medium/large graphs.
#'  `"cytoscape"` uses an externally calculated Cytoscape layout.
#' @param data_in_log Boolean. Whether or not the provided gene expression data is logged to the base of 2.
#' @export
hc_set_global_settings <- function(hc,
                                   organism = "human",
                                   control_keyword = "none",
                                   variable_of_interest = "merged",
                                   min_nodes_number_for_network = 50,
                                   min_nodes_number_for_cluster = 50,
                                   range_GFC = 2.0,
                                   layout_algorithm = "layout_with_stress",
                                   data_in_log = TRUE) {
  .hc_set_global_settings_impl(
    hc = hc,
    organism = organism,
    control_keyword = control_keyword,
    variable_of_interest = variable_of_interest,
    min_nodes_number_for_network = min_nodes_number_for_network,
    min_nodes_number_for_cluster = min_nodes_number_for_cluster,
    range_GFC = range_GFC,
    layout_algorithm = layout_algorithm,
    data_in_log = data_in_log
  )
}

#' Set layer settings (S4 API)
#'
#' @noRd
# Internal implementation shared by S4 and legacy entry points.
.hc_recycle_layer_arg <- function(x, n_layers, arg_name) {
  if (base::length(x) == 1 && n_layers > 1) {
    return(base::rep(x, n_layers))
  }
  if (base::length(x) != n_layers) {
    stop("`", arg_name, "` must have length 1 or match the number of layers (", n_layers, ").")
  }
  x
}

.hc_set_layer_settings_impl <- function(hc,
                                        top_var,
                                        min_corr = 0.7,
                                        range_cutoff_length,
                                        print_distribution_plots = FALSE) {
  if (!inherits(hc, "HCoCenaExperiment")) {
    stop("`hc` must be a `HCoCenaExperiment`.")
  }
  if (base::nrow(hc@config@layer) == 0 || !("layer_id" %in% base::colnames(hc@config@layer))) {
    stop("No layers found. Run `hc_define_layers()` before `hc_set_layer_settings()`.")
  }

  n_layers <- base::nrow(hc@config@layer)
  top_var <- .hc_recycle_layer_arg(top_var, n_layers, "top_var")
  min_corr <- .hc_recycle_layer_arg(min_corr, n_layers, "min_corr")
  range_cutoff_length <- .hc_recycle_layer_arg(range_cutoff_length, n_layers, "range_cutoff_length")
  print_distribution_plots <- .hc_recycle_layer_arg(
    print_distribution_plots,
    n_layers,
    "print_distribution_plots"
  )

  rows <- vector("list", n_layers)
  for (i in base::seq_len(n_layers)) {
    row <- list()
    for (nm in base::colnames(hc@config@layer)) {
      row[[nm]] <- hc@config@layer[[nm]][[i]]
    }
    row[["top_var"]] <- top_var[[i]]
    row[["min_corr"]] <- min_corr[[i]]
    row[["range_cutoff_length"]] <- range_cutoff_length[[i]]
    row[["print_distribution_plots"]] <- print_distribution_plots[[i]]
    rows[[i]] <- row
  }

  # Compare as numbers: these are stored with whatever type they were supplied
  # with, so identical() would call integer 40 and double 40 a change and throw
  # away a correct analysis.
  relevant <- c("top_var", "min_corr", "range_cutoff_length")
  settings_digest <- function(cfg) {
    keys <- base::intersect(relevant, base::colnames(cfg))
    if (base::length(keys) == 0) {
      return(NULL)
    }
    base::unlist(base::lapply(keys, function(k) base::as.numeric(cfg[[k]])))
  }
  before <- tryCatch(settings_digest(hc@config@layer), error = function(e) NULL)
  hc@config@layer <- .hc_rows_to_data_frame(rows)
  methods::validObject(hc)
  after <- tryCatch(settings_digest(hc@config@layer), error = function(e) NULL)
  if (!isTRUE(base::all.equal(before, after))) {
    hc <- .hc_invalidate_from(hc, from = "data")
  }
  hc
}
#
#' Define Layer Settings
#'
#' Receives all settings that are dataset specific and not globally valid.
#' @param hc A `HCoCenaExperiment`.
#' @examples
#' extdir <- paste0(
#'   normalizePath(system.file("extdata", package = "hcocena"), winslash = "/"),
#'   "/"
#' )
#' outdir <- file.path(tempdir(), "hcocena-layer-settings")
#' dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
#' hc <- hc_init()
#' hc <- hc_set_paths(
#'   hc,
#'   dir_count_data = extdir,
#'   dir_annotation = extdir,
#'   dir_reference_files = extdir,
#'   dir_output = outdir
#' )
#' hc <- hc_define_layers(
#'   hc,
#'   data_sets = list(
#'     Layer1 = c("toy_layer1_counts.tsv", "toy_layer1_anno.tsv"),
#'     Layer2 = c("toy_layer2_counts.tsv", "toy_layer2_anno.tsv")
#'   )
#' )
#' hc <- hc_read_data(
#'   hc,
#'   sep_counts = "\t",
#'   sep_anno = "\t",
#'   gene_symbol_col = "SYMBOL",
#'   sample_col = "SampleID",
#'   count_has_rn = FALSE,
#'   anno_has_rn = FALSE
#' )
#' hc <- hc_set_layer_settings(
#'   hc,
#'   top_var = c(2, 2),
#'   min_corr = 0.1,
#'   range_cutoff_length = c(2, 2),
#'   print_distribution_plots = FALSE
#' )
#' as.data.frame(methods::slot(hc_config(hc), "layer"))
#' @return Updated `HCoCenaExperiment`.
#' @param top_var A vector with the length equal to the number of datasets/layers. Each entry of the vector is either "all" or an integer.
#' 	Defines the number of most variable genes to be extracted per dataset. The order in the vector corresponds to the order in which the datasets have been declared in define_layers().
#' @param min_corr A vector of floats with the length equal to the number of datasets/layers. To construct a meaningful co-expression network for each layer, correlation cut-offs must be determined for every dataset that mark the lower boundary for the correlation
#' 	of two genes in order for their co-expression to be represented as an edge in the network. To facilitate choosing these cut-offs, a series of parameters will be calculated for a defined number of different cut-offs.
#' 	The number of cut-offs for which these parameters are calculated is determined by 'range_cutoff_length'.
#' 	The range from which these possible cut-offs are taken is on the lower end restricted by 'min_corr' and on the upper end by the maximum correlation calculated between any two genes in the dataset. Default is 0.7.
#' @param range_cutoff_length A vector of integers with the length equal to the number of datasets/layers. Details see "min_corr".
#' @param print_distribution_plots A vector of Booleans with the length equal to the number of datasets/layers. Whether or not to print the degree distribution plots for all tested cut-offs to pdf files.
#' 	The number of plots per data set will therefore be equal to the 'range_cutoff_length' parameter you have set.
#' 	Given the potential size, this should only be set to TRUE, if 'range_cutoff_length' is small or if one wishes to thoroughly analyse how the degree distribution changes in detail for differing cut-offs.
#' 	Default is FALSE.
#' @export
hc_set_layer_settings <- function(hc,
                                  top_var,
                                  min_corr = 0.7,
                                  range_cutoff_length,
                                  print_distribution_plots = FALSE) {
  .hc_set_layer_settings_impl(
    hc = hc,
    top_var = top_var,
    min_corr = min_corr,
    range_cutoff_length = range_cutoff_length,
    print_distribution_plots = print_distribution_plots
  )
}

#' Register supplementary references (S4 API)
#'
#' @noRd
# Internal implementation shared by S4 and legacy entry points.
.hc_set_supp_files_impl <- function(hc,
                                    Tf = NULL,
                                    Hallmark = NULL,
                                    Go = NULL,
                                    Kegg = NULL,
                                    Reactome = NULL,
                                    ...) {
  if (!inherits(hc, "HCoCenaExperiment")) {
    stop("`hc` must be a `HCoCenaExperiment`.")
  }

  hc@references@registry <- .hc_reference_registry(
    base::list(Tf = Tf, Hallmark = Hallmark, Go = Go, Kegg = Kegg, Reactome = Reactome, ...)
  )
  methods::validObject(hc)
  hc
}
#
#' Set Supplementary Files
#'
#' Receives the file names of the supplementary files
#' @param hc A `HCoCenaExperiment`.
#' @param ... Additional supplementary files passed through to `set_supp_files()`.
#' @examples
#' hc <- hc_init()
#' hc <- hc_set_supp_files(
#'   hc,
#'   Hallmark = "hallmark.gmt",
#'   Go = "go.gmt"
#' )
#' sort(names(hcocena:::as_hcobject(hc)$supplement))
#' @return Updated `HCoCenaExperiment`.
#' @param Tf A file name. The file must contain a list of at least two columns. One column should be titled with the organism from which the datasets originate and contain gene names, the other column must always be the last column and contain the information if the gene is a transcription factor ("TF"),
#' 	a co-factor ("Co_factor"), a chromatin remodelling protein ("Chromatin_remodeller") or a ribonucleic acid binding protein ("RNBP"). The name of the last column may vary. An exmemplary file for mouse and human is found in the Reference Files Folder in the repository.
#' @param Hallmark A file name for a .gmt MSigDB hallmark gene set file. You can find one in the reference file folder in the repository.
#' @param Go A file name for a .gmt gene ontology file. You can find one in the reference file folder in the repository.
#' @param Kegg A file name for a .gmt KEGG database file. You can find one in the reference file folder in the repository.
#' @param Reactome A file name for a .gmt Reactomte database file. You can find one in the reference file folder in the repository.
#' @export
hc_set_supp_files <- function(hc, Tf = NULL, Hallmark = NULL, Go = NULL, Kegg = NULL, Reactome = NULL, ...) {
  .hc_set_supp_files_impl(
    hc = hc,
    Tf = Tf,
    Hallmark = Hallmark,
    Go = Go,
    Kegg = Kegg,
    Reactome = Reactome,
    ...
  )
}

#' Load supplementary references (S4 API)
#'
#' @noRd
# Internal implementation shared by S4 and legacy entry points.
.hc_read_supplementary_target_name <- function(name) {
  if (base::identical(name, "Tf")) {
    return("TF")
  }
  if (name %in% c("Hallmark", "Go", "Kegg", "Reactome")) {
    return(name)
  }
  stringr::str_to_title(name)
}

.hc_read_supplementary_one <- function(name, source, dir_reference_files) {
  source_path <- .hc_resolve_source_path(dir_reference_files, source)
  if (!base::file.exists(source_path)) {
    stop("Supplementary source `", source, "` for `", name, "` does not exist.")
  }

  if (base::identical(name, "Tf")) {
    return(utils::read.delim(source_path, header = TRUE, check.names = FALSE))
  }
  if (name %in% c("Hallmark", "Go", "Kegg", "Reactome")) {
    return(clusterProfiler::read.gmt(source_path))
  }
  if (base::grepl("\\.csv$", source, ignore.case = TRUE)) {
    return(utils::read.csv(source_path, header = TRUE, stringsAsFactors = FALSE, quote = ""))
  }
  if (base::grepl("\\.gmt$", source, ignore.case = TRUE)) {
    return(clusterProfiler::read.gmt(source_path))
  }

  warning(
    "invalid input format of database: ", name, ". Valid inputs are: .csv and .gmt files!",
    call. = FALSE
  )
  NULL
}

.hc_read_supplementary_impl <- function(hc) {
  if (!inherits(hc, "HCoCenaExperiment")) {
    stop("`hc` must be a `HCoCenaExperiment`.")
  }

  registry <- hc@references@registry
  if (!(base::nrow(registry) > 0 && base::all(c("name", "source") %in% base::colnames(registry)))) {
    return(hc)
  }

  paths <- .hc_row_to_list(hc@config@paths)
  dir_reference_files <- paths[["dir_reference_files"]]
  if (base::is.null(dir_reference_files) || .hc_legacy_dir_is_false(dir_reference_files)) {
    stop("`dir_reference_files` is not set. Please run `hc_set_paths()` first.")
  }

  loaded <- as.list(hc@references@data)
  for (i in base::seq_len(base::nrow(registry))) {
    name <- as.character(registry$name[[i]])
    source <- as.character(registry$source[[i]])
    target_name <- .hc_read_supplementary_target_name(name)
    value <- .hc_read_supplementary_one(name, source, dir_reference_files)
    if (!base::is.null(value)) {
      loaded[[target_name]] <- value
    }
  }

  hc@references@data <- S4Vectors::SimpleList(loaded)
  methods::validObject(hc)
  hc
}
#
#' Read Supplementary Data
#'
#' Function to read and collect the supplementary files set in set_supp_files().
#' @param hc A `HCoCenaExperiment`.
#' @examples
#' refdir <- file.path(tempdir(), "hcocena-reference-files")
#' dir.create(refdir, recursive = TRUE, showWarnings = FALSE)
#' writeLines("ToySet\\tdescription\\tG1\\tG2", file.path(refdir, "hallmark.gmt"))
#' hc <- hc_init()
#' hc <- hc_set_paths(
#'   hc,
#'   dir_count_data = FALSE,
#'   dir_annotation = FALSE,
#'   dir_reference_files = paste0(normalizePath(refdir, winslash = "/"), "/"),
#'   dir_output = tempdir()
#' )
#' hc <- hc_set_supp_files(hc, Hallmark = "hallmark.gmt")
#' hc <- hc_read_supplementary(hc)
#' names(hcocena:::as_hcobject(hc)$supplementary_data)
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_read_supplementary <- function(hc) {
  .hc_read_supplementary_impl(hc)
}

#' Run expression analysis part I (S4 API)
#'
#' @noRd
# Internal implementation shared by S4 and legacy entry points.
.hc_run_expression_analysis_1_driver <- function(padj,
                                                 export,
                                                 import,
                                                 bayes,
                                                 prior,
                                                 alpha,
                                                 corr_method,
                                                 corr_backend) {
  for (x in base::seq_len(base::length(hcobject[["layers"]]))) {
    .hc_set_bridge_hcobject_slot(
      c("layer_specific_outputs", base::paste0("set", x), "part1"),
      run_expression_analysis_1_body(
        x = x,
        bayes = bayes,
        prior = prior,
        alpha = alpha,
        padj = padj,
        export = export,
        import = import,
        corr_method = corr_method,
        corr_backend = corr_backend
      )
    )
  }
}

.hc_run_expression_analysis_1_impl <- function(hc,
                                               padj = "none",
                                               export = FALSE,
                                               import = NULL,
                                               bayes = FALSE,
                                               prior = 2,
                                               alpha = 0.5,
                                               corr_method = "pearson",
                                               corr_backend = "auto") {
  if (!inherits(hc, "HCoCenaExperiment")) {
    stop("`hc` must be a `HCoCenaExperiment`.")
  }
  if (!base::is.character(corr_method) || base::length(corr_method) != 1L ||
    base::is.na(corr_method) || !corr_method %in% c("pearson", "spearman")) {
    stop("Parameter 'corr_method' must be either 'pearson' or 'spearman'.")
  }
  if (!base::is.character(corr_backend) || base::length(corr_backend) != 1L ||
    base::is.na(corr_backend) || !corr_backend %in% c("auto", "rcorr")) {
    stop("Parameter 'corr_backend' must be either 'auto' or 'rcorr'.")
  }
  if (!(base::nrow(hc@config@layer) > 0 && "layer_id" %in% base::colnames(hc@config@layer))) {
    stop("No layers found. Run `hc_define_layers()` before `hc_run_expression_analysis_1()`.")
  }

  .hc_run_driver(
    hc = hc,
    fun = .hc_run_expression_analysis_1_driver,
    padj = padj,
    export = export,
    import = import,
    bayes = bayes,
    prior = prior,
    alpha = alpha,
    corr_method = corr_method,
    corr_backend = corr_backend
  )
}
#
#' Run First Part Of The Gene Expression Analysis
#'
#' This function executes the frist part of the data processing procedure. It leads up to choosing the correlation cut-off for each layer.
#' 	All datasets will be filtered for their most variant genes as defined in the layer-specific settings.
#' 	After this filtering step, the pair-wise correlation coefficients for all pairs of genes are calculated.
#' 	Correlations that are negative or that have an associated p-value higher than 0.05 are immediately discarded.
#' 	Next, a set of statistics will be calculated for the set range of cut-off values that aim to facilitate the cut-off choice.
#' 	This includes determining the number of graph components resulting from creating a network when cutting the data with the respective cut-off,
#' 	as well as the number of nodes and edges this network comprises.
#' 	The last parameter that is evaluated is the R^2-value of the data to a linear regression through the logged degree distribution for the given network.
#' @param hc A `HCoCenaExperiment`.
#' @return Updated `HCoCenaExperiment`.
#' @param padj A String. Defines the method to be used for p-value adjustment. Valid values are "none" (default) or values for "method" in stats::p.adjust.
#' @param export A Boolean. If TRUE, correlation values and p-values will be exported.
#' 	This can save time if you plan on re-running the analysis since computing pari-wise correlations is a bottleneck of the analysis. Default is FALSE.
#' @param import A list. Each slot in the list corresponds to one of the layers (datasets) and is either a vector of two strings (1. path to file holding the correlation matrix and
#'  2. Path to the file holding the p-value matrix) or NA. A list slot is set to NA if for that layer you do not want to import a pre-calculated correlation matrix.
#'  The files do not necessarily have to be exported from a previous run, but can have any kind of origin (created with a different program or method).
#'  For compatibility it is only important, that it is a whitespace-separated text (.txt) file containing a symmetric, numeric matrix.
#'  The first line has to be gene names, there must be no row names, since the first line will be used for column names and row names.
#'  Also, you must provide a matrix with correlation values AND a matrix with corresponding p-values, where cells in the matrices correspond to each other (only a correlation matrix will not be sufficient).
#'  Default is NULL.
#' @param bayes Sanchez-Taltavull et al. (2016) suggest superiority of Bayesian correlation analysis to Pearson correlation in some cases.
#' 	Therefore, the Pearson correlation values can be weighted with Bayesian correlation values. To do so, set the "bayes"-parameter to TRUE. Default is FALSE, using only Pearson correlations.
#' @param alpha A numeric value from 0 to 1. Allows to adjust the strength of the Bayes weighting: For alpha = 0 the Pearson correlation values remain unaltered, for alpha = 1 the Pearson correlation value and the Bayesian correlation value contribute equally to the final correlation.
#' @param prior An integer, either 2 or 3, using prior 2 or 3 for the Bayes weighting as described in "Bayesian correlation analysis for sequence count data" by Sanchez-Taltavull et al. (2016).
#' @param corr_method Correlation method to use. Supported values are "pearson"
#'   (default) and "spearman".
#' @param corr_backend Correlation backend. `"auto"` (default) uses a fast
#'   cross-product path (BLAS matmul plus analytic t p-values) that is
#'   numerically identical to `Hmisc::rcorr` but ~13-25x faster, and
#'   automatically falls back to `Hmisc::rcorr` when the expression matrix
#'   contains `NA`s (to preserve pairwise-complete semantics). `"rcorr"` forces
#'   the original `Hmisc::rcorr` computation.
#' @examples
#' hc <- hc_example_data("prepared")
#' hc <- hc_run_expression_analysis_1(hc)
#' @export
hc_run_expression_analysis_1 <- function(hc,
                                         padj = "none",
                                         export = FALSE,
                                         import = NULL,
                                         bayes = FALSE,
                                         prior = 2,
                                         alpha = 0.5,
                                         corr_method = "pearson",
                                         corr_backend = "auto") {
  .hc_run_expression_analysis_1_impl(
    hc = hc,
    padj = padj,
    export = export,
    import = import,
    bayes = bayes,
    prior = prior,
    alpha = alpha,
    corr_method = corr_method,
    corr_backend = corr_backend
  )
}

.hc_set_cutoff_driver <- function(cutoff_vector) {
  .hc_set_bridge_hcobject_slot("cutoff_vec", cutoff_vector)
}

#' Set network cutoffs (S4 API)
#'
#' @param hc A `HCoCenaExperiment`.
#' @param cutoff_vector A vector of correlation cutoff values, one per layer,
#'   in the order the layers were declared in [hc_define_layers()]. Ignored
#'   when `auto = TRUE`.
#' @param auto Logical; if `TRUE`, automatically select cutoffs from available
#'   tuning outputs in this priority:
#'   1) `hc@satellite$cutoff_tuning$applied_cutoff_vector`,
#'   2) `hc@satellite$cutoff_tuning$recommended_cutoff_vector`,
#'   3) `hc@satellite$auto_tune$cutoff$recommended_cutoff_vector`,
#'   4) internal simple auto cutoff extraction from layer results,
#'   5) existing `hc@config@layer$cutoff`,
#'   6) user-provided `cutoff_vector`,
#'   7) `fallback_cutoff`.
#' @param fallback_cutoff Numeric fallback cutoff used for still-missing layers
#'   (for both `auto = TRUE` and `auto = FALSE`).
#' @param verbose Logical; print additional source details.
#'   The final applied cutoff vector is always printed.
#' @return Updated `HCoCenaExperiment`.
#' @examples
#' hc <- hc_example_data("after_part1")
#' hc <- hc_set_cutoff(hc, cutoff_vector = c(0.3, 0.3))
#' @export
hc_set_cutoff <- function(hc,
                          cutoff_vector = base::c(),
                          auto = FALSE,
                          fallback_cutoff = 0.982,
                          verbose = TRUE) {
  if (!inherits(hc, "HCoCenaExperiment")) {
    stop("`hc` must be a `HCoCenaExperiment`.")
  }
  if (!is.logical(auto) || length(auto) != 1 || is.na(auto)) {
    stop("`auto` must be TRUE or FALSE.")
  }
  if (!is.numeric(fallback_cutoff) || length(fallback_cutoff) != 1 || !is.finite(fallback_cutoff)) {
    stop("`fallback_cutoff` must be a finite numeric scalar.")
  }
  if (!is.logical(verbose) || length(verbose) != 1 || is.na(verbose)) {
    stop("`verbose` must be TRUE or FALSE.")
  }

  layer_ids <- character(0)
  if (base::nrow(hc@config@layer) > 0 && "layer_id" %in% base::colnames(hc@config@layer)) {
    layer_ids <- as.character(hc@config@layer$layer_id)
  } else if (base::length(hc@layer_results) > 0) {
    layer_ids <- base::names(hc@layer_results)
  }

  n_layers <- base::length(layer_ids)
  if (n_layers == 0 && base::nrow(hc@config@layer) > 0) {
    n_layers <- base::nrow(hc@config@layer)
    layer_ids <- paste0("set", seq_len(n_layers))
  }
  if (n_layers == 0 && base::length(hc@layer_results) > 0) {
    n_layers <- base::length(hc@layer_results)
    layer_ids <- if (!is.null(base::names(hc@layer_results))) {
      base::names(hc@layer_results)
    } else {
      paste0("set", seq_len(n_layers))
    }
  }
  if (n_layers == 0 && base::length(cutoff_vector) > 0) {
    n_layers <- base::length(cutoff_vector)
    layer_ids <- paste0("set", seq_len(n_layers))
  }
  if (n_layers == 0) {
    stop("No layers found to set cutoffs for.")
  }

  # Layers can be addressed by the name the user gave them in
  # hc_define_layers() or by the internal id. Both are accepted; a name that is
  # ambiguous between the two namespaces is an error rather than a guess.
  display_names <- tryCatch(
    as.character(hc@config@layer$layer_name),
    error = function(e) character(0)
  )
  if (length(display_names) != n_layers || anyNA(display_names)) {
    display_names <- character(0)
  }

  resolve_layer_name <- function(nm) {
    hits <- unique(c(which(layer_ids == nm),
                     if (length(display_names)) which(display_names == nm) else integer(0)))
    if (length(hits) == 0) {
      stop(
        "`cutoff_vector` names a layer that does not exist: `", nm, "`. ",
        "Known layers: ", paste(unique(c(display_names, layer_ids)), collapse = ", "), ".",
        call. = FALSE
      )
    }
    if (length(hits) > 1) {
      stop(
        "`cutoff_vector` name `", nm, "` is ambiguous between layers.",
        call. = FALSE
      )
    }
    hits[[1]]
  }

  check_range <- function(v, where) {
    bad <- v[is.finite(v) & (v < -1 | v > 1)]
    if (length(bad) > 0) {
      stop(
        "Correlation cutoffs must lie in [-1, 1]; ", where, " contains ",
        paste(unique(bad), collapse = ", "), ".",
        call. = FALSE
      )
    }
    invisible(NULL)
  }

  align_vec <- function(x, strict = FALSE) {
    out <- rep(NA_real_, n_layers)
    if (length(x) == 0) {
      return(out)
    }
    xv <- .hc_as_numeric_safely(x)
    if (isTRUE(strict)) {
      check_range(xv, "`cutoff_vector`")
    }

    nms <- names(x)
    has_names <- !is.null(nms) && any(nzchar(nms))

    if (has_names) {
      if (!all(nzchar(nms))) {
        if (isTRUE(strict)) {
          stop(
            "`cutoff_vector` mixes named and unnamed entries. Name every ",
            "element or none.",
            call. = FALSE
          )
        }
        return(out)
      }
      if (anyDuplicated(nms) > 0) {
        if (isTRUE(strict)) {
          stop(
            "`cutoff_vector` names a layer more than once: ",
            paste(unique(nms[duplicated(nms)]), collapse = ", "), ".",
            call. = FALSE
          )
        }
        return(out)
      }
      for (i in seq_along(nms)) {
        pos <- if (isTRUE(strict)) {
          resolve_layer_name(nms[[i]])
        } else {
          hits <- unique(c(which(layer_ids == nms[[i]]),
                           if (length(display_names)) which(display_names == nms[[i]]) else integer(0)))
          if (length(hits) == 1) hits[[1]] else NA_integer_
        }
        if (!is.na(pos)) {
          out[[pos]] <- xv[[i]]
        }
      }
      # Named input addresses exactly the layers it names. Positions left over
      # stay NA and are filled from the next source in the precedence list.
      return(out)
    }

    if (isTRUE(strict) && length(xv) != 1L && length(xv) != n_layers) {
      stop(
        "`cutoff_vector` has ", length(xv), " values for ", n_layers,
        " layers. Supply one value, one per layer, or a named vector.",
        call. = FALSE
      )
    }
    if (length(xv) == 1L) {
      out[] <- xv[[1]]
      return(out)
    }
    n_copy <- min(n_layers, length(xv))
    out[seq_len(n_copy)] <- xv[seq_len(n_copy)]
    out
  }

  cutoff_state <- new.env(parent = emptyenv())
  cutoff_state$chosen <- rep(NA_real_, n_layers)
  cutoff_state$source_per_layer <- rep(NA_character_, n_layers)

  fill_missing <- function(vec, source_name, strict = FALSE) {
    aligned <- align_vec(vec, strict = strict)
    idx <- which(!is.finite(cutoff_state$chosen) & is.finite(aligned))
    if (length(idx) > 0) {
      cutoff_state$chosen[idx] <- aligned[idx]
      cutoff_state$source_per_layer[idx] <- source_name
    }
    invisible(NULL)
  }

  if (isTRUE(auto)) {
    sat <- as.list(hc@satellite)
    if ("cutoff_tuning" %in% names(sat) && is.list(sat[["cutoff_tuning"]])) {
      ct <- sat[["cutoff_tuning"]]
      if ("applied_cutoff_vector" %in% names(ct)) {
        fill_missing(ct[["applied_cutoff_vector"]], "cutoff_tuning.applied")
      }
      if ("recommended_cutoff_vector" %in% names(ct)) {
        fill_missing(ct[["recommended_cutoff_vector"]], "cutoff_tuning.recommended")
      }
    }
    if ("auto_tune" %in% names(sat) &&
      is.list(sat[["auto_tune"]]) &&
      "cutoff" %in% names(sat[["auto_tune"]]) &&
      is.list(sat[["auto_tune"]][["cutoff"]]) &&
      "recommended_cutoff_vector" %in% names(sat[["auto_tune"]][["cutoff"]])) {
      fill_missing(sat[["auto_tune"]][["cutoff"]][["recommended_cutoff_vector"]], "auto_tune.recommended")
    }

    if (base::exists(".hc_auto_collect_cutoffs", mode = "function", inherits = TRUE)) {
      simple_auto <- tryCatch(.hc_auto_collect_cutoffs(hc), error = function(e) numeric(0))
      fill_missing(simple_auto, "simple_auto")
    }

    if (base::nrow(hc@config@layer) > 0 && "cutoff" %in% base::colnames(hc@config@layer)) {
      fill_missing(hc@config@layer$cutoff, "config.layer.cutoff")
    }

    fill_missing(cutoff_vector, "user.cutoff_vector", strict = TRUE)
  } else {
    fill_missing(cutoff_vector, "user.cutoff_vector", strict = TRUE)
    if (base::nrow(hc@config@layer) > 0 && "cutoff" %in% base::colnames(hc@config@layer)) {
      fill_missing(hc@config@layer$cutoff, "config.layer.cutoff")
    }
  }

  chosen <- cutoff_state$chosen
  source_per_layer <- cutoff_state$source_per_layer
  miss <- which(!is.finite(chosen))
  if (length(miss) > 0) {
    chosen[miss] <- as.numeric(fallback_cutoff)
    source_per_layer[miss] <- "fallback"
  }

  message(
    "hc_set_cutoff(auto=", if (isTRUE(auto)) "TRUE" else "FALSE", "): applying cutoffs = ",
    paste0(format(chosen, digits = 4), collapse = ", ")
  )
  if (isTRUE(verbose)) {
    src_df <- data.frame(
      layer_id = layer_ids,
      source = as.character(source_per_layer),
      stringsAsFactors = FALSE
    )
    message(
      "hc_set_cutoff(): source per layer -> ",
      paste0(src_df$layer_id, ":", src_df$source, collapse = ", ")
    )
  }

  previous_cutoffs <- tryCatch(as.numeric(hc@config@layer$cutoff), error = function(e) NULL)

  hc <- .hc_run_driver(
    hc = hc,
    fun = .hc_set_cutoff_driver,
    cutoff_vector = chosen
  )

  # A different cutoff means a different filtered network, so part2 onwards is
  # obsolete. The correlations in part1 do not depend on the cutoff and stay.
  if (!isTRUE(all.equal(previous_cutoffs, as.numeric(chosen)))) {
    hc <- .hc_invalidate_from(hc, from = "cutoff")
  }

  sat2 <- as.list(hc@satellite)
  sat2[["cutoff_selection"]] <- list(
    created_at = as.character(base::Sys.time()),
    mode = if (isTRUE(auto)) "auto" else "manual",
    fallback_cutoff = as.numeric(fallback_cutoff),
    applied_cutoff_vector = chosen,
    layer_source = data.frame(
      layer_id = layer_ids,
      cutoff = as.numeric(chosen),
      source = as.character(source_per_layer),
      stringsAsFactors = FALSE
    )
  )
  hc@satellite <- S4Vectors::SimpleList(sat2)
  methods::validObject(hc)
  hc
}

#' Run expression analysis part II (S4 API)
#'
#' @noRd
# Internal implementation shared by S4 and legacy entry points.
.hc_run_expression_analysis_2_driver <- function(grouping_v,
                                                 plot_HM,
                                                 method,
                                                 additional_anno,
                                                 cols) {
  for (x in base::seq_len(base::length(hcobject[["layers"]]))) {
    extra_anno <- if (base::is.null(additional_anno) || base::length(additional_anno) < x) {
      NULL
    } else {
      additional_anno[[x]]
    }

    run_expression_analysis_2_body(
      x = x,
      grouping_v = grouping_v,
      plot_HM = plot_HM,
      method = method,
      additional_anno = extra_anno,
      title = hcobject[["layers_names"]][x],
      cols = cols
    )
  }
}

.hc_run_expression_analysis_2_impl <- function(hc,
                                               grouping_v = NULL,
                                               plot_HM = TRUE,
                                               method = "complete",
                                               additional_anno = NULL,
                                               cols = NULL) {
  if (!inherits(hc, "HCoCenaExperiment")) {
    stop("`hc` must be a `HCoCenaExperiment`.")
  }
  if (!(base::nrow(hc@config@layer) > 0 && "layer_id" %in% base::colnames(hc@config@layer))) {
    stop("No layers found. Run `hc_define_layers()` before `hc_run_expression_analysis_2()`.")
  }

  .hc_run_driver(
    hc = hc,
    fun = .hc_run_expression_analysis_2_driver,
    grouping_v = grouping_v,
    plot_HM = plot_HM,
    method = method,
    additional_anno = additional_anno,
    cols = cols
  )
}
#
#' Run Second Part Of The Gene Expression Analysis
#'
#' This function plots a heatmap for the network genes in each data layer and computes the Group-Fold-Changes for all genes per layer.
#' @param hc A `HCoCenaExperiment`.
#' @return Updated `HCoCenaExperiment`.
#' @param grouping_v A string giving a column name present in all annotation files, if this variable shall be used for grouping the samles isntead of the variable of interest. Default is NULL.
#' @param plot_HM A Boolean. Whether or not to plot the heatmap (for networks with many genes this may be very demanding for your computer if you are running the analysis locally). Default is TRUE.
#' @param method The method used for clustering the heatmap in the pheatmap function. Default is "complete".
#' @param additional_anno A list, with one slot per data set. A slot contains a vector of column names from that data set's annotation file that you wish to annotate with.
#' 	If for some of the data sets you don't wish any further annotation, you can set the corresponding list slot to NULL. Default is NULL.
#' @param cols A named list of color vectors. The list names need to match the chosen annotation column names. Default is NULL which uses implemented colors.
#' @examples
#' hc <- hc_example_data("after_part1")
#' hc <- hc_run_expression_analysis_2(hc, plot_HM = FALSE)
#' @export
hc_run_expression_analysis_2 <- function(hc,
                                         grouping_v = NULL,
                                         plot_HM = TRUE,
                                         method = "complete",
                                         additional_anno = NULL,
                                         cols = NULL) {
  .hc_run_expression_analysis_2_impl(
    hc = hc,
    grouping_v = grouping_v,
    plot_HM = plot_HM,
    method = method,
    additional_anno = additional_anno,
    cols = cols
  )
}

#' Build integrated network (S4 API)
#'
#' @noRd
# Internal implementation shared by S4 and legacy entry points.
.hc_build_integrated_network_driver <- function(mode,
                                                with,
                                                multi_edges,
                                                GFC_when_missing) {
  if (mode == "u") {
    message("Intergrating network based on union.")
    get_union()
  } else if (mode == "i") {
    message("Intergrating network based on intersection.")
    if (is.null(with)) {
      stop("The 'with' parameter must be specified.")
    }
    get_intersection(with = with)
  } else {
    stop("No valid choice of 'mode'. Must be either 'u' for integration by union or 'i' for integration by intersection.")
  }

  merged_net <- igraph::graph_from_data_frame(
    hcobject[["integrated_output"]][["combined_edgelist"]],
    directed = FALSE
  )
  merged_net <- igraph::simplify(
    merged_net,
    edge.attr.comb = list(weight = multi_edges, "ignore")
  )
  new_edgelist <- base::cbind(
    igraph::get.edgelist(merged_net),
    base::round(igraph::E(merged_net)$weight, 7)
  ) %>%
    base::as.data.frame()
  base::colnames(new_edgelist) <- base::colnames(
    hcobject[["integrated_output"]][["combined_edgelist"]]
  )
  .hc_set_bridge_hcobject_slot(c("integrated_output", "combined_edgelist"), new_edgelist)
  .hc_set_bridge_hcobject_slot(c("integrated_output", "merged_net"), merged_net)
  .hc_set_bridge_hcobject_slot(
    c("integrated_output", "GFC_all_layers"),
    merge_GFCs(GFC_when_missing = GFC_when_missing)
  )

  # Rebuilding the integrated graph invalidates prior cluster-dependent outputs.
  .hc_set_bridge_hcobject_slot(c("integrated_output", "cluster_calc"), list())
  .hc_set_bridge_hcobject_slot(c("integrated_output", "enrichments"), NULL)
  .hc_set_bridge_hcobject_slot(c("integrated_output", "upstream_inference"), NULL)
  .hc_set_bridge_hcobject_slot(c("integrated_output", "knowledge_network"), NULL)
  .hc_set_bridge_hcobject_slot(c("satellite_outputs", "enrichments"), NULL)
  .hc_set_bridge_hcobject_slot(c("satellite_outputs", "upstream_inference"), NULL)
  .hc_set_bridge_hcobject_slot(c("satellite_outputs", "knowledge_network"), NULL)
  .hc_set_bridge_hcobject_slot(c("satellite_outputs", "labelled_network"), NULL)
  .hc_set_bridge_hcobject_slot(c("satellite_outputs", "network_col_by_module"), NULL)
}
#
#' Network integration
#'
#' The previously constructed layer-specific networks are being integrated.
#' @param hc A `HCoCenaExperiment`.
#' @return Updated `HCoCenaExperiment`.
#' @param mode A string, either "u", if network integration is to be done by union (default), or "i", if integration is to be done by intersection.
#'  For details please refer to the information pages provided in the repository's Wiki.
#' @param multi_edges One of "min", "mean" or "max" resulting in the simplification of the multigraph by using the minimum, the mean or the maximum edge weight among the multiple edges, respectively.
#'  Multiple edges occur when the edge was present in more than one datset.
#' @param GFC_when_missing The value to substitute missing data in the case where some genes were not measured in all but only some of the datasets.
#' @param with Either an integer giving the number of the dataset to be used as reference (e.g., 1) or the name given to the layer.
#'  Can be ignored when integration is done by union.
#' @examples
#' hc <- hc_example_data("after_part2")
#' hc <- hc_build_integrated_network(hc, mode = "u")
#' @export
hc_build_integrated_network <- function(hc,
                                        mode = "u",
                                        with = NULL,
                                        multi_edges = "min",
                                        GFC_when_missing = NULL) {
  if (is.null(GFC_when_missing) &&
    base::nrow(hc@config@global) > 0 &&
    "range_GFC" %in% base::colnames(hc@config@global)) {
    GFC_when_missing <- -hc@config@global$range_GFC[[1]]
  }
  if (is.null(GFC_when_missing)) {
    GFC_when_missing <- -2.0
  }
  .hc_run_driver(
    hc = hc,
    fun = .hc_build_integrated_network_driver,
    mode = mode,
    with = with,
    multi_edges = multi_edges,
    GFC_when_missing = GFC_when_missing
  )
}

#' Cluster integrated network (S4 API)
#'
#' @noRd
# Internal implementation shared by S4 and legacy entry points.
.hc_cluster_calculation_driver <- function(cluster_algo,
                                           no_of_iterations,
                                           resolution,
                                           partition_type,
                                           max_cluster_count_per_gene,
                                           return_result) {
  .hc_set_bridge_hcobject_slot(c("global_settings", "chosen_clustering_algo"), cluster_algo)

  if (cluster_algo == "cluster_leiden") {
    .hc_set_bridge_hcobject_slot(
      c("integrated_output", "cluster_calc", "cluster_information"),
      leiden_clustering(
        g = hcobject[["integrated_output"]][["merged_net"]],
        num_it = no_of_iterations,
        resolution = resolution,
        partition_type = partition_type
      )
    )
    .hc_set_bridge_hcobject_slot(c("integrated_output", "cluster_calc", "labelled_network"), NULL)
    .hc_set_bridge_hcobject_slot(c("integrated_output", "cluster_calc", "network_col_by_module"), NULL)
    .hc_set_bridge_hcobject_slot(c("satellite_outputs", "labelled_network"), NULL)
    .hc_set_bridge_hcobject_slot(c("satellite_outputs", "network_col_by_module"), NULL)
    return(invisible(NULL))
  }

  cluster_algo_list <- c(
    "cluster_label_prop",
    "cluster_fast_greedy",
    "cluster_louvain",
    "cluster_infomap",
    "cluster_walktrap",
    "cluster_leiden"
  )

  if (cluster_algo == "auto") {
    algos_to_use <- cluster_algo_list
  } else {
    algos_to_use <- cluster_algo
  }

  cluster_run <- .hc_with_seed(168575L, {
    if (cluster_algo == "auto") {
      message("Testing clustering algorithms: ", base::paste(algos_to_use, collapse = ", "))
      df_modularity_score <- base::do.call(
        "rbind",
        base::lapply(algos_to_use, function(algo_now) {
          cluster_calculation_internal(
            graph_obj = hcobject[["integrated_output"]][["merged_net"]],
            algo = algo_now,
            case = "test",
            resolution = resolution,
            partition_type = partition_type
          )
        })
      )

      cluster_algo_used <- df_modularity_score %>%
        dplyr::filter(modularity_score == base::max(modularity_score)) %>%
        dplyr::select(cluster_algorithm) %>%
        base::as.character()

      message(cluster_algo_used, " will be used based on the highest modularity score.")
    } else {
      cluster_algo_used <- cluster_algo
      message(cluster_algo_used, " will be used based on your input.")
    }

    # Deterministic algorithms return the identical partition every time, so
    # extra replicates only cost runtime and can never change the vote.
    iterations_to_run <- no_of_iterations
    if (cluster_algo_used %in% .hc_deterministic_cluster_algos() &&
      no_of_iterations > 1) {
      message(
        cluster_algo_used, " is deterministic; running 1 iteration instead of ",
        no_of_iterations, " (replicates would be identical)."
      )
      iterations_to_run <- 1L
    }

    gene_which_cluster <- base::do.call(
      "cbind",
      base::lapply(base::seq_len(iterations_to_run), function(iteration_now) {
        cluster_calculation_internal(
          graph_obj = hcobject[["integrated_output"]][["merged_net"]],
          algo = cluster_algo_used,
          case = "best",
          resolution = resolution,
          partition_type = partition_type,
          it = no_of_iterations,
          seed_offset = iteration_now
        )
      })
    )

    list(
      cluster_algo_used = cluster_algo_used,
      gene_which_cluster = gene_which_cluster
    )
  })
  cluster_algo_used <- cluster_run$cluster_algo_used
  gene_which_cluster <- cluster_run$gene_which_cluster

  if (base::ncol(gene_which_cluster) > 1) {
    # Align each replicate's community labels to the first one before voting.
    # Without this the vote compares arbitrary label permutations instead of
    # actual assignments, and sends the majority of genes to cluster 0 even when
    # the partitions agree: on a 726-gene test network cluster_louvain lost
    # 68% of genes this way, versus 20% once the labels are matched.
    gene_which_cluster <- base::apply(gene_which_cluster, 2, base::as.character)
    reference_partition <- gene_which_cluster[, 1]
    for (j in base::seq_len(base::ncol(gene_which_cluster))[-1]) {
      gene_which_cluster[, j] <- .hc_match_partition_labels(
        gene_which_cluster[, j],
        reference_partition
      )
    }

    gene_cluster_ident <- base::apply(gene_which_cluster, 1, function(x) {
      if (base::length(base::unique(x)) > max_cluster_count_per_gene) {
        0
      } else {
        base::names(base::which(base::table(x) == base::max(base::table(x))))[1]
      }
    })
  } else {
    gene_cluster_ident <- gene_which_cluster[, 1]
  }

  white_genes_clustercounts <- base::as.integer(
    base::length(base::grep(gene_cluster_ident, pattern = "\\b0\\b"))
  )

  message(
    white_genes_clustercounts,
    " genes were assigned to more than ",
    max_cluster_count_per_gene,
    " cluster(s). These genes are assigned to Cluster 0 (white) and will be left out of the network and further analyses."
  )

  cluster_data <- base::data.frame(
    genes = igraph::vertex_attr(hcobject[["integrated_output"]][["merged_net"]], "name"),
    clusters = base::paste0("Cluster ", gene_cluster_ident),
    stringsAsFactors = FALSE
  )

  dfk <- cluster_data %>%
    dplyr::count(clusters, genes) %>%
    dplyr::group_by(clusters) %>%
    dplyr::summarise(
      gene_no = base::sum(n),
      gene_n = base::paste0(genes, collapse = ",")
    ) %>%
    dplyr::mutate(
      cluster_included = base::ifelse(
        gene_no >= hcobject[["global_settings"]][["min_nodes_number_for_cluster"]],
        "yes",
        "no"
      ),
      color = "white"
    )

  color.cluster <- get_cluster_colours()
  plot_clusters <- ggplot2::ggplot(
    data = dfk[dfk$cluster_included == "yes" & dfk$clusters != "Cluster 0", ],
    ggplot2::aes(x = clusters)
  ) +
    ggplot2::geom_bar(ggplot2::aes(fill = clusters)) +
    ggplot2::scale_fill_manual(values = color.cluster)

  plot_clust <- ggplot2::ggplot_build(plot_clusters)
  dfk[dfk$cluster_included == "yes" & dfk$clusters != "Cluster 0", "color"] <- plot_clust$data[[1]][["fill"]]

  white_genes_clustersize <- base::as.integer(
    dfk %>%
      dplyr::filter(cluster_included == "no") %>%
      dplyr::summarise(n = base::sum(gene_no)) %>%
      purrr::map(1)
  )

  message(
    white_genes_clustersize,
    " genes were assigned to clusters with a smaller size than the defined minimal cluster size of ",
    hcobject[["global_settings"]][["min_nodes_number_for_cluster"]],
    " genes per cluster. These genes will also be assigned to Cluster 0 (white) and left out of the network and further analyses."
  )

  dfk_allinfo <- base::do.call(
    "rbind",
    base::lapply(
      base::seq_len(base::nrow(dfk)),
      gfc_mean_clustergene,
      cluster_df = dfk,
      gfc_dat = hcobject[["integrated_output"]][["GFC_all_layers"]]
    )
  )
  dfk_allinfo$vertexsize <- base::ifelse(dfk_allinfo$cluster_included == "yes", 3, 1)

  if (isTRUE(return_result)) {
    return(dfk_allinfo)
  }

  .hc_set_bridge_hcobject_slot(c("integrated_output", "cluster_calc", "cluster_information"), dfk_allinfo)
  .hc_set_bridge_hcobject_slot(c("integrated_output", "cluster_calc", "labelled_network"), NULL)
  .hc_set_bridge_hcobject_slot(c("integrated_output", "cluster_calc", "network_col_by_module"), NULL)
  .hc_set_bridge_hcobject_slot(c("satellite_outputs", "labelled_network"), NULL)
  .hc_set_bridge_hcobject_slot(c("satellite_outputs", "network_col_by_module"), NULL)
  invisible(NULL)
}

.hc_cluster_calculation_impl <- function(hc,
                                         cluster_algo = "cluster_leiden",
                                         no_of_iterations = 2,
                                         resolution = 0.1,
                                         partition_type = "RBConfigurationVertexPartition",
                                         max_cluster_count_per_gene = 1,
                                         return_result = FALSE) {
  out <- .hc_run_driver_capture(
    hc = hc,
    fun = .hc_cluster_calculation_driver,
    cluster_algo = cluster_algo,
    no_of_iterations = no_of_iterations,
    resolution = resolution,
    partition_type = partition_type,
    max_cluster_count_per_gene = max_cluster_count_per_gene,
    return_result = return_result
  )

  if (isTRUE(return_result) && !base::is.null(out$result)) {
    return(list(hc = out$hc, result = out$result))
  }

  out$hc
}
#
#' Cluster Calculation
#'
#' The function offers several community detection algorithms to identify dense regions in the co-expression network.
#'  These dense regions represent collections of highly co-expressed genes likely to form a functional group.
#' @param hc A `HCoCenaExperiment`.
#' @return Updated `HCoCenaExperiment`. If `return_result = TRUE` and the
#'   clustering backend returns a cluster table, that table is returned
#'   instead.
#' @param cluster_algo The clustering algorithm to be used. The choice is between "cluster_leiden" (default), "cluster_louvain", "cluster_label_prop", "cluster_fast_greedy",
#'  "cluster_infomap", "cluster_walktrap" and "auto" (in which case all are tested and the one with the highest modularity is chosen).
#' @param no_of_iterations Some of the algorithms are iterative (e.g. Leiden Algorithm). Set here, how many iterations should be performed.
#'  For information on which other algorithms are iterative, please refer to their documentation in the igraph or leidenbase package. Default is 2.
#' @param max_cluster_count_per_gene The maximum number of different clusters a
#'  gene is allowed to be associated with during the different iterations
#'  before it is marked as indecisive and removed. Default is 1.
#' @param resolution The cluster resolution if the cluster algorithm is set to "cluster_leiden". Default is 0.1. Higher values result in more clusters and vice versa.
#' @param partition_type Name of the partition type. Select from 'CPMVertexPartition', 'ModularityVertexPartition', 'RBConfigurationVertexPartition' and 'RBERVertexPartition'. Default is 'RBConfigurationVertexPartition'.
#' @param return_result Logical. If `TRUE`, return the cluster table instead of
#'  storing it in `hcobject`.
#' @examples
#' hc <- hc_example_data("clustered")
#' hc <- hc_cluster_calculation(hc, cluster_algo = "cluster_leiden")
#' @export
hc_cluster_calculation <- function(hc,
                                   cluster_algo = "cluster_leiden",
                                   no_of_iterations = 2,
                                   resolution = 0.1,
                                   partition_type = "RBConfigurationVertexPartition",
                                   max_cluster_count_per_gene = 1,
                                   return_result = FALSE) {
  out <- .hc_cluster_calculation_impl(
    hc = hc,
    cluster_algo = cluster_algo,
    no_of_iterations = no_of_iterations,
    resolution = resolution,
    partition_type = partition_type,
    max_cluster_count_per_gene = max_cluster_count_per_gene,
    return_result = return_result
  )

  if (base::is.list(out) &&
    "hc" %in% base::names(out) &&
    inherits(out[["hc"]], "HCoCenaExperiment")) {
    return(out[["result"]])
  }

  out
}

#' Merge modules based on module-heatmap similarity (S4 API)
#'
#' @noRd
.hc_merge_clusters_impl <- function(hc,
                                    k = "auto",
                                    save = TRUE,
                                    method = "complete",
                                    k_min = 2,
                                    k_max = NULL,
                                    auto_parsimony_penalty = 1e-04,
                                    verbose = TRUE) {
  out <- .hc_run_driver_capture(
    hc = hc,
    fun = .hc_merge_clusters_driver,
    k = k,
    save = save,
    method = method,
    k_min = k_min,
    k_max = k_max,
    auto_parsimony_penalty = auto_parsimony_penalty,
    verbose = verbose
  )

  if (!base::is.null(out$result)) {
    return(list(hc = out$hc, result = out$result))
  }

  out$hc
}
#
#' Merge similar modules
#'
#' Cuts the dendrogram of the module heatmap to merge modules with similar
#' expression patterns. Run with `save = FALSE` first to preview a given `k`;
#' `save = TRUE` overwrites the current module assignment and cannot be undone.
#'
#' @param hc A `HCoCenaExperiment`.
#' @param k Either an integer (the number of merged modules) or `"auto"`
#'   (default) to pick a data-driven value from the average silhouette score
#'   across candidate cuts.
#' @param save Logical. If `FALSE`, only preview the merge. If `TRUE`, replace
#'   the current module assignment. Default is `TRUE`.
#' @param method Agglomeration method for the dendrogram. Default `"complete"`.
#' @param k_min Smallest candidate `k` considered when `k = "auto"`. Default 2.
#' @param k_max Largest candidate `k` considered when `k = "auto"`. Defaults to
#'   `min(10, number of modules - 1)`.
#' @param auto_parsimony_penalty Small penalty applied in auto mode to prefer
#'   fewer modules when silhouette scores are nearly tied. Default `1e-04`.
#' @param verbose Logical. Print the selected `k` and selection diagnostics.
#' @return Updated `HCoCenaExperiment`.
#' @examples
#' hc <- hc_example_data("clustered")
#' hc <- hc_plot_cluster_heatmap(hc, file_name = FALSE)
#' hc <- hc_merge_clusters(hc, k = 2, save = FALSE)
#' @export
hc_merge_clusters <- function(hc,
                              k = "auto",
                              save = TRUE,
                              method = "complete",
                              k_min = 2,
                              k_max = NULL,
                              auto_parsimony_penalty = 1e-04,
                              verbose = TRUE) {
  out <- .hc_merge_clusters_impl(
    hc = hc,
    k = k,
    save = save,
    method = method,
    k_min = k_min,
    k_max = k_max,
    auto_parsimony_penalty = auto_parsimony_penalty,
    verbose = verbose
  )

  if (base::is.list(out) &&
    "hc" %in% base::names(out) &&
    inherits(out[["hc"]], "HCoCenaExperiment")) {
    return(out[["hc"]])
  }

  out
}

#' Split one or multiple modules into submodules (S4 API)
#'
#' Re-clusters genes within selected modules and replaces those modules by
#' submodules (with color shades of the parent module).
#'
#' @param hc A `HCoCenaExperiment`.
#' @param modules Character/numeric vector of modules to split. Accepts module
#'   labels (e.g. `"M3"`), module colors, or module indices.
#' @return Updated `HCoCenaExperiment`.
#' @param cluster_algo Clustering algorithm for within-module splitting.
#'  One of `"cluster_leiden"` (default), `"cluster_louvain"`,
#'  `"cluster_fast_greedy"`, `"cluster_infomap"`, `"cluster_walktrap"`,
#'  `"cluster_label_prop"` or `"auto"`.
#' @param no_of_iterations Number of Leiden iterations (used only for Leiden).
#' @param resolution Leiden resolution (used only for Leiden). Use either one
#'  positive value for all selected modules or one positive value per selected
#'  module in the same order as `modules`. Named vectors may use module labels
#'  or module colors.
#' @param resolution_grid Optional numeric vector of candidate resolutions to
#'  test before splitting. For each candidate, hCoCena reports how many
#'  submodules would be retained after size filtering.
#' @param resolution_test_only Logical; if `TRUE`, only run the resolution test
#'  (when `resolution_grid` is set) and do not apply any split.
#' @param partition_type Leiden partition type (used only for Leiden).
#' @param seed Random seed used for deterministic clustering.
#' @param drop_small_submodules Logical; if `TRUE` (default), submodules with
#'  fewer than `min_submodule_size` genes are dropped from the module
#'  annotation (their genes become unassigned/white in network visualizations).
#' @param min_submodule_size Optional minimum size for retained split
#'  submodules. If `NULL`, uses `global_settings$min_nodes_number_for_cluster`.
#' @param min_module_size Optional alias for `min_submodule_size`. If set, it
#'  takes precedence.
#' @param verbose Logical; print progress messages.
#' @examples
#' hc <- hc_example_data("clustered")
#' hc <- hc_plot_cluster_heatmap(hc, file_name = FALSE)
#' hc <- hc_split_modules(hc, modules = "M1", resolution = 1)
#' @export
hc_split_modules <- function(hc,
                             modules,
                             cluster_algo = "cluster_leiden",
                             no_of_iterations = 2,
                             resolution = 0.1,
                             resolution_grid = NULL,
                             resolution_test_only = FALSE,
                             partition_type = "RBConfigurationVertexPartition",
                             seed = 168575,
                             drop_small_submodules = TRUE,
                             min_submodule_size = NULL,
                             min_module_size = NULL,
                             verbose = TRUE) {
  .hc_run_driver(
    hc = hc,
    fun = .hc_split_modules_driver,
    modules = modules,
    cluster_algo = cluster_algo,
    no_of_iterations = no_of_iterations,
    resolution = resolution,
    resolution_grid = resolution_grid,
    resolution_test_only = resolution_test_only,
    partition_type = partition_type,
    seed = seed,
    drop_small_submodules = drop_small_submodules,
    min_submodule_size = min_submodule_size,
    min_module_size = min_module_size,
    verbose = verbose
  )
}

#' Undo module splitting (S4 API)
#'
#' Restores module structure from split history.
#'
#' @param hc A `HCoCenaExperiment`.
#' @return Updated `HCoCenaExperiment`.
#' @param which Either `"last"` (undo one split step) or `"all"` (restore the
#'  original pre-split cluster state).
#' @param verbose Logical; print progress messages.
#' @examples
#' hc <- hc_example_data("clustered")
#' hc <- hc_plot_cluster_heatmap(hc, file_name = FALSE)
#' hc <- hc_split_modules(hc, modules = "M1", resolution = 1)
#' hc <- hc_unsplit_modules(hc, which = "all")
#' @export
hc_unsplit_modules <- function(hc,
                               which = c("last", "all"),
                               verbose = TRUE) {
  which <- base::match.arg(which)
  .hc_run_driver(
    hc = hc,
    fun = .hc_unsplit_modules_driver,
    which = which,
    verbose = verbose
  )
}

# Internal helpers for S4 functional-enrichment panel redraw.
# @noRd
.hc_functional_enrichment_entries <- function(hc) {
  sat <- as.list(hc@satellite)
  enrich <- sat[["enrichments"]]
  if (is.null(enrich)) {
    sat_keys <- names(sat)
    if (!is.null(sat_keys) && any(grepl("^top_", sat_keys))) {
      enrich <- sat
    }
  }
  if (is.null(enrich)) {
    return(list())
  }
  enrich_names <- names(enrich)
  enrich <- as.list(enrich)
  if (is.null(names(enrich)) && !is.null(enrich_names) && length(enrich_names) == length(enrich)) {
    names(enrich) <- enrich_names
  }
  if (length(enrich) == 0) {
    return(list())
  }
  enrich
}

.hc_functional_enrichment_panel_keys <- function(enrich) {
  if (is.null(enrich) || length(enrich) == 0) {
    return(character(0))
  }

  panel_keys <- names(enrich)
  panel_keys <- panel_keys[grepl("^top_", panel_keys)]
  if (length(panel_keys) == 0) {
    return(character(0))
  }
  panel_keys <- unique(panel_keys)

  requested_order <- enrich[["panel_order"]]
  if (!is.null(requested_order)) {
    requested_order <- base::as.character(requested_order)
    requested_order <- requested_order[!base::is.na(requested_order) & base::nzchar(requested_order)]
    requested_order <- requested_order[requested_order %in% panel_keys]
    panel_keys <- unique(base::c(requested_order, panel_keys))
  }

  single_db_keys <- panel_keys[!panel_keys %in% c("top_all_dbs", "top_all_dbs_mixed")]
  combined_keys <- intersect(c("top_all_dbs", "top_all_dbs_mixed"), panel_keys)
  panel_keys <- unique(c(single_db_keys, combined_keys))

  panel_keys
}

.hc_functional_enrichment_has_draw_object <- function(entry) {
  entry <- tryCatch(as.list(entry), error = function(e) NULL)
  if (is.null(entry)) {
    return(FALSE)
  }
  for (nm in c("p", "hc_heatmap", "enrichment_plot")) {
    if (nm %in% names(entry) && !is.null(entry[[nm]])) {
      return(TRUE)
    }
  }
  FALSE
}

.hc_functional_enrichment_draw_object <- function(entry, heatmap_side = "left") {
  if (is.null(entry)) {
    return(NULL)
  }
  entry <- tryCatch(as.list(entry), error = function(e) NULL)
  if (is.null(entry)) {
    return(NULL)
  }

  clone_obj <- function(x) {
    .hc_safe_deep_clone(x, context = "stored enrichment panel object")
  }

  draw_obj <- NULL
  has_components <- (
    "hc_heatmap" %in% names(entry) &&
      "enrichment_plot" %in% names(entry)
  )
  if (isTRUE(has_components)) {
    hc_ht <- clone_obj(entry[["hc_heatmap"]])
    enr_ht <- clone_obj(entry[["enrichment_plot"]])
    if (inherits(hc_ht, "Heatmap") && inherits(enr_ht, "Heatmap")) {
      draw_obj <- if (identical(heatmap_side, "right")) {
        enr_ht + hc_ht
      } else {
        hc_ht + enr_ht
      }
    }
  }
  if (is.null(draw_obj)) {
    p <- clone_obj(entry[["p"]])
    if (inherits(p, "HeatmapList")) {
      draw_obj <- p
    }
  }
  draw_obj
}

.hc_functional_enrichment_panel_title <- function(panel_key) {
  if (is.null(panel_key) || length(panel_key) != 1 || is.na(panel_key)) {
    return(NULL)
  }
  panel_key <- as.character(panel_key[[1]])
  if (identical(panel_key, "top_all_dbs")) {
    return("Combined enrichment (all DBs)")
  }
  if (identical(panel_key, "top_all_dbs_mixed")) {
    return("Combined enrichment (all DBs, mixed)")
  }
  db_name <- sub("^top_", "", panel_key)
  if (!nzchar(db_name) || identical(db_name, panel_key)) {
    return(NULL)
  }
  db_name <- gsub("_", " ", db_name, fixed = TRUE)
  paste0(db_name, " enrichment")
}

.hc_functional_enrichment_panel_target <- function(panel_key) {
  if (is.null(panel_key) || length(panel_key) != 1 || is.na(panel_key)) {
    return(NULL)
  }
  panel_key <- as.character(panel_key[[1]])
  if (identical(panel_key, "top_all_dbs")) {
    return("enrichment_all")
  }
  if (identical(panel_key, "top_all_dbs_mixed")) {
    return("enrichment_all_mixed")
  }
  "enrichment"
}

.hc_draw_functional_enrichment_panel_title <- function(panel_key, fontsize = 16) {
  panel_title <- .hc_functional_enrichment_panel_title(panel_key)
  panel_target <- .hc_functional_enrichment_panel_target(panel_key)
  if (is.null(panel_title) || !nzchar(panel_title) || is.null(panel_target)) {
    return(invisible(NULL))
  }

  panel_key_chr <- as.character(panel_key[[1]])
  title_drawn <- FALSE
  is_combined_panel <- panel_key_chr %in% c("top_all_dbs", "top_all_dbs_mixed")
  slice_candidates <- if (identical(panel_key_chr, "top_all_dbs")) c(2L, 1L) else 1L
  title_offset_mm <- if (isTRUE(is_combined_panel)) {
    max(8, 0.50 * as.numeric(fontsize))
  } else {
    max(4, 0.35 * as.numeric(fontsize))
  }
  components <- try(ComplexHeatmap::list_components(), silent = TRUE)
  if (inherits(components, "try-error")) {
    components <- character(0)
  }

  if (isTRUE(is_combined_panel) && length(components) > 0) {
    body_pattern <- paste0("^", panel_target, "_heatmap_body_[0-9]+_[0-9]+$")
    body_components <- components[grepl(body_pattern, components)]
    if (length(body_components) >= 1) {
      draw_res <- try(
        {
          left_edges <- numeric(0)
          right_edges <- numeric(0)
          top_edges <- numeric(0)

          for (vp_name in body_components) {
            grid::seekViewport(vp_name)
            loc_left <- grid::deviceLoc(
              x = grid::unit(0, "npc"),
              y = grid::unit(1, "npc")
            )
            loc_right <- grid::deviceLoc(
              x = grid::unit(1, "npc"),
              y = grid::unit(1, "npc")
            )
            grid::upViewport(0)
            left_edges <- c(left_edges, grid::convertX(loc_left[["x"]], "inches", valueOnly = TRUE))
            right_edges <- c(right_edges, grid::convertX(loc_right[["x"]], "inches", valueOnly = TRUE))
            top_edges <- c(top_edges, grid::convertY(loc_left[["y"]], "inches", valueOnly = TRUE))
          }

          loc_x <- grid::unit((min(left_edges) + max(right_edges)) / 2, "inches")
          loc_y <- grid::unit(max(top_edges), "inches")
          grid::grid.text(
            label = panel_title,
            x = loc_x,
            y = loc_y + grid::unit(title_offset_mm, "mm"),
            just = c("center", "bottom"),
            gp = grid::gpar(fontsize = fontsize, fontface = "bold")
          )
        },
        silent = TRUE
      )
      try(grid::upViewport(0), silent = TRUE)
      if (!inherits(draw_res, "try-error")) {
        return(invisible(NULL))
      }
    }
  }

  for (slice_idx in slice_candidates) {
    body_vp_name <- paste0(panel_target, "_heatmap_body_", slice_idx, "_1")
    if (length(components) > 0 && !(body_vp_name %in% components)) {
      next
    }
    draw_res <- try(
      {
        grid::seekViewport(body_vp_name)
        loc <- grid::deviceLoc(
          x = grid::unit(0.5, "npc"),
          y = grid::unit(1, "npc")
        )
        grid::upViewport(0)
        loc_x <- loc[["x"]]
        loc_y <- loc[["y"]]
        grid::grid.text(
          label = panel_title,
          x = loc_x,
          y = loc_y + grid::unit(title_offset_mm, "mm"),
          just = c("center", "bottom"),
          gp = grid::gpar(fontsize = fontsize, fontface = "bold")
        )
      },
      silent = TRUE
    )
    try(grid::upViewport(0), silent = TRUE)
    if (!inherits(draw_res, "try-error")) {
      title_drawn <- TRUE
      break
    }
  }
  if (!title_drawn) {
    return(invisible(NULL))
  }
  invisible(NULL)
}

.hc_draw_functional_enrichment_panels <- function(hc, heatmap_side = "left", record_history = FALSE) {
  enrich <- .hc_functional_enrichment_entries(hc)
  panel_keys <- .hc_functional_enrichment_panel_keys(enrich)
  if (length(panel_keys) == 0) {
    return(invisible(NULL))
  }

  if (isTRUE(record_history) && grDevices::dev.cur() > 1) {
    try(grDevices::dev.control(displaylist = "enable"), silent = TRUE)
  }

  draw_errors <- character(0)
  for (k in panel_keys) {
    draw_obj <- .hc_functional_enrichment_draw_object(enrich[[k]], heatmap_side = heatmap_side)
    if (is.null(draw_obj)) {
      next
    }

    draw_res <- try(
      ComplexHeatmap::draw(
        draw_obj,
        newpage = TRUE,
        merge_legends = TRUE,
        show_annotation_legend = TRUE,
        show_heatmap_legend = TRUE
      ),
      silent = TRUE
    )
    if (inherits(draw_res, "try-error")) {
      draw_errors <- c(draw_errors, k)
      next
    }
    try(.hc_draw_functional_enrichment_panel_title(k), silent = TRUE)
    if (isTRUE(record_history) && grDevices::dev.cur() > 1) {
      try(grDevices::recordPlot(), silent = TRUE)
    }
  }

  if (length(draw_errors) > 0) {
    warning(
      "Could not redraw enrichment panel(s): ",
      paste(unique(draw_errors), collapse = ", "),
      call. = FALSE
    )
  }
  invisible(NULL)
}

.hc_emit_functional_enrichment_knitr_panels <- function(hc,
                                                        heatmap_side = "left",
                                                        panel_keys = NULL,
                                                        res = 200) {
  if (!isTRUE(getOption("knitr.in.progress")) &&
    !isTRUE(getOption("rstudio.notebook.executing"))) {
    return(invisible(NULL))
  }
  if (requireNamespace("knitr", quietly = TRUE)) {
    old_fig_keep <- try(knitr::opts_current$get("fig.keep"), silent = TRUE)
    if (!inherits(old_fig_keep, "try-error")) {
      on.exit(
        try(knitr::opts_current$set(fig.keep = old_fig_keep), silent = TRUE),
        add = TRUE
      )
    }
    try(knitr::opts_current$set(fig.keep = "all"), silent = TRUE)
  }

  enrich <- .hc_functional_enrichment_entries(hc)
  available_panel_keys <- .hc_functional_enrichment_panel_keys(enrich)
  combined_panel_keys <- intersect(c("top_all_dbs", "top_all_dbs_mixed"), available_panel_keys)
  if (is.null(panel_keys)) {
    panel_keys <- available_panel_keys
  } else {
    panel_keys <- base::as.character(panel_keys)
    panel_keys <- panel_keys[!base::is.na(panel_keys) & base::nzchar(panel_keys)]
    panel_keys <- panel_keys[panel_keys %in% available_panel_keys]
  }
  panel_keys <- unique(panel_keys)
  panel_keys <- c(
    panel_keys[!panel_keys %in% combined_panel_keys],
    intersect(combined_panel_keys, panel_keys),
    setdiff(panel_keys[panel_keys %in% combined_panel_keys], combined_panel_keys)
  )
  if (length(panel_keys) == 0) {
    return(invisible(NULL))
  }
  panel_keys <- c(
    panel_keys[!panel_keys %in% c("top_all_dbs", "top_all_dbs_mixed")],
    intersect(c("top_all_dbs", "top_all_dbs_mixed"), panel_keys)
  )

  tmp_plot_dir <- tempfile(pattern = "hc_functional_enrichment_panels_")
  dir.create(tmp_plot_dir, recursive = TRUE, showWarnings = FALSE)
  png_paths <- character(0)
  failed_panels <- character(0)

  for (idx in base::seq_along(panel_keys)) {
    k <- panel_keys[[idx]]
    draw_obj <- .hc_functional_enrichment_draw_object(enrich[[k]], heatmap_side = heatmap_side)
    if (is.null(draw_obj)) {
      next
    }

    is_combined_panel <- k %in% c("top_all_dbs", "top_all_dbs_mixed")
    panel_width_in <- if (is_combined_panel) 18 else 14
    panel_height_in <- if (is_combined_panel) 10 else 9
    safe_key <- gsub("[^A-Za-z0-9._-]+", "_", k)
    png_file <- base::file.path(tmp_plot_dir, sprintf("%03d_%s.png", idx, safe_key))
    opened <- FALSE
    open_res <- try(
      {
        grDevices::png(
          filename = png_file,
          width = panel_width_in,
          height = panel_height_in,
          units = "in",
          res = res
        )
        opened <- TRUE
      },
      silent = TRUE
    )
    if (!isTRUE(opened) || inherits(open_res, "try-error")) {
      failed_panels <- c(failed_panels, k)
      next
    }

    draw_ok <- FALSE
    draw_res <- try(
      {
        ComplexHeatmap::draw(
          draw_obj,
          newpage = TRUE,
          merge_legends = TRUE,
          show_annotation_legend = TRUE,
          show_heatmap_legend = TRUE
        )
        .hc_draw_functional_enrichment_panel_title(k)
        draw_ok <- TRUE
      },
      silent = TRUE
    )
    try(grDevices::dev.off(), silent = TRUE)

    if (inherits(draw_res, "try-error") || !isTRUE(draw_ok) || !file.exists(png_file)) {
      failed_panels <- c(failed_panels, k)
      if (file.exists(png_file)) {
        try(unlink(png_file, force = TRUE), silent = TRUE)
      }
      next
    }
    png_paths <- c(png_paths, png_file)
  }

  rendered_paths <- NULL
  if (length(png_paths) > 0) {
    if (requireNamespace("knitr", quietly = TRUE)) {
      rendered_paths <- knitr::include_graphics(png_paths)
    } else {
      warning(
        "Package `knitr` is required for notebook panel rendering.",
        call. = FALSE
      )
    }
  }

  if (length(failed_panels) > 0) {
    warning(
      "Could not render enrichment panel(s) for knitr output: ",
      paste(unique(failed_panels), collapse = ", "),
      call. = FALSE
    )
  }

  if (!base::is.null(rendered_paths)) {
    return(rendered_paths)
  }

  invisible(NULL)
}

#' Draw saved functional-enrichment panels from the current object
#'
#' This is a convenience plotting helper for notebook/interactive usage. It
#' draws the enrichment panels already stored in `hc@satellite$enrichments`
#' after `hc_functional_enrichment()`, without recomputing enrichment.
#'
#' @param hc A `HCoCenaExperiment`.
#' @param panels Optional character vector of panel names to draw.
#'   Default (`NULL`) draws all available panels in standard order.
#'   Accepted values include `"top_Hallmark"`, `"top_Kegg"`,
#'   `"top_all_dbs"`, `"top_all_dbs_mixed"` (or without `"top_"` prefix).
#' @param heatmap_side One of `"left"` (default) or `"right"`.
#' @return Invisibly returns `hc`.
#' @examples
#' hc <- hc_example_data("clustered")
#' gmt <- tempfile(fileext = ".gmt")
#' writeLines(
#'   c(
#'     paste(c("T cells", "-", "CD3D", "CD3E", "CD3G", "CD2", "CD28",
#'             "LCK", "ZAP70", "IL7R", "CD7", "TRAC"), collapse = "\t"),
#'     paste(c("B cells", "-", "CD19", "MS4A1", "CD79A", "CD79B", "BLNK",
#'             "PAX5", "CR2", "FCRL1", "TNFRSF13B", "VPREB3"), collapse = "\t"),
#'     paste(c("Monocytes", "-", "CD14", "LYZ", "FCN1", "VCAN", "S100A8",
#'             "S100A9", "CSF1R", "ITGAM", "CD68", "FCGR3A"), collapse = "\t"),
#'     paste(c("NK cells", "-", "NKG7", "GNLY", "KLRD1", "KLRF1", "PRF1",
#'             "GZMB", "NCR1", "KLRC1", "FGFBP2", "SPON2"), collapse = "\t")
#'   ),
#'   gmt
#' )
#' hc <- hc_functional_enrichment(
#'   hc,
#'   gene_sets = character(0),
#'   custom_gmt_files = c(CellTypes = gmt)
#' )
#' p <- hc_plot_enrichment_panels(hc)
#' @export
hc_plot_enrichment_panels <- function(hc,
                                      panels = NULL,
                                      heatmap_side = "left") {
  if (!inherits(hc, "HCoCenaExperiment")) {
    stop("`hc` must be a `HCoCenaExperiment`.")
  }

  heatmap_side <- base::match.arg(heatmap_side, choices = c("left", "right"))

  enrich <- .hc_functional_enrichment_entries(hc)
  if (length(enrich) == 0) {
    stop(
      "No enrichment panels found in `hc@satellite$enrichments`. ",
      "Run `hc_functional_enrichment()` first."
    )
  }

  available_panel_keys <- .hc_functional_enrichment_panel_keys(enrich)
  if (length(available_panel_keys) == 0) {
    stop(
      "No drawable enrichment panels found. ",
      "Expected entries like `top_Hallmark`, `top_all_dbs`, ... in ",
      "`hc@satellite$enrichments`."
    )
  }

  combined_panel_keys <- intersect(c("top_all_dbs", "top_all_dbs_mixed"), available_panel_keys)
  single_panel_keys <- available_panel_keys[!available_panel_keys %in% combined_panel_keys]
  panel_keys <- c(single_panel_keys, combined_panel_keys)
  if (!is.null(panels)) {
    panels <- base::as.character(panels)
    panels <- panels[!base::is.na(panels) & base::nzchar(panels)]
    panels <- ifelse(grepl("^top_", panels), panels, paste0("top_", panels))
    missing_panels <- panels[!panels %in% available_panel_keys]
    if (length(missing_panels) > 0) {
      warning(
        "Skipping unknown enrichment panel(s): ",
        paste(unique(missing_panels), collapse = ", "),
        call. = FALSE
      )
    }
    requested <- panels[panels %in% available_panel_keys]
    if (length(requested) == 0) {
      stop(
        "None of the requested panels are available. Available panels: ",
        paste(available_panel_keys, collapse = ", ")
      )
    }

    panel_keys <- unique(requested)
  }

  panel_keys <- unique(panel_keys)
  drawable_mask <- base::vapply(
    panel_keys,
    function(k) .hc_functional_enrichment_has_draw_object(enrich[[k]]),
    FUN.VALUE = base::logical(1)
  )
  if (!base::any(drawable_mask)) {
    stop(
      "Selected enrichment panels were not stored as drawable objects in `hc`. ",
      "This usually happens in memory-saving mode for multi-database enrichment runs. ",
      "Rerun `hc_functional_enrichment(..., store_panel_objects = 'always')` if you need redraw from the object."
    )
  }
  message(
    "Drawing enrichment panels in order: ",
    paste(panel_keys, collapse = ", ")
  )

  in_notebook <- isTRUE(getOption("knitr.in.progress")) ||
    isTRUE(getOption("rstudio.notebook.executing"))
  if (isTRUE(in_notebook)) {
    for (k in panel_keys) {
      .hc_emit_functional_enrichment_knitr_panels(
        hc = hc,
        heatmap_side = heatmap_side,
        panel_keys = k
      )
    }
    return(invisible(hc))
  }

  draw_errors <- character(0)
  for (k in panel_keys) {
    draw_obj <- .hc_functional_enrichment_draw_object(enrich[[k]], heatmap_side = heatmap_side)
    if (is.null(draw_obj)) {
      draw_errors <- c(draw_errors, k)
      next
    }

    draw_res <- try(
      {
        grid::grid.newpage()
        ComplexHeatmap::draw(
          draw_obj,
          newpage = FALSE,
          merge_legends = TRUE,
          show_annotation_legend = TRUE,
          show_heatmap_legend = TRUE
        )
      },
      silent = TRUE
    )
    if (inherits(draw_res, "try-error")) {
      draw_errors <- c(draw_errors, k)
      next
    }
    try(.hc_draw_functional_enrichment_panel_title(k), silent = TRUE)
  }

  if (length(draw_errors) > 0) {
    warning(
      "Could not draw enrichment panel(s): ",
      paste(unique(draw_errors), collapse = ", "),
      call. = FALSE
    )
  }

  invisible(hc)
}

#' @noRd
# Internal implementation shared by the S4 entry point and the driver.
.hc_functional_enrichment_impl <- function(hc, ...) {
  out <- .hc_run_driver_capture(
    hc = hc,
    fun = ".hc_functional_enrichment_driver",
    ...
  )
  list(hc = out$hc, result = out$result)
}
#
#' Functional enrichment (S4 API)
#'
#' @param hc A `HCoCenaExperiment`.
#' @param ... Further arguments passed through to the underlying
#'   implementation. Every supported argument is documented on this page;
#'   `...` only carries the tail of the argument list.
#' @return Updated `HCoCenaExperiment`.
#' @rdname hc_functional_enrichment
#' @section Significance columns:
#' The result tables carry `pvalue` (raw), `p_adjusted` (corrected with the
#' method given in `padj`, and what the `qval` filter and the ranking use),
#' `padj_method` (which correction that was) and `q_storey` (clusterProfiler's
#' Storey q-value, `NA` where pi0 could not be fitted - which is common for
#' small gene-set collections). `qvalue` is kept as an alias of `p_adjusted`
#' so existing scripts and exports keep working.
#' @param universe Background gene set for the hypergeometric test. `"all_genes"`
#'  (default) uses every gene measured in any layer, asking whether a module is
#'  enriched relative to the transcriptome. `"network"` uses only the genes that
#'  entered the integrated network -- the genes that could have landed in a
#'  module at all -- asking whether a module is enriched relative to the other
#'  modules. `"network"` removes the bias that the top-variance selection and
#'  the correlation cutoff introduce, and typically returns fewer but more
#'  module-specific terms.
#' @param gene_sets A vector. The names of databases enrichment should be performed for. Choose one or multiple of "Go", "Kegg", "Hallmark", and/or "Reactome".
#'  Available databases depend on supplement files previously set
#'  Default is "Hallmark".
#' @param custom_gmt_files Optional custom GMT file(s) to include directly in
#'  enrichment. Accepts a character vector or named list of file paths.
#'  Unnamed entries are auto-labeled (`CustomEnrichment1`, ...). Paths can be
#'  absolute/relative or file names inside `dir_reference_files`.
#' @param top Integer. The number of most strongly enriched terms to return per cluster. Default is 5.
#' @param clusters Either "all" (default) or a vector of clusters as strings. Defines for which clusters to perform the enrichment.
#' @param padj Method to use for multiple testing correction. Can be one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".  Default is "BH" (Benjamini-Hochberg).
#'  The reported `qvalue` column holds the p-values adjusted with this method,
#'  and `qval` is applied to it. (Up to and including 0.99.7 the `qvalue` column
#'  carried `clusterProfiler`'s Storey q-value, which is computed independently
#'  of `padj`, so this argument did not affect which terms were called
#'  significant. Result tables from earlier versions can therefore differ.)
#' @param qval Upper threshold for the adjusted p-value. Default is 0.05.
#' @param consistent_terms Logical. If `TRUE` (default), use the union of the
#'  top `top` enriched terms across the selected modules for each database and
#'  also show those terms in other modules whenever they are significantly
#'  enriched there. If `FALSE`, show the top enriched terms per module
#'  separately.
#' @param heatmap_side Position of the hCoCena heatmap in the combined output.
#'  Choose one of "left" (default) or "right".
#' @param heatmap_cluster_rows A Boolean whether or not to cluster rows in the hCoCena heatmap.
#' @param cluster_columns A Boolean whether or not to cluster columns in the hCoCena heatmap.
#' @param heatmap_cluster_columns Legacy alias for `cluster_columns`.
#' @param heatmap_show_row_dend A Boolean whether to show the row dendrogram when `heatmap_cluster_rows = TRUE`.
#' @param heatmap_show_column_dend A Boolean whether to show the column dendrogram when `cluster_columns = TRUE`.
#' @param col_order Optional character vector overriding the hCoCena
#'   heatmap column order for this enrichment plot only. If `NULL` (default),
#'   the column order from the main module heatmap is reused when available.
#' @param heatmap_col_order Legacy alias for `col_order`.
#' @param heatmap_order Optional character vector specifying module order in the hCoCena heatmap.
#'  Entries can be module colors (e.g. "turquoise") or module labels from the main heatmap
#'  (e.g. "M1", "M2", ...). Modules not listed are appended afterwards.
#' @param heatmap_module_label_mode Controls labels inside module color boxes of the hCoCena heatmap.
#'  One of "same" or "none". "same" reuses the prefix from
#'  `hc_plot_cluster_heatmap()`
#'  and reindexes modules consecutively (M1, M2, ...) after final module ordering.
#' @param heatmap_show_gene_counts A Boolean whether or not to show gene counts per module
#'  in the right annotation of the hCoCena heatmap. Default is FALSE.
#' @param heatmap_column_label_fontsize Optional numeric font size for hCoCena
#'  heatmap column labels in enrichment plots. If NULL (default), uses automatic sizing.
#' @param heatmap_module_label_fontsize Optional numeric font size for module labels
#'  inside module color boxes (`M1`, `M2`, ...). If NULL (default), uses automatic sizing.
#' @param legend_fontsize Optional numeric base font size for enrichment-related legends
#'  (GFC legend and enrichment significance legend). If NULL (default), uses automatic sizing.
#' @param enrichment_label_fontsize Optional numeric font size for enrichment term labels.
#'  If NULL (default), uses automatic sizing.
#' @param enrichment_db_header_fontsize Optional numeric font size for database
#'  headers in combined all-DB enrichment plots (e.g. "Go", "Kegg", "Hallmark").
#'  If NULL (default), uses automatic sizing.
#' @param enrichment_label_wrap Logical. If TRUE, wraps enrichment term labels using
#'  `enrichment_label_wrap_width`. Default is FALSE.
#' @param enrichment_label_wrap_width Integer wrap width used when
#'  `enrichment_label_wrap = TRUE`. Default is 30.
#' @param gfc_scale_limits Optional numeric vector controlling the module-heatmap
#'  color scale limits used in enrichment plots. Provide one positive number
#'  (`x` -> `c(-x, x)`) or two numbers (`c(min, max)`). If NULL, uses stored
#'  limits from the latest main module heatmap when available, otherwise
#'  falls back to `c(-range_GFC, range_GFC)`.
#' @param pdf_width Optional numeric width (inches) for enrichment PDFs.
#'  If NULL (default), width is auto-estimated from content.
#' @param pdf_height Optional numeric height (inches) for enrichment PDFs.
#'  If NULL (default), height is auto-estimated from content.
#' @param pdf_pointsize Numeric base pointsize used for PDF devices.
#'  Default is 11.
#' @param store_panel_objects One of `"auto"` (default), `"always"`, or `"never"`.
#'  Controls whether heavy heatmap/panel objects are stored inside `hc` for later
#'  redraw with `hc_plot_enrichment_panels()`. `"auto"` stores them only for a
#'  single selected database; multi-database runs keep only tables to save memory.
#' @examples
#' hc <- hc_example_data("clustered")
#' gmt <- tempfile(fileext = ".gmt")
#' writeLines(
#'   c(
#'     paste(c("T cells", "-", "CD3D", "CD3E", "CD3G", "CD2", "CD28",
#'             "LCK", "ZAP70", "IL7R", "CD7", "TRAC"), collapse = "\t"),
#'     paste(c("B cells", "-", "CD19", "MS4A1", "CD79A", "CD79B", "BLNK",
#'             "PAX5", "CR2", "FCRL1", "TNFRSF13B", "VPREB3"), collapse = "\t"),
#'     paste(c("Monocytes", "-", "CD14", "LYZ", "FCN1", "VCAN", "S100A8",
#'             "S100A9", "CSF1R", "ITGAM", "CD68", "FCGR3A"), collapse = "\t"),
#'     paste(c("NK cells", "-", "NKG7", "GNLY", "KLRD1", "KLRF1", "PRF1",
#'             "GZMB", "NCR1", "KLRC1", "FGFBP2", "SPON2"), collapse = "\t")
#'   ),
#'   gmt
#' )
#' hc <- hc_functional_enrichment(
#'   hc,
#'   gene_sets = character(0),
#'   custom_gmt_files = c(CellTypes = gmt)
#' )
#' @export
hc_functional_enrichment <- function(hc,
                                     gene_sets = c("Go", "Kegg", "Hallmark", "Reactome"),
                                     custom_gmt_files = NULL,
                                     top = 5,
                                     clusters = c("all"),
                                     padj = "BH",
                                     qval = 0.05,
                                     universe = c("all_genes", "network"),
                                     consistent_terms = TRUE,
                                     heatmap_side = "left",
                                     heatmap_cluster_rows = FALSE,
                                     cluster_columns = FALSE,
                                     heatmap_cluster_columns = NULL,
                                     heatmap_show_row_dend = FALSE,
                                     heatmap_show_column_dend = FALSE,
                                     col_order = NULL,
                                     heatmap_col_order = NULL,
                                     heatmap_order = NULL,
                                     heatmap_module_label_mode = "same",
                                     heatmap_show_gene_counts = FALSE,
                                     heatmap_column_label_fontsize = NULL,
                                     heatmap_module_label_fontsize = NULL,
                                     legend_fontsize = NULL,
                                     enrichment_label_fontsize = NULL,
                                     enrichment_db_header_fontsize = NULL,
                                     enrichment_label_wrap = FALSE,
                                     enrichment_label_wrap_width = 30,
                                     gfc_scale_limits = NULL,
                                     store_panel_objects = c("auto", "always", "never"),
                                     pdf_width = NULL,
                                     pdf_height = NULL,
                                     pdf_pointsize = 11,
                                     ...) {
  .hc_functional_enrichment_impl(
    hc = hc,
    gene_sets = gene_sets,
    custom_gmt_files = custom_gmt_files,
    top = top,
    clusters = clusters,
    padj = padj,
    qval = qval,
    universe = universe,
    consistent_terms = consistent_terms,
    heatmap_side = heatmap_side,
    heatmap_cluster_rows = heatmap_cluster_rows,
    cluster_columns = cluster_columns,
    heatmap_cluster_columns = heatmap_cluster_columns,
    heatmap_show_row_dend = heatmap_show_row_dend,
    heatmap_show_column_dend = heatmap_show_column_dend,
    col_order = col_order,
    heatmap_col_order = heatmap_col_order,
    heatmap_order = heatmap_order,
    heatmap_module_label_mode = heatmap_module_label_mode,
    heatmap_show_gene_counts = heatmap_show_gene_counts,
    heatmap_column_label_fontsize = heatmap_column_label_fontsize,
    heatmap_module_label_fontsize = heatmap_module_label_fontsize,
    legend_fontsize = legend_fontsize,
    enrichment_label_fontsize = enrichment_label_fontsize,
    enrichment_db_header_fontsize = enrichment_db_header_fontsize,
    enrichment_label_wrap = enrichment_label_wrap,
    enrichment_label_wrap_width = enrichment_label_wrap_width,
    gfc_scale_limits = gfc_scale_limits,
    store_panel_objects = store_panel_objects,
    pdf_width = pdf_width,
    pdf_height = pdf_height,
    pdf_pointsize = pdf_pointsize,
    ...
  )[["hc"]]
}

#' @noRd
# Internal implementation shared by the S4 entry point and the driver.
.hc_upstream_inference_impl <- function(hc, ...) {
  out <- .hc_run_driver_capture(
    hc = hc,
    fun = ".hc_upstream_inference_driver",
    ...
  )
  list(hc = out$hc, result = out$result)
}
#
#' Upstream regulator/pathway inference (S4 API)
#'
#' @param hc A `HCoCenaExperiment`.
#' @return Updated `HCoCenaExperiment`.
#' @rdname hc_upstream_inference
#' @param resources Character vector of upstream resources to use.
#'   Allowed values are `"TF"` and `"Pathway"`. Default is both.
#' @param top Integer. Number of top significant regulators/pathways to keep
#'   per module and resource in the selected summary. Default is 5.
#' @param clusters Either `"all"` (default) or a character vector of module
#'   colors to process.
#' @param padj Multiple-testing correction method passed to
#'   [stats::p.adjust()]. Default is `"BH"`.
#' @param qval Adjusted p-value threshold for significance. Default is 0.05.
#' @param tf_confidence Character vector of DoRothEA confidence levels to keep.
#'   Default is `c("A", "B", "C")`.
#' @param minsize Minimum number of targets required per source in
#'   `decoupleR::run_ulm()`. Default is 5.
#' @param method Inference method name used via `decoupleR::run_<method>`.
#'   Currently only `"ulm"` is supported. Default is `"ulm"`.
#' @param activity_input Character scalar selecting the matrix used for
#'   decoupleR activity inference:
#'   `"gfc"` (default) uses `integrated_output$GFC_all_layers`,
#'   `"fc"` uses user-defined pairwise fold-changes from `fc_comparisons`,
#'   `"expression"` uses layer-wise mean expression values (anti-log transformed
#'   when `data_in_log = TRUE`) across samples.
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
#' @examples
#' if (requireNamespace("decoupleR", quietly = TRUE)) {
#'   hc <- hc_example_data("clustered")
#'   gmt <- system.file("extdata", "toy_celltype_markers.gmt",
#'                      package = "hcocena")
#'   # Pathway activity from a local GMT; the DoRothEA TF priors
#'   # (`resources = "TF"`) need the `dorothea` package.
#'   hc <- hc_upstream_inference(
#'     hc,
#'     resources = "Pathway",
#'     custom_pathway_gmt = c(CellTypes = gmt),
#'     plot = FALSE,
#'     save_pdf = FALSE
#'   )
#'   head(hc_satellite(hc, "upstream_inference")$selected_upstream_all)
#' }
#' @export
hc_upstream_inference <- function(hc,
                                  resources = c("TF", "Pathway"),
                                  top = 5,
                                  clusters = c("all"),
                                  padj = "BH",
                                  qval = 0.05,
                                  tf_confidence = c("A", "B", "C"),
                                  minsize = 5,
                                  method = "ulm",
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
  .hc_upstream_inference_impl(
    hc = hc,
    resources = resources,
    top = top,
    clusters = clusters,
    padj = padj,
    qval = qval,
    tf_confidence = tf_confidence,
    minsize = minsize,
    method = method,
    activity_input = activity_input,
    fc_comparisons = fc_comparisons,
    custom_pathway_gmt = custom_pathway_gmt,
    heatmap_side = heatmap_side,
    cluster_columns = cluster_columns,
    heatmap_cluster_columns = heatmap_cluster_columns,
    col_order = col_order,
    heatmap_col_order = heatmap_col_order,
    gfc_scale_limits = gfc_scale_limits,
    plot = plot,
    save_pdf = save_pdf,
    pdf_width = pdf_width,
    pdf_height = pdf_height,
    pdf_pointsize = pdf_pointsize,
    plot_per_comparison = plot_per_comparison,
    consistent_terms = consistent_terms,
    overall_plot_scale = overall_plot_scale
  )[["hc"]]
}

#' @noRd
# Internal implementation shared by the S4 entry point and the driver.
.hc_celltype_annotation_impl <- function(hc, ...) {
  out <- .hc_run_driver_capture(
    hc = hc,
    fun = ".hc_celltype_annotation_driver",
    ...
  )
  list(hc = out$hc, result = out$result)
}
#
#' Module cell-type annotation from Enrichr (S4 API)
#'
#' @param hc A `HCoCenaExperiment`.
#' @param ... Further arguments passed through to the underlying
#'   implementation. Every supported argument is documented on this page;
#'   `...` only carries the tail of the argument list.
#' @return Updated `HCoCenaExperiment`.
#' @rdname hc_celltype_annotation
#' @param databases Character vector of Enrichr library names.
#' @param custom_gmt_files Optional custom GMT file(s) to include in the
#'   cell-type annotation. Accepts a character vector or named list of file
#'   paths. Paths can be absolute/relative or file names inside
#'   `dir_reference_files`.
#' @param clusters Either `"all"` (default) or a character vector of module
#'   IDs/colors.
#' @param mode Either `"coarse"` (broad classes) or `"fine"` (specific terms).
#' @param top Number of selected categories per module.
#' @param qval Maximum adjusted p-value (`qvalue`) for significant terms.
#' @param padj Multiple-testing correction method for [clusterProfiler::enricher()].
#' @param min_term_genes Minimum number of genes required per Enrichr term.
#' @param min_gs_size Minimum gene-set size used in [clusterProfiler::enricher()]
#'   (`minGSSize`).
#' @param max_gs_size Maximum gene-set size used in [clusterProfiler::enricher()]
#'   (`maxGSSize`).
#' @param annotation_slot Backward-compatibility option when only one database
#'   is used. With multiple databases, one slot per DB is always written
#'   (`enriched_per_cluster_<db>`), and previous DB slots are reset on each run.
#' @param slot_suffix Optional character suffix appended to generated annotation
#'   slots (for example `"decoupler"` -> `enriched_per_cluster_<db>_decoupler`).
#'   Useful to keep multiple annotation runs side by side.
#' @param clear_previous_slots Logical. If `TRUE` (default), previous
#'   `enriched_per_cluster*` slots are removed before writing new results.
#'   Set to `FALSE` to keep existing slots.
#' @param coarse_map Optional named character vector with regex rules for
#'   `mode = "coarse"`. Names are output class labels.
#' @param coarse_include_other Logical. Keep unmatched terms as `"Other"` in
#'   `mode = "coarse"`.
#' @param refresh_db Logical. If `TRUE`, re-download Enrichr metadata/libraries.
#' @param export_excel Logical. If `TRUE`, write summary tables to Excel.
#' @param excel_file Excel file name in the configured output folder.
#' @param plot_heatmap Logical. If `TRUE`, run [hc_plot_cluster_heatmap()] after
#'   updating annotation slots.
#' @param heatmap_file_name File name used when `plot_heatmap = TRUE`.
#' @examples
#' hc <- hc_example_data("clustered")
#' gmt <- system.file("extdata", "toy_celltype_markers.gmt", package = "hcocena")
#' # A local marker GMT instead of the default Enrichr libraries, so that the
#' # example runs offline:
#' hc <- hc_celltype_annotation(
#'   hc,
#'   databases = character(0),
#'   custom_gmt_files = c(CellTypes = gmt),
#'   export_excel = FALSE
#' )
#' hc_satellite(hc, "celltype_annotation")$selected_celltypes
#' @export
hc_celltype_annotation <- function(hc,
                                   databases = c("Descartes_Cell_Types_and_Tissue_2021", "Human_Gene_Atlas"),
                                   custom_gmt_files = NULL,
                                   clusters = c("all"),
                                   mode = c("coarse", "fine"),
                                   top = 3,
                                   qval = 0.1,
                                   padj = "BH",
                                   min_term_genes = 5,
                                   min_gs_size = 10,
                                   max_gs_size = 5000,
                                   annotation_slot = c("auto", "enriched_per_cluster", "enriched_per_cluster2"),
                                   slot_suffix = NULL,
                                   clear_previous_slots = TRUE,
                                   coarse_map = NULL,
                                   coarse_include_other = TRUE,
                                   refresh_db = FALSE,
                                   export_excel = TRUE,
                                   excel_file = "Module_Celltype_Annotation.xlsx",
                                   plot_heatmap = FALSE,
                                   heatmap_file_name = "Heatmap_modules_celltype_annotation.pdf",
                                   ...) {
  .hc_celltype_annotation_impl(
    hc = hc,
    databases = databases,
    custom_gmt_files = custom_gmt_files,
    clusters = clusters,
    mode = mode,
    top = top,
    qval = qval,
    padj = padj,
    min_term_genes = min_term_genes,
    min_gs_size = min_gs_size,
    max_gs_size = max_gs_size,
    annotation_slot = annotation_slot,
    slot_suffix = slot_suffix,
    clear_previous_slots = clear_previous_slots,
    coarse_map = coarse_map,
    coarse_include_other = coarse_include_other,
    refresh_db = refresh_db,
    export_excel = export_excel,
    excel_file = excel_file,
    plot_heatmap = plot_heatmap,
    heatmap_file_name = heatmap_file_name,
    ...
  )[["hc"]]
}

#' @noRd
# Internal implementation shared by the S4 entry point and the driver.
.hc_celltype_activity_decoupler_impl <- function(hc, ...) {
  out <- .hc_run_driver_capture(
    hc = hc,
    fun = ".hc_celltype_activity_decoupler_driver",
    ...
  )
  list(hc = out$hc, result = out$result)
}
#
#' Module cell-type activity from Enrichr markers via decoupleR (S4 API)
#'
#' @param hc A `HCoCenaExperiment`.
#' @param ... Further arguments passed through to the underlying
#'   implementation. Every supported argument is documented on this page;
#'   `...` only carries the tail of the argument list.
#' @return Updated `HCoCenaExperiment`.
#' @rdname hc_celltype_activity_decoupler
#' @param databases Character vector of Enrichr library names.
#' @param custom_gmt_files Optional custom GMT file(s) to include as additional
#' marker resources. Accepts a character vector or named list of file paths.
#' @param clusters Either `"all"` (default) or a character vector of module IDs/colors.
#' @param mode Either `"coarse"` (default) or `"fine"`.
#' @param top Number of selected marker activities per module.
#' @param qval Adjusted p-value cutoff.
#' @param padj Multiple-testing correction method.
#' @param activity_input One of `"gfc"` (default), `"fc"`, `"expression"`.
#' @param fc_comparisons Vector like `c("A_vs_B", "C_vs_B")` when `activity_input = "fc"`.
#' @param method decoupleR method suffix (e.g. `"ulm"`).
#' @param minsize Minimum target size passed to decoupleR.
#' @param min_term_genes Minimum genes per Enrichr term.
#' @param annotation_slot Legacy slot override for single-database runs.
#' @param slot_suffix Optional slot suffix; default `"decoupler"`.
#' @param clear_previous_slots If TRUE, clears old `enriched_per_cluster*` slots first.
#' @param coarse_map Optional named regex map used in coarse mode.
#' @param coarse_include_other Keep `"Other"` class in coarse mode.
#' @param refresh_db Refresh Enrichr cache.
#' @param export_excel Write summary workbook.
#' @param excel_file Excel filename.
#' @param plot_heatmap Replot cluster heatmap with dynamic slots.
#' @param heatmap_file_name Heatmap filename when `plot_heatmap = TRUE`.
#' @examples
#' if (requireNamespace("decoupleR", quietly = TRUE)) {
#'   hc <- hc_example_data("clustered")
#'   gmt <- system.file("extdata", "toy_celltype_markers.gmt",
#'                      package = "hcocena")
#'   hc <- hc_celltype_activity_decoupler(
#'     hc,
#'     databases = character(0),
#'     custom_gmt_files = c(CellTypes = gmt),
#'     export_excel = FALSE
#'   )
#'   names(hc_satellite(hc, "celltype_activity_decoupler"))
#' }
#' @export
hc_celltype_activity_decoupler <- function(hc,
                                           databases = c("Descartes_Cell_Types_and_Tissue_2021", "Human_Gene_Atlas"),
                                           custom_gmt_files = NULL,
                                           clusters = c("all"),
                                           mode = c("coarse", "fine"),
                                           top = 3,
                                           qval = 0.1,
                                           padj = "BH",
                                           activity_input = c("gfc", "fc", "expression"),
                                           fc_comparisons = NULL,
                                           method = "ulm",
                                           minsize = 5,
                                           min_term_genes = 5,
                                           annotation_slot = c("auto", "enriched_per_cluster", "enriched_per_cluster2"),
                                           slot_suffix = "decoupler",
                                           clear_previous_slots = FALSE,
                                           coarse_map = NULL,
                                           coarse_include_other = TRUE,
                                           refresh_db = FALSE,
                                           export_excel = TRUE,
                                           excel_file = "Module_Celltype_Activity_decoupler.xlsx",
                                           plot_heatmap = FALSE,
                                           heatmap_file_name = "Heatmap_modules_celltype_activity_decoupler.pdf",
                                           ...) {
  .hc_celltype_activity_decoupler_impl(
    hc = hc,
    databases = databases,
    custom_gmt_files = custom_gmt_files,
    clusters = clusters,
    mode = mode,
    top = top,
    qval = qval,
    padj = padj,
    activity_input = activity_input,
    fc_comparisons = fc_comparisons,
    method = method,
    minsize = minsize,
    min_term_genes = min_term_genes,
    annotation_slot = annotation_slot,
    slot_suffix = slot_suffix,
    clear_previous_slots = clear_previous_slots,
    coarse_map = coarse_map,
    coarse_include_other = coarse_include_other,
    refresh_db = refresh_db,
    export_excel = export_excel,
    excel_file = excel_file,
    plot_heatmap = plot_heatmap,
    heatmap_file_name = heatmap_file_name,
    ...
  )[["hc"]]
}

#' @noRd
# Internal implementation shared by the S4 entry point and the driver.
.hc_plot_enrichment_upstream_network_impl <- function(hc, ...) {
  out <- .hc_run_driver_capture(
    hc = hc,
    fun = ".hc_plot_enrichment_upstream_network_driver",
    ...
  )
  list(hc = out$hc, result = out$result)
}
#
#' Plot module knowledge network from enrichment + upstream inference (S4 API)
#'
#' @param hc A `HCoCenaExperiment`.
#' @return Updated `HCoCenaExperiment`.
#' @rdname hc_plot_enrichment_upstream_network
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
#' @examples
#' if (requireNamespace("decoupleR", quietly = TRUE)) {
#'   hc <- hc_example_data("clustered")
#'   gmt <- system.file("extdata", "toy_celltype_markers.gmt",
#'                      package = "hcocena")
#'   hc <- hc_functional_enrichment(
#'     hc,
#'     gene_sets = character(0),
#'     custom_gmt_files = c(CellTypes = gmt),
#'     universe = "network"
#'   )
#'   hc <- hc_upstream_inference(
#'     hc,
#'     resources = "Pathway",
#'     custom_pathway_gmt = c(CellTypes = gmt),
#'     plot = FALSE,
#'     save_pdf = FALSE
#'   )
#'   hc <- hc_plot_enrichment_upstream_network(hc, save_pdf = FALSE)
#' }
#' @export
hc_plot_enrichment_upstream_network <- function(hc,
                                                enrichment_mode = "selected",
                                                upstream_mode = "selected",
                                                clusters = c("all"),
                                                max_enrichment_per_module = NULL,
                                                max_upstream_per_module = NULL,
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
  .hc_plot_enrichment_upstream_network_impl(
    hc = hc,
    enrichment_mode = enrichment_mode,
    upstream_mode = upstream_mode,
    clusters = clusters,
    max_enrichment_per_module = max_enrichment_per_module,
    max_upstream_per_module = max_upstream_per_module,
    label_mode = label_mode,
    show_plot = show_plot,
    save_pdf = save_pdf,
    pdf_name = pdf_name,
    gfc_scale_limits = gfc_scale_limits,
    col_order = col_order,
    heatmap_col_order = heatmap_col_order,
    cluster_columns = cluster_columns,
    heatmap_cluster_columns = heatmap_cluster_columns,
    pdf_width = pdf_width,
    pdf_height = pdf_height,
    pdf_pointsize = pdf_pointsize,
    overall_plot_scale = overall_plot_scale
  )[["hc"]]
}

#' Run `check_dirs()` with S4 state synchronization
#'
#' @noRd
# Internal implementation shared by S4 and legacy entry points.
.hc_check_dirs_impl <- function(hc, create_output_dir = TRUE) {
  if (!inherits(hc, "HCoCenaExperiment")) {
    stop("`hc` must be a `HCoCenaExperiment`.")
  }

  paths <- .hc_row_to_list(hc@config@paths)
  if (base::length(paths) == 0) {
    return(hc)
  }

  for (nm in base::names(paths)) {
    current_dir <- paths[[nm]]
    if (nm == "dir_output" &&
      !base::identical(current_dir, FALSE) &&
      isTRUE(create_output_dir) &&
      !base::dir.exists(current_dir)) {
      base::dir.create(current_dir, recursive = TRUE, showWarnings = FALSE)
      message("Created missing output directory: ", current_dir)
    }
    paths[[nm]] <- fix_dir(current_dir)
  }

  hc@config@paths <- .hc_to_data_frame(paths)
  methods::validObject(hc)
  hc
}
#
#' Fixes Directories
#'
#' Iteratively calls fix_dir() on the provided working directory paths to fix them if necessary or throw an error if they are invalid.
#' @param hc A `HCoCenaExperiment`.
#' @param create_output_dir Boolean. If TRUE and `dir_output` is missing, create it.
#' @return Updated `HCoCenaExperiment`.
#' @examples
#' hc <- hc_example_data("prepared")
#' hc <- hc_set_paths(
#'   hc,
#'   dir_count_data = FALSE,
#'   dir_annotation = FALSE,
#'   dir_reference_files = FALSE,
#'   dir_output = tempdir()
#' )
#' hc <- hc_check_dirs(hc)
#' @export
hc_check_dirs <- function(hc, create_output_dir = TRUE) {
  .hc_check_dirs_impl(hc = hc, create_output_dir = create_output_dir)
}

#' Run `init_save_folder()` with S4 state synchronization
#'
#' @noRd
# Internal implementation shared by S4 and legacy entry points.
.hc_init_save_folder_impl <- function(hc, name, use_output_dir = FALSE) {
  if (!inherits(hc, "HCoCenaExperiment")) {
    stop("`hc` must be a `HCoCenaExperiment`.")
  }
  if (isTRUE(use_output_dir)) {
    name <- ""
  }
  if (base::is.null(name) || base::length(name) != 1 || base::is.na(name)) {
    stop("`name` must be a non-NA character scalar. Use \"\" to skip a subfolder.")
  }

  paths <- .hc_row_to_list(hc@config@paths)
  out_dir <- paths[["dir_output"]]
  if (base::is.null(out_dir) || base::identical(out_dir, FALSE) || !base::nzchar(base::as.character(out_dir))) {
    stop("`dir_output` is not set. Please run `hc_set_paths()` first.")
  }

  if (!base::dir.exists(out_dir)) {
    base::dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    message("Created output directory: ", out_dir)
  }

  save_folder <- base::as.character(name)
  if (base::identical(save_folder, "")) {
    message("Using output directory directly (no additional save subfolder): ", out_dir)
  } else {
    target_dir <- base::file.path(out_dir, save_folder)
    if (!base::dir.exists(target_dir)) {
      base::dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
      message("Created save folder: ", target_dir)
    } else {
      message("Using existing save folder: ", target_dir)
    }
  }

  global_cfg <- .hc_row_to_list(hc@config@global)
  global_cfg[["save_folder"]] <- name
  hc@config@global <- .hc_to_data_frame(global_cfg)
  methods::validObject(hc)
  hc
}
#
#' Creates a Save Folder
#'
#' A folder with the given name is created in the output directory. All analysis outputs will be saved to this folder.
#' @param hc A `HCoCenaExperiment`.
#' @param name Folder name. Use `""` to write directly into `dir_output`.
#' @param use_output_dir Boolean. If TRUE, ignore `name` and use `dir_output` directly.
#' @return Updated `HCoCenaExperiment`.
#' @examples
#' hc <- hc_example_data("prepared")
#' hc <- hc_init_save_folder(hc, name = "")
#' @export
hc_init_save_folder <- function(hc, name, use_output_dir = FALSE) {
  .hc_init_save_folder_impl(hc = hc, name = name, use_output_dir = use_output_dir)
}

#' Plot cut-off diagnostics (S4 API)
#'
#' @noRd
#' @param hc A `HCoCenaExperiment`.
#' @param ... Additional arguments for cut-off plotting.
#' @examples
#' hc <- hc_example_data("after_part1")
#' hc <- hc_plot_cutoffs(hc)
#' @return Updated `HCoCenaExperiment`.
# Internal implementation shared by S4 and legacy entry points.
.hc_plot_cutoffs_impl <- function(hc, ...) {
  out <- .hc_run_driver_capture(
    hc = hc,
    fun = ".hc_plot_cutoffs_driver",
    ...
  )
  list(hc = out$hc, result = out$result)
}
#
#' plot cutoffs
#'
#' @param hc A `HCoCenaExperiment`.
#' @param interactive Logical. If `TRUE` (default), render the cutoff
#'   statistics as an interactive plotly widget; if `FALSE`, as a static plot.
#' @param hline Named list of horizontal reference lines to draw, with the
#'   entries `R.squared`, `no_edges`, `no_nodes` and `no_networks`. `NULL`
#'   entries (the default) draw no line for that panel.
#' @return Updated `HCoCenaExperiment`, invisibly; called for the cut-off
#'   diagnostic plot it draws.
#' @examples
#' hc <- hc_example_data("after_part1")
#' hc <- hc_plot_cutoffs(hc, interactive = FALSE)
#' @export
hc_plot_cutoffs <- function(hc,
                            interactive = TRUE,
                            hline = list(
                              R.squared = NULL, no_edges = NULL,
                              no_nodes = NULL, no_networks = NULL
                            )) {
  .hc_plot_cutoffs_impl(hc = hc, interactive = interactive, hline = hline)[["hc"]]
}

#' Plot degree distributions (S4 API)
#'
#' @noRd
#' @param hc A `HCoCenaExperiment`.
#' @examples
#' hc <- hc_example_data("after_part1")
#' hc <- hc_plot_deg_dist(hc)
#' @return Updated `HCoCenaExperiment`.
#' @rdname hc_plot_deg_dist
# Internal implementation shared by S4 and legacy entry points.
.hc_plot_deg_dist_impl <- function(hc) {
  out <- .hc_run_driver_capture(
    hc = hc,
    fun = ".hc_plot_deg_dist_driver"
  )
  list(hc = out$hc, result = out$result)
}
#
#' Plot degree distributions (S4 API)
#'
#' Plots the logged degree distribution and its linear fit for each layer
#' at the chosen correlation cutoff, to check how well the network follows
#' a scale-free topology.
#'
#' @param hc A `HCoCenaExperiment`.
#' @return Updated `HCoCenaExperiment`.
#' @examples
#' hc <- hc_example_data("after_part1")
#' hc <- hc_plot_deg_dist(hc)
#' @export
hc_plot_deg_dist <- function(hc) {
  .hc_plot_deg_dist_impl(hc = hc)[["hc"]]
}

#' Plot module heatmap (S4 API)
#'
#' @param hc A `HCoCenaExperiment`.
#' @param file_name Optional export file name for the module heatmap.
#'   Defaults to `"Heatmap_modules.pdf"`. Use `FALSE` to skip file export.
#' @param ... Additional plotting arguments forwarded to the heatmap backend,
#'   including `smart_column_gaps`, `column_gap_by`, and `column_gap_mm`.
#' @examples
#' hc <- hc_example_data("clustered")
#' hc <- hc_plot_cluster_heatmap(hc, file_name = FALSE)
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_plot_cluster_heatmap <- function(hc, file_name = "Heatmap_modules.pdf", ...) {
  if (!inherits(hc, "HCoCenaExperiment")) {
    stop("`hc` must be a `HCoCenaExperiment`.")
  }
  dot_args <- list(...)
  if (!("module_label_numbering" %in% names(dot_args)) &&
    .hc_module_label_map_has_split_labels(hc@integration@cluster[["module_label_map"]])) {
    dot_args[["module_label_numbering"]] <- "preserve_existing"
  }

  legacy_envo <- .hc_bridge_state_env()
  legacy_state <- .hc_bind_bridge_hcobject(
    .hc_as_bridge_object_for_cluster_plot(hc),
    envo = legacy_envo
  )
  on.exit(.hc_restore_bridge_hcobject(legacy_state), add = TRUE)

  old_opt <- getOption("hcocena.suppress_legacy_warning", FALSE)
  options(hcocena.suppress_legacy_warning = TRUE)
  on.exit(options(hcocena.suppress_legacy_warning = old_opt), add = TRUE)

  base::do.call(.hc_plot_cluster_heatmap_driver, c(list(file_name = file_name), dot_args))
  .hc_update_hc_from_cluster_plot(
    hc = hc,
    hcobject = base::get("hcobject", envir = legacy_envo, inherits = FALSE)
  )
}

#' change grouping parameter
#'
#' Recalculate grouped GFCs and replot module heatmap (S4 API)
#'
#' Uses the existing module definitions from the current object and only
#' recalculates grouped GFCs for a different grouping variable. No integrated
#' network rebuild and no module re-clustering is performed.
#'
#' @param hc A `HCoCenaExperiment`.
#' @param group_by Grouping column to use (must be present in all annotation tables).
#' @param col_order Optional heatmap column order.
#' @param cluster_columns Whether to cluster heatmap columns. Default is
#'   `FALSE`, so the stored main hCoCena column order is reused when available.
#' @param row_order Optional heatmap row order.
#' @param cluster_rows Whether to cluster heatmap rows.
#' @examples
#' hc <- hc_example_data("clustered")
#' hc <- hc_change_grouping_parameter(hc, group_by = "batch")
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_change_grouping_parameter <- function(hc,
                                         group_by,
                                         col_order = NULL,
                                         cluster_columns = FALSE,
                                         row_order = NULL,
                                         cluster_rows = TRUE) {
  .hc_run_driver(
    hc = hc,
    fun = .hc_change_grouping_parameter_driver,
    group_by = group_by,
    col_order = col_order,
    cluster_columns = cluster_columns,
    row_order = row_order,
    cluster_rows = cluster_rows
  )
}

#' Plot integrated network (S4 API)
#'
#' @param hc A `HCoCenaExperiment`.
#' @param layout Optional pre-computed layout matrix. If `NULL` (default), the
#'   stored layout is reused or a new one is computed.
#' @param layout_algorithm Optional igraph layout function name, e.g.
#'   `"layout_with_fr"`. Overrides the layout set in the global settings.
#' @param gene_labels Optional character vector of genes to label in the plot.
#' @param save Logical. Write the network to PDF. Default is `TRUE`.
#' @param store_plot Logical. Keep the plot object in `hc`. Default is `FALSE`.
#' @param label_offset Distance between a node and its label. Default is 50.
#' @examples
#' hc <- hc_example_data("clustered")
#' hc <- hc_plot_integrated_network(hc)
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_plot_integrated_network <- function(hc, layout = NULL,
                                       layout_algorithm = NULL,
                                       gene_labels = NULL, save = TRUE,
                                       store_plot = FALSE,
                                       label_offset = 50) {
  .hc_run_driver(
    hc = hc, fun = .hc_plot_integrated_network_driver,
    layout = layout, layout_algorithm = layout_algorithm,
    gene_labels = gene_labels, save = save, store_plot = store_plot,
    label_offset = label_offset
  )
}

#' plot gfc network
#'
#' Plot network colored by GFC (S4 API)
#'
#' @param hc A `HCoCenaExperiment`.
#' @examples
#' hc <- hc_example_data("clustered")
#' hc <- hc_plot_gfc_network(hc)
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_plot_gfc_network <- function(hc) {
  .hc_run_driver(hc = hc, fun = .hc_plot_GFC_network_driver)
}

#' TF enrichment per module (S4 API)
#'
#' @param hc A `HCoCenaExperiment`.
#' @param clusters Either "all" (default) or a vector of module colours/labels
#'   for which the TF enrichment should be performed.
#' @param topTF Integer. Number of top ranking transcription factors to return
#'   per module. Default is 5.
#' @param topTarget Integer. Number of top ranking targets to return per
#'   transcription factor. Default is 5.
#' @return Updated `HCoCenaExperiment`.
#' @examples
#' \donttest{
#' # Queries the ChEA3 web service.
#' hc <- hc_example_data("clustered")
#' hc <- hc_tf_overrep_module(hc, topTF = 3, topTarget = 3)
#' }
#' @export
hc_tf_overrep_module <- function(hc, clusters = "all", topTF = 5,
                                 topTarget = 5) {
  .hc_run_driver(
    hc = hc, fun = .hc_TF_overrep_module_driver,
    clusters = clusters, topTF = topTF, topTarget = topTarget
  )
}

#' TF enrichment network-wide (S4 API)
#'
#' @param hc A `HCoCenaExperiment`.
#' @param topTF Integer. Number of transcription factors with the highest
#'   number of enriched targets network-wide. Default is 100.
#' @param topTarget Integer. Number of top enriched targets to return per
#'   transcription factor. Default is 30.
#' @return Updated `HCoCenaExperiment`.
#' @examples
#' \donttest{
#' # Queries the ChEA3 web service.
#' hc <- hc_example_data("clustered")
#' hc <- hc_tf_overrep_network(hc, topTF = 10, topTarget = 5)
#' names(hc_satellite(hc, "tf_network_targets"))
#' }
#' @export
hc_tf_overrep_network <- function(hc, topTF = 100, topTarget = 30) {
  .hc_run_driver(
    hc = hc, fun = .hc_TF_overrep_network_driver,
    topTF = topTF, topTarget = topTarget
  )
}

#' Check TF targets (S4 API)
#'
#' @param hc A `HCoCenaExperiment`.
#' @param TF Transcription factor symbol.
#' @return Updated `HCoCenaExperiment`.
#' @examples
#' \donttest{
#' # Needs the ChEA3 results of hc_tf_overrep_network(), which queries the
#' # ChEA3 web service.
#' hc <- hc_example_data("clustered")
#' hc <- hc_tf_overrep_network(hc, topTF = 10, topTarget = 5)
#' tfs <- names(hc_satellite(hc, "tf_network_targets"))
#' if (length(tfs) > 0) {
#'   hc <- hc_check_tf(hc, TF = tfs[[1]])
#' }
#' }
#' @export
hc_check_tf <- function(hc, TF) {
  .hc_run_driver(hc = hc, fun = .hc_check_tf_driver, TF = TF)
}

#' write session info
#'
#' Write session info (S4 API)
#'
#' @param hc A `HCoCenaExperiment`.
#' @examples
#' hc <- hc_example_data("prepared")
#' hc <- hc_write_session_info(hc)
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_write_session_info <- function(hc) {
  .hc_run_driver(hc = hc, fun = .hc_write_session_info_driver)
}

#' suggest topvar
#'
#' Suggest top variable genes (S4 API)
#'
#' @param hc A `HCoCenaExperiment`.
#' @details `hc_suggest_topvar()` operates on the counts stored in `hc`, i.e.
#'   after `hc_read_data()` preprocessing. Because `hc_read_data()` removes
#'   zero-variance genes and drops non-numeric helper columns from object-based
#'   count inputs, suggested inflection points can differ slightly from older
#'   releases on the same raw input even though the `hc_suggest_topvar()`
#'   heuristic itself is unchanged.
#' @examples
#' hc <- hc_example_data("prepared")
#' hc <- hc_suggest_topvar(hc)
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_suggest_topvar <- function(hc) {
  .hc_run_driver(hc = hc, fun = .hc_suggest_topvar_driver)
}

#' plot sample distributions
#'
#' Plot sample distributions (S4 API)
#'
#' @param hc A `HCoCenaExperiment`.
#' @param plot_type Either "boxplot" (default) or "freqdist" for a
#'   distribution-focused view.
#' @param log_2 Logical. Log2-transform the values before plotting.
#'   Default is `TRUE`.
#' @param plot Logical. Draw the plot. Default is `TRUE`.
#' @examples
#' hc <- hc_example_data("prepared")
#' hc <- hc_plot_sample_distributions(hc)
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_plot_sample_distributions <- function(hc, plot_type = "boxplot",
                                         log_2 = TRUE, plot = TRUE) {
  .hc_run_driver(
    hc = hc, fun = .hc_plot_sample_distributions_driver,
    plot_type = plot_type, log_2 = log_2, plot = plot
  )
}

#' pca
#'
#' PCA plotting (S4 API)
#'
#' @param hc A `HCoCenaExperiment`.
#' @param which Which gene set to run the PCA on: "all" (default, all genes),
#'   "topvar" (the top-variable genes of each layer) or "network_genes" (only
#'   genes that ended up in the integrated network).
#' @param color_by Optional annotation column name used to colour the samples.
#'   If `NULL` (default), the variable of interest from the global settings is
#'   used.
#' @param ellipses Logical. Whether to draw group confidence ellipses.
#'   Default is `FALSE`.
#' @param cols Optional named vector of colours for the groups. If `NULL`
#'   (default), a built-in palette is used.
#' @examples
#' hc <- hc_example_data("prepared")
#' hc <- hc_pca(hc)
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_pca <- function(hc, which = "all", color_by = NULL, ellipses = FALSE,
                   cols = NULL) {
  .hc_run_driver(
    hc = hc, fun = .hc_PCA_driver,
    which = which, color_by = color_by, ellipses = ellipses, cols = cols
  )
}

#' meta plot
#'
#' Meta-data plotting (S4 API)
#'
#' @param hc A `HCoCenaExperiment`.
#' @param set An integer. Number of the dataset/layer to inspect.
#' @param group_col Annotation column used to group the samples, e.g. the
#'   timepoint or condition.
#' @param meta_col Annotation column whose distribution is plotted.
#' @param type Either "cat" (default) for categorical `meta_col` or "num" for
#'   continuous ones.
#' @param cols Optional named vector of colours for the groups.
#' @examples
#' hc <- hc_example_data("prepared")
#' hc <- hc_meta_plot(hc, set = 1, group_col = "group", meta_col = "batch")
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_meta_plot <- function(hc, set, group_col = NULL, meta_col = NULL,
                         type = "cat", cols = NULL) {
  .hc_run_driver(
    hc = hc, fun = .hc_meta_plot_driver,
    set = set, group_col = group_col, meta_col = meta_col,
    type = type, cols = cols
  )
}

#' export clusters
#'
#' Export clusters (S4 API)
#'
#' @param hc A `HCoCenaExperiment`.
#' @examples
#' hc <- hc_example_data("clustered")
#' hc <- hc_export_clusters(hc)
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_export_clusters <- function(hc) {
  .hc_run_driver(hc = hc, fun = .hc_export_clusters_driver)
}

#' get module scores
#'
#' Module scores (S4 API)
#'
#' @param hc A `HCoCenaExperiment`.
#' @param save Logical. Write the module-score box plot to PDF.
#'   Default is `TRUE`.
#' @examples
#' hc <- hc_example_data("clustered")
#' hc <- hc_get_module_scores(hc)
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_get_module_scores <- function(hc, save = TRUE) {
  .hc_run_driver(hc = hc, fun = .hc_get_module_scores_driver, save = save)
}

#' algo alluvial
#'
#' Alluvial comparison plots (S4 API)
#'
#' @param hc A `HCoCenaExperiment`.
#' @examples
#' hc <- hc_example_data("clustered")
#' hc <- hc_algo_alluvial(hc)
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_algo_alluvial <- function(hc) {
  .hc_run_driver(hc = hc, fun = .hc_algo_alluvial_driver)
}

#' pca algo compare
#'
#' PCA algorithm comparison (S4 API)
#'
#' @param hc A `HCoCenaExperiment`.
#' @param gtc Optional gene-to-cluster table to compare against. If `NULL`
#'   (default), the current module assignment is used.
#' @param algo Optional name of the clustering algorithm the `gtc` came from,
#'   used for the plot title.
#' @param cols Optional named vector of colours for the sample groups.
#' @examples
#' hc <- hc_example_data("clustered")
#' hc <- hc_pca_algo_compare(hc)
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_pca_algo_compare <- function(hc, gtc = NULL, algo = NULL, cols = NULL) {
  .hc_run_driver(
    hc = hc, fun = .hc_PCA_algo_compare_driver,
    gtc = gtc, algo = algo, cols = cols
  )
}

#' update clustering algorithm
#'
#' Update clustering algorithm (S4 API)
#'
#' @param hc A `HCoCenaExperiment`.
#' @param new_algo Name of the clustering algorithm to switch to, e.g.
#'   `"cluster_louvain"`. Mutually exclusive with `gtc`.
#' @param gtc Optional externally supplied gene-to-cluster table to adopt
#'   instead of re-running an algorithm.
#' @examples
#' hc <- hc_example_data("clustered")
#' hc <- hc_update_clustering_algorithm(hc, new_algo = "cluster_louvain")
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_update_clustering_algorithm <- function(hc, new_algo = NULL, gtc = NULL) {
  .hc_run_driver(
    hc = hc, fun = .hc_update_clustering_algorithm_driver,
    new_algo = new_algo, gtc = gtc
  )
}

#' Export network to local folder (S4 API)
#'
#' @param hc A `HCoCenaExperiment`.
#' @param file Target folder for the exported network files. If omitted, the
#'   current output/save folder configured in `hc` is used.
#' @return Updated `HCoCenaExperiment`.
#' @examples
#' hc <- hc_example_data("clustered")
#' hc <- hc_export_to_local_folder(hc, file = tempdir())
#' @export
hc_export_to_local_folder <- function(hc, file) {
  # `file`'s default in the driver is built from the live `hcobject`, which only
  # exists inside the driver's frame - so forward it only when supplied instead
  # of copying the default up here.
  if (missing(file)) {
    .hc_run_driver(hc = hc, fun = .hc_export_to_local_folder_driver)
  } else {
    .hc_run_driver(hc = hc, fun = .hc_export_to_local_folder_driver, file = file)
  }
}

#' Import Cytoscape layout from local folder (S4 API)
#'
#' @param hc A `HCoCenaExperiment`.
#' @param file Path to the `network_layout.csv` exported from Cytoscape. If
#'   omitted, the current output/save folder is used.
#' @return Updated `HCoCenaExperiment`.
#' @examples
#' hc <- hc_example_data("clustered")
#' # A layout CSV as written by hc_import_layout_from_cytoscape(): one row per
#' # node, named by gene, with x and y coordinates.
#' g <- hc_graph(hc)
#' xy <- igraph::layout_with_fr(g)
#' rownames(xy) <- igraph::V(g)$name
#' layout_file <- file.path(tempdir(), "network_layout.csv")
#' utils::write.csv(xy, layout_file)
#' hc <- hc_import_layout_from_local_folder(hc, file = layout_file)
#' @export
hc_import_layout_from_local_folder <- function(hc, file) {
  # the driver's default is built from the live `hcobject`, so only forward
  # `file` when the caller actually supplied one
  if (missing(file)) {
    .hc_run_driver(hc = hc, fun = .hc_import_layout_from_local_folder_driver)
  } else {
    .hc_run_driver(
      hc = hc, fun = .hc_import_layout_from_local_folder_driver, file = file
    )
  }
}

#' Export network to Cytoscape (S4 API)
#'
#' @param hc A `HCoCenaExperiment`.
#' @param name Name under which the network is created in Cytoscape.
#'   Default is "my igraph".
#' @param docker_container Logical. Set to `TRUE` when running hCoCena inside a
#'   Docker container that talks to a Cytoscape instance on the host.
#'   Default is `FALSE`.
#' @return Updated `HCoCenaExperiment`.
#' @examples
#' \donttest{
#' # Needs a running Cytoscape instance reachable through RCy3.
#' cytoscape_up <- requireNamespace("RCy3", quietly = TRUE) &&
#'   isTRUE(tryCatch({
#'     RCy3::cytoscapePing()
#'     TRUE
#'   }, error = function(e) FALSE))
#' if (cytoscape_up) {
#'   hc <- hc_example_data("clustered")
#'   hc <- hc_export_to_cytoscape(hc, name = "hcocena example")
#' }
#' }
#' @export
hc_export_to_cytoscape <- function(hc, name = "my igraph",
                                   docker_container = FALSE) {
  .hc_run_driver(
    hc = hc, fun = .hc_export_to_cytoscape_driver,
    name = name, docker_container = docker_container
  )
}

#' Import layout from Cytoscape (S4 API)
#'
#' @param hc A `HCoCenaExperiment`.
#' @return Updated `HCoCenaExperiment`.
#' @examples
#' \donttest{
#' # Needs a running Cytoscape instance reachable through RCy3.
#' cytoscape_up <- requireNamespace("RCy3", quietly = TRUE) &&
#'   isTRUE(tryCatch({
#'     RCy3::cytoscapePing()
#'     TRUE
#'   }, error = function(e) FALSE))
#' if (cytoscape_up) {
#'   hc <- hc_example_data("clustered")
#'   hc <- hc_export_to_cytoscape(hc, name = "hcocena example")
#'   # ... arrange the network in Cytoscape, then:
#'   hc <- hc_import_layout_from_cytoscape(hc)
#' }
#' }
#' @export
hc_import_layout_from_cytoscape <- function(hc) {
  .hc_run_driver(hc = hc, fun = .hc_import_layout_from_cytoscape_driver)
}

#' Hub detection (S4 API)
#'
#' Determines hub genes per module from a combined ranking of weighted degree,
#' closeness and betweenness centrality.
#'
#' @param hc A `HCoCenaExperiment`.
#' @param clusters Either "all" (default) or a vector of module colours for
#'   which hub detection should be performed.
#' @param top Integer. Number of top-ranked genes per module to report as hubs.
#'   Default is 10.
#' @param tree_layout Logical. Draw the per-module hub network with a tree
#'   layout. Default is `FALSE`.
#' @param TF_only Either `FALSE` (default, consider all genes), `"all"` (only
#'   genes from the transcription-factor reference file), or one gene category
#'   from that file's last column.
#' @param save Logical. Write the hub network and expression heatmap to PDF.
#'   Default is `FALSE`.
#' @param plot Logical. Draw the per-module network. Default is `FALSE`.
#' @examples
#' hc <- hc_example_data("clustered")
#' hc <- hc_find_hubs(hc)
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_find_hubs <- function(hc, clusters = c("all"), top = 10,
                         tree_layout = FALSE, TF_only = FALSE,
                         save = FALSE, plot = FALSE) {
  .hc_run_driver(
    hc = hc, fun = .hc_find_hubs_driver,
    clusters = clusters, top = top, tree_layout = tree_layout,
    TF_only = TF_only, save = save, plot = plot
  )
}

#' Visualize gene expression (S4 API)
#'
#' Plots the mean expression per condition for a set of genes as a heatmap.
#'
#' @param hc A `HCoCenaExperiment`.
#' @param genes Character vector of gene symbols to plot.
#' @param name Optional name used for the plot title and output file. If `NULL`
#'   (default), a name is derived automatically.
#' @param width Plot width in inches. Default is 15.
#' @param height Plot height in inches. Default is 10.
#' @param save Logical. Write the heatmap to PDF. Default is `TRUE`.
#' @examples
#' hc <- hc_example_data("clustered")
#' hc <- hc_visualize_gene_expression(hc, genes = "CD3D")
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_visualize_gene_expression <- function(hc, genes, name = NULL, width = 15,
                                         height = 10, save = TRUE) {
  .hc_run_driver(
    hc = hc, fun = .hc_visualize_gene_expression_driver,
    genes = genes, name = name, width = width, height = height, save = save
  )
}

#' Highlight gene set in network (S4 API)
#'
#' Re-plots the integrated network with a gene set of interest highlighted.
#'
#' @param hc A `HCoCenaExperiment`.
#' @param gene_set Character vector of gene symbols to highlight.
#' @param name Optional name used for the plot title and output file.
#' @param col Colour used for the highlighted genes. Default is `"black"`.
#' @param label_offset Distance between a node and its label. Default is 3.
#' @param plot Logical. Draw the network. Default is `TRUE`.
#' @param save Logical. Write the network to PDF. Default is `TRUE`.
#' @examples
#' hc <- hc_example_data("clustered")
#' hc <- hc_highlight_geneset(hc, gene_set = c("CD3D", "CD19"))
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_highlight_geneset <- function(hc, gene_set, name = NULL, col = "black",
                                 label_offset = 3, plot = TRUE, save = TRUE) {
  .hc_run_driver(
    hc = hc, fun = .hc_highlight_geneset_driver,
    gene_set = gene_set, name = name, col = col,
    label_offset = label_offset, plot = plot, save = save
  )
}

#' Highlight single cluster in network (S4 API)
#'
#' Re-plots the integrated network with one module emphasised by colour and
#' node size.
#'
#' @param hc A `HCoCenaExperiment`.
#' @param cluster The module to highlight, given as its colour or module label.
#' @examples
#' hc <- hc_example_data("clustered")
#' hc <- hc_colour_single_cluster(hc, cluster = "gold")
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_colour_single_cluster <- function(hc, cluster) {
  .hc_run_driver(
    hc = hc, fun = .hc_colour_single_cluster_driver, cluster = cluster
  )
}

#' Add categorical module heatmap annotations (S4 API)
#'
#' @param hc A `HCoCenaExperiment`.
#' @param variables Character vector of categorical annotation column names to
#'   add as heatmap column annotations.
#' @param variable_label Optional display label(s) for `variables`. If `NULL`
#'   (default), the column names themselves are used.
#' @param type Either "abs" (default) for absolute counts or "rel" for
#'   relative proportions.
#' @return Updated `HCoCenaExperiment`.
#' @examples
#' hc <- hc_example_data("clustered")
#' hc <- hc_col_anno_categorical(hc, variables = "batch")
#' @export
hc_col_anno_categorical <- function(hc, variables, variable_label = NULL,
                                    type = "abs") {
  .hc_run_driver(
    hc = hc, fun = .hc_col_anno_categorical_driver,
    variables = variables, variable_label = variable_label, type = type
  )
}

#' Correlate categorical metadata with modules (S4 API)
#'
#' @param hc A `HCoCenaExperiment`.
#' @param meta A single string. Name of the categorical annotation column to
#'   correlate with the module expression patterns.
#' @param set An integer. Number of the dataset/layer to use.
#' @param p_val Maximum adjusted p-value for a correlation to count as
#'   significant. Default is 0.05. Non-significant cells are shown in grey.
#' @param padj Multiple-testing correction method passed to
#'   [stats::p.adjust()]. Default is "BH".
#' @return Updated `HCoCenaExperiment`.
#' @examples
#' hc <- hc_example_data("clustered")
#' hc_meta_correlation_cat(hc, set = 1, meta = "batch")
#' @export
hc_meta_correlation_cat <- function(hc, meta, set, p_val = 0.05,
                                    padj = "BH") {
  .hc_run_driver(
    hc = hc, fun = .hc_meta_correlation_cat_driver,
    meta = meta, set = set, p_val = p_val, padj = padj
  )
}

#' Test module differences between conditions (S4 API)
#'
#' @param hc A `HCoCenaExperiment`.
#' @return Updated `HCoCenaExperiment`.
#' @param set Integer vector of layer indices or `"all"` (default).
#' @param condition_col Column name in annotation files defining conditions.
#'   If `NULL`, uses `hcobject[["global_settings"]][["voi"]]`.
#' @param donor_col Optional donor identifier column in annotation files.
#'   Required if `run_lmm = TRUE`.
#' @param time_col Optional time column in annotation files (used by LMM when
#'   available and `lmm_include_time = TRUE`).
#' @param run_wilcox Logical; run Wilcoxon/Kruskal module-wise test.
#' @param run_limma Logical; run limma pairwise contrasts.
#' @param run_lmm Logical; run linear mixed model (`nlme::lme`) per module.
#' @param lmm_include_time Logical; include `time_col` as fixed effect in LMM
#'   when available.
#' @param limma_reference Optional reference condition for one-vs-reference
#'   limma contrasts. If `NULL`, all pairwise contrasts are tested.
#' @param limma_trend Logical; passed to `limma::eBayes(trend = ...)`.
#' @param standardize_modules Logical; z-score each module across the samples
#'   of a layer before testing. Raw module means carry the module's absolute
#'   expression level, so limma's variance moderation - which borrows
#'   information across modules - is dominated by the highly expressed ones,
#'   and the reported effect sizes are not comparable between modules.
#'
#'   The standardization is applied per layer, so with several layers pooled
#'   (`set = "all"`) it also removes the between-layer offset. With a single
#'   layer the rank-based Wilcoxon/Kruskal p-values are unchanged, because
#'   z-scoring is monotone within a row; pooled across layers they do change,
#'   because the layers are brought onto a common scale first. That is usually
#'   what you want when the layers are different tissues or assays, since
#'   neither the limma design (`~ 0 + condition`) nor the LMM
#'   (`value ~ condition`, random `~ 1 | donor`) carries a layer term.
#'
#'   Default `FALSE` keeps the previous behaviour.
#' @param padj Multiple-testing correction method (passed to `p.adjust`/limma).
#' @param export_excel Logical; write result tables to Excel in save folder.
#' @param excel_file File name for the Excel export.
#' @param slot_name Satellite slot name for storing results.
#' @examples
#' hc <- hc_example_data("clustered")
#' hc <- hc_module_condition_significance(
#'   hc,
#'   condition_col = "group",
#'   run_limma = FALSE,
#'   export_excel = FALSE
#' )
#' @export
hc_module_condition_significance <- function(hc,
                                             set = "all",
                                             condition_col = NULL,
                                             donor_col = NULL,
                                             time_col = NULL,
                                             run_wilcox = TRUE,
                                             run_limma = TRUE,
                                             run_lmm = FALSE,
                                             lmm_include_time = TRUE,
                                             limma_reference = NULL,
                                             limma_trend = TRUE,
                                             standardize_modules = FALSE,
                                             padj = "BH",
                                             export_excel = TRUE,
                                             excel_file = "Module_condition_significance.xlsx",
                                             slot_name = "module_condition_significance") {
  .hc_run_driver(
    hc = hc, fun = .hc_module_condition_significance_driver,
    set = set, condition_col = condition_col, donor_col = donor_col,
    time_col = time_col, run_wilcox = run_wilcox, run_limma = run_limma,
    run_lmm = run_lmm, lmm_include_time = lmm_include_time,
    limma_reference = limma_reference, limma_trend = limma_trend,
    standardize_modules = standardize_modules,
    padj = padj, export_excel = export_excel, excel_file = excel_file,
    slot_name = slot_name
  )
}
