# Internal helpers for user-provided GMT gene sets.

.hc_normalize_custom_gmt_files <- function(custom_gmt_files,
                                           default_prefix = "CustomGMT") {
  if (base::is.null(custom_gmt_files)) {
    return(stats::setNames(base::character(0), base::character(0)))
  }

  files <- custom_gmt_files
  if (base::is.list(files)) {
    files <- base::unlist(files, recursive = TRUE, use.names = TRUE)
  }
  # as.character() drops names, and the names are the database labels the
  # user asked for; carry them over explicitly.
  files_nm <- base::names(files)
  files <- base::as.character(files)
  base::names(files) <- files_nm
  if (base::length(files) == 0) {
    return(stats::setNames(base::character(0), base::character(0)))
  }

  keep <- !base::is.na(files) & base::trimws(files) != ""
  files <- base::trimws(files[keep])
  if (base::length(files) == 0) {
    return(stats::setNames(base::character(0), base::character(0)))
  }

  nm <- base::names(files)
  if (base::is.null(nm)) {
    nm <- base::rep("", base::length(files))
  }
  nm <- base::trimws(base::as.character(nm))
  fallback <- base::paste0(default_prefix, base::seq_along(files))
  nm[nm == "" | base::is.na(nm)] <- fallback[nm == "" | base::is.na(nm)]
  nm <- gsub("\\s+", "_", nm, perl = TRUE)
  nm <- gsub("[^[:alnum:]_\\-.]", "_", nm, perl = TRUE)
  nm <- base::make.unique(nm, sep = "_")

  stats::setNames(files, nm)
}

.hc_resolve_custom_gmt_path <- function(path) {
  path <- base::as.character(path[[1]])
  if (base::is.na(path) || base::trimws(path) == "") {
    stop("Invalid custom GMT path.")
  }
  path <- base::trimws(path)

  if (base::file.exists(path)) {
    return(base::normalizePath(path, winslash = "/", mustWork = TRUE))
  }

  ref_dir <- tryCatch(
    hcobject[["working_directory"]][["dir_reference_files"]],
    error = function(e) NULL
  )
  if (!base::is.null(ref_dir) && !identical(ref_dir, FALSE)) {
    ref_dir <- base::as.character(ref_dir[[1]])
    if (!base::is.na(ref_dir) && base::trimws(ref_dir) != "") {
      candidate <- base::file.path(ref_dir, path)
      if (base::file.exists(candidate)) {
        return(base::normalizePath(candidate, winslash = "/", mustWork = TRUE))
      }
    }
  }

  stop(
    "Custom GMT file not found: ", path,
    ". Use an absolute/relative path or place the file in `dir_reference_files`."
  )
}

.hc_load_custom_gmt_term2gene <- function(custom_gmt_files,
                                          default_prefix = "CustomGMT",
                                          uppercase_genes = FALSE,
                                          add_database_prefix = FALSE,
                                          title_case_names = FALSE) {
  files <- .hc_normalize_custom_gmt_files(
    custom_gmt_files = custom_gmt_files,
    default_prefix = default_prefix
  )
  if (base::length(files) == 0) {
    return(list())
  }
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    stop("Package `clusterProfiler` is required to read custom GMT files.")
  }

  db_names <- base::names(files)
  if (isTRUE(title_case_names)) {
    db_names <- stringr::str_to_title(db_names)
  }
  db_names <- base::make.unique(db_names, sep = "_")

  out <- vector("list", base::length(files))
  base::names(out) <- db_names
  resolved_paths <- stats::setNames(base::character(base::length(files)), db_names)

  for (i in base::seq_along(files)) {
    db_nm <- db_names[[i]]
    gmt_path <- .hc_resolve_custom_gmt_path(files[[i]])
    resolved_paths[[db_nm]] <- gmt_path

    t2g <- clusterProfiler::read.gmt(gmt_path)
    t2g <- t2g %>% base::as.data.frame(stringsAsFactors = FALSE)
    if (!all(c("term", "gene") %in% base::colnames(t2g))) {
      stop("Custom GMT does not contain expected `term`/`gene` columns: ", gmt_path)
    }
    t2g <- t2g[, c("term", "gene"), drop = FALSE]
    t2g$term <- base::trimws(base::as.character(t2g$term))
    t2g$gene <- base::trimws(base::as.character(t2g$gene))
    t2g <- t2g[
      !base::is.na(t2g$term) & t2g$term != "" &
        !base::is.na(t2g$gene) & t2g$gene != "", ,
      drop = FALSE
    ]
    if (isTRUE(uppercase_genes)) {
      t2g$gene <- base::toupper(t2g$gene)
    }
    t2g <- base::unique(t2g)

    if (isTRUE(add_database_prefix)) {
      term_raw <- t2g$term
      t2g <- base::data.frame(
        term = base::paste0("[", db_nm, "] ", term_raw),
        term_raw = term_raw,
        database = db_nm,
        gene = t2g$gene,
        stringsAsFactors = FALSE
      )
      t2g <- base::unique(t2g)
    }

    out[[db_nm]] <- t2g
  }

  base::attr(out, "paths") <- resolved_paths
  out
}

.hc_load_custom_gmt_pathway_network <- function(custom_pathway_gmt,
                                                default_prefix = "CustomPathway") {
  t2g_list <- .hc_load_custom_gmt_term2gene(
    custom_gmt_files = custom_pathway_gmt,
    default_prefix = default_prefix,
    uppercase_genes = TRUE,
    add_database_prefix = TRUE,
    title_case_names = FALSE
  )
  if (base::length(t2g_list) == 0) {
    out <- base::data.frame(
      source = base::character(0),
      target = base::character(0),
      mor = base::numeric(0),
      stringsAsFactors = FALSE
    )
    base::attr(out, "paths") <- stats::setNames(base::character(0), base::character(0))
    base::attr(out, "databases") <- base::character(0)
    return(out)
  }

  rows <- base::lapply(base::names(t2g_list), function(db_nm) {
    t2g <- t2g_list[[db_nm]]
    if (base::is.null(t2g) || base::nrow(t2g) == 0) {
      return(NULL)
    }
    base::data.frame(
      source = base::as.character(t2g$term),
      target = base::as.character(t2g$gene),
      mor = 1,
      stringsAsFactors = FALSE
    )
  })
  rows <- rows[!base::vapply(rows, base::is.null, FUN.VALUE = base::logical(1))]
  out <- if (base::length(rows) == 0) {
    base::data.frame(
      source = base::character(0),
      target = base::character(0),
      mor = base::numeric(0),
      stringsAsFactors = FALSE
    )
  } else {
    base::do.call(base::rbind, rows)
  }
  out <- base::unique(out)
  base::rownames(out) <- NULL
  base::attr(out, "paths") <- base::attr(t2g_list, "paths")
  base::attr(out, "databases") <- base::names(t2g_list)
  out
}
