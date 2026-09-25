#' List Enrichr cell-type related databases
#'
#' Lists Enrichr libraries and optionally filters for cell/tissue related names.
#'
#' @param pattern Optional regular expression applied to database names.
#' @param include_all Logical. If `FALSE` (default), applies a built-in
#'   cell-type/tissue oriented filter before `pattern`.
#' @param ignore_case Logical. Passed to [base::grepl()].
#'
#' @return A data frame with Enrichr library metadata.
#' @examples
#' hc_list_celltype_databases()
#' @export
hc_list_celltype_databases <- function(pattern = NULL,
                                       include_all = FALSE,
                                       ignore_case = TRUE) {
  dbs <- .hc_enrichr_db_table(refresh = FALSE)
  if (base::nrow(dbs) == 0) {
    return(dbs)
  }

  out <- dbs
  if (!isTRUE(include_all)) {
    default_pattern <- paste0(
      "(cell|atlas|immune|leukocyte|hematopo|blood|pbmc|",
      "descartes|panglao|tabula|single[[:space:]_-]*cell|",
      "lymph|myeloid|monocyte|macrophage|t[[:space:]_-]*cell|b[[:space:]_-]*cell|nk)"
    )
    keep <- base::grepl(default_pattern, out$database, ignore.case = TRUE, perl = TRUE)
    out <- out[keep, , drop = FALSE]
  }
  if (!base::is.null(pattern) && base::nzchar(base::as.character(pattern[[1]]))) {
    keep <- base::grepl(
      pattern = base::as.character(pattern[[1]]),
      x = out$database,
      ignore.case = isTRUE(ignore_case),
      perl = TRUE
    )
    out <- out[keep, , drop = FALSE]
  }
  out <- out[base::order(out$database), , drop = FALSE]
  base::rownames(out) <- NULL
  out
}

# Internal Enrichr cache.
.hc_enrichr_cache_env <- local({
  cache <- new.env(parent = emptyenv())
  cache$db_table <- NULL
  cache$libraries <- list()
  function() cache
})

.hc_enrichr_db_table <- function(refresh = FALSE, timeout_sec = 60) {
  cache <- .hc_enrichr_cache_env()
  if (!isTRUE(refresh) && !base::is.null(cache$db_table)) {
    return(cache$db_table)
  }

  req <- tryCatch(
    httr::GET("https://maayanlab.cloud/Enrichr/datasetStatistics", httr::timeout(timeout_sec)),
    error = function(e) NULL
  )
  if (base::is.null(req) || httr::status_code(req) >= 300) {
    warning("Could not fetch Enrichr database statistics.")
    out <- base::data.frame(
      database = base::character(0),
      category_id = base::integer(0),
      num_terms = base::integer(0),
      genes_per_term = base::numeric(0),
      gene_coverage = base::numeric(0),
      stringsAsFactors = FALSE
    )
    cache$db_table <- out
    return(out)
  }

  payload <- tryCatch(
    jsonlite::fromJSON(httr::content(req, as = "text", encoding = "UTF-8")),
    error = function(e) NULL
  )
  stats <- NULL
  if (!base::is.null(payload) && "statistics" %in% base::names(payload)) {
    stats <- payload$statistics
  }
  if (base::is.null(stats) || base::nrow(stats) == 0) {
    out <- base::data.frame(
      database = base::character(0),
      category_id = base::integer(0),
      num_terms = base::integer(0),
      genes_per_term = base::numeric(0),
      gene_coverage = base::numeric(0),
      stringsAsFactors = FALSE
    )
    cache$db_table <- out
    return(out)
  }

  out <- base::data.frame(
    database = base::as.character(stats$libraryName),
    category_id = .hc_as_integer_safely(stats$categoryId),
    num_terms = .hc_as_integer_safely(stats$numTerms),
    genes_per_term = .hc_as_numeric_safely(stats$genesPerTerm),
    gene_coverage = .hc_as_numeric_safely(stats$geneCoverage),
    stringsAsFactors = FALSE
  )
  out <- out[!base::is.na(out$database) & out$database != "", , drop = FALSE]
  out <- out[base::order(out$database), , drop = FALSE]
  base::rownames(out) <- NULL
  cache$db_table <- out
  out
}

.hc_enrichr_token_to_gene <- function(token) {
  token <- base::trimws(base::as.character(token))
  if (!base::nzchar(token)) {
    return(NA_character_)
  }
  if (base::startsWith(token, "http://") || base::startsWith(token, "https://")) {
    return(NA_character_)
  }
  gene <- sub(",.*$", "", token, perl = TRUE)
  gene <- base::toupper(base::trimws(gene))
  if (!base::nzchar(gene)) {
    return(NA_character_)
  }
  gene
}

.hc_enrichr_fetch_library_term2gene <- function(database,
                                                refresh = FALSE,
                                                timeout_sec = 120,
                                                min_genes_per_term = 3) {
  if (!base::is.character(database) || base::length(database) != 1 || base::is.na(database) || !base::nzchar(database)) {
    stop("`database` must be a single non-empty character string.")
  }
  if (!base::is.numeric(min_genes_per_term) || base::length(min_genes_per_term) != 1 || base::is.na(min_genes_per_term) || min_genes_per_term < 1) {
    stop("`min_genes_per_term` must be a positive integer.")
  }
  min_genes_per_term <- base::as.integer(min_genes_per_term)

  cache <- .hc_enrichr_cache_env()
  cache_key <- base::as.character(database)
  if (!isTRUE(refresh) &&
    !base::is.null(cache$libraries[[cache_key]]) &&
    base::is.list(cache$libraries[[cache_key]]) &&
    !base::is.null(cache$libraries[[cache_key]]$term2gene)) {
    return(cache$libraries[[cache_key]])
  }

  req <- tryCatch(
    httr::GET(
      "https://maayanlab.cloud/Enrichr/geneSetLibrary",
      query = list(mode = "text", libraryName = database),
      httr::timeout(timeout_sec)
    ),
    error = function(e) NULL
  )
  if (base::is.null(req) || httr::status_code(req) >= 300) {
    warning("Could not download Enrichr library: ", database)
    out <- list(
      database = database,
      term2gene = base::data.frame(
        term = base::character(0),
        term_raw = base::character(0),
        database = base::character(0),
        gene = base::character(0),
        stringsAsFactors = FALSE
      )
    )
    cache$libraries[[cache_key]] <- out
    return(out)
  }

  txt <- httr::content(req, as = "text", encoding = "UTF-8")
  lines <- unlist(base::strsplit(txt, "\n", fixed = TRUE))
  lines <- lines[base::nzchar(base::trimws(lines))]

  rows <- vector("list", base::length(lines))
  idx <- 0L
  for (ln in lines) {
    fields <- unlist(base::strsplit(ln, "\t", fixed = TRUE))
    if (base::length(fields) < 2) {
      next
    }
    term_raw <- base::trimws(base::as.character(fields[[1]]))
    if (!base::nzchar(term_raw)) {
      next
    }
    term <- base::paste0("[", database, "] ", term_raw)
    tokens <- fields[-1]
    genes <- base::vapply(tokens, .hc_enrichr_token_to_gene, FUN.VALUE = base::character(1))
    genes <- genes[!base::is.na(genes) & genes != ""]
    genes <- base::unique(genes)
    if (base::length(genes) < min_genes_per_term) {
      next
    }
    idx <- idx + 1L
    rows[[idx]] <- base::data.frame(
      term = term,
      term_raw = term_raw,
      database = database,
      gene = genes,
      stringsAsFactors = FALSE
    )
  }
  rows <- rows[base::seq_len(idx)]
  t2g <- if (base::length(rows) == 0) {
    base::data.frame(
      term = base::character(0),
      term_raw = base::character(0),
      database = base::character(0),
      gene = base::character(0),
      stringsAsFactors = FALSE
    )
  } else {
    base::do.call(base::rbind, rows)
  }
  if (base::nrow(t2g) > 0) {
    t2g <- t2g[!base::is.na(t2g$gene) & t2g$gene != "", , drop = FALSE]
    t2g <- base::unique(t2g)
  }

  out <- list(database = database, term2gene = t2g)
  cache$libraries[[cache_key]] <- out
  out
}

#' Preview terms and genes from one Enrichr database
#'
#' Downloads and previews terms from a selected Enrichr library to help choose
#' tissue/cell-type databases.
#'
#' @param database Enrichr library name.
#' @param pattern Backward-compatible alias for `tissue_pattern`.
#' @param n Backward-compatible alias for `max_terms`.
#' @param tissue_pattern Optional regular expression used to filter term names.
#' @param max_terms Maximum number of terms to return. Default is 80.
#' @param max_genes_per_term Maximum number of genes shown in `genes_preview`
#'   per term. Default is 25.
#' @param include_genes Logical. If `TRUE`, include `genes_preview`.
#'
#' @return A data frame with terms and gene counts.
#' @examples
#' hc_preview_celltype_database("PanglaoDB")
#' @export
hc_preview_celltype_database <- function(database,
                                         pattern = NULL,
                                         n = NULL,
                                         tissue_pattern = NULL,
                                         max_terms = 80,
                                         max_genes_per_term = 25,
                                         include_genes = TRUE) {
  if (!base::is.character(database) || base::length(database) != 1 || base::is.na(database) || !base::nzchar(database)) {
    stop("`database` must be a single non-empty character string.")
  }
  if (!base::is.null(pattern)) {
    if (!base::is.null(tissue_pattern)) {
      stop("Use either `pattern` or `tissue_pattern`, not both.")
    }
    tissue_pattern <- pattern
  }
  if (!base::is.null(n)) {
    if (!base::is.null(max_terms) && !identical(max_terms, 80)) {
      stop("Use either `n` or `max_terms`, not both.")
    }
    max_terms <- n
  }
  if (!base::is.numeric(max_terms) || base::length(max_terms) != 1 || base::is.na(max_terms) || max_terms < 1) {
    stop("`max_terms` must be a positive integer.")
  }
  if (!base::is.numeric(max_genes_per_term) || base::length(max_genes_per_term) != 1 || base::is.na(max_genes_per_term) || max_genes_per_term < 1) {
    stop("`max_genes_per_term` must be a positive integer.")
  }
  if (!base::is.logical(include_genes) || base::length(include_genes) != 1 || base::is.na(include_genes)) {
    stop("`include_genes` must be TRUE or FALSE.")
  }

  lib <- .hc_enrichr_fetch_library_term2gene(database = database, refresh = FALSE)
  t2g <- lib$term2gene
  if (base::nrow(t2g) == 0) {
    return(
      base::data.frame(
        database = base::character(0),
        term = base::character(0),
        n_genes = base::integer(0),
        stringsAsFactors = FALSE
      )
    )
  }

  by_term <- base::split(t2g$gene, t2g$term_raw)
  out <- base::data.frame(
    database = database,
    term = base::names(by_term),
    n_genes = base::as.integer(base::lengths(by_term)),
    stringsAsFactors = FALSE
  )

  if (!base::is.null(tissue_pattern) && base::nzchar(base::as.character(tissue_pattern[[1]]))) {
    keep <- base::grepl(
      pattern = base::as.character(tissue_pattern[[1]]),
      x = out$term,
      ignore.case = TRUE,
      perl = TRUE
    )
    out <- out[keep, , drop = FALSE]
    by_term <- by_term[out$term]
  }

  out <- out[base::order(-out$n_genes, out$term), , drop = FALSE]
  if (isTRUE(include_genes) && base::nrow(out) > 0) {
    preview <- base::vapply(
      out$term,
      function(tt) {
        genes <- base::as.character(by_term[[tt]])
        genes <- genes[!base::is.na(genes) & genes != ""]
        genes <- base::unique(genes)
        if (base::length(genes) > max_genes_per_term) {
          genes <- base::c(genes[base::seq_len(max_genes_per_term)], "...")
        }
        base::paste0(genes, collapse = ", ")
      },
      FUN.VALUE = base::character(1)
    )
    out$genes_preview <- preview
  }

  if (base::nrow(out) > max_terms) {
    out <- out[base::seq_len(max_terms), , drop = FALSE]
  }
  base::rownames(out) <- NULL
  out
}

.hc_default_coarse_celltype_map <- function() {
  c(
    "Monocyte/Macrophage" = "(monocyte|macrophage|myeloid|cd14\\+?\\s*mon|cd16\\+?\\s*mon|kupffer)",
    "Neutrophil" = "(neutrophil|pmn|granulocyte)",
    "Dendritic cell" = "(dendritic|dc\\b|plasmacytoid\\s*dendritic|cDC|pDC)",
    "T cell" = "(\\bt\\s*cell\\b|cd4\\+?\\s*t|cd8\\+?\\s*t|th1|th2|th17|treg|cytotoxic\\s*t)",
    "B cell/Plasma" = "(\\bb\\s*cell\\b|plasma\\s*cell|plasmablast|memory\\s*b)",
    "NK cell/ILC" = "(\\bnk\\b|natural\\s*killer|innate\\s*lymphoid|\\bilc\\b)",
    "Endothelial" = "(endothelial|vascular\\s*endotheli|venous|arterial\\s*endotheli)",
    "Epithelial" = "(epithelial|epithelium|keratinocyte)",
    "Fibroblast/Stromal" = "(fibroblast|stromal|mesenchymal)",
    "Smooth muscle/Pericyte" = "(smooth\\s*muscle|pericyte|myofibroblast)",
    "Hepatocyte" = "(hepatocyte|liver\\s*cell)",
    "Neuron/Glia" = "(neuron|astrocyte|oligodendro|microglia|glia)",
    "Other immune" = "(immune|leukocyte|lymphocyte)"
  )
}

.hc_assign_coarse_celltype <- function(term, coarse_map) {
  term <- base::as.character(term)
  if (base::length(term) == 0) {
    return(base::character(0))
  }
  out <- base::rep("Other", base::length(term))
  nms <- base::names(coarse_map)
  for (i in base::seq_along(coarse_map)) {
    patt <- coarse_map[[i]]
    lbl <- nms[[i]]
    hit <- base::grepl(patt, term, ignore.case = TRUE, perl = TRUE)
    out[hit & out == "Other"] <- lbl
  }
  out
}

.hc_extract_cluster_genes <- function(cluster_info, cluster_id) {
  genes <- dplyr::filter(cluster_info, color == cluster_id) %>%
    dplyr::pull(., "gene_n") %>%
    base::strsplit(., split = ",") %>%
    base::unlist(.)
  genes <- base::toupper(base::trimws(base::as.character(genes)))
  genes <- genes[!base::is.na(genes) & genes != ""]
  base::unique(genes)
}

.hc_find_output_dir <- function() {
  out_dir <- tryCatch(hcobject[["working_directory"]][["dir_output"]], error = function(e) NULL)
  save_folder <- tryCatch(hcobject[["global_settings"]][["save_folder"]], error = function(e) NULL)
  if (base::is.null(out_dir) || base::is.null(save_folder)) {
    return(NULL)
  }
  out_dir <- base::as.character(out_dir[[1]])
  save_folder <- base::as.character(save_folder[[1]])
  if (base::is.na(out_dir) || !base::nzchar(out_dir) || base::is.na(save_folder)) {
    return(NULL)
  }
  normalizePath(
    file.path(out_dir, save_folder),
    winslash = "/",
    mustWork = FALSE
  )
}

.hc_sanitize_db_token <- function(x) {
  tok <- base::gsub("[^A-Za-z0-9]+", "_", base::as.character(x), perl = TRUE)
  tok <- base::gsub("^_+|_+$", "", tok, perl = TRUE)
  tok <- base::tolower(tok)
  if (!base::nzchar(tok)) {
    tok <- "db"
  }
  tok
}

#' Automatic module cell-type annotation from Enrichr resources
#'
#' Performs module-wise enrichment against one or more Enrichr cell/tissue
#' libraries and stores the resulting composition in the same format that
#' `hc_plot_cluster_heatmap()` uses for right-side stacked-bar annotations.
#'
#' This enables a quick "module deconvolution"-style readout, e.g.
#' "module M4 consists of 60% monocyte-like and 40% neutrophil-like signal"
#' based on selected Enrichr term hits.
#'
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
#' @param ... Additional arguments passed to [hc_plot_cluster_heatmap()].
#'
#' @return Invisibly returns a list with selected/significant results and a
#'   `annotation_slot_map` (database -> slot).
#' @noRd
.hc_celltype_annotation_driver <- function(
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
  ...
) {

  mode <- base::match.arg(mode)
  annotation_slot <- base::match.arg(annotation_slot)

  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    stop("Package `clusterProfiler` is required for `.hc_celltype_annotation_driver()`.")
  }
  if (base::is.null(databases)) {
    databases <- base::character(0)
  }
  if (!base::is.character(databases)) {
    stop("`databases` must be a character vector.")
  }
  databases <- base::unique(base::trimws(base::as.character(databases)))
  databases <- databases[!base::is.na(databases) & databases != ""]
  if (!base::is.numeric(top) || base::length(top) != 1 || base::is.na(top) || top < 1) {
    stop("`top` must be a positive integer.")
  }
  top <- base::as.integer(top)
  if (!base::is.numeric(qval) || base::length(qval) != 1 || base::is.na(qval) || qval <= 0 || qval > 1) {
    stop("`qval` must be a numeric scalar in (0, 1].")
  }
  if (!base::is.numeric(min_term_genes) || base::length(min_term_genes) != 1 || base::is.na(min_term_genes) || min_term_genes < 1) {
    stop("`min_term_genes` must be a positive integer.")
  }
  min_term_genes <- base::as.integer(min_term_genes)
  if (!base::is.numeric(min_gs_size) || base::length(min_gs_size) != 1 || base::is.na(min_gs_size) || min_gs_size < 1) {
    stop("`min_gs_size` must be a positive integer.")
  }
  if (!base::is.numeric(max_gs_size) || base::length(max_gs_size) != 1 || base::is.na(max_gs_size) || max_gs_size < min_gs_size) {
    stop("`max_gs_size` must be >= `min_gs_size`.")
  }
  min_gs_size <- base::as.integer(min_gs_size)
  max_gs_size <- base::as.integer(max_gs_size)
  if (!base::is.logical(coarse_include_other) || base::length(coarse_include_other) != 1 || base::is.na(coarse_include_other)) {
    stop("`coarse_include_other` must be TRUE or FALSE.")
  }
  if (!base::is.logical(refresh_db) || base::length(refresh_db) != 1 || base::is.na(refresh_db)) {
    stop("`refresh_db` must be TRUE or FALSE.")
  }
  if (!base::is.logical(export_excel) || base::length(export_excel) != 1 || base::is.na(export_excel)) {
    stop("`export_excel` must be TRUE or FALSE.")
  }
  if (!base::is.logical(plot_heatmap) || base::length(plot_heatmap) != 1 || base::is.na(plot_heatmap)) {
    stop("`plot_heatmap` must be TRUE or FALSE.")
  }
  if (!base::is.null(slot_suffix)) {
    if (!base::is.character(slot_suffix) || base::length(slot_suffix) != 1) {
      stop("`slot_suffix` must be NULL or a single character string.")
    }
    slot_suffix <- base::trimws(slot_suffix)
    if (!base::nzchar(slot_suffix)) {
      stop("`slot_suffix` must be non-empty when provided.")
    }
    slot_suffix <- .hc_sanitize_db_token(slot_suffix)
    if (!base::nzchar(slot_suffix)) {
      stop("`slot_suffix` contains no usable characters after sanitization.")
    }
  }
  if (!base::is.logical(clear_previous_slots) ||
    base::length(clear_previous_slots) != 1 ||
    base::is.na(clear_previous_slots)) {
    stop("`clear_previous_slots` must be TRUE or FALSE.")
  }

  cluster_info <- tryCatch(hcobject[["integrated_output"]][["cluster_calc"]][["cluster_information"]], error = function(e) NULL)
  if (base::is.null(cluster_info) || base::nrow(cluster_info) == 0) {
    stop("No cluster information found. Run `hc_cluster_calculation()` first.")
  }
  if (!all(c("color", "gene_n") %in% base::colnames(cluster_info))) {
    stop("`cluster_information` is missing required columns (`color`, `gene_n`).")
  }

  all_clusters <- base::unique(base::as.character(cluster_info$color))
  all_clusters <- all_clusters[!base::is.na(all_clusters) & all_clusters != "white"]
  if (base::length(all_clusters) == 0) {
    stop("No non-white modules available.")
  }

  if (base::length(clusters) == 1 && identical(clusters[[1]], "all")) {
    clusters <- all_clusters
  } else {
    clusters <- base::as.character(clusters)
    missing_clusters <- base::setdiff(clusters, all_clusters)
    if (base::length(missing_clusters) > 0) {
      warning("Ignoring unknown cluster(s): ", base::paste(missing_clusters, collapse = ", "))
      clusters <- clusters[clusters %in% all_clusters]
    }
    if (base::length(clusters) == 0) {
      stop("No valid clusters selected.")
    }
  }

  if (mode == "coarse") {
    if (base::is.null(coarse_map)) {
      coarse_map <- .hc_default_coarse_celltype_map()
    }
    if (!base::is.character(coarse_map) || base::is.null(base::names(coarse_map)) || any(base::names(coarse_map) == "")) {
      stop("`coarse_map` must be a named character vector.")
    }
  }

  custom_t2g <- .hc_load_custom_gmt_term2gene(
    custom_gmt_files = custom_gmt_files,
    default_prefix = "CustomCellType",
    uppercase_genes = TRUE,
    add_database_prefix = TRUE,
    title_case_names = FALSE
  )
  custom_gmt_paths <- base::attr(custom_t2g, "paths")

  if (base::length(databases) > 0 && base::length(custom_t2g) > 0) {
    overlap_db <- base::intersect(databases, base::names(custom_t2g))
    if (base::length(overlap_db) > 0) {
      stop(
        "Overlapping Enrichr and custom GMT database names: ",
        base::paste(overlap_db, collapse = ", "),
        ". Rename custom GMT entries (names(custom_gmt_files))."
      )
    }
  }

  db_t2g <- list()
  if (base::length(databases) > 0) {
    .hc_enrichr_db_table(refresh = isTRUE(refresh_db))
    db_t2g <- lapply(
      databases,
      function(db) {
        .hc_enrichr_fetch_library_term2gene(
          database = db,
          refresh = isTRUE(refresh_db),
          min_genes_per_term = min_term_genes
        )$term2gene
      }
    )
    names(db_t2g) <- databases
  }
  if (base::length(custom_t2g) > 0) {
    db_t2g <- base::c(db_t2g, custom_t2g)
  }
  db_names <- base::names(db_t2g)
  if (base::length(db_names) == 0) {
    stop("No databases available. Provide `databases` and/or `custom_gmt_files`.")
  }

  combined_t2g <- base::do.call(base::rbind, db_t2g)
  if (base::is.null(combined_t2g) || base::nrow(combined_t2g) == 0) {
    stop("No usable TERM2GENE entries found for selected Enrichr/custom databases.")
  }
  combined_t2g <- base::unique(combined_t2g)
  term_meta <- base::unique(combined_t2g[, c("term", "term_raw", "database"), drop = FALSE])

  # Universe/background genes from all loaded count matrices.
  data_names <- tryCatch(base::names(hcobject[["data"]]), error = function(e) base::character(0))
  count_names <- data_names[base::grepl("_counts$", data_names)]
  universe <- base::unique(base::unlist(
    lapply(
      count_names,
      function(nm) {
        rn <- base::rownames(hcobject[["data"]][[nm]])
        base::toupper(base::trimws(base::as.character(rn)))
      }
    )
  ))
  universe <- universe[!base::is.na(universe) & universe != ""]
  if (base::length(universe) == 0) {
    universe <- base::unique(base::unlist(lapply(clusters, function(cl) .hc_extract_cluster_genes(cluster_info, cl))))
  }

  selected_by_db_nested <- stats::setNames(vector("list", base::length(db_names)), db_names)
  significant_by_db_nested <- stats::setNames(vector("list", base::length(db_names)), db_names)

  for (cl in clusters) {
    genes <- .hc_extract_cluster_genes(cluster_info, cl)
    if (base::length(genes) == 0) {
      next
    }
    enr <- tryCatch(
      clusterProfiler::enricher(
        genes,
        TERM2GENE = combined_t2g[, c("term", "gene"), drop = FALSE],
        pvalueCutoff = 1,
        qvalueCutoff = 1,
        pAdjustMethod = padj,
        universe = universe,
        minGSSize = min_gs_size,
        maxGSSize = max_gs_size
      ),
      error = function(e) NULL
    )
    if (base::is.null(enr) || base::is.null(enr@result) || base::nrow(enr@result) == 0) {
      next
    }

    sig <- enr@result
    # Storey's q-value is NA whenever the p-value distribution is too
    # degenerate to estimate pi0 (e.g. a single testable term), which would
    # silently drop every hit. Fall back to the `padj`-adjusted p-value there,
    # as `hc_functional_enrichment()` does.
    q_na <- base::is.na(.hc_as_numeric_safely(sig$qvalue))
    sig$qvalue[q_na] <- sig$p.adjust[q_na]
    sig <- sig[!base::is.na(sig$qvalue) & sig$qvalue <= qval, , drop = FALSE]
    if (base::nrow(sig) == 0) {
      next
    }
    term_col <- if ("Description" %in% base::colnames(sig)) "Description" else if ("ID" %in% base::colnames(sig)) "ID" else NULL
    if (base::is.null(term_col)) {
      next
    }

    sig$term <- base::as.character(sig[[term_col]])
    idx_meta <- base::match(sig$term, term_meta$term)
    sig$database <- term_meta$database[idx_meta]
    sig$term_raw <- term_meta$term_raw[idx_meta]
    sig$database[base::is.na(sig$database)] <- "Unknown"
    sig$term_raw[base::is.na(sig$term_raw)] <- sig$term[base::is.na(sig$term_raw)]
    sig$cluster <- cl
    sig$qvalue <- .hc_as_numeric_safely(sig$qvalue)
    sig$p.adjust <- .hc_as_numeric_safely(sig$p.adjust)
    sig$Count <- .hc_as_numeric_safely(sig$Count)
    sig$Count[!base::is.finite(sig$Count)] <- 0
    sig <- sig[base::order(sig$qvalue, sig$p.adjust, -sig$Count, sig$term), , drop = FALSE]

    for (db_nm in db_names) {
      sig_db <- sig[sig$database == db_nm, , drop = FALSE]
      if (base::nrow(sig_db) == 0) {
        next
      }

      sig_out_db <- sig_db[, c("cluster", "database", "term_raw", "term", "Count", "qvalue", "p.adjust", "geneID"), drop = FALSE]

      if (mode == "fine") {
        sel_db <- sig_db[base::seq_len(base::min(top, base::nrow(sig_db))), , drop = FALSE]
        sel_db$count <- sel_db$Count
        sel_db$cell_type <- sel_db$term
        sel_db$support_terms <- sel_db$term_raw
        sel_db <- sel_db[, c("cluster", "cell_type", "count", "qvalue", "p.adjust", "database", "support_terms"), drop = FALSE]
      } else {
        sig_db$class <- .hc_assign_coarse_celltype(sig_db$term_raw, coarse_map)
        if (!isTRUE(coarse_include_other)) {
          sig_db <- sig_db[sig_db$class != "Other", , drop = FALSE]
        }
        if (base::nrow(sig_db) == 0) {
          next
        }

        by_class <- base::split(base::seq_len(base::nrow(sig_db)), sig_db$class)
        class_rows <- lapply(base::names(by_class), function(cls) {
          idx <- by_class[[cls]]
          genes_hit <- base::unique(base::unlist(base::strsplit(base::paste(sig_db$geneID[idx], collapse = "/"), "/", fixed = TRUE)))
          genes_hit <- genes_hit[!base::is.na(genes_hit) & genes_hit != ""]
          base::data.frame(
            cluster = cl,
            cell_type = cls,
            count = base::length(genes_hit),
            qvalue = base::min(sig_db$qvalue[idx], na.rm = TRUE),
            p.adjust = base::min(sig_db$p.adjust[idx], na.rm = TRUE),
            database = db_nm,
            support_terms = base::paste(base::unique(sig_db$term_raw[idx]), collapse = " | "),
            stringsAsFactors = FALSE
          )
        })
        sel_db <- base::do.call(base::rbind, class_rows)
        if (base::is.null(sel_db) || base::nrow(sel_db) == 0) {
          next
        }
        sel_db <- sel_db[base::order(sel_db$qvalue, sel_db$p.adjust, -sel_db$count, sel_db$cell_type), , drop = FALSE]
        sel_db <- sel_db[base::seq_len(base::min(top, base::nrow(sel_db))), , drop = FALSE]
      }

      if (base::nrow(sel_db) == 0) {
        next
      }

      sel_db$rank <- base::seq_len(base::nrow(sel_db))
      hits <- base::sum(sel_db$count, na.rm = TRUE)
      sel_db$hits <- hits
      sel_db$pct <- if (hits > 0) {
        (sel_db$count / hits) * 100
      } else {
        0
      }

      selected_by_db_nested[[db_nm]][[cl]] <- sel_db
      significant_by_db_nested[[db_nm]][[cl]] <- sig_out_db
    }
  }

  selected_by_db <- lapply(
    selected_by_db_nested,
    function(x) {
      if (base::length(x) == 0) {
        return(NULL)
      }
      base::do.call(base::rbind, x)
    }
  )
  significant_by_db <- lapply(
    significant_by_db_nested,
    function(x) {
      if (base::length(x) == 0) {
        return(NULL)
      }
      base::do.call(base::rbind, x)
    }
  )

  non_empty_selected <- selected_by_db[!base::vapply(selected_by_db, base::is.null, FUN.VALUE = logical(1))]
  if (base::length(non_empty_selected) == 0) {
    stop("No significant cell-type terms found for selected modules and settings.")
  }

  selected_all <- base::do.call(base::rbind, non_empty_selected)
  selected_all <- selected_all[, c("cluster", "cell_type", "count", "pct", "hits", "rank", "qvalue", "p.adjust", "database", "support_terms"), drop = FALSE]

  non_empty_significant <- significant_by_db[!base::vapply(significant_by_db, base::is.null, FUN.VALUE = logical(1))]
  significant_all <- if (base::length(non_empty_significant) > 0) {
    base::do.call(base::rbind, non_empty_significant)
  } else {
    base::data.frame()
  }

  celltype_stats <- stats::aggregate(
    cbind(count, qvalue) ~ database + cell_type,
    data = selected_all[, c("database", "cell_type", "count", "qvalue"), drop = FALSE],
    FUN = function(x) c(sum = base::sum(x, na.rm = TRUE), min = base::min(x, na.rm = TRUE))
  )
  celltype_stats <- base::data.frame(
    database = celltype_stats$database,
    cell_type = celltype_stats$cell_type,
    total_count = celltype_stats$count[, "sum"],
    best_q = celltype_stats$qvalue[, "min"],
    stringsAsFactors = FALSE
  )
  n_clusters_by_type <- stats::aggregate(cluster ~ database + cell_type, data = selected_all, FUN = function(x) base::length(base::unique(x)))
  base::colnames(n_clusters_by_type)[3] <- "n_clusters"
  celltype_stats <- dplyr::left_join(celltype_stats, n_clusters_by_type, by = c("database", "cell_type"))
  celltype_stats <- celltype_stats[base::order(celltype_stats$database, celltype_stats$best_q, -celltype_stats$n_clusters, -celltype_stats$total_count, celltype_stats$cell_type), , drop = FALSE]

  categories_by_db <- list()
  celltype_order_by_db <- list()
  for (db_nm in db_names) {
    sel_db <- selected_by_db[[db_nm]]
    if (base::is.null(sel_db) || base::nrow(sel_db) == 0) {
      categories_by_db[[db_nm]] <- base::data.frame(
        count = base::numeric(0),
        cell_type = base::character(0),
        cluster = base::character(0),
        hits = base::numeric(0),
        stringsAsFactors = FALSE
      )
      celltype_order_by_db[[db_nm]] <- base::data.frame(
        cell_type = base::character(0),
        total_count = base::numeric(0),
        best_q = base::numeric(0),
        n_clusters = base::integer(0),
        stringsAsFactors = FALSE
      )
      next
    }

    celltype_stats_db <- stats::aggregate(
      cbind(count, qvalue) ~ cell_type,
      data = sel_db[, c("cell_type", "count", "qvalue"), drop = FALSE],
      FUN = function(x) c(sum = base::sum(x, na.rm = TRUE), min = base::min(x, na.rm = TRUE))
    )
    celltype_stats_db <- base::data.frame(
      cell_type = celltype_stats_db$cell_type,
      total_count = celltype_stats_db$count[, "sum"],
      best_q = celltype_stats_db$qvalue[, "min"],
      stringsAsFactors = FALSE
    )
    n_clusters_db <- stats::aggregate(cluster ~ cell_type, data = sel_db, FUN = function(x) base::length(base::unique(x)))
    base::colnames(n_clusters_db)[2] <- "n_clusters"
    celltype_stats_db <- dplyr::left_join(celltype_stats_db, n_clusters_db, by = "cell_type")
    celltype_stats_db <- celltype_stats_db[base::order(celltype_stats_db$best_q, -celltype_stats_db$n_clusters, -celltype_stats_db$total_count, celltype_stats_db$cell_type), , drop = FALSE]
    cell_types_db <- celltype_stats_db$cell_type

    grid_db <- base::expand.grid(
      cluster = clusters,
      cell_type = cell_types_db,
      stringsAsFactors = FALSE
    )
    categories_db <- dplyr::left_join(
      grid_db,
      sel_db[, c("cluster", "cell_type", "pct"), drop = FALSE],
      by = c("cluster", "cell_type")
    )
    categories_db$pct[base::is.na(categories_db$pct)] <- 0
    hits_map_db <- stats::setNames(
      tapply(sel_db$count, sel_db$cluster, base::sum),
      base::names(tapply(sel_db$count, sel_db$cluster, base::sum))
    )
    categories_db$hits <- hits_map_db[categories_db$cluster]
    categories_db$hits[base::is.na(categories_db$hits)] <- 0
    categories_db$count <- categories_db$pct
    categories_db <- categories_db[, c("count", "cell_type", "cluster", "hits"), drop = FALSE]
    categories_db$cluster <- base::factor(categories_db$cluster, levels = clusters)
    categories_db$cell_type <- base::factor(categories_db$cell_type, levels = cell_types_db)
    categories_db <- categories_db[base::order(categories_db$cluster, categories_db$cell_type), , drop = FALSE]
    categories_db$cluster <- base::as.character(categories_db$cluster)
    categories_db$cell_type <- base::as.character(categories_db$cell_type)
    base::rownames(categories_db) <- NULL

    categories_by_db[[db_nm]] <- categories_db
    celltype_order_by_db[[db_nm]] <- celltype_stats_db
  }

  categories <- base::do.call(
    base::rbind,
    lapply(
      base::names(categories_by_db),
      function(db_nm) {
        x <- categories_by_db[[db_nm]]
        if (base::nrow(x) == 0) {
          return(NULL)
        }
        x$database <- db_nm
        x <- x[, c("database", "count", "cell_type", "cluster", "hits"), drop = FALSE]
        x
      }
    )
  )
  if (base::is.null(categories)) {
    categories <- base::data.frame(
      database = base::character(0),
      count = base::numeric(0),
      cell_type = base::character(0),
      cluster = base::character(0),
      hits = base::numeric(0),
      stringsAsFactors = FALSE
    )
  }

  db_tokens <- base::make.unique(
    base::vapply(db_names, .hc_sanitize_db_token, FUN.VALUE = base::character(1)),
    sep = "_"
  )
  db_slot_names <- base::paste0("enriched_per_cluster_", db_tokens)
  if (!base::is.null(slot_suffix) && base::nzchar(slot_suffix)) {
    db_slot_names <- base::paste0(db_slot_names, "_", slot_suffix)
  }
  base::names(db_slot_names) <- db_names
  if (base::length(db_names) == 1 && annotation_slot != "auto") {
    db_slot_names[[1]] <- annotation_slot
    if (!base::is.null(slot_suffix) && base::nzchar(slot_suffix)) {
      db_slot_names[[1]] <- base::paste0(db_slot_names[[1]], "_", slot_suffix)
    }
  }

  sat_out <- hcobject[["satellite_outputs"]]
  if (base::is.null(sat_out) || !base::is.list(sat_out)) {
    sat_out <- list()
  }
  if (isTRUE(clear_previous_slots)) {
    sat_names <- base::names(sat_out)
    old_dynamic_slots <- sat_names[base::grepl("^enriched_per_cluster_", sat_names)]
    for (nm in old_dynamic_slots) {
      sat_out[[nm]] <- NULL
    }
    # Reset legacy slots as well to avoid stale mixed annotations after reruns.
    for (nm in c("enriched_per_cluster", "enriched_per_cluster2")) {
      if (nm %in% sat_names) {
        sat_out[[nm]] <- NULL
      }
    }
  }

  hidden_common <- list(
    source = "celltype_annotation",
    databases = db_names,
    custom_gmt_files = custom_gmt_paths,
    clusters = clusters,
    mode = mode,
    top = top,
    qval = qval,
    padj = padj,
    slot_suffix = slot_suffix,
    clear_previous_slots = clear_previous_slots
  )
  base::attr(categories, "hidden") <- hidden_common

  for (db_nm in db_names) {
    slot_nm <- db_slot_names[[db_nm]]
    if (base::is.na(slot_nm) || !base::nzchar(slot_nm)) {
      slot_nm <- base::paste0("enriched_per_cluster_", .hc_sanitize_db_token(db_nm))
    }
    categories_db <- categories_by_db[[db_nm]]
    if (base::is.null(categories_db)) {
      categories_db <- base::data.frame(
        count = base::numeric(0),
        cell_type = base::character(0),
        cluster = base::character(0),
        hits = base::numeric(0),
        stringsAsFactors = FALSE
      )
    }
    hidden_db <- base::c(hidden_common, list(database = db_nm, label = db_nm, slot = slot_nm))
    base::attr(categories_db, "hidden") <- hidden_db
    slot_payload <- list(categories_per_cluster = categories_db)
    base::attr(slot_payload[["categories_per_cluster"]], "hidden") <- hidden_db
    sat_out[[slot_nm]] <- slot_payload
  }

  out <- list(
    categories_per_cluster = categories,
    categories_per_cluster_by_db = categories_by_db,
    selected_celltypes = selected_all,
    selected_celltypes_by_db = selected_by_db,
    significant_terms = significant_all,
    significant_terms_by_db = significant_by_db,
    celltype_order = celltype_stats,
    celltype_order_by_db = celltype_order_by_db,
    annotation_slots = base::unname(db_slot_names),
    annotation_slot_map = base::data.frame(
      database = db_names,
      slot = base::unname(db_slot_names),
      stringsAsFactors = FALSE
    ),
    parameters = list(
      databases = db_names,
      clusters = clusters,
      mode = mode,
      top = top,
      qval = qval,
      padj = padj,
      min_term_genes = min_term_genes,
      min_gs_size = min_gs_size,
      max_gs_size = max_gs_size,
      custom_gmt_files = custom_gmt_paths,
      annotation_slot = annotation_slot,
      annotation_slots = base::unname(db_slot_names),
      slot_suffix = slot_suffix,
      clear_previous_slots = clear_previous_slots,
      coarse_include_other = coarse_include_other
    )
  )
  base::attr(out$categories_per_cluster, "hidden") <- hidden_common

  sat_out[["celltype_annotation"]] <- out
  .hc_set_bridge_hcobject_slot("satellite_outputs", sat_out)

  if (isTRUE(export_excel)) {
    out_dir <- .hc_find_output_dir()
    if (!base::is.null(out_dir)) {
      if (!dir.exists(out_dir)) {
        dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      }
      excel_path <- file.path(out_dir, excel_file)
      make_sheet_name <- function(prefix, db_name, used_names) {
        base_name <- base::paste0(prefix, "_", .hc_sanitize_db_token(db_name))
        base_name <- base::substr(base_name, 1, 31)
        candidate <- base_name
        k <- 1
        while (candidate %in% used_names) {
          suffix <- base::paste0("_", k)
          candidate <- base::paste0(
            base::substr(base_name, 1, base::max(1, 31 - base::nchar(suffix))),
            suffix
          )
          k <- k + 1
        }
        candidate
      }
      xlsx_tables <- list(
        selected_celltypes = selected_all,
        significant_terms = significant_all,
        categories_per_cluster = categories,
        celltype_order = celltype_stats
      )
      used_sheet_names <- base::names(xlsx_tables)
      for (db_nm in db_names) {
        sel_db <- selected_by_db[[db_nm]]
        sig_db <- significant_by_db[[db_nm]]
        cat_db <- categories_by_db[[db_nm]]
        ord_db <- celltype_order_by_db[[db_nm]]
        if (base::is.null(sel_db)) {
          sel_db <- base::data.frame()
        }
        if (base::is.null(sig_db)) {
          sig_db <- base::data.frame()
        }
        if (base::is.null(cat_db)) {
          cat_db <- base::data.frame()
        }
        if (base::is.null(ord_db)) {
          ord_db <- base::data.frame()
        }
        nm_sel <- make_sheet_name("selected", db_nm, used_sheet_names)
        xlsx_tables[[nm_sel]] <- sel_db
        used_sheet_names <- base::c(used_sheet_names, nm_sel)
        nm_sig <- make_sheet_name("signif", db_nm, used_sheet_names)
        xlsx_tables[[nm_sig]] <- sig_db
        used_sheet_names <- base::c(used_sheet_names, nm_sig)
        nm_cat <- make_sheet_name("cat", db_nm, used_sheet_names)
        xlsx_tables[[nm_cat]] <- cat_db
        used_sheet_names <- base::c(used_sheet_names, nm_cat)
        nm_ord <- make_sheet_name("order", db_nm, used_sheet_names)
        xlsx_tables[[nm_ord]] <- ord_db
        used_sheet_names <- base::c(used_sheet_names, nm_ord)
      }
      tryCatch(
        .hc_write_xlsx_atomic(
          x = xlsx_tables,
          file = excel_path,
          overwrite = TRUE
        ),
        error = function(e) warning("Could not write cell-type annotation Excel: ", conditionMessage(e))
      )
    } else {
      warning("Output folder not configured; skipping Excel export.")
    }
  }

  if (isTRUE(plot_heatmap)) {
    hm_args <- list(
      file_name = heatmap_file_name,
      gene_count_mode = "none",
      module_label_preset = "compact",
      include_dynamic_enrichment_slots = TRUE,
      celltype_bar_show_dominant = TRUE,
      celltype_bar_top_n = top,
      celltype_bar_include_other = TRUE
    )

    cluster_calc <- tryCatch(
      hcobject[["integrated_output"]][["cluster_calc"]],
      error = function(e) NULL
    )
    if (base::is.list(cluster_calc) && base::length(cluster_calc) > 0) {
      if ("module_label_mode" %in% base::names(cluster_calc) &&
        !base::is.null(cluster_calc[["module_label_mode"]])) {
        hm_args$module_label_mode <- cluster_calc[["module_label_mode"]]
      }
      if ("module_label_numbering" %in% base::names(cluster_calc) &&
        !base::is.null(cluster_calc[["module_label_numbering"]])) {
        hm_args$module_label_numbering <- cluster_calc[["module_label_numbering"]]
      }
      if ("module_prefix" %in% base::names(cluster_calc) &&
        !base::is.null(cluster_calc[["module_prefix"]])) {
        hm_args$module_prefix <- cluster_calc[["module_prefix"]]
      }
      if ("module_label_fontsize" %in% base::names(cluster_calc) &&
        base::is.numeric(cluster_calc[["module_label_fontsize"]]) &&
        base::length(cluster_calc[["module_label_fontsize"]]) == 1 &&
        base::is.finite(cluster_calc[["module_label_fontsize"]]) &&
        cluster_calc[["module_label_fontsize"]] > 0) {
        hm_args$module_label_fontsize <- cluster_calc[["module_label_fontsize"]]
      }
      if ("module_label_pt_size" %in% base::names(cluster_calc) &&
        base::is.numeric(cluster_calc[["module_label_pt_size"]]) &&
        base::length(cluster_calc[["module_label_pt_size"]]) == 1 &&
        base::is.finite(cluster_calc[["module_label_pt_size"]]) &&
        cluster_calc[["module_label_pt_size"]] > 0) {
        hm_args$module_label_pt_size <- cluster_calc[["module_label_pt_size"]]
      }
      if ("module_box_width_cm" %in% base::names(cluster_calc) &&
        base::is.numeric(cluster_calc[["module_box_width_cm"]]) &&
        base::length(cluster_calc[["module_box_width_cm"]]) == 1 &&
        base::is.finite(cluster_calc[["module_box_width_cm"]]) &&
        cluster_calc[["module_box_width_cm"]] > 0) {
        hm_args$module_box_width_cm <- cluster_calc[["module_box_width_cm"]]
      }
    }

    extra <- list(...)
    if (base::length(extra) > 0) {
      hm_args <- utils::modifyList(hm_args, extra)
    }
    base::do.call(.hc_plot_cluster_heatmap_driver, hm_args)
  }

  invisible(out)
}


