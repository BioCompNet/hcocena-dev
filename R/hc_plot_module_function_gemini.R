#' Plot AI-assisted module function summaries
#'
#' Draws a compact overview of stored AI-assisted module function summaries with
#' a colored module box on the left and the inferred overarching function on the
#' right.
#'
#' @param hc An `HCoCenaExperiment` with stored results from
#'   `hc_module_function_llm()`.
#' @param slot_name Satellite slot name used for storage. Default is
#'   `"llm_module_function"`.
#' @param modules Optional character vector to subset modules.
#' @param fields Character vector selecting which LLM fields to plot. Supported
#'   values are `"general_processes"`, `"contextual_state"`,
#'   `"key_regulators"`, `"rag_general_processes"`,
#'   `"rag_contextual_state"`, and `"rag_key_regulators"`. The aliases
#'   `"general_processes_rag"`, `"contextual_state_rag"`, and
#'   `"key_regulators_rag"` are also accepted. Defaults to the three baseline
#'   fields.
#' @param max_chars Maximum number of characters shown per term. Default is
#'   `90`.
#' @param text_size Numeric text size passed to `ggplot2::geom_text()`.
#' @param with_heatmap Logical. If `TRUE`, reuse the stored hCoCena cluster
#'   heatmap and place the LLM interpretations as a right-side annotation.
#'   Falls back to a text-only plot if no heatmap cache is available.
#' @param col_order Optional character vector overriding the hCoCena
#'   heatmap column order for this LLM plot only. If `NULL` (default), the
#'   column order from the main module heatmap is reused when available.
#' @param heatmap_col_order Legacy alias for `col_order`.
#' @param cluster_columns Logical. If `FALSE` (default), reuse the
#'   column order from the main hCoCena heatmap when available. If `TRUE`,
#'   cluster the columns for this LLM plot instead.
#' @param heatmap_cluster_columns Legacy alias for `cluster_columns`.
#' @param heatmap_rel_width Relative width of the heatmap panel when
#'   `with_heatmap = TRUE`.
#' @param text_rel_width Relative width of the text panel when
#'   `with_heatmap = TRUE`.
#' @param title Plot title. Default is `"AI-assisted module interpretation"`.
#' @param module_label_fontsize Optional numeric fontsize for module labels
#'   inside the colored module boxes. If `NULL`, a compact
#'   `hc_plot_cluster_heatmap()`-like default is used.
#' @param module_label_pt_size Optional numeric glyph size for module labels
#'   inside the colored module boxes. Uses `snpc` units. If `NULL`, a compact
#'   `hc_plot_cluster_heatmap()`-like default is used.
#' @param module_box_width_cm Optional numeric width (cm) for the module color
#'   boxes. If `NULL`, a compact `hc_plot_cluster_heatmap()`-like default is
#'   used.
#' @param save Logical. If `TRUE` (default), export each selected field plot as
#'   both PDF and PNG into the configured hCoCena output directory.
#' @param file_stem Base file name used for exported plots. Field names and, for
#'   single-module selections, the module label are appended automatically.
#'   Default is `"LLM_module_function"`.
#' @param pdf_width Optional numeric export width in inches. If `NULL`, a
#'   field-appropriate default is chosen.
#' @param pdf_height Optional numeric export height in inches. If `NULL`, a
#'   module-count dependent default is chosen.
#' @param dpi Numeric export DPI for PNG output. Default is `300`.
#' @param ... Used by the backward-compatible wrapper alias
#'   `hc_plot_module_function_gemini()`.
#'
#' @return A single plot object when one field is selected, otherwise a named
#'   list of plot objects.
#' @export
#'
#' @examples
#' # A stand-in for the result of hc_module_function_llm(), which needs an LLM
#' # endpoint:
#' hc <- hc_init()
#' methods::slot(hc, "satellite")[["llm_module_function"]] <- list(
#'   module_1 = list(status = "ok")
#' )
#' methods::slot(hc, "satellite")[["llm_module_function_summary"]] <- data.frame(
#'   module = "module_1",
#'   module_color = "steelblue",
#'   general_processes = "Interferon signaling",
#'   contextual_state = "Acute antiviral activation",
#'   key_regulators = "STAT1, IRF7",
#'   stringsAsFactors = FALSE
#' )
#' p <- hc_plot_module_function_llm(hc, with_heatmap = FALSE, save = FALSE)
#' print(p)
hc_plot_module_function_llm <- function(hc,
                                        slot_name = "llm_module_function",
                                        modules = NULL,
                                        fields = c("general_processes", "contextual_state", "key_regulators"),
                                        max_chars = 90,
                                        text_size = 4,
                                        with_heatmap = TRUE,
                                        col_order = NULL,
                                        heatmap_col_order = NULL,
                                        cluster_columns = FALSE,
                                        heatmap_cluster_columns = NULL,
                                        heatmap_rel_width = 1.7,
                                        text_rel_width = 1.25,
                                        title = "AI-assisted module interpretation",
                                        module_label_fontsize = NULL,
                                        module_label_pt_size = NULL,
                                        module_box_width_cm = NULL,
                                        save = TRUE,
                                        file_stem = "LLM_module_function",
                                        pdf_width = NULL,
                                        pdf_height = NULL,
                                        dpi = 300) {
  if (missing(hc) || is.null(hc)) {
    stop("`hc` must be provided.")
  }
  if (!is.numeric(max_chars) || length(max_chars) != 1 || is.na(max_chars) || max_chars < 20) {
    stop("`max_chars` must be a single number >= 20.")
  }
  if (!is.numeric(text_size) || length(text_size) != 1 || is.na(text_size) || text_size <= 0) {
    stop("`text_size` must be a single positive number.")
  }
  if (!is.logical(with_heatmap) || length(with_heatmap) != 1 || is.na(with_heatmap)) {
    stop("`with_heatmap` must be TRUE or FALSE.")
  }
  col_order <- .hc_resolve_col_order_alias(
    col_order = col_order,
    heatmap_col_order = heatmap_col_order,
    col_order_missing = missing(col_order),
    heatmap_col_order_missing = missing(heatmap_col_order),
    context = "hc_plot_module_function_llm()"
  )
  cluster_columns <- .hc_resolve_cluster_columns_alias(
    cluster_columns = cluster_columns,
    heatmap_cluster_columns = heatmap_cluster_columns,
    cluster_columns_missing = missing(cluster_columns),
    heatmap_cluster_columns_missing = missing(heatmap_cluster_columns),
    context = "hc_plot_module_function_llm()"
  )
  if (!is.logical(cluster_columns) || length(cluster_columns) != 1 || is.na(cluster_columns)) {
    stop("`cluster_columns` must be TRUE or FALSE.")
  }
  if (!is.numeric(heatmap_rel_width) || length(heatmap_rel_width) != 1 || is.na(heatmap_rel_width) || heatmap_rel_width <= 0) {
    stop("`heatmap_rel_width` must be a single positive number.")
  }
  if (!is.numeric(text_rel_width) || length(text_rel_width) != 1 || is.na(text_rel_width) || text_rel_width <= 0) {
    stop("`text_rel_width` must be a single positive number.")
  }
  if (!is.null(module_label_fontsize) &&
    (!is.numeric(module_label_fontsize) || length(module_label_fontsize) != 1 || is.na(module_label_fontsize) || module_label_fontsize <= 0)) {
    stop("`module_label_fontsize` must be NULL or a single positive number.")
  }
  if (!is.null(module_label_pt_size) &&
    (!is.numeric(module_label_pt_size) || length(module_label_pt_size) != 1 || is.na(module_label_pt_size) || module_label_pt_size <= 0)) {
    stop("`module_label_pt_size` must be NULL or a single positive number.")
  }
  if (!is.null(module_box_width_cm) &&
    (!is.numeric(module_box_width_cm) || length(module_box_width_cm) != 1 || is.na(module_box_width_cm) || module_box_width_cm <= 0)) {
    stop("`module_box_width_cm` must be NULL or a single positive number.")
  }
  if (!is.logical(save) || length(save) != 1 || is.na(save)) {
    stop("`save` must be TRUE or FALSE.")
  }
  if (!is.character(file_stem) || length(file_stem) != 1 || !base::nzchar(file_stem)) {
    stop("`file_stem` must be a single non-empty character string.")
  }
  if (!is.null(pdf_width) && (!is.numeric(pdf_width) || length(pdf_width) != 1 || is.na(pdf_width) || pdf_width <= 0)) {
    stop("`pdf_width` must be NULL or a single positive number.")
  }
  if (!is.null(pdf_height) && (!is.numeric(pdf_height) || length(pdf_height) != 1 || is.na(pdf_height) || pdf_height <= 0)) {
    stop("`pdf_height` must be NULL or a single positive number.")
  }
  if (!is.numeric(dpi) || length(dpi) != 1 || is.na(dpi) || dpi <= 0) {
    stop("`dpi` must be a single positive number.")
  }
  field_map <- .hc_llm_plot_field_map()
  fields <- unique(as.character(fields))
  fields <- fields[!is.na(fields) & fields != ""]
  valid_fields <- names(field_map$source)
  if (length(fields) == 0 || !all(fields %in% valid_fields)) {
    stop("`fields` must contain one or more of: ", paste(valid_fields, collapse = ", "), ".")
  }

  sat <- tryCatch(base::as.list(hc@satellite), error = function(e) list())
  stored_results <- sat[[slot_name]]
  if (is.null(stored_results) || !is.list(stored_results) || length(stored_results) == 0) {
    stop(
      "No stored module interpretations found in `hc@satellite$", slot_name,
      "`. Run `hc_module_function_llm()` first."
    )
  }

  summary_tbl <- sat[[paste0(slot_name, "_summary")]]
  if (is.null(summary_tbl) || !is.data.frame(summary_tbl) || nrow(summary_tbl) == 0) {
    summary_tbl <- .hc_llm_summary_from_results(results = stored_results, hc = hc)
  }
  if (nrow(summary_tbl) == 0) {
    stop("No module interpretation summaries are available for plotting.")
  }

  if (!is.null(modules)) {
    modules <- as.character(modules)
    summary_tbl <- summary_tbl[summary_tbl$module %in% modules, , drop = FALSE]
  }
  if (nrow(summary_tbl) == 0) {
    stop("No matching modules found for plotting.")
  }

  heatmap_info <- if (isTRUE(with_heatmap)) .hc_llm_heatmap_info(hc) else NULL
  summary_tbl <- .hc_llm_reorder_summary_for_display(summary_tbl = summary_tbl, heatmap_info = heatmap_info)
  summary_tbl$module <- as.character(summary_tbl$module)
  summary_tbl$module_color <- as.character(summary_tbl$module_color)
  summary_tbl$module_color[is.na(summary_tbl$module_color) | summary_tbl$module_color == ""] <- "grey70"
  needed_summary_fields <- unique(unname(field_map$source[fields]))
  for (nm in needed_summary_fields) {
    if (!nm %in% colnames(summary_tbl)) {
      summary_tbl[[nm]] <- NA_character_
    }
    summary_tbl[[nm]] <- as.character(summary_tbl[[nm]])
  }
  summary_tbl$text_color <- vapply(summary_tbl$module_color, .hc_llm_darken_color, FUN.VALUE = character(1))
  summary_tbl$label_color <- vapply(summary_tbl$module_color, .hc_llm_contrast_text_color, FUN.VALUE = character(1))
  summary_tbl$module_factor <- factor(summary_tbl$module, levels = rev(summary_tbl$module))

  out <- lapply(fields, function(field_nm) {
    field_tbl <- summary_tbl
    source_nm <- unname(field_map$source[[field_nm]])
    field_tbl$term_plot <- vapply(
      ifelse(
        is.na(field_tbl[[source_nm]]) | field_tbl[[source_nm]] == "",
        "No interpretation available.",
        field_tbl[[source_nm]]
      ),
      .hc_llm_prepare_display_title,
      FUN.VALUE = character(1),
      max_chars = as.integer(max_chars[[1]])
    )

    plot_title <- paste0(title, ": ", field_map$title[[field_nm]])

    if (isTRUE(with_heatmap) && !is.null(heatmap_info) && isTRUE(heatmap_info$draw_supported)) {
      combined_grob <- .hc_llm_capture_combined_heatmap_grob(
        heatmap_info = heatmap_info,
        summary_tbl = field_tbl,
        max_chars = max_chars,
        text_size = text_size,
        module_label_fontsize = module_label_fontsize,
        module_label_pt_size = module_label_pt_size,
        module_box_width_cm = module_box_width_cm,
        heatmap_col_order = col_order,
        heatmap_cluster_columns = cluster_columns
      )
      return(.hc_llm_add_title_grob(combined_grob, title = plot_title, text_size = text_size))
    }

    ggplot2::ggplot(field_tbl, ggplot2::aes(y = module_factor)) +
      ggplot2::geom_tile(
        ggplot2::aes(x = 1, fill = module_color),
        width = 0.34,
        height = 0.78,
        show.legend = FALSE
      ) +
      ggplot2::geom_text(
        ggplot2::aes(x = 1, label = module, color = label_color),
        fontface = "bold",
        size = text_size * 0.85,
        show.legend = FALSE
      ) +
      ggplot2::geom_segment(
        ggplot2::aes(x = 1.22, xend = 1.34, yend = module_factor, color = text_color),
        linewidth = 0.6,
        show.legend = FALSE
      ) +
      ggplot2::geom_text(
        ggplot2::aes(x = 1.38, label = term_plot, color = text_color),
        hjust = 0,
        size = text_size,
        show.legend = FALSE
      ) +
      ggplot2::scale_fill_identity() +
      ggplot2::scale_color_identity() +
      ggplot2::coord_cartesian(xlim = c(0.75, 3.6), clip = "off") +
      ggplot2::scale_y_discrete(expand = ggplot2::expansion(mult = c(0.02, 0.02))) +
      ggplot2::labs(title = plot_title) +
      ggplot2::theme_void() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", hjust = 0, size = text_size * 4),
        plot.margin = ggplot2::margin(t = 10, r = 180, b = 10, l = 4)
      )
  })
  names(out) <- fields

  if (isTRUE(save) && length(out) > 0) {
    output_dir <- .hc_llm_plot_output_dir(hc)
    export_width <- if (is.null(pdf_width)) {
      if (isTRUE(with_heatmap)) 14 else 11
    } else {
      as.numeric(pdf_width[[1]])
    }
    export_height <- if (is.null(pdf_height)) {
      max(4.5, min(18, 2.4 + (0.46 * nrow(summary_tbl))))
    } else {
      as.numeric(pdf_height[[1]])
    }
    module_token <- .hc_llm_plot_module_token(summary_tbl)
    export_files <- lapply(names(out), function(field_nm) {
      export_file <- base::file.path(
        output_dir,
        base::paste0(
          .hc_export_sanitize_stem(file_stem, default = "LLM_module_function"),
          "_",
          module_token,
          "_",
          .hc_export_sanitize_stem(field_nm, default = "field"),
          ".pdf"
        )
      )
      .hc_export_single_page_plot(
        file = export_file,
        width = export_width,
        height = export_height,
        pointsize = 11,
        res = dpi,
        draw_fun = function() {
          .hc_display_object(out[[field_nm]])
        }
      )
    })
    names(export_files) <- names(out)
    attr(out, "output_files") <- export_files
  }

  if (length(out) == 1) {
    single_plot <- out[[1]]
    single_exports <- attr(out, "output_files", exact = TRUE)
    if (!is.null(single_exports)) {
      attr(single_plot, "output_files") <- single_exports[[1]]
    }
    return(single_plot)
  }
  out
}

#' @rdname hc_plot_module_function_llm
#' @export
hc_plot_module_function_gemini <- function(...) {
  hc_plot_module_function_llm(...)
}

.hc_llm_plot_field_map <- function() {
  source <- c(
    general_processes = "general_processes",
    contextual_state = "contextual_state",
    key_regulators = "key_regulators",
    rag_general_processes = "rag_general_processes",
    rag_contextual_state = "rag_contextual_state",
    rag_key_regulators = "rag_key_regulators",
    general_processes_rag = "rag_general_processes",
    contextual_state_rag = "rag_contextual_state",
    key_regulators_rag = "rag_key_regulators",
    enrichment_general_processes = "enrichment_general_processes",
    enrichment_contextual_state = "enrichment_contextual_state",
    enrichment_key_regulators = "enrichment_key_regulators",
    general_processes_enrichment = "enrichment_general_processes",
    contextual_state_enrichment = "enrichment_contextual_state",
    key_regulators_enrichment = "enrichment_key_regulators"
  )
  title <- c(
    general_processes = "General processes",
    contextual_state = "Contextual state",
    key_regulators = "Key regulators",
    rag_general_processes = "General processes with RAG",
    rag_contextual_state = "Contextual state with RAG",
    rag_key_regulators = "Key regulators with RAG",
    general_processes_rag = "General processes with RAG",
    contextual_state_rag = "Contextual state with RAG",
    key_regulators_rag = "Key regulators with RAG",
    enrichment_general_processes = "General processes with enrichment",
    enrichment_contextual_state = "Contextual state with enrichment",
    enrichment_key_regulators = "Key regulators with enrichment",
    general_processes_enrichment = "General processes with enrichment",
    contextual_state_enrichment = "Contextual state with enrichment",
    key_regulators_enrichment = "Key regulators with enrichment"
  )
  list(source = source, title = title)
}

.hc_llm_plot_output_dir <- function(hc) {
  paths <- tryCatch(.hc_row_to_list(hc@config@paths), error = function(e) list())
  global_cfg <- tryCatch(.hc_row_to_list(hc@config@global), error = function(e) list())
  out_dir <- paths[["dir_output"]]
  if (is.null(out_dir) || length(out_dir) == 0) {
    stop("No output directory configured in `hc@config@paths$dir_output`.")
  }
  out_dir <- as.character(out_dir[[1]])
  if (is.na(out_dir) || !base::nzchar(out_dir)) {
    stop("No output directory configured in `hc@config@paths$dir_output`.")
  }
  save_folder <- global_cfg[["save_folder"]]
  save_folder <- if (is.null(save_folder) || length(save_folder) == 0) {
    ""
  } else {
    as.character(save_folder[[1]])
  }
  if (is.na(save_folder) || save_folder %in% c("", "FALSE", "false")) {
    target_dir <- out_dir
  } else {
    target_dir <- base::file.path(out_dir, save_folder)
  }
  if (!base::dir.exists(target_dir)) {
    base::dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
  }
  target_dir
}

.hc_llm_plot_module_token <- function(summary_tbl) {
  mods <- unique(as.character(summary_tbl$module))
  mods <- mods[!is.na(mods) & mods != ""]
  if (length(mods) == 1) {
    return(.hc_export_sanitize_stem(mods[[1]], default = "module"))
  }
  if (length(mods) == 0) {
    return("modules")
  }
  base::paste0(length(mods), "_modules")
}

.hc_llm_title_wrap_width <- function(available_width_in = NA_real_) {
  width_in <- .hc_first_numeric_value(available_width_in[[1]])
  if (!base::is.finite(width_in) || width_in <= 0) {
    return(42L)
  }
  as.integer(base::max(32, base::min(72, base::floor(width_in * 7.5))))
}

.hc_llm_add_title_grob <- function(grob, title, text_size) {
  class(grob) <- unique(c("hc_llm_heatmap_plot", class(grob)))
  attr(grob, "llm_title") <- title
  attr(grob, "llm_text_size") <- text_size
  grob
}

.hc_llm_top_align_grob <- function(grob) {
  if (is.null(grob) || !inherits(grob, "grob")) {
    return(grob)
  }

  grob_vp <- tryCatch(grob$childrenvp[[1]]$parent, error = function(e) NULL)
  if (is.null(grob_vp)) {
    return(grob)
  }

  grob_vp$y <- grid::unit(1, "npc")
  grob_vp$justification <- c(0.5, 1)
  grob_vp$valid.just <- c(0.5, 1)
  grob$childrenvp[[1]]$parent <- grob_vp
  grob
}

#' @export
print.hc_llm_heatmap_plot <- function(x, ...) {
  title <- attr(x, "llm_title", exact = TRUE)
  text_size <- attr(x, "llm_text_size", exact = TRUE)
  inner_grob <- attr(x, "llm_inner_grob", exact = TRUE)

  if (is.null(inner_grob)) {
    inner_grob <- if (inherits(x, "gtable") && length(x$grobs) >= 2) x$grobs[[2]] else x
  }
  inner_grob <- .hc_llm_top_align_grob(inner_grob)
  if (is.null(title) && inherits(x, "gtable") && length(x$grobs) >= 1) {
    title <- tryCatch(x$grobs[[1]]$label, error = function(e) NULL)
  }
  if (is.null(text_size) || !is.numeric(text_size) || length(text_size) != 1 || is.na(text_size)) {
    text_size <- 4
  }

  title_wrapped <- NULL
  title_height_lines <- NULL
  if (!is.null(title) && is.character(title) && base::nzchar(title)) {
    device_width_in <- tryCatch(grDevices::dev.size("in")[[1]], error = function(e) NA_real_)
    title_wrapped <- stringr::str_wrap(
      title,
      width = .hc_llm_title_wrap_width(device_width_in)
    )
    title_line_count <- base::length(base::strsplit(title_wrapped, "\n", fixed = TRUE)[[1]])
    title_height_lines <- grid::unit(
      base::max(2.2, (title_line_count * 1.15) + 0.8),
      "lines"
    )
  }

  grid::grid.newpage()
  if (!is.null(title_wrapped) && !is.null(title_height_lines)) {
    grid::pushViewport(
      grid::viewport(
        layout = grid::grid.layout(
          nrow = 2,
          ncol = 1,
          heights = grid::unit.c(
            title_height_lines,
            grid::unit(1, "null")
          )
        ),
        clip = "off"
      )
    )
    grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 1, clip = "off"))
    grid::grid.text(
      label = title_wrapped,
      x = 0.5,
      y = 0.5,
      just = c("center", "center"),
      gp = grid::gpar(
        fontface = "bold",
        fontsize = base::max(10.5, text_size * 4),
        lineheight = 1.05
      )
    )
    grid::popViewport()
    grid::pushViewport(
      grid::viewport(
        layout.pos.row = 2,
        layout.pos.col = 1,
        x = 0.5,
        y = 1,
        width = 1,
        height = 1,
        just = c("center", "top"),
        clip = "off"
      )
    )
    grid::grid.draw(inner_grob)
    grid::popViewport()
    grid::popViewport()
  } else {
    grid::grid.draw(inner_grob)
  }
  invisible(x)
}

.hc_llm_darken_color <- function(col, factor = 0.35) {
  rgb_mat <- tryCatch(grDevices::col2rgb(col) / 255, error = function(e) NULL)
  if (is.null(rgb_mat)) {
    return("#222222")
  }
  rgb_new <- pmax(0, rgb_mat[, 1] * (1 - factor))
  grDevices::rgb(rgb_new[[1]], rgb_new[[2]], rgb_new[[3]])
}

.hc_llm_contrast_text_color <- function(col) {
  rgb_mat <- tryCatch(grDevices::col2rgb(col) / 255, error = function(e) NULL)
  if (is.null(rgb_mat)) {
    return("white")
  }
  luminance <- 0.299 * rgb_mat[1, 1] + 0.587 * rgb_mat[2, 1] + 0.114 * rgb_mat[3, 1]
  if (luminance > 0.65) "black" else "white"
}

.hc_llm_reorder_summary_for_display <- function(summary_tbl, heatmap_info = NULL) {
  if (is.null(heatmap_info) || is.null(heatmap_info$module_order) || length(heatmap_info$module_order) == 0) {
    return(summary_tbl)
  }
  keep <- heatmap_info$module_order[heatmap_info$module_order %in% summary_tbl$module]
  extras <- base::setdiff(summary_tbl$module, keep)
  ord <- base::match(base::c(keep, extras), summary_tbl$module)
  ord <- ord[!base::is.na(ord)]
  summary_tbl[ord, , drop = FALSE]
}

.hc_llm_heatmap_info <- function(hc) {
  cache_info <- .hc_heatmap_cache_info(hc@integration@cluster)
  raw_heatmap_obj <- cache_info$raw_heatmap_obj
  heatmap_obj <- cache_info$heatmap_obj
  mat <- cache_info$matrix
  if (is.null(mat) || is.null(rownames(mat))) {
    return(NULL)
  }

  row_ids <- cache_info$row_order
  if (is.null(row_ids) || length(row_ids) == 0) {
    row_ids <- rownames(mat)
  }

  module_order <- row_ids
  label_map <- tryCatch(hc@integration@cluster[["module_label_map"]], error = function(e) NULL)
  if (!is.null(label_map) && length(label_map) > 0) {
    label_map <- as.character(label_map)
    map_names <- tryCatch(names(hc@integration@cluster[["module_label_map"]]), error = function(e) NULL)
    if (!is.null(map_names) && length(map_names) == length(label_map)) {
      names(label_map) <- as.character(map_names)
    }
    hit <- label_map[row_ids]
    if (length(hit) == length(row_ids) && any(!is.na(hit))) {
      module_order <- as.character(hit)
      module_order[is.na(module_order) | module_order == ""] <- row_ids[is.na(module_order) | module_order == ""]
    }
  }

  list(
    matrix = mat,
    raw_heatmap_obj = raw_heatmap_obj,
    heatmap_obj = heatmap_obj,
    draw_supported = !is.null(mat) && nrow(mat) > 0 && ncol(mat) > 0,
    hcobject = tryCatch(.hc_as_bridge_object_for_cluster_plot(hc), error = function(e) NULL),
    col_order = cache_info$col_order,
    col_labels_display = tryCatch(base::as.character(hc@integration@cluster[["heatmap_column_labels_display"]]), error = function(e) NULL),
    row_ids = row_ids,
    module_order = module_order,
    module_by_row = module_order,
    stored_module_label_fontsize = tryCatch(as.numeric(hc@integration@cluster[["module_label_fontsize"]]), error = function(e) NA_real_),
    stored_module_label_pt_size = tryCatch(as.numeric(hc@integration@cluster[["module_label_pt_size"]]), error = function(e) NA_real_),
    stored_module_box_width_cm = tryCatch(as.numeric(hc@integration@cluster[["module_box_width_cm"]]), error = function(e) NA_real_),
    stored_heatmap_cell_size_mm = tryCatch(as.numeric(hc@integration@cluster[["heatmap_cell_size_mm"]]), error = function(e) NA_real_),
    stored_gfc_colors = tryCatch(base::as.character(hc@integration@cluster[["gfc_colors"]]), error = function(e) NULL),
    stored_gfc_scale_limits = tryCatch(.hc_as_numeric_safely(hc@integration@cluster[["gfc_scale_limits"]]), error = function(e) NULL),
    stored_overall_plot_scale = tryCatch(as.numeric(hc@integration@cluster[["overall_plot_scale"]]), error = function(e) NA_real_)
  )
}

.hc_llm_capture_heatmap_grob <- function(heatmap_obj) {
  clone <- .hc_safe_deep_clone(heatmap_obj, context = "LLM module function heatmap")
  grid::grid.grabExpr(
    ComplexHeatmap::draw(
      clone,
      newpage = FALSE,
      merge_legends = TRUE,
      show_annotation_legend = TRUE,
      show_heatmap_legend = TRUE
    )
  )
}

.hc_llm_extract_heatmap_matrix <- function(heatmap_obj) {
  if (inherits(heatmap_obj, "Heatmap")) {
    return(tryCatch(heatmap_obj@matrix, error = function(e) NULL))
  }
  if (inherits(heatmap_obj, "HeatmapList")) {
    return(tryCatch(heatmap_obj@ht_list[[1]]@matrix, error = function(e) NULL))
  }
  NULL
}

.hc_llm_align_display_labels <- function(source_ids, display_labels, target_ids) {
  source_ids <- base::as.character(source_ids)
  target_ids <- base::as.character(target_ids)
  display_labels <- base::as.character(display_labels)
  if (base::length(source_ids) == 0 || base::length(display_labels) != base::length(source_ids)) {
    return(target_ids)
  }
  if (base::identical(source_ids, target_ids)) {
    return(display_labels)
  }

  source_occ <- stats::ave(base::seq_along(source_ids), source_ids, FUN = base::seq_along)
  target_occ <- stats::ave(base::seq_along(target_ids), target_ids, FUN = base::seq_along)
  out <- base::vapply(
    base::seq_along(target_ids),
    function(i) {
      hit <- base::which(source_ids == target_ids[[i]] & source_occ == target_occ[[i]])
      if (base::length(hit) > 0) {
        display_labels[hit[[1]]]
      } else {
        target_ids[[i]]
      }
    },
    FUN.VALUE = base::character(1)
  )
  out
}

.hc_llm_normalize_gfc_scale_limits <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  x <- .hc_as_numeric_safely(x)
  if (base::length(x) == 1) {
    if (!base::is.finite(x) || x <= 0) {
      return(NULL)
    }
    return(c(-base::abs(x), base::abs(x)))
  }
  if (base::length(x) != 2 || base::any(!base::is.finite(x))) {
    return(NULL)
  }
  x <- base::sort(x)
  if (base::identical(x[[1]], x[[2]])) {
    return(NULL)
  }
  x
}

.hc_llm_default_heatmap_cell_size_mm <- function(n_heat_rows,
                                                 n_heat_cols,
                                                 overall_plot_scale = 1) {
  cell_size_mm <- 5.4
  if (n_heat_rows > 20) {
    cell_size_mm <- 4.9
  }
  if (n_heat_rows > 30) {
    cell_size_mm <- 4.3
  }
  if (n_heat_rows > 45) {
    cell_size_mm <- 3.8
  }
  if (n_heat_cols > 10) {
    cell_size_mm <- base::min(cell_size_mm, 4.6)
  }
  if (n_heat_cols <= 4 && n_heat_rows <= 24) {
    min_body_w_mm <- if (n_heat_cols <= 3) 30 else 34
    min_body_h_mm <- if (n_heat_rows <= 12) 90 else 108
    max_cell_mm <- if (n_heat_rows <= 12) 10 else 8
    boosted_cell_mm <- base::max(
      min_body_w_mm / base::max(1, n_heat_cols),
      min_body_h_mm / base::max(1, n_heat_rows)
    )
    cell_size_mm <- base::max(cell_size_mm, base::min(max_cell_mm, boosted_cell_mm))
  }
  cell_size_mm * overall_plot_scale
}

.hc_llm_cluster_heatmap_style_defaults <- function(module_labels,
                                                   n_heat_rows,
                                                   n_heat_cols,
                                                   overall_plot_scale = 1) {
  max_label_chars <- if (base::length(module_labels) == 0) {
    1
  } else {
    base::max(base::nchar(module_labels), na.rm = TRUE)
  }
  module_label_fontsize <- 4.2
  module_box_width_cm <- base::max(0.56, base::min(4.8, (max_label_chars * 0.09) + 0.15))
  module_label_pt_size <- 0.90
  cell_size_mm <- .hc_llm_default_heatmap_cell_size_mm(
    n_heat_rows = n_heat_rows,
    n_heat_cols = n_heat_cols,
    overall_plot_scale = overall_plot_scale
  )
  list(
    module_label_fontsize = module_label_fontsize,
    module_label_pt_size = module_label_pt_size,
    module_box_width_cm = module_box_width_cm,
    cell_size_mm = cell_size_mm
  )
}

.hc_llm_resolve_heatmap_style <- function(heatmap_info,
                                          module_labels,
                                          n_heat_rows,
                                          n_heat_cols,
                                          mat_use,
                                          module_label_fontsize = NULL,
                                          module_label_pt_size = NULL,
                                          module_box_width_cm = NULL) {
  stored_module_label_fontsize <- .hc_first_numeric_value(heatmap_info$stored_module_label_fontsize[[1]])
  if (!base::is.finite(stored_module_label_fontsize) || stored_module_label_fontsize <= 0) {
    stored_module_label_fontsize <- NULL
  }
  stored_module_label_pt_size <- .hc_first_numeric_value(heatmap_info$stored_module_label_pt_size[[1]])
  if (!base::is.finite(stored_module_label_pt_size) || stored_module_label_pt_size <= 0) {
    stored_module_label_pt_size <- NULL
  }
  stored_module_box_width_cm <- .hc_first_numeric_value(heatmap_info$stored_module_box_width_cm[[1]])
  if (!base::is.finite(stored_module_box_width_cm) || stored_module_box_width_cm <= 0) {
    stored_module_box_width_cm <- NULL
  }
  stored_heatmap_cell_size_mm <- .hc_first_numeric_value(heatmap_info$stored_heatmap_cell_size_mm[[1]])
  if (!base::is.finite(stored_heatmap_cell_size_mm) || stored_heatmap_cell_size_mm <= 0) {
    stored_heatmap_cell_size_mm <- NULL
  }
  stored_overall_plot_scale <- .hc_first_numeric_value(heatmap_info$stored_overall_plot_scale[[1]])
  if (!base::is.finite(stored_overall_plot_scale) || stored_overall_plot_scale <= 0) {
    stored_overall_plot_scale <- 1
  }
  stored_overall_plot_scale <- base::max(0.5, base::min(3, stored_overall_plot_scale))

  stored_gfc_colors <- tryCatch(base::as.character(heatmap_info$stored_gfc_colors), error = function(e) NULL)
  if (is.null(stored_gfc_colors) || base::length(stored_gfc_colors) < 2 ||
    any(base::is.na(stored_gfc_colors)) || any(stored_gfc_colors == "")) {
    stored_gfc_colors <- .hc_default_gfc_colors()
  }
  stored_gfc_scale_limits <- .hc_llm_normalize_gfc_scale_limits(heatmap_info$stored_gfc_scale_limits)
  if (is.null(stored_gfc_scale_limits)) {
    fallback_lim <- .hc_max_finite(base::abs(mat_use))
    if (!base::is.finite(fallback_lim) || fallback_lim <= 0) {
      fallback_lim <- 2
    }
    stored_gfc_scale_limits <- c(-base::abs(fallback_lim), base::abs(fallback_lim))
  }

  compact_defaults <- .hc_llm_cluster_heatmap_style_defaults(
    module_labels = module_labels,
    n_heat_rows = n_heat_rows,
    n_heat_cols = n_heat_cols,
    overall_plot_scale = stored_overall_plot_scale
  )

  resolved_module_label_fontsize <- if (!is.null(module_label_fontsize)) {
    as.numeric(module_label_fontsize[[1]])
  } else {
    compact_defaults$module_label_fontsize
  }
  resolved_module_label_pt_size <- if (!is.null(module_label_pt_size)) {
    as.numeric(module_label_pt_size[[1]])
  } else {
    compact_defaults$module_label_pt_size
  }
  resolved_module_box_width_cm <- if (!is.null(module_box_width_cm)) {
    as.numeric(module_box_width_cm[[1]])
  } else {
    compact_defaults$module_box_width_cm
  }
  cell_size_mm <- if (!is.null(stored_heatmap_cell_size_mm)) {
    stored_heatmap_cell_size_mm
  } else {
    compact_defaults$cell_size_mm
  }

  list(
    module_label_fontsize = resolved_module_label_fontsize,
    module_label_pt_size = resolved_module_label_pt_size,
    module_box_width_cm = resolved_module_box_width_cm,
    cell_size_mm = cell_size_mm,
    gfc_colors = stored_gfc_colors,
    gfc_scale_limits = stored_gfc_scale_limits,
    overall_plot_scale = stored_overall_plot_scale
  )
}

.hc_llm_capture_combined_heatmap_grob <- function(heatmap_info,
                                                  summary_tbl,
                                                  max_chars,
                                                  text_size,
                                                  module_label_fontsize = NULL,
                                                  module_label_pt_size = NULL,
                                                  module_box_width_cm = NULL,
                                                  heatmap_col_order = NULL,
                                                  heatmap_cluster_columns = FALSE) {
  source_obj <- if (!is.null(heatmap_info$heatmap_obj)) heatmap_info$heatmap_obj else heatmap_info$raw_heatmap_obj
  mat <- heatmap_info$matrix
  if (is.null(mat) && !is.null(source_obj)) {
    mat <- .hc_llm_extract_heatmap_matrix(source_obj)
  }
  if (is.null(mat) || is.null(rownames(mat)) || nrow(mat) == 0) {
    stop("No valid heatmap matrix available for combined LLM visualization.")
  }

  row_ids <- as.character(heatmap_info$row_ids)
  module_by_row <- as.character(heatmap_info$module_by_row)
  keep_row <- module_by_row %in% as.character(summary_tbl$module)
  row_ids <- row_ids[keep_row]
  module_by_row <- module_by_row[keep_row]
  if (length(row_ids) == 0) {
    stop("No overlapping modules between LLM summaries and stored heatmap.")
  }

  summary_tbl$module <- as.character(summary_tbl$module)
  if (!"term_plot" %in% colnames(summary_tbl)) {
    if (!"short_title" %in% colnames(summary_tbl)) {
      summary_tbl$short_title <- NA_character_
    }
    summary_tbl$term_plot <- vapply(
      ifelse(
        is.na(summary_tbl$short_title) | summary_tbl$short_title == "",
        as.character(summary_tbl$overarching_function),
        as.character(summary_tbl$short_title)
      ),
      .hc_llm_prepare_display_title,
      FUN.VALUE = character(1),
      max_chars = max_chars
    )
  } else {
    summary_tbl$term_plot <- vapply(
      as.character(summary_tbl$term_plot),
      .hc_llm_prepare_display_title,
      FUN.VALUE = character(1),
      max_chars = max_chars
    )
  }
  summary_tbl$text_color <- vapply(
    as.character(summary_tbl$module_color),
    .hc_llm_darken_color,
    FUN.VALUE = character(1)
  )
  summary_tbl$module_color <- as.character(summary_tbl$module_color)
  summary_tbl$module_color[is.na(summary_tbl$module_color) | summary_tbl$module_color == ""] <- "grey70"
  summary_tbl$label_color <- vapply(summary_tbl$module_color, .hc_llm_contrast_text_color, FUN.VALUE = character(1))

  fallback_col_ids <- if (!is.null(source_obj)) .hc_llm_extract_column_ids(source_obj, mat) else base::colnames(mat)
  prepared_cols <- .hc_prepare_plot_heatmap_columns(
    mat = mat[row_ids, , drop = FALSE],
    cluster_columns = heatmap_cluster_columns,
    plot_order = heatmap_col_order,
    main_order = heatmap_info$col_order,
    fallback_order = fallback_col_ids,
    context = "LLM module function heatmap"
  )
  mat_use <- prepared_cols$mat
  column_labels_display <- .hc_llm_align_display_labels(
    source_ids = base::colnames(mat),
    display_labels = heatmap_info$col_labels_display,
    target_ids = base::colnames(mat_use)
  )

  module_lookup <- stats::setNames(summary_tbl$term_plot, summary_tbl$module)
  color_lookup <- stats::setNames(summary_tbl$text_color, summary_tbl$module)
  fill_lookup <- stats::setNames(summary_tbl$module_color, summary_tbl$module)
  label_lookup <- stats::setNames(summary_tbl$label_color, summary_tbl$module)

  term_labels <- unname(module_lookup[module_by_row])
  term_cols <- unname(color_lookup[module_by_row])
  fill_cols <- unname(fill_lookup[module_by_row])
  label_cols <- unname(label_lookup[module_by_row])

  term_labels[is.na(term_labels)] <- ""
  term_cols[is.na(term_cols) | term_labels == ""] <- "#00000000"
  fill_cols[is.na(fill_cols) | fill_cols == ""] <- "grey70"
  label_cols[is.na(label_cols) | label_cols == ""] <- "black"

  max_width_chars <- if (length(term_labels) == 0) 0 else max(nchar(term_labels), na.rm = TRUE)
  text_width_cm <- max(7.5, min(18, (max_width_chars * 0.14) + 0.8))

  module_levels <- unique(module_by_row)
  module_col_map <- stats::setNames(fill_lookup[module_levels], module_levels)
  module_col_map[is.na(module_col_map) | module_col_map == ""] <- "grey70"
  style <- .hc_llm_resolve_heatmap_style(
    heatmap_info = heatmap_info,
    module_labels = module_by_row,
    n_heat_rows = nrow(mat_use),
    n_heat_cols = ncol(mat_use),
    mat_use = mat_use,
    module_label_fontsize = module_label_fontsize,
    module_label_pt_size = module_label_pt_size,
    module_box_width_cm = module_box_width_cm
  )
  label_fontsize <- style$module_label_fontsize
  module_box_width_cm <- style$module_box_width_cm
  module_pt_size <- style$module_label_pt_size
  cell_mm <- style$cell_size_mm
  overall_plot_scale <- style$overall_plot_scale
  column_gap_spec <- .hc_heatmap_column_gap_spec(
    hcobject = heatmap_info$hcobject,
    cols = base::colnames(mat_use),
    cluster_columns = heatmap_cluster_columns,
    gap_mm = 0.6 * overall_plot_scale
  )
  heatmap_body_w_mm <- max(18, (ncol(mat_use) * cell_mm) + column_gap_spec$total_gap_mm)
  heatmap_body_h_mm <- max(20, nrow(mat_use) * cell_mm)
  shared_heatmap_line_lwd <- 0.5
  module_box_border_gp <- grid::gpar(col = "black", lwd = shared_heatmap_line_lwd)
  heatmap_cell_border_gp <- grid::gpar(col = "black", lwd = shared_heatmap_line_lwd)
  dendrogram_line_gp <- grid::gpar(col = "black", lwd = shared_heatmap_line_lwd)
  max_col_chars <- if (base::length(column_labels_display) > 0) {
    max(nchar(column_labels_display), na.rm = TRUE)
  } else {
    10
  }
  shared_column_name_max_cm <- max(
    2.8,
    min(
      8.0,
      (max_col_chars * (10 * overall_plot_scale) * 0.022) + 0.6
    )
  )
  row_dend_width_mm <- if (nrow(mat_use) > 1) {
    max(8, min(14, nrow(mat_use) * 0.75))
  } else {
    6
  }
  row_dend_width_mm <- row_dend_width_mm * overall_plot_scale
  column_dend_height_mm <- if (ncol(mat_use) > 1) {
    max(10, min(18, ncol(mat_use) * 2.2))
  } else {
    6
  }
  column_dend_height_mm <- column_dend_height_mm * overall_plot_scale
  gfc_palette <- grDevices::colorRampPalette(style$gfc_colors)(51)
  gfc_col_fun <- circlize::colorRamp2(
    seq(style$gfc_scale_limits[1], style$gfc_scale_limits[2], length.out = length(gfc_palette)),
    gfc_palette
  )
  gfc_scale_breaks <- pretty(style$gfc_scale_limits, n = 5)
  gfc_scale_breaks <- gfc_scale_breaks[
    gfc_scale_breaks >= style$gfc_scale_limits[1] - .Machine$double.eps^0.5 &
      gfc_scale_breaks <= style$gfc_scale_limits[2] + .Machine$double.eps^0.5
  ]
  if (!any(base::abs(gfc_scale_breaks) < .Machine$double.eps^0.5)) {
    gfc_scale_breaks <- base::sort(base::unique(base::c(gfc_scale_breaks, 0)))
  }
  if (base::length(gfc_scale_breaks) < 3) {
    gfc_scale_breaks <- base::seq(style$gfc_scale_limits[1], style$gfc_scale_limits[2], length.out = 5)
  }
  gfc_scale_labels <- base::formatC(gfc_scale_breaks, format = "fg", digits = 3)
  gfc_scale_labels <- base::trimws(gfc_scale_labels)
  gfc_label_width <- base::max(base::nchar(gfc_scale_labels), na.rm = TRUE)
  gfc_scale_labels <- base::format(gfc_scale_labels, width = gfc_label_width, justify = "right")
  legend_height_mm <- max(24, 4.5 * max(1, length(gfc_scale_breaks))) * overall_plot_scale
  heat_legend_param <- list(
    title = "GFC",
    at = gfc_scale_breaks,
    labels = gfc_scale_labels,
    title_gp = grid::gpar(fontsize = 7.6 * overall_plot_scale, fontface = "bold"),
    labels_gp = grid::gpar(fontsize = 6.6 * overall_plot_scale),
    grid_width = grid::unit(3.2 * overall_plot_scale, "mm"),
    legend_height = grid::unit(legend_height_mm, "mm")
  )
  module_label_fit <- .hc_module_label_fit_pt(
    module_label_pt_size = module_pt_size,
    module_box_width_cm = module_box_width_cm,
    module_labels_display = module_by_row,
    n_heat_rows = nrow(mat_use),
    cell_size_mm = cell_mm,
    module_label_fontsize = label_fontsize,
    use_fontsize_request = !base::is.null(module_label_fontsize) &&
      base::is.null(module_label_pt_size),
    fontface = "bold"
  )
  module_label_pt_size_unit <- if (base::length(module_by_row) > 0 &&
    base::length(module_label_fit$pt_size) == base::length(module_by_row)) {
    grid::unit(module_label_fit$pt_size, "pt")
  } else {
    grid::unit(module_pt_size, "snpc")
  }
  module_label_fontsize_draw <- .hc_module_label_effective_fontsize(
    module_label_fit = module_label_fit,
    fallback_fontsize = label_fontsize,
    fallback_pt_size = module_pt_size
  )

  module_box_anno <- .hc_module_label_box_annotation(
    values = module_by_row,
    colors = module_col_map,
    labels = module_by_row,
    label_color = "white",
    label_fontsize_pt = module_label_fit$base_pt_size,
    fontface = "bold",
    width_cm = module_box_width_cm,
    border_gp = module_box_border_gp,
    which = "row"
  )

  right_anno <- ComplexHeatmap::HeatmapAnnotation(
    modules = module_box_anno,
    llm = ComplexHeatmap::anno_text(
      term_labels,
      which = "row",
      just = "left",
      location = 0,
      gp = grid::gpar(col = term_cols, fontsize = max(7.2, text_size * 2.25)),
      width = grid::unit(text_width_cm, "cm")
    ),
    which = "row",
    show_legend = FALSE,
    show_annotation_name = FALSE,
    gap = grid::unit(2, "mm")
  )

  full_row_match <- identical(as.character(row_ids), as.character(heatmap_info$row_ids))
  row_dend <- if (isTRUE(full_row_match) && !is.null(source_obj)) {
    .hc_llm_extract_dendrogram(source_obj, which = "row")
  } else {
    NULL
  }
  col_dend <- prepared_cols$col_dend
  if (is.null(col_dend) && is.null(heatmap_col_order) && !is.null(source_obj)) {
    col_dend <- .hc_llm_extract_dendrogram(source_obj, which = "column")
  }
  row_dend <- .hc_llm_validate_axis_dendrogram(row_dend, axis_ids = row_ids)
  col_dend <- .hc_llm_validate_axis_dendrogram(col_dend, axis_ids = base::colnames(mat_use))

  combined_ht_args <- list(
    matrix = mat_use,
    name = "GFC",
    col = gfc_col_fun,
    clustering_distance_rows = "euclidean",
    clustering_distance_columns = "euclidean",
    clustering_method_rows = "complete",
    clustering_method_columns = "complete",
    cluster_rows = if (is.null(row_dend)) FALSE else row_dend,
    cluster_columns = if (is.null(col_dend)) FALSE else col_dend,
    row_dend_reorder = FALSE,
    show_row_names = FALSE,
    show_heatmap_legend = TRUE,
    show_row_dend = !is.null(row_dend),
    show_column_dend = !is.null(col_dend),
    right_annotation = right_anno,
    width = grid::unit(heatmap_body_w_mm, "mm"),
    height = grid::unit(heatmap_body_h_mm, "mm"),
    row_dend_width = grid::unit(row_dend_width_mm, "mm"),
    column_dend_height = grid::unit(column_dend_height_mm, "mm"),
    row_dend_gp = dendrogram_line_gp,
    column_dend_gp = dendrogram_line_gp,
    column_names_max_height = grid::unit(shared_column_name_max_cm, "cm"),
    rect_gp = heatmap_cell_border_gp,
    row_names_side = "right",
    column_names_side = "bottom",
    column_names_rot = 90,
    column_labels = column_labels_display,
    column_names_centered = FALSE,
    column_names_gp = grid::gpar(fontsize = 10 * overall_plot_scale),
    heatmap_legend_param = heat_legend_param
  )
  combined_ht_args <- .hc_heatmap_add_column_gap_args(
    combined_ht_args,
    column_gap_spec = column_gap_spec,
    title_gp = grid::gpar(fontsize = 10 * overall_plot_scale, fontface = "bold")
  )
  combined_ht <- do.call(ComplexHeatmap::Heatmap, combined_ht_args)
  right_pad_mm <- max(28, 6 + (gfc_label_width * 2.2))
  capture_width_mm <- base::max(
    120,
    row_dend_width_mm +
      heatmap_body_w_mm +
      (module_box_width_cm * 10) +
      (text_width_cm * 10) +
      right_pad_mm +
      38
  )
  capture_height_mm <- base::max(
    70,
    column_dend_height_mm +
      heatmap_body_h_mm +
      (shared_column_name_max_cm * 10) +
      30
  )
  grid::grid.grabExpr(
    ComplexHeatmap::draw(
      combined_ht,
      newpage = FALSE,
      merge_legends = TRUE,
      show_annotation_legend = FALSE,
      show_heatmap_legend = TRUE,
      heatmap_legend_side = "right",
      annotation_legend_side = "right",
      padding = grid::unit(c(8, 18, 18, right_pad_mm) * overall_plot_scale, "mm")
    ),
    width = capture_width_mm / 25.4,
    height = capture_height_mm / 25.4
  )
}

.hc_llm_attach_right_annotation <- function(heatmap_obj, annotation, heatmap_info = NULL) {
  pin_heatmap_size <- function(ht) {
    if (is.null(heatmap_info) || is.null(heatmap_info$matrix)) {
      return(ht)
    }
    cell_mm_raw <- heatmap_info$stored_heatmap_cell_size_mm
    if (is.null(cell_mm_raw) || length(cell_mm_raw) == 0) {
      return(ht)
    }
    cell_mm <- .hc_first_numeric_value(cell_mm_raw[[1]])
    if (!is.finite(cell_mm) || cell_mm <= 0) {
      return(ht)
    }
    n_cols <- ncol(heatmap_info$matrix)
    n_rows <- nrow(heatmap_info$matrix)
    if (!is.finite(n_cols) || !is.finite(n_rows) || n_cols <= 0 || n_rows <= 0) {
      return(ht)
    }
    ht@matrix_param$width <- grid::unit(n_cols * cell_mm, "mm")
    ht@matrix_param$height <- grid::unit(n_rows * cell_mm, "mm")
    ht
  }

  append_annotation <- function(ht) {
    ht <- pin_heatmap_size(ht)
    existing <- tryCatch(ht@right_annotation, error = function(e) NULL)
    if (is.null(existing)) {
      ht@right_annotation <- annotation
    } else {
      ht@right_annotation <- c(existing, annotation)
    }
    ht
  }

  cloned_obj <- .hc_safe_deep_clone(
    heatmap_obj,
    context = "LLM module function heatmap reuse"
  )

  if (inherits(cloned_obj, "Heatmap")) {
    return(append_annotation(cloned_obj))
  }

  if (inherits(cloned_obj, "HeatmapList")) {
    heatmap_idx <- which(vapply(cloned_obj@ht_list, inherits, logical(1), "Heatmap"))
    if (length(heatmap_idx) == 0) {
      stop("No heatmap component available for LLM heatmap reuse.")
    }
    cloned_obj@ht_list[[heatmap_idx[[1]]]] <- append_annotation(cloned_obj@ht_list[[heatmap_idx[[1]]]])
    return(cloned_obj)
  }

  stop("Unsupported cached heatmap object for LLM heatmap reuse.")
}

.hc_llm_extract_column_ids <- function(heatmap_obj, mat) {
  col_ids <- colnames(mat)
  ord <- suppressWarnings(
    tryCatch(ComplexHeatmap::column_order(heatmap_obj), error = function(e) NULL)
  )
  if (is.list(ord) && length(ord) > 0) {
    ord <- ord[[1]]
  }
  if (is.numeric(ord) && length(ord) == ncol(mat)) {
    col_ids <- col_ids[ord]
  } else if (is.character(ord) && length(ord) > 0) {
    col_ids <- ord[ord %in% col_ids]
  }
  col_ids
}

.hc_llm_extract_dendrogram <- function(heatmap_obj, which = c("row", "column")) {
  which <- match.arg(which)
  fn <- if (identical(which, "row")) ComplexHeatmap::row_dend else ComplexHeatmap::column_dend
  out <- suppressWarnings(tryCatch(fn(heatmap_obj), error = function(e) NULL))
  if (is.list(out) && length(out) > 0) {
    out <- out[[1]]
  }
  if (inherits(out, "dendrogram") || inherits(out, "hclust")) {
    return(out)
  }
  NULL
}

.hc_llm_dendrogram_leaf_count <- function(dend) {
  if (is.null(dend)) {
    return(NA_integer_)
  }
  if (inherits(dend, "hclust")) {
    n <- tryCatch(base::length(dend$order), error = function(e) NA_integer_)
    return(as.integer(n))
  }
  if (inherits(dend, "dendrogram")) {
    n <- tryCatch(base::length(stats::order.dendrogram(dend)), error = function(e) NA_integer_)
    return(as.integer(n))
  }
  NA_integer_
}

.hc_llm_validate_axis_dendrogram <- function(dend, axis_ids = NULL) {
  if (is.null(dend)) {
    return(NULL)
  }
  expected_n <- base::length(axis_ids)
  if (!base::is.numeric(expected_n) || expected_n <= 0) {
    return(NULL)
  }
  dend_n <- .hc_llm_dendrogram_leaf_count(dend)
  if (!base::is.finite(dend_n) || dend_n != expected_n) {
    return(NULL)
  }

  dend_labels <- tryCatch(base::as.character(labels(dend)), error = function(e) NULL)
  if (!is.null(dend_labels) && base::length(dend_labels) > 0) {
    dend_labels <- dend_labels[!base::is.na(dend_labels) & base::nzchar(dend_labels)]
    axis_ids_chr <- base::as.character(axis_ids)
    if (base::length(dend_labels) > 0 &&
      !base::setequal(base::unique(dend_labels), base::unique(axis_ids_chr))) {
      return(NULL)
    }
  }
  dend
}

.hc_llm_prepare_display_term <- function(term, max_chars = 90) {
  term <- as.character(term[[1]])
  if (!nzchar(term) || is.na(term)) {
    return("No interpretation available.")
  }
  term <- .hc_llm_clean_text(term)
  term <- stringr::str_squish(term)
  term <- sub("^This module primarily reflects\\s+", "", term, ignore.case = TRUE)
  term <- sub("^This module reflects\\s+", "", term, ignore.case = TRUE)
  term <- sub("^This module primarily captures\\s+", "", term, ignore.case = TRUE)
  term <- sub("^This module captures\\s+", "", term, ignore.case = TRUE)
  term <- sub("^This module is characterized by\\s+", "", term, ignore.case = TRUE)
  term <- sub("^This module represents\\s+", "", term, ignore.case = TRUE)
  term <- sub("\\.$", "", term)
  if (is.infinite(max_chars) || is.na(max_chars) || max_chars <= 0) {
    return(term)
  }
  stringr::str_trunc(term, width = max_chars)
}

.hc_llm_prepare_display_title <- function(term, max_chars = 60) {
  term <- .hc_llm_short_title_fallback(term = term, max_chars = max_chars)
  term <- .hc_llm_clean_text(term)
  term <- stringr::str_squish(term)
  if (!nzchar(term) || is.na(term)) {
    return("No interpretation available")
  }
  stringr::str_trunc(term, width = max_chars)
}
