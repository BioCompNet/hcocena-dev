test_that("regression: utils no longer references working_director", {
  src <- paste(deparse(get("run_expression_analysis_1_body", asNamespace("hcocena"))), collapse = "\n")
  expect_false(grepl("working_director", src, fixed = TRUE))
})


test_that("regression: rho support is removed from the correlation path", {
  src <- paste(deparse(get("pwcorr", asNamespace("hcocena"))), collapse = "\n")
  expect_false(grepl(".hc_require_namespace(\"propr\"", src, fixed = TRUE))
  expect_false(grepl("getExportedValue(\"propr\", \"propr\")", src, fixed = TRUE))
  expect_false(grepl("corr_method == 'rho'", src, fixed = TRUE))
})


test_that("regression: hub centrality helpers return finite legacy-style outputs", {
  dc_fun <- get("weighted_DC", asNamespace("hcocena"))
  cc_fun <- get("weighted_CC", asNamespace("hcocena"))
  bc_fun <- get("weighted_BC", asNamespace("hcocena"))
  combined_fun <- get("combined_centrality", asNamespace("hcocena"))

  g <- igraph::graph_from_edgelist(
    matrix(c("A", "B",
             "B", "C",
             "A", "C"),
           ncol = 2, byrow = TRUE),
    directed = FALSE
  )
  igraph::E(g)$weight <- c(0.8, 0.6, 0.4)

  dc <- dc_fun(g)
  cc <- cc_fun(g)
  bc <- bc_fun(g)
  combined <- combined_fun(g)

  expect_equal(length(dc), 3)
  expect_equal(length(cc), 3)
  expect_equal(length(bc), 3)
  expect_equal(sum(dc), 2)
  expect_true(all(is.finite(dc)))
  expect_true(all(is.finite(cc)))
  expect_true(all(is.finite(bc)))
  expect_true(any(bc > 0))
  expect_identical(sort(unique(combined$node)), c("A", "B", "C"))
  expect_identical(sort(unique(rownames(combined))), c("A", "B", "C"))

  g_disc <- igraph::graph_from_edgelist(
    matrix(c("A", "B",
             "C", "D"),
            ncol = 2, byrow = TRUE),
    directed = FALSE
  )
  igraph::E(g_disc)$weight <- c(0.7, 0.9)
  expect_true(all(is.finite(cc_fun(g_disc))))
  expect_true(all(is.finite(suppressWarnings(bc_fun(g_disc)))))
})


test_that("regression: longitudinal module labels preserve split suffixes", {
  label_fun <- get(".hc_normalize_longitudinal_module_labels", asNamespace("hcocena"))
  structured_fun <- get(".hc_is_structured_module_label", asNamespace("hcocena"))
  map_fun <- get(".hc_resolve_module_label_map_for_colors", asNamespace("hcocena"))
  color_order_fun <- get(".hc_module_colors_in_cluster_order", asNamespace("hcocena"))
  natural_fun <- get(".hc_natural_module_order", asNamespace("hcocena"))
  reorder_lookup_fun <- get(".hc_reorder_module_lookup_natural", asNamespace("hcocena"))

  expect_true(all(structured_fun(c("M1", "M2", "M1.1", "M1.2"))))
  expect_equal(
    label_fun(
      module_labels = c("M1.1", "M1.2", "M2"),
      module_colors = c("red", "blue", "green")
    ),
    c("M1.1", "M1.2", "M2")
  )
  expect_equal(
    label_fun(
      module_labels = c("red", "blue"),
      module_colors = c("red", "blue")
    ),
    c("M1", "M2")
  )
  expect_equal(
    map_fun(
      label_map = stats::setNames(c("M1.1", "M1.2"), c("red", "blue")),
      module_colors = c("red", "blue")
    ),
    stats::setNames(c("M1.1", "M1.2"), c("red", "blue"))
  )
  expect_equal(
    map_fun(
      label_map = stats::setNames(c("red", "blue"), c("M1.1", "M1.2")),
      module_colors = c("red", "blue")
    ),
    stats::setNames(c("M1.1", "M1.2"), c("red", "blue"))
  )
  expect_equal(
    color_order_fun(
      cluster_info = data.frame(
        color = c("gold_b", "gold_c", "gold_a", "gold_a"),
        gene_n = c("g1", "g2", "g3", "g4"),
        stringsAsFactors = FALSE
      ),
      module_genes = list(gold_a = "g3", gold_b = "g1", gold_c = "g2")
    ),
    c("gold_b", "gold_c", "gold_a")
  )
  expect_equal(
    natural_fun(c("M4", "M10", "M8", "M5", "M7", "M3", "M9", "M6", "M1.2", "M1.3", "M1.1", "M2")),
    c("M1.1", "M1.2", "M1.3", "M2", "M3", "M4", "M5", "M6", "M7", "M8", "M9", "M10")
  )
  expect_equal(
    reorder_lookup_fun(
      data.frame(
        module = c("M4", "M10", "M8", "M5", "M7", "M3", "M9", "M6", "M1.2", "M1.3", "M1.1", "M2"),
        module_color = seq_len(12),
        stringsAsFactors = FALSE
      )
    )$module,
    c("M1.1", "M1.2", "M1.3", "M2", "M3", "M4", "M5", "M6", "M7", "M8", "M9", "M10")
  )
})


test_that("regression: heatmap gets subtle smart column gaps only at useful group boundaries", {
  gap_fun <- get(".hc_heatmap_column_gap_spec", asNamespace("hcocena"))
  ggplot_layout_fun <- get(".hc_heatmap_ggplot_column_layout", asNamespace("hcocena"))

  hc_stub <- list(
    layers = list(set1 = TRUE, set2 = TRUE),
    layers_names = c("RNA", "PROT"),
    data = list(
      set1_anno = data.frame(Group = c("T1", "T2", "T3"), Phase = c("early", "early", "late"), stringsAsFactors = FALSE),
      set2_anno = data.frame(Group = c("T1", "T2", "T3"), Phase = c("early", "early", "late"), stringsAsFactors = FALSE)
    ),
    global_settings = list(voi = "Group"),
    layer_specific_outputs = NULL
  )

  default_off <- gap_fun(
    hcobject = hc_stub,
    cols = c("T1", "T2", "T3", "T1", "T2", "T3"),
    cluster_columns = FALSE,
    gap_mm = 0.6
  )
  expect_null(default_off$column_split)
  expect_equal(default_off$total_gap_mm, 0)

  by_layer <- gap_fun(
    hcobject = hc_stub,
    cols = c("T1", "T2", "T3", "T1", "T2", "T3"),
    cluster_columns = FALSE,
    gap_mm = 0.6,
    enabled = TRUE
  )
  expect_identical(by_layer$source, "layer_name")
  expect_equal(base::as.integer(by_layer$column_split), c(1L, 1L, 1L, 2L, 2L, 2L))
  expect_equal(base::levels(by_layer$column_split), c("RNA", "PROT"))
  expect_equal(by_layer$slice_count, 2L)
  expect_equal(by_layer$slice_titles, c("RNA", "PROT"))
  expect_gt(by_layer$total_gap_mm, 0)
  ggplot_layout <- ggplot_layout_fun(
    cols = c("T1", "T2", "T3", "T1", "T2", "T3"),
    column_gap_spec = by_layer,
    default_cell_mm = 5
  )
  expect_equal(ggplot_layout$slice_df$title, c("RNA", "PROT"))
  expect_gt(ggplot_layout$x[[4]] - ggplot_layout$x[[3]], 1)

  by_prefix <- gap_fun(
    hcobject = hc_stub,
    cols = c("MC1_T1_RNA", "MC1_T1_PROT", "MC2_T1_RNA", "MC2_T1_PROT"),
    cluster_columns = FALSE,
    gap_mm = 0.6,
    enabled = TRUE
  )
  expect_identical(by_prefix$source, "prefix_before_layer")
  expect_equal(base::as.integer(by_prefix$column_split), c(1L, 1L, 2L, 2L))
  expect_equal(by_prefix$slice_titles, c("", ""))

  by_metadata <- gap_fun(
    hcobject = hc_stub,
    cols = c("T1", "T2", "T1", "T2", "T3", "T3"),
    cluster_columns = FALSE,
    gap_mm = 0.6,
    metadata_column = "Phase"
  )
  expect_identical(by_metadata$source, "metadata:Phase")
  expect_equal(base::as.integer(by_metadata$column_split), c(1L, 1L, 1L, 1L, 2L, 2L))
  expect_equal(by_metadata$slice_titles, c("early", "late"))

  by_metadata_singletons <- gap_fun(
    hcobject = hc_stub,
    cols = c("T1", "T3"),
    cluster_columns = FALSE,
    gap_mm = 0.6,
    metadata_column = "Phase"
  )
  expect_identical(by_metadata_singletons$source, "metadata:Phase")
  expect_equal(base::as.integer(by_metadata_singletons$column_split), c(1L, 2L))
  expect_equal(by_metadata_singletons$slice_titles, c("early", "late"))

  clustered <- gap_fun(
    hcobject = hc_stub,
    cols = c("T1", "T2", "T3", "T1", "T2", "T3"),
    cluster_columns = TRUE,
    gap_mm = 0.6,
    enabled = TRUE
  )
  expect_null(clustered$column_split)
  expect_equal(clustered$total_gap_mm, 0)
})


test_that("regression: additional heatmap paths reuse layer gap and layer-title logic", {
  fun_enrich_path <- test_path("..", "..", "R", "hc_functional_enrichment.R")
  knowledge_path <- test_path("..", "..", "R", "hc_plot_enrichment_upstream_network.R")
  if (!file.exists(fun_enrich_path) || !file.exists(knowledge_path)) {
    skip("Source files are not available in the installed-package test context.")
  }

  fun_enrich_src <- paste(
    readLines(fun_enrich_path, warn = FALSE),
    collapse = "\n"
  )
  llm_src <- paste(deparse(get(".hc_llm_capture_combined_heatmap_grob", asNamespace("hcocena"))), collapse = "\n")
  upstream_src <- paste(deparse(get(".hc_ui_build_upstream_combined_heatmap", asNamespace("hcocena"))), collapse = "\n")
  knowledge_src <- paste(
    readLines(knowledge_path, warn = FALSE),
    collapse = "\n"
  )

  expect_true(grepl(".hc_heatmap_column_gap_spec", fun_enrich_src, fixed = TRUE))
  expect_true(grepl(".hc_heatmap_add_column_gap_args", fun_enrich_src, fixed = TRUE))
  expect_true(grepl(".hc_heatmap_column_gap_spec", llm_src, fixed = TRUE))
  expect_true(grepl(".hc_heatmap_add_column_gap_args", llm_src, fixed = TRUE))
  expect_true(grepl(".hc_heatmap_column_gap_spec", upstream_src, fixed = TRUE))
  expect_true(grepl(".hc_heatmap_add_column_gap_args", upstream_src, fixed = TRUE))
  expect_true(grepl(".hc_heatmap_column_gap_spec", knowledge_src, fixed = TRUE))
  expect_true(grepl(".hc_heatmap_ggplot_column_layout", knowledge_src, fixed = TRUE))
})


test_that("regression: regrouped heatmap uses the modern heatmap renderer and compact gene-count text mode", {
  regroup_src <- paste(deparse(get(".hc_change_grouping_parameter_driver", asNamespace("hcocena"))), collapse = "\n")

  expect_true(grepl("plot_cluster_heatmap_new", regroup_src, fixed = TRUE))
  expect_true(grepl("gene_count_mode = \"text\"", regroup_src, fixed = TRUE))
  expect_true(grepl("c(\"global_settings\", \"voi\")", regroup_src, fixed = TRUE))
  expect_true(grepl("\"regrouped\"", regroup_src, fixed = TRUE))
  expect_false(grepl("replot_cluster_heatmap", regroup_src, fixed = TRUE))
})


test_that("regression: longitudinal step1 can loop all layers with unique slots and file prefixes", {
  hc <- hc_init()
  se1 <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(1, nrow = 1, ncol = 1, dimnames = list("g1", "s1"))),
    colData = S4Vectors::DataFrame(PatID = "p1", Timepoint_rough_num = "1", row.names = "s1")
  )
  se2 <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(2, nrow = 1, ncol = 1, dimnames = list("g2", "s2"))),
    colData = S4Vectors::DataFrame(PatID = "p2", Timepoint_rough_num = "1", row.names = "s2")
  )
  hc@mae <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = S4Vectors::SimpleList(set1 = se1, set2 = se2)
  )
  hc@config@layer <- S4Vectors::DataFrame(
    layer_id = c("set1", "set2"),
    layer_name = c("RNA layer", "Protein layer")
  )

  means_calls <- list()
  cluster_calls <- list()
  cap_calls <- list()

  testthat::local_mocked_bindings(
    hc_longitudinal_module_means = function(hc,
                                            donor_col,
                                            time_col,
                                            layer,
                                            group_col,
                                            use_module_labels,
                                            time_levels,
                                            impute_missing,
                                            slot_name,
                                            value_label = NULL) {
      means_calls[[length(means_calls) + 1L]] <<- list(
        layer = layer,
        slot_name = slot_name
      )
      sat <- as.list(hc@satellite)
      sat[[slot_name]] <- list(layer_id = layer)
      hc@satellite <- S4Vectors::SimpleList(sat)
      hc
    },
    .hc_run_direct_step1_exact = function(hc, means_slot, output_slot, ...) {
      sat <- as.list(hc@satellite)
      sat[[output_slot]] <- list(
        module_cluster_score = paste0("score_", output_slot),
        module_cluster_best_k = 2L
      )
      hc@satellite <- S4Vectors::SimpleList(sat)
      hc
    },
    hc_plot_longitudinal_module_means = function(hc, slot_name, save_pdf, file_prefix, ...) {
      cluster_calls[[length(cluster_calls) + 1L]] <<- list(
        type = "means",
        slot_name = slot_name,
        file_prefix = file_prefix
      )
      list(module_means = paste(slot_name, file_prefix, sep = "::"))
    },
    hc_plot_longitudinal_module_clusters = function(hc, slot_name, save_pdf, file_prefix, ...) {
      cluster_calls[[length(cluster_calls) + 1L]] <<- list(
        type = "clusters",
        slot_name = slot_name,
        file_prefix = file_prefix
      )
      list(
        module_cluster_waves = paste(slot_name, file_prefix, sep = "::"),
        module_cluster_heatmap = list(
          data = data.frame(donor = factor("d1", levels = "d1"), stringsAsFactors = TRUE)
        )
      )
    },
    hc_plot_longitudinal_cap = function(hc, slot_name, save_pdf, file_prefix, show_values, donor_order = NULL, ...) {
      cap_calls[[length(cap_calls) + 1L]] <<- list(
        slot_name = slot_name,
        file_prefix = file_prefix,
        donor_order = donor_order
      )
      list(cap_heatmap = paste(slot_name, file_prefix, sep = "::"))
    },
    .package = "hcocena"
  )

  out <- hcocena::hc_longitudinal_step1_module_donor(
    hc,
    donor_col = "PatID",
    time_col = "Timepoint_rough_num",
    layer = "all",
    time_levels = c("1"),
    k = 2,
    method = "kmeans",
    nstart = 1,
    cap_runs = 1,
    impute = FALSE,
    ntree = 10,
    min_cluster_fraction = 0.1,
    score_method = "calinski_harabasz",
    scale_features = FALSE,
    seed = 42
  )

  expect_named(out$plots, c("set1", "set2"))
  expect_named(out$diagnostics, c("set1", "set2"))
  expect_equal(out$layer_info$means_slot, c("longitudinal_module_means_set1", "longitudinal_module_means_set2"))
  expect_equal(out$layer_info$output_slot, c("longitudinal_endotypes_set1", "longitudinal_endotypes_set2"))
  expect_equal(vapply(means_calls, `[[`, character(1), "layer"), c("set1", "set2"))
  expect_equal(vapply(means_calls, `[[`, character(1), "slot_name"), c("longitudinal_module_means_set1", "longitudinal_module_means_set2"))
  expect_true(any(vapply(cluster_calls, function(x) identical(x$file_prefix, "Longitudinal_ModuleMeans_RNA_layer"), logical(1))))
  expect_true(any(vapply(cluster_calls, function(x) identical(x$file_prefix, "Longitudinal_ModuleClusters_Protein_layer"), logical(1))))
  expect_equal(vapply(cap_calls, `[[`, character(1), "file_prefix"), c("Longitudinal_CAP_RNA_layer", "Longitudinal_CAP_Protein_layer"))
  expect_equal(out$diagnostics$set1$module_cluster_score, "score_longitudinal_endotypes_set1")
  expect_equal(out$diagnostics$set2$module_cluster_best_k, 2L)
})


test_that("regression: longitudinal direct step1 can loop all layers with unique slots and file prefixes", {
  hc <- hc_init()
  se1 <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(1, nrow = 1, ncol = 1, dimnames = list("g1", "s1"))),
    colData = S4Vectors::DataFrame(PatID = "p1", Timepoint_rough_num = "1", row.names = "s1")
  )
  se2 <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(2, nrow = 1, ncol = 1, dimnames = list("g2", "s2"))),
    colData = S4Vectors::DataFrame(PatID = "p2", Timepoint_rough_num = "1", row.names = "s2")
  )
  hc@mae <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = S4Vectors::SimpleList(set1 = se1, set2 = se2)
  )
  hc@config@layer <- S4Vectors::DataFrame(
    layer_id = c("set1", "set2"),
    layer_name = c("RNA layer", "Protein layer")
  )

  means_calls <- list()
  cluster_calls <- list()
  cap_calls <- list()

  testthat::local_mocked_bindings(
    hc_longitudinal_module_means = function(hc,
                                            donor_col,
                                            time_col,
                                            layer,
                                            group_col,
                                            use_module_labels,
                                            time_levels,
                                            impute_missing,
                                            slot_name,
                                            value_label = NULL) {
      means_calls[[length(means_calls) + 1L]] <<- list(
        layer = layer,
        slot_name = slot_name
      )
      sat <- as.list(hc@satellite)
      sat[[slot_name]] <- list(layer_id = layer)
      hc@satellite <- S4Vectors::SimpleList(sat)
      hc
    },
    .hc_run_direct_step1_exact = function(hc, means_slot, output_slot, ...) {
      sat <- as.list(hc@satellite)
      sat[[output_slot]] <- list(
        module_cluster_score = paste0("direct_score_", output_slot),
        module_cluster_best_k = 2L,
        direct_module_filtering = data.frame(module = "M1", donors_used = 2, stringsAsFactors = FALSE)
      )
      hc@satellite <- S4Vectors::SimpleList(sat)
      hc
    },
    hc_plot_longitudinal_module_means = function(hc, slot_name, save_pdf, file_prefix, ...) {
      cluster_calls[[length(cluster_calls) + 1L]] <<- list(
        type = "means",
        slot_name = slot_name,
        file_prefix = file_prefix
      )
      list(module_means = paste(slot_name, file_prefix, sep = "::"))
    },
    hc_plot_longitudinal_module_clusters = function(hc, slot_name, save_pdf, file_prefix, ...) {
      cluster_calls[[length(cluster_calls) + 1L]] <<- list(
        type = "clusters",
        slot_name = slot_name,
        file_prefix = file_prefix
      )
      list(
        module_cluster_waves = paste(slot_name, file_prefix, sep = "::"),
        module_cluster_heatmap = list(
          data = data.frame(donor = factor("d1", levels = "d1"), stringsAsFactors = TRUE)
        )
      )
    },
    hc_plot_longitudinal_cap = function(hc, slot_name, save_pdf, file_prefix, show_values, donor_order = NULL, ...) {
      cap_calls[[length(cap_calls) + 1L]] <<- list(
        slot_name = slot_name,
        file_prefix = file_prefix,
        donor_order = donor_order
      )
      list(cap_heatmap = paste(slot_name, file_prefix, sep = "::"))
    },
    .package = "hcocena"
  )

  out <- hcocena::hc_longitudinal_step1_module_donor_direct(
    hc,
    donor_col = "PatID",
    time_col = "Timepoint_rough_num",
    layer = "all",
    time_levels = c("1"),
    k = 2,
    method = "kmeans",
    nstart = 5,
    cap_runs = 3,
    impute = FALSE,
    ntree = 10,
    min_cluster_fraction = 0.1,
    score_method = "calinski_harabasz",
    seed = 42
  )

  expect_named(out$plots, c("set1", "set2"))
  expect_named(out$diagnostics, c("set1", "set2"))
  expect_equal(out$layer_info$means_slot, c("longitudinal_module_means_direct_set1", "longitudinal_module_means_direct_set2"))
  expect_equal(out$layer_info$output_slot, c("longitudinal_endotypes_direct_set1", "longitudinal_endotypes_direct_set2"))
  expect_equal(vapply(means_calls, `[[`, character(1), "layer"), c("set1", "set2"))
  expect_equal(vapply(means_calls, `[[`, character(1), "slot_name"), c("longitudinal_module_means_direct_set1", "longitudinal_module_means_direct_set2"))
  expect_true(any(vapply(cluster_calls, function(x) identical(x$file_prefix, "Longitudinal_Direct_ModuleMeans_RNA_layer"), logical(1))))
  expect_true(any(vapply(cluster_calls, function(x) identical(x$file_prefix, "Longitudinal_Direct_ModuleClusters_Protein_layer"), logical(1))))
  expect_equal(vapply(cap_calls, `[[`, character(1), "file_prefix"), c("Longitudinal_Direct_CAP_RNA_layer", "Longitudinal_Direct_CAP_Protein_layer"))
  expect_equal(out$diagnostics$set1$module_cluster_score, "direct_score_longitudinal_endotypes_direct_set1")
  expect_equal(out$diagnostics$set2$module_cluster_best_k, 2L)
})


test_that("regression: longitudinal step1 keeps per-layer nesting for explicit layer='all' with one layer", {
  hc <- hc_init()
  se1 <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(1, nrow = 1, ncol = 1, dimnames = list("g1", "s1"))),
    colData = S4Vectors::DataFrame(PatID = "p1", Timepoint_rough_num = "1", row.names = "s1")
  )
  hc@mae <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = S4Vectors::SimpleList(set1 = se1)
  )
  hc@config@layer <- S4Vectors::DataFrame(
    layer_id = "set1",
    layer_name = "RNA layer"
  )

  testthat::local_mocked_bindings(
    hc_longitudinal_module_means = function(hc,
                                            donor_col,
                                            time_col,
                                            layer,
                                            group_col,
                                            use_module_labels,
                                            time_levels,
                                            impute_missing,
                                            slot_name,
                                            value_label = NULL) {
      sat <- as.list(hc@satellite)
      sat[[slot_name]] <- list(layer_id = layer)
      hc@satellite <- S4Vectors::SimpleList(sat)
      hc
    },
    .hc_run_direct_step1_exact = function(hc, means_slot, output_slot, ...) {
      sat <- as.list(hc@satellite)
      sat[[output_slot]] <- list(
        module_cluster_score = paste0("score_", output_slot),
        module_cluster_best_k = 2L,
        direct_module_filtering = data.frame(module = "M1", donors_used = 2, stringsAsFactors = FALSE)
      )
      hc@satellite <- S4Vectors::SimpleList(sat)
      hc
    },
    hc_plot_longitudinal_module_means = function(hc, slot_name, save_pdf, file_prefix, ...) {
      list(module_means = paste(slot_name, file_prefix, sep = "::"))
    },
    hc_plot_longitudinal_module_clusters = function(hc, slot_name, save_pdf, file_prefix, ...) {
      list(
        module_cluster_waves = paste(slot_name, file_prefix, sep = "::"),
        module_cluster_heatmap = list(
          data = data.frame(donor = factor("d1", levels = "d1"), stringsAsFactors = TRUE)
        )
      )
    },
    hc_plot_longitudinal_cap = function(hc, slot_name, save_pdf, file_prefix, show_values, donor_order = NULL, ...) {
      list(cap_heatmap = paste(slot_name, file_prefix, sep = "::"))
    },
    .package = "hcocena"
  )

  out <- hcocena::hc_longitudinal_step1_module_donor(
    hc,
    donor_col = "PatID",
    time_col = "Timepoint_rough_num",
    layer = "all",
    time_levels = c("1"),
    k = 2,
    method = "kmeans",
    nstart = 1,
    cap_runs = 1,
    impute = FALSE,
    ntree = 10,
    min_cluster_fraction = 0.1,
    score_method = "calinski_harabasz",
    scale_features = FALSE,
    seed = 42
  )

  expect_named(out$plots, "set1")
  expect_named(out$diagnostics, "set1")
  expect_equal(out$layer_info$means_slot, "longitudinal_module_means_set1")
  expect_equal(out$layer_info$output_slot, "longitudinal_endotypes_set1")
  expect_equal(out$plots$set1$module_means_waves, "longitudinal_module_means_set1::Longitudinal_ModuleMeans_RNA_layer")
  expect_true(is.list(out$plots$set1$module_cluster_heatmap))
  expect_equal(out$plots$set1$cap_heatmap, "longitudinal_endotypes_set1::Longitudinal_CAP_RNA_layer")
})


test_that("regression: longitudinal direct workflow nests the three quick steps", {
  hc <- hc_init()

  testthat::local_mocked_bindings(
    hc_longitudinal_step1_module_donor_direct = function(hc, ..., output_slot = "longitudinal_endotypes_direct") {
      sat <- as.list(hc@satellite)
      sat[[output_slot]] <- list(cap_matrix = matrix(1, nrow = 1, dimnames = list("d1", "M1__1")))
      hc@satellite <- S4Vectors::SimpleList(sat)
      list(
        hc = hc,
        plots = list(module_cluster_heatmap = "step1_plot"),
        diagnostics = list(module_cluster_best_k = data.frame(module = "M1", best_k = 2)),
        layer_info = data.frame(layer_id = "set1", output_slot = output_slot, stringsAsFactors = FALSE)
      )
    },
    hc_longitudinal_step2_meta_clustering = function(hc, slot_name, ...) {
      list(
        hc = hc,
        plots = list(pca = paste0(slot_name, "::step2")),
        diagnostics = list(cross_tab = paste0(slot_name, "::diag")),
        slot_info = data.frame(slot_name = slot_name, stringsAsFactors = FALSE)
      )
    },
    hc_longitudinal_step3_meta_module_trajectories = function(hc, slot_name, ...) {
      list(
        hc = hc,
        plots = list(meta_module_waves = paste0(slot_name, "::step3")),
        slot_info = data.frame(slot_name = slot_name, stringsAsFactors = FALSE)
      )
    },
    .package = "hcocena"
  )

  out <- hcocena::hc_longitudinal_workflow_direct(
    hc,
    donor_col = "PatID",
    time_col = "Timepoint_rough_num",
    time_levels = c("1")
  )

  expect_equal(out$plots$step1_module_donor$module_cluster_heatmap, "step1_plot")
  expect_equal(out$plots$step2_meta_clustering$pca, "longitudinal_endotypes_direct::step2")
  expect_equal(out$plots$step3_meta_module_trajectories$meta_module_waves, "longitudinal_endotypes_direct::step3")
  expect_equal(out$slot_info$step2$slot_name, "longitudinal_endotypes_direct")
  expect_equal(out$slot_info$step3$slot_name, "longitudinal_endotypes_direct")
})


test_that("regression: direct longitudinal step1 exact produces module clusters and CAP", {
  hc <- hc_init()
  sat <- list(
    longitudinal_module_means_direct = list(
      donor_col = "Subject",
      time_col = "Time_token",
      group_col = NULL,
      time_levels = c("T1", "T2", "T3", "T4", "T5", "T6"),
      module_lookup = data.frame(
        module = c("M1", "M2"),
        module_color = c("#1f77b4", "#d62728"),
        stringsAsFactors = FALSE
      ),
      donor_time_module = data.frame(
        donor = rep(paste0("D", 1:6), each = 12),
        module = rep(rep(c("M1", "M2"), each = 6), times = 6),
        time = rep(c("T1", "T2", "T3", "T4", "T5", "T6"), times = 12),
        value = c(
          1.0, 2.0, 3.0, 4.0, 4.4, 4.7, 4.5, 3.7, 2.6, 1.5, 1.0, 0.7,
          1.1, 2.1, 3.0, 4.2, 4.1, 4.4, 4.3, 3.8, 2.8, 1.7, 1.2, 0.8,
          0.9, 2.0, 3.2, 4.1, 4.0, 4.6, 4.7, 3.9, 2.7, 1.6, 1.1, 0.9,
          4.2, 3.1, 2.0, 1.0, 0.8, 0.6, 0.8, 1.9, 3.0, 4.0, 4.4, 4.8,
          4.0, 3.0, 1.9, 1.2, 1.0, 0.8, 1.0, 2.1, 3.1, 4.1, 4.5, 4.9,
          4.1, 2.9, 2.1, 1.1, 0.9, 0.7, 0.9, 2.0, 3.2, 4.2, 4.6, 5.0
        ),
        stringsAsFactors = FALSE
      ),
      donor_feature_matrix = matrix(
        0,
        nrow = 6,
        ncol = 12,
        dimnames = list(
          paste0("D", 1:6),
          paste0(rep(c("M1", "M2"), each = 6), "__", rep(c("T1", "T2", "T3", "T4", "T5", "T6"), times = 2))
        )
      ),
      value_label = "Scaled mean VST",
      impute_missing = "none"
    )
  )
  hc@satellite <- S4Vectors::SimpleList(sat)

  hc2 <- hcocena:::.hc_run_direct_step1_exact(
    hc,
    means_slot = "longitudinal_module_means_direct",
    output_slot = "longitudinal_endotypes_direct",
    k = 2:3,
    method = "kmeans",
    nstart = 5,
    cap_runs = 4,
    impute = FALSE,
    min_cluster_fraction = 0.1,
    score_method = "calinski_harabasz",
    seed = 42
  )

  obj <- hc2@satellite[["longitudinal_endotypes_direct"]]
  expect_equal(as.integer(obj$module_cluster_best_k$best_k), c(2L, 2L))
  expect_equal(dim(obj$cap_matrix), c(6L, 4L))
  expect_equal(nrow(obj$module_cluster_assignments), 12L)
  expect_equal(nrow(obj$module_cluster_score_table), 4L)
  expect_true(all(rownames(obj$cap_matrix) == paste0("D", 1:6)))
})


test_that("regression: longitudinal step2 can loop suffixed step1 slots with unique prefixes", {
  hc <- hc_init()
  hc@config@layer <- S4Vectors::DataFrame(
    layer_id = c("set1", "set2"),
    layer_name = c("RNA layer", "Protein layer")
  )
  hc@satellite <- S4Vectors::SimpleList(
    longitudinal_endotypes_set1 = list(
      cap_matrix = matrix(
        c(1, 2, 3, 4),
        nrow = 2,
        dimnames = list(c("d1", "d2"), c("m1", "m2"))
      ),
      source_slot = "longitudinal_module_means_set1"
    ),
    longitudinal_endotypes_set2 = list(
      cap_matrix = matrix(
        c(5, 6, 7, 8),
        nrow = 2,
        dimnames = list(c("d3", "d4"), c("m1", "m2"))
      ),
      source_slot = "longitudinal_module_means_set2"
    )
  )

  run_calls <- list()
  plot_calls <- list()

  testthat::local_mocked_bindings(
    .hc_run_longitudinal_step2_graph = function(hc, slot_name, ...) {
      run_calls[[length(run_calls) + 1L]] <<- slot_name
      sat <- as.list(hc@satellite)
      sat[[slot_name]]$meta_cluster <- data.frame(
        donor = c("d1", "d2"),
        meta_cluster = c("MC1", "MC2"),
        stringsAsFactors = FALSE
      )
      sat[[slot_name]]$meta_method_comparison <- data.frame(
        method = "graph_leiden",
        stringsAsFactors = FALSE
      )
      sat[[slot_name]]$meta_score_table <- data.frame(
        k = 2L,
        stringsAsFactors = FALSE
      )
      hc@satellite <- S4Vectors::SimpleList(sat)
      hc
    },
    hc_plot_longitudinal_meta_embeddings = function(hc,
                                                    slot_name,
                                                    save_pdf,
                                                    file_prefix,
                                                    show_endotype_crosstab,
                                                    show_cluster_labels,
                                                    save_tables,
                                                    table_format,
                                                    table_detail) {
      plot_calls[[length(plot_calls) + 1L]] <<- list(
        slot_name = slot_name,
        file_prefix = file_prefix
      )
      list(
        pca = paste(slot_name, "pca", sep = "::"),
        umap = paste(slot_name, "umap", sep = "::"),
        cross_tab = paste(slot_name, "cross", sep = "::"),
        tables = list(Method_Clusters = paste(slot_name, "table", sep = "::"))
      )
    },
    .package = "hcocena"
  )

  out <- hcocena::hc_longitudinal_step2_meta_clustering(
    hc,
    slot_name = "longitudinal_endotypes",
    dimensions = 4,
    graph_method = "knn",
    knn_method = "annoy",
    graph_k = 7,
    resolution = 0.3,
    leiden_method = "RBConfigurationVertexPartition"
  )

  expect_named(out$plots, c("set1", "set2"))
  expect_named(out$diagnostics, c("set1", "set2"))
  expect_equal(run_calls, list("longitudinal_endotypes_set1", "longitudinal_endotypes_set2"))
  expect_equal(vapply(plot_calls, `[[`, character(1), "file_prefix"), c("Longitudinal_Meta_RNA_layer", "Longitudinal_Meta_Protein_layer"))
  expect_equal(out$slot_info$slot_name, c("longitudinal_endotypes_set1", "longitudinal_endotypes_set2"))
  expect_equal(out$slot_info$layer_id, c("set1", "set2"))
  expect_equal(out$plots$set1$pca, "longitudinal_endotypes_set1::pca")
  expect_equal(out$diagnostics$set2$cross_tab, "longitudinal_endotypes_set2::cross")
})


test_that("regression: longitudinal step2 prefers suffixed family slots over stale unsuffixed slot", {
  hc <- hc_init()
  hc@satellite <- S4Vectors::SimpleList(
    longitudinal_endotypes = list(
      cap_matrix = matrix(
        c(9, 9, 9, 9),
        nrow = 2,
        dimnames = list(c("dx1", "dx2"), c("m1", "m2"))
      )
    ),
    longitudinal_endotypes_set1 = list(
      cap_matrix = matrix(
        c(1, 2, 3, 4),
        nrow = 2,
        dimnames = list(c("d1", "d2"), c("m1", "m2"))
      ),
      source_slot = "longitudinal_module_means_set1"
    ),
    longitudinal_endotypes_set2 = list(
      cap_matrix = matrix(
        c(5, 6, 7, 8),
        nrow = 2,
        dimnames = list(c("d3", "d4"), c("m1", "m2"))
      ),
      source_slot = "longitudinal_module_means_set2"
    )
  )

  resolved <- hcocena:::.hc_longitudinal_step2_target_slots(hc, slot_name = "longitudinal_endotypes")
  expect_equal(resolved, c("longitudinal_endotypes_set1", "longitudinal_endotypes_set2"))
})


test_that("regression: longitudinal step2 keeps per-slot nesting for family slot_name with one matching slot", {
  hc <- hc_init()
  hc@config@layer <- S4Vectors::DataFrame(
    layer_id = "set1",
    layer_name = "RNA layer"
  )
  hc@satellite <- S4Vectors::SimpleList(list(
    longitudinal_endotypes_set1 = list(
      cap_matrix = matrix(
        c(1, 2, 3, 4),
        nrow = 2,
        dimnames = list(c("d1", "d2"), c("m1", "m2"))
      ),
      source_slot = "longitudinal_module_means_set1"
    )
  ))

  testthat::local_mocked_bindings(
    .hc_run_longitudinal_step2_graph = function(hc, slot_name, ...) {
      sat <- as.list(hc@satellite)
      sat[[slot_name]]$meta_cluster <- data.frame(donor = "d1", meta_cluster = "MC1", stringsAsFactors = FALSE)
      sat[[slot_name]]$meta_method_comparison <- data.frame(method = "graph_leiden", stringsAsFactors = FALSE)
      sat[[slot_name]]$meta_score_table <- data.frame(k = 2L, stringsAsFactors = FALSE)
      hc@satellite <- S4Vectors::SimpleList(sat)
      hc
    },
    hc_plot_longitudinal_meta_embeddings = function(hc,
                                                    slot_name,
                                                    save_pdf,
                                                    file_prefix,
                                                    show_endotype_crosstab,
                                                    show_cluster_labels,
                                                    save_tables,
                                                    table_format,
                                                    table_detail) {
      list(
        pca = paste(slot_name, "pca", sep = "::"),
        umap = paste(slot_name, "umap", sep = "::"),
        cross_tab = paste(slot_name, "cross", sep = "::"),
        tables = list(Method_Clusters = paste(slot_name, "table", sep = "::"))
      )
    },
    .package = "hcocena"
  )

  out <- hcocena::hc_longitudinal_step2_meta_clustering(
    hc,
    slot_name = "longitudinal_endotypes",
    dimensions = 4,
    graph_method = "knn",
    knn_method = "annoy",
    graph_k = 7,
    resolution = 0.3,
    leiden_method = "RBConfigurationVertexPartition"
  )

  expect_named(out$plots, "set1")
  expect_named(out$diagnostics, "set1")
  expect_equal(out$slot_info$slot_name, "longitudinal_endotypes_set1")
  expect_equal(out$slot_info$layer_id, "set1")
  expect_equal(out$plots$set1$pca, "longitudinal_endotypes_set1::pca")
  expect_equal(out$plots$set1$umap, "longitudinal_endotypes_set1::umap")
})


test_that("regression: longitudinal step3 can loop suffixed step2 slots with unique prefixes", {
  hc <- hc_init()
  hc@config@layer <- S4Vectors::DataFrame(
    layer_id = c("set1", "set2"),
    layer_name = c("RNA layer", "Protein layer")
  )
  hc@satellite <- S4Vectors::SimpleList(
    longitudinal_endotypes = list(meta_cluster = data.frame(donor = "stale", stringsAsFactors = FALSE)),
    longitudinal_endotypes_set1 = list(
      meta_cluster = data.frame(donor = c("d1", "d2"), meta_cluster = c("MC1", "MC2"), stringsAsFactors = FALSE),
      source_slot = "longitudinal_module_means_set1"
    ),
    longitudinal_endotypes_set2 = list(
      meta_cluster = data.frame(donor = c("d3", "d4"), meta_cluster = c("MC1", "MC2"), stringsAsFactors = FALSE),
      source_slot = "longitudinal_module_means_set2"
    )
  )

  plot_calls <- list()

  testthat::local_mocked_bindings(
    hc_plot_longitudinal_meta_module_waves = function(hc,
                                                      slot_name,
                                                      save_pdf,
                                                      file_prefix,
                                                      facet_ncol,
                                                      free_y,
                                                      square_panels,
                                                      value_mode,
                                                      value_range,
                                                      save_width,
                                                      save_height) {
      plot_calls[[length(plot_calls) + 1L]] <<- list(
        slot_name = slot_name,
        file_prefix = file_prefix,
        value_mode = value_mode
      )
      list(meta_module_waves = paste(slot_name, file_prefix, sep = "::"))
    },
    .package = "hcocena"
  )

  out <- hcocena::hc_longitudinal_step3_meta_module_trajectories(
    hc,
    slot_name = "longitudinal_endotypes",
    facet_ncol = 4,
    free_y = FALSE,
    square_panels = TRUE,
    value_mode = "scaled_mean_vst",
    value_range = c(-2, 2)
  )

  expect_named(out$plots, c("set1", "set2"))
  expect_equal(vapply(plot_calls, `[[`, character(1), "slot_name"), c("longitudinal_endotypes_set1", "longitudinal_endotypes_set2"))
  expect_equal(vapply(plot_calls, `[[`, character(1), "file_prefix"), c("Longitudinal_Meta_ModuleWaves_RNA_layer", "Longitudinal_Meta_ModuleWaves_Protein_layer"))
  expect_equal(out$slot_info$slot_name, c("longitudinal_endotypes_set1", "longitudinal_endotypes_set2"))
  expect_equal(out$plots$set1$meta_module_waves, "longitudinal_endotypes_set1::Longitudinal_Meta_ModuleWaves_RNA_layer")
  expect_equal(out$plots$set2$meta_module_waves, "longitudinal_endotypes_set2::Longitudinal_Meta_ModuleWaves_Protein_layer")
})


test_that("regression: longitudinal step3 keeps per-slot nesting for family slot_name with one matching slot", {
  hc <- hc_init()
  hc@config@layer <- S4Vectors::DataFrame(
    layer_id = "set1",
    layer_name = "RNA layer"
  )
  hc@satellite <- S4Vectors::SimpleList(list(
    longitudinal_endotypes_set1 = list(
      meta_cluster = data.frame(donor = c("d1", "d2"), meta_cluster = c("MC1", "MC2"), stringsAsFactors = FALSE),
      source_slot = "longitudinal_module_means_set1"
    )
  ))

  testthat::local_mocked_bindings(
    hc_plot_longitudinal_meta_module_waves = function(hc,
                                                      slot_name,
                                                      save_pdf,
                                                      file_prefix,
                                                      facet_ncol,
                                                      free_y,
                                                      square_panels,
                                                      value_mode,
                                                      value_range,
                                                      save_width,
                                                      save_height) {
      list(meta_module_waves = paste(slot_name, file_prefix, sep = "::"))
    },
    .package = "hcocena"
  )

  out <- hcocena::hc_longitudinal_step3_meta_module_trajectories(
    hc,
    slot_name = "longitudinal_endotypes",
    facet_ncol = 4,
    free_y = FALSE,
    square_panels = TRUE,
    value_mode = "scaled_mean_vst",
    value_range = c(-2, 2)
  )

  expect_named(out$plots, "set1")
  expect_equal(out$slot_info$slot_name, "longitudinal_endotypes_set1")
  expect_equal(out$slot_info$layer_id, "set1")
  expect_equal(out$plots$set1$meta_module_waves, "longitudinal_endotypes_set1::Longitudinal_Meta_ModuleWaves_RNA_layer")
})


test_that("regression: longitudinal enrichment meta waves resolve family slot_name with one matching slot", {
  hc <- hc_init()
  hc@config@layer <- S4Vectors::DataFrame(
    layer_id = "set1",
    layer_name = "RNA layer"
  )
  hc@satellite <- S4Vectors::SimpleList(list(
    longitudinal_endotypes_set1 = list(
      meta_cluster = data.frame(
        donor = c("d1", "d2"),
        meta_cluster = c("MC1", "MC1"),
        stringsAsFactors = FALSE
      ),
      source_slot = "longitudinal_module_means_set1"
    )
  ))

  testthat::local_mocked_bindings(
    hc_plot_longitudinal_enrichment_waves = function(hc, slot_name, ...) {
      expect_equal(slot_name, "longitudinal_endotypes_set1")
      list(
        plots = list(Hallmark = "base_plot"),
        top_terms = list(Hallmark = data.frame(
          module = "M1",
          module_color = "red",
          term = "IFN signaling",
          rank = 1,
          qvalue = 0.01,
          stringsAsFactors = FALSE
        )),
        donor_trajectories = list(Hallmark = data.frame(
          donor = c("d1", "d2"),
          module = c("M1", "M1"),
          term = c("IFN signaling", "IFN signaling"),
          time = c("T1", "T1"),
          score = c(1, 2),
          stringsAsFactors = FALSE
        )),
        mean_trajectories = list(),
        score_method_used = list(Hallmark = "rank_mean")
      )
    },
    .package = "hcocena"
  )

  out <- hcocena::hc_plot_longitudinal_enrichment_meta_waves(
    hc,
    slot_name = "longitudinal_endotypes",
    databases = "Hallmark",
    top = 1,
    score_method = "rank_mean",
    show_donor_lines = FALSE,
    save_pdf = FALSE,
    export_excel = FALSE,
    donor_col = "PatID",
    time_col = "Timepoint_rough_num",
    time_levels = "T1"
  )

  expect_named(out$plots, "Hallmark")
  expect_s3_class(out$plots$Hallmark, "ggplot")
  expect_named(out$top_terms, "Hallmark")
})


test_that("regression: meta-time grouping uses per-layer longitudinal slots and returns regrouped heatmap order", {
  hc <- hc_init()
  se1 <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(c(1, 2), nrow = 1, dimnames = list("g1", c("s1", "s2")))),
    colData = S4Vectors::DataFrame(
      PatID = c("p1", "p2"),
      Timepoint_rough = c("T1", "T2"),
      row.names = c("s1", "s2")
    )
  )
  se2 <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(c(3, 4), nrow = 1, dimnames = list("g2", c("s3", "s4")))),
    colData = S4Vectors::DataFrame(
      PatID = c("p3", "p4"),
      Timepoint_rough = c("T1", "T2"),
      row.names = c("s3", "s4")
    )
  )
  hc@mae <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = S4Vectors::SimpleList(set1 = se1, set2 = se2)
  )
  hc@config@layer <- S4Vectors::DataFrame(
    layer_id = c("set1", "set2"),
    layer_name = c("pretm", "pretm2")
  )
  hc@satellite <- S4Vectors::SimpleList(
    longitudinal_endotypes = list(meta_cluster = data.frame(donor = "stale", meta_cluster = "MC0", stringsAsFactors = FALSE)),
    longitudinal_endotypes_set1 = list(
      meta_cluster = data.frame(donor = c("p1", "p2"), meta_cluster = c("MC1", "MC2"), stringsAsFactors = FALSE),
      source_slot = "longitudinal_module_means_set1"
    ),
    longitudinal_endotypes_set2 = list(
      meta_cluster = data.frame(donor = c("p3", "p4"), meta_cluster = c("MC3", "MC4"), stringsAsFactors = FALSE),
      source_slot = "longitudinal_module_means_set2"
    )
  )

  hc <- hcocena::hc_add_meta_time_grouping(
    hc,
    donor_col = "PatID",
    time_col = "Timepoint_rough",
    grouping_col = "meta_cluster_time",
    slot_name = "longitudinal_endotypes"
  )

  anno1 <- as.data.frame(SummarizedExperiment::colData(MultiAssayExperiment::experiments(hc@mae)[["set1"]]), stringsAsFactors = FALSE)
  anno2 <- as.data.frame(SummarizedExperiment::colData(MultiAssayExperiment::experiments(hc@mae)[["set2"]]), stringsAsFactors = FALSE)
  expect_equal(as.character(anno1$meta_cluster_time), c("MC1__T1", "MC2__T2"))
  expect_equal(as.character(anno2$meta_cluster_time), c("MC3__T1", "MC4__T2"))
  expect_equal(
    hcocena::hc_get_meta_time_col_order(hc, slot_name = "longitudinal_endotypes", layer = 1),
    c("MC1__T1_pretm", "MC2__T2_pretm")
  )
  expect_equal(
    hcocena::hc_get_meta_time_col_order(hc, slot_name = "longitudinal_endotypes"),
    c("MC1__T1_pretm", "MC2__T2_pretm", "MC3__T1_pretm2", "MC4__T2_pretm2")
  )
})


test_that("regression: safe deep clone falls back instead of aborting", {
  clone_fun <- get(".hc_safe_deep_clone", asNamespace("hcocena"))

  x <- list(a = 1, b = list(2))
  y <- clone_fun(x, context = "test object")
  y$b[[1]] <- 3
  expect_equal(x$b[[1]], 2)

  env <- new.env(parent = emptyenv())
  env$a <- 1
  testthat::local_mocked_bindings(
    serialize = function(...) stop("mock oom"),
    .package = "base"
  )
  expect_warning(
    out <- clone_fun(env, context = "test object"),
    "Could not deep-clone test object"
  )
  expect_identical(out, env)

  expect_warning(
    out_big <- clone_fun(x, context = "large test object", max_bytes = 1),
    "Skipping deep clone of large test object"
  )
  expect_identical(out_big, x)
})


test_that("regression: cutoff stats match legacy igraph behavior on small graphs", {
  rsquaredfun <- get("rsquaredfun", asNamespace("hcocena"))
  union_find <- get(".hc_union_find_components", asNamespace("hcocena"))

  graph_df <- data.frame(
    V1 = c("g1", "g1", "g2", "g4"),
    V2 = c("g2", "g3", "g3", "g5"),
    rval = c(0.91, 0.88, 0.9, 0.95),
    stringsAsFactors = FALSE
  )

  ref_graph <- igraph::graph_from_data_frame(graph_df, directed = FALSE, vertices = NULL)
  ref_components <- igraph::components(ref_graph)
  keep_components <- which(ref_components$csize >= 3)
  gene_to_comp <- data.frame(
    gene = names(ref_components$membership),
    component = ref_components$membership,
    stringsAsFactors = FALSE
  )
  nodes_to_remove <- dplyr::filter(gene_to_comp, !component %in% keep_components) %>%
    dplyr::pull(gene)
  ref_graph <- igraph::delete_vertices(ref_graph, nodes_to_remove)
  ref_degree <- igraph::degree(ref_graph, mode = "all")
  ref_dd <- igraph::degree_distribution(ref_graph, mode = "all", cumulative = FALSE)
  ref_prob <- ref_dd[-1]
  ref_deg_vals <- seq_len(max(ref_degree))
  nonzero <- which(ref_prob != 0)
  ref_prob <- ref_prob[nonzero]
  ref_deg_vals <- ref_deg_vals[nonzero]
  ref_r2 <- summary(stats::lm(log(ref_prob) ~ log(ref_deg_vals)))$r.squared

  uf <- union_find(from_idx = c(1L, 1L, 2L, 4L), to_idx = c(2L, 3L, 3L, 5L), n_vertices = 5L)
  expect_equal(sort(uf$csize), c(2L, 3L))
  expect_equal(uf$csize[uf$membership], c(3L, 3L, 3L, 2L, 2L))

  out <- rsquaredfun(
    graph_df = graph_df,
    cutoff = 0.88,
    print.all.plots = FALSE,
    min_nodes = 3,
    x = 1
  )

  expect_equal(out$no_of_networks, 1)
  expect_equal(out$no_nodes, 3)
  expect_equal(out$no_edges, 3)
  expect_equal(out$degree[[1]], ref_deg_vals)
  expect_equal(out$Probs[[1]], ref_prob)
  expect_equal(out$R.squared[[1]], ref_r2)
})


test_that("regression: longitudinal top terms are compacted after filtering", {
  rank_lookup <- get(".hc_longitudinal_rank_lookup", asNamespace("hcocena"))
  panel_label <- get(".hc_longitudinal_panel_label", asNamespace("hcocena"))

  top_df <- data.frame(
    module = c("M1", "M1", "M1", "M2", "M2"),
    term = c("KEGG_T1", "KEGG_T2", "KEGG_T3", "HALLMARK_U1", "HALLMARK_U2"),
    module_color = c("red", "red", "red", "blue", "blue"),
    rank = c(1, 2, 3, 1, 2),
    qvalue = c(0.01, 0.02, 0.03, 0.01, 0.02),
    stringsAsFactors = FALSE
  )
  valid_terms <- data.frame(
    module = c("M1", "M1", "M2"),
    term = c("KEGG_T1", "KEGG_T3", "HALLMARK_U2"),
    stringsAsFactors = FALSE
  )

  out <- rank_lookup(
    top_df = top_df,
    wrap_term = function(x) base::gsub("_", " ", x, fixed = TRUE),
    valid_terms = valid_terms
  )

  expect_equal(out[out$module == "M1", "term"], c("KEGG_T1", "KEGG_T3"))
  expect_equal(out[out$module == "M1", "term_rank"], c(1L, 2L))
  expect_equal(out[out$module == "M1", "term_rank_label"], c("Top 1", "Top 2"))
  expect_false(any(base::is.na(out$term_label)))
  expect_false(any(!base::nzchar(base::trimws(out$term_label))))

  fallback <- panel_label(module = NA_character_, rank_label = NA_character_, term_display = NA_character_)
  expect_identical(fallback, "Module: Top\nTerm unavailable")
})


test_that("regression: longitudinal enrichment uses observed module label/color pairs", {
  harmonize <- get(".hc_harmonize_module_lookup_with_enrichment", asNamespace("hcocena"))

  module_lookup <- data.frame(
    module = c("M1", "M2", "M3"),
    module_color = c("gold", "lightblue", "wheat"),
    stringsAsFactors = FALSE
  )
  enrich_df <- data.frame(
    cluster = c("gold", "darkgreen", "cyan", "wheat"),
    module_label = c("M1", "M7", "M10", ""),
    stringsAsFactors = FALSE
  )

  out <- harmonize(module_lookup = module_lookup, enrich_df = enrich_df)

  expect_equal(out$module_color[match("M1", out$module)], "gold")
  expect_equal(out$module_color[match("M7", out$module)], "darkgreen")
  expect_equal(out$module_color[match("M10", out$module)], "cyan")
  expect_equal(out$module[match("wheat", out$module_color)], "M3")
})


test_that("regression: longitudinal plotting drops panel labels outside final facet levels", {
  align_panels <- get(".hc_align_longitudinal_panel_factors", asNamespace("hcocena"))

  df <- data.frame(
    panel_label = c("M1: Top 1\nA", "M1: Top 2\nB", NA, "NA"),
    term_display = c("A", "B", "A", "B"),
    value = 1:4,
    stringsAsFactors = FALSE
  )

  out <- align_panels(
    df = df,
    panel_levels = c("M1: Top 1\nA"),
    term_display_levels = c("A")
  )

  expect_equal(base::as.character(out$panel_label), "M1: Top 1\nA")
  expect_equal(base::as.character(out$term_display), "A")
  expect_equal(out$value, 1L)
})


test_that("regression: longitudinal enrichment meta waves can flatten multi-slot outputs", {
  hc <- hc_init()
  hc@config@layer <- S4Vectors::DataFrame(
    layer_id = c("set1", "set2"),
    layer_name = c("RNA layer", "Protein layer")
  )
  hc@satellite <- S4Vectors::SimpleList(
    longitudinal_endotypes = list(meta_cluster = data.frame(donor = "stale", meta_cluster = "MC0", stringsAsFactors = FALSE)),
    longitudinal_endotypes_set1 = list(
      meta_cluster = data.frame(donor = c("d1", "d2"), meta_cluster = c("MC1", "MC2"), stringsAsFactors = FALSE),
      source_slot = "longitudinal_module_means_set1"
    ),
    longitudinal_endotypes_set2 = list(
      meta_cluster = data.frame(donor = c("d3", "d4"), meta_cluster = c("MC1", "MC2"), stringsAsFactors = FALSE),
      source_slot = "longitudinal_module_means_set2"
    )
  )

  enrichment_calls <- list()

  testthat::local_mocked_bindings(
    hc_plot_longitudinal_enrichment_waves = function(hc,
                                                     slot_name,
                                                     databases,
                                                     top,
                                                     custom_terms,
                                                     term_match,
                                                     enrichment_table,
                                                     score_method,
                                                     score_scale,
                                                     layer,
                                                     donor_col,
                                                     time_col,
                                                     time_levels,
                                                     min_term_genes,
                                                     impute_missing,
                                                     qvalue_max,
                                                     show_donor_lines,
                                                     facet_ncol,
                                                     free_y,
                                                     save_pdf,
                                                     file_prefix,
                                                     save_width,
                                                     save_height,
                                                     export_excel) {
      enrichment_calls[[length(enrichment_calls) + 1L]] <<- list(
        slot_name = slot_name,
        file_prefix = file_prefix
      )
      donor_ids <- if (grepl("set1$", slot_name)) c("d1", "d2") else c("d3", "d4")
      list(
        plots = list(),
        top_terms = list(
          Hallmark = data.frame(
            module = "M1",
            term = "TERM_A",
            module_color = "gold",
            rank = 1,
            qvalue = 0.01,
            stringsAsFactors = FALSE
          )
        ),
        donor_trajectories = list(
          Hallmark = data.frame(
            donor = rep(donor_ids, each = 2),
            module = "M1",
            term = "TERM_A",
            time = rep(c("1", "2"), times = length(donor_ids)),
            score = c(0.1, 0.2, 0.3, 0.4),
            stringsAsFactors = FALSE
          )
        ),
        mean_trajectories = list(),
        score_method_used = list(Hallmark = score_method)
      )
    },
    .package = "hcocena"
  )

  out <- hcocena::hc_plot_longitudinal_enrichment_meta_waves(
    hc,
    slot_name = "longitudinal_endotypes",
    databases = "Hallmark",
    top = 1,
    show_donor_lines = FALSE,
    score_method = "ssgsea",
    score_scale = "z",
    save_pdf = FALSE,
    export_excel = FALSE
  )

  expect_equal(vapply(enrichment_calls, `[[`, character(1), "slot_name"), c("longitudinal_endotypes_set1", "longitudinal_endotypes_set2"))
  expect_equal(vapply(enrichment_calls, `[[`, character(1), "file_prefix"), c("Longitudinal_Enrichment_MetaWaves_RNA_layer", "Longitudinal_Enrichment_MetaWaves_Protein_layer"))
  expect_named(out$results_by_slot, c("set1", "set2"))
  expect_true(all(c("set1__Hallmark", "set2__Hallmark") %in% names(out$plots)))
  expect_true(inherits(out$plots$set1__Hallmark, "ggplot"))
  expect_equal(out$slot_info$slot_name, c("longitudinal_endotypes_set1", "longitudinal_endotypes_set2"))
})


test_that("regression: emmeans summaries are coerced to numeric safely", {
  normalize_emm <- get(".hc_normalize_emmeans_summary", asNamespace("hcocena"))

  df <- data.frame(
    response = c("0.4", "0.7"),
    lower.CL = c("nonEst", "0.2"),
    upper.CL = c("0.9", "1.1"),
    stringsAsFactors = FALSE
  )

  out <- normalize_emm(df)

  expect_equal(out$emmean, c(0.4, 0.7))
  expect_true(is.na(out$lower.CL[[1]]))
  expect_equal(out$lower.CL[[2]], 0.2)
  expect_equal(out$upper.CL, c(0.9, 1.1))

  out$lower.CL <- suppressWarnings(as.numeric(out$lower.CL))
  out$upper.CL <- suppressWarnings(as.numeric(out$upper.CL))
  ci_missing <- !is.finite(out$lower.CL) | !is.finite(out$upper.CL)
  out$lower.CL[ci_missing] <- out$emmean[ci_missing]
  out$upper.CL[ci_missing] <- out$emmean[ci_missing]

  expect_equal(out$lower.CL[[1]], out$emmean[[1]])
  expect_equal(out$upper.CL[[1]], out$emmean[[1]])
})


test_that("regression: heatmap col_order ignores missing columns safely", {
  resolve_cols <- get(".hc_resolve_heatmap_col_order", asNamespace("hcocena"))

  expect_warning(
    out <- resolve_cols(
      mat_cols = c("MC1_T1", "MC1_T2", "MC2_T1"),
      requested_order = c("old_group", "MC2_T1"),
      context = "regrouped cluster heatmap"
    ),
    "Ignoring 1 `col_order` entries"
  )

  expect_equal(out, c("MC2_T1", "MC1_T1", "MC1_T2"))
})


test_that("regression: enrichment panel storage defaults to memory-saving for multi-db runs", {
  resolve_mode <- get(".hc_resolve_panel_storage_mode", asNamespace("hcocena"))
  has_draw_obj <- get(".hc_functional_enrichment_has_draw_object", asNamespace("hcocena"))

  expect_identical(resolve_mode("auto", n_databases = 1L), "always")
  expect_identical(resolve_mode("auto", n_databases = 3L), "never")
  expect_identical(resolve_mode("always", n_databases = 3L), "always")

  expect_true(has_draw_obj(list(p = structure(list(), class = "HeatmapList"))))
  expect_false(has_draw_obj(list(
    p = NULL,
    hc_heatmap = NULL,
    enrichment_plot = NULL,
    panel_objects_stored = FALSE
  )))
})


test_that("regression: large result stores are no longer mirrored across legacy slots", {
  fun_enrich_src <- paste(deparse(get(".hc_functional_enrichment_driver", asNamespace("hcocena"))), collapse = "\n")
  heatmap_src <- paste(deparse(get(".hc_plot_cluster_heatmap_driver", asNamespace("hcocena"))), collapse = "\n")
  heatmap_new_src <- paste(deparse(get("plot_cluster_heatmap_new", asNamespace("hcocena"))), collapse = "\n")
  upstream_src <- paste(deparse(get(".hc_upstream_inference_driver", asNamespace("hcocena"))), collapse = "\n")
  knowledge_src <- paste(deparse(get(".hc_plot_enrichment_upstream_network_driver", asNamespace("hcocena"))), collapse = "\n")
  network_plot_src <- paste(
    deparse(get(".hc_plot_integrated_network_driver", asNamespace("hcocena"))),
    collapse = "\n"
  )
  gfc_network_src <- paste(
    deparse(get(".hc_plot_GFC_network_driver", asNamespace("hcocena"))),
    collapse = "\n"
  )

  expect_false(grepl(
    'hcobject[["satellite_outputs"]][["enrichments"]] <<- hcobject[["integrated_output"]][["enrichments"]]',
    fun_enrich_src,
    fixed = TRUE
  ))
  expect_false(grepl(
    'hcobject[["integrated_output"]][["enrichments"]][["all_enrichments_all_dbs"]] <<-',
    fun_enrich_src,
    fixed = TRUE
  ))
  expect_false(grepl(
    'hcobject[["integrated_output"]][["cluster_calc"]][["module_gene_list"]] <<-',
    heatmap_src,
    fixed = TRUE
  ))
  expect_true(grepl(
    '\\.hc_set_bridge_hcobject_slot\\(c\\("integrated_output", "cluster_calc",\\s*"heatmap_matrix"\\)',
    heatmap_new_src,
    perl = TRUE
  ))
  expect_true(grepl(
    '\\.hc_set_bridge_hcobject_slot\\(c\\("integrated_output", "cluster_calc",\\s*"heatmap_row_order"\\)',
    heatmap_new_src,
    perl = TRUE
  ))
  expect_true(grepl(
    '\\.hc_set_bridge_hcobject_slot\\(c\\("integrated_output", "cluster_calc",\\s*"heatmap_column_order"\\)',
    heatmap_new_src,
    perl = TRUE
  ))
  expect_false(grepl(
    'hcobject[["integrated_output"]][["upstream_inference"]] <<- output',
    upstream_src,
    fixed = TRUE
  ))
  expect_false(grepl(
    'hcobject[["integrated_output"]][["knowledge_network"]] <<- output',
    knowledge_src,
    fixed = TRUE
  ))
  expect_false(grepl(
    'hcobject[["integrated_output"]][["cluster_calc"]][["network_col_by_module"]] <<- network',
    network_plot_src,
    fixed = TRUE
  ))
  expect_false(grepl(
    'hcobject[["integrated_output"]][["cluster_calc"]][["labelled_network"]] <<- network2',
    network_plot_src,
    fixed = TRUE
  ))
  expect_true(grepl(
    'if (isTRUE(store_plot))',
    network_plot_src,
    fixed = TRUE
  ))
  expect_true(grepl(
    'sat[["network_col_by_module"]]',
    gfc_network_src,
    fixed = TRUE
  ))
  expect_true(grepl(
    'network <- hcobject[["integrated_output"]][["merged_net"]]',
    gfc_network_src,
    fixed = TRUE
  ))
})


test_that("regression: duplicate GFC condition names survive S4-to-legacy conversion", {
  to_base_df <- get(".hc_to_base_data_frame_preserve_names", asNamespace("hcocena"))
  cluster_plot_hco <- get(".hc_as_bridge_object_for_cluster_plot", asNamespace("hcocena"))

  dup_df <- data.frame(
    T1 = c(1, 2),
    T1 = c(3, 4),
    Gene = c("g1", "g2"),
    check.names = FALSE
  )
  expect_identical(colnames(to_base_df(S4Vectors::DataFrame(dup_df, check.names = FALSE))), c("T1", "T1", "Gene"))

  hc <- hc_init()
  hc@integration@gfc <- S4Vectors::DataFrame(dup_df, check.names = FALSE)

  legacy_full <- hcocena:::as_hcobject(hc)
  expect_identical(colnames(legacy_full$integrated_output$GFC_all_layers), c("T1", "T1", "Gene"))

  legacy_plot <- cluster_plot_hco(hc)
  expect_identical(colnames(legacy_plot$integrated_output$GFC_all_layers), c("T1", "T1", "Gene"))
})


test_that("regression: duplicate GFC condition names no longer break cluster summaries", {
  cond_names_fun <- get(".hc_gfc_condition_names", asNamespace("hcocena"))
  colmeans_fun <- get(".hc_gfc_colmeans_for_genes", asNamespace("hcocena"))
  cluster_mean_fun <- get("gfc_mean_clustergene", asNamespace("hcocena"))

  gfc_df <- data.frame(
    T1 = c(1, 2),
    T2 = c(3, 4),
    T3 = c(5, 6),
    T1 = c(7, 8),
    T2 = c(9, 10),
    T3 = c(11, 12),
    Gene = c("g1", "g2"),
    check.names = FALSE
  )

  expect_identical(cond_names_fun(gfc_df), c("T1", "T2", "T3", "T1", "T2", "T3"))
  expect_equal(
    unname(colmeans_fun(gfc_df, genes = "g1")),
    c(1, 3, 5, 7, 9, 11)
  )

  out <- cluster_mean_fun(
    rownum = 1,
    cluster_df = data.frame(gene_n = "g1,g2", stringsAsFactors = FALSE),
    gfc_dat = gfc_df
  )
  expect_identical(out$conditions, "T1#T2#T3#T1#T2#T3")
  expect_identical(out$grp_means, "1.5,3.5,5.5,7.5,9.5,11.5")
})


test_that("split_modules resolves numeric inputs by current heatmap row order", {
  resolve_fun <- get(".hc_resolve_modules_for_split", asNamespace("hcocena"))
  order_fun <- get(".hc_split_module_order", asNamespace("hcocena"))

  module_label_map <- c(red = "M-red", blue = "M-blue", green = "M-green")
  module_order <- order_fun(
    cluster_calc = list(heatmap_row_order = c("blue", "green", "red")),
    available_colors = c("red", "blue", "green")
  )

  resolved <- resolve_fun(
    modules = c(1, 3),
    available_colors = c("red", "blue", "green"),
    module_label_map = module_label_map,
    module_order = module_order
  )

  expect_equal(module_order, c("blue", "green", "red"))
  expect_equal(resolved$target_colors, c("blue", "red"))
  expect_equal(resolved$resolved_labels, c("M-blue", "M-red"))
  expect_equal(resolved$resolution_table$resolved_index, c(1L, 3L))
  expect_equal(resolved$resolution_table$resolved_color, c("blue", "red"))

  bad <- resolve_fun(
    modules = 1.5,
    available_colors = c("red", "blue", "green"),
    module_label_map = module_label_map,
    module_order = module_order
  )
  expect_equal(bad$resolution_table$status, "not_found")
})


test_that("split_modules maps resolution values per requested module", {
  resolution_fun <- get(".hc_split_resolution_by_module", asNamespace("hcocena"))
  target_colors <- c("blue", "red", "green", "gold")
  resolved_labels <- c("M2", "M3", "M4", "M8")
  resolution_table <- data.frame(
    input = resolved_labels,
    resolved_color = target_colors,
    resolved_label = resolved_labels,
    status = "ok",
    stringsAsFactors = FALSE
  )

  expect_equal(
    resolution_fun(0.8, target_colors, resolved_labels, resolution_table),
    stats::setNames(rep(0.8, 4), target_colors)
  )
  expect_equal(
    resolution_fun(c(0.8, 0.8, 0.9, 0.3), target_colors, resolved_labels, resolution_table),
    stats::setNames(c(0.8, 0.8, 0.9, 0.3), target_colors)
  )
  expect_equal(
    resolution_fun(
      c(M2 = 0.8, M3 = 0.8, green = 0.9, M8 = 0.3),
      target_colors,
      resolved_labels,
      resolution_table
    ),
    stats::setNames(c(0.8, 0.8, 0.9, 0.3), target_colors)
  )

  expect_error(
    resolution_fun(c(0.8, 0.9), target_colors, resolved_labels, resolution_table),
    "one value per resolved module"
  )
  expect_error(
    resolution_fun(c(M2 = 0.8, M9 = 0.9), target_colors, resolved_labels, resolution_table),
    "Unmatched: M9"
  )
  expect_error(
    resolution_fun(c(M2 = 0.8, M3 = 0.8), target_colors, resolved_labels, resolution_table),
    "Missing"
  )
})


test_that("split_modules normalizes label maps and preserves duplicate GFC columns in children", {
  normalize_fun <- get(".hc_normalize_module_label_map_for_split", asNamespace("hcocena"))
  child_fun <- get(".hc_build_child_cluster_rows", asNamespace("hcocena"))

  expect_equal(
    normalize_fun(c(M1 = "red", M2 = "blue"), available_colors = c("red", "blue")),
    c(red = "M1", blue = "M2")
  )

  template_row <- data.frame(
    clusters = "M1",
    gene_no = 4L,
    gene_n = "g1,g2,g3,g4",
    cluster_included = "yes",
    color = "red",
    conditions = "",
    grp_means = "",
    vertexsize = 3,
    stringsAsFactors = FALSE
  )
  membership <- c(g1 = 1L, g2 = 1L, g3 = 2L, g4 = 2L)
  child_rows <- child_fun(
    template_row = template_row,
    membership = membership,
    member_levels = c("1", "2"),
    child_colors = c("#111111", "#222222"),
    child_labels = c("#111111" = "M1.1", "#222222" = "M1.2"),
    gfc_all = data.frame(
      T1 = c(1, 3, 10, 12),
      T1 = c(5, 7, 20, 22),
      Gene = c("g1", "g2", "g3", "g4"),
      check.names = FALSE
    )
  )

  expect_identical(child_rows$conditions, c("T1#T1", "T1#T1"))
  expect_identical(child_rows$grp_means, c("2,6", "11,21"))
})


test_that("split_modules resolution_test_only without grid does not split", {
  hc <- methods::new("HCoCenaExperiment")

  expect_error(
    hcocena::hc_split_modules(hc, modules = "M1", resolution_test_only = TRUE),
    "requires `resolution_grid`"
  )
})


test_that("split-module labels trigger preserve-existing heatmap numbering", {
  has_split_labels <- get(".hc_module_label_map_has_split_labels", asNamespace("hcocena"))

  expect_true(has_split_labels(c(turquoise = "M10.1", red = "M2")))
  expect_false(has_split_labels(c(turquoise = "M10", red = "M2")))
  # Repeated splits append further numeric suffixes and must still be detected.
  expect_true(has_split_labels(c(turquoise = "M10.1.2", red = "M2")))
  expect_false(has_split_labels(NULL))
  expect_false(has_split_labels(character()))
})


test_that("cluster heatmap view resolves module order and full GFC scale", {
  resolve_modules <- get(".hc_heatmap_view_resolve_modules", asNamespace("hcocena"))
  resolve_scale <- get(".hc_heatmap_view_resolve_gfc_scale_limits", asNamespace("hcocena"))
  resolve_row_order <- get(".hc_resolve_cluster_heatmap_row_order", asNamespace("hcocena"))

  cluster_calc <- list(
    cluster_information = data.frame(
      color = c("red", "blue", "green"),
      cluster_included = c("yes", "yes", "yes"),
      stringsAsFactors = FALSE
    ),
    heatmap_matrix = matrix(
      c(-5, 1, 3, 2, -1, 4),
      nrow = 3,
      dimnames = list(c("red", "blue", "green"), c("T1", "T2"))
    ),
    heatmap_row_order = c("blue", "red", "green"),
    module_label_map = c(red = "M1", blue = "M2", green = "M3")
  )

  resolved <- resolve_modules(
    modules = c("M3", 1, "red"),
    cluster_calc = cluster_calc
  )
  expect_equal(resolved$target_colors, c("green", "blue", "red"))
  expect_equal(resolved$resolved_labels, c("M3", "M2", "M1"))
  expect_equal(
    resolve_row_order(
      row_order = c("M3", "1", "red"),
      cluster_calc = cluster_calc,
      available_colors = c("red", "blue", "green")
    ),
    c("green", "blue", "red")
  )
  expect_error(
    resolve_row_order(
      row_order = "M9",
      cluster_calc = cluster_calc,
      available_colors = c("red", "blue", "green")
    ),
    "Unknown entries in `row_order`"
  )

  hc <- methods::new("HCoCenaExperiment")
  hc@integration@cluster <- S4Vectors::SimpleList(cluster_calc)
  expect_equal(resolve_scale(hc, cluster_calc), c(-5, 5))
  cluster_calc[["gfc_scale_limits"]] <- c(-3, 3)
  expect_equal(resolve_scale(hc, cluster_calc), c(-3, 3))
  expect_true("modules" %in% names(formals(hcocena::hc_plot_cluster_heatmap_view)))
  expect_true("write_module_tables" %in% names(formals(get("plot_cluster_heatmap_new", asNamespace("hcocena")))))
})


test_that("cluster colour palette is unique (merge_clusters keys modules by colour)", {
  palette <- hcocena:::get_cluster_colours()
  expect_false(any(duplicated(palette)))
  expect_true(all(nzchar(palette)))
})


test_that(".hc_with_seed handles NULL / empty / NA seeds without erroring", {
  with_seed <- get(".hc_with_seed", asNamespace("hcocena"))

  expect_identical(with_seed(NULL, 1 + 1), 2)
  expect_identical(with_seed(integer(0), 1 + 1), 2)
  expect_identical(with_seed(NA, 1 + 1), 2)

  # A real seed still makes the draw reproducible.
  a <- with_seed(42, stats::runif(1))
  b <- with_seed(42, stats::runif(1))
  expect_identical(a, b)
})


test_that(".hc_longitudinal_group_singletons leaves all-singleton input unchanged", {
  group_singletons <- get(".hc_longitudinal_group_singletons", asNamespace("hcocena"))

  ids <- c(a = "1", b = "2", c = "3")
  expect_warning(
    out <- group_singletons(ids, snn = NULL, group.singletons = TRUE, verbose = FALSE),
    NA
  )
  expect_identical(out, ids)
})


test_that("regression: duplicate heatmap condition labels get layer prefixes and keep axis multiplicity", {
  display_fun <- get(".hc_gfc_display_col_labels", asNamespace("hcocena"))
  count_fun <- get(".hc_gfc_display_count_labels", asNamespace("hcocena"))
  width_scale_fun <- get(".hc_gfc_duplicate_condition_width_scale", asNamespace("hcocena"))
  norm_order_fun <- get(".hc_normalize_heatmap_axis_order", asNamespace("hcocena"))

  hcobject <- hcocena:::.hc_default_object()
  hcobject$layers <- stats::setNames(list(character(0), character(0)), c("set1", "set2"))
  hcobject$layers_names <- c("RNA", "PROT")
  hcobject$global_settings$voi <- "group"
  hcobject$data$set1_anno <- data.frame(
    group = c("T1", "T1", "T2", "T3"),
    stringsAsFactors = FALSE,
    row.names = paste0("s", 1:4)
  )
  hcobject$data$set2_anno <- data.frame(
    group = c("T1", "T2", "T2", "T3"),
    stringsAsFactors = FALSE,
    row.names = paste0("p", 1:4)
  )
  hcobject$layer_specific_outputs$set1 <- list(
    part2 = list(
      GFC_all_genes = data.frame(
        T1 = 1,
        T2 = 2,
        T3 = 3,
        Gene = "g1",
        check.names = FALSE
      )
    )
  )
  hcobject$layer_specific_outputs$set2 <- list(
    part2 = list(
      GFC_all_genes = data.frame(
        T1 = 4,
        T2 = 5,
        T3 = 6,
        Gene = "g1",
        check.names = FALSE
      )
    )
  )

  raw_cols <- c("T1", "T2", "T3", "T1", "T2", "T3")
  expect_identical(
    display_fun(hcobject, raw_cols),
    c("RNA: T1", "RNA: T2", "RNA: T3", "PROT: T1", "PROT: T2", "PROT: T3")
  )
  expect_identical(
    count_fun(hcobject, raw_cols),
    c("RNA: T1  [2]", "RNA: T2  [1]", "RNA: T3  [1]", "PROT: T1  [1]", "PROT: T2  [2]", "PROT: T3  [1]")
  )
  expect_equal(width_scale_fun(hcobject, raw_cols), 1.12)
  expect_equal(width_scale_fun(hcobject, c("T1", "T2", "T3")), 1)
  expect_identical(
    norm_order_fun(raw_cols, raw_cols),
    raw_cols
  )
})


test_that("regression: enrichment-related heatmap legends use standard font settings", {
  fun_enrich_src <- paste(deparse(get(".hc_functional_enrichment_driver", asNamespace("hcocena"))), collapse = "\n")
  upstream_src <- paste(deparse(get(".hc_upstream_inference_driver", asNamespace("hcocena"))), collapse = "\n")
  replot_src <- paste(deparse(get("replot_cluster_heatmap", asNamespace("hcocena"))), collapse = "\n")

  expect_false(grepl('fontfamily = "mono"', fun_enrich_src, fixed = TRUE))
  expect_false(grepl('fontfamily = "mono"', upstream_src, fixed = TRUE))
  expect_false(grepl('fontfamily = "mono"', replot_src, fixed = TRUE))
})


test_that("regression: enrichment plot body borders use the same thin line style as the main heatmap", {
  fun_enrich_path <- test_path("..", "..", "R", "hc_functional_enrichment.R")
  if (!file.exists(fun_enrich_path)) {
    skip("Source files are not available in the installed-package test context.")
  }

  fun_enrich_src <- paste(
    readLines(fun_enrich_path, warn = FALSE),
    collapse = "\n"
  )

  expect_true(grepl('shared_heatmap_line_lwd <- 0\\.5', fun_enrich_src))
  expect_true(grepl('rect_gp = grid::gpar\\(col = "black",\\s*lwd = shared_heatmap_line_lwd\\)', fun_enrich_src))
  expect_true(grepl('panel_border_gp <- grid::gpar\\(col = "black",\\s*fill = NA,\\s*lwd = shared_heatmap_line_lwd\\)', fun_enrich_src))
  expect_true(grepl('panel_border_gp_all <- grid::gpar\\(col = "black",\\s*fill = NA,\\s*lwd = shared_heatmap_line_lwd\\)', fun_enrich_src))
  expect_true(grepl('panel_border_slices_all <- base::seq_len\\(base::nlevels\\(term_db_levels\\)\\)', fun_enrich_src))
  expect_true(grepl('for \\(slice_idx in slices_use\\)', fun_enrich_src))
  expect_true(grepl('decorate_heatmap_body\\(nm, slice = slice_idx, \\{', fun_enrich_src))
  expect_true(grepl('panel_border_slices = panel_border_slices_all', fun_enrich_src, fixed = TRUE))
  expect_false(grepl('border_targets <- base::unique(base::c("GFC", "enrichment"', fun_enrich_src, fixed = TRUE))
  expect_false(grepl('border_targets <- base::unique(base::c("GFC", "enrichment_all"', fun_enrich_src, fixed = TRUE))
})


test_that("regression: enrichment module labels and all-db headers avoid visual overlap", {
  fun_enrich_path <- test_path("..", "..", "R", "hc_functional_enrichment.R")
  if (!file.exists(fun_enrich_path)) {
    skip("Source files are not available in the installed-package test context.")
  }

  fun_enrich_src <- paste(
    readLines(fun_enrich_path, warn = FALSE),
    collapse = "\n"
  )

  expect_true(grepl("module_label_fit <- .hc_module_label_fit_pt", fun_enrich_src, fixed = TRUE))
  expect_true(grepl("module_box_anno <- .hc_module_label_box_annotation", fun_enrich_src, fixed = TRUE))
  expect_true(grepl("label_fontsize_pt = module_label_fit$base_pt_size", fun_enrich_src, fixed = TRUE))
  expect_true(grepl("panel_title_header_offset_mm_all <-", fun_enrich_src, fixed = TRUE))
  expect_true(grepl("panel_title_offset_mm = panel_title_header_offset_mm_all", fun_enrich_src, fixed = TRUE))
})


test_that("regression: auxiliary heatmap module labels use the shared box-fit guard", {
  llm_path <- test_path("..", "..", "R", "hc_plot_llm_enrichment.R")
  upstream_path <- test_path("..", "..", "R", "hc_upstream_inference.R")
  api_path <- test_path("..", "..", "R", "hc_api.R")
  if (!file.exists(llm_path) || !file.exists(upstream_path) || !file.exists(api_path)) {
    skip("Source files are not available in the installed-package test context.")
  }

  llm_src <- paste(readLines(llm_path, warn = FALSE), collapse = "\n")
  upstream_src <- paste(readLines(upstream_path, warn = FALSE), collapse = "\n")
  api_src <- paste(readLines(api_path, warn = FALSE), collapse = "\n")

  expect_true(grepl("module_label_fit <- .hc_module_label_fit_pt", llm_src, fixed = TRUE))
  expect_true(grepl("module_box_anno <- .hc_module_label_box_annotation", llm_src, fixed = TRUE))
  expect_true(grepl("label_fontsize_pt = module_label_fit$base_pt_size", llm_src, fixed = TRUE))
  expect_true(grepl("module_label_fit <- .hc_module_label_fit_pt", upstream_src, fixed = TRUE))
  expect_true(grepl("module_box_anno <- .hc_module_label_box_annotation", upstream_src, fixed = TRUE))
  expect_true(grepl("label_fontsize_pt = module_label_fit$base_pt_size", upstream_src, fixed = TRUE))
  expect_true(grepl("is_combined_panel <- panel_key_chr %in% c(\"top_all_dbs\", \"top_all_dbs_mixed\")", api_src, fixed = TRUE))
  expect_true(grepl("max(8, 0.50 * as.numeric(fontsize))", api_src, fixed = TRUE))
})


test_that("regression: lightweight heatmap cache works without ComplexHeatmap object", {
  cache_info <- get(".hc_heatmap_cache_info", asNamespace("hcocena"))
  select_col_order <- get(".hc_select_heatmap_col_order", asNamespace("hcocena"))
  llm_heatmap_info <- get(".hc_llm_heatmap_info", asNamespace("hcocena"))
  llm_capture <- get(".hc_llm_capture_combined_heatmap_grob", asNamespace("hcocena"))
  llm_style <- get(".hc_llm_resolve_heatmap_style", asNamespace("hcocena"))
  llm_title_wrap_width <- get(".hc_llm_title_wrap_width", asNamespace("hcocena"))
  plot_heatmap <- get(".hc_plot_cluster_heatmap_driver", asNamespace("hcocena"))
  plot_heatmap_new <- get("plot_cluster_heatmap_new", asNamespace("hcocena"))
  plot_network <- get(".hc_plot_integrated_network_driver", asNamespace("hcocena"))
  fun_enrich <- get(".hc_functional_enrichment_driver", asNamespace("hcocena"))
  up_inf <- get(".hc_upstream_inference_driver", asNamespace("hcocena"))
  knowledge_plot <- get(".hc_plot_enrichment_upstream_network_driver", asNamespace("hcocena"))
  llm_plot <- get("hc_plot_llm_enrichment", asNamespace("hcocena"))

  cluster_calc <- list(
    heatmap_matrix = matrix(
      c(1, 2, 3, 4),
      nrow = 2,
      dimnames = list(c("M1", "M2"), c("T1", "T2"))
    ),
    heatmap_row_order = c("M2", "M1"),
    heatmap_column_order = c("T2", "T1")
  )

  info <- cache_info(cluster_calc)

  expect_equal(info$row_order, c("M2", "M1"))
  expect_equal(info$col_order, c("T2", "T1"))
  expect_equal(base::rownames(info$matrix), c("M1", "M2"))
  expect_equal(
    select_col_order(
      available_cols = c("A", "B", "C"),
      plot_order = c("C"),
      main_order = c("B", "A"),
      fallback_order = c("A", "C")
    ),
    c("C", "A", "B")
  )
  expect_equal(
    select_col_order(
      available_cols = c("A", "B", "C"),
      main_order = c("B", "A"),
      fallback_order = c("C")
    ),
    c("B", "A", "C")
  )
  expect_equal(
    select_col_order(
      available_cols = c("A", "B", "C"),
      main_order = c("Z"),
      fallback_order = c("C", "A")
    ),
    c("C", "A", "B")
  )
  expect_identical(formals(plot_heatmap)$return_HM, FALSE)
  expect_identical(formals(plot_heatmap_new)$return_HM, FALSE)
  expect_identical(formals(plot_network)$store_plot, FALSE)
  expect_true("col_order" %in% names(formals(fun_enrich)))
  expect_true("consistent_terms" %in% names(formals(fun_enrich)))
  expect_true("col_order" %in% names(formals(up_inf)))
  expect_true("col_order" %in% names(formals(knowledge_plot)))
  expect_true("col_order" %in% names(formals(llm_plot)))
  expect_true("cluster_columns" %in% names(formals(fun_enrich)))
  expect_true("cluster_columns" %in% names(formals(up_inf)))
  expect_true("cluster_columns" %in% names(formals(knowledge_plot)))
  expect_true("cluster_columns" %in% names(formals(llm_plot)))
  expect_true("heatmap_col_order" %in% names(formals(fun_enrich)))
  expect_true("heatmap_col_order" %in% names(formals(up_inf)))
  expect_true("heatmap_col_order" %in% names(formals(knowledge_plot)))
  expect_true("heatmap_col_order" %in% names(formals(llm_plot)))
  expect_true("heatmap_cluster_columns" %in% names(formals(fun_enrich)))
  expect_true("heatmap_cluster_columns" %in% names(formals(up_inf)))
  expect_true("heatmap_cluster_columns" %in% names(formals(knowledge_plot)))
  expect_true("heatmap_cluster_columns" %in% names(formals(llm_plot)))
  expect_true("module_label_fontsize" %in% names(formals(llm_plot)))
  expect_true("module_label_pt_size" %in% names(formals(llm_plot)))
  expect_true("module_box_width_cm" %in% names(formals(llm_plot)))

  hc <- methods::new("HCoCenaExperiment")
  hc@integration@cluster <- S4Vectors::SimpleList(
    heatmap_matrix = matrix(
      c(-1, 0.5, 1, -0.25),
      nrow = 2,
      dimnames = list(c("red", "blue"), c("T1", "T2"))
    ),
    heatmap_cluster_raw = ComplexHeatmap::add_heatmap(
      ComplexHeatmap::Heatmap(
        matrix(
          c(-1, 0.5, 1, -0.25),
          nrow = 2,
          dimnames = list(c("red", "blue"), c("T1", "T2"))
        ),
        name = "GFC",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_row_names = FALSE,
        width = grid::unit(90, "mm"),
        height = grid::unit(24, "mm")
      ),
      ComplexHeatmap::columnAnnotation(
        groups = ComplexHeatmap::anno_text(c("T1", "T2"))
      ),
      direction = "vertical"
    ),
    heatmap_row_order = c("red", "blue"),
    heatmap_column_order = c("T2", "T1"),
    module_label_map = c(red = "M1", blue = "M2"),
    module_label_fontsize = 9,
    module_label_pt_size = 0.3,
    module_box_width_cm = 0.9,
    heatmap_cell_size_mm = 4.4,
    gfc_colors = c("#112233", "#f7f7f7", "#cc3311"),
    gfc_scale_limits = c(-3, 3),
    overall_plot_scale = 1.25
  )
  info2 <- llm_heatmap_info(hc)
  expect_true(isTRUE(info2$draw_supported))
  expect_equal(info2$col_order, c("T2", "T1"))
  expect_equal(info2$module_order, c("M1", "M2"))
  expect_equal(info2$stored_heatmap_cell_size_mm, 4.4)
  expect_equal(info2$stored_gfc_colors, c("#112233", "#f7f7f7", "#cc3311"))
  expect_equal(info2$stored_gfc_scale_limits, c(-3, 3))
  expect_equal(info2$stored_overall_plot_scale, 1.25)

  style <- llm_style(
    heatmap_info = info2,
    module_labels = c("M1", "M2"),
    n_heat_rows = 2,
    n_heat_cols = 2,
    mat_use = matrix(c(-1, 0.5, 1, -0.25), nrow = 2)
  )
  expect_equal(style$module_label_fontsize, 4.2)
  expect_equal(style$module_label_pt_size, 0.90)
  expect_equal(style$module_box_width_cm, 0.56)
  expect_equal(style$cell_size_mm, 4.4)
  expect_equal(style$gfc_colors, c("#112233", "#f7f7f7", "#cc3311"))
  expect_equal(style$gfc_scale_limits, c(-3, 3))
  expect_equal(style$overall_plot_scale, 1.25)
  expect_equal(llm_title_wrap_width(NA_real_), 42L)
  expect_gte(llm_title_wrap_width(6), 32L)
  expect_lte(llm_title_wrap_width(6), 72L)

  style_override <- llm_style(
    heatmap_info = info2,
    module_labels = c("M1", "M2"),
    n_heat_rows = 2,
    n_heat_cols = 2,
    mat_use = matrix(c(-1, 0.5, 1, -0.25), nrow = 2),
    module_label_fontsize = 4.2,
    module_label_pt_size = 0.22,
    module_box_width_cm = 0.58
  )
  expect_equal(style_override$module_label_fontsize, 4.2)
  expect_equal(style_override$module_label_pt_size, 0.22)
  expect_equal(style_override$module_box_width_cm, 0.58)

  summary_tbl <- data.frame(
    module = c("M1", "M2"),
    module_color = c("red", "blue"),
    term_plot = c("alpha process", "beta process"),
    text_color = c("#111111", "#222222"),
    label_color = c("white", "white"),
    stringsAsFactors = FALSE
  )
  grob <- llm_capture(
    heatmap_info = info2,
    summary_tbl = summary_tbl,
    max_chars = 90,
    text_size = 4
  )
  expect_s3_class(grob, "grob")

  hc@satellite <- S4Vectors::SimpleList(list(
    llm_enrichment = list(module_1 = list(status = "ok")),
    llm_enrichment_summary = data.frame(
      module = c("M1", "M2"),
      module_color = c("red", "blue"),
      general_processes = c("alpha process", "beta process"),
      contextual_state = c("state a", "state b"),
      key_regulators = c("reg a", "reg b"),
      stringsAsFactors = FALSE
    )
  ))
  p_heat <- llm_plot(
    hc,
    fields = "general_processes",
    save = FALSE
  )
  expect_s3_class(p_heat, "hc_llm_heatmap_plot")

  out_dir <- file.path(tempdir(), paste0("hc_llm_heatmap_export_", Sys.getpid()))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  hc@config@paths <- S4Vectors::DataFrame(dir_output = out_dir)
  hc@config@global <- S4Vectors::DataFrame(save_folder = "exports")
  p_heat_export <- llm_plot(
    hc,
    fields = "general_processes",
    save = TRUE,
    file_stem = "llm_heatmap_test"
  )
  heatmap_exports <- attr(p_heat_export, "output_files", exact = TRUE)
  expect_true(file.exists(heatmap_exports$pdf))
  expect_true(file.exists(heatmap_exports$png))
  if (requireNamespace("png", quietly = TRUE)) {
    read_png <- getExportedValue("png", "readPNG")
    png_img <- read_png(heatmap_exports$png)
    rgb_img <- png_img[, , seq_len(min(3L, dim(png_img)[[3L]])), drop = FALSE]
    expect_gt(stats::sd(as.numeric(rgb_img)), 0.02)
    expect_gt(mean(rgb_img), 0.1)
    expect_lt(mean(rgb_img), 0.99)
  }
})


test_that("regression: llm plot export writes pdf and png into configured output dir", {
  hc <- methods::new("HCoCenaExperiment")
  out_dir <- file.path(tempdir(), paste0("hc_llm_plot_export_", Sys.getpid()))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  hc@config@paths <- S4Vectors::DataFrame(dir_output = out_dir)
  hc@config@global <- S4Vectors::DataFrame(save_folder = "exports")
  hc@satellite <- S4Vectors::SimpleList(list(
    llm_enrichment = list(module_1 = list(status = "ok")),
    llm_enrichment_summary = data.frame(
      module = "M1",
      module_color = "steelblue",
      general_processes = "Interferon signaling",
      contextual_state = "Acute antiviral activation",
      key_regulators = "STAT1 / IRF7",
      stringsAsFactors = FALSE
    )
  ))

  p <- hcocena::hc_plot_llm_enrichment(
    hc,
    with_heatmap = FALSE,
    fields = "general_processes",
    save = TRUE,
    file_stem = "llm_test"
  )

  export_files <- attr(p, "output_files", exact = TRUE)
  expect_true(file.exists(export_files$pdf))
  expect_true(file.exists(export_files$png))
  expect_match(export_files$pdf, "exports")
  expect_match(export_files$pdf, "llm_test_M1_general_processes\\.pdf$")

  if (requireNamespace("png", quietly = TRUE)) {
    read_png <- getExportedValue("png", "readPNG")
    png_img <- read_png(export_files$png)
    rgb_img <- png_img[, , seq_len(min(3L, dim(png_img)[[3L]])), drop = FALSE]
    border_rgb <- c(
      as.numeric(rgb_img[1L, , , drop = FALSE]),
      as.numeric(rgb_img[dim(rgb_img)[[1L]], , , drop = FALSE]),
      as.numeric(rgb_img[, 1L, , drop = FALSE]),
      as.numeric(rgb_img[, dim(rgb_img)[[2L]], , drop = FALSE])
    )
    expect_gt(mean(border_rgb), 0.85)
    expect_gt(mean(rgb_img), 0.4)
  }
})


test_that("regression: llm heatmap print reserves a separate title row", {
  dummy <- grid::rectGrob(gp = grid::gpar(fill = "grey80", col = NA))
  class(dummy) <- unique(c("hc_llm_heatmap_plot", class(dummy)))
  attr(dummy, "llm_title") <- "AI-assisted module interpretation: General processes"
  attr(dummy, "llm_text_size") <- 4

  captured <- grid::grid.grabExpr(print(dummy))
  title_child <- captured$children[[1]]
  plot_child <- captured$children[[2]]

  expect_s3_class(title_child, "text")
  expect_match(title_child$label, "AI-assisted module interpretation")
  expect_false(is.null(title_child$vp))
  expect_false(is.null(plot_child$vp))
  expect_false(identical(title_child$vp, plot_child$vp))
})


test_that("regression: llm heatmap plot body is top-aligned under the title", {
  hc <- methods::new("HCoCenaExperiment")
  hc@integration@cluster <- S4Vectors::SimpleList(
    heatmap_matrix = matrix(
      c(-1, 0.5, 1, -0.25),
      nrow = 2,
      dimnames = list(c("red", "blue"), c("T1", "T2"))
    ),
    heatmap_cluster_raw = ComplexHeatmap::add_heatmap(
      ComplexHeatmap::Heatmap(
        matrix(
          c(-1, 0.5, 1, -0.25),
          nrow = 2,
          dimnames = list(c("red", "blue"), c("T1", "T2"))
        ),
        name = "GFC",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_row_names = FALSE,
        width = grid::unit(90, "mm"),
        height = grid::unit(24, "mm")
      ),
      ComplexHeatmap::columnAnnotation(
        groups = ComplexHeatmap::anno_text(c("T1", "T2"))
      ),
      direction = "vertical"
    ),
    heatmap_row_order = c("red", "blue"),
    heatmap_column_order = c("T2", "T1"),
    module_label_map = c(red = "M1", blue = "M2"),
    module_label_fontsize = 9,
    module_label_pt_size = 0.3,
    module_box_width_cm = 0.9,
    heatmap_cell_size_mm = 4.4,
    gfc_colors = c("#112233", "#f7f7f7", "#cc3311"),
    gfc_scale_limits = c(-3, 3),
    overall_plot_scale = 1.25
  )
  hc@satellite <- S4Vectors::SimpleList(list(
    llm_enrichment = list(module_1 = list(status = "ok")),
    llm_enrichment_summary = data.frame(
      module = c("M1", "M2"),
      module_color = c("red", "blue"),
      general_processes = c("alpha process", "beta process"),
      contextual_state = c("state a", "state b"),
      key_regulators = c("reg a", "reg b"),
      stringsAsFactors = FALSE
    )
  ))

  p <- hcocena::hc_plot_llm_enrichment(
    hc,
    fields = "general_processes",
    save = FALSE
  )

  captured <- grid::grid.grabExpr(print(p))
  plot_child <- captured$children[[2]]
  plot_parent_vp <- plot_child$childrenvp[[1]]$parent

  expect_equal(grid::convertY(plot_parent_vp$y, "npc", valueOnly = TRUE), 1)
  expect_equal(plot_parent_vp$justification, c(0.5, 1))
})


test_that("regression: llm heatmap plot ignores cached dendrograms with mismatched row count", {
  llm_capture <- get(".hc_llm_capture_combined_heatmap_grob", asNamespace("hcocena"))

  source_mat <- matrix(
    c(-1, 0.5,
      1, -0.25,
      0.2, 0.8),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("red", "blue", "green"), c("T1", "T2"))
  )
  heatmap_info <- list(
    matrix = source_mat[c("red", "blue"), , drop = FALSE],
    raw_heatmap_obj = ComplexHeatmap::Heatmap(
      source_mat,
      name = "GFC",
      cluster_rows = TRUE,
      cluster_columns = FALSE,
      show_row_names = FALSE
    ),
    heatmap_obj = NULL,
    draw_supported = TRUE,
    col_order = c("T1", "T2"),
    row_ids = c("red", "blue"),
    module_order = c("M1", "M2"),
    module_by_row = c("M1", "M2"),
    stored_module_label_fontsize = 5,
    stored_module_label_pt_size = 0.3,
    stored_module_box_width_cm = 0.9,
    stored_heatmap_cell_size_mm = 4.4,
    stored_gfc_colors = c("#112233", "#f7f7f7", "#cc3311"),
    stored_gfc_scale_limits = c(-3, 3),
    stored_overall_plot_scale = 1
  )
  summary_tbl <- data.frame(
    module = c("M1", "M2"),
    module_color = c("red", "blue"),
    term_plot = c("alpha process", "beta process"),
    text_color = c("#111111", "#222222"),
    label_color = c("white", "white"),
    stringsAsFactors = FALSE
  )

  expect_s3_class(
    llm_capture(
      heatmap_info = heatmap_info,
      summary_tbl = summary_tbl,
      max_chars = 90,
      text_size = 4
    ),
    "grob"
  )
})


test_that("regression: split labels and significance suffixes expand module boxes only when needed", {
  draw_width_fun <- get(".hc_module_label_draw_width_cm", asNamespace("hcocena"))

  expect_equal(
    draw_width_fun(
      module_box_width_cm = 0.62,
      module_labels_display = c("M1", "M2"),
      user_set_module_box_width_cm = FALSE,
      module_sig_integrated = FALSE,
      max_sig_stars = 0
    ),
    0.62
  )
  expect_equal(
    draw_width_fun(
      module_box_width_cm = 0.62,
      module_labels_display = c("M1.2", "M2"),
      user_set_module_box_width_cm = FALSE,
      module_sig_integrated = FALSE,
      max_sig_stars = 0
    ),
    0.94
  )
  expect_equal(
    draw_width_fun(
      module_box_width_cm = 0.62,
      module_labels_display = c("M1**", "M2"),
      user_set_module_box_width_cm = FALSE,
      module_sig_integrated = TRUE,
      max_sig_stars = 2
    ),
    0.94
  )
  expect_equal(
    draw_width_fun(
      module_box_width_cm = 0.62,
      module_labels_display = c("M1.2**", "M2"),
      user_set_module_box_width_cm = FALSE,
      module_sig_integrated = TRUE,
      max_sig_stars = 2
    ),
    1.32
  )
  expect_equal(
    draw_width_fun(
      module_box_width_cm = 0.62,
      module_labels_display = c("M1.2", "M2"),
      user_set_module_box_width_cm = TRUE,
      module_sig_integrated = FALSE,
      max_sig_stars = 0
    ),
    0.62
  )
})


test_that("regression: module label auto-fit uses one safe size for all labels", {
  fit_fun <- get(".hc_module_label_fit_pt", asNamespace("hcocena"))
  effective_fontsize_fun <- get(".hc_module_label_effective_fontsize", asNamespace("hcocena"))
  width_fun <- get(".hc_module_label_text_width_cm", asNamespace("hcocena"))

  fit <- fit_fun(
    module_label_pt_size = 0.55,
    module_box_width_cm = 0.9,
    module_labels_display = c("M1", "M4.12***"),
    n_heat_rows = 13,
    cell_size_mm = 5.4
  )
  fitted_widths <- width_fun(c("M1", "M4.12***"), fontsize_pt = fit$pt_size)

  expect_true(fit$shrunk)
  expect_true(fit$width_limited || fit$height_limited)
  expect_equal(fit$pt_size[[1]], fit$pt_size[[2]])
  expect_lt(fit$pt_size[[1]], fit$base_pt_size[[1]])
  expect_true(all(fitted_widths <= fit$available_width_cm + 1e-6))
  expect_lte(max(fit$pt_size), fit$available_height_pt + 1e-6)

  short_fit <- fit_fun(
    module_label_pt_size = 0.22,
    module_box_width_cm = 0.9,
    module_labels_display = c("M1", "M2"),
    n_heat_rows = 20,
    cell_size_mm = 5.4
  )
  expect_false(short_fit$shrunk)
  expect_equal(short_fit$pt_size[[1]], short_fit$pt_size[[2]])

  fontsize_fit <- fit_fun(
    module_label_pt_size = 0.22,
    module_box_width_cm = 1.2,
    module_labels_display = "M1",
    n_heat_rows = 10,
    cell_size_mm = 5.4,
    module_label_fontsize = 11,
    use_fontsize_request = TRUE
  )
  expect_equal(fontsize_fit$base_pt_size[[1]], 11)

  fill_fit <- fit_fun(
    module_label_pt_size = 0.95,
    module_box_width_cm = 0.80,
    module_labels_display = c("M1", "M2", "M3"),
    n_heat_rows = 5,
    cell_size_mm = 7
  )
  fill_widths <- width_fun(c("M1", "M2", "M3"), fontsize_pt = fill_fit$pt_size)
  expect_equal(fill_fit$pt_size[[1]], fill_fit$pt_size[[2]])
  expect_gt(fill_fit$pt_size[[1]], 10)
  expect_true(all(fill_widths <= fill_fit$available_width_cm + 1e-6))
  expect_equal(effective_fontsize_fun(fill_fit, fallback_fontsize = 5), fill_fit$pt_size[[1]])

  tiny_box_fit <- fit_fun(
    module_label_pt_size = 0.9,
    module_box_width_cm = 0.08,
    module_labels_display = c("M1", "M123456789"),
    n_heat_rows = 6,
    cell_size_mm = 5.4,
    module_label_fontsize = 80,
    use_fontsize_request = TRUE
  )
  tiny_box_widths <- width_fun(c("M1", "M123456789"), fontsize_pt = tiny_box_fit$pt_size)
  expect_equal(tiny_box_fit$pt_size[[1]], tiny_box_fit$pt_size[[2]])
  expect_true(tiny_box_fit$below_preferred_min)
  expect_true(all(tiny_box_widths <= tiny_box_fit$available_width_cm + 1e-6))
})


test_that("regression: current-device heatmap view scales down when device is smaller than export size", {
  screen_fit_fun <- get(".hc_heatmap_screen_fit_scale", asNamespace("hcocena"))

  expect_equal(
    screen_fit_fun(
      total_width_mm = 120,
      total_height_mm = 100,
      device_size_mm = c(240, 180)
    ),
    1
  )
  expect_lt(
    screen_fit_fun(
      total_width_mm = 260,
      total_height_mm = 220,
      device_size_mm = c(180, 140)
    ),
    1
  )
  expect_gt(
    screen_fit_fun(
      total_width_mm = 260,
      total_height_mm = 220,
      device_size_mm = c(180, 140)
    ),
    0
  )
})


test_that("regression: split suffix detection differs between unsplit and split module labels", {
  labels_unsplit <- c("M1", "M2", "M10")
  labels_split <- c("M1.2", "M1.3", "M2")

  expect_false(any(grepl("\\.[0-9]+", labels_unsplit)))
  expect_true(any(grepl("\\.[0-9]+", labels_split)))
})


test_that("regression: .hc_run_driver resolves functions from the legacy env", {
  run_legacy <- get(".hc_run_driver", asNamespace("hcocena"))
  hc <- methods::new("HCoCenaExperiment")
  legacy_env <- new.env(parent = baseenv())
  legacy_env$target_env <- legacy_env

  fun <- function(gene_sets = "Hallmark", heatmap_col_order = NULL) {
    base::assign("hcobject", list(used = heatmap_col_order), envir = target_env)
    invisible(NULL)
  }
  environment(fun) <- legacy_env
  legacy_env$functional_enrichment <- fun
  legacy_env$hcobject <- list(initial = TRUE)

  global_fun <- function(gene_sets = "Hallmark") {
    stop("wrong global function selected")
  }
  base::assign("functional_enrichment", global_fun, envir = .GlobalEnv)
  on.exit(base::rm("functional_enrichment", envir = .GlobalEnv), add = TRUE)

  testthat::local_mocked_bindings(
    as_hcobject = function(hc) list(initial = TRUE),
    as_hcocena = function(x) x,
    .hc_bind_bridge_hcobject = function(hcobject, envo = legacy_env) {
      base::assign("hcobject", hcobject, envir = envo)
      list(envo = envo, had_existing = FALSE, old_hcobject = NULL, binding_locked = FALSE)
    },
    .hc_restore_bridge_hcobject = function(state) invisible(NULL),
    .package = "hcocena"
  )

  out <- run_legacy(
    hc,
    fun = "functional_enrichment",
    gene_sets = "Hallmark",
    heatmap_col_order = c("T1", "T2"),
    envo = legacy_env
  )

  expect_equal(out$used, c("T1", "T2"))
})


test_that("regression: vllm think blocks are stripped before JSON parsing", {
  strip_think <- get(".hc_llm_strip_think_blocks", asNamespace("hcocena"))
  strip_fences <- get(".hc_llm_strip_json_fences", asNamespace("hcocena"))

  raw_text <- paste(
    "<think>",
    "First inspect the genes and reason privately.",
    "</think>",
    "```json",
    '{"general_processes":"interferon signaling","contextual_state":"interferon-high inflammatory state","key_regulators":"STAT1 / IRF7 / IRF9"}',
    "```",
    sep = "\n"
  )

  cleaned <- strip_fences(strip_think(raw_text))

  expect_false(grepl("<think>", cleaned, fixed = TRUE))
  expect_false(grepl("</think>", cleaned, fixed = TRUE))
  expect_true(grepl('"general_processes":"interferon signaling"', cleaned, fixed = TRUE))
})


test_that("regression: claude provider is accepted and resolves api key/model settings", {
  resolve_api_key <- get(".hc_llm_resolve_api_key", asNamespace("hcocena"))
  default_model <- get(".hc_llm_default_model", asNamespace("hcocena"))

  expect_true("claude" %in% eval(formals(hcocena::hc_llm_enrichment)$provider))
  expect_identical(default_model(NULL, "claude"), "claude-sonnet-4-6")

  withr::local_envvar(c(ANTHROPIC_API_KEY = "test-anthropic-key"))
  expect_identical(
    resolve_api_key(api_key = NULL, llm = "claude"),
    "test-anthropic-key"
  )
})

test_that("regression: only the hc_llm_enrichment interface is exported", {
  exports <- getNamespaceExports("hcocena")

  expect_true(all(c("hc_llm_enrichment", "hc_llm_import", "hc_plot_llm_enrichment") %in% exports))
  expect_false(any(c(
    "hc_module_function_llm", "hc_module_function_gemini",
    "hc_module_function_vllm", "hc_module_function_claude",
    "hc_plot_module_function_llm", "hc_plot_module_function_gemini"
  ) %in% exports))
})

test_that("regression: llm request helpers use ellmer backends", {
  gemini_src <- paste(deparse(get(".hc_llm_request_gemini", asNamespace("hcocena"))), collapse = "\n")
  claude_src <- paste(deparse(get(".hc_llm_request_claude", asNamespace("hcocena"))), collapse = "\n")
  openai_src <- paste(deparse(get(".hc_llm_request_openai", asNamespace("hcocena"))), collapse = "\n")
  vllm_src <- paste(deparse(get(".hc_llm_request_vllm", asNamespace("hcocena"))), collapse = "\n")

  expect_true(grepl("ellmer::chat_google_gemini", gemini_src, fixed = TRUE))
  expect_true(grepl("ellmer::chat_anthropic", claude_src, fixed = TRUE))
  expect_true(grepl("ellmer::chat_openai", openai_src, fixed = TRUE))
  expect_true(grepl("ellmer::chat_vllm", vllm_src, fixed = TRUE))
  expect_true(grepl("credentials = .hc_llm_api_key_credentials", gemini_src, fixed = TRUE))
  expect_true(grepl("credentials = .hc_llm_api_key_credentials", claude_src, fixed = TRUE))
  expect_true(grepl("credentials = .hc_llm_api_key_credentials", openai_src, fixed = TRUE))
  expect_true(grepl("credentials = .hc_llm_api_key_credentials", vllm_src, fixed = TRUE))
  expect_true(grepl("run_request\\(include_temperature = FALSE\\)", gemini_src))
  expect_true(grepl("max_tokens = 8000", vllm_src, fixed = TRUE))
  expect_true(grepl("chat_template_kwargs = list", vllm_src, fixed = TRUE))
  expect_true(grepl("enable_thinking = FALSE", vllm_src, fixed = TRUE))
})


test_that("regression: llm summary builder is robust for error-only results", {
  summary_fun <- get(".hc_llm_summary_from_results", asNamespace("hcocena"))
  err_res <- list(
    M1 = list(
      label = "M1",
      module = "M1",
      llm = "gemini",
      model = "gemini-3.1-pro-preview",
      gene_count_input = 10L,
      gene_count_sent = 0L,
      truncated = FALSE,
      status = "error",
      error_message = "HTTP 400",
      response = list(
        general_processes = NA_character_,
        contextual_state = NA_character_,
        key_regulators = "HTTP 400"
      ),
      timestamp = "2026-03-12 12:00:00"
    )
  )

  out <- summary_fun(err_res, hc = NULL)

  expect_equal(nrow(out), 1)
  expect_equal(out$module, "M1")
  expect_equal(out$status, "error")
  expect_equal(out$error_message, "HTTP 400")
  expect_equal(out$gene_count_sent, 0L)
})


test_that("regression: llm module order uses natural ordering for module='all'", {
  natural_order <- get(".hc_llm_natural_module_order", asNamespace("hcocena"))
  module_order <- get(".hc_llm_module_order", asNamespace("hcocena"))
  resolve_modules <- get(".hc_llm_resolve_modules", asNamespace("hcocena"))

  expect_equal(
    natural_order(c("M1", "M10", "M2", "M3")),
    c("M1", "M2", "M3", "M10")
  )

  hc <- methods::new("HCoCenaExperiment")
  hc@satellite <- S4Vectors::SimpleList(list(
    module_gene_list = data.frame(
      module = c("M1", "M10", "M2", "M3"),
      genes = c("A", "B", "C", "D"),
      stringsAsFactors = FALSE
    )
  ))

  expect_equal(
    module_order(hc),
    c("M1", "M2", "M3", "M10")
  )
  expect_equal(
    resolve_modules(hc, module = "all"),
    c("M1", "M2", "M3", "M10")
  )
})


test_that("regression: plot heatmaps default to main order unless clustering is enabled", {
  prepare_cols <- get(".hc_prepare_plot_heatmap_columns", asNamespace("hcocena"))

  mat <- matrix(
    c(1, 0, 2,
      2, 1, 0,
      3, 2, 1),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("M1", "M2", "M3"), c("T3", "T1", "T2"))
  )

  ordered <- prepare_cols(
    mat = mat,
    cluster_columns = FALSE,
    plot_order = NULL,
    main_order = c("T1", "T2", "T3"),
    fallback_order = c("T3", "T2", "T1"),
    context = "test heatmap"
  )
  expect_equal(colnames(ordered$mat), c("T1", "T2", "T3"))
  expect_null(ordered$col_dend)

  clustered <- prepare_cols(
    mat = mat,
    cluster_columns = TRUE,
    plot_order = c("T1", "T2", "T3"),
    main_order = c("T1", "T2", "T3"),
    fallback_order = c("T3", "T2", "T1"),
    context = "test heatmap"
  )
  expect_setequal(colnames(clustered$mat), c("T1", "T2", "T3"))
  expect_true(inherits(clustered$col_dend, "dendrogram"))

  expect_true("cluster_columns" %in% names(formals(hcocena:::.hc_upstream_inference_driver)))
  expect_true("cluster_columns" %in% names(formals(hcocena:::.hc_plot_enrichment_upstream_network_driver)))
  expect_true("cluster_columns" %in% names(formals(hcocena::hc_upstream_inference)))
  expect_true("cluster_columns" %in% names(formals(hcocena::hc_plot_enrichment_upstream_network)))
  expect_true("cluster_columns" %in% names(formals(hcocena::hc_plot_llm_enrichment)))
  expect_true("heatmap_cluster_columns" %in% names(formals(hcocena:::.hc_upstream_inference_driver)))
  expect_true("heatmap_cluster_columns" %in% names(formals(hcocena:::.hc_plot_enrichment_upstream_network_driver)))
  expect_true("heatmap_cluster_columns" %in% names(formals(hcocena::hc_upstream_inference)))
  expect_true("heatmap_cluster_columns" %in% names(formals(hcocena::hc_plot_enrichment_upstream_network)))
  expect_true("heatmap_cluster_columns" %in% names(formals(hcocena::hc_plot_llm_enrichment)))
  expect_identical(formals(hcocena:::.hc_plot_cluster_heatmap_driver)$cluster_columns, FALSE)
  expect_identical(formals(hcocena:::plot_cluster_heatmap_new)$cluster_columns, FALSE)
  expect_identical(formals(hcocena:::.hc_plot_cluster_heatmap_driver)$smart_column_gaps, FALSE)
  expect_identical(formals(hcocena:::plot_cluster_heatmap_new)$smart_column_gaps, FALSE)
  expect_true("column_gap_by" %in% names(formals(hcocena:::plot_cluster_heatmap_new)))
  expect_identical(formals(hcocena:::plot_cluster_heatmap_new)$column_gap_mm, 0.6)
  expect_identical(formals(hcocena:::replot_cluster_heatmap)$cluster_columns, FALSE)
  expect_identical(formals(hcocena:::.hc_change_grouping_parameter_driver)$cluster_columns, FALSE)
  expect_identical(formals(hcocena::hc_change_grouping_parameter)$cluster_columns, FALSE)
})

test_that("regression: duplicate heatmap column names keep layer-specific values when reordered", {
  prepare_cols <- get(".hc_prepare_plot_heatmap_columns", asNamespace("hcocena"))

  mat <- matrix(
    c(0.5, 1.5, -0.5, -1.5,
      0.2, 1.2, -0.2, -1.2),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("M1", "M2"), c("Moderate", "Severe", "Moderate", "Severe"))
  )

  ordered <- prepare_cols(
    mat = mat,
    cluster_columns = FALSE,
    plot_order = c("Moderate", "Severe", "Moderate", "Severe"),
    main_order = NULL,
    fallback_order = NULL,
    context = "duplicate condition heatmap"
  )

  expect_equal(as.numeric(ordered$mat["M1", ]), c(0.5, 1.5, -0.5, -1.5))
  expect_equal(as.numeric(ordered$mat["M2", ]), c(0.2, 1.2, -0.2, -1.2))
  expect_equal(colnames(ordered$mat), c("Moderate", "Severe", "Moderate", "Severe"))
})

test_that("regression: heatmap API alias helpers prefer new names and keep legacy aliases working", {
  resolve_col_order_alias <- get(".hc_resolve_col_order_alias", asNamespace("hcocena"))
  resolve_cluster_columns_alias <- get(".hc_resolve_cluster_columns_alias", asNamespace("hcocena"))

  expect_equal(
    resolve_col_order_alias(
      col_order = c("T2", "T1"),
      heatmap_col_order = NULL,
      col_order_missing = FALSE,
      heatmap_col_order_missing = TRUE
    ),
    c("T2", "T1")
  )
  expect_equal(
    resolve_col_order_alias(
      col_order = NULL,
      heatmap_col_order = c("T2", "T1"),
      col_order_missing = TRUE,
      heatmap_col_order_missing = FALSE
    ),
    c("T2", "T1")
  )
  expect_error(
    resolve_col_order_alias(
      col_order = c("T1", "T2"),
      heatmap_col_order = c("T2", "T1"),
      col_order_missing = FALSE,
      heatmap_col_order_missing = FALSE,
      context = "test"
    ),
    "Use either `col_order` or legacy `heatmap_col_order`"
  )

  expect_identical(
    resolve_cluster_columns_alias(
      cluster_columns = TRUE,
      heatmap_cluster_columns = NULL,
      cluster_columns_missing = FALSE,
      heatmap_cluster_columns_missing = TRUE
    ),
    TRUE
  )
  expect_identical(
    resolve_cluster_columns_alias(
      cluster_columns = FALSE,
      heatmap_cluster_columns = TRUE,
      cluster_columns_missing = TRUE,
      heatmap_cluster_columns_missing = FALSE
    ),
    TRUE
  )
  expect_error(
    resolve_cluster_columns_alias(
      cluster_columns = FALSE,
      heatmap_cluster_columns = TRUE,
      cluster_columns_missing = FALSE,
      heatmap_cluster_columns_missing = FALSE,
      context = "test"
    ),
    "Use either `cluster_columns` or legacy `heatmap_cluster_columns`"
  )
})


test_that("regression: celltype annotation matrix shows module labels instead of raw colors", {
  hc <- methods::new("HCoCenaExperiment")
  hc@integration@cluster <- S4Vectors::SimpleList()
  hc@integration@cluster[["module_prefix"]] <- "M"
  hc@integration@cluster[["module_label_map"]] <- c(yellow = "M1", blue = "M2")
  hc@satellite <- S4Vectors::SimpleList()
  hc@satellite[["celltype_annotation"]] <- list(
    selected_celltypes = base::data.frame(
      cluster = c("yellow", "blue"),
      cell_type = c("B cell", "T cell"),
      pct = c(55, 72),
      count = c(11, 9),
      stringsAsFactors = FALSE
    )
  )

  out <- hcocena::hc_plot_celltype_annotation_matrix(hc, return_data = TRUE)

  expect_equal(base::rownames(out$matrix), c("M1", "M2"))
  expect_setequal(base::unique(base::as.character(out$long_data$module_color)), c("yellow", "blue"))
  expect_true(inherits(out$plot, "ggplot"))
})

test_that("regression: longitudinal step2 and step3 accept prior step result lists", {
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE) ||
      !requireNamespace("SummarizedExperiment", quietly = TRUE) ||
      !requireNamespace("MultiAssayExperiment", quietly = TRUE) ||
      !requireNamespace("S4Vectors", quietly = TRUE)) {
    skip("longitudinal quick wrappers need Bioconductor core packages")
  }

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(c(1, 2), nrow = 1, dimnames = list("g1", c("s1", "s2")))),
    colData = S4Vectors::DataFrame(
      Subject = c("d1", "d2"),
      Time_token = c("T1", "T1"),
      row.names = c("s1", "s2")
    )
  )
  mae <- MultiAssayExperiment::MultiAssayExperiment(experiments = list(set1 = se))
  hc <- methods::new("HCoCenaExperiment")
  hc@mae <- mae
  hc@config <- methods::new("HCoCenaConfig")
  hc@satellite <- S4Vectors::SimpleList(list(
    longitudinal_endotypes = list(
      cap_matrix = matrix(
        c(1, 0, 0, 1),
        nrow = 2,
        dimnames = list(c("d1", "d2"), c("M1__1", "M1__2"))
      ),
      module_cluster_matrix = matrix(
        c(1, 2),
        nrow = 2,
        dimnames = list(c("d1", "d2"), c("M1"))
      )
    )
  ))

  with_mocked_bindings(
    .hc_run_longitudinal_step2_graph = function(hc, slot_name, ...) {
      sat <- as.list(hc@satellite)
      sat[[slot_name]]$meta_cluster <- data.frame(
        donor = c("d1", "d2"),
        meta_cluster = c("MC1", "MC2"),
        stringsAsFactors = FALSE
      )
      sat[[slot_name]]$meta_method_comparison <- data.frame(method = "knn", stringsAsFactors = FALSE)
      sat[[slot_name]]$meta_score_table <- data.frame(k = 2, score = 1, stringsAsFactors = FALSE)
      hc@satellite <- S4Vectors::SimpleList(sat)
      hc
    },
    hc_plot_longitudinal_meta_embeddings = function(hc, ...) {
      list(pca = "pca_plot", umap = "umap_plot", cross_tab = "cross_tab", tables = list())
    },
    hc_plot_longitudinal_meta_module_waves = function(hc, ...) {
      list(meta_module_waves = "meta_waves_plot")
    },
    {
      step1_res <- list(hc = hc, plots = list(), diagnostics = list())
      step2_res <- hcocena::hc_longitudinal_step2_meta_clustering(step1_res)
      expect_s4_class(step2_res$hc, "HCoCenaExperiment")
      expect_named(step2_res$plots, "longitudinal_endotypes")
      expect_equal(step2_res$plots$longitudinal_endotypes$pca, "pca_plot")

      step3_res <- hcocena::hc_longitudinal_step3_meta_module_trajectories(step2_res)
      expect_s4_class(step3_res$hc, "HCoCenaExperiment")
      expect_named(step3_res$plots, "longitudinal_endotypes")
      expect_equal(step3_res$plots$longitudinal_endotypes$meta_module_waves, "meta_waves_plot")
    }
  )
})

test_that("regression: longitudinal step1 accepts prior step result lists", {
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE) ||
      !requireNamespace("SummarizedExperiment", quietly = TRUE) ||
      !requireNamespace("MultiAssayExperiment", quietly = TRUE) ||
      !requireNamespace("S4Vectors", quietly = TRUE)) {
    skip("longitudinal quick wrappers need Bioconductor core packages")
  }

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(c(1, 2), nrow = 1, dimnames = list("g1", c("s1", "s2")))),
    colData = S4Vectors::DataFrame(
      Subject = c("d1", "d2"),
      Time_token = c("T1", "T2"),
      row.names = c("s1", "s2")
    )
  )
  mae <- MultiAssayExperiment::MultiAssayExperiment(experiments = list(set1 = se))
  hc <- methods::new("HCoCenaExperiment")
  hc@mae <- mae
  hc@config <- methods::new("HCoCenaConfig")
  hc@satellite <- S4Vectors::SimpleList()

  with_mocked_bindings(
    .hc_longitudinal_step1_run_single_direct = function(hc, ...) {
      list(hc = hc, plots = list(ok = TRUE), diagnostics = list(ok = TRUE))
    },
    {
      wrapped <- list(hc = hc, plots = list(), diagnostics = list())
      out <- hcocena::hc_longitudinal_step1_module_donor(
        wrapped,
        donor_col = "Subject",
        time_col = "Time_token",
        time_levels = c("T1", "T2"),
        method = "kmeans",
        nstart = 1,
        cap_runs = 1,
        score_method = "calinski_harabasz",
        scale_features = FALSE
      )
      expect_s4_class(out$hc, "HCoCenaExperiment")
      expect_true(isTRUE(out$plots$ok))
    }
  )
})

test_that("regression: htmlwidgets display inline during HTML knitting", {
  skip_if_not_installed("htmlwidgets")
  skip_if_not_installed("knitr")

  display_fun <- get(".hc_display_object", asNamespace("hcocena"))
  widget <- htmlwidgets::createWidget(
    name = "hcocena-test-widget",
    x = list(value = 1),
    package = "htmlwidgets"
  )

  old_options <- options(knitr.in.progress = TRUE)
  old_knit <- knitr::opts_knit$get()
  old_current <- knitr::opts_current$get()
  on.exit({
    options(old_options)
    do.call(knitr::opts_knit$set, old_knit)
    do.call(knitr::opts_current$set, old_current)
  }, add = TRUE)

  knitr::opts_knit$set(rmarkdown.pandoc.to = "html")
  knitr::opts_current$set(out.width.px = "100%", out.height.px = "400px")

  out <- capture.output(display_fun(widget))

  expect_true(any(grepl("hcocena-test-widget html-widget", out, fixed = TRUE)))
  expect_true(length(knitr::knit_meta()) > 0)
})

test_that("regression: htmlwidgets dispatch through the print generic for RStudio inline rendering", {
  display_src_path <- test_path("..", "..", "R", "hc_display_helpers.R")
  if (!file.exists(display_src_path)) {
    skip("Source file is not available in the installed-package test context.")
  }
  display_src <- paste(readLines(display_src_path, warn = FALSE), collapse = "\n")

  # RStudio's notebook renders htmlwidgets inline only when they go through the
  # `print` generic (this works even from inside a function). Calling
  # print.htmlwidget() directly (e.g. via getFromNamespace) bypasses that hook
  # and the widget lands in the Viewer pane instead of inline below the chunk.
  expect_true(grepl("print(x)", display_src, fixed = TRUE))
  expect_false(grepl("getFromNamespace(\"print.htmlwidget\"", display_src, fixed = TRUE))
})

test_that("regression: htmlwidgets emit raw HTML only when explicitly opted in", {
  skip_if_not_installed("htmlwidgets")
  skip_if_not_installed("repr")

  display_fun <- get(".hc_display_object", asNamespace("hcocena"))
  widget <- htmlwidgets::createWidget(
    name = "hcocena-test-widget",
    x = list(value = 1),
    package = "htmlwidgets"
  )

  old_options <- options(
    rstudio.notebook.executing = TRUE,
    viewer = function(...) stop("viewer should not be used", call. = FALSE),
    hcocena.htmlwidget_display = "auto",
    hcocena.htmlwidget_emit_html = TRUE
  )
  on.exit(options(old_options), add = TRUE)

  out <- capture.output(display_fun(widget))

  expect_true(any(grepl("hcocena-test-widget html-widget", out, fixed = TRUE)))
  expect_true(any(grepl("<!doctype html>", tolower(out), fixed = TRUE)))
})
