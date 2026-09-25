# Provenance for the example fixtures in inst/extdata.
#
# Builds the four small HCoCenaExperiment objects that the man-page examples
# run against, capturing the workflow at prepared / after-part1 / after-part2 /
# clustered. Run from the package root with the native Rscript:
#
#   Rscript inst/scripts/make-fixtures.R
#
# Design notes
# ------------
# * 16 donors x 2 timepoints = 32 samples per layer. The donor count is what
#   the longitudinal meta-clustering needs: its UMAP step requires more points
#   than neighbours, so a handful of donors is not enough.
# * Genes are simulated in four correlated blocks so that clustering actually
#   finds several modules on such a small matrix. Each block is a set of real
#   marker symbols for one blood cell type, so the resulting modules are
#   recognizable to the cell-type annotation and enrichment functions rather
#   than being unmatchable placeholder names.
# * Two stored plots are dropped from every layer before saving:
#   `cutoff_calc_out$dd_plot_calculated_optimal` (a ggplot, ~16 MB per layer in
#   memory) and `heatmap_out$heatmap` (a ComplexHeatmap, ~0.7 MB). Both are
#   written by the pipeline and never read back - hc_plot_deg_dist() rebuilds
#   the degree-distribution plot from `cutoff_stats`, and nothing consumes the
#   stored Heatmap at all. Keeping them costs ~1.6 MB of tarball for nothing.
# * Noise sd 1.2 puts within-block correlations near 0.8 rather than 0.99. At
#   0.99 the whole degree distribution changes only in the last thousandth of
#   the cutoff grid and the automatic cutoff can collapse a layer to a single
#   edge.
# * saveRDS() cannot write these objects: DelayedArray's
#   .S4_object_contains_out_of_memory_data hook walks the S4 slots recursively
#   and overflows the node stack. serialize() to a gzfile writes the same
#   format without that hook.

suppressPackageStartupMessages({
  library(MultiAssayExperiment)
  library(SummarizedExperiment)
})
pkgload::load_all(".", quiet = TRUE)

set.seed(42)
extdata <- file.path("inst", "extdata")
out_dir <- tempfile("hc-fixtures-")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

N_DONORS <- 16L
TIMEPOINTS <- c("T1", "T2")
# Four blocks of ten markers: T cells, B cells, monocytes, NK cells.
GENE_BLOCKS <- list(
  T_cell = c("CD3D", "CD3E", "CD3G", "CD2", "CD28",
             "LCK", "ZAP70", "IL7R", "CD7", "TRAC"),
  B_cell = c("CD19", "MS4A1", "CD79A", "CD79B", "BLNK",
             "PAX5", "CR2", "FCRL1", "TNFRSF13B", "VPREB3"),
  Monocyte = c("CD14", "LYZ", "FCN1", "VCAN", "S100A8",
               "S100A9", "CSF1R", "ITGAM", "CD68", "FCGR3A"),
  NK_cell = c("NKG7", "GNLY", "KLRD1", "KLRF1", "PRF1",
              "GZMB", "NCR1", "KLRC1", "FGFBP2", "SPON2")
)
GENES <- unlist(GENE_BLOCKS, use.names = FALSE)
BLOCK_OF <- rep(seq_along(GENE_BLOCKS), lengths(GENE_BLOCKS))
N_GENES <- length(GENES)
N_BLOCKS <- length(GENE_BLOCKS)
NOISE_SD <- 1.2
# Offset large enough that no simulated count comes out negative: several
# functions reject negative expression as mis-set `data_in_log`.
BASELINE <- 25

make_layer <- function(prefix, seed) {
  set.seed(seed)
  samples <- paste0(prefix, "_", rep(sprintf("D%02d", seq_len(N_DONORS)),
                                     each = length(TIMEPOINTS)),
                    "_", rep(TIMEPOINTS, times = N_DONORS))
  n <- length(samples)
  latent <- matrix(stats::rnorm(N_BLOCKS * n), nrow = N_BLOCKS)
  counts <- t(vapply(seq_len(N_GENES), function(g) {
    latent[BLOCK_OF[[g]], ] * 3 + stats::rnorm(n, sd = NOISE_SD) + BASELINE
  }, numeric(n)))
  dimnames(counts) <- list(GENES, samples)
  stopifnot(all(counts > 0))

  anno <- data.frame(
    SampleID = samples,
    # Three groups, not two: hc_meta_correlation_cat() correlates a variable
    # across the groups of the variable of interest, which needs at least
    # three points to be defined.
    group = rep(c("control", "mild", "severe"), length.out = n),
    batch = rep(c("B1", "B2"), each = 2, length.out = n),
    donor = rep(sprintf("D%02d", seq_len(N_DONORS)), each = length(TIMEPOINTS)),
    timepoint = rep(TIMEPOINTS, times = N_DONORS),
    # one numeric covariate, so hc_meta_correlation_num() has something to
    # correlate module scores against
    age = rep(round(stats::runif(N_DONORS, 20, 70)), each = length(TIMEPOINTS)),
    row.names = samples,
    stringsAsFactors = FALSE
  )
  list(counts = counts, anno = anno)
}

l1 <- make_layer("L1", 1)
l2 <- make_layer("L2", 2)
set1_counts <- l1$counts; set1_anno <- l1$anno
set2_counts <- l2$counts; set2_anno <- l2$anno

# The same two layers as tab-separated files. The vignette reads these rather
# than the serialized objects, so that it demonstrates the import path a user
# actually takes - and so that vignette and fixtures show one dataset, not two.
write_tsv_layer <- function(counts, anno, idx) {
  cnt <- data.frame(SYMBOL = rownames(counts), counts,
                    check.names = FALSE, stringsAsFactors = FALSE)
  utils::write.table(cnt, file.path(extdata, sprintf("toy_layer%d_counts.tsv", idx)),
                     sep = "	", quote = FALSE, row.names = FALSE)
  utils::write.table(anno, file.path(extdata, sprintf("toy_layer%d_anno.tsv", idx)),
                     sep = "	", quote = FALSE, row.names = FALSE)
}
write_tsv_layer(round(set1_counts, 3), set1_anno, 1)
write_tsv_layer(round(set2_counts, 3), set2_anno, 2)

# ---- drop the unread ggplot before writing -------------------------------
strip_plots <- function(hc) {
  for (nm in names(hc@layer_results)) {
    lr <- hc@layer_results[[nm]]
    p1 <- as.list(lr@part1)
    if (!is.null(p1$cutoff_calc_out)) {
      p1$cutoff_calc_out$dd_plot_calculated_optimal <- NULL
    }
    lr@part1 <- S4Vectors::SimpleList(p1)
    p2 <- as.list(lr@part2)
    if (!is.null(p2$heatmap_out)) p2$heatmap_out$heatmap <- NULL
    lr@part2 <- S4Vectors::SimpleList(p2)
    hc@layer_results[[nm]] <- lr
  }
  hc
}

write_fixture <- function(hc, name) {
  hc <- strip_plots(hc)
  path <- file.path(extdata, name)
  con <- gzfile(path, "wb")
  on.exit(close(con), add = TRUE)
  serialize(hc, con)
  invisible(path)
}

# ---- build the workflow --------------------------------------------------
hc <- hc_init()
hc <- hc_set_paths(hc, dir_count_data = FALSE, dir_annotation = FALSE,
                   dir_reference_files = FALSE, dir_output = out_dir)
hc <- hc_check_dirs(hc)
hc <- hc_init_save_folder(hc, name = "")
# The list names become the layer names used in file names and plot titles;
# the internal layer_results keys stay set1/set2 regardless.
hc <- hc_define_layers(hc, data_sets = list(
  Layer1 = c("set1_counts", "set1_anno"),
  Layer2 = c("set2_counts", "set2_anno")
))
hc <- hc_read_data(hc, count_has_rn = TRUE, anno_has_rn = TRUE,
                   auto_setup_output = FALSE)
hc <- hc_set_global_settings(
  hc, organism = "human", control_keyword = "none",
  variable_of_interest = "group",
  min_nodes_number_for_network = 1, min_nodes_number_for_cluster = 1,
  range_GFC = 2, layout_algorithm = "layout_with_stress", data_in_log = FALSE
)
hc <- hc_set_layer_settings(hc, top_var = c(N_GENES, N_GENES),
                            min_corr = c(0.5, 0.5),
                            range_cutoff_length = c(50, 50),
                            print_distribution_plots = c(FALSE, FALSE))
write_fixture(hc, "hc_prepared.rds")

hc <- hc_run_expression_analysis_1(hc, corr_method = "pearson")
hc <- hc_set_cutoff(hc, auto = TRUE)
write_fixture(hc, "hc_after_part1.rds")

hc <- hc_run_expression_analysis_2(hc, plot_HM = FALSE)
write_fixture(hc, "hc_after_part2.rds")

hc <- hc_build_integrated_network(hc, mode = "u", multi_edges = "min")
hc <- hc_cluster_calculation(hc, cluster_algo = "cluster_leiden",
                             no_of_iterations = 2, resolution = 1)
write_fixture(hc, "hc_clustered.rds")

# The four marker panels as a GMT, for the enrichment, cell-type and upstream
# examples. Each panel has exactly the ten genes of its simulated block.
gmt_names <- c(T_cell = "T cells", B_cell = "B cells",
               Monocyte = "Monocytes", NK_cell = "NK cells")
writeLines(
  vapply(names(GENE_BLOCKS), function(b) {
    paste(c(gmt_names[[b]], "-", GENE_BLOCKS[[b]]), collapse = "\t")
  }, character(1)),
  file.path(extdata, "toy_celltype_markers.gmt")
)

for (f in c("toy_layer1_counts.tsv", "toy_layer1_anno.tsv",
            "toy_layer2_counts.tsv", "toy_layer2_anno.tsv",
            "toy_celltype_markers.gmt")) {
  message(sprintf("%-22s %5.0f KB", f, file.size(file.path(extdata, f)) / 1024))
}
for (f in c("hc_prepared.rds", "hc_after_part1.rds",
            "hc_after_part2.rds", "hc_clustered.rds")) {
  p <- file.path(extdata, f)
  chk <- readRDS(p)
  stopifnot(inherits(chk, "HCoCenaExperiment"))
  message(sprintf("%-22s %5.0f KB", f, file.size(p) / 1024))
}
ci <- as.data.frame(as.list(readRDS(file.path(extdata, "hc_clustered.rds"))@
                              integration@cluster)[["cluster_information"]])
message("Module im clustered-Fixture: ", nrow(ci))
