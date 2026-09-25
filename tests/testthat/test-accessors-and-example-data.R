## Accessors that replace direct slot access, the example-data loader, and
## fixes that surfaced while writing runnable examples for them.

toy_gmt <- function() {
  path <- system.file("extdata", "toy_celltype_markers.gmt", package = "hcocena")
  skip_if(!nzchar(path), "The toy marker GMT is unavailable.")
  path
}

quietly <- function(expr) suppressMessages(suppressWarnings(expr))

test_that("hc_example_data() points the output into a fresh temporary directory", {
  hc <- hc_example_data("clustered")
  out <- as.character(as.data.frame(hc_config(hc)@paths)$dir_output)
  tmp <- normalizePath(tempdir(), winslash = "/")

  expect_s4_class(hc, "HCoCenaExperiment")
  expect_true(dir.exists(out))
  expect_true(startsWith(normalizePath(out, winslash = "/"), tmp))
  expect_match(out, "/$")
  expect_error(hc_example_data("no_such_stage"))
})

test_that("hc_graph() and hc_satellite() return the stored objects", {
  hc <- hc_example_data("clustered")

  expect_identical(hc_graph(hc), hc@integration@graph)
  expect_true(inherits(hc_graph(hc), "igraph"))
  expect_identical(hc_satellite(hc), hc@satellite)
  expect_identical(
    hc_satellite(hc, "cutoff_selection"),
    hc@satellite[["cutoff_selection"]]
  )
  expect_null(hc_satellite(hc, "no_such_entry"))
  expect_error(hc_satellite(hc, c("a", "b")), "single character string")
  expect_null(hc_graph(hc_init()))
})

test_that("names given to custom GMT files become the database labels", {
  files <- hcocena:::.hc_normalize_custom_gmt_files(
    c(CellTypes = "a.gmt", "b.gmt"),
    default_prefix = "CustomEnrichment"
  )
  expect_identical(names(files), c("CellTypes", "CustomEnrichment2"))

  files <- hcocena:::.hc_normalize_custom_gmt_files(list(Markers = "a.gmt"))
  expect_identical(names(files), "Markers")
})

test_that("cell-type annotation keeps hits whose Storey q-value is undefined", {
  skip_if_not_installed("clusterProfiler")
  # Each module overlaps exactly one marker set, so every test sees a single
  # p-value and qvalue::pi0est() cannot fit - the q-value comes back NA.
  hc <- hc_example_data("clustered")
  hc <- quietly(hc_celltype_annotation(
    hc,
    databases = character(0),
    custom_gmt_files = c(CellTypes = toy_gmt()),
    mode = "fine",
    export_excel = FALSE
  ))
  sel <- as.data.frame(hc_satellite(hc, "celltype_annotation")$selected_celltypes)

  expect_gt(nrow(sel), 0)
  expect_setequal(
    unique(sel$cell_type),
    paste("[CellTypes]", c("T cells", "B cells", "Monocytes", "NK cells"))
  )
})

test_that("upstream inference keeps same-named groups of different layers apart", {
  skip_if_not_installed("decoupleR")
  # Both example layers use the groups control / mild / severe.
  hc <- hc_example_data("clustered")
  hc <- quietly(hc_upstream_inference(
    hc,
    resources = "Pathway",
    custom_pathway_gmt = c(CellTypes = toy_gmt()),
    plot = FALSE,
    save_pdf = FALSE
  ))
  res <- hc_satellite(hc, "upstream_inference")
  conds <- unique(as.character(as.data.frame(res$all_upstream_by_condition)$condition))

  expect_length(conds, 6)
  expect_false(anyDuplicated(conds) > 0)
})

test_that("hc_check_tf() says what to run first", {
  hc <- hc_example_data("clustered")
  expect_error(quietly(hc_check_tf(hc, TF = "PAX5")), "hc_tf_overrep_network")
})
