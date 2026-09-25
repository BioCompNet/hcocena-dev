#' Load an example `HCoCenaExperiment`
#'
#' Returns one of the small example objects shipped in `inst/extdata`, taken
#' at a given stage of the standard workflow, so that individual steps can be
#' tried without running everything before them. The data are two simulated
#' layers of 40 marker genes for four blood cell types over 32 samples; see
#' `inst/scripts/make-fixtures.R` for the generating code.
#'
#' The output directory stored in the object is replaced by a fresh
#' directory inside [tempdir()], so functions that write plots or tables
#' never write outside the session's temporary directory.
#'
#' @param stage Workflow stage to load:
#'   * `"prepared"`: data read and all settings made;
#'   * `"after_part1"`: correlations computed and the cutoff chosen
#'     ([hc_run_expression_analysis_1()], [hc_set_cutoff()]);
#'   * `"after_part2"`: group fold changes computed
#'     ([hc_run_expression_analysis_2()]);
#'   * `"clustered"`: integrated network built and modules detected
#'     ([hc_build_integrated_network()], [hc_cluster_calculation()]).
#' @param dir_output Output directory to store in the object. Defaults to a
#'   new directory inside [tempdir()].
#' @examples
#' hc <- hc_example_data("clustered")
#' hc
#' table(hc_gene_to_cluster(hc)$color)
#' @return A `HCoCenaExperiment`.
#' @export
hc_example_data <- function(stage = c("prepared", "after_part1",
                                      "after_part2", "clustered"),
                            dir_output = NULL) {
  stage <- base::match.arg(stage)
  # Unprefixed on purpose: pkgload shims `system.file()` for development
  # loads, and `base::system.file()` would bypass that shim.
  path <- system.file(
    "extdata", base::paste0("hc_", stage, ".rds"),
    package = "hcocena", mustWork = TRUE
  )
  hc <- base::readRDS(path)

  if (base::is.null(dir_output)) {
    dir_output <- base::tempfile(pattern = "hcocena-example-")
  }
  dir_output <- .hc_normalize_path_value(dir_output, "dir_output")
  base::dir.create(dir_output, recursive = TRUE, showWarnings = FALSE)

  paths <- hc@config@paths
  paths$dir_output <- .hc_path_with_trailing_slash(dir_output)
  hc@config@paths <- paths
  methods::validObject(hc)
  hc
}
