# hcocena

[![R-CMD-check](https://github.com/BioCompNet/hcocena/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/BioCompNet/hcocena/actions/workflows/R-CMD-check.yaml)

`hcocena` is an R package for horizontal integration and downstream analysis of
transcriptomics datasets. It provides a modern S4 workflow built around
`HCoCenaExperiment` for reproducible network-centric transcriptomics analyses.

![hcocena overview](.github/assets/hcocena-overview.jpg)

The package supports both multi-layer integration, such as RNA-seq plus array
data, and single-layer analyses using the same API. The focus is a
module-centric workflow: from data import and correlation-based network
construction to clustering, heatmaps, functional enrichment, upstream
inference, cell-type annotation, longitudinal analysis, and optional
LLM-assisted module interpretation.

## What hcocena provides

- S4-first workflow with `HCoCenaExperiment`, `MultiAssayExperiment`, and
  `SummarizedExperiment`
- Correlation cutoff tuning and automatic cutoff selection helpers
- Clustering, integrated network construction, module splitting, and hCoCena
  heatmaps
- Functional enrichment across multiple databases with export helpers
- Upstream inference with DoRothEA and PROGENy via `decoupleR`
- Cell-type annotation helpers and reference-data preview utilities
- Longitudinal module and endotype analyses
- A Docker workflow with bundled `reference_files` for a ready-to-run setup

## Repository structure

- Package source is at the repository root and follows a Bioconductor-style
  layout
- Docker support lives in [`docker/`](docker), including bundled
  `reference_files`
- GitHub-only workflow notebooks are kept in [`inst/scripts/workflows/`](inst/scripts/workflows/)
- CI for package checks is defined in
  [`.github/workflows/bioc-check.yaml`](.github/workflows/bioc-check.yaml)

## Contributors

- [Waqar Hanif](https://github.com/waqarhanif-biocode)

## Installation

After Bioconductor acceptance:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("hcocena")
```

For local development or pre-submission testing from a checkout:

```r
install.packages("remotes")
remotes::install_local(".", dependencies = TRUE, upgrade = "never")
```

## Docker

To build a ready-to-use RStudio image from this repository:

```bash
docker build -f docker/Dockerfile -t hcocena .
docker run --rm -p 8787:8787 -e PASSWORD=hcocena hcocena
```

The container prepares a workspace at `/home/rstudio/hcocena` and includes:

- the local `hcocena` installation
- bundled `reference_files/` with pathway, GO, hallmark, TF, and immune helper references
- visible workflow notebooks under `/home/rstudio/hcocena/inst/scripts/workflows/`
  including `hcocena_main.Rmd` and `hcocena_satellite.Rmd`
- preinstalled optional packages for common workflows, including
  `DESeq2`, `limma`, `sva`, `edgeR`, `tximport`, `apeglm`, `ashr`,
  `EnhancedVolcano`, `pheatmap`, `org.Hs.eg.db`, `org.Mm.eg.db`,
  `CALIBERrfimpute`, `RCy3`, `SpatialExperiment`, and `GSVA`
- empty `count_data`, `annotation_data`, and `output` directories

See [`docker/README.md`](docker/README.md) for the Docker-specific notes.

## Minimal S4 workflow

```r
library(hcocena)

hc <- hc_init()
hc <- hc_set_paths(
  hc,
  dir_count_data = "/path/to/counts/",
  dir_annotation = "/path/to/annotation/",
  dir_reference_files = "/path/to/reference/",
  dir_output = tempdir()
)
hc <- hc_define_layers(
  hc,
  data_sets = list(
    Layer1 = c("counts.tsv", "anno.tsv")
  )
)
hc <- hc_read_data(
  hc,
  gene_symbol_col = "SYMBOL",
  sample_col = "SampleID",
  count_has_rn = FALSE,
  anno_has_rn = FALSE
)
hc <- hc_run_expression_analysis_1(hc, export = FALSE)
hc <- hc_plot_cutoffs(hc, interactive = FALSE)
```

For a reproducible package-based example, install `hcocena` and run:

```r
browseVignettes("hcocena")
```

The package ships toy data and prepared example objects in `inst/extdata` to
support documentation, testing, and manual smoke tests.

## Real-data regression checks

Local real-data checks live behind an explicit opt-in runner. By default, the
runner looks for the STAR protocol data in `../data` relative to this repository
and writes ignored outputs to `realdata-output/`.

`quick` uses the real Array/RNA-seq data with the top 2000 variable genes per
layer. It runs integration, clustering, module heatmaps, module splitting,
Hallmark/KEGG enrichment, and a module-label font-size probe. `full` keeps the
same checks but uses the full imported gene set and broader enrichment defaults.

```bash
Rscript scripts/run_realdata_regression.R --mode quick
Rscript scripts/run_realdata_regression.R --mode quick --update-reference
Rscript scripts/run_realdata_regression.R --mode full
```

Set `HCOCENA_REALDATA_DIR` to point at another data directory. Normal package
tests skip the real-data run; enable it explicitly with:

```bash
HCOCENA_RUN_REALDATA=true HCOCENA_REALDATA_MODE=quick Rscript -e "testthat::test_file('tests/testthat/test-realdata-regression.R')"
```

PowerShell:

```powershell
$env:HCOCENA_RUN_REALDATA = "true"; $env:HCOCENA_REALDATA_MODE = "quick"; Rscript -e "testthat::test_file('tests/testthat/test-realdata-regression.R')"
```

The exported CSV/JSON artifacts are designed as golden-master inputs for later
R/Python comparisons: QC metrics, edge lists, GFC matrices, module tables,
split diagnostics, enrichment tables, heatmap matrices, and a module-label
font-size probe. Each run also writes a `visual_check_report_<mode>.pdf` with
the generated plot variants, non-default plot parameters, and image pages for
manual inspection, including the pre-split heatmap, post-split heatmap, and
combined enrichment plot.

## Analysis templates

The longer walkthroughs the package is normally driven from are installed with
it, so they are available from an installed copy and not only from a clone:

- `inst/scripts/workflows/hcocena_main.Rmd` -- full analysis, import to enrichment
- `inst/scripts/workflows/hcocena_satellite.Rmd` -- optional downstream analyses

```r
dir(system.file("scripts", "workflows", package = "hcocena"))

file.copy(
  system.file("scripts", "workflows", "hcocena_main.Rmd", package = "hcocena"),
  "hcocena_main.Rmd"
)
```

They are templates, not reproducible documents: they point at your own count
and annotation files, and several steps contact external services (ChEA3,
Enrichr, DoRothEA, Cytoscape, LLM providers). The built vignette
(`vignette("hcocena-s4-workflow")`) is the runnable short version.

One practical note: when using longitudinal imputation with
`impute_method = "rfcont"`, attach `CALIBERrfimpute` in the session first:

```r
library(CALIBERrfimpute)
```

## Documentation and references

- Method paper: Oestreich et al. (2022), *Bioinformatics*,
  https://doi.org/10.1093/bioinformatics/btac589
- Step-by-step protocol: Holsten et al. (2024), *STAR Protocols*,
  https://doi.org/10.1016/j.xpro.2024.102922

### Coming from the papers

Both papers use the earlier hCoCena interface, in which the functions worked
on a global `hcobject` and had no `hc_` prefix (repository
`MarieOestreich/hCoCena`, Docker image `mo126/hcocena`). In this package every
step takes and returns an `HCoCenaExperiment` (`hc <- hc_step(hc, ...)`), and
the current Docker image is `therealtomek/hcocena`. The functions map as
follows:

| In the papers | In this package |
|---|---|
| `init_wd()`, `check_dirs()`, `init_save_folder()` | `hc_init()`, `hc_set_paths()`, `hc_check_dirs()`, `hc_init_save_folder()` |
| `define_layers()`, `read_data()` | `hc_define_layers()`, `hc_read_data()` |
| `set_supp_files()`, `read_supplementary()` | `hc_set_supp_files()`, `hc_read_supplementary()` |
| `set_global_settings()`, `set_layer_settings()` | `hc_set_global_settings()`, `hc_set_layer_settings()` |
| `suggest_topvar()` | `hc_suggest_topvar()` |
| `run_expression_analysis_1()`, `plot_cutoffs()`, `set_cutoff()`, `plot_deg_dist()` | `hc_run_expression_analysis_1()`, `hc_plot_cutoffs()`, `hc_set_cutoff()`, `hc_plot_deg_dist()` |
| `run_expression_analysis_2()` | `hc_run_expression_analysis_2()` |
| `build_integrated_network()`, `plot_integrated_network()` | `hc_build_integrated_network()`, `hc_plot_integrated_network()` |
| `cluster_calculation()`, `plot_cluster_heatmap()` | `hc_cluster_calculation()`, `hc_plot_cluster_heatmap()` |
| `functional_enrichment()` | `hc_functional_enrichment()` |
| `TF_overrep_module()`, `TF_overrep_network()`, `check_tf()` | `hc_tf_overrep_module()`, `hc_tf_overrep_network()`, `hc_check_tf()` |
| `find_hubs()` | `hc_find_hubs()` |
| `change_grouping_parameter()`, `cut_hclust()` | `hc_change_grouping_parameter()`, `hc_cut_hclust()` |
| `write_session_info()` | `hc_write_session_info()` |

Results that the papers read from `hcobject` are available through
accessors such as `hc_clusters()`, `hc_gene_to_cluster()`, `hc_graph()` and
`hc_satellite()`.

## Citation

Please cite the method paper, and the protocol where you followed it:

Oestreich M, Holsten L, Agrawal S, Dahm K, Koch P, Jin H, Becker M, Ulas T
(2022). "hCoCena: horizontal integration and analysis of transcriptomics
datasets." *Bioinformatics* 38(20):4727-4734.
doi:10.1093/bioinformatics/btac589

Holsten L, Dahm K, Oestreich M, Becker M, Ulas T (2024). "hCoCena: A toolbox
for network-based co-expression analysis and horizontal integration of
transcriptomic datasets." *STAR Protocols* 5(1):102922.
doi:10.1016/j.xpro.2024.102922

In R: `citation("hcocena")`.
