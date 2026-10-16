# hcocena Docker image

Public RStudio-based image for `hcocena`, the R package for horizontal
integration and downstream analysis of transcriptomics datasets.

## What is included

- the `hcocena` package preinstalled in the container, with its vignette:
  `vignette("hcocena-s4-workflow")` runs the whole analysis on bundled toy data
- a project workspace in your RStudio home, set up at the first start (also
  when a platform mounts its own home directory):
  - `README_hcocena.md`
  - `reference_files/` with pathway, GO, hallmark, TF, and immune helper
    references, shared by all projects
  - `projects/example_project/` with the workflow notebooks
    `01_hcocena_main.Rmd` and `02_hcocena_satellite.Rmd` and empty
    `count_data/`, `annotation_data/`, `project_data/` and `output/` folders;
    the notebook paths are relative to the project folder
- further projects: copy `example_project/` anywhere, e.g. into a mounted
  folder; the reference files are found from any location
- the unmodified templates for use outside the container, at
  `system.file("scripts", "workflows", package = "hcocena")`
- RNA-seq and differential-expression packages including `DESeq2`, `limma`,
  `sva`, `edgeR`, `tximport`, `apeglm`, `ashr`, `EnhancedVolcano`,
  `pheatmap`, `org.Hs.eg.db`, and `org.Mm.eg.db`

## Recommended tags

- `latest` for the current image and quick-start commands
- `0.99.12` for a pinned, reproducible setup

Older tags still kept for older reproducible runs:

- `0.99.10`
- `0.99.7`
- `1.100`
- `1.99`
- `1.98`
- `1.97`
- `1.96`
- `1.95`
- `1.94`
- `1.9`
- `1.28`
- `1.1.2`

## Quick start

Pull the current image:

```bash
docker pull therealtomek/hcocena:latest
```

Run RStudio Server:

```bash
docker run --rm -p 8787:8787 -e PASSWORD=hcocena therealtomek/hcocena:latest
```

For reproducible runs, replace `latest` with a pinned tag such as `0.99.12`.

Then open:

```text
http://localhost:8787
```

Login:

- user: `rstudio`
- password: the value passed in `PASSWORD`

## Project links

- Package (Bioconductor submission): https://github.com/BioCompNet/hcocena
- Docker setup and development repository: https://github.com/BioCompNet/hcocena-dev
- Docker Hub tags: https://hub.docker.com/r/therealtomek/hcocena/tags

## Notes

The image is intended as a ready-to-use environment for exploring the package
and running the main and satellite workflows in RStudio.
