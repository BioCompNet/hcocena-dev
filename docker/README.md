# Docker

Build the image from the repository root:

```bash
docker build -f docker/Dockerfile -t hcocena .
```

The image includes:

- the local `hcocena` package installation (the source tree is removed after
  the install, so users only see their own home)
- a project workspace in the home of the RStudio user:
  - `~/README_workspace.md`
  - `~/reference_files/` with pathway, GO, hallmark, TF, and immune helper
    references, shared by all hCoCena projects
  - `~/projects/hcocena_example/` with the TB-neo-specific
    `01_hcocena_main.Rmd` and `02_hcocena_satellite.Rmd` and empty
    `count_data/`, `annotation_data/`, `project_data/` and `output/` folders
  - `~/projects/deseq2_example/` with the RNA-seq DESeq2 workflow from
    `docker/deseq2_workflow/` (workflow Rmd, functions file, cheat sheet);
    its `project_dir` is set to the notebook folder during the build
- further projects: copy an example project anywhere (e.g. into a mounted
  folder); the notebook paths are relative to the project folder
- preinstalled optional packages for common workflows, including
  `graphlayouts`, `ellmer`, `DESeq2`, `limma`, `sva`, `edgeR`, `tximport`,
  `apeglm`, `ashr`, `EnhancedVolcano`, `pheatmap`, `org.Hs.eg.db`,
  `org.Mm.eg.db`, `CALIBERrfimpute`, `decoupleR`, `OmnipathR`, `progeny`,
  `RCy3`, `SpatialExperiment`, and `GSVA`, and the optional packages of the
  DESeq2 workflow (`IHW`, `msigdbr` with `msigdbrdata`, `fgsea`, `UpSetR`, `vsn`,
  `multcomp`, `lme4`)

The DESeq2 workflow is not part of this repository (`docker/deseq2_workflow/`
is git-ignored). Before building, copy the lab's current version (workflow Rmd,
functions file, cheat sheet) into `docker/deseq2_workflow/`; the build fails
without it.

The workspace templates live in `/usr/local/share/hcocena`. The home is filled
from them at container start (`/etc/cont-init.d/99_hcocena_workspace`) and at
the first R session (`Rprofile.site`), so the workspace also appears when a
platform mounts its own home directory. Each part (hcocena, deseq2) is set up
once per home (listed in `~/.hcocena_workspace`), so a part added in a newer
image appears once, and existing files are never overwritten.

RStudio Server is exposed on port `8787` by the base image.

## Updating the Docker Hub description

When publishing a new public Docker tag, update the overview file from the
repository root with:

```bash
Rscript docker/update_dockerhub_overview.R <new-tag>
```

This updates [`DOCKERHUB_OVERVIEW.md`](DOCKERHUB_OVERVIEW.md) so the new tag
becomes the recommended pinned version and moves the previous recommended tag
into the older reproducibility list. The `docker pull` / `docker run` quick
start examples intentionally stay on `latest`, so the Docker Hub page still
points to the current image even if the pinned-version text is not copied right
away.

If Docker Hub is not syncing this file automatically, copy the updated markdown
to the repository Overview field there right after pushing the new tag.
