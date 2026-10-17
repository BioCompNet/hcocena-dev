# Workspace

```text
~/
├── README_workspace.md        this file
├── reference_files/           hCoCena references (pathways, GO, Hallmark, TFs, immune cells)
└── projects/
    ├── hcocena_example/
    │   ├── 01_hcocena_main.Rmd       hCoCena core workflow
    │   ├── 02_hcocena_satellite.Rmd  optional analyses
    │   ├── count_data/               expression files
    │   ├── annotation_data/          sample annotation files
    │   ├── project_data/             other inputs, e.g. .rds objects
    │   └── output/                   results
    └── deseq2_example/
        ├── RNAseq_DESeq2_Workflow_V5_1.Rmd    DESeq2 workflow
        ├── RNAseq_DESeq2_Functions_V5_1.R     its functions (sourced by the workflow)
        └── RNAseq_DESeq2_V5_1_CheatSheet.pdf  overview of all steps
```

Each project keeps its data and results in its own folder, and the notebooks
use paths relative to that folder.

- **hCoCena:** the reference files come from `~/reference_files/` (or from the
  image's copy in `/usr/local/share/hcocena/reference_files/`), wherever the
  project lies.
- **DESeq2:** put the sample table, the counts (or kallisto / STAR output) and
  optionally the annotation and MultiQC table into the project folder and set
  the file names in `create_rnaseq_analysis()`. The results are written to the
  `output_dir` given there, inside the project folder.

## Start a new project

Copy an example project (for example into a mounted project folder) and rename
it. A copied DESeq2 project needs the functions file next to the workflow.

## Help

- hCoCena vignette: `vignette("hcocena-s4-workflow")`
- function help, e.g. `?hcocena::hc_define_layers`
- unmodified hCoCena notebook templates:
  `system.file("scripts", "workflows", package = "hcocena")`
- DESeq2 workflow: the cheat sheet in the project and `check_requirements()`
- hCoCena source code and issues: https://github.com/BioCompNet/hcocena
