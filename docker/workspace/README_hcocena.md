# hCoCena workspace

```text
~/
├── README_hcocena.md        this file
├── reference_files/         shared references (pathways, GO, Hallmark, TFs, immune cells)
└── projects/
    └── example_project/
        ├── 01_hcocena_main.Rmd       core workflow
        ├── 02_hcocena_satellite.Rmd  optional analyses
        ├── count_data/               expression files
        ├── annotation_data/          sample annotation files
        ├── project_data/             other inputs, e.g. .rds objects
        └── output/                   results
```

Each project keeps its data and results in its own folder. The notebooks use
paths relative to the project folder, so a copied or renamed project keeps
working. All projects share `~/reference_files/`.

## Start a new project

Copy `example_project/` (for example into a mounted project folder) and rename
it. The notebooks find their data next to themselves; the reference files come
from `~/reference_files/` (or from the image's copy in
`/usr/local/share/hcocena/reference_files/`), wherever the project lies.

## Help

- package vignette: `vignette("hcocena-s4-workflow")`
- function help, e.g. `?hcocena::hc_define_layers`
- unmodified notebook templates:
  `system.file("scripts", "workflows", package = "hcocena")`
- source code and issues: https://github.com/BioCompNet/hcocena
