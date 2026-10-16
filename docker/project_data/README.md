Project-specific TB-neo input files are intentionally not bundled into the
public Docker image.

Before running `01_hcocena_main.Rmd` of this project, place these files into
this `project_data/` folder:

- `sample_table_230530.rds`
- `vst_anno_log_230530.rds`
- `exclusion_vector_lm_full_240625.rds`

The notebook reads them from here (`project_data/` next to the notebook).
