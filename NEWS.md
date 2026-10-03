# hcocena 0.99.12

## One simple function for LLM-based module interpretation

- `hc_llm_enrichment()` replaces `hc_module_function_llm()`,
  `hc_module_function_gemini()` and `hc_module_function_vllm()`. It needs only
  the feature sets (modules of `hc` or any named list of genes/proteins), an
  optional system prompt, an optional context and the provider:
  `"claude"`, `"openai"`, `"gemini"`, `"local"` (vLLM or any
  OpenAI-compatible server) or `"manual"`. API keys are read from the usual
  environment variables.
- `provider = "manual"` works without an API key: the complete request is
  copied to the clipboard (or written to a file where no clipboard exists,
  e.g. in Docker), pasted into the web chat of any provider, and the answer is
  read back with the new `hc_llm_import()`. The answer is checked against the
  request id; code fences, surrounding text and typographic quotes are
  tolerated. Long requests are split into parts, and the next part is copied
  automatically after each import. Imported results are stored exactly like
  API results.
- The default system prompt is neutral (genes, proteins, metabolites) instead
  of being written for monocyte transcriptomics; `system_prompt =` replaces it.
- Each set is interpreted in separate requests so that every field sees only
  its own input: `general_processes` from the features alone (no context),
  `contextual_state` and `key_regulators` with the context, and the new
  `grounded_processes` (+ `supporting_evidence`) with the context and the
  significant terms of `hc_functional_enrichment()` (`evidence =`; optionally
  also `"upstream"` regulators and `"rag"` literature passages from a DoRAG
  server). `hc_plot_llm_enrichment()` draws all four plots; without stored
  enrichment the evidence request is left out. In the web-chat route this
  means one prompt per request type, each pasted into a new chat.
- `hc_plot_module_function_llm()` is now `hc_plot_llm_enrichment()` and reads
  the slot `"llm_enrichment"`. The alias `hc_plot_module_function_gemini()` is
  gone.
- Removed with the old function: DoRAG literature retrieval and the
  enrichment-grounded interpretation level.
- `clipr` (Suggests) is used for clipboard access.
- The prompt asks for `general_processes` as short phrases of at most about
  90 characters, without naming or listing individual genes or proteins.
- `hc_plot_llm_enrichment()` leaves gene/protein lists in parentheses out of
  the plot and wraps long terms onto at most two lines, never inside
  parentheses. Text that still does not fit is cut after a whole process (at a
  `" / "`), not mid-phrase. The heatmap cells stay square.
- Models that reject the `temperature` parameter (e.g. Claude 5) are retried
  once without it.
- Greek letters and typographic dashes are spelled out in plot text, which
  the PNG devices on Windows could not draw.

## Upstream inference rebuilt

- `hc_upstream_inference()` now links modules and regulators by
  over-representation: a hypergeometric test of each regulator's targets
  among the module genes, against all network genes, with BH correction over
  all module x regulator tests (`method = "ora"`). Before, decoupleR ULM was
  run on the genes of one module only, so targets were compared with genes
  that move the same way; the smallest p-value over all conditions was used
  without correction, and scores were averaged over GFC conditions, where they
  cancel out.
- Per condition, the result reports the signed score of the module's targets
  (mode of regulation x GFC/FC), the share of targets that agree
  (`consistency`) and the regulator activity from `decoupleR::run_ulm()` over
  all network genes (`activity`, `activity_qvalue`). The per-condition
  heatmaps mark regulators that are linked to the module and active in the
  condition.
- TFs that are themselves genes of the module are reported
  (`regulator_in_module`, `regulator_cor`) and framed in the heatmaps.
- Regulators whose module targets repeat a better-ranked one are flagged in
  `redundant_with` and left out of the plots (`collapse_redundant`,
  `redundancy_jaccard`). Links need at least `min_overlap = 3` targets.
- TF regulons come from CollecTRI when `OmnipathR` can load it and from
  DoRothEA otherwise (`tf_resource`); the resource used is stored. PROGENy uses
  the top `progeny_top = 100` genes per pathway; the fallback used before took
  the whole model (all genes) as pathway targets.
- `method = "ulm"` and `activity_input = "expression"` were removed; FC values
  are no longer clipped at `range_GFC`.

## Module knowledge network redrawn

- `hc_plot_enrichment_upstream_network()` now reads causally from left to
  right: module heatmap | upstream regulators -> modules -> enriched terms,
  drawn as one figure with square heatmap cells and module boxes.
- Arrows show how a regulator acts on its targets in the module (new column
  `regulation` of `hc_upstream_inference()`: activating, repressing, mixed),
  instead of the direction of change in the strongest condition.
- Regulators that are themselves genes of a module are filled with that
  module's colour; TFs and PROGENy pathways have different symbols.
- Enrichment terms whose module genes lie mostly in a better term are left
  out (`collapse_redundant_terms`). On the per-module pages, dashed lines
  link a regulator to the terms over-represented among its targets
  (`link_min_share`).

# hcocena 0.99.11

- `citation("hcocena")` also lists the STAR Protocols paper (Holsten et al.
  2024) next to the method paper (Oestreich et al. 2022).
- The vignette ends with a "Citing hcocena" section.
- The README maps the function names used in both papers (e.g.
  `functional_enrichment()`, `hcobject`) to their `hc_` counterparts.

# hcocena 0.99.10

- Condition messages are built from their parts instead of with `paste()`,
  and three closures no longer assign to their enclosing function with
  `<<-`; behaviour is unchanged.
- `KEGG` added to the biocViews.

# hcocena 0.99.9

- No code changes. Version bump to re-run the Bioconductor build: the
  previous BiocCheck run stopped before checking the package, because
  downloading a resource from bioconductor.org timed out (HTTP 504).

# hcocena 0.99.8

## New accessors and example data

- `hc_graph()` returns the integrated network and `hc_satellite()` the stored
  downstream results (all of them, or one entry by name), so that examples,
  the vignette and user code no longer need to reach into slots with `@`.
- `hc_example_data()` loads the example objects shipped in `inst/extdata` at
  a given workflow stage and points their output directory at a fresh
  directory inside `tempdir()`. The stored objects carried the temporary
  directory of the session that built them, so examples wrote their files
  outside the current session's temporary directory (and, on other systems,
  into a relative path below the working directory).
- `inst/extdata/toy_celltype_markers.gmt` holds the four marker panels of the
  example data, for offline examples of the enrichment, cell-type and
  upstream functions.

## Removed

- `hc_get_reference_data()` is gone, and with it the `ExperimentHub`
  suggestion. It was a placeholder for a planned `hcocenaData` package that
  does not exist, so there was nothing it could fetch.

## Bug fixes

- Names given to `custom_gmt_files` / `custom_pathway_gmt`
  (e.g. `c(CellTypes = "markers.gmt")`) are used as database labels, as
  documented. They were dropped, so every custom file was labelled
  `CustomEnrichment1`, `CustomCellType1`, ... regardless.
- `hc_celltype_annotation()` falls back to the `padj`-adjusted p-value where
  clusterProfiler's Storey q-value is `NA` (pi0 cannot be estimated, e.g. when
  a module is tested against a single term), as `hc_functional_enrichment()`
  already did. Such hits were silently dropped, down to "No significant
  cell-type terms found".
- `hc_upstream_inference()` no longer fails with "factor level is duplicated"
  when several layers share group names. Those GFC columns are labelled with
  their layer (as in the module heatmap) and stay separate conditions.
- `hc_check_tf()` works through the S4 API: the targets found by
  `hc_tf_overrep_network()` are now kept in the object
  (`hc_satellite(hc, "tf_network_targets")`); they were lost when the result
  was written back, so `hc_check_tf()` could never find them. It now also says
  when `hc_tf_overrep_network()` has to be run first, or which TFs are
  available.

## Documentation

- Every exported function has an example. Those that need a web service
  (ChEA3, Cytoscape, LLM providers) are wrapped in `\donttest{}` and guarded
  so that they do nothing when the service is unavailable.
- The vignette has an installation section and uses accessors throughout.

# hcocena 0.99.7

## Reliability of the analysis state

- Results are discarded when the inputs or parameters they were computed from
  change. `hc_read_data()`, `hc_set_layer_settings()` and `hc_set_cutoff()`
  drop everything downstream along data -> correlation -> layer network/GFC ->
  integration -> modules -> downstream analyses, but only when something
  actually changed, so a repeated call with identical arguments keeps the
  analysis. Previously a finished object survived having all its counts
  replaced, and setting a new cutoff left the object reporting a cutoff its
  network had not been built with.
- Downstream analyses report a missing module assignment instead of failing
  inside `apply()`.
- `hc_set_cutoff()` resolves named vectors by name only, against both the
  layer names set in `hc_define_layers()` and the internal ids. Names were
  matched against the internal ids alone, so naming a layer as the user sees it
  fell through to positional assignment - with the cutoffs swapped whenever the
  named order differed. A named value no longer spills onto layers it did not
  name. Unknown or duplicate names, mixed named/unnamed input, values outside
  `[-1, 1]` and a length that is neither 1 nor the layer count are now errors.
- Imported correlation and p-value matrices are aligned by gene name instead of
  by position, so files sorted differently no longer attach p-values to the
  wrong gene pairs, and are checked for squareness, numeric content, unique
  gene names and plausible value ranges. A single layer is handled like several
  layers: `NA` means "do not import" everywhere.
- Enrichment results name each statistic for what it is: `p_adjusted` (the
  requested correction, and what the filter uses), `padj_method`, and
  `q_storey` (clusterProfiler's Storey q-value, `NA` where pi0 could not be
  fitted). `qvalue` remains as a documented alias of `p_adjusted`.
- S4 calls no longer overwrite a variable called `hcobject` in the user's
  workspace. The bridge mirrored its whole legacy object into `.GlobalEnv`
  whenever it found one there, a leftover from the removed legacy API.
- Validation rejects a non-`igraph` object in `integration@graph`, cutoffs
  outside `[-1, 1]`, and duplicate layer ids.
- CI checks Windows as well as Linux, fails on warnings, and runs BiocCheck.
  The README badge pointed at a workflow file that did not exist.

A full analysis captured before and after these changes - cutoffs, edge list
and weights, module membership, GFCs, enrichment, module significance, GFC
tables, hub genes - is identical in every field.

## Methodological options

- `hc_functional_enrichment(universe = )` chooses the background of the
  hypergeometric test. `"all_genes"` (default, unchanged) tests against every
  measured gene; `"network"` tests against the genes that entered the
  integrated network - the genes that could have landed in a module at all -
  which removes the bias introduced by the top-variance selection and the
  correlation cutoff.
- `hc_module_condition_significance(standardize_modules = TRUE)` z-scores each
  module across the samples of a layer before testing, so modules are compared
  on a common scale rather than on their absolute expression. Because it is
  applied per layer it also removes the between-layer offset when layers are
  pooled - neither the limma design nor the LMM carries a layer term. Default
  `FALSE`.
- `hc_module_function_llm(use_enrichment = TRUE)` adds a further module
  interpretation grounded in the module's significant enrichment terms,
  alongside the existing gene-only and RAG interpretations. The enrichment
  terms are a test on exactly the gene list being interpreted, whereas RAG
  passages are retrieved by similarity, so the prompt states that ranking
  explicitly. Results are stored in `enrichment_response` and are available to
  `hc_plot_module_function_llm()` as `enrichment_general_processes`,
  `enrichment_contextual_state` and `enrichment_key_regulators`. Terms can also
  be supplied directly via `enrichment_terms`, which is the only way to combine
  this with `genes =`.

## Examples and fixtures

- Rebuilt the bundled example fixtures with a 16-donor, two-timepoint
  longitudinal design, so the meta-clustering functions can be demonstrated on
  real output rather than only described. `inst/scripts/make-fixtures.R`
  documents how they are generated.
- Dropped two stored-but-never-read plot objects (`dd_plot_calculated_optimal`
  and the layer heatmap) from the fixtures. Together with the new design this
  takes the source tarball from 2.8 MB to 1.2 MB.
- Added runnable examples to 17 further help pages, covering the longitudinal
  meta-clustering chain, the direct workflow, `hc_meta_correlation_num()`, the
  cell-type database helpers, and the `hc_sample_regrouping` plot method.

## Fixes

- The Cytoscape helpers now say what is wrong. With Cytoscape closed, RCy3
  receives an empty CyREST response and dies inside its own parsing with
  "$ operator is invalid for atomic vectors", which gave the user nothing to
  act on. `hc_export_to_cytoscape()` and `hc_import_layout_from_cytoscape()`
  ping Cytoscape first and explain how to make it reachable, and
  `hc_import_layout_from_local_folder()` names the file it expected instead of
  failing with a bare connection error.
- `hc_longitudinal_workflow_direct()` failed immediately with `'arg' must be of
  length 1`. It forwards `cap_na_impute` explicitly, but the receiving formal
  defaults to the already-matched `na_impute`, so `match.arg()` saw a single
  choice. The three affected call sites now pass `choices` explicitly.
- `hc_meta_correlation_cat()` correlates across groups of the variable of
  interest, so it needs at least three of them. With two it fell through to a
  cryptic `cor.test()` error; it now reports the requirement and the groups it
  found.
- `hc_upstream_inference()`, `hc_plot_enrichment_upstream_network()` and
  `hc_celltype_activity_decoupler()` silently used only the first layer.
  `GFC_all_layers` repeats its condition columns once per layer, and these
  functions selected the value columns with `setdiff()` on the column *names*,
  which collapses the duplicates - so with two layers sharing their group
  names, half the matrix was dropped without any message. They now select by
  position, as `.hc_gfc_value_col_idx()` already did. Upstream-regulator and
  cell-type activity results computed on multi-layer objects whose layers share
  group names should be recomputed.

## Export reliability

- Build and validate XLSX workbooks on R's local temporary filesystem before
  publishing them to synchronized or bind-mounted output directories.
- Verify staged XLSX transfers byte-for-byte and atomically replace existing
  outputs without exposing partially written workbooks.
- Route both table-based exports and updated `Hub_genes.xlsx` workbooks through
  the same local-staging path.

# hcocena 0.99.6

## Enrichment defaults

- Made module-heatmap column gaps opt-in via `smart_column_gaps`, with
  `column_gap_by` for explicit metadata-based splits and `column_gap_mm` for
  gap size control.
- Switched functional-enrichment defaults to consistent term selection across
  modules and wrappers.
- Added optional DoRAG retrieval support to `hc_module_function_llm()` so LLM
  module summaries can be grounded in retrieved passages and stored citations.
- RAG runs now preserve the normal context-aware interpretation in `response`
  and store a separate literature-supported interpretation in `rag_response`.
- Added `rag_connect_timeout_sec` and `rag_continue_on_error` to make DoRAG
  retrieval robust to unreachable or slow RAG servers.
- Extended `hc_plot_module_function_llm()` to plot separate RAG fields such as
  `rag_contextual_state` (also available as `contextual_state_rag`).
- Added common RNA-seq differential-expression packages to the Docker image,
  including `DESeq2`, `limma`, `sva`, `edgeR`, and supporting visualization,
  shrinkage, import, and organism annotation packages.
- Allowed `hc_split_modules()` to use one Leiden `resolution` value per
  selected module.
- Clarified that `hc_read_data()` now removes zero-variance genes and drops
  non-numeric helper columns from object-based count inputs, which can shift
  `hc_suggest_topvar()` inflection points slightly compared with older
  releases.
- Refreshed the Docker release metadata for the next public image tag.
- Rendered interactive cutoff plots inline during HTML/R Markdown knitting
  instead of opening them in the RStudio Viewer.

# hcocena 0.99.5

## Bioconductor readiness and Docker release

- Finalized the S4/legacy bridge cleanup, including removal of remaining
  package-level `<<-` usage from the active R sources.
- Hardened regression coverage for the updated heatmap and auto-tuning paths.
- Refined package formatting and documentation metadata ahead of submission.
- Refreshed the public Docker release metadata for the next image tag.

# hcocena 0.99.1

## Bioconductor preparation

- Aligned the package version with Bioconductor pre-submission conventions.
- Simplified `DESCRIPTION` metadata for the first Bioconductor submission.
- Added a package-level `README.md` and `inst/CITATION`.
- Expanded the workflow and migration vignettes to use `BiocStyle` and
  reproducible examples based on `inst/extdata`.
- Added ignore rules for local build artifacts and reduced the default branch to
  package source for submission.
- Added compatibility fixes for longitudinal `rfcont` imputation and clarified
  that this workflow requires `library(CALIBERrfimpute)` in the active session.
- Made LLM-related examples safe for package checks and non-interactive builds.
