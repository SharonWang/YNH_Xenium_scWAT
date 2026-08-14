# scWAT Xenium Extended QC Diagnostics Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Extend the existing section and slide-summary R notebooks with direct cycle-alarm evidence, candidate affected-gene diagnostics, subset-versus-full ranking, spatial QC clustering/artifact review, and within-mouse technical concordance.

**Architecture:** Reusable calculations remain in `R/source.R`; notebook cells only validate inputs, orchestrate functions, display bounded results, and write/reload artifacts. Section notebooks produce section-local alarm, gene-quality, and spatial summaries. The slide notebook combines exactly four validated section bundles and answers the four scientific questions without changing the existing Phase 0-2 readiness gates.

**Tech Stack:** R 4.3+, Matrix, jsonlite, ggplot2, arrow for HPC transcript-Parquet aggregation, RANN for scalable HPC nearest-neighbour calculations, Jupyter IRkernel, PowerShell local smoke-test launcher, Bash/Slurm HPC launcher.

## Global Constraints

- scWAT only; colon is outside this task.
- Local project and temporary files must remain below `D:/Xiaonan/CODEX_projects/Yanan_Xenium`; no project or temporary file may be written to `C:`.
- HPC root is `/dssg/home/acct-svetoslav_chakarov/svetoslav_chakarov/Lab_members/Yanan_Hu`.
- Local execution uses only the deterministic subset below `adipose_analysis/subset_input`.
- Full transcript diagnostics require projected/lazy Arrow aggregation; never collect the complete transcript table.
- Exact poor-cycle identity remains `CYCLE_IDENTITY_UNRESOLVED_REQUIRES_10X`.
- Candidate genes remain `CANDIDATE_NOT_CONFIRMED`.
- Coordinate statistics must not label hotspots as folds or tears; morphology review is required.
- Mouse is the biological replicate; section is a technical processing unit.
- No new diagnostic may delete cells or overwrite raw Xenium inputs.
- Existing Phase 0-2 readiness gates remain authoritative and unchanged by advisory diagnostics.
- Use fixed seed `20260814` for permutation diagnostics.

## File Structure

- Create `config/extended_qc_defaults.tsv`: versioned thresholds and computational settings.
- Create `config/subset_qc_reference.tsv`: fixed validated subset review-rate reference.
- Create `tests/test_extended_qc.R`: focused pure-R tests for every new diagnostic behavior.
- Modify `R/source.R`: reusable configuration, alarm, gene-quality, spatial, ranking, concordance, plotting, and artifact functions.
- Modify `scripts/render_notebooks.py`: generate the new notebook chunks and parameters.
- Regenerate `notebooks/01_section_phase0_2_QC.ipynb` and `notebooks/02_slide_QC_summary.ipynb` from the renderer.
- Modify `scripts/execute_local_subset.ps1`: run new tests and verify explicit local skip states/artifacts.
- Modify `shell/run_notebook_qc_hpc.sh` and `slurm/scwat_notebook_qc.sbatch`: full-data dependency/input preflight and diagnostic mode.
- Modify `README.md`: exact inputs, outputs, local mode, HPC chunks, and interpretation limits.
- Modify `reports/2026-08-14_scwat_qc_summary/artifact.json` and rebuild `report.html`: distinguish tested local calculations from HPC-only pending diagnostics.
- Modify `D:/Xiaonan/CODEX_projects/Yanan_Xenium/PROJECT_RESEARCH_PLAN_AND_PROGRESS.md`: checkpoint and result history.

---

### Task 1: Versioned diagnostic configuration and execution-mode contract

**Files:**
- Create: `config/extended_qc_defaults.tsv`
- Create: `config/subset_qc_reference.tsv`
- Create: `tests/test_extended_qc.R`
- Modify: `R/source.R`

**Interfaces:**
- Consumes: project-root-confined config paths, requested mode `AUTO|LOCAL_SUBSET|FULL_HPC`, section input directory.
- Produces: `read_extended_qc_config(path) -> named list`; `resolve_extended_qc_mode(requested_mode, region_dir) -> character(1)`; `extended_qc_preflight(mode, region_dir, config) -> data.frame`.

- [ ] **Step 1: Write failing configuration and mode tests**

Add literal fixtures to `tests/test_extended_qc.R`:

```r
source(file.path(repo_root, "R", "source.R"))
cfg <- read_extended_qc_config(file.path(repo_root, "config", "extended_qc_defaults.tsv"))
stopifnot(cfg$qv_threshold == 20, cfg$spatial_k == 15L, cfg$permutations == 999L)
stopifnot(resolve_extended_qc_mode("AUTO", subset_region_dir) == "LOCAL_SUBSET")
stopifnot(resolve_extended_qc_mode("AUTO", full_region_dir) == "FULL_HPC")
expect_error(resolve_extended_qc_mode("BAD", subset_region_dir), "AUTO")
preflight <- extended_qc_preflight("LOCAL_SUBSET", subset_region_dir, cfg)
stopifnot(preflight$status[preflight$check == "transcripts_parquet"] == "SKIP_ALLOWED")
```

The test fixture creates empty `transcripts.parquet` only for `full_region_dir`; it does not write outside `tempdir()`, which the launcher forces below D:.

- [ ] **Step 2: Run the new test and observe the expected failure**

Run:

```powershell
$env:TMPDIR='D:\Xiaonan\CODEX_projects\Yanan_Xenium\adipose_analysis\tmp'
$env:TEMP=$env:TMPDIR; $env:TMP=$env:TMPDIR
& 'C:\Program Files\R\R-4.3.3\bin\Rscript.exe' 'D:\Xiaonan\CODEX_projects\Yanan_Xenium\adipose_analysis\YNH_Xenium_scWAT\tests\test_extended_qc.R'
```

Expected: FAIL because `read_extended_qc_config()` is not defined.

- [ ] **Step 3: Add the two configuration tables**

`extended_qc_defaults.tsv` contains typed values for:

```text
qv_threshold=20; spatial_k=15; grid_size_um=100; permutations=999;
dense_quantile=0.90; min_bin_cells=20; hotspot_fdr=0.05;
candidate_abs_log2_ratio=0.50; candidate_abs_q20_delta=0.05;
concordance_gene_spearman=0.90; concordance_review_rate_difference=0.05;
concordance_count_ratio_lower=0.67; concordance_count_ratio_upper=1.50;
concordance_feature_ratio_lower=0.80; concordance_feature_ratio_upper=1.25;
seed=20260814
```

`subset_qc_reference.tsv` contains literal validated values:

```text
Region_1  500  37  0.074  1
Region_2  500  16  0.032  3
Region_3  500  27  0.054  2
Region_4  500  11  0.022  4
```

- [ ] **Step 4: Implement minimal config/mode/preflight functions**

Implementation requirements:

```r
read_extended_qc_config <- function(path) {
  x <- utils::read.delim(path, check.names = FALSE)
  # require key/value/type; reject duplicates, blanks, unknown types
  # coerce integer, numeric, logical, or character and return named list
}

resolve_extended_qc_mode <- function(requested_mode = "AUTO", region_dir) {
  requested_mode <- toupper(requested_mode)
  if (!requested_mode %in% c("AUTO", "LOCAL_SUBSET", "FULL_HPC")) stop(...)
  if (requested_mode != "AUTO") return(requested_mode)
  if (file.exists(file.path(region_dir, "transcripts.parquet"))) "FULL_HPC" else "LOCAL_SUBSET"
}
```

`extended_qc_preflight()` returns one row per input/package with `check`, `required`, `available`, `status`, and `details`. In full mode, missing `transcripts.parquet`, `arrow`, or `RANN` is `FAIL`; in subset mode it is `SKIP_ALLOWED`.

- [ ] **Step 5: Run focused and regression tests**

Expected: both commands exit 0.

```powershell
& $Rscript "$RepoRoot\tests\test_extended_qc.R"
& $Rscript "$RepoRoot\tests\test_source.R"
```

- [ ] **Step 6: Commit Task 1**

```bash
git add config/extended_qc_defaults.tsv config/subset_qc_reference.tsv tests/test_extended_qc.R R/source.R
git commit -m "feat: add extended QC configuration contract"
```

### Task 2: Direct cycle-alarm evidence and bounded per-gene quality summaries

**Files:**
- Modify: `tests/test_extended_qc.R`
- Modify: `R/source.R`

**Interfaces:**
- Consumes: existing alarm table; sparse gene-expression matrix/features; bounded transcript-quality aggregate; signature table.
- Produces: `build_cycle_alarm_evidence(alarms, region_id) -> data.frame`; `summarise_gene_matrix_qc(counts, features, region_id, gene_sets) -> data.frame`; `summarise_transcript_quality_table(transcripts, region_id, qv_threshold) -> data.frame`; `summarise_transcript_quality_arrow(path, region_id, qv_threshold) -> data.frame`; `combine_gene_quality(matrix_qc, transcript_qc) -> data.frame`.

- [ ] **Step 1: Write failing tests with hand-calculated expected values**

```r
evidence <- build_cycle_alarm_evidence(alarm_fixture, "Region_1")
stopifnot(evidence$evidence_status == "DIRECT_EVIDENCE")
stopifnot(evidence$cycle_identity_status == "CYCLE_IDENTITY_UNRESOLVED_REQUIRES_10X")

tx <- data.frame(feature_name=c("A","A","B"), qv=c(30,10,25), codeword_index=c(1,1,2))
txq <- summarise_transcript_quality_table(tx, "Region_1", 20)
stopifnot(txq$transcript_rows[txq$gene=="A"] == 2L)
stopifnot(txq$mean_qv[txq$gene=="A"] == 20)
stopifnot(txq$fraction_q20[txq$gene=="A"] == 0.5)
```

Add a no-alarm fixture and ambiguous-schema fixture. The schema test must fail with a message naming missing gene/QV fields.

- [ ] **Step 2: Run and observe failure because functions are absent**

Run `Rscript tests/test_extended_qc.R` with D:-local temp. Expected: FAIL at `build_cycle_alarm_evidence`.

- [ ] **Step 3: Implement pure-R summaries and schema resolver**

The schema resolver accepts `feature_name|gene|target_name`, `qv|quality_value`, and optional `codeword_index`. It rejects multiple matching aliases. `build_cycle_alarm_evidence()` always returns one bounded row, including a no-alarm row.

- [ ] **Step 4: Implement the Arrow adapter without full-table collection**

```r
summarise_transcript_quality_arrow <- function(path, region_id, qv_threshold = 20) {
  require_package("arrow")
  dataset <- arrow::open_dataset(path, format = "parquet")
  fields <- resolve_transcript_schema(names(dataset$schema))
  # select only gene, QV, and optional codeword; group and aggregate in Arrow
  # collect only the per-gene result
}
```

Record `reader="arrow_projected_aggregate"` and source row count. Never expose a code path that calls `collect(dataset)` before grouping.

- [ ] **Step 5: Run tests**

Pure-R tests must pass locally. Arrow integration is conditionally tested only when `requireNamespace("arrow", quietly=TRUE)`; otherwise the test asserts the full-mode preflight failure and subset skip record.

- [ ] **Step 6: Commit Task 2**

```bash
git add tests/test_extended_qc.R R/source.R
git commit -m "feat: add Xenium alarm and gene quality diagnostics"
```

### Task 3: Scalable spatial edge, density, clustering, and hotspot diagnostics

**Files:**
- Modify: `tests/test_extended_qc.R`
- Modify: `R/source.R`

**Interfaces:**
- Consumes: cell table with `region_id`, `cell_id`, `x_centroid`, `y_centroid`, `qc_review_flag`; config; seed.
- Produces: `calculate_knn_density(cells, k, mode) -> numeric`; `assign_spatial_grid(cells, grid_size_um) -> data.frame`; `summarise_spatial_enrichment(annotated_cells) -> data.frame`; `test_spatial_flag_clustering(annotated_cells, k, permutations, seed) -> data.frame`; `find_spatial_qc_hotspots(annotated_cells, permutations, min_bin_cells, fdr, seed) -> data.frame`; `plot_extended_spatial_qc(...) -> named list`.

- [ ] **Step 1: Write failing deterministic spatial tests**

Create a literal 6x6 grid with review flags concentrated in one corner. Assert:

```r
ann <- assign_spatial_grid(spatial_fixture, grid_size_um=10)
stopifnot(sum(ann$edge_proxy) == 20L)
enrich <- summarise_spatial_enrichment(ann)
stopifnot(enrich$flagged[enrich$class=="edge"] == 4L)
hot <- find_spatial_qc_hotspots(ann, permutations=199L, min_bin_cells=4L, fdr=0.10, seed=20260814L)
stopifnot(any(hot$hotspot_status == "MORPHOLOGY_REVIEW_REQUIRED"))
```

Add a randomized/no-flag fixture that returns `NOT_ESTIMABLE` rather than dividing by zero.

- [ ] **Step 2: Run and observe the expected undefined-function failure**

Run `Rscript tests/test_extended_qc.R`. Expected: FAIL at `assign_spatial_grid`.

- [ ] **Step 3: Implement deterministic grid/edge and enrichment functions**

Grid indices are `floor((coordinate-min_coordinate)/grid_size_um)`. An occupied bin is an edge proxy when any of its eight neighboring bins is unoccupied. Risk ratios use a 0.5 continuity correction and always retain numerator/denominator counts.

- [ ] **Step 4: Implement local-density and clustering functions**

For `n <= 2000` in local subset mode, use a chunked base-R distance fallback. For full mode, require `RANN::nn2()` and reject a base-R full distance matrix. Use a symmetric kNN adjacency for the global binary-flag statistic. Permute the flag labels within section using the fixed seed.

- [ ] **Step 5: Implement grid hotspot permutation and plots**

For each eligible bin, compare observed flagged count/rate with within-section permuted labels, calculate empirical p-values `(1 + exceedances)/(1 + permutations)`, adjust with BH, and require both adjusted p-value and positive rate difference. Output bounding boxes in microns and `MORPHOLOGY_REVIEW_REQUIRED`.

- [ ] **Step 6: Run focused and regression tests**

Expected: deterministic repeated runs are identical; all existing tests remain green.

- [ ] **Step 7: Commit Task 3**

```bash
git add tests/test_extended_qc.R R/source.R
git commit -m "feat: add spatial Xenium QC diagnostics"
```

### Task 4: Candidate affected-gene evidence tiers and subset/full ranking

**Files:**
- Modify: `tests/test_extended_qc.R`
- Modify: `R/source.R`

**Interfaces:**
- Consumes: four-section gene-quality tables; verified manifest; extended config; subset reference; full QC summary.
- Produces: `rank_candidate_cycle_genes(gene_quality, manifest, config) -> data.frame`; `compare_subset_full_qc(full_summary, subset_reference) -> list(ranking, agreement)`; candidate/ranking plots.

- [ ] **Step 1: Write failing candidate-tier tests**

Use literal gene metrics where:

- gene A is depleted and has lower Q20 in Region 4 versus Region 3, expecting Tier A;
- gene B has concordant evidence in Regions 1 and 2, expecting Tier B;
- gene C has depletion only, expecting Tier C;
- gene D lacks transcript metrics, expecting Unranked/insufficient.

Assert every row includes `candidate_status == "CANDIDATE_NOT_CONFIRMED"` and `exact_cycle_status == "REQUIRES_10X_DIAGNOSTICS"`.

- [ ] **Step 2: Write failing subset/full ranking test**

```r
full <- data.frame(region_id=paste0("Region_",1:4), input_cells=1000L,
                   review_flagged=c(80L,20L,60L,40L))
cmp <- compare_subset_full_qc(full, subset_reference_fixture)
stopifnot(identical(cmp$ranking$full_rank, c(1L,4L,2L,3L)))
stopifnot(is.numeric(cmp$agreement$spearman_rho))
stopifnot(cmp$agreement$interpretation == "DESCRIPTIVE_FOUR_SECTIONS")
```

- [ ] **Step 3: Run and observe failure at missing functions**

- [ ] **Step 4: Implement candidate ranking with explicit comparison types**

Calculate per-gene log2 CP10k ratio, detection-fraction difference, and Q20 difference against Region 3. Mark Region 4 versus Region 3 as `WITHIN_MOUSE_TECHNICAL_PAIR`; mark Regions 1/2 versus Region 3 as `CROSS_MOUSE_DIAGNOSTIC_REFERENCE`. Preserve raw metrics and threshold columns in every output row.

- [ ] **Step 5: Implement subset/full comparison and plots**

Use deterministic descending rank with `region_id` as tie-breaker. Calculate Spearman and Kendall coefficients without presenting p-values as biological evidence.

- [ ] **Step 6: Run all tests and commit Task 4**

```bash
git add tests/test_extended_qc.R R/source.R
git commit -m "feat: rank cycle candidates and full-data QC burden"
```

### Task 5: Within-mouse technical concordance

**Files:**
- Modify: `tests/test_extended_qc.R`
- Modify: `R/source.R`

**Interfaces:**
- Consumes: verified manifest, section QC summary, combined cell metadata, combined gene-quality table, config.
- Produces: `build_section_pairs(manifest) -> data.frame`; `calculate_within_mouse_concordance(...) -> list(summary, genes)`; concordance plots.

- [ ] **Step 1: Write failing pair-construction tests**

Assert exact pairs:

```r
pairs <- build_section_pairs(manifest_fixture)
stopifnot(pairs$section_a == c("62308","62310"))
stopifnot(pairs$section_b == c("62309","62311"))
stopifnot(pairs$region_a == c("Region_1","Region_3"))
stopifnot(pairs$region_b == c("Region_2","Region_4"))
```

Reject mice with other than two sections using `NOT_ESTIMABLE` rows rather than silently choosing a pair.

- [ ] **Step 2: Write failing metric tests with literal ratios and flags**

Construct one concordant pair and one review pair. Assert the review-rate difference, median ratios, gene Spearman correlations, each advisory criterion, and final `CONCORDANT|REVIEW|NOT_ESTIMABLE` status.

- [ ] **Step 3: Run and observe failure at `build_section_pairs`**

- [ ] **Step 4: Implement bounded quantile distances and gene correlations**

Use quantiles `seq(0.01,0.99,0.01)` for transcript count, feature count, area, and control fraction. Report mean absolute quantile difference scaled by the pooled median. Use complete-case Spearman correlations for log1p CP10k and detection fractions; require at least 20 complete panel genes.

- [ ] **Step 5: Implement advisory criteria and plots**

Write one Boolean/result column per criterion. Overall `CONCORDANT` requires every estimable advisory criterion to pass; any failure is `REVIEW`; insufficient data is `NOT_ESTIMABLE`. State that thresholds are advisory and do not update readiness.

- [ ] **Step 6: Run all tests and commit Task 5**

```bash
git add tests/test_extended_qc.R R/source.R
git commit -m "feat: assess within-mouse section concordance"
```

### Task 6: Extended section artifact contract and section notebook chunks

**Files:**
- Modify: `tests/test_source.R`
- Modify: `tests/test_extended_qc.R`
- Modify: `R/source.R`
- Modify: `scripts/render_notebooks.py`
- Regenerate: `notebooks/01_section_phase0_2_QC.ipynb`

**Interfaces:**
- Consumes: Task 1-5 functions and current Phase 0-2 section object.
- Produces: extended section tables/figures, `extended_qc_status.tsv`, and reload validation.

- [ ] **Step 1: Write failing artifact-contract tests**

Update a temporary section fixture and assert `extended_section_required_artifacts(region_id, mode)` requires:

```text
extended_qc_preflight.tsv
cycle_alarm_evidence.tsv
gene_transcript_quality.tsv
spatial_qc_global.tsv
spatial_qc_edge_density.tsv
spatial_qc_hotspots.tsv
spatial_qc_cell_annotations.tsv.gz
spatial_manual_review_manifest.tsv
extended_qc_status.tsv
figures/<region>_extended_spatial_qc.pdf
```

In local subset mode, `gene_transcript_quality.tsv` contains valid matrix metrics plus transcript status `NOT_RUN_LOCAL_SUBSET`. In full mode, this skip state is rejected.

- [ ] **Step 2: Run tests and observe missing-contract failure**

- [ ] **Step 3: Implement writer/reloader/validator functions**

Add `write_extended_section_artifacts(...)`, `read_extended_section_artifacts(...)`, and `validate_extended_section_artifacts(...)`. Use existing path guards and atomic bounded tables. Preserve Phase 0-2 artifacts unchanged.

- [ ] **Step 4: Add section-notebook parameters and chunks**

Add parameters:

```r
EXTENDED_QC_MODE <- "AUTO"
EXTENDED_QC_CONFIG_PATH <- file.path(PIPELINE_REPO,"config","extended_qc_defaults.tsv")
```

Add chunks in order: extended preflight; direct alarm evidence; matrix/transcript gene quality; spatial diagnostics; extended figures; artifact write/reload. Each chunk prints its inputs, outputs, mode, and caveats.

- [ ] **Step 5: Regenerate and validate the section notebook**

```powershell
& $Python "$RepoRoot\scripts\render_notebooks.py" --build-section
& $Python "$RepoRoot\scripts\render_notebooks.py" --validate "$RepoRoot\notebooks\01_section_phase0_2_QC.ipynb" --type section
```

Expected: validation passes; notebook contains no `C:` path.

- [ ] **Step 6: Run tests and commit Task 6**

```bash
git add R/source.R tests/test_source.R tests/test_extended_qc.R scripts/render_notebooks.py notebooks/01_section_phase0_2_QC.ipynb
git commit -m "feat: extend section QC notebook diagnostics"
```

### Task 7: Slide aggregation, four-question outputs, and final QC notebook

**Files:**
- Modify: `tests/test_extended_qc.R`
- Modify: `R/source.R`
- Modify: `scripts/render_notebooks.py`
- Regenerate: `notebooks/02_slide_QC_summary.ipynb`

**Interfaces:**
- Consumes: exactly four validated extended section bundles, subset reference, verified manifest, config.
- Produces: combined gene/spatial tables, candidate ranking, subset/full comparison, concordance tables/figures, extended slide status.

- [ ] **Step 1: Write failing four-section aggregation tests**

Create four minimal extended section fixtures. Assert that missing, duplicate, or mixed-mode sections fail. Assert combined row counts and exact Region 1-4 order.

- [ ] **Step 2: Run and observe missing aggregation-function failure**

- [ ] **Step 3: Implement extended slide reader/writer**

Add `read_extended_slide_qc_outputs()`, `summarise_extended_slide_qc()`, `write_extended_slide_qc_artifacts()`, and `validate_extended_slide_qc_artifacts()`. Required combined outputs:

```text
combined_cycle_alarm_evidence.tsv
combined_gene_transcript_quality.tsv
candidate_cycle_affected_genes.tsv
subset_full_qc_ranking.tsv
subset_full_qc_rank_agreement.tsv
combined_spatial_qc.tsv
combined_spatial_hotspots.tsv
within_mouse_section_concordance.tsv
within_mouse_gene_concordance.tsv
extended_slide_qc_status.tsv
figures/scwat_extended_qc_diagnostics.pdf
```

- [ ] **Step 4: Add one notebook section per scientific question**

The summary notebook order is: coverage/reload; alarm evidence and candidate genes; subset/full ranking; spatial clustering/manual morphology review; within-mouse concordance; unchanged readiness; outputs/reload checks.

- [ ] **Step 5: Regenerate/validate the summary notebook and run tests**

Expected: both notebook validators and all R tests pass.

- [ ] **Step 6: Commit Task 7**

```bash
git add R/source.R tests/test_extended_qc.R scripts/render_notebooks.py notebooks/02_slide_QC_summary.ipynb
git commit -m "feat: summarize extended slide QC diagnostics"
```

### Task 8: Local subset execution with explicit HPC-only skip states

**Files:**
- Modify: `scripts/execute_local_subset.ps1`
- Modify: `README.md`
- Generated outside Git: `adipose_analysis/scwat_qc_outputs/local_extended_qc_test/`

**Interfaces:**
- Consumes: local subset and repository manifest/config.
- Produces: four executed section notebooks, slide summary, and validated local skip/status records.

- [ ] **Step 1: Add the extended test and mode parameters to the local runner**

Run `tests/test_extended_qc.R` after `tests/test_source.R`. Inject `EXTENDED_QC_MODE=LOCAL_SUBSET`. After each section, assert extended artifacts reload and transcript status is `NOT_RUN_LOCAL_SUBSET` when Parquet is absent.

- [ ] **Step 2: Execute all four local section notebooks and summary**

```powershell
& 'C:\Windows\System32\WindowsPowerShell\v1.0\powershell.exe' -NoProfile -ExecutionPolicy Bypass `
  -File "$RepoRoot\scripts\execute_local_subset.ps1" `
  -ProjectRoot 'D:\Xiaonan\CODEX_projects\Yanan_Xenium' `
  -RunLabel 'local_extended_qc_test'
```

- [ ] **Step 3: Reconcile unchanged core results**

Assert 2,000 cells, core-pass counts `493,500,497,499`, review counts `37,16,27,11`, metadata PASS in all sections, readiness `HOLD,HOLD,PASS,HOLD`, and zero deleted cells.

- [ ] **Step 4: Validate new local results**

Assert spatial tables exist for all sections, full-data ranking is explicitly unavailable in local mode where appropriate, candidate-cycle gene confirmation is never claimed, and cycle identity always requires 10x diagnostics.

- [ ] **Step 5: Update README local inputs/outputs and commit Task 8**

```bash
git add scripts/execute_local_subset.ps1 README.md
git commit -m "test: validate extended QC on local subsets"
```

### Task 9: HPC full-data preflight and copy/paste execution chunks

**Files:**
- Modify: `shell/run_notebook_qc_hpc.sh`
- Modify: `slurm/scwat_notebook_qc.sbatch`
- Modify: `README.md`

**Interfaces:**
- Consumes: full `${PROJECT_ROOT}/adipose_data`, verified manifest/config, R/Jupyter environment with arrow/RANN.
- Produces: full extended run and exact checkpoint commands.

- [ ] **Step 1: Add full-mode preflight**

Require `transcripts.parquet` for every Region 1-4 directory and R packages `IRkernel`, `Matrix`, `jsonlite`, `ggplot2`, `arrow`, `RANN`. Print projected input sizes and free space. Inject `EXTENDED_QC_MODE=FULL_HPC`.

- [ ] **Step 2: Add post-section checks**

After each executed notebook, call `validate_extended_section_artifacts(..., mode="FULL_HPC")` and reject any `NOT_RUN_LOCAL_SUBSET` transcript status.

- [ ] **Step 3: Add post-summary checks**

Reload the extended slide RDS/tables; assert four regions, two mouse pairs, nonempty gene-quality output, direct alarm pattern Regions 1/2/4, and `REQUIRES_10X_DIAGNOSTICS` cycle status.

- [ ] **Step 4: Document copy/paste HPC chunks**

README chunks must separately show: environment variables; package/input preflight; one-section interactive test; four-section full run; summary-only rerun; post-run validation; exact input/output paths.

- [ ] **Step 5: Validate shell syntax on HPC and commit**

Local Windows records `bash -n` as pending if Bash is unavailable. On HPC run:

```bash
bash -n shell/run_notebook_qc_hpc.sh
bash -n slurm/scwat_notebook_qc.sbatch
```

Commit:

```bash
git add shell/run_notebook_qc_hpc.sh slurm/scwat_notebook_qc.sbatch README.md
git commit -m "feat: run extended Xenium QC on HPC"
```

### Task 10: Report, progress record, and final verification

**Files:**
- Modify: `reports/2026-08-14_scwat_qc_summary/artifact.json`
- Regenerate: `reports/2026-08-14_scwat_qc_summary/report.html`
- Modify: `README.md`
- Modify outside repo: `D:/Xiaonan/CODEX_projects/Yanan_Xenium/PROJECT_RESEARCH_PLAN_AND_PROGRESS.md`

**Interfaces:**
- Consumes: verified local smoke-test results and current HPC availability state.
- Produces: updated portable report, documented pending HPC-only questions, clean repository.

- [ ] **Step 1: Update the progress record after every completed checkpoint**

Record test evidence, local skip states, output root, unchanged Phase 0-2 QC results, and remaining full-HPC/10x requirements.

- [ ] **Step 2: Update report source paths and findings**

The report distinguishes:

- directly observed alarm evidence;
- locally tested spatial/calculation functions;
- HPC-only transcript/candidate/full-distribution results not yet run;
- exact cycle identity requiring 10x diagnostics;
- advisory concordance thresholds from universal readiness gates.

- [ ] **Step 3: Rebuild and validate the portable report**

Use the packaged report builder with `TMPDIR/TMP/TEMP` set below D:. Expected: schema and package PASS; disclose structural-only browser verification if Chromium remains unavailable.

- [ ] **Step 4: Run the complete verification suite**

```powershell
& $Rscript "$RepoRoot\tests\test_source.R"
& $Rscript "$RepoRoot\tests\test_extended_qc.R"
& $Python "$RepoRoot\scripts\render_notebooks.py" --validate "$RepoRoot\notebooks\01_section_phase0_2_QC.ipynb" --type section
& $Python "$RepoRoot\scripts\render_notebooks.py" --validate "$RepoRoot\notebooks\02_slide_QC_summary.ipynb" --type summary
git -C $RepoRoot diff --check
git -C $RepoRoot status --short
```

Freshly reload every local section and slide artifact and reconcile exact counts/statuses before any completion claim.

- [ ] **Step 5: Commit Task 10**

```bash
git add reports/2026-08-14_scwat_qc_summary/artifact.json reports/2026-08-14_scwat_qc_summary/report.html README.md
git commit -m "docs: report extended Xenium QC diagnostics"
```

## Execution checkpoints

1. **Checkpoint A:** Tasks 1-2 — configuration, modes, direct alarms, gene-quality summaries.
2. **Checkpoint B:** Tasks 3-5 — spatial diagnostics, candidate/ranking logic, within-mouse concordance.
3. **Checkpoint C:** Tasks 6-7 — both notebooks and full artifact contracts.
4. **Checkpoint D:** Task 8 — complete local subset smoke test.
5. **Checkpoint E:** Tasks 9-10 — HPC scripts, report, progress record, and final verification.

Each checkpoint requires fresh passing tests, `git diff --check`, a progress-record update, and an intentional commit before continuing.
