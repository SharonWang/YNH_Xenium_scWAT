# YNH Xenium scWAT: Phase 0-2 QC notebooks

This repository contains the reusable scWAT Xenium QC code. Run one parameterized R notebook independently for each of the four sections, then run the slide-level summary notebook.

## Files

- `R/source.R`: general Xenium path, metadata, integrity, panel, sparse-import, QC, aggregation, artifact, and plotting functions.
- `notebooks/01_section_phase0_2_QC.ipynb`: reusable combined Phase 0-2 source template.
- `notebooks/01_section_phase0_2_QC_Region1.ipynb` through `01_section_phase0_2_QC_Region4.ipynb`: committed, parameter-locked notebooks to run and review one section at a time.
- `notebooks/02_slide_QC_summary.ipynb`: the sole slide-level QC summary/report, with final decision tables and Cell-inspired figures inline.
- `config/eos_gene_sets.tsv`: 100 unique expected genes: 7 common, 47 short-lived, and 46 long-lived.
- `tests/test_source.R`: reusable-function tests.
- `tests/test_extended_qc.R`: alarm, gene-quality, spatial, ranking, concordance, and artifact-contract tests.
- `config/extended_qc_defaults.tsv`: versioned extended-QC thresholds and deterministic seed.
- `config/subset_qc_reference.tsv`: fixed 500-cell subset review-burden reference.
- `shell/run_notebook_qc_hpc.sh` and `slurm/scwat_notebook_qc.sbatch`: full-data HPC execution.

Raw data, generated outputs, temporary files, R libraries, logs, and executed notebooks are not stored in Git.

## Section notebook inputs

For `REGION_ID=Region_1`, `Region_2`, `Region_3`, or `Region_4`, the notebook discovers exactly one matching directory below `INPUT_ROOT`. Each directory must contain:

- `experiment.xenium`
- `metrics_summary.csv`
- `analysis_summary.html`
- `gene_panel.json`
- `cells.csv.gz`
- `cell_feature_matrix/features.tsv.gz`
- `cell_feature_matrix/barcodes.tsv.gz`
- `cell_feature_matrix/matrix.mtx.gz`

`FULL_HPC` extended QC additionally requires `transcripts.parquet` in each section directory and the R packages `arrow`, `dplyr`, and `RANN`. `LOCAL_SUBSET` does not fabricate transcript-QV results: it writes `NOT_RUN_LOCAL_SUBSET` explicitly.

If `METADATA_PATH` is empty, deterministic synthetic metadata is used for testing and biological interpretation remains blocked.

## Section outputs

Each section writes below:

`${PROJECT_ROOT}/adipose_analysis/scwat_qc_outputs/<RUN_LABEL>/sections/<REGION_ID>/`

Core outputs include configuration/environment records, manifest, inventory, integrity table, feature types, panel reconciliation, alarms, QC thresholds/summary, gzipped cell-QC metadata, a raw-count-preserving sparse RDS, readiness gates, session information, and section PDF/300-dpi PNG figures.

Extended outputs include `extended_qc_preflight.tsv`, direct `cycle_alarm_evidence.tsv`, `gene_transcript_quality.tsv`, global/edge-density/hotspot spatial tables, gzipped cell-level spatial annotations, a morphology-review manifest, `extended_qc_status.tsv`, and `<REGION_ID>_extended_spatial_qc.pdf`. Exact cycle identity is never inferred from these tables and requires 10x diagnostics.

QC is calculated separately for each section. No cell is automatically removed. Droplet-specific mitochondrial and `scDblFinder` filters are intentionally not used for image-based Xenium segmentation.

## Slide-summary outputs

The summary notebook requires all four section bundles and writes below `<RUN_ROOT>/slide_summary/`:

- combined QC summary, thresholds, alarms, readiness, and gzipped cell metadata;
- `slide_qc_summary.rds`;
- `slide_qc_status.tsv`;
- a multi-page PDF and six 300-dpi PNG figures;
- `sessionInfo.txt`.

The extended slide contract additionally writes combined alarm and gene-quality tables, the complete `candidate_cycle_affected_genes.tsv`, the filtered `candidate_cycle_affected_genes_affected_only.tsv`, subset/full ranking and agreement tables, combined spatial diagnostics/hotspots, within-mouse section/gene concordance tables, `extended_slide_qc_status.tsv`, and `figures/scwat_extended_qc_diagnostics.pdf`. The filtered table contains only section-gene comparisons crossing the prespecified depletion and/or Q20-loss thresholds. Candidate genes remain labelled `CANDIDATE_NOT_CONFIRMED`; spatial hotspot labels require morphology/image review; concordance thresholds are advisory and do not change readiness.

## Local subset validation

Local execution uses only `adipose_analysis/subset_input/adipose_data` and writes by default to `adipose_analysis/scwat_qc_outputs/local_extended_qc_test`. All temporary files are forced below the D: project. The optional `-RepoRoot` argument allows verification from an isolated D:-local worktree.

```powershell
& 'C:\Windows\System32\WindowsPowerShell\v1.0\powershell.exe' -NoProfile -ExecutionPolicy Bypass `
  -File 'D:\Xiaonan\CODEX_projects\Yanan_Xenium\adipose_analysis\YNH_Xenium_scWAT\scripts\execute_local_subset.ps1' `
  -ProjectRoot 'D:\Xiaonan\CODEX_projects\Yanan_Xenium' `
  -RunLabel 'local_extended_qc_test' `
  -RegionId 'ALL'
```

Use `-RegionId Region_1`, `Region_2`, `Region_3`, or `Region_4` to execute only one named section notebook. Use `-RegionId SUMMARY` after all four section bundles exist under the same run label. `ALL` runs the four sections and then the summary.

The local computer does not have IRkernel/Jupyter notebook packages. Therefore, local validation evaluates the R cells sequentially in one clean R environment and saves executed notebook JSON. HPC uses the standard registered Jupyter R kernel.

## HPC preflight and submission

Required HPC software: `Rscript`, `python3`, `jupyter`, registered kernelspec `ir`, and R packages `IRkernel`, `Matrix`, `jsonlite`, `ggplot2`, `arrow`, `dplyr`, and `RANN`. The runner installs nothing. Each Region 1-4 directory must contain `transcripts.parquet`; Arrow projects and aggregates this input before collection and never loads the full transcript table into R memory.

### Chunk 1 - Environment and D/HPC-local paths

```bash
PROJECT_ROOT=/dssg/home/acct-svetoslav_chakarov/svetoslav_chakarov/Lab_members/Yanan_Hu/YNH_Xenium
export PROJECT_ROOT
export PIPELINE_REPO="${PROJECT_ROOT}/adipose_analysis/YNH_Xenium_scWAT"
export INPUT_ROOT="${PROJECT_ROOT}/adipose_data"
export RUN_LABEL="full_notebook_qc_v2"
export METADATA_PATH="${PIPELINE_REPO}/config/scwat_sample_manifest.tsv"
export RUN_ROOT="${PROJECT_ROOT}/adipose_analysis/scwat_qc_outputs/${RUN_LABEL}"
export TMPDIR="${PROJECT_ROOT}/adipose_analysis/tmp"
export TMP="${TMPDIR}"
export TEMP="${TMPDIR}"
export R_LIBS_USER="${PROJECT_ROOT}/adipose_analysis/R_libs"
mkdir -p "${TMPDIR}" "${R_LIBS_USER}" "${PROJECT_ROOT}/adipose_analysis/scwat_qc_logs" "${RUN_ROOT}/executed_notebooks"
```

Inputs: full Xenium sections below `${INPUT_ROOT}`, verified manifest, repository notebooks/config. Outputs: `${RUN_ROOT}`; temporary/R-library/log files remain below `${PROJECT_ROOT}/adipose_analysis`.

### Chunk 2 - Package, kernel, transcript, size, and shell preflight

```bash
command -v Rscript
command -v python3
command -v jupyter
jupyter kernelspec list
Rscript -e 'p <- c("IRkernel","Matrix","jsonlite","ggplot2","arrow","dplyr","RANN"); ok <- vapply(p, requireNamespace, logical(1), quietly=TRUE); print(data.frame(package=p, available=ok)); stopifnot(all(ok))'
for region in Region_1 Region_2 Region_3 Region_4; do
  region_dir=$(find "${INPUT_ROOT}" -mindepth 1 -maxdepth 1 -type d -name "*__${region}__*" -print -quit)
  test -n "${region_dir}"
  test -f "${region_dir}/transcripts.parquet"
  du -sh "${region_dir}/transcripts.parquet" "${region_dir}/cell_feature_matrix" "${region_dir}/cells.csv.gz"
done
df -h "${PROJECT_ROOT}"
bash -n "${PIPELINE_REPO}/shell/run_notebook_qc_hpc.sh"
bash -n "${PIPELINE_REPO}/slurm/scwat_notebook_qc.sbatch"
```

### Chunk 3 - Run and inspect one section notebook

Set `REGION_ID` to run one committed section notebook. Inputs are the matching Xenium section directory and `transcripts.parquet`. Outputs are `${RUN_ROOT}/sections/<Region_ID>/` and `${RUN_ROOT}/executed_notebooks/<Region_ID>.executed.ipynb`.

```bash
REGION_ID=Region_1 bash "${PIPELINE_REPO}/shell/run_notebook_qc_hpc.sh"
```

After reviewing Region 1, repeat with the same `RUN_LABEL`:

```bash
REGION_ID=Region_2 bash "${PIPELINE_REPO}/shell/run_notebook_qc_hpc.sh"
REGION_ID=Region_3 bash "${PIPELINE_REPO}/shell/run_notebook_qc_hpc.sh"
REGION_ID=Region_4 bash "${PIPELINE_REPO}/shell/run_notebook_qc_hpc.sh"
```

### Chunk 4 - Run all sections and summary in one job

`ALL` executes the four named notebooks in order and then the summary notebook.

```bash
REGION_ID=ALL bash "${PIPELINE_REPO}/shell/run_notebook_qc_hpc.sh"
```

Alternatively submit the configured Slurm wrapper:

```bash
sbatch "${PIPELINE_REPO}/slurm/scwat_notebook_qc.sbatch"
```

### Chunk 5 - Final summary notebook only

Use `SUMMARY` only after all four full-HPC section bundles exist under the same `RUN_LABEL`. Inputs are `${RUN_ROOT}/sections/Region_1` through `Region_4`. Outputs are `${RUN_ROOT}/slide_summary/` and `${RUN_ROOT}/executed_notebooks/slide_summary.executed.ipynb`.

```bash
REGION_ID=SUMMARY bash "${PIPELINE_REPO}/shell/run_notebook_qc_hpc.sh"
```

### Chunk 6 - Post-run validation and output locations

```bash
Rscript -e "source('${PIPELINE_REPO}/R/source.R'); for(r in paste0('Region_',1:4)) stopifnot(validate_extended_section_artifacts(file.path('${RUN_ROOT}','sections',r),r,'FULL_HPC',stop_on_error=TRUE)); stopifnot(validate_extended_slide_qc_artifacts('${RUN_ROOT}',stop_on_error=TRUE)); stopifnot(validate_evidence_only_qc_artifacts('${RUN_ROOT}',stop_on_error=TRUE))"
find "${RUN_ROOT}/slide_summary" -maxdepth 2 -type f -print | sort
find "${RUN_ROOT}/downstream_inputs" -maxdepth 1 -type f -print | sort
```

The evidence-only downstream contract is in `${RUN_ROOT}/slide_summary/cell_downstream_masks.tsv.gz`, `section_downstream_decision.tsv`, `gene_downstream_decision.tsv`, `eos_gene_decision_summary.tsv`, `hotspot_sensitivity_decision.tsv`, and `evidence_only_qc_release.tsv`. Final region inputs are `${RUN_ROOT}/downstream_inputs/Region_1.downstream_input.rds` through `Region_4.downstream_input.rds`, indexed by `downstream_input_manifest.tsv`.

## Evidence-only downstream decisions

- Region 1: `PRIMARY_CONDITIONAL`; Region 2: `PRIMARY_CONDITIONAL`; Region 3: `PRIMARY`; Region 4: `SENSITIVITY_ONLY`.
- `primary_include` excludes core-QC failures, segmentation multiplets, and high-control cells. `strict_include` excludes every review-flagged cell. `hotspot_sensitivity_include` additionally excludes Region 3 morphology-review hotspot cells without removing them from primary analysis.
- All genes remain in `RAW_COMPLETE_PANEL`. The primary feature eligibility field is `PROVISIONAL_PRIMARY_FEATURES`; its zero-alarm subset is `CONSERVATIVE_NO_SIGNAL_DETECTED`. Genes recurring in at least two alarm-positive sections are `TECHNICAL_RISK_SENSITIVITY_ONLY` and cannot define primary clusters.
- Region 4 never contributes to reference discovery or primary gene-level results. Its downstream bundle requires mapping to the finalized Region 1-3 reference and an `Uncertain` label for insufficient-confidence assignments.
- Phase 0-2 produces inputs rather than PCA, clusters, Eos states, or Region 4 mapping. Therefore the overall primary release remains `PENDING_DOWNSTREAM_ANALYSIS` until the automated downstream stability gates are supplied.

## Current scientific gate

The metadata-aware subset rerun reproduces 500 cells per section and 2,000 total cells. The panel reconciles 100/100. Metadata passes in all four sections. Regions 1, 2, and 4 contain poor-quality-cycle ERROR alarms and remain on `HOLD`; Region 3 is now `PASS`. Overall slide readiness remains `HOLD`, so cross-section biological interpretation is still blocked.

The applied study design is: sections 62308 and 62309 are from Mouse 1; sections 62310 and 62311 are from Mouse 2. All samples are untreated WT scWAT from normal 8-week-old mice. Left/right is not a design factor. Section is the technical processing unit and mouse is the biological replicate.

The current reader-facing QC report is `notebooks/02_slide_QC_summary.ipynb`. It reports direct alarm evidence, clearly labelled candidate-gene evidence, spatial and within-mouse diagnostics, section decisions, and the exact full-HPC/10x evidence still required.

### Local validation evidence (2026-08-14)

- Reusable and extended R tests exited successfully.
- The extended section notebook executed 14 R code cells independently for each of Regions 1-4.
- The extended slide-summary notebook executed 11 R code cells after all core and extended section artifact checks passed.
- Core-pass counts were 493, 500, 497, and 499; review-flag counts were 37, 16, 27, and 11.
- All four sparse RDS objects reloaded with 500 cells; the combined summary reloaded with 2,000 cells and four unique regions.
- Source and executed notebook JSON contracts passed structural validation.
- All four regenerated manifests contain `VERIFIED_USER_SUPPLIED`; all four metadata gates are `PASS`; Region 3 overall readiness changed from `PENDING` to `PASS`.
- Local transcript-QV and subset-versus-full comparisons are explicitly `NOT_RUN_LOCAL_SUBSET`; no exact-cycle or confirmed affected-gene claim is made locally.
- Direct poor-cycle alarms remain present in Regions 1, 2, and 4. Local spatial tables and two advisory within-mouse section comparisons were generated and reload-validated.
- Full-data HPC execution and `bash -n` remain pending because Bash/Jupyter/IRkernel are unavailable on the local Windows test environment.
