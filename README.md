# YNH Xenium scWAT: Phase 0-2 QC notebooks

This repository contains the reusable scWAT Xenium QC code. Run one parameterized R notebook independently for each of the four sections, then run the slide-level summary notebook.

## Files

- `R/source.R`: general Xenium path, metadata, integrity, panel, sparse-import, QC, aggregation, artifact, and plotting functions.
- `notebooks/01_section_phase0_2_QC.ipynb`: combined Phase 0-2 notebook for one section.
- `notebooks/02_slide_QC_summary.ipynb`: final four-section QC tables and Cell-inspired figures.
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
  -RunLabel 'local_extended_qc_test'
```

The local computer does not have IRkernel/Jupyter notebook packages. Therefore, local validation evaluates the R cells sequentially in one clean R environment and saves executed notebook JSON. HPC uses the standard registered Jupyter R kernel.

## HPC preflight and submission

Required HPC software: `Rscript`, `python3`, `jupyter`, registered kernelspec `ir`, and R packages `IRkernel`, `Matrix`, `jsonlite`, `ggplot2`, `arrow`, `dplyr`, and `RANN`. The runner installs nothing. Each Region 1-4 directory must contain `transcripts.parquet`; Arrow projects and aggregates this input before collection and never loads the full transcript table into R memory.

### Chunk 1 - Environment and D/HPC-local paths

```bash
PROJECT_ROOT=/dssg/home/acct-svetoslav_chakarov/svetoslav_chakarov/Lab_members/Yanan_Hu/YNH_Xenium
export PROJECT_ROOT
export PIPELINE_REPO="${PROJECT_ROOT}/adipose_analysis/YNH_Xenium_scWAT"
export INPUT_ROOT="${PROJECT_ROOT}/adipose_data"
export RUN_LABEL="full_extended_qc_v1"
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

### Chunk 3 - One-section interactive checkpoint

This tests Region 1 first and writes its executed notebook plus section bundle under `${RUN_ROOT}`.

```bash
python3 "${PIPELINE_REPO}/scripts/render_notebooks.py" --inject \
  "${PIPELINE_REPO}/notebooks/01_section_phase0_2_QC.ipynb" \
  "${RUN_ROOT}/executed_notebooks/Region_1.executed.ipynb" \
  --set "PROJECT_ROOT=${PROJECT_ROOT}" --set "PIPELINE_REPO=${PIPELINE_REPO}" \
  --set "INPUT_ROOT=${INPUT_ROOT}" --set "REGION_ID=Region_1" --set "RUN_LABEL=${RUN_LABEL}" \
  --set "METADATA_PATH=${METADATA_PATH}" --set "EXPECTED_SECTION_COUNT=4L" \
  --set "SEED=20260814L" --set "STRICT_MODE=FALSE" --set "EXTENDED_QC_MODE=FULL_HPC" \
  --set "EXTENDED_QC_CONFIG_PATH=${PIPELINE_REPO}/config/extended_qc_defaults.tsv"
jupyter nbconvert --execute --to notebook --inplace \
  --ExecutePreprocessor.kernel_name=ir --ExecutePreprocessor.timeout=-1 \
  "${RUN_ROOT}/executed_notebooks/Region_1.executed.ipynb"
Rscript -e "source('${PIPELINE_REPO}/R/source.R'); stopifnot(validate_extended_section_artifacts('${RUN_ROOT}/sections/Region_1','Region_1','FULL_HPC',stop_on_error=TRUE))"
```

### Chunk 4 - Full four-section run

The runner repeats/overwrites Region 1 within the same run label, then processes Regions 2-4 and the slide summary.

```bash
bash "${PIPELINE_REPO}/shell/run_notebook_qc_hpc.sh"
```

Alternatively submit the configured Slurm wrapper:

```bash
sbatch "${PIPELINE_REPO}/slurm/scwat_notebook_qc.sbatch"
```

### Chunk 5 - Summary-only rerun

Use this only after all four full-HPC section bundles validate.

```bash
python3 "${PIPELINE_REPO}/scripts/render_notebooks.py" --inject \
  "${PIPELINE_REPO}/notebooks/02_slide_QC_summary.ipynb" \
  "${RUN_ROOT}/executed_notebooks/slide_summary.executed.ipynb" \
  --set "PROJECT_ROOT=${PROJECT_ROOT}" --set "PIPELINE_REPO=${PIPELINE_REPO}" \
  --set "RUN_LABEL=${RUN_LABEL}" --set "EXPECTED_SECTION_COUNT=4L" \
  --set "METADATA_PATH=${METADATA_PATH}" \
  --set "EXTENDED_QC_CONFIG_PATH=${PIPELINE_REPO}/config/extended_qc_defaults.tsv" \
  --set "SUBSET_REFERENCE_PATH=${PIPELINE_REPO}/config/subset_qc_reference.tsv"
jupyter nbconvert --execute --to notebook --inplace \
  --ExecutePreprocessor.kernel_name=ir --ExecutePreprocessor.timeout=-1 \
  "${RUN_ROOT}/executed_notebooks/slide_summary.executed.ipynb"
```

### Chunk 6 - Post-run validation and output locations

```bash
Rscript -e "source('${PIPELINE_REPO}/R/source.R'); for(r in paste0('Region_',1:4)) stopifnot(validate_extended_section_artifacts(file.path('${RUN_ROOT}','sections',r),r,'FULL_HPC',stop_on_error=TRUE)); stopifnot(validate_extended_slide_qc_artifacts('${RUN_ROOT}',stop_on_error=TRUE))"
find "${RUN_ROOT}/slide_summary" -maxdepth 2 -type f -print | sort
```

Primary answers are in `${RUN_ROOT}/slide_summary/candidate_cycle_affected_genes.tsv`, `subset_full_qc_ranking.tsv`, `combined_spatial_qc.tsv`, `combined_spatial_hotspots.tsv`, `within_mouse_section_concordance.tsv`, and `figures/scwat_extended_qc_diagnostics.pdf`. Exact cycle identity still requires the separate 10x diagnostic output.

## Current scientific gate

The metadata-aware subset rerun reproduces 500 cells per section and 2,000 total cells. The panel reconciles 100/100. Metadata passes in all four sections. Regions 1, 2, and 4 contain poor-quality-cycle ERROR alarms and remain on `HOLD`; Region 3 is now `PASS`. Overall slide readiness remains `HOLD`, so cross-section biological interpretation is still blocked.

The applied study design is: sections 62308 and 62309 are from Mouse 1; sections 62310 and 62311 are from Mouse 2. All samples are untreated WT scWAT from normal 8-week-old mice. Left/right is not a design factor. Section is the technical processing unit and mouse is the biological replicate.

The updated portable technical QC report is at `reports/2026-08-14_scwat_qc_summary/report.html`; its canonical data/provenance specification is `reports/2026-08-14_scwat_qc_summary/artifact.json`. It reports verified local alarm/spatial/concordance evidence, a clearly labelled depletion-only candidate preview, and the exact full-HPC/10x evidence still required.

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
