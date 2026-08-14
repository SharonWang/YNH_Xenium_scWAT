# YNH Xenium scWAT: Phase 0-2 QC notebooks

This repository contains the reusable scWAT Xenium QC code. Run one parameterized R notebook independently for each of the four sections, then run the slide-level summary notebook.

## Files

- `R/source.R`: general Xenium path, metadata, integrity, panel, sparse-import, QC, aggregation, artifact, and plotting functions.
- `notebooks/01_section_phase0_2_QC.ipynb`: combined Phase 0-2 notebook for one section.
- `notebooks/02_slide_QC_summary.ipynb`: final four-section QC tables and Cell-inspired figures.
- `config/eos_gene_sets.tsv`: 100 unique expected genes: 7 common, 47 short-lived, and 46 long-lived.
- `tests/test_source.R`: reusable-function tests.
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

If `METADATA_PATH` is empty, deterministic synthetic metadata is used for testing and biological interpretation remains blocked.

## Section outputs

Each section writes below:

`${PROJECT_ROOT}/adipose_analysis/scwat_qc_outputs/<RUN_LABEL>/sections/<REGION_ID>/`

Outputs include configuration/environment records, manifest, inventory, integrity table, feature types, panel reconciliation, alarms, QC thresholds/summary, gzipped cell-QC metadata, a raw-count-preserving sparse RDS, readiness gates, session information, and section PDF/300-dpi PNG figures.

QC is calculated separately for each section. No cell is automatically removed. Droplet-specific mitochondrial and `scDblFinder` filters are intentionally not used for image-based Xenium segmentation.

## Slide-summary outputs

The summary notebook requires all four section bundles and writes below `<RUN_ROOT>/slide_summary/`:

- combined QC summary, thresholds, alarms, readiness, and gzipped cell metadata;
- `slide_qc_summary.rds`;
- `slide_qc_status.tsv`;
- a multi-page PDF and six 300-dpi PNG figures;
- `sessionInfo.txt`.

## Local subset validation

Local execution uses only `adipose_analysis/subset_input/adipose_data` and writes to `adipose_analysis/scwat_qc_outputs/local_notebook_test`. All temporary files are forced below the D: project.

```powershell
& 'C:\Windows\System32\WindowsPowerShell\v1.0\powershell.exe' -NoProfile -ExecutionPolicy Bypass `
  -File 'D:\Xiaonan\CODEX_projects\Yanan_Xenium\adipose_analysis\YNH_Xenium_scWAT\scripts\execute_local_subset.ps1'
```

The local computer does not have IRkernel/Jupyter notebook packages. Therefore, local validation evaluates the R cells sequentially in one clean R environment and saves executed notebook JSON. HPC uses the standard registered Jupyter R kernel.

## HPC preflight and submission

Required HPC software: `Rscript`, `python3`, `jupyter`, registered kernelspec `ir`, and R packages `IRkernel`, `Matrix`, `jsonlite`, `ggplot2`. The runner installs nothing; load/install these using the site-approved environment before submission.

```bash
PROJECT_ROOT=/dssg/home/acct-svetoslav_chakarov/svetoslav_chakarov/Lab_members/Yanan_Hu
mkdir -p "${PROJECT_ROOT}/adipose_analysis/scwat_qc_logs"
bash -n "${PROJECT_ROOT}/adipose_analysis/YNH_Xenium_scWAT/shell/run_notebook_qc_hpc.sh"
bash -n "${PROJECT_ROOT}/adipose_analysis/YNH_Xenium_scWAT/slurm/scwat_notebook_qc.sbatch"
sbatch "${PROJECT_ROOT}/adipose_analysis/YNH_Xenium_scWAT/slurm/scwat_notebook_qc.sbatch"
```

To use real metadata later:

```bash
export METADATA_PATH="${PROJECT_ROOT}/adipose_analysis/metadata/scwat_sample_manifest.tsv"
export RUN_LABEL="full_notebook_qc_real_metadata_v1"
bash "${PROJECT_ROOT}/adipose_analysis/YNH_Xenium_scWAT/shell/run_notebook_qc_hpc.sh"
```

## Current scientific gate

The validated subset reproduces 500 cells per section and 2,000 total cells. The panel reconciles 100/100. Regions 1, 2, and 4 contain poor-quality-cycle ERROR alarms and are on `HOLD`; Region 3 is `PENDING` because metadata is synthetic. Overall slide readiness is `HOLD`, so diagnostic QC is valid but biological interpretation is not.
