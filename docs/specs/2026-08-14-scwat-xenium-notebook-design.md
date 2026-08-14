# scWAT Xenium Phase 0-2 Notebook Design

**Date:** 2026-08-14  
**Scope:** scWAT only  
**Canonical repository:** `D:/Xiaonan/CODEX_projects/Yanan_Xenium/adipose_analysis/YNH_Xenium_scWAT`  
**Local data root:** `D:/Xiaonan/CODEX_projects/Yanan_Xenium`  
**HPC project root:** `/dssg/home/acct-svetoslav_chakarov/svetoslav_chakarov/Lab_members/Yanan_Hu`

## Objective

Provide one parameterized R-kernel Jupyter notebook that runs Xenium Phases 0-2 independently for each of the four scWAT sections, followed by one slide-level QC summary notebook. General Xenium logic lives in one reusable `source.R` file inside this repository.

## Deliverables

1. `R/source.R`: reusable, tissue-agnostic functions for path safety, provenance, integrity, panel reconciliation, sparse import, cell QC, aggregation, artifact validation, and plotting.
2. `notebooks/01_section_phase0_2_QC.ipynb`: parameterized section notebook. Set `REGION_ID` to `Region_1`, `Region_2`, `Region_3`, or `Region_4` and run once per section.
3. `notebooks/02_slide_QC_summary.ipynb`: validates and summarizes all four completed section runs.
4. `tests/test_source.R`: unit tests for reusable functions.
5. `shell/run_notebook_qc_hpc.sh` and `slurm/scwat_notebook_qc.sbatch`: HPC execution wrappers.
6. `README.md`: exact inputs, outputs, manual cells, and HPC commands.

The repository contains code and documentation only. Data, generated outputs, R libraries, logs, caches, and temporary files stay below the external project roots and are excluded from Git.

## Section notebook parameters

- `PROJECT_ROOT`: defaults to the supplied HPC project root.
- `PIPELINE_REPO`: defaults to `${PROJECT_ROOT}/adipose_analysis/YNH_Xenium_scWAT`.
- `INPUT_ROOT`: defaults to `${PROJECT_ROOT}/adipose_data`.
- `REGION_ID`: one of `Region_1` through `Region_4`.
- `RUN_LABEL`: stable output namespace.
- `METADATA_PATH`: optional; empty creates deterministic synthetic placeholder metadata.
- `EXPECTED_SECTION_COUNT`: exactly `4`.
- `SEED`: `20260814`.
- `STRICT_MODE`: defaults to `FALSE`, allowing diagnostic outputs when a biological gate is on `HOLD`.

For local validation, `PROJECT_ROOT` remains `D:/Xiaonan/CODEX_projects/Yanan_Xenium` and `INPUT_ROOT` is the deterministic subset at `adipose_analysis/subset_input/adipose_data`. No project or temporary files may be written to C:.

## Section inputs

Each run discovers exactly one directory matching the chosen region below `INPUT_ROOT` and validates:

- `experiment.xenium`;
- `metrics_summary.csv`;
- `analysis_summary.html`;
- `gene_panel.json`;
- `cells.csv.gz`;
- `cell_feature_matrix/features.tsv.gz`;
- `cell_feature_matrix/barcodes.tsv.gz`;
- `cell_feature_matrix/matrix.mtx.gz`.

The expected Eos panel table is `config/eos_gene_sets.tsv` in the repository.

## Section workflow

The notebook runs top-to-bottom:

1. validate project, repository, input, output, and temporary paths;
2. verify the R kernel and required packages;
3. source `R/source.R`;
4. discover and validate the selected section;
5. create or validate its one-row metadata record;
6. inventory required files with sizes, timestamps, and checksums;
7. reconcile matrix dimensions, features, barcodes, and cell rows;
8. extract Xenium analysis alarms;
9. reconcile installed and expected panel genes;
10. import raw gene-expression counts sparsely;
11. calculate section-specific QC thresholds and review flags;
12. create section-level QC plots;
13. preserve raw sparse counts and all cells in an RDS;
14. write readiness gates, output checks, and `sessionInfo()`.

## Section outputs

Each run writes only below:

`${PROJECT_ROOT}/adipose_analysis/scwat_qc_outputs/<RUN_LABEL>/sections/<REGION_ID>/`

Required files:

- `configuration.tsv`;
- `section_manifest.tsv`;
- `environment_preflight.tsv`;
- `file_inventory.tsv`;
- `integrity_summary.tsv`;
- `feature_type_summary.tsv`;
- `panel_reconciliation.tsv`;
- `analysis_alerts.tsv`;
- `qc_thresholds.tsv`;
- `qc_summary.tsv`;
- `cell_qc_metadata.tsv.gz`;
- `<REGION_ID>.phase0_2_qc.rds`;
- `section_readiness_gates.tsv`;
- section PDF and 300-dpi PNG figures;
- `sessionInfo.txt`.

Missing or corrupt required inputs stop the run. Xenium `ERROR` alarms set biological readiness to `HOLD` but do not suppress diagnostic outputs when `STRICT_MODE=FALSE`.

## Reusable function contract

`R/source.R` groups functions by responsibility while remaining one requested source file:

- runtime and path validation;
- deterministic metadata contracts;
- section discovery and required-file inventory;
- Matrix Market, barcode, feature, and cell-table integrity checks;
- sparse Xenium import;
- analysis-alarm extraction;
- panel reconciliation;
- robust per-section cell QC;
- section artifact writing and reload validation;
- four-section completeness validation;
- slide-level aggregation;
- Cell-inspired plotting helpers.

Functions accept explicit arguments and return R objects. They do not rely on notebook globals. Functions that write files require a previously validated output root. Raw Xenium inputs are never modified.

## QC policy

- QC is computed independently for each section.
- Count and detected-feature limits use robust median/MAD intervals with quantile fallback.
- Nucleus absence, multiple nuclei, cell-area outliers, high control fraction, and segmentation-multiplet evidence remain separate flags.
- No cell is automatically deleted.
- Raw sparse gene-expression counts remain unchanged in saved objects.
- Droplet-specific mitochondrial filters and `scDblFinder` are not used for image-based Xenium segmentation.
- Synthetic metadata blocks biological interpretation.
- Poor-quality-cycle errors remain a `HOLD` gate pending cycle/codeword/gene-level review.

## Slide summary notebook

The summary notebook requires exactly four completed, uniquely identified section directories. It combines QC summaries, thresholds, gates, alarms, and cell metadata; compares distributions without treating cells as biological replicates; writes slide-level TSV/RDS outputs; and creates cross-section figures. It performs no mouse/side inference while metadata is synthetic.

## Figure style

Figures use a restrained Cell-inspired scientific style:

- white background, thin dark axes, compact typography, and clear units;
- one consistent colorblind-safe palette for Regions 1-4;
- violin/box/point distribution summaries;
- directly labeled QC pass/review proportions;
- equal-coordinate spatial centroid plots with review flags overlaid;
- PDF and 300-dpi PNG outputs;
- no decorative gradients, 3D effects, or truncated axes that alter interpretation.

## Reproducibility and HPC execution

- Both notebooks declare an R kernel.
- Required R packages are `IRkernel`, `Matrix`, `jsonlite`, and `ggplot2`.
- `TMPDIR`, `TMP`, `TEMP`, and `R_LIBS_USER` are set below the project root before notebook execution.
- The HPC runner executes the section notebook four times with parameter injection, then the slide summary notebook.
- Executed notebook copies and logs are stored under the external run directory, not Git.

## Validation criteria

- unit tests pass for reusable calculation and contract functions;
- notebook JSON/schema validation passes;
- the section notebook executes top-to-bottom for all four 500-cell local subsets;
- the summary notebook executes across those four completed outputs;
- all four section RDS objects and the slide summary reload successfully;
- no local output or temporary path escapes the D: project;
- HPC performs `bash -n`, kernel/package preflight, and full-data runtime validation.

Local testing cannot certify HPC module, kernel, scheduler, or memory configuration; those remain explicit HPC validation steps.
