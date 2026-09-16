# YNH Xenium scWAT

## Project abstract

This project investigates the cellular composition and spatial organization of normal mouse subcutaneous white adipose tissue (scWAT) using the 10x Genomics Xenium platform. The dataset contains four tissue sections collected from two untreated, 8-week-old wild-type mice. Sections 62308 and 62309 are technical sections from Mouse 1; sections 62310 and 62311 are technical sections from Mouse 2. Left/right is not used as an experimental factor. The mouse is the biological replicate, whereas each tissue section is analysed separately as a technical and spatial unit.

| Region | Xenium section | Mouse | Genotype | Treatment | Age | Analysis role |
|---|---:|---|---|---|---:|---|
| Region 1 | 62308 | Mouse 1 | WT | None | 8 weeks | `PRIMARY_CONDITIONAL` |
| Region 2 | 62309 | Mouse 1 | WT | None | 8 weeks | `PRIMARY_CONDITIONAL` |
| Region 3 | 62310 | Mouse 2 | WT | None | 8 weeks | `PRIMARY`; initial anchor reference |
| Region 4 | 62311 | Mouse 2 | WT | None | 8 weeks | `SENSITIVITY_ONLY`; mapping-only |

The Xenium panel contains 479 genes. The analysis first performs technical and cell-level QC independently for every section, then summarizes QC across the slide. Downstream analyses use Region 3 as the initial anchor: Regions 1 and 2 are admitted only if they agree with the anchor, while Region 4 is mapped to the finalized reference without influencing reference construction, PCA, integration, clustering, marker selection, or reference labels.

The biological focus is Eosinophil identification and the spatial comparison of short-lived-like and long-lived-like Eosinophil states. The panel includes a prespecified 100-gene Eosinophil set comprising 7 common genes, 47 short-lived genes, and 46 long-lived genes. Eosinophil identity must be validated before Eosinophil-state, nearest-neighbour, or ligand–receptor results are interpreted.

The verified sample metadata are stored in [`config/scwat_sample_manifest.tsv`](config/scwat_sample_manifest.tsv). Raw Xenium files, large analysis objects, executed notebooks, and generated results are not stored in this Git repository.

## Notebooks and source code

### Initial QC notebooks

These notebooks perform Phase 0–2 QC. Run the four region notebooks separately and then run the slide-level summary.

| Notebook | Description |
|---|---|
| [`notebooks/01_QC_Region1.ipynb`](notebooks/01_QC_Region1.ipynb) | Configuration, file-integrity checks, panel reconciliation, cell-level QC, spatial review flags, and output validation for Region 1. |
| [`notebooks/01_QC_Region2.ipynb`](notebooks/01_QC_Region2.ipynb) | The same region-specific QC workflow for Region 2. |
| [`notebooks/01_QC_Region3.ipynb`](notebooks/01_QC_Region3.ipynb) | The same workflow for Region 3, which is subsequently used as the anchor region. |
| [`notebooks/01_QC_Region4.ipynb`](notebooks/01_QC_Region4.ipynb) | The same workflow for Region 4, which remains sensitivity/mapping-only. |
| [`notebooks/02_slide_QC_summary.ipynb`](notebooks/02_slide_QC_summary.ipynb) | Combines the four section-level QC bundles, summarizes slide-level QC statistics and plots, and creates the downstream cell masks and section/gene decisions. |

The fixed core cell-QC limits used by these notebooks are stored in [`config/fixed_cell_qc_thresholds.tsv`](config/fixed_cell_qc_thresholds.tsv):

```text
5 < nFeature_Xenium < 200
10 < nCount_Xenium < 1000
```

### Full-panel region notebooks

The active B2 workflow analyses all 479 panel genes. For every region, it is divided into three notebooks that must be run in the order shown below.

| Notebook pattern | Regions | Description |
|---|---|---|
| `notebooks/B2_RegionX_all_QCpass_479.ipynb` | 1–4 | Analyses all cells passing `primary_include_revised`, performs normalization, PCA, clustering and annotation, and defines the largest non-noise DBSCAN lymph-node domain. It writes the frozen domain manifest used by the next two notebooks. |
| `notebooks/B2_RegionX_adipose_only_479.ipynb` | 1–4 | Excludes cells in the frozen lymph-node domain and repeats the analysis for the adipose-tissue compartment. |
| `notebooks/B2_RegionX_lymph_node_only_479.ipynb` | 1–4 | Restricts the analysis to the frozen lymph-node domain and stops cleanly if the domain contains too few cells. |

Replace `X` with `1`, `2`, `3`, or `4`. For example, the Region 3 notebooks are:

- [`notebooks/B2_Region3_all_QCpass_479.ipynb`](notebooks/B2_Region3_all_QCpass_479.ipynb)
- [`notebooks/B2_Region3_adipose_only_479.ipynb`](notebooks/B2_Region3_adipose_only_479.ipynb)
- [`notebooks/B2_Region3_lymph_node_only_479.ipynb`](notebooks/B2_Region3_lymph_node_only_479.ipynb)

The twelve active B2 notebooks are generated from [`scripts/build_b2_split_notebooks.py`](scripts/build_b2_split_notebooks.py). Changes shared across regions or tissue branches should be made in the builder and regenerated consistently.

The files `notebooks/B2_Region1_primary_479.ipynb` through `notebooks/B2_Region4_primary_479.ipynb` are older monolithic development notebooks retained for provenance; they are not the current canonical execution path.

### Reference and marker notebooks

| Notebook | Description |
|---|---|
| [`notebooks/CellType_Markers.ipynb`](notebooks/CellType_Markers.ipynb) | Reviews and organizes marker genes used for scWAT cell-type annotation. |
| [`notebooks/RefCheck_Inhouse_scWAT.ipynb`](notebooks/RefCheck_Inhouse_scWAT.ipynb) | Inspects the in-house scWAT reference, its metadata, cell labels, and suitability for reference-based annotation. |
| [`notebooks/RefCheck_WangScience2025.ipynb`](notebooks/RefCheck_WangScience2025.ipynb) | Inspects the Wang Science 2025 reference, including age groups, cell labels, Eosinophil representation, and compatibility with the Xenium panel. |

### Reusable R source

| File | Description |
|---|---|
| [`R/source.R`](R/source.R) | Active reusable functions called by the QC and downstream notebooks. It contains path and environment checks, Xenium import, mask validation, plotting, clustering diagnostics, cell annotation, Eosinophil scoring, spatial-neighbour analysis, reference mapping, CellChat adapters, safe table export, and artifact validation. Functions include roxygen-style documentation describing their inputs and outputs. |
| [`R/source_bk.R`](R/source_bk.R) | Archived or obsolete functions that are no longer required by the active notebooks. This file is retained for provenance and is not sourced by the current workflow. |

### Configuration files

| File | Description |
|---|---|
| [`config/scwat_sample_manifest.tsv`](config/scwat_sample_manifest.tsv) | Verified mouse, section, condition, and analysis-role metadata. |
| [`config/fixed_cell_qc_thresholds.tsv`](config/fixed_cell_qc_thresholds.tsv) | Versioned cell-level QC thresholds. |
| [`config/extended_qc_defaults.tsv`](config/extended_qc_defaults.tsv) | Thresholds and deterministic settings for extended QC. |
| [`config/eos_gene_sets.tsv`](config/eos_gene_sets.tsv) | Prespecified common, short-lived, and long-lived Eosinophil gene sets. |
| [`config/subset_qc_reference.tsv`](config/subset_qc_reference.tsv) | Fixed reference values used to check the local QC subset. |

### Pipeline and validation code

| Location | Description |
|---|---|
| [`scripts/`](scripts/) | Notebook builders, local notebook execution, rendering, and source-code audit utilities. |
| [`shell/run_notebook_qc_hpc.sh`](shell/run_notebook_qc_hpc.sh) | Executes a selected region QC notebook, all regions, or the slide summary on HPC. |
| [`slurm/scwat_notebook_qc.sbatch`](slurm/scwat_notebook_qc.sbatch) | Slurm submission wrapper for the HPC QC workflow. |
| [`tests/`](tests/) | R and Python tests for reusable functions, fixed QC masks, notebook contracts, B2 notebook generation, Eosinophil helpers, tissue-domain helpers, and local subset smoke tests. |
| [`docs/`](docs/) | Method designs, implementation plans, and validation audits. |
| [`SCWAT_QUERY_AND_RESEARCH_PLAN.md`](SCWAT_QUERY_AND_RESEARCH_PLAN.md) | Detailed project decisions and chronological research/implementation record. |

Repository: <https://github.com/SharonWang/YNH_Xenium_scWAT>
