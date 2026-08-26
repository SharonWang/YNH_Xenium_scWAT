#!/usr/bin/env python3
"""Build output-free scWAT initial-QC notebooks from the Colon reader flow."""

from __future__ import annotations

import json
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
NOTEBOOK_DIR = REPO / "notebooks"


def split_source(text: str) -> list[str]:
    lines = text.strip("\n").splitlines(keepends=True)
    if lines and not lines[-1].endswith("\n"):
        lines[-1] += "\n"
    return lines


def cell(cell_type: str, source: str, index: int) -> dict:
    result = {
        "cell_type": cell_type,
        "id": f"cell-{index:03d}",
        "metadata": {},
        "source": split_source(source),
    }
    if cell_type == "code":
        result.update(execution_count=None, outputs=[])
    return result


def notebook(cells: list[tuple[str, str]]) -> dict:
    return {
        "cells": [cell(kind, source, index) for index, (kind, source) in enumerate(cells, 1)],
        "metadata": {
            "kernelspec": {"display_name": "R", "language": "R", "name": "ir"},
            "language_info": {"name": "R", "mimetype": "text/x-r-source", "file_extension": ".r"},
        },
        "nbformat": 4,
        "nbformat_minor": 5,
    }


def region_cells(region: str) -> list[tuple[str, str]]:
    return [
        ("markdown", f"""# scWAT Xenium QC - {region}

## Goal

Run reproducible, non-destructive technical QC for **{region}**. This notebook is one of four structurally identical region notebooks; its region identifier is locked below. Full-data results must be produced on HPC. Local runs use a deterministic small subset only and cannot establish transcript-level or final slide readiness."""),
        ("markdown", """## Setup

The cell below defines HPC defaults and accepts explicit environment overrides for local subset validation. It refuses writable paths outside `adipose_analysis`. The pipeline installs nothing; missing packages cause a clear preflight stop."""),
        ("code", f'''REGION_ID <- "{region}"
EXECUTION_MODE <- toupper(Sys.getenv("SCWAT_QC_MODE", "FULL_HPC"))
stopifnot(EXECUTION_MODE %in% c("FULL_HPC", "LOCAL_SUBSET"))
PROJECT_ROOT <- Sys.getenv("SCWAT_PROJECT_ROOT", "/dssg/home/acct-svetoslav_chakarov/svetoslav_chakarov/Lab_members/Yanan_Hu/YNH_Xenium")
PIPELINE_REPO <- Sys.getenv("SCWAT_PIPELINE_REPO", file.path(PROJECT_ROOT, "adipose_analysis", "YNH_Xenium_scWAT"))
INPUT_ROOT <- Sys.getenv("SCWAT_INPUT_ROOT", file.path(PROJECT_ROOT, "adipose_data"))
RUN_LABEL <- Sys.getenv("SCWAT_RUN_LABEL", "scwat_qc_hpc")
RUN_ROOT <- Sys.getenv("SCWAT_RUN_ROOT", file.path(PROJECT_ROOT, "adipose_analysis", "scwat_qc_outputs", RUN_LABEL))
TEMP_ROOT <- Sys.getenv("SCWAT_TEMP_ROOT", file.path(PROJECT_ROOT, "adipose_analysis", "tmp", RUN_LABEL))
dir.create(TEMP_ROOT, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(TMPDIR = TEMP_ROOT, TMP = TEMP_ROOT, TEMP = TEMP_ROOT)
source(file.path(PIPELINE_REPO, "R", "source.R"))
set.seed(20260814L)'''),
        ("markdown", """### Key assumptions

All four regions use one installed Xenium panel. Mouse is the biological replicate; Regions 1–2 are sections from Mouse 1 and Regions 3–4 are sections from Mouse 2. No left/right or treatment factor is defined. This notebook performs technical QC only. The primary cell cohort is prespecified as `5 < nFeature_Xenium < 200` and `10 < nCount_Xenium < 1000`. Segmentation, nucleus, area, control, and spatial findings remain separate review flags."""),
        ("code", '''runtime <- validate_runtime_paths(project_root = PROJECT_ROOT, input_root = INPUT_ROOT, output_root = RUN_ROOT, temp_root = TEMP_ROOT)
regions <- discover_xenium_sections(INPUT_ROOT, expected_section_count = 4L)
stopifnot(identical(regions$region_id, expected_scwat_regions()))
manifest <- utils::read.delim(file.path(PIPELINE_REPO, "config", "scwat_sample_manifest.tsv"), check.names = FALSE)
manifest_check <- validate_sample_manifest(manifest, expected_scwat_regions())
stopifnot(isTRUE(manifest_check$valid))
fixed_thresholds <- read_fixed_cell_qc_thresholds(file.path(PIPELINE_REPO, "config", "fixed_cell_qc_thresholds.tsv"))
manifest[manifest$region_id == REGION_ID, , drop = FALSE]'''),
        ("markdown", """## Inputs and integrity

Discover the region directory from its identifier, inventory required files, and verify matrix/metadata alignment before loading counts. Raw data are read-only. Full transcript Parquet is not collected into memory; later full-HPC transcript summaries must use Arrow-backed aggregation."""),
        ("code", '''region_record <- discover_one_section(INPUT_ROOT, REGION_ID)
region_dir <- region_record$region_dir[[1L]]
inventory <- inventory_section_files(region_dir, REGION_ID, calculate_md5 = FALSE)
integrity <- validate_section_integrity(region_dir, REGION_ID)
stopifnot(all(inventory$exists), isTRUE(integrity$dimension_match[[1L]]))
bundle <- import_xenium_mex(region_dir)
stopifnot(inherits(bundle$counts, "sparseMatrix"), identical(colnames(bundle$counts), bundle$cells$cell_id))'''),
        ("markdown", """## Cell QC

Compute targeted-panel complexity, control burden, and segmentation review fields without deleting cells. Exact boundary values 5/200 features and 10/1000 counts fail. `primary_include` contains only the fixed rule; `strict_include` additionally excludes all review flags."""),
        ("code", '''cell_qc <- calculate_xenium_cell_qc(bundle$counts, bundle$cells, REGION_ID, fixed_thresholds)
masks <- build_cell_downstream_masks(cell_qc$cell_metadata, provenance = paste(RUN_LABEL, EXECUTION_MODE, sep = "::"))
stopifnot(nrow(masks) == ncol(bundle$counts), all(!masks$strict_include | masks$primary_include))
table(primary_include = masks$primary_include, strict_include = masks$strict_include)'''),
        ("markdown", """## Outputs and checks

Write an auditable region bundle beneath the run root: sparse raw counts, all cell metadata and masks, thresholds, inventory/integrity, run configuration, readiness evidence, and session information. The source notebook stays output-free; executed copies belong under the run directory."""),
        ("code", '''section_output_dir <- file.path(RUN_ROOT, "sections", REGION_ID)
result <- write_scwat_region_qc_bundle(
  project_root = PROJECT_ROOT, output_dir = section_output_dir, region_id = REGION_ID,
  run_label = RUN_LABEL, execution_mode = EXECUTION_MODE, manifest = manifest,
  inventory = inventory, integrity = integrity, bundle = bundle, cell_qc = cell_qc, masks = masks
)
stopifnot(validate_scwat_region_qc_bundle(section_output_dir, REGION_ID))
result'''),
        ("markdown", """## Next steps

Review region tables and figures, especially assignment/control burden, boundary counts, segmentation flags, and spatial patterns. A `PASS` here is technical evidence only. After all four regions finish under the same run label and mode, execute `02_slide_QC_summary.ipynb`. Full-data HPC output—not a local subset—determines final readiness before the separate Region 3 anchor workflow."""),
    ]


def summary_cells() -> list[tuple[str, str]]:
    return [
        ("markdown", """# scWAT Xenium slide QC summary

## Goal

Combine exactly four validated region bundles from one run label and execution mode. The summary compares two sections within each of `Mouse_1` and `Mouse_2` without treating cells as biological replicates and without inventing a left/right factor."""),
        ("markdown", """## Setup

HPC paths are defaults. Local validation injects D-drive paths and `LOCAL_SUBSET`. No package installation or writes outside `adipose_analysis` are permitted."""),
        ("code", '''EXECUTION_MODE <- toupper(Sys.getenv("SCWAT_QC_MODE", "FULL_HPC"))
PROJECT_ROOT <- Sys.getenv("SCWAT_PROJECT_ROOT", "/dssg/home/acct-svetoslav_chakarov/svetoslav_chakarov/Lab_members/Yanan_Hu/YNH_Xenium")
PIPELINE_REPO <- Sys.getenv("SCWAT_PIPELINE_REPO", file.path(PROJECT_ROOT, "adipose_analysis", "YNH_Xenium_scWAT"))
RUN_LABEL <- Sys.getenv("SCWAT_RUN_LABEL", "scwat_qc_hpc")
RUN_ROOT <- Sys.getenv("SCWAT_RUN_ROOT", file.path(PROJECT_ROOT, "adipose_analysis", "scwat_qc_outputs", RUN_LABEL))
TEMP_ROOT <- Sys.getenv("SCWAT_TEMP_ROOT", file.path(PROJECT_ROOT, "adipose_analysis", "tmp", RUN_LABEL))
dir.create(TEMP_ROOT, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(TMPDIR = TEMP_ROOT, TMP = TEMP_ROOT, TEMP = TEMP_ROOT)
source(file.path(PIPELINE_REPO, "R", "source.R"))
regions <- expected_scwat_regions()'''),
        ("markdown", """## Inputs and integrity

Require Region_1 through Region_4, reject missing or duplicate regions and mixed run labels or modes, then reload every saved sparse object and table."""),
        ("code", '''manifest <- utils::read.delim(file.path(PIPELINE_REPO, "config", "scwat_sample_manifest.tsv"), check.names = FALSE)
manifest_check <- validate_sample_manifest(manifest, regions)
stopifnot(isTRUE(manifest_check$valid))
slide_data <- read_scwat_slide_qc_outputs(RUN_ROOT, expected_regions = regions, run_label = RUN_LABEL, execution_mode = EXECUTION_MODE)
stopifnot(identical(slide_data$coverage$region_id, regions))'''),
        ("markdown", """## Cell QC

Summarize fixed-bound pass fractions and review flags. `primary_include` is never redefined at slide level. Region readiness is derived from current evidence rather than inherited anchor or sensitivity labels."""),
        ("code", '''slide_summary <- summarise_slide_qc(slide_data)
mouse_summary <- summarise_scwat_mouse_sections(slide_summary$section_summary, manifest)
mouse_summary$within_mouse
mouse_summary$mouse_summary'''),
        ("markdown", """## Outputs and checks

Save combined tables, descriptive within-mouse and mouse-level summaries, readiness gates, figures, session information, and a reloadable slide object."""),
        ("code", '''summary_output_dir <- file.path(RUN_ROOT, "slide_summary")
result <- write_scwat_slide_qc_bundle(PROJECT_ROOT, summary_output_dir, RUN_LABEL, EXECUTION_MODE, slide_data, slide_summary, mouse_summary)
stopifnot(validate_scwat_slide_qc_bundle(summary_output_dir, expected_regions = regions))
result'''),
        ("markdown", """## Next steps

Use full-HPC evidence for the final initial-QC decision. With two mice, mouse and section comparisons are descriptive. Only after these gates pass should Region 3 be built as the anchor, Regions 1-2 be tested for admission, and Region 4 be handled as mapping-only sensitivity data."""),
    ]


def write_notebook(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")


def main() -> None:
    for index in range(1, 5):
        write_notebook(NOTEBOOK_DIR / f"01_QC_Region{index}.ipynb", notebook(region_cells(f"Region_{index}")))
    write_notebook(NOTEBOOK_DIR / "02_slide_QC_summary.ipynb", notebook(summary_cells()))


if __name__ == "__main__":
    main()
