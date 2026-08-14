#!/usr/bin/env python3
"""Build, inject parameters into, and structurally validate the scWAT R notebooks."""

import argparse
import json
import re
from pathlib import Path


def markdown(text):
    return {"cell_type": "markdown", "metadata": {}, "source": text.splitlines(True)}


def code(text, tags=None):
    metadata = {} if not tags else {"tags": tags}
    return {"cell_type": "code", "execution_count": None, "metadata": metadata, "outputs": [], "source": text.splitlines(True)}


def notebook(cells):
    return {
        "cells": cells,
        "metadata": {
            "kernelspec": {"display_name": "R", "language": "R", "name": "ir"},
            "language_info": {"codemirror_mode": "r", "file_extension": ".r", "mimetype": "text/x-r-source", "name": "R", "pygments_lexer": "r", "version": "4.3"},
        },
        "nbformat": 4,
        "nbformat_minor": 5,
    }


def section_cells(incomplete=False):
    cells = [
        markdown("# scWAT Xenium section QC: Phases 0-2\n\nRun this notebook once per section. It retains every cell and preserves raw sparse counts."),
        markdown("## Goal\n\nValidate one Xenium section, reconcile its panel, import sparse counts, calculate section-specific QC flags, and write a reload-validated artifact bundle."),
        markdown("## Setup\n\n### Parameters\n\nChange `REGION_ID` for manual execution. Launchers inject the same parameters without editing the source notebook."),
        code(
            'PROJECT_ROOT <- "/dssg/home/acct-svetoslav_chakarov/svetoslav_chakarov/Lab_members/Yanan_Hu"\n'
            'PIPELINE_REPO <- file.path(PROJECT_ROOT, "adipose_analysis", "YNH_Xenium_scWAT")\n'
            'INPUT_ROOT <- file.path(PROJECT_ROOT, "adipose_data")\n'
            'REGION_ID <- "Region_1"\n'
            'RUN_LABEL <- "full_notebook_qc_v1"\n'
            'METADATA_PATH <- ""\n'
            'EXPECTED_SECTION_COUNT <- 4L\n'
            'SEED <- 20260814L\n'
            'STRICT_MODE <- FALSE\n',
            tags=["parameters"],
        ),
        code(
            'OUTPUT_ROOT <- file.path(PROJECT_ROOT, "adipose_analysis", "scwat_qc_outputs", RUN_LABEL)\n'
            'SECTION_OUTPUT_DIR <- file.path(OUTPUT_ROOT, "sections", REGION_ID)\n'
            'source(file.path(PIPELINE_REPO, "R", "source.R"))\n'
            'for (package in c("Matrix", "jsonlite", "ggplot2")) require_package(package)\n'
            'validate_runtime_paths(PROJECT_ROOT, INPUT_ROOT, SECTION_OUTPUT_DIR, tempdir())\n'
            'dir.create(SECTION_OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)\n'
            'stopifnot(REGION_ID %in% paste0("Region_", seq_len(EXPECTED_SECTION_COUNT)))\n'
            'cat("Section:", REGION_ID, "\\nOutput:", SECTION_OUTPUT_DIR, "\\n")\n'
        ),
        markdown("## Inputs\n\nThe selected raw/subset section directory is discovered by region ID. All raw files remain read-only."),
        code(
            'all_sections <- discover_xenium_sections(INPUT_ROOT, EXPECTED_SECTION_COUNT)\n'
            'selected_section <- discover_one_section(INPUT_ROOT, REGION_ID)\n'
            'selected_section\n'
        ),
        markdown("## Phase 0 - Configuration and metadata contract"),
        code(
            'if (nzchar(METADATA_PATH)) {\n'
            '  assert_path_within(PROJECT_ROOT, METADATA_PATH)\n'
            '  full_manifest <- utils::read.delim(METADATA_PATH, check.names = FALSE)\n'
            '} else {\n'
            '  full_manifest <- create_synthetic_manifest(all_sections$region_id, SEED)\n'
            '}\n'
            'manifest_check <- validate_sample_manifest(full_manifest, all_sections$region_id)\n'
            'if (!manifest_check$valid) stop(paste(manifest_check$issues, collapse = "; "))\n'
            'section_manifest <- full_manifest[full_manifest$region_id == REGION_ID, , drop = FALSE]\n'
            'section_manifest$region_dir <- selected_section$region_dir\n'
            'configuration <- data.frame(\n'
            '  key = c("project_root", "pipeline_repo", "input_root", "output_root", "region_id", "run_label", "seed", "strict_mode"),\n'
            '  value = c(PROJECT_ROOT, PIPELINE_REPO, INPUT_ROOT, OUTPUT_ROOT, REGION_ID, RUN_LABEL, SEED, STRICT_MODE)\n'
            ')\n'
            'environment <- data.frame(\n'
            '  item = c("R_version", "platform", "tempdir", "Matrix", "jsonlite", "ggplot2"),\n'
            '  value = c(R.version.string, R.version$platform, tempdir(), as.character(packageVersion("Matrix")), as.character(packageVersion("jsonlite")), as.character(packageVersion("ggplot2")))\n'
            ')\n'
            'section_manifest\n'
        ),
        markdown("## Phase 1 - Provenance, integrity, panel reconciliation, and QC gate"),
        code(
            'region_dir <- selected_section$region_dir[[1]]\n'
            'inventory <- inventory_section_files(region_dir, REGION_ID, calculate_md5 = TRUE)\n'
            'if (!all(inventory$exists)) stop(paste("Missing required files:", paste(inventory$relative_path[!inventory$exists], collapse = ", ")))\n'
            'integrity <- validate_section_integrity(region_dir, REGION_ID)\n'
            'alarms <- extract_analysis_alarms(file.path(region_dir, "analysis_summary.html"))\n'
            'signature <- utils::read.delim(file.path(PIPELINE_REPO, "config", "eos_gene_sets.tsv"), check.names = FALSE)\n'
            'installed_genes <- read_custom_panel_genes(file.path(region_dir, "gene_panel.json"))\n'
            'gene_sets <- setNames(signature$gene_set, signature$gene)\n'
            'panel_reconciliation <- reconcile_panel(signature$gene, installed_genes, gene_sets)\n'
            'features_preview <- read_xenium_features(file.path(region_dir, "cell_feature_matrix", "features.tsv.gz"))\n'
            'feature_type_summary <- as.data.frame(table(features_preview$feature_type), stringsAsFactors = FALSE)\n'
            'names(feature_type_summary) <- c("feature_type", "n_features")\n'
            'list(integrity = integrity, alarms = alarms, panel = table(panel_reconciliation$status))\n'
        ),
        markdown("## Phase 2 - Sparse import and section-specific cell QC"),
        code(
            'xenium <- import_xenium_mex(region_dir)\n'
            'qc <- calculate_xenium_cell_qc(xenium$counts, xenium$cells, REGION_ID)\n'
            'metadata_columns <- intersect(c("mouse_id", "side", "section_id", "biological_replicate_id", "metadata_status", "do_not_interpret"), names(section_manifest))\n'
            'for (column in metadata_columns) qc$cell_metadata[[column]] <- section_manifest[[column]][[1]]\n'
            'qc$summary\n'
        ),
        code(
            'section_plots <- plot_section_qc(qc$cell_metadata, REGION_ID)\n'
            'print(section_plots$counts)\n'
            'print(section_plots$features)\n'
            'print(section_plots$area)\n'
            'print(section_plots$spatial)\n'
        ),
        markdown("## Checks\n\nWrite every required artifact, then reload the saved sparse object and verify dimensions and cell alignment."),
        code(
            'artifact_paths <- write_section_artifacts(\n'
            '  project_root = PROJECT_ROOT, output_dir = SECTION_OUTPUT_DIR, region_id = REGION_ID,\n'
            '  configuration = configuration, manifest = section_manifest, environment = environment,\n'
            '  inventory = inventory, integrity = integrity, feature_type_summary = feature_type_summary,\n'
            '  panel_reconciliation = panel_reconciliation, alarms = alarms, qc = qc,\n'
            '  counts = xenium$counts, features = xenium$features, strict_mode = STRICT_MODE\n'
            ')\n'
            'stopifnot(validate_section_artifacts(SECTION_OUTPUT_DIR, REGION_ID))\n'
            'readiness <- utils::read.delim(file.path(SECTION_OUTPUT_DIR, "section_readiness_gates.tsv"), check.names = FALSE)\n'
            'readiness\n'
        ),
        markdown("## Outputs\n\nAll outputs are section-specific. A `HOLD` or `PENDING` gate blocks biological interpretation but preserves diagnostic QC."),
        code(
            'data.frame(artifact = basename(artifact_paths), path = artifact_paths)[seq_len(min(length(artifact_paths), 20L)), , drop = FALSE]\n'
            'cat("Completed", REGION_ID, "with", nrow(qc$cell_metadata), "cells; zero cells deleted.\\n")\n'
        ),
    ]
    if incomplete:
        cells = [cell for cell in cells if not (cell["cell_type"] == "markdown" and "## Checks" in "".join(cell["source"]))]
    return cells


def summary_cells(incomplete=False):
    cells = [
        markdown("# scWAT Xenium slide-level QC summary\n\nAggregate the four independently processed sections without treating cells as biological replicates."),
        markdown("## Goal\n\nVerify complete Region 1-4 coverage, compare QC distributions and gates, and write Cell-inspired slide-level figures."),
        markdown("## Setup\n\n### Parameters"),
        code(
            'PROJECT_ROOT <- "/dssg/home/acct-svetoslav_chakarov/svetoslav_chakarov/Lab_members/Yanan_Hu"\n'
            'PIPELINE_REPO <- file.path(PROJECT_ROOT, "adipose_analysis", "YNH_Xenium_scWAT")\n'
            'RUN_LABEL <- "full_notebook_qc_v1"\n'
            'EXPECTED_SECTION_COUNT <- 4L\n',
            tags=["parameters"],
        ),
        code(
            'RUN_ROOT <- file.path(PROJECT_ROOT, "adipose_analysis", "scwat_qc_outputs", RUN_LABEL)\n'
            'source(file.path(PIPELINE_REPO, "R", "source.R"))\n'
            'for (package in c("Matrix", "jsonlite", "ggplot2")) require_package(package)\n'
            'assert_path_within(PROJECT_ROOT, RUN_ROOT)\n'
            'assert_path_within(PROJECT_ROOT, tempdir())\n'
            'stopifnot(EXPECTED_SECTION_COUNT == 4L)\n'
            'cat("Slide QC run root:", RUN_ROOT, "\\n")\n'
        ),
        markdown("## Inputs\n\nExactly four independently completed section bundles are required."),
        markdown("## Completeness Checks"),
        code(
            'coverage <- validate_four_section_outputs(RUN_ROOT, paste0("Region_", seq_len(EXPECTED_SECTION_COUNT)))\n'
            'coverage\n'
        ),
        markdown("## QC Results"),
        code(
            'slide_data <- read_slide_qc_outputs(RUN_ROOT, coverage$region_id)\n'
            'slide_summary <- summarise_slide_qc(slide_data)\n'
            'slide_summary$section_summary\n'
        ),
        markdown("## Cell-style Figures\n\nColors are fixed across sections; distributions are descriptive and do not imply cell-level biological replication."),
        code(
            'slide_plots <- plot_slide_qc(slide_data, slide_summary)\n'
            'for (plot in slide_plots) print(plot)\n'
        ),
        markdown("## Readiness\n\nThe worst section gate determines slide readiness. Synthetic metadata and unresolved imaging errors block biology."),
        code(
            'slide_summary$readiness\n'
            'cat("Overall slide QC status:", slide_summary$overall_status, "\\n")\n'
        ),
        markdown("## Outputs"),
        code(
            'slide_artifacts <- write_slide_qc_artifacts(PROJECT_ROOT, RUN_ROOT, slide_data, slide_summary, slide_plots)\n'
            'saved_summary <- readRDS(file.path(RUN_ROOT, "slide_summary", "slide_qc_summary.rds"))\n'
            'stopifnot(nrow(saved_summary$data$coverage) == 4L)\n'
            'stopifnot(length(unique(saved_summary$data$cell_metadata$region_id)) == 4L)\n'
            'data.frame(artifact = basename(slide_artifacts), path = slide_artifacts)\n'
        ),
    ]
    if incomplete:
        cells = [cell for cell in cells if not (cell["cell_type"] == "markdown" and "## Readiness" in "".join(cell["source"]))]
    return cells


def validate_notebook(path, notebook_type="section"):
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    errors = []
    if data.get("nbformat") != 4 or not isinstance(data.get("cells"), list): errors.append("invalid nbformat structure")
    if data.get("metadata", {}).get("kernelspec", {}).get("name") != "ir": errors.append("R kernelspec is required")
    text = "\n".join("".join(cell.get("source", [])) for cell in data.get("cells", []))
    if re.search(r"(?i)(?:^|[\"'])C:[/\\]", text, flags=re.MULTILINE): errors.append("notebook contains a C: path")
    parameter_cells = [c for c in data.get("cells", []) if "parameters" in c.get("metadata", {}).get("tags", [])]
    if len(parameter_cells) != 1: errors.append("exactly one tagged parameters cell is required")
    if notebook_type == "section":
        required_parameters = ["PROJECT_ROOT", "PIPELINE_REPO", "INPUT_ROOT", "REGION_ID", "RUN_LABEL", "METADATA_PATH", "EXPECTED_SECTION_COUNT", "SEED", "STRICT_MODE"]
        required_sections = ["## Goal", "## Setup", "## Inputs", "## Phase 0", "## Phase 1", "## Phase 2", "## Checks", "## Outputs"]
        for value in required_parameters + required_sections:
            if value not in text: errors.append(f"missing required section/parameter: {value}")
    if notebook_type == "summary":
        required_parameters = ["PROJECT_ROOT", "PIPELINE_REPO", "RUN_LABEL", "EXPECTED_SECTION_COUNT"]
        required_sections = ["## Goal", "## Setup", "## Inputs", "## Completeness Checks", "## QC Results", "## Cell-style Figures", "## Readiness", "## Outputs"]
        for value in required_parameters + required_sections:
            if value not in text: errors.append(f"missing required section/parameter: {value}")
    if errors: raise ValueError("; ".join(errors))
    return True


def inject_parameters(source_path, output_path, assignments):
    data = json.loads(Path(source_path).read_text(encoding="utf-8"))
    cells = [c for c in data["cells"] if "parameters" in c.get("metadata", {}).get("tags", [])]
    if len(cells) != 1: raise ValueError("Cannot inject parameters without exactly one tagged cell")
    text = "".join(cells[0]["source"])
    for key, value in assignments.items():
        if value in ("TRUE", "FALSE") or re.fullmatch(r"[0-9]+L?", value): rendered = value
        else: rendered = json.dumps(value)
        text, count = re.subn(rf"^{re.escape(key)}\s*<-.*$", f"{key} <- {rendered}", text, flags=re.MULTILINE)
        if count != 1: raise ValueError(f"Parameter not found exactly once: {key}")
    cells[0]["source"] = text.splitlines(True)
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    Path(output_path).write_text(json.dumps(data, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--build-section", action="store_true")
    parser.add_argument("--build-summary", action="store_true")
    parser.add_argument("--incomplete", action="store_true")
    parser.add_argument("--validate")
    parser.add_argument("--type", default="section", choices=["section", "summary"])
    parser.add_argument("--inject", nargs=2, metavar=("SOURCE", "OUTPUT"))
    parser.add_argument("--set", action="append", default=[])
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    section_path = repo / "notebooks" / "01_section_phase0_2_QC.ipynb"
    summary_path = repo / "notebooks" / "02_slide_QC_summary.ipynb"
    if args.build_section:
        section_path.write_text(json.dumps(notebook(section_cells(args.incomplete)), indent=1, ensure_ascii=False) + "\n", encoding="utf-8")
        print(section_path)
    if args.build_summary:
        summary_path.write_text(json.dumps(notebook(summary_cells(args.incomplete)), indent=1, ensure_ascii=False) + "\n", encoding="utf-8")
        print(summary_path)
    if args.inject:
        values = dict(item.split("=", 1) for item in args.set)
        inject_parameters(args.inject[0], args.inject[1], values)
    if args.validate:
        validate_notebook(args.validate, args.type)
        print(f"Notebook validation passed: {args.validate}")


if __name__ == "__main__":
    main()
