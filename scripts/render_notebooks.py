#!/usr/bin/env python3
"""Build, inject parameters into, and structurally validate the scWAT R notebooks."""

import argparse
import json
import re
from pathlib import Path


REGION_NOTEBOOKS = {
    "Region_1": "01_section_phase0_2_QC_Region1.ipynb",
    "Region_2": "01_section_phase0_2_QC_Region2.ipynb",
    "Region_3": "01_section_phase0_2_QC_Region3.ipynb",
    "Region_4": "01_section_phase0_2_QC_Region4.ipynb",
}


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
        markdown("# scWAT Xenium section QC: Phases 0-2\n\nRun this notebook once per section. It retains every cell, preserves raw sparse counts, and adds advisory alarm/gene/spatial diagnostics."),
        markdown("## Goal\n\nValidate one Xenium section, reconcile its panel, import sparse counts, calculate section-specific QC flags, and write reload-validated core and extended artifact bundles."),
        markdown("## Setup\n\n### Parameters\n\nChange `REGION_ID` for manual execution. Launchers inject the same parameters without editing the source notebook."),
        code(
            'PROJECT_ROOT <- "/dssg/home/acct-svetoslav_chakarov/svetoslav_chakarov/Lab_members/Yanan_Hu/YNH_Xenium"\n'
            'PIPELINE_REPO <- file.path(PROJECT_ROOT, "adipose_analysis", "YNH_Xenium_scWAT")\n'
            'INPUT_ROOT <- file.path(PROJECT_ROOT, "adipose_data")\n'
            'REGION_ID <- "Region_1"\n'
            'RUN_LABEL <- "full_notebook_qc_v2"\n'
            'METADATA_PATH <- file.path(PIPELINE_REPO, "config", "scwat_sample_manifest.tsv")\n'
            'EXPECTED_SECTION_COUNT <- 4L\n'
            'SEED <- 20260814L\n'
            'STRICT_MODE <- FALSE\n'
            'EXTENDED_QC_MODE <- "AUTO"\n'
            'EXTENDED_QC_CONFIG_PATH <- file.path(PIPELINE_REPO, "config", "extended_qc_defaults.tsv")\n',
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
            'metadata_columns <- intersect(c("mouse_id", "side", "section_id", "biological_replicate_id", "genotype", "treatment", "condition", "age_weeks", "metadata_status", "do_not_interpret"), names(section_manifest))\n'
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
        markdown("## Extended QC preflight and direct alarm evidence\n\nInputs are the section directory, versioned extended-QC settings, and existing Xenium alarm records. The exact affected cycle remains unresolved without 10x diagnostics."),
        code(
            'extended_config <- read_extended_qc_config(EXTENDED_QC_CONFIG_PATH)\n'
            'extended_mode <- resolve_extended_qc_mode(EXTENDED_QC_MODE, region_dir)\n'
            'extended_preflight <- extended_qc_preflight(extended_mode, region_dir, extended_config)\n'
            'if (extended_mode == "FULL_HPC" && any(extended_preflight$status == "FAIL")) {\n'
            '  stop(paste("FULL_HPC extended preflight failed:", paste(extended_preflight$check[extended_preflight$status == "FAIL"], collapse = ", ")))\n'
            '}\n'
            'cycle_alarm_evidence <- build_cycle_alarm_evidence(alarms, REGION_ID)\n'
            'cat("Extended mode:", extended_mode, "\\nInput:", region_dir, "\\nOutput:", SECTION_OUTPUT_DIR, "\\n")\n'
            'extended_preflight\n'
            'cycle_alarm_evidence\n'
        ),
        markdown("## Extended per-gene quality diagnostics\n\nMatrix metrics are calculated for every panel feature. Full-HPC mode lazily aggregates `transcripts.parquet`; local subset mode records `NOT_RUN_LOCAL_SUBSET` and does not imply acceptable transcript quality."),
        code(
            'matrix_gene_quality <- summarise_gene_matrix_qc(xenium$counts, REGION_ID, gene_sets)\n'
            'if (extended_mode == "FULL_HPC") {\n'
            '  transcript_gene_quality <- summarise_transcript_quality_arrow(file.path(region_dir, "transcripts.parquet"), REGION_ID, extended_config$qv_threshold)\n'
            '  gene_quality <- combine_gene_quality(matrix_gene_quality, transcript_gene_quality)\n'
            '} else {\n'
            '  gene_quality <- matrix_gene_quality\n'
            '  gene_quality$transcript_rows <- NA_integer_\n'
            '  gene_quality$mean_qv <- NA_real_\n'
            '  gene_quality$fraction_q20 <- NA_real_\n'
            '  gene_quality$represented_codewords <- NA_integer_\n'
            '  gene_quality$transcript_status <- "NOT_RUN_LOCAL_SUBSET"\n'
            '}\n'
            'gene_quality[seq_len(min(12L, nrow(gene_quality))), , drop = FALSE]\n'
        ),
        markdown("## Extended spatial diagnostics\n\nCoordinates test global clustering, an occupied-grid tissue-edge proxy, a kNN dense-aggregate proxy, and candidate hotspot bins. These labels are coordinate evidence only: folds, tears, and other morphology require image review."),
        code(
            'spatial_cells <- assign_spatial_grid(qc$cell_metadata, extended_config$grid_size_um)\n'
            'spatial_k <- min(as.integer(extended_config$spatial_k), nrow(spatial_cells) - 1L)\n'
            'spatial_cells$local_density <- calculate_knn_density(spatial_cells, spatial_k, extended_mode)\n'
            'density_cutoff <- stats::quantile(spatial_cells$local_density, extended_config$dense_quantile, na.rm = TRUE)\n'
            'spatial_cells$dense_aggregate <- spatial_cells$local_density >= density_cutoff\n'
            'spatial_global <- test_spatial_flag_clustering(spatial_cells, spatial_k, extended_config$permutations, extended_config$seed, extended_mode)\n'
            'spatial_global$region_id <- REGION_ID\n'
            'spatial_edge_density <- summarise_spatial_enrichment(spatial_cells)\n'
            'spatial_edge_density$region_id <- REGION_ID\n'
            'spatial_hotspots <- find_spatial_qc_hotspots(spatial_cells, extended_config$permutations, extended_config$min_bin_cells, extended_config$hotspot_fdr, extended_config$seed)\n'
            'spatial_hotspots$region_id <- rep(REGION_ID, nrow(spatial_hotspots))\n'
            'manual_review_manifest <- spatial_hotspots[spatial_hotspots$hotspot_status == "MORPHOLOGY_REVIEW_REQUIRED", , drop = FALSE]\n'
            'if (!nrow(manual_review_manifest)) manual_review_manifest$review_note <- character() else manual_review_manifest$review_note <- "Inspect morphology/image for edge, fold, tear, or dense aggregate context"\n'
            'list(global = spatial_global, enrichment = spatial_edge_density, candidate_hotspots = manual_review_manifest)\n'
        ),
        markdown("## Extended spatial figures"),
        code(
            'extended_plots <- plot_extended_spatial_qc(spatial_cells, spatial_edge_density, spatial_hotspots, REGION_ID)\n'
            'for (plot in extended_plots) print(plot)\n'
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
        code(
            'extended_artifact_paths <- write_extended_section_artifacts(\n'
            '  project_root = PROJECT_ROOT, output_dir = SECTION_OUTPUT_DIR, region_id = REGION_ID, mode = extended_mode,\n'
            '  preflight = extended_preflight, cycle_alarm_evidence = cycle_alarm_evidence, gene_quality = gene_quality,\n'
            '  spatial_global = spatial_global, spatial_edge_density = spatial_edge_density, spatial_hotspots = spatial_hotspots,\n'
            '  spatial_cells = spatial_cells, manual_review_manifest = manual_review_manifest, plots = extended_plots\n'
            ')\n'
            'stopifnot(validate_extended_section_artifacts(SECTION_OUTPUT_DIR, REGION_ID, extended_mode))\n'
            'extended_reload <- read_extended_section_artifacts(SECTION_OUTPUT_DIR, REGION_ID, extended_mode)\n'
            'stopifnot(nrow(extended_reload$spatial_cells) == nrow(qc$cell_metadata))\n'
            'extended_reload$status\n'
        ),
        markdown("## Outputs\n\nAll outputs are section-specific. A `HOLD` or `PENDING` gate blocks biological interpretation but preserves diagnostic QC."),
        code(
            'data.frame(artifact = basename(artifact_paths), path = artifact_paths)[seq_len(min(length(artifact_paths), 20L)), , drop = FALSE]\n'
            'data.frame(artifact = basename(extended_artifact_paths), path = extended_artifact_paths)\n'
            'cat("Completed", REGION_ID, "with", nrow(qc$cell_metadata), "cells; zero cells deleted. Extended mode:", extended_mode, "\\n")\n'
        ),
    ]
    if incomplete:
        cells = [cell for cell in cells if not (cell["cell_type"] == "markdown" and "## Checks" in "".join(cell["source"]))]
    return cells


def summary_cells(incomplete=False):
    cells = [
        markdown("# scWAT Xenium slide-level QC summary\n\nThis notebook is the sole reader-facing QC report for the four independently processed scWAT sections. Sections are technical units; the two mice are the biological units."),
        markdown("## Setup\n\n### Parameters\n\nInputs: four completed section bundles below `RUN_ROOT`, the verified sample manifest, and versioned QC settings. Outputs: combined tables and Cell-style figures below `${RUN_ROOT}/slide_summary/`."),
        code(
            'PROJECT_ROOT <- "/dssg/home/acct-svetoslav_chakarov/svetoslav_chakarov/Lab_members/Yanan_Hu/YNH_Xenium"\n'
            'PIPELINE_REPO <- file.path(PROJECT_ROOT, "adipose_analysis", "YNH_Xenium_scWAT")\n'
            'RUN_LABEL <- "full_notebook_qc_v2"\n'
            'EXPECTED_SECTION_COUNT <- 4L\n'
            'METADATA_PATH <- file.path(PIPELINE_REPO, "config", "scwat_sample_manifest.tsv")\n'
            'EXTENDED_QC_CONFIG_PATH <- file.path(PIPELINE_REPO, "config", "extended_qc_defaults.tsv")\n'
            'SUBSET_REFERENCE_PATH <- file.path(PIPELINE_REPO, "config", "subset_qc_reference.tsv")\n',
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
        markdown("## Inputs and validation\n\nExactly four independently completed core and extended section bundles are required. The verified manifest maps 62308/62309 to Mouse 1 and 62310/62311 to Mouse 2. Left/right is not an analysis factor."),
        code(
            'coverage <- validate_four_section_outputs(RUN_ROOT, paste0("Region_", seq_len(EXPECTED_SECTION_COUNT)))\n'
            'extended_coverage <- validate_four_extended_section_outputs(RUN_ROOT, coverage$region_id)\n'
            'stopifnot(identical(coverage$region_id, extended_coverage$region_id))\n'
            'slide_data <- read_slide_qc_outputs(RUN_ROOT, coverage$region_id)\n'
            'slide_summary <- summarise_slide_qc(slide_data)\n'
            'extended_slide_data <- read_extended_slide_qc_outputs(RUN_ROOT, coverage$region_id)\n'
            'manifest <- utils::read.delim(METADATA_PATH, check.names = FALSE)\n'
            'extended_config <- read_extended_qc_config(EXTENDED_QC_CONFIG_PATH)\n'
            'subset_reference <- utils::read.delim(SUBSET_REFERENCE_PATH, check.names = FALSE)\n'
            'extended_slide_summary <- summarise_extended_slide_qc(extended_slide_data, slide_summary$section_summary, manifest, extended_config, subset_reference)\n'
            'stopifnot(length(unique(slide_data$cell_metadata$region_id)) == 4L)\n'
            'list(core = coverage, extended = extended_coverage, cells = nrow(slide_data$cell_metadata))\n'
        ),
        markdown("## TL;DR and QC decision\n\nThe table below is the review entry point. Direct poor-cycle alarms are available for Regions 1, 2, and 4. Exact cycle/channel/codeword identity is `REQUIRES_10X_DIAGNOSTICS`. Region 3 requires morphology review; Region 4 additionally requires segmentation/cell-area review."),
        code(
            'qc_decision <- merge(slide_summary$section_summary, slide_summary$readiness[, c("region_id", "status")], by = "region_id", all.x = TRUE, sort = FALSE)\n'
            'qc_decision <- qc_decision[match(coverage$region_id, qc_decision$region_id), , drop = FALSE]\n'
            'direct_alarm_regions <- unique(extended_slide_data$cycle_alarm_evidence$region_id[extended_slide_data$cycle_alarm_evidence$evidence_status == "DIRECT_EVIDENCE"])\n'
            'qc_decision$direct_poor_cycle_alarm <- qc_decision$region_id %in% direct_alarm_regions\n'
            'qc_decision$reported_qc_status <- ifelse(qc_decision$region_id == "Region_3" & qc_decision$status == "PASS", "CONDITIONAL_PASS", qc_decision$status)\n'
            'qc_decision$required_next_action <- c("Obtain 10x poor-cycle diagnostics", "Obtain 10x poor-cycle diagnostics", "Review morphology in FDR-positive hotspot bins", "Obtain 10x diagnostics and review segmentation/cell area")[match(qc_decision$region_id, paste0("Region_", 1:4))]\n'
            'qc_decision[, c("region_id", "input_cells", "core_pass_fraction", "review_fraction", "direct_poor_cycle_alarm", "reported_qc_status", "required_next_action")]\n'
        ),
        markdown("## Core QC distributions\n\nThese are descriptive section-level and cell-level QC summaries. No cells are automatically deleted, and cell-level distributions do not create biological replication."),
        code(
            'slide_summary$section_summary\n'
            'stopifnot(all(slide_summary$section_summary$cells_deleted == 0L))\n'
            'slide_plots <- plot_slide_qc(slide_data, slide_summary)\n'
            'for (plot in slide_plots) print(plot)\n'
        ),
        markdown("## Alarm evidence and candidate genes\n\nThe alarm table reports directly available evidence. The candidate list uses cross-section transcript depletion/quality patterns and remains `CANDIDATE_NOT_CONFIRMED`; it is not a definitive cycle-to-gene map. Exact cycle identity is `REQUIRES_10X_DIAGNOSTICS` and must come from 10x diagnostic output."),
        code(
            'extended_slide_data$cycle_alarm_evidence\n'
            'candidate_display <- extended_slide_summary$candidates[extended_slide_summary$candidates$section_candidate_flag, , drop = FALSE]\n'
            'gene_sets <- unique(extended_slide_data$gene_quality[, c("gene", "gene_set")])\n'
            'candidate_display <- merge(candidate_display, gene_sets, by = "gene", all.x = TRUE, sort = FALSE)\n'
            'candidate_counts <- if (nrow(candidate_display)) aggregate(section_candidate_flag ~ region_id + section_evidence_status, candidate_display, sum) else data.frame(region_id = character(), section_evidence_status = character(), section_candidate_flag = integer())\n'
            'candidate_counts\n'
            'if (nrow(candidate_display)) with(candidate_display, table(region_id, gene_set, useNA = "ifany"))\n'
            'candidate_display <- candidate_display[order(candidate_display$region_id, candidate_display$evidence_tier, candidate_display$log2_count_ratio, candidate_display$q20_difference), , drop = FALSE]\n'
            'candidate_display[seq_len(min(30L, nrow(candidate_display))), , drop = FALSE]\n'
        ),
        markdown("## Subset versus full-data burden\n\nThe comparison is descriptive across four technical sections. `NOT_RUN_LOCAL_SUBSET` means the full-data ranking remains an HPC checkpoint; rank correlations across only four sections are not biological evidence."),
        code(
            'extended_slide_summary$ranking\n'
            'extended_slide_summary$rank_agreement\n'
        ),
        markdown("## Spatial QC\n\nGlobal kNN clustering, tissue-edge proxies, dense-cell proxies, and candidate hotspot bins are coordinate-based diagnostics. A hotspot is only `MORPHOLOGY_REVIEW_REQUIRED`; folds, tears, tissue edges, and aggregates require image review."),
        code(
            'extended_slide_data$spatial_global\n'
            'extended_slide_data$spatial_edge_density\n'
            'extended_slide_data$manual_review_manifest\n'
        ),
        markdown("## Within-mouse concordance\n\nThe verified technical pairs are 62308/62309 for Mouse 1 and 62310/62311 for Mouse 2. Thresholds are advisory. The Region 3/4 cell-area contrast is reviewed separately because it was not part of the original concordance gate."),
        code(
            'extended_slide_summary$concordance$summary\n'
        ),
        markdown("## Diagnostic Cell-style figures\n\nColors and scales are consistent across sections where scientifically appropriate. These plots support review rather than biological inference."),
        code(
            'extended_slide_plots <- plot_extended_slide_qc(extended_slide_data, extended_slide_summary)\n'
            'for (plot in extended_slide_plots) print(plot)\n'
        ),
        markdown("## Final QC decision and next actions\n\nRegions 1 and 2 remain `HOLD` pending 10x diagnostics. Region 3 is a conditional pass pending morphology review. Region 4 remains `HOLD` pending 10x diagnostics plus segmentation/cell-area review. Comparative gene-level analysis remains on hold until the affected-gene decision is frozen."),
        code(
            'qc_decision[, c("region_id", "reported_qc_status", "required_next_action")]\n'
            'extended_slide_summary$status\n'
            'cat("Overall slide QC status:", slide_summary$overall_status, "\\n")\n'
        ),
        markdown("## Outputs and reload checks\n\nAll tables, figures, and serialized objects are written below `${RUN_ROOT}/slide_summary/`. The notebook then reloads and validates the complete artifact contract."),
        code(
            'slide_artifacts <- write_slide_qc_artifacts(PROJECT_ROOT, RUN_ROOT, slide_data, slide_summary, slide_plots)\n'
            'extended_slide_artifacts <- write_extended_slide_qc_artifacts(PROJECT_ROOT, RUN_ROOT, extended_slide_data, extended_slide_summary, extended_slide_plots)\n'
            'saved_summary <- readRDS(file.path(RUN_ROOT, "slide_summary", "slide_qc_summary.rds"))\n'
            'stopifnot(nrow(saved_summary$data$coverage) == 4L)\n'
            'stopifnot(length(unique(saved_summary$data$cell_metadata$region_id)) == 4L)\n'
            'stopifnot(validate_extended_slide_qc_artifacts(RUN_ROOT, stop_on_error = TRUE))\n'
            'list(core = data.frame(artifact = basename(slide_artifacts), path = slide_artifacts),\n'
            '     extended = data.frame(artifact = basename(extended_slide_artifacts), path = extended_slide_artifacts))\n'
        ),
    ]
    if incomplete:
        cells = [cell for cell in cells if not (cell["cell_type"] == "markdown" and "## Final QC decision" in "".join(cell["source"]))]
    return cells


def validate_notebook(path, notebook_type="section", expected_region_id=None):
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    errors = []
    if data.get("nbformat") != 4 or not isinstance(data.get("cells"), list): errors.append("invalid nbformat structure")
    if data.get("metadata", {}).get("kernelspec", {}).get("name") != "ir": errors.append("R kernelspec is required")
    text = "\n".join("".join(cell.get("source", [])) for cell in data.get("cells", []))
    if re.search(r"(?i)(?:^|[\"'])C:[/\\]", text, flags=re.MULTILINE): errors.append("notebook contains a C: path")
    parameter_cells = [c for c in data.get("cells", []) if "parameters" in c.get("metadata", {}).get("tags", [])]
    if len(parameter_cells) != 1: errors.append("exactly one tagged parameters cell is required")
    if notebook_type == "section":
        required_parameters = ["PROJECT_ROOT", "PIPELINE_REPO", "INPUT_ROOT", "REGION_ID", "RUN_LABEL", "METADATA_PATH", "EXPECTED_SECTION_COUNT", "SEED", "STRICT_MODE", "EXTENDED_QC_MODE", "EXTENDED_QC_CONFIG_PATH"]
        required_sections = ["## Goal", "## Setup", "## Inputs", "## Phase 0", "## Phase 1", "## Phase 2", "## Extended QC preflight", "## Extended per-gene", "## Extended spatial diagnostics", "## Checks", "## Outputs"]
        for value in required_parameters + required_sections:
            if value not in text: errors.append(f"missing required section/parameter: {value}")
        if expected_region_id is not None:
            parameter_text = "".join(parameter_cells[0].get("source", [])) if len(parameter_cells) == 1 else ""
            match = re.search(r'^REGION_ID\s*<-\s*"(Region_[1-4])"$', parameter_text, flags=re.MULTILINE)
            if match is None or match.group(1) != expected_region_id:
                errors.append(f"expected fixed REGION_ID {expected_region_id}")
    if notebook_type == "summary":
        required_parameters = ["PROJECT_ROOT", "PIPELINE_REPO", "RUN_LABEL", "EXPECTED_SECTION_COUNT", "METADATA_PATH", "EXTENDED_QC_CONFIG_PATH", "SUBSET_REFERENCE_PATH"]
        required_sections = ["## Setup", "## Inputs and validation", "## TL;DR and QC decision", "## Core QC distributions", "## Alarm evidence and candidate genes", "## Subset versus full-data burden", "## Spatial QC", "## Within-mouse concordance", "## Diagnostic Cell-style figures", "## Final QC decision and next actions", "## Outputs and reload checks"]
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


def build_region_notebooks(repo):
    source = repo / "notebooks" / "01_section_phase0_2_QC.ipynb"
    built = []
    for region_id, filename in REGION_NOTEBOOKS.items():
        output = repo / "notebooks" / filename
        inject_parameters(source, output, {"REGION_ID": region_id})
        validate_notebook(output, "section", expected_region_id=region_id)
        built.append(output)
    return built


def region_notebook_path(repo, region_id):
    if region_id not in REGION_NOTEBOOKS:
        raise ValueError(f"Unknown region: {region_id}")
    path = repo / "notebooks" / REGION_NOTEBOOKS[region_id]
    if not path.is_file():
        raise ValueError(f"Committed region notebook is missing: {path}")
    validate_notebook(path, "section", expected_region_id=region_id)
    return path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--build-section", action="store_true")
    parser.add_argument("--build-regions", action="store_true")
    parser.add_argument("--build-summary", action="store_true")
    parser.add_argument("--incomplete", action="store_true")
    parser.add_argument("--validate")
    parser.add_argument("--region-notebook")
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
    if args.build_regions:
        for path in build_region_notebooks(repo):
            print(path)
    if args.build_summary:
        summary_path.write_text(json.dumps(notebook(summary_cells(args.incomplete)), indent=1, ensure_ascii=False) + "\n", encoding="utf-8")
        print(summary_path)
    if args.inject:
        values = dict(item.split("=", 1) for item in args.set)
        inject_parameters(args.inject[0], args.inject[1], values)
    if args.validate:
        validate_notebook(args.validate, args.type)
        print(f"Notebook validation passed: {args.validate}")
    if args.region_notebook:
        try:
            print(region_notebook_path(repo, args.region_notebook))
        except ValueError as error:
            parser.error(str(error))


if __name__ == "__main__":
    main()
