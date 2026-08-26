options(stringsAsFactors = FALSE)

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (!length(script_arg)) stop("Run this maintenance script with Rscript.", call. = FALSE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = TRUE)
repo_root <- dirname(dirname(script_path))
source_path <- file.path(repo_root, "R", "source.R")
backup_path <- file.path(repo_root, "R", "source_bk.R")

# This explicit allow-list was derived from all code cells in the committed
# notebooks plus active R/Python/PowerShell/shell tests and support scripts.
# Functions listed here have no active consumer and are archived, not deleted.
archive_names <- c(
  "annotate_seurat_clusters",
  "apply_manual_admission_review",
  "assert_downstream_model_environment",
  "assess_mapped_marker_coherence",
  "assess_region3_morphology_review",
  "assign_spatial_mapping_folds",
  "build_eligible_consensus",
  "build_seurat_reference",
  "calibrate_mapping_thresholds",
  "calibrate_region3_mapping",
  "classify_mapping_uncertainty",
  "cluster_marker_evidence",
  "downstream_model_dependencies",
  "downstream_model_preflight",
  "evaluate_section_admission",
  "finalize_downstream_release",
  "find_primary_markers",
  "fit_region3_anchor_branches",
  "get_seurat_normalized_data",
  "identify_eosinophils_independently",
  "join_seurat_layers_if_needed",
  "label_agreement_on_major",
  "map_and_evaluate_conditional_region",
  "map_query_to_frozen_reference",
  "map_region4_sensitivity",
  "match_cluster_labels",
  "parse_downstream_cli",
  "plot_annotation_overlap_heatmap",
  "plot_spatial_discrete_overlay",
  "prediction_score_columns",
  "prepare_seurat_query",
  "primary_model_feature_policy",
  "read_canonical_marker_config",
  "read_downstream_handoff",
  "read_downstream_reference_config",
  "read_region3_anchor_artifacts",
  "reference_label_centroids",
  "region3_anchor_branch_definitions",
  "require_downstream_argument",
  "run_eosinophil_robustness",
  "score_standardized_gene_set",
  "validate_anchor_branches",
  "validate_downstream_handoff",
  "validate_stage_files",
  "write_consensus_artifacts",
  "write_mapping_calibration_artifacts",
  "write_region_admission_artifacts",
  "write_region3_anchor_artifacts"
)

source_lines <- readLines(source_path, warn = FALSE)
definition_pattern <- "^([A-Za-z][A-Za-z0-9._]*)[[:space:]]*<-[[:space:]]*function[[:space:]]*\\(.*$"
definition_lines <- grep(definition_pattern, source_lines)
definition_names <- sub(definition_pattern, "\\1", source_lines[definition_lines])
missing_archive <- setdiff(archive_names, definition_names)
if (length(missing_archive)) {
  stop(sprintf("Archive functions not found: %s", paste(missing_archive, collapse = ", ")), call. = FALSE)
}

definition_ends <- c(definition_lines[-1L] - 1L, length(source_lines))
archive_indexes <- match(archive_names, definition_names)
remove_line <- rep(FALSE, length(source_lines))
archive_blocks <- vector("list", length(archive_indexes))
for (block_index in seq_along(archive_indexes)) {
  definition_index <- archive_indexes[[block_index]]
  line_index <- definition_lines[[definition_index]]:definition_ends[[definition_index]]
  remove_line[line_index] <- TRUE
  archive_blocks[[block_index]] <- source_lines[line_index]
}

active_lines <- source_lines[!remove_line]
active_lines <- active_lines[!grepl("Scientifically revised downstream reference workflow|Single canonical implementation; superseded checkpoint definitions were removed", active_lines, fixed = FALSE)]
while (length(active_lines) && !nzchar(tail(active_lines, 1L))) active_lines <- head(active_lines, -1L)

backup_header <- c(
  "# ARCHIVAL SOURCE — NOT LOADED BY ACTIVE QC NOTEBOOKS",
  "#",
  "# Purpose: preserve functions with no active notebook, test, or support-script",
  "# consumer as of 2026-08-26. This file is retained for provenance and possible",
  "# later recovery; it must not be sourced by the Phase 0-2 QC notebooks.",
  "#",
  "# Inputs: source R/source.R first if an archived function is intentionally",
  "# reactivated, because archived workflows may call shared active helpers.",
  "# Output: function definitions only; sourcing this file writes no artifacts.",
  ""
)
archive_lines <- c(backup_header, unlist(Map(
  function(name, block) c(sprintf("# Archived function: %s", name), block, ""),
  archive_names,
  archive_blocks,
  USE.NAMES = FALSE
)))

writeLines(active_lines, source_path, useBytes = TRUE)
writeLines(archive_lines, backup_path, useBytes = TRUE)
cat(sprintf("Kept %d active functions; archived %d functions.\n", length(definition_names) - length(archive_names), length(archive_names)))
