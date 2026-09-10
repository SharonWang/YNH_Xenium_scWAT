options(stringsAsFactors = FALSE)

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (!length(script_arg)) stop("Run this test with Rscript.", call. = FALSE)
test_path <- normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = TRUE)
repo_root <- dirname(dirname(test_path))

source_path <- file.path(repo_root, "R", "source.R")
backup_path <- file.path(repo_root, "R", "source_bk.R")
stopifnot(file.exists(source_path), file.exists(backup_path))

active_environment <- new.env(parent = baseenv())
sys.source(source_path, envir = active_environment)
required_active_functions <- c(
  "read_fixed_cell_qc_thresholds", "apply_fixed_primary_bounds",
  "calculate_xenium_cell_qc", "build_cell_downstream_masks",
  "expected_scwat_regions", "write_scwat_region_qc_bundle",
  "validate_scwat_region_qc_bundle", "read_scwat_slide_qc_outputs",
  "summarise_slide_qc", "summarise_scwat_mouse_sections",
  "write_scwat_slide_qc_bundle", "validate_scwat_slide_qc_bundle",
  "region_bundle_to_spatial_seurat", "refine_eosinophil_identity",
  "score_eosinophil_likeness", "plot_eos_with_celltypes",
  "validate_runtime_paths", "derive_primary_include_revised",
  "build_tissue_branch_manifest", "select_tissue_branch_ids",
  "validate_disjoint_spatial_pools", "validate_neighbour_distances",
  "adjusted_rand_index", "summarise_cluster_stability",
  "summarise_pca_qc_correlations", "derive_lymph_node_domain",
  "summarise_ln_boundary_sensitivity", "write_validated_seurat_checkpoint"
)
stopifnot(all(vapply(
  required_active_functions,
  exists,
  logical(1),
  envir = active_environment,
  mode = "function",
  inherits = FALSE
)))

archived_functions <- c(
  "fit_region3_anchor_branches", "map_and_evaluate_conditional_region",
  "build_eligible_consensus", "map_region4_sensitivity",
  "finalize_downstream_release", "write_section_artifacts",
  "write_extended_section_artifacts", "write_evidence_only_qc_artifacts"
)
stopifnot(!any(vapply(
  archived_functions,
  exists,
  logical(1),
  envir = active_environment,
  mode = "function",
  inherits = FALSE
)))

archive_environment <- new.env(parent = active_environment)
sys.source(backup_path, envir = archive_environment)
stopifnot(all(vapply(
  archived_functions,
  exists,
  logical(1),
  envir = archive_environment,
  mode = "function",
  inherits = FALSE
)))

# Every active public/helper function must carry an immediately adjacent
# roxygen block that states its return contract. This catches regression to the
# former generic Purpose/Inputs/Output comments.
source_lines <- readLines(source_path, warn = FALSE)
definition_lines <- grep("^[A-Za-z][A-Za-z0-9._]*[[:space:]]*<-[[:space:]]*function", source_lines)
documentation_ok <- vapply(definition_lines, function(line_number) {
  cursor <- line_number - 1L
  while (cursor > 0L && !nzchar(trimws(source_lines[[cursor]]))) cursor <- cursor - 1L
  end <- cursor
  while (cursor > 0L && grepl("^[[:space:]]*#'", source_lines[[cursor]])) cursor <- cursor - 1L
  block <- if (end >= cursor + 1L) source_lines[(cursor + 1L):end] else character()
  length(block) > 0L && any(grepl("@return", block, fixed = TRUE))
}, logical(1))
stopifnot(all(documentation_ok))

cat("Active and archival source contracts passed.\n")
