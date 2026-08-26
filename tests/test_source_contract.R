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
  "write_section_artifacts", "summarise_slide_qc",
  "region_bundle_to_spatial_seurat", "refine_eosinophil_identity",
  "score_eosinophil_likeness", "plot_eos_with_celltypes",
  "summarise_transcript_quality_table", "validate_runtime_paths"
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
  "finalize_downstream_release"
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

cat("Active and archival source contracts passed.\n")
