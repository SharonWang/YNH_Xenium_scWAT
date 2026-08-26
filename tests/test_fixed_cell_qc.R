options(stringsAsFactors = FALSE)

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (!length(script_arg)) stop("Run this test with Rscript.", call. = FALSE)
test_path <- normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = TRUE)
repo_root <- dirname(dirname(test_path))

source(file.path(repo_root, "R", "source.R"))

threshold_path <- file.path(repo_root, "config", "fixed_cell_qc_thresholds.tsv")
stopifnot(file.exists(threshold_path))
fixed_thresholds <- read_fixed_cell_qc_thresholds(threshold_path)

# Break caught: replacing strict inequalities with inclusive comparisons would
# incorrectly retain cells exactly on one of the four approved boundaries.
boundary_counts <- c(11, 11, 999, 999, 10, 1000)
boundary_features <- c(5, 6, 199, 200, 6, 6)
expected_core_pass <- c(FALSE, TRUE, TRUE, FALSE, FALSE, FALSE)
stopifnot(identical(
  apply_fixed_primary_bounds(boundary_features, boundary_counts, fixed_thresholds),
  expected_core_pass
))

require_package("Matrix")
count_matrix <- Matrix::Matrix(0, nrow = 200L, ncol = length(boundary_counts), sparse = TRUE)
for (cell_index in seq_along(boundary_counts)) {
  feature_count <- boundary_features[[cell_index]]
  count_matrix[seq_len(feature_count), cell_index] <- 1
  count_matrix[1L, cell_index] <- count_matrix[1L, cell_index] + boundary_counts[[cell_index]] - feature_count
}
colnames(count_matrix) <- paste0("boundary_cell_", seq_along(boundary_counts))
rownames(count_matrix) <- paste0("Gene_", seq_len(nrow(count_matrix)))

cell_metadata <- data.frame(
  cell_id = colnames(count_matrix),
  total_counts = boundary_counts,
  control_probe_counts = 0,
  genomic_control_counts = 0,
  control_codeword_counts = 0,
  cell_area = 100,
  nucleus_count = 1,
  stringsAsFactors = FALSE
)

qc <- calculate_xenium_cell_qc(
  count_matrix,
  cell_metadata,
  region_id = "Region_1",
  fixed_thresholds = fixed_thresholds
)
stopifnot(identical(qc$cell_metadata$qc_core_pass, expected_core_pass))
stopifnot(identical(
  qc$thresholds$method[match(c("nCount_Xenium", "nFeature_Xenium"), qc$thresholds$metric)],
  rep("fixed_exclusive_user_approved", 2L)
))

# Break caught: adding segmentation/high-control exclusions back into the
# primary mask would violate the approved Colon-consistent cohort definition.
mask_fixture <- data.frame(
  region_id = rep("Region_3", 4L),
  cell_id = paste0("mask_cell_", seq_len(4L)),
  qc_core_pass = c(TRUE, TRUE, FALSE, TRUE),
  segmentation_multiplet_flag = c(FALSE, TRUE, FALSE, FALSE),
  high_control_flag = c(FALSE, FALSE, FALSE, TRUE),
  qc_review_flag = c(FALSE, TRUE, TRUE, TRUE),
  stringsAsFactors = FALSE
)
masks <- build_cell_downstream_masks(mask_fixture, provenance = "fixed_boundary_unit_test")
stopifnot(identical(masks$primary_include, c(TRUE, TRUE, FALSE, TRUE)))
stopifnot(identical(masks$strict_include, c(TRUE, FALSE, FALSE, FALSE)))
stopifnot(all(!masks$strict_include | masks$primary_include))

section_decision <- build_one_section_downstream_decision(
  masks,
  run_label = "fixed_boundary_unit_test",
  execution_mode = "LOCAL_SUBSET",
  provenance = "fixed_boundary_unit_test",
  generated_utc = "2026-08-26 UTC"
)
stopifnot(identical(
  section_decision$primary_mask_rule,
  "nFeature_Xenium > 5 AND nFeature_Xenium < 200 AND nCount_Xenium > 10 AND nCount_Xenium < 1000"
))
stopifnot(identical(section_decision$strict_mask_rule, "primary_include AND NOT qc_review_flag"))

cat("Fixed cell-QC boundary and mask tests passed.\n")
