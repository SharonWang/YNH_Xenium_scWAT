#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
test_path <- normalizePath(sub("^--file=", "", script_arg[[1L]]), winslash = "/", mustWork = TRUE)
repo_root <- dirname(dirname(test_path))
source(file.path(repo_root, "R", "source.R"))

expect_error <- function(expr, pattern) {
  error <- tryCatch({ force(expr); NULL }, error = identity)
  stopifnot(inherits(error, "error"), grepl(pattern, conditionMessage(error), ignore.case = TRUE))
}

test_root <- file.path(tempdir(), "scwat_source_contract")
unlink(test_root, recursive = TRUE, force = TRUE)
dir.create(test_root, recursive = TRUE)
on.exit(unlink(test_root, recursive = TRUE, force = TRUE), add = TRUE)

stopifnot(assert_path_within(test_root, file.path(test_root, "adipose_analysis")))
expect_error(assert_path_within(test_root, "C:/unsafe"), "outside")

input_root <- file.path(test_root, "adipose_data")
dir.create(input_root)
for (index in 1:4) dir.create(file.path(input_root, sprintf("output-XETG__Region_%d__unit", index)))
sections <- discover_xenium_sections(input_root, 4L)
stopifnot(identical(sections$region_id, expected_scwat_regions()))
stopifnot(identical(discover_one_section(input_root, "Region_3")$region_id, "Region_3"))

manifest <- utils::read.delim(file.path(repo_root, "config", "scwat_sample_manifest.tsv"), check.names = FALSE)
stopifnot(validate_sample_manifest(manifest, expected_scwat_regions())$valid)
stopifnot(identical(as.integer(table(manifest$mouse_id)), c(2L, 2L)))

# A real sparse count fixture exercises fixed QC, non-destructive masks, and
# mouse-aware summary behavior without mocking any pipeline boundary.
counts <- Matrix::Matrix(matrix(c(
  2, 0, 1,
  1, 1, 0,
  0, 1, 1,
  0, 0, 1,
  0, 0, 1,
  0, 0, 1
), nrow = 6L), sparse = TRUE)
rownames(counts) <- paste0("Gene", 1:6)
colnames(counts) <- paste0("cell", 1:3)
cells <- data.frame(
  cell_id = colnames(counts), x_centroid = 1:3, y_centroid = 3:1,
  transcript_counts = Matrix::colSums(counts), control_probe_counts = 0,
  genomic_control_counts = 0, control_codeword_counts = 0,
  total_counts = Matrix::colSums(counts), cell_area = c(50, 60, 70),
  nucleus_count = 1L, segmentation_method = "nucleus_expansion",
  stringsAsFactors = FALSE
)
fixed <- read_fixed_cell_qc_thresholds(file.path(repo_root, "config", "fixed_cell_qc_thresholds.tsv"))
qc <- calculate_xenium_cell_qc(counts, cells, "Region_1", fixed)
masks <- build_cell_downstream_masks(qc$cell_metadata, provenance = "unit::LOCAL_SUBSET")
stopifnot(
  nrow(masks) == ncol(counts),
  identical(masks$primary_include, qc$cell_metadata$qc_core_pass),
  all(!masks$strict_include | masks$primary_include),
  all(!masks$hotspot_sensitivity_include | masks$primary_include)
)

section_summary <- data.frame(
  region_id = expected_scwat_regions(), input_cells = c(100L, 110L, 120L, 130L),
  core_qc_pass = c(90L, 99L, 108L, 117L), stringsAsFactors = FALSE
)
mouse <- summarise_scwat_mouse_sections(section_summary, manifest)
stopifnot(nrow(mouse$within_mouse) == 4L, nrow(mouse$mouse_summary) == 2L)

plots <- plot_section_qc(masks, "Region_1")
stopifnot(identical(names(plots), c("counts", "features", "area", "spatial")))
stopifnot(inherits(plots$spatial$coordinates, "CoordCartesian"), identical(plots$spatial$coordinates$ratio, 1))

cat("All active scWAT Xenium source tests passed.\n")
