options(stringsAsFactors = FALSE)

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (!length(script_arg)) stop("Run this test with Rscript.", call. = FALSE)
test_path <- normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = TRUE)
repo_root <- dirname(dirname(test_path))

panel_path <- file.path(repo_root, "config", "eos_gene_sets.tsv")
stopifnot(file.exists(panel_path))
panel <- utils::read.delim(panel_path, check.names = FALSE)
stopifnot(identical(names(panel), c("gene", "gene_set")))
stopifnot(nrow(panel) == 100L, length(unique(panel$gene)) == 100L)
stopifnot(!anyNA(panel), !any(trimws(as.matrix(panel)) == ""))
group_sizes <- table(panel$gene_set)
stopifnot(
  identical(as.integer(group_sizes[c("common", "short_lived", "long_lived")]), c(7L, 47L, 46L))
)

source_path <- file.path(repo_root, "R", "source.R")
stopifnot(file.exists(source_path))
source(source_path)

expect_error <- function(expr, pattern = NULL) {
  error <- tryCatch({ force(expr); NULL }, error = identity)
  stopifnot(inherits(error, "error"))
  if (!is.null(pattern)) stopifnot(grepl(pattern, conditionMessage(error), ignore.case = TRUE))
  invisible(error)
}

test_root <- file.path(tempdir(), "scwat_source_tests")
unlink(test_root, recursive = TRUE, force = TRUE)
dir.create(test_root, recursive = TRUE)
on.exit(unlink(test_root, recursive = TRUE, force = TRUE), add = TRUE)

# Path safety.
stopifnot(assert_path_within(test_root, file.path(test_root, "outputs")))
expect_error(assert_path_within(test_root, "C:/unsafe_output"), "outside")

# Four-section discovery and one-section selection.
input_root <- file.path(test_root, "adipose_data")
dir.create(input_root)
section_names <- sprintf("output-XETG__Region_%d__20260814", 4:1)
invisible(vapply(file.path(input_root, section_names), dir.create, logical(1)))
sections <- discover_xenium_sections(input_root, expected_section_count = 4L)
stopifnot(identical(sections$region_id, paste0("Region_", 1:4)))
stopifnot(identical(discover_one_section(input_root, "Region_3")$region_id, "Region_3"))
expect_error(discover_xenium_sections(input_root, expected_section_count = 3L), "Expected 3")
expect_error(discover_one_section(input_root, "Region_5"), "exactly one")

# Deterministic placeholder metadata.
m1 <- create_synthetic_manifest(paste0("Region_", 1:4), seed = 20260814L)
m2 <- create_synthetic_manifest(paste0("Region_", 1:4), seed = 20260814L)
stopifnot(identical(m1, m2), all(m1$do_not_interpret))
stopifnot(identical(as.integer(sort(table(m1$mouse_id))), c(2L, 2L)))
stopifnot(identical(as.integer(sort(table(m1$side))), c(2L, 2L)))
stopifnot(validate_sample_manifest(m1, paste0("Region_", 1:4))$valid)

# Required-file inventory reports absent files without mutating inputs.
region_dir <- sections$region_dir[sections$region_id == "Region_1"]
inventory <- inventory_section_files(region_dir, "Region_1", calculate_md5 = FALSE)
stopifnot(nrow(inventory) == 8L, !any(inventory$exists))

# Minimal internally aligned Xenium bundle.
matrix_dir <- file.path(region_dir, "cell_feature_matrix")
dir.create(matrix_dir)
write_gz_lines <- function(lines, path) {
  con <- gzfile(path, "wt"); on.exit(close(con), add = TRUE); writeLines(lines, con)
}
write_gz_lines(c("gene1\tGene1\tGene Expression", "ctrl1\tCtrl1\tNegative Control Probe"), file.path(matrix_dir, "features.tsv.gz"))
write_gz_lines(c("cell1", "cell2"), file.path(matrix_dir, "barcodes.tsv.gz"))
write_gz_lines(c("%%MatrixMarket matrix coordinate integer general", "%", "2 2 2", "1 1 3", "1 2 4"), file.path(matrix_dir, "matrix.mtx.gz"))
cells <- data.frame(cell_id = c("cell1", "cell2"), transcript_counts = c(3,4))
con <- gzfile(file.path(region_dir, "cells.csv.gz"), "wt"); utils::write.csv(cells, con, row.names = FALSE); close(con)
integrity <- validate_section_integrity(region_dir, "Region_1")
stopifnot(integrity$dimension_match, integrity$matrix_cells == 2L, integrity$matrix_features == 2L)

write_gz_lines("cell1", file.path(matrix_dir, "barcodes.tsv.gz"))
expect_error(validate_section_integrity(region_dir, "Region_1"), "integrity")

# Xenium alarms, panel reconciliation, sparse import, and cell QC.
write_gz_lines(c("cell1", "cell2"), file.path(matrix_dir, "barcodes.tsv.gz"))
cells <- data.frame(
  cell_id = c("cell1", "cell2"), x_centroid = c(1, 2), y_centroid = c(3, 4),
  transcript_counts = c(3, 4), control_probe_counts = c(0, 1),
  genomic_control_counts = c(0, 0), control_codeword_counts = c(0, 0),
  total_counts = c(3, 5), cell_area = c(20, 80), nucleus_count = c(1, 2),
  segmentation_method = c("nucleus_expansion", "nucleus_expansion")
)
con <- gzfile(file.path(region_dir, "cells.csv.gz"), "wt"); utils::write.csv(cells, con, row.names = FALSE); close(con)

alarm_html <- paste0(
  '<html>"alarms":{"alarms":[{"raw_value":true,"formatted_value":"true",',
  '"raised":true,"title":"Poor cycles","message":"Review cycles",',
  '"level":"ERROR","id":"poor_quality_cycles_detected"}]},"sample":{}</html>'
)
writeLines(alarm_html, file.path(region_dir, "analysis_summary.html"))
alarms <- extract_analysis_alarms(file.path(region_dir, "analysis_summary.html"))
stopifnot(nrow(alarms) == 1L, alarms$level[[1]] == "ERROR", alarms$id[[1]] == "poor_quality_cycles_detected")

panel_json <- paste0(
  '{"payload":{"targets":[',
  '{"source":{"category":"current"},"type":{"descriptor":"gene","data":{"name":"Gene1"}}},',
  '{"source":{"category":"current"},"type":{"descriptor":"gene","data":{"name":"ExtraGene"}}}',
  ']}}'
)
writeLines(panel_json, file.path(region_dir, "gene_panel.json"))
installed <- read_custom_panel_genes(file.path(region_dir, "gene_panel.json"))
panel_check <- reconcile_panel(c("Gene1", "Gene2"), installed)
stopifnot(sum(panel_check$status == "MATCH") == 1L, sum(panel_check$status == "MISSING") == 1L, sum(panel_check$status == "EXTRA") == 1L)

imported <- import_xenium_mex(region_dir)
stopifnot(inherits(imported$counts, "sparseMatrix"), identical(dim(imported$counts), c(1L, 2L)))
stopifnot(identical(colnames(imported$counts), imported$cells$cell_id))

qc <- calculate_xenium_cell_qc(imported$counts, imported$cells, "Region_1")
stopifnot(nrow(qc$cell_metadata) == 2L, nrow(qc$thresholds) == 4L, nrow(qc$summary) == 1L)
stopifnot(qc$cell_metadata$multiple_nuclei_flag[[2]], qc$cell_metadata$segmentation_multiplet_flag[[2]])
stopifnot(all(c("qc_core_pass", "qc_review_flag", "high_control_flag") %in% names(qc$cell_metadata)))

cat("All reusable scWAT Xenium function tests passed.\n")
