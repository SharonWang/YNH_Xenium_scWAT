options(stringsAsFactors = FALSE)

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (!length(script_arg)) stop("Run this test with Rscript.", call. = FALSE)
test_path <- normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = TRUE)
repo_root <- dirname(dirname(test_path))
source(file.path(repo_root, "R", "source.R"))

expect_error <- function(expr, pattern = NULL) {
  error <- tryCatch({ force(expr); NULL }, error = identity)
  stopifnot(inherits(error, "error"))
  if (!is.null(pattern)) stopifnot(grepl(pattern, conditionMessage(error), ignore.case = TRUE))
  invisible(error)
}

# A broken config parser would silently change scientific thresholds.
config_path <- file.path(repo_root, "config", "extended_qc_defaults.tsv")
config <- read_extended_qc_config(config_path)
stopifnot(config$qv_threshold == 20, config$spatial_k == 15L, config$permutations == 999L)
stopifnot(config$seed == 20260814L, is.integer(config$seed), is.numeric(config$dense_quantile))

bad_config <- utils::read.delim(config_path, check.names = FALSE)
bad_config <- rbind(bad_config, bad_config[1, , drop = FALSE])
bad_path <- file.path(tempdir(), "duplicate_extended_qc_config.tsv")
utils::write.table(bad_config, bad_path, sep = "\t", quote = FALSE, row.names = FALSE)
expect_error(read_extended_qc_config(bad_path), "duplicate")

# AUTO must select full mode only when the large transcript input exists.
mode_root <- file.path(tempdir(), "extended_qc_modes")
unlink(mode_root, recursive = TRUE, force = TRUE)
subset_region_dir <- file.path(mode_root, "subset")
full_region_dir <- file.path(mode_root, "full")
dir.create(subset_region_dir, recursive = TRUE)
dir.create(full_region_dir, recursive = TRUE)
invisible(file.create(file.path(full_region_dir, "transcripts.parquet")))
stopifnot(resolve_extended_qc_mode("AUTO", subset_region_dir) == "LOCAL_SUBSET")
stopifnot(resolve_extended_qc_mode("AUTO", full_region_dir) == "FULL_HPC")
stopifnot(resolve_extended_qc_mode("local_subset", full_region_dir) == "LOCAL_SUBSET")
expect_error(resolve_extended_qc_mode("BAD", subset_region_dir), "AUTO")

# Local mode may skip large inputs but must report the skip explicitly.
preflight <- extended_qc_preflight("LOCAL_SUBSET", subset_region_dir, config)
stopifnot(preflight$status[preflight$check == "transcripts_parquet"] == "SKIP_ALLOWED")
stopifnot(preflight$status[preflight$check == "arrow"] == "SKIP_ALLOWED")
full_preflight <- extended_qc_preflight("FULL_HPC", full_region_dir, config)
if (!requireNamespace("arrow", quietly = TRUE)) stopifnot(full_preflight$status[full_preflight$check == "arrow"] == "FAIL")
if (!requireNamespace("RANN", quietly = TRUE)) stopifnot(full_preflight$status[full_preflight$check == "RANN"] == "FAIL")
stopifnot("dplyr" %in% full_preflight$check)
if (!requireNamespace("dplyr", quietly = TRUE)) stopifnot(full_preflight$status[full_preflight$check == "dplyr"] == "FAIL")

# Alarm evidence must never imply that the exact cycle or gene mapping is known.
alarm_fixture <- data.frame(
  raw_value = TRUE, formatted_value = "True", raised = TRUE,
  title = "Poor quality imaging cycles detected", message = "Diagnostic message",
  level = "ERROR", id = "poor_quality_cycles_detected", stringsAsFactors = FALSE
)
evidence <- build_cycle_alarm_evidence(alarm_fixture, "Region_1")
stopifnot(evidence$evidence_status == "DIRECT_EVIDENCE")
stopifnot(evidence$cycle_identity_status == "CYCLE_IDENTITY_UNRESOLVED_REQUIRES_10X")
stopifnot(evidence$gene_effect_status == "GENE_EFFECT_UNCONFIRMED")
no_evidence <- build_cycle_alarm_evidence(alarm_fixture[0, ], "Region_3")
stopifnot(no_evidence$evidence_status == "NO_ALARM_REPORTED", nrow(no_evidence) == 1L)

# Transcript summaries use literal hand-calculated QV and codeword results.
transcript_fixture <- data.frame(
  feature_name = c("A", "A", "B"), qv = c(30, 10, 25), codeword_index = c(1L, 1L, 2L),
  stringsAsFactors = FALSE
)
schema <- resolve_transcript_schema(names(transcript_fixture))
stopifnot(schema$gene == "feature_name", schema$qv == "qv", schema$codeword == "codeword_index")
expect_error(resolve_transcript_schema(c("gene", "feature_name", "qv")), "ambiguous")
expect_error(resolve_transcript_schema(c("feature_name", "x")), "QV")
transcript_qc <- summarise_transcript_quality_table(transcript_fixture, "Region_1", 20)
stopifnot(transcript_qc$transcript_rows[transcript_qc$gene == "A"] == 2L)
stopifnot(transcript_qc$mean_qv[transcript_qc$gene == "A"] == 20)
stopifnot(transcript_qc$fraction_q20[transcript_qc$gene == "A"] == 0.5)
stopifnot(transcript_qc$represented_codewords[transcript_qc$gene == "A"] == 1L)

# Matrix summaries preserve all panel genes, including zero-count genes.
require_package("Matrix")
count_fixture <- Matrix::Matrix(matrix(c(2, 0, 3, 0, 0, 0), nrow = 3L, byrow = TRUE), sparse = TRUE)
rownames(count_fixture) <- c("A", "B", "C")
colnames(count_fixture) <- c("cell1", "cell2")
gene_sets <- c(A = "common", B = "short_lived", C = "long_lived")
matrix_qc <- summarise_gene_matrix_qc(count_fixture, "Region_1", gene_sets)
stopifnot(matrix_qc$raw_counts[matrix_qc$gene == "A"] == 2)
stopifnot(matrix_qc$detection_fraction[matrix_qc$gene == "B"] == 0.5)
stopifnot(matrix_qc$raw_counts[matrix_qc$gene == "C"] == 0)
combined_qc <- combine_gene_quality(matrix_qc, transcript_qc)
stopifnot(nrow(combined_qc) == 3L, combined_qc$transcript_status[combined_qc$gene == "C"] == "NO_TRANSCRIPTS")

cat("All extended scWAT Xenium QC tests passed.\n")
