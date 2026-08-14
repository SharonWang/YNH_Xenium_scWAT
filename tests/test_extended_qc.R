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
file.create(file.path(full_region_dir, "transcripts.parquet"))
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

cat("Extended QC configuration and mode tests passed.\n")
