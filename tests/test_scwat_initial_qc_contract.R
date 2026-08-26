#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
test_path <- normalizePath(sub("^--file=", "", script_arg[[1L]]), winslash = "/", mustWork = TRUE)
repo_root <- dirname(dirname(test_path))
source(file.path(repo_root, "R", "source.R"))

expect_error <- function(expr, pattern) {
  message <- tryCatch({ force(expr); NA_character_ }, error = function(e) conditionMessage(e))
  stopifnot(!is.na(message), grepl(pattern, message, ignore.case = TRUE))
}

# Break caught: a six-region or inherited anchor/sensitivity contract must not
# enter the four-region initial-QC pipeline.
regions <- paste0("Region_", 1:4)
stopifnot(identical(expected_scwat_regions(), regions))
stopifnot(identical(names(section_palette()), regions))

index <- data.frame(
  region_id = regions,
  run_label = rep("unit_run", 4L),
  execution_mode = rep("LOCAL_SUBSET", 4L),
  section_output_dir = file.path("sections", regions),
  stringsAsFactors = FALSE
)
validated <- validate_scwat_section_bundle_index(index, regions)
stopifnot(identical(validated$region_id, regions))
expect_error(validate_scwat_section_bundle_index(index[-4L, ], regions), "exactly four")
duplicate <- index; duplicate$region_id[[4L]] <- "Region_3"
expect_error(validate_scwat_section_bundle_index(duplicate, regions), "duplicate|exactly four")
mixed_run <- index; mixed_run$run_label[[4L]] <- "other"
expect_error(validate_scwat_section_bundle_index(mixed_run, regions), "run label")
mixed_mode <- index; mixed_mode$execution_mode[[4L]] <- "FULL_HPC"
expect_error(validate_scwat_section_bundle_index(mixed_mode, regions), "execution mode")

# Break caught: the slide summary must respect mouse as the replicate and must
# not invent a left/right or anatomical-position variable for these sections.
manifest <- utils::read.delim(file.path(repo_root, "config", "scwat_sample_manifest.tsv"), check.names = FALSE)
section_summary <- data.frame(
  region_id = regions,
  input_cells = c(100L, 110L, 120L, 130L),
  core_qc_pass = c(90L, 99L, 108L, 117L),
  stringsAsFactors = FALSE
)
mouse_summary <- summarise_scwat_mouse_sections(section_summary, manifest)
stopifnot(
  identical(mouse_summary$region_summary$region_id, regions),
  nrow(mouse_summary$within_mouse) == 4L,
  nrow(mouse_summary$mouse_summary) == 2L,
  identical(mouse_summary$mouse_summary$mouse_id, c("Mouse_1", "Mouse_2")),
  all(mouse_summary$mouse_summary$n_sections == 2L),
  !any(c("side", "position") %in% names(mouse_summary$region_summary))
)

gates <- expand.grid(region_id = regions, gate = "overall", stringsAsFactors = FALSE)
gates$status <- c("PASS", "WARN", "HOLD", "PASS")
readiness <- derive_scwat_region_readiness(gates, regions)
stopifnot(
  identical(readiness$region_id, regions),
  identical(readiness$status, gates$status),
  !any(readiness$status %in% c("PRIMARY", "PRIMARY_CONDITIONAL", "SENSITIVITY_ONLY"))
)

cat("Four-region scWAT initial-QC contracts passed.\n")
