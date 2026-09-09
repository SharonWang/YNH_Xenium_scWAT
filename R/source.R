# scWAT Xenium reusable QC functions
#
# Active contract:
# - This file is sourced by the four region QC notebooks, the slide summary,
#   and active exploratory Region 3 notebooks.
# - Raw matrices/objects are never overwritten; QC decisions are added as
#   metadata, summaries, masks, or separately written artifacts.
# - Functions validate required columns, alignment, path containment, and
#   output provenance before returning or writing results.
# - Functions with no active notebook/test/support consumer are preserved in
#   R/source_bk.R and are intentionally not loaded by the QC notebooks.
#
# Documentation convention: every function states its purpose, required and
# optional inputs, and output/write behavior immediately above its definition.

options(stringsAsFactors = FALSE)

#' Use a fallback for null, empty, missing, or blank values.
#'
#' @param x Candidate value.
#' @param y Fallback returned when `x` has no usable first value.
#' @return `x` when usable; otherwise `y`.
`%||%` <- function(x, y) {
  if (is.null(x) || !length(x)) return(y)
  first <- x[[1]]
  missing_scalar <- is.atomic(first) && length(first) == 1L && (is.na(first) || !nzchar(as.character(first)))
  if (missing_scalar) y else x
}

#' Canonical path.
#'
#' @param path Required `path` input; validated before computation.
#' @return A deterministic derived value or annotated copy; the supplied raw data are not overwritten.
canonical_path <- function(path) {
  value <- normalizePath(path, winslash = "/", mustWork = FALSE)
  value <- sub("/+$", "", value)
  if (.Platform$OS.type == "windows") value <- tolower(value)
  value
}

#' Assert path within.
#'
#' @param project_root Required `project_root` input; validated before computation.
#' @param candidate Required `candidate` input; validated before computation.
#' @return Validated input or `TRUE`; invalid input stops with an informative error.
assert_path_within <- function(project_root, candidate) {
  root <- canonical_path(project_root)
  path <- canonical_path(candidate)
  valid <- identical(path, root) || startsWith(path, paste0(root, "/"))
  if (!valid) stop(sprintf("Unsafe path outside project root: %s", path), call. = FALSE)
  invisible(TRUE)
}

#' Validate runtime paths.
#'
#' @param project_root Required `project_root` input; validated before computation.
#' @param input_root Required `input_root` input; validated before computation.
#' @param output_root Required `output_root` input; validated before computation.
#' @param temp_root Optional `temp_root` input with the default shown in the function signature.
#' @return Validated input or `TRUE`; invalid input stops with an informative error.
validate_runtime_paths <- function(project_root, input_root, output_root, temp_root = tempdir()) {
  if (!dir.exists(project_root)) stop(sprintf("Project root does not exist: %s", project_root), call. = FALSE)
  if (!dir.exists(input_root)) stop(sprintf("Input root does not exist: %s", input_root), call. = FALSE)
  for (path in c(input_root, output_root, temp_root)) assert_path_within(project_root, path)
  invisible(TRUE)
}

#' Discover xenium sections.
#'
#' @param input_root Required `input_root` input; validated before computation.
#' @param expected_section_count Optional `expected_section_count` input with the default shown in the function signature.
#' @return A deterministic derived value or annotated copy; the supplied raw data are not overwritten.
discover_xenium_sections <- function(input_root, expected_section_count = 4L) {
  if (!dir.exists(input_root)) stop(sprintf("Input root does not exist: %s", input_root), call. = FALSE)
  dirs <- list.dirs(input_root, recursive = FALSE, full.names = TRUE)
  matched <- grepl("__Region_[0-9]+__", basename(dirs))
  dirs <- dirs[matched]
  region_number <- as.integer(sub(".*__Region_([0-9]+)__.*", "\\1", basename(dirs)))
  if (length(dirs) != as.integer(expected_section_count)) {
    stop(sprintf("Expected %d Xenium section directories under %s but found %d.", as.integer(expected_section_count), input_root, length(dirs)), call. = FALSE)
  }
  if (anyDuplicated(region_number)) stop("Duplicate Xenium Region identifiers detected.", call. = FALSE)
  order_index <- order(region_number)
  data.frame(
    region_id = paste0("Region_", region_number[order_index]),
    region_number = region_number[order_index],
    region_dir = normalizePath(dirs[order_index], winslash = "/", mustWork = TRUE),
    stringsAsFactors = FALSE
  )
}

#' Discover one section.
#'
#' @param input_root Required `input_root` input; validated before computation.
#' @param region_id Required `region_id` input; validated before computation.
#' @return A deterministic derived value or annotated copy; the supplied raw data are not overwritten.
discover_one_section <- function(input_root, region_id) {
  if (!grepl("^Region_[1-9][0-9]*$", region_id)) stop("region_id must use Region_<integer> format.", call. = FALSE)
  dirs <- list.dirs(input_root, recursive = FALSE, full.names = TRUE)
  matched <- dirs[grepl(sprintf("__%s__", region_id), basename(dirs), fixed = TRUE)]
  if (length(matched) != 1L) stop(sprintf("Expected exactly one directory for %s but found %d.", region_id, length(matched)), call. = FALSE)
  data.frame(
    region_id = region_id,
    region_number = as.integer(sub("^Region_", "", region_id)),
    region_dir = normalizePath(matched, winslash = "/", mustWork = TRUE),
    stringsAsFactors = FALSE
  )
}
#' Validate sample manifest.
#'
#' @param manifest Required `manifest` input; validated before computation.
#' @param expected_regions Required `expected_regions` input; validated before computation.
#' @return Validated input or `TRUE`; invalid input stops with an informative error.
validate_sample_manifest <- function(manifest, expected_regions) {
  required <- c("tissue", "region_id", "mouse_id", "section_id", "biological_replicate_id", "technical_replicate_id", "metadata_status")
  issues <- character()
  missing <- setdiff(required, names(manifest))
  if (length(missing)) issues <- c(issues, paste0("missing_columns:", paste(missing, collapse = ",")))
  if (!length(missing)) {
    if (anyDuplicated(manifest$region_id)) issues <- c(issues, "duplicate_region_id")
    if (anyDuplicated(manifest$section_id)) issues <- c(issues, "duplicate_section_id")
    if (!setequal(manifest$region_id, expected_regions)) issues <- c(issues, "region_set_mismatch")
    if ("side" %in% names(manifest) && any(!manifest$side %in% c("Left", "Right", "Not_applicable", "Unknown"))) issues <- c(issues, "invalid_side")
    values <- as.matrix(manifest[, required, drop = FALSE])
    if (anyNA(values) || any(trimws(values) == "")) issues <- c(issues, "missing_required_values")
  }
  list(valid = !length(issues), issues = issues)
}

#' Xenium required files.
#'
#' @return A deterministic derived value or annotated copy; the supplied raw data are not overwritten.
xenium_required_files <- function() {
  c(
    "experiment.xenium", "metrics_summary.csv", "analysis_summary.html", "gene_panel.json",
    "cells.csv.gz", "cell_feature_matrix/features.tsv.gz",
    "cell_feature_matrix/barcodes.tsv.gz", "cell_feature_matrix/matrix.mtx.gz"
  )
}

#' Inventory section files.
#'
#' @param region_dir Required `region_dir` input; validated before computation.
#' @param region_id Required `region_id` input; validated before computation.
#' @param calculate_md5 Optional `calculate_md5` input with the default shown in the function signature.
#' @return A deterministic derived value or annotated copy; the supplied raw data are not overwritten.
inventory_section_files <- function(region_dir, region_id, calculate_md5 = TRUE) {
  relative <- xenium_required_files()
  paths <- file.path(region_dir, relative)
  exists <- file.exists(paths)
  info <- file.info(paths)
  hashes <- rep(NA_character_, length(paths))
  if (isTRUE(calculate_md5) && any(exists)) hashes[exists] <- unname(tools::md5sum(paths[exists]))
  data.frame(
    region_id = region_id, relative_path = relative, exists = exists,
    size_bytes = ifelse(exists, info$size, NA_real_),
    modified_utc = ifelse(exists, format(info$mtime, tz = "UTC", usetz = TRUE), NA_character_),
    md5 = hashes, stringsAsFactors = FALSE
  )
}

#' Read mtx dimensions.
#'
#' @param path Required `path` input; validated before computation.
#' @return Parsed and validated R data with input row or matrix alignment preserved.
read_mtx_dimensions <- function(path) {
  if (!file.exists(path)) stop(sprintf("Matrix file not found: %s", path), call. = FALSE)
  con <- gzfile(path, "rt"); on.exit(close(con), add = TRUE)
  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (!length(line)) stop("Matrix Market dimension line is absent.", call. = FALSE)
    if (!startsWith(line, "%")) break
  }
  values <- suppressWarnings(as.numeric(strsplit(trimws(line), "[[:space:]]+")[[1]]))
  if (length(values) != 3L || anyNA(values)) stop("Invalid Matrix Market dimensions.", call. = FALSE)
  setNames(as.integer(values), c("features", "cells", "nonzero"))
}

#' Read gz rows.
#'
#' @param path Required `path` input; validated before computation.
#' @param header Optional `header` input with the default shown in the function signature.
#' @return Parsed and validated R data with input row or matrix alignment preserved.
read_gz_rows <- function(path, header = FALSE) {
  utils::read.delim(gzfile(path), header = header, quote = "", check.names = FALSE, stringsAsFactors = FALSE)
}

#' Validate section integrity.
#'
#' @param region_dir Required `region_dir` input; validated before computation.
#' @param region_id Required `region_id` input; validated before computation.
#' @return Validated input or `TRUE`; invalid input stops with an informative error.
validate_section_integrity <- function(region_dir, region_id) {
  matrix_dir <- file.path(region_dir, "cell_feature_matrix")
  dimensions <- read_mtx_dimensions(file.path(matrix_dir, "matrix.mtx.gz"))
  feature_rows <- nrow(read_gz_rows(file.path(matrix_dir, "features.tsv.gz")))
  barcode_rows <- length(scan(gzfile(file.path(matrix_dir, "barcodes.tsv.gz")), what = character(), quiet = TRUE))
  cells <- utils::read.csv(gzfile(file.path(region_dir, "cells.csv.gz")), stringsAsFactors = FALSE, check.names = FALSE)
  duplicate_cell_ids <- if ("cell_id" %in% names(cells)) sum(duplicated(cells$cell_id)) else NA_integer_
  dimension_match <- dimensions[["features"]] == feature_rows && dimensions[["cells"]] == barcode_rows && barcode_rows == nrow(cells) && identical(duplicate_cell_ids, 0L)
  out <- data.frame(
    region_id = region_id, matrix_features = dimensions[["features"]], matrix_cells = dimensions[["cells"]],
    matrix_nonzero = dimensions[["nonzero"]], feature_rows = feature_rows, barcode_rows = barcode_rows,
    cell_rows = nrow(cells), duplicate_cell_ids = duplicate_cell_ids, dimension_match = dimension_match,
    stringsAsFactors = FALSE
  )
  if (!dimension_match) stop(sprintf("Xenium integrity mismatch for %s.", region_id), call. = FALSE)
  out
}

#' Flatten optional-package tables to a rectangular TSV-safe data frame
#'
#' Some R packages return nominal data frames containing matrix or list-valued
#' columns. Base `write.table()` cannot reliably construct column names for
#' these objects and can fail with a dimnames/array-extent error. Matrix/data
#' frame columns are expanded with their parent column name; list cells are
#' serialized to a deterministic ` | `-separated character value.
#'
#' @param x Object coercible to a data frame.
#'
#' @return A rectangular data frame containing only atomic columns and the same
#'   number and order of rows as `x`.
flatten_tsv_table <- function(x) {
  table <- if (is.data.frame(x)) x else as.data.frame(x, stringsAsFactors = FALSE)
  if (!ncol(table)) return(table)
  n_rows <- nrow(table)
  pieces <- lapply(seq_along(table), function(index) {
    column <- table[[index]]
    parent <- names(table)[[index]]
    if (is.matrix(column) || is.data.frame(column)) {
      nested_value <- if (is.matrix(column)) unclass(column) else column
      expanded <- as.data.frame(nested_value, stringsAsFactors = FALSE, check.names = FALSE)
      if (nrow(expanded) != n_rows) {
        stop("Nested TSV column has incompatible row count: ", parent, call. = FALSE)
      }
      child_names <- names(expanded)
      if (is.null(child_names) || any(!nzchar(child_names))) {
        child_names <- paste0("V", seq_len(ncol(expanded)))
      }
      names(expanded) <- paste(parent, child_names, sep = ".")
      return(expanded)
    }
    if (is.list(column)) {
      serialize_cell <- function(value) {
        if (is.null(value) || !length(value)) return(NA_character_)
        if (is.atomic(value)) return(paste(as.character(value), collapse = " | "))
        paste(capture.output(dput(value)), collapse = " ")
      }
      column <- vapply(column, serialize_cell, character(1))
    }
    output <- data.frame(column, stringsAsFactors = FALSE, check.names = FALSE)
    names(output) <- parent
    output
  })
  output <- do.call(cbind, pieces)
  names(output) <- make.unique(names(output), sep = ".")
  rownames(output) <- rownames(table)
  output
}

#' Write tsv.
#'
#' @param x Required `x` input; validated before computation.
#' @param path Required `path` input; validated before computation.
#' @param project_root Required `project_root` input; validated before computation.
#' @return Written artifact path(s) or an auditable one-row write summary.
write_tsv <- function(x, path, project_root) {
  assert_path_within(project_root, path)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  table <- flatten_tsv_table(x)
  tryCatch(
    utils::write.table(table, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"),
    error = function(error) {
      stop("Failed to write TSV '", path, "': ", conditionMessage(error), call. = FALSE)
    }
  )
  invisible(path)
}

#' Require package.
#'
#' @param package Required `package` input; validated before computation.
#' @return A deterministic derived value or annotated copy; the supplied raw data are not overwritten.
require_package <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) stop(sprintf("Required R package '%s' is unavailable.", package), call. = FALSE)
}
#' Empty alarm table.
#'
#' @return A deterministic derived value or annotated copy; the supplied raw data are not overwritten.
empty_alarm_table <- function() {
  data.frame(
    raw_value = logical(), formatted_value = character(), raised = logical(),
    title = character(), message = character(), level = character(), id = character(),
    stringsAsFactors = FALSE
  )
}
#' Read xenium features.
#'
#' @param path Required `path` input; validated before computation.
#' @return Parsed and validated R data with input row or matrix alignment preserved.
read_xenium_features <- function(path) {
  features <- read_gz_rows(path, header = FALSE)
  if (ncol(features) != 3L) stop(sprintf("Expected three columns in %s", path), call. = FALSE)
  names(features) <- c("feature_id", "feature_name", "feature_type")
  features
}

#' Read xenium barcodes.
#'
#' @param path Required `path` input; validated before computation.
#' @return Parsed and validated R data with input row or matrix alignment preserved.
read_xenium_barcodes <- function(path) {
  scan(gzfile(path), what = character(), quiet = TRUE)
}

#' Import xenium mex.
#'
#' @param region_dir Required `region_dir` input; validated before computation.
#' @return Parsed and validated R data with input row or matrix alignment preserved.
import_xenium_mex <- function(region_dir) {
  require_package("Matrix")
  matrix_dir <- file.path(region_dir, "cell_feature_matrix")
  features_path <- file.path(matrix_dir, "features.tsv.gz")
  barcodes_path <- file.path(matrix_dir, "barcodes.tsv.gz")
  matrix_path <- file.path(matrix_dir, "matrix.mtx.gz")
  cells_path <- file.path(region_dir, "cells.csv.gz")
  required <- c(features_path, barcodes_path, matrix_path, cells_path)
  missing <- required[!file.exists(required)]
  if (length(missing)) stop(sprintf("Missing Xenium import inputs: %s", paste(missing, collapse = ", ")), call. = FALSE)
  features <- read_xenium_features(features_path)
  barcodes <- read_xenium_barcodes(barcodes_path)
  raw_matrix <- methods::as(Matrix::readMM(gzfile(matrix_path)), "CsparseMatrix")
  if (!identical(dim(raw_matrix), c(nrow(features), length(barcodes)))) stop("Matrix dimensions do not match features/barcodes.", call. = FALSE)
  rownames(raw_matrix) <- make.unique(features$feature_name)
  colnames(raw_matrix) <- barcodes
  cells <- utils::read.csv(gzfile(cells_path), stringsAsFactors = FALSE, check.names = FALSE)
  if (!"cell_id" %in% names(cells) || anyDuplicated(cells$cell_id)) stop("Cell metadata requires unique cell_id values.", call. = FALSE)
  order_index <- match(barcodes, cells$cell_id)
  if (anyNA(order_index)) stop("Matrix barcodes and cell metadata are not aligned.", call. = FALSE)
  cells <- cells[order_index, , drop = FALSE]
  rownames(cells) <- cells$cell_id
  gene_rows <- features$feature_type == "Gene Expression"
  counts <- methods::as(raw_matrix[gene_rows, , drop = FALSE], "CsparseMatrix")
  list(
    counts = counts, cells = cells, features = features,
    feature_type_summary = as.data.frame(table(features$feature_type), stringsAsFactors = FALSE),
    raw_matrix_dimensions = c(features = nrow(features), cells = length(barcodes), nonzero = length(raw_matrix@x))
  )
}

#' Safe quantile.
#'
#' @param x Required `x` input; validated before computation.
#' @param probability Required `probability` input; validated before computation.
#' @return A deterministic derived value or annotated copy; the supplied raw data are not overwritten.
safe_quantile <- function(x, probability) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  unname(stats::quantile(x, probability, names = FALSE, na.rm = TRUE, type = 7))
}

#' Robust interval.
#'
#' @param x Required `x` input; validated before computation.
#' @param lower_mads Optional `lower_mads` input with the default shown in the function signature.
#' @param upper_mads Optional `upper_mads` input with the default shown in the function signature.
#' @param floor_value Optional `floor_value` input with the default shown in the function signature.
#' @return A deterministic derived value or annotated copy; the supplied raw data are not overwritten.
robust_interval <- function(x, lower_mads = 3, upper_mads = 5, floor_value = -Inf) {
  x <- x[is.finite(x)]
  if (!length(x)) stop("Cannot calculate a robust interval from empty data.", call. = FALSE)
  median_value <- stats::median(x)
  mad_value <- stats::mad(x, center = median_value, constant = 1.4826)
  if (!is.finite(mad_value) || mad_value == 0) {
    lower <- safe_quantile(x, 0.01); upper <- safe_quantile(x, 0.99); method <- "q01_q99_fallback"
  } else {
    lower <- median_value - lower_mads * mad_value
    upper <- median_value + upper_mads * mad_value
    method <- sprintf("median_minus_%gMAD_plus_%gMAD", lower_mads, upper_mads)
  }
  c(lower = max(floor_value, lower), upper = upper, median = median_value, mad = mad_value, method = method)
}

#' Read and validate the prespecified Xenium cell-complexity thresholds.
#'
#' @param path Existing TSV with columns metric, bound, value, inclusive, and
#'   purpose. Exactly one lower and upper row is required for nFeature_Xenium
#'   and nCount_Xenium; this pipeline requires all four bounds to be exclusive.
#' @return Named list with numeric feature_lower, feature_upper, count_lower,
#'   and count_upper values used by calculate_xenium_cell_qc().
read_fixed_cell_qc_thresholds <- function(path) {
  if (!file.exists(path)) stop(sprintf("Fixed cell-QC threshold file not found: %s", path), call. = FALSE)
  thresholds <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  required <- c("metric", "bound", "value", "inclusive", "purpose")
  missing <- setdiff(required, names(thresholds))
  if (length(missing)) stop(sprintf("Fixed cell-QC thresholds missing columns: %s", paste(missing, collapse = ", ")), call. = FALSE)
  expected_keys <- c(
    "nFeature_Xenium::lower", "nFeature_Xenium::upper",
    "nCount_Xenium::lower", "nCount_Xenium::upper"
  )
  keys <- paste(thresholds$metric, thresholds$bound, sep = "::")
  if (nrow(thresholds) != 4L || anyDuplicated(keys) || !setequal(keys, expected_keys)) {
    stop("Fixed cell-QC thresholds must contain exactly one lower and upper row for nFeature_Xenium and nCount_Xenium.", call. = FALSE)
  }
  inclusive <- toupper(trimws(as.character(thresholds$inclusive)))
  if (any(inclusive != "FALSE")) stop("The approved Xenium primary bounds must all be exclusive.", call. = FALSE)
  values <- suppressWarnings(as.numeric(thresholds$value))
  if (anyNA(values) || any(!is.finite(values))) stop("Fixed cell-QC threshold values must be finite numbers.", call. = FALSE)
  value_for <- function(metric, bound) values[keys == paste(metric, bound, sep = "::")][[1L]]
  out <- list(
    feature_lower = value_for("nFeature_Xenium", "lower"),
    feature_upper = value_for("nFeature_Xenium", "upper"),
    count_lower = value_for("nCount_Xenium", "lower"),
    count_upper = value_for("nCount_Xenium", "upper")
  )
  if (out$feature_lower >= out$feature_upper || out$count_lower >= out$count_upper) {
    stop("Each fixed cell-QC lower bound must be smaller than its upper bound.", call. = FALSE)
  }
  out
}

#' Apply the approved fixed, open-interval Xenium primary bounds.
#'
#' @param n_feature Numeric detected-gene counts per cell.
#' @param n_count Numeric transcript counts per cell, aligned to n_feature.
#' @param fixed_thresholds Named list returned by
#'   read_fixed_cell_qc_thresholds().
#' @return Logical vector; TRUE only when both metrics lie strictly inside
#'   their configured lower and upper bounds.
apply_fixed_primary_bounds <- function(n_feature, n_count, fixed_thresholds) {
  required <- c("feature_lower", "feature_upper", "count_lower", "count_upper")
  if (!is.list(fixed_thresholds) || length(setdiff(required, names(fixed_thresholds)))) {
    stop("fixed_thresholds must be returned by read_fixed_cell_qc_thresholds().", call. = FALSE)
  }
  if (length(n_feature) != length(n_count)) stop("n_feature and n_count must have equal length.", call. = FALSE)
  n_feature > fixed_thresholds$feature_lower & n_feature < fixed_thresholds$feature_upper &
    n_count > fixed_thresholds$count_lower & n_count < fixed_thresholds$count_upper
}

#' Calculate per-cell Xenium QC metrics without deleting or modifying raw data.
#'
#' @param counts Sparse gene-by-cell raw-count matrix with cell IDs as columns.
#' @param cells Cell metadata aligned exactly to the count-matrix columns.
#' @param region_id Single Xenium region identifier.
#' @param fixed_thresholds Validated fixed-threshold list from
#'   read_fixed_cell_qc_thresholds().
#' @return List with cell_metadata (metrics and flags), thresholds (auditable
#'   cutoffs/methods), and summary (one-row cell counts for the region).
calculate_xenium_cell_qc <- function(counts, cells, region_id, fixed_thresholds) {
  require_package("Matrix")
  required <- c("cell_id", "total_counts", "control_probe_counts", "genomic_control_counts", "control_codeword_counts", "cell_area", "nucleus_count")
  missing <- setdiff(required, names(cells))
  if (length(missing)) stop(sprintf("Cell metadata missing QC columns: %s", paste(missing, collapse = ",")), call. = FALSE)
  if (!identical(colnames(counts), cells$cell_id)) stop("Count columns and cell metadata are not aligned.", call. = FALSE)
  n_count <- as.numeric(Matrix::colSums(counts))
  n_feature <- as.numeric(Matrix::colSums(counts > 0))
  count_descriptive <- robust_interval(n_count, 3, 5, 1)
  feature_descriptive <- robust_interval(n_feature, 3, 5, 1)
  area_bounds <- robust_interval(cells$cell_area, 5, 5, 0)
  control_count <- cells$control_probe_counts + cells$genomic_control_counts + cells$control_codeword_counts
  control_fraction <- ifelse(cells$total_counts > 0, control_count / cells$total_counts, 0)
  control_upper <- max(0.05, safe_quantile(control_fraction, 0.995))
  qc_core_pass <- apply_fixed_primary_bounds(n_feature, n_count, fixed_thresholds)
  nucleus_missing <- cells$nucleus_count == 0
  multiple_nuclei <- cells$nucleus_count > 1
  area_outlier <- cells$cell_area < as.numeric(area_bounds["lower"]) | cells$cell_area > as.numeric(area_bounds["upper"])
  high_control <- control_fraction > control_upper
  # Robust high-tail limits remain diagnostic evidence for a possible
  # segmentation multiplet. They do not define qc_core_pass or primary_include.
  high_complexity <- n_count > as.numeric(count_descriptive["upper"]) |
    n_feature > as.numeric(feature_descriptive["upper"])
  segmentation_multiplet <- multiple_nuclei | (high_complexity & cells$cell_area > as.numeric(area_bounds["upper"]))
  out <- cells
  out$region_id <- region_id; out$nCount_Xenium <- n_count; out$nFeature_Xenium <- n_feature
  out$control_fraction_cell <- control_fraction; out$nucleus_missing_flag <- nucleus_missing
  out$multiple_nuclei_flag <- multiple_nuclei; out$cell_area_outlier_flag <- area_outlier
  out$high_control_flag <- high_control; out$segmentation_multiplet_flag <- segmentation_multiplet
  out$qc_core_pass <- qc_core_pass
  out$qc_review_flag <- nucleus_missing | segmentation_multiplet | area_outlier | high_control | !qc_core_pass
  thresholds <- rbind(
    data.frame(metric = "nCount_Xenium", lower = fixed_thresholds$count_lower, upper = fixed_thresholds$count_upper, value = stats::median(n_count), method = "fixed_exclusive_user_approved"),
    data.frame(metric = "nFeature_Xenium", lower = fixed_thresholds$feature_lower, upper = fixed_thresholds$feature_upper, value = stats::median(n_feature), method = "fixed_exclusive_user_approved"),
    data.frame(metric = "cell_area", lower = as.numeric(area_bounds["lower"]), upper = as.numeric(area_bounds["upper"]), value = as.numeric(area_bounds["median"]), method = area_bounds["method"]),
    data.frame(metric = "control_fraction_cell", lower = 0, upper = control_upper, value = stats::median(control_fraction), method = "max_5pct_or_q99.5")
  )
  thresholds$region_id <- region_id
  thresholds <- thresholds[, c("region_id", "metric", "lower", "upper", "value", "method")]
  summary <- data.frame(
    region_id = region_id, input_cells = nrow(out), core_qc_pass = sum(out$qc_core_pass),
    core_qc_fail = sum(!out$qc_core_pass), review_flagged = sum(out$qc_review_flag),
    nucleus_missing = sum(out$nucleus_missing_flag), multiple_nuclei = sum(out$multiple_nuclei_flag),
    segmentation_multiplet = sum(out$segmentation_multiplet_flag), area_outlier = sum(out$cell_area_outlier_flag),
    high_control = sum(out$high_control_flag), cells_deleted = 0L, stringsAsFactors = FALSE
  )
  list(cell_metadata = out, thresholds = thresholds, summary = summary)
}

#' Overall readiness.
#'
#' @param status Required `status` input; validated before computation.
#' @return A deterministic derived value or annotated copy; the supplied raw data are not overwritten.
overall_readiness <- function(status) {
  status <- toupper(as.character(status))
  if (any(status %in% c("FAIL", "HOLD", "BLOCKED"))) return("HOLD")
  if (any(status %in% c("PENDING", "WARN"))) return("PENDING")
  "PASS"
}
#' Section palette.
#'
#' @return A deterministic derived value or annotated copy; the supplied raw data are not overwritten.
section_palette <- function() {
  c(Region_1 = "#3C5488", Region_2 = "#00A087", Region_3 = "#E64B35", Region_4 = "#F39B7F")
}

#' Return the canonical ordered scWAT Xenium region identifiers.
#'
#' @return A character vector containing `Region_1` through `Region_4` in
#'   acquisition order.
#' @examples
#' expected_scwat_regions()
expected_scwat_regions <- function() paste0("Region_", seq_len(4L))

#' Validate the index for one coherent four-region initial-QC run.
#'
#' The function rejects missing or duplicated regions and prevents the slide
#' summary from combining section bundles generated under different run labels
#' or execution modes.
#'
#' @param index Data frame with one row per section bundle.
#' @param expected_regions Ordered region identifiers. Defaults to
#'   [expected_scwat_regions()].
#' @return `index`, reordered to `expected_regions`.
#' @export
validate_scwat_section_bundle_index <- function(index, expected_regions = expected_scwat_regions()) {
  required <- c("region_id", "run_label", "execution_mode", "section_output_dir")
  missing <- setdiff(required, names(index))
  if (length(missing)) {
    stop(sprintf("Section bundle index missing columns: %s", paste(missing, collapse = ", ")), call. = FALSE)
  }
  if (nrow(index) != length(expected_regions) ||
      !setequal(index$region_id, expected_regions) || anyDuplicated(index$region_id)) {
    stop(sprintf(
      "Expected exactly four unique scWAT section bundles (%s); duplicates and omissions are invalid.",
      paste(expected_regions, collapse = ", ")
    ), call. = FALSE)
  }
  if (length(unique(index$run_label)) != 1L) {
    stop("All scWAT section bundles must have one run label.", call. = FALSE)
  }
  if (length(unique(index$execution_mode)) != 1L) {
    stop("All scWAT section bundles must have one execution mode.", call. = FALSE)
  }
  index[match(expected_regions, index$region_id), , drop = FALSE]
}

#' Derive four-region initial-QC readiness from current technical gates.
#'
#' This is deliberately independent of the later Region 3 anchor admission
#' workflow. Values such as `PRIMARY`, `PRIMARY_CONDITIONAL`, and
#' `SENSITIVITY_ONLY` are downstream decisions and are not generated here.
#'
#' @param gates Data frame containing `region_id`, `gate`, and `status`.
#' @param expected_regions Ordered region identifiers.
#' @return One row per region with technical status and a logical indicator of
#'   whether biological interpretation is currently allowed.
#' @export
derive_scwat_region_readiness <- function(gates, expected_regions = expected_scwat_regions()) {
  required <- c("region_id", "gate", "status")
  missing <- setdiff(required, names(gates))
  if (length(missing)) {
    stop(sprintf("Readiness gates missing columns: %s", paste(missing, collapse = ", ")), call. = FALSE)
  }
  overall <- gates[gates$gate == "overall", c("region_id", "status"), drop = FALSE]
  if (nrow(overall) != length(expected_regions) || anyDuplicated(overall$region_id) ||
      !setequal(overall$region_id, expected_regions)) {
    stop("Exactly one overall readiness gate is required for each of the four scWAT regions.", call. = FALSE)
  }
  overall <- overall[match(expected_regions, overall$region_id), , drop = FALSE]
  overall$biological_interpretation_allowed <- overall$status == "PASS"
  rownames(overall) <- NULL
  overall
}

#' Summarise the two scWAT sections within each mouse.
#'
#' Cells are never treated as biological replicates. The output provides
#' section-level values and descriptive mouse-level means only. It does not
#' create a left/right or anatomical-position factor because the supplied
#' metadata explicitly defines neither.
#'
#' @param section_summary Region-level QC summary with `region_id`,
#'   `input_cells`, and `core_qc_pass`.
#' @param manifest scWAT sample manifest with `region_id`, `mouse_id`, and
#'   `section_id`.
#' @return A named list containing `region_summary`, `within_mouse`, and
#'   `mouse_summary` data frames.
#' @export
summarise_scwat_mouse_sections <- function(section_summary, manifest) {
  needed_summary <- c("region_id", "input_cells", "core_qc_pass")
  needed_manifest <- c("region_id", "mouse_id", "section_id")
  if (length(setdiff(needed_summary, names(section_summary)))) {
    stop("Section summary lacks region_id/input_cells/core_qc_pass.", call. = FALSE)
  }
  if (length(setdiff(needed_manifest, names(manifest)))) {
    stop("Manifest lacks region_id/mouse_id/section_id.", call. = FALSE)
  }
  regions <- expected_scwat_regions()
  if (!setequal(section_summary$region_id, regions) || !setequal(manifest$region_id, regions)) {
    stop("Section summary and manifest must both contain all four scWAT regions.", call. = FALSE)
  }
  region_summary <- merge(
    manifest[, needed_manifest, drop = FALSE], section_summary,
    by = "region_id", sort = FALSE
  )
  region_summary <- region_summary[match(regions, region_summary$region_id), , drop = FALSE]
  region_summary$core_pass_fraction <- ifelse(
    region_summary$input_cells > 0,
    region_summary$core_qc_pass / region_summary$input_cells,
    NA_real_
  )
  within_mouse <- region_summary[order(region_summary$mouse_id, region_summary$region_id), , drop = FALSE]
  mouse_summary <- stats::aggregate(
    region_summary[, c("input_cells", "core_qc_pass", "core_pass_fraction"), drop = FALSE],
    by = list(mouse_id = region_summary$mouse_id), FUN = mean, na.rm = TRUE
  )
  mouse_summary$n_sections <- as.integer(table(region_summary$mouse_id)[mouse_summary$mouse_id])
  mouse_summary <- mouse_summary[order(mouse_summary$mouse_id), , drop = FALSE]
  rownames(region_summary) <- rownames(within_mouse) <- rownames(mouse_summary) <- NULL
  list(region_summary = region_summary, within_mouse = within_mouse, mouse_summary = mouse_summary)
}

#' Section downstream status.
#'
#' @param region_id Required `region_id` input; validated before computation.
#' @return A deterministic derived value or annotated copy; the supplied raw data are not overwritten.
section_downstream_status <- function(region_id) {
  status <- c(
    Region_1 = "PRIMARY_CONDITIONAL",
    Region_2 = "PRIMARY_CONDITIONAL",
    Region_3 = "PRIMARY",
    Region_4 = "SENSITIVITY_ONLY"
  )
  region_id <- as.character(region_id)
  unknown <- setdiff(unique(region_id), names(status))
  if (length(unknown)) {
    stop(sprintf("Unknown scWAT region: %s", paste(unknown, collapse = ", ")), call. = FALSE)
  }
  unname(status[region_id])
}

#' Build provenance-preserving primary and sensitivity cell masks.
#'
#' @param cell_metadata Per-cell QC metadata from calculate_xenium_cell_qc().
#' @param spatial_hotspots Optional hotspot table with grid_id and
#'   hotspot_status; used only for the Region 3 hotspot sensitivity mask.
#' @param provenance Non-empty run identifier written to every output row.
#' @return Original cell metadata plus primary_include, strict_include,
#'   hotspot_sensitivity_include, section status, rule text, and provenance.
build_cell_downstream_masks <- function(cell_metadata, spatial_hotspots = data.frame(), provenance) {
  required <- c(
    "region_id", "cell_id", "qc_core_pass", "segmentation_multiplet_flag",
    "high_control_flag", "qc_review_flag"
  )
  missing <- setdiff(required, names(cell_metadata))
  if (length(missing)) {
    stop(sprintf("Cell metadata lacks downstream-mask fields: %s", paste(missing, collapse = ", ")), call. = FALSE)
  }
  if (anyDuplicated(paste(cell_metadata$region_id, cell_metadata$cell_id, sep = "::"))) {
    stop("Cell downstream-mask keys must be unique within region.", call. = FALSE)
  }
  if (length(provenance) != 1L || is.na(provenance) || !nzchar(provenance)) {
    stop("A non-empty provenance value is required for downstream masks.", call. = FALSE)
  }
  out <- cell_metadata
  # Primary is the prespecified fixed feature/count cohort. Segmentation,
  # control, nucleus, and area evidence remains visible in review fields and is
  # excluded by strict_include rather than silently redefining primary_include.
  out$primary_include <- as.logical(out$qc_core_pass)
  out$strict_include <- out$primary_include & !as.logical(out$qc_review_flag)
  hotspot_ids <- character()
  if (nrow(spatial_hotspots)) {
    hotspot_required <- c("grid_id", "hotspot_status")
    hotspot_missing <- setdiff(hotspot_required, names(spatial_hotspots))
    if (length(hotspot_missing)) stop("Spatial hotspots lack grid_id or hotspot_status.", call. = FALSE)
    hotspot_ids <- unique(as.character(spatial_hotspots$grid_id[
      spatial_hotspots$hotspot_status == "MORPHOLOGY_REVIEW_REQUIRED"
    ]))
  }
  out$hotspot_review_cell <- FALSE
  if ("grid_id" %in% names(out)) {
    out$hotspot_review_cell <- out$region_id == "Region_3" & as.character(out$grid_id) %in% hotspot_ids
  }
  out$hotspot_sensitivity_include <- out$primary_include & !out$hotspot_review_cell
  out$section_status <- section_downstream_status(out$region_id)
  out$mask_rule_primary <- "nFeature_Xenium > 5 AND nFeature_Xenium < 200 AND nCount_Xenium > 10 AND nCount_Xenium < 1000"
  out$mask_rule_strict <- "primary_include AND NOT qc_review_flag"
  out$mask_rule_hotspot_sensitivity <- "primary_include AND NOT Region_3 morphology-review hotspot cell"
  out$provenance <- provenance
  out
}

#' Cell style theme.
#'
#' @param base_size Optional `base_size` input with the default shown in the function signature.
#' @return A deterministic derived value or annotated copy; the supplied raw data are not overwritten.
cell_style_theme <- function(base_size = 14) {
  require_package("ggplot2")
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      axis.line = ggplot2::element_line(linewidth = 0.35, colour = "#222222"),
      axis.ticks = ggplot2::element_line(linewidth = 0.3, colour = "#222222"),
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 1),
      plot.subtitle = ggplot2::element_text(colour = "#4D4D4D", size = base_size - 1),
      legend.title = ggplot2::element_text(face = "bold"),
      strip.background = ggplot2::element_blank(), strip.text = ggplot2::element_text(face = "bold")
    )
}

#' Return the biological display order used for scWAT cell types
#'
#' @param labels Optional observed labels. When supplied, only observed known
#'   labels are returned and previously unseen labels are appended in their
#'   first-observed order.
#'
#' @return A character vector of cell-type labels in biological display order.
scwat_cell_type_order <- function(labels = NULL) {
  canonical <- c(
    "Adipocyte", "MatureAdip",
    "ASC", "APC", "Fibroblast", "Stromal_Fibroblast",
    "VSMC", "Pericyte", "Mural",
    "Capillary_EC", "Venous_EC", "Lymphatic_EC", "Endothelial",
    "Schwann", "Neural", "Mesothelial", "Epithelial",
    "LYVE1_resident_Mac", "Scavenging_macrophage", "Inflammatory_Mac",
    "TREM2_LAM", "Macrophage", "Monocyte", "Myeloid",
    "DC", "cDC", "cDC1", "cDC2", "CCR7_migratory_DC", "Neutrophil",
    "Eosinophil", "Mast cell",
    "ILC", "ILC2", "NK", "T", "γδ T", "T cell", "B", "B cell",
    "Plasma", "Uncertain"
  )
  if (is.null(labels)) {
    return(canonical)
  }
  observed <- unique(as.character(labels[!is.na(labels) & nzchar(as.character(labels))]))
  c(intersect(canonical, observed), setdiff(observed, canonical))
}

#' Apply the scWAT biological order without changing label values
#'
#' @param labels Character or factor cell-type labels.
#'
#' @return A factor containing the original values and biologically ordered
#'   levels. Missing values remain missing.
apply_scwat_cell_type_order <- function(labels) {
  factor(as.character(labels), levels = scwat_cell_type_order(labels))
}

#' Return stable macaron colours for scWAT categories
#'
#' @param labels Optional observed labels. If omitted, colours for the complete
#'   canonical scWAT order are returned. Unknown labels receive deterministic
#'   fallback colours after the canonical palette.
#'
#' @return A named character vector of six-digit hexadecimal colours.
cell_macaron_palette <- function(labels = NULL) {
  canonical <- scwat_cell_type_order()
  roots <- c(
    "#F2B8A2", "#F6D7A7", "#E7C6A5", "#D7BDE2", "#C9B4D9",
    "#B9CDE5", "#AFC6E9", "#A9D6C8", "#BFD8A8", "#D7E8B2",
    "#B7DDD2", "#9FD3C7", "#A8DADC", "#B8D8E8", "#C7C5E8",
    "#D8C4E6", "#E6C7D5", "#E8B4B8", "#D99C9C", "#E7AAA2",
    "#D6A59A", "#E3B7A0", "#D8B58A", "#DCC98F", "#C6D59B",
    "#B3D1A4", "#A6CDB4", "#A4CEC6", "#B7D9D0", "#F2C48D",
    "#E89A8F", "#E5B2C5", "#C9B6DF", "#B7BEE0", "#A9C5E6",
    "#9FC8D8", "#B5D6C6", "#C7DDAF", "#D9DEA8", "#E7D2A8",
    "#E6C2A5", "#D8D8D8"
  )
  names(roots) <- canonical
  if (is.null(labels)) {
    return(roots)
  }
  observed <- unique(as.character(labels[!is.na(labels) & nzchar(as.character(labels))]))
  output <- roots[intersect(canonical, observed)]
  unknown <- setdiff(observed, canonical)
  if (length(unknown)) {
    fallback <- grDevices::hcl.colors(length(unknown), palette = "Pastel 1")
    names(fallback) <- unknown
    output <- c(output, fallback)
  }
  output[observed]
}

#' Apply the shared cell-style theme to a ggplot-compatible object
#'
#' @param plot A ggplot or patchwork-compatible plot object.
#' @param base_size Base font size.
#' @param legend_position ggplot legend-position setting.
#'
#' @return The styled plot object; the input data are not modified.
style_cell_plot <- function(plot, base_size = 12, legend_position = "right") {
  require_package("ggplot2")
  plot +
    cell_style_theme(base_size = base_size) +
    ggplot2::theme(
      legend.position = legend_position,
      plot.margin = ggplot2::margin(8, 12, 8, 8)
    )
}

#' Plot focal categories above smaller contextual points
#'
#' @param data Data frame containing coordinates and a grouping field.
#' @param x_col,y_col Numeric coordinate column names.
#' @param group_col Categorical grouping column name.
#' @param target_labels Groups to draw in the final, larger point layer.
#' @param palette Named colours covering all observed groups.
#' @param background_size,target_size Point sizes for context and focal layers.
#' @param background_alpha,target_alpha Point opacity for context and focal layers.
#' @param fixed_coordinates Whether to use an equal-aspect coordinate system.
#' @param title Optional title. Subtitles are deliberately unsupported.
#' @param legend_title Optional legend title.
#'
#' @return A cell-style ggplot with contextual points in layer one and target
#'   points in layer two, guaranteeing that focal cells are drawn on top.
plot_target_overlay <- function(
    data,
    x_col,
    y_col,
    group_col,
    target_labels,
    palette = NULL,
    background_size = 0.35,
    target_size = 1.1,
    background_alpha = 0.45,
    target_alpha = 0.95,
    fixed_coordinates = FALSE,
    title = NULL,
    legend_title = NULL
) {
  require_package("ggplot2")
  plot_data <- as.data.frame(data, stringsAsFactors = FALSE)
  required <- c(x_col, y_col, group_col)
  missing <- setdiff(required, names(plot_data))
  if (length(missing)) stop("Overlay plot data missing: ", paste(missing, collapse = ", "), call. = FALSE)
  if (any(!is.finite(as.numeric(plot_data[[x_col]]))) || any(!is.finite(as.numeric(plot_data[[y_col]])))) {
    stop("Overlay plot coordinates must be finite.", call. = FALSE)
  }
  groups <- unique(as.character(plot_data[[group_col]]))
  groups <- groups[!is.na(groups) & nzchar(groups)]
  targets <- intersect(as.character(target_labels), groups)
  if (is.null(palette)) palette <- cell_macaron_palette(groups)
  if (is.null(names(palette)) || !all(groups %in% names(palette))) {
    stop("palette must be named and cover every observed group.", call. = FALSE)
  }
  plot_data[[group_col]] <- factor(as.character(plot_data[[group_col]]), levels = groups)
  is_target <- as.character(plot_data[[group_col]]) %in% targets
  background <- plot_data[!is_target, , drop = FALSE]
  focal <- plot_data[is_target, , drop = FALSE]
  plot <- ggplot2::ggplot() +
    ggplot2::geom_point(
      data = background,
      ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]], colour = .data[[group_col]]),
      size = background_size, alpha = background_alpha
    ) +
    ggplot2::geom_point(
      data = focal,
      ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]], colour = .data[[group_col]]),
      size = target_size, alpha = target_alpha
    ) +
    ggplot2::scale_colour_manual(values = palette[groups], drop = FALSE) +
    ggplot2::labs(title = title, x = NULL, y = NULL, colour = legend_title)
  if (isTRUE(fixed_coordinates)) plot <- plot + ggplot2::coord_fixed()
  style_cell_plot(plot)
}

#' Save a cell-style plot as matched PNG and PDF artifacts
#'
#' @param plot A ggplot or patchwork-compatible plot object.
#' @param stem Filename stem without an extension.
#' @param output_dir Existing or creatable output directory below
#'   `project_root`.
#' @param project_root Absolute project root used for path-safety validation.
#' @param width,height Plot dimensions in inches.
#' @param dpi PNG resolution in dots per inch.
#'
#' @return A named character vector containing the validated PNG and PDF paths.
#'   The supplied plot and analysis objects are not modified.
save_cell_plot <- function(
    plot,
    stem,
    output_dir,
    project_root,
    width = 12,
    height = 7,
    dpi = 300
) {
  require_package("ggplot2")
  stem <- as.character(stem)[[1L]]
  if (is.na(stem) || !grepl("^[A-Za-z0-9_.-]+$", stem)) {
    stop("stem must contain only letters, numbers, dot, underscore or hyphen.", call. = FALSE)
  }
  paths <- c(
    png = file.path(output_dir, paste0(stem, ".png")),
    pdf = file.path(output_dir, paste0(stem, ".pdf"))
  )
  invisible(lapply(paths, function(path) assert_path_within(project_root, path)))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(paths[["png"]], plot = plot, width = width, height = height, dpi = dpi, bg = "white")
  ggplot2::ggsave(paths[["pdf"]], plot = plot, width = width, height = height, device = grDevices::cairo_pdf, bg = "white")
  paths
}

#' Order canonical-marker rows and build matching DotPlot feature groups
#'
#' @param marker_df Data frame containing `Gene_Symbol` and
#'   `CellType_subtype`.
#' @param available_genes Character vector of genes present in the assay.
#'
#' @return A list containing the ordered marker table, a named feature list and
#'   the observed biological cell-type order.
order_marker_features <- function(marker_df, available_genes) {
  required <- c("Gene_Symbol", "CellType_subtype")
  missing <- setdiff(required, names(marker_df))
  if (length(missing)) {
    stop("marker_df missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  marker_table <- as.data.frame(marker_df, stringsAsFactors = FALSE)
  marker_table <- marker_table[
    marker_table$Gene_Symbol %in% as.character(available_genes),
    , drop = FALSE
  ]
  marker_table <- marker_table[!duplicated(marker_table$Gene_Symbol), , drop = FALSE]
  subtype_order <- scwat_cell_type_order(marker_table$CellType_subtype)
  marker_table$CellType_subtype <- factor(
    as.character(marker_table$CellType_subtype),
    levels = subtype_order
  )
  marker_table <- marker_table[
    order(marker_table$CellType_subtype, seq_len(nrow(marker_table))),
    , drop = FALSE
  ]
  rownames(marker_table) <- NULL
  feature_groups <- split(
    as.character(marker_table$Gene_Symbol),
    factor(as.character(marker_table$CellType_subtype), levels = subtype_order),
    drop = TRUE
  )
  list(
    marker_table = marker_table,
    feature_groups = feature_groups,
    cell_type_order = subtype_order
  )
}

#' Plot section qc.
#'
#' @param cell_metadata Required `cell_metadata` input; validated before computation.
#' @param region_id Required `region_id` input; validated before computation.
#' @return A `ggplot` object or named list of plots; input objects are not modified.
plot_section_qc <- function(cell_metadata, region_id) {
  require_package("ggplot2")
  colour <- unname(section_palette()[region_id])
  if (is.na(colour)) stop(sprintf("No section color configured for %s.", region_id), call. = FALSE)
  common_hist <- function(metric, label, title) {
    ggplot2::ggplot(cell_metadata, ggplot2::aes(x = .data[[metric]])) +
      ggplot2::geom_histogram(bins = 45, fill = colour, colour = "white", linewidth = 0.15) +
      ggplot2::labs(title = title, subtitle = region_id, x = label, y = "Cells") + cell_style_theme()
  }
  status <- factor(ifelse(cell_metadata$qc_review_flag, "Review", "Pass"), levels = c("Pass", "Review"))
  spatial_data <- transform(cell_metadata, QC_status = status)
  list(
    counts = common_hist("nCount_Xenium", "Gene-expression transcripts per cell", "Transcript-count distribution"),
    features = common_hist("nFeature_Xenium", "Genes detected per cell", "Detected-feature distribution"),
    area = common_hist("cell_area", expression("Cell area ("*mu*"m"^2*")"), "Cell-area distribution"),
    spatial = ggplot2::ggplot(
  spatial_data,
  ggplot2::aes(x = x_centroid, y = y_centroid)
) +
  ggplot2::geom_point(
    data = spatial_data[spatial_data$QC_status == "Pass", ],
    ggplot2::aes(colour = QC_status),
    size = 0.35,
    alpha = 0.75
  ) +
  ggplot2::geom_point(
    data = spatial_data[spatial_data$QC_status == "Review", ],
    ggplot2::aes(colour = QC_status),
    size = 0.35,
    alpha = 0.9
  ) +
  ggplot2::scale_colour_manual(
    values = c(
      Pass = "#BDBDBD",
      Review = "#D73027"
    ),
    drop = FALSE
  ) +
  ggplot2::coord_fixed() +
  ggplot2::labs(
    title = "Spatial QC review map",
    subtitle = region_id,
    x = "X centroid",
    y = "Y centroid",
    colour = "QC"
  ) +
  cell_style_theme()
  )
}

#' Save section plots.
#'
#' @param plots Required `plots` input; validated before computation.
#' @param figure_dir Required `figure_dir` input; validated before computation.
#' @param region_id Required `region_id` input; validated before computation.
#' @param project_root Required `project_root` input; validated before computation.
#' @return Written artifact path(s) or an auditable one-row write summary.
save_section_plots <- function(plots, figure_dir, region_id, project_root) {
  require_package("ggplot2")
  assert_path_within(project_root, figure_dir)
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  pdf_path <- file.path(figure_dir, paste0(region_id, "_qc_overview.pdf"))
  grDevices::pdf(pdf_path, width = 7, height = 5, onefile = TRUE)
  for (plot in plots) print(plot)
  grDevices::dev.off()
  png_paths <- vapply(names(plots), function(name) {
    path <- file.path(figure_dir, paste0(region_id, "_", name, ".png"))
    ggplot2::ggsave(path, plots[[name]], width = 7, height = 5, units = "in", dpi = 300, bg = "white")
    path
  }, character(1))
  c(pdf_path, unname(png_paths))
}

#' Write gz tsv.
#'
#' @param x Required `x` input; validated before computation.
#' @param path Required `path` input; validated before computation.
#' @param project_root Required `project_root` input; validated before computation.
#' @return Written artifact path(s) or an auditable one-row write summary.
write_gz_tsv <- function(x, path, project_root) {
  assert_path_within(project_root, path)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  con <- gzfile(path, "wt"); on.exit(close(con), add = TRUE)
  table <- flatten_tsv_table(x)
  tryCatch(
    utils::write.table(table, con, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"),
    error = function(error) {
      stop("Failed to write gzipped TSV '", path, "': ", conditionMessage(error), call. = FALSE)
    }
  )
  invisible(path)
}

#' Write one reloadable, non-destructive scWAT region QC bundle.
#'
#' The raw sparse count matrix and every segmented cell are retained. Cell
#' filtering is represented only by `primary_include`, `strict_include`, and
#' `hotspot_sensitivity_include` columns. Initial technical readiness is saved
#' separately from later Region 3 anchor and mapping decisions.
#'
#' @param project_root Absolute Xenium project root. Every output path must be
#'   below this directory.
#' @param output_dir Region output directory below `adipose_analysis_B2`.
#' @param region_id One value from [expected_scwat_regions()].
#' @param run_label Non-empty identifier shared by all four regions.
#' @param execution_mode Either `FULL_HPC` or `LOCAL_SUBSET`.
#' @param manifest Four-row scWAT sample manifest.
#' @param inventory File-inventory data frame from [inventory_section_files()].
#' @param integrity Matrix-alignment result from [validate_section_integrity()].
#' @param bundle Imported Xenium list containing sparse `counts`, `cells`, and
#'   `features`.
#' @param cell_qc Result from [calculate_xenium_cell_qc()].
#' @param masks Cell table returned by [build_cell_downstream_masks()].
#' @return A one-row data frame reporting the output directory and retained
#'   cell counts. The function also writes tables, plots, session information,
#'   and `<region_id>.phase0_2_qc.rds`.
#' @export
write_scwat_region_qc_bundle <- function(project_root, output_dir, region_id, run_label,
                                         execution_mode, manifest, inventory, integrity,
                                         bundle, cell_qc, masks) {
  require_package("Matrix")
  assert_path_within(project_root, output_dir)
  normalized_output <- canonical_path(output_dir)
  if (!grepl("/adipose_analysis_B2/", normalized_output, fixed = TRUE)) {
    stop("scWAT QC outputs must be below adipose_analysis_B2.", call. = FALSE)
  }
  if (grepl("colon_analysis", normalized_output, fixed = TRUE)) {
    stop("scWAT QC outputs cannot use colon_analysis.", call. = FALSE)
  }
  if (!region_id %in% expected_scwat_regions()) {
    stop("Invalid scWAT region identifier.", call. = FALSE)
  }
  if (!inherits(bundle$counts, "sparseMatrix") || ncol(bundle$counts) != nrow(masks)) {
    stop("Sparse counts and cell masks are not aligned.", call. = FALSE)
  }
  if (!identical(colnames(bundle$counts), as.character(masks$cell_id))) {
    stop("Sparse count columns and mask cell IDs differ.", call. = FALSE)
  }
  if (!identical(as.character(manifest$region_id), expected_scwat_regions())) {
    stop("The complete ordered four-region scWAT manifest is required.", call. = FALSE)
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  masks$initial_qc_status <- "QC_EVIDENCE_PENDING"
  configuration <- data.frame(
    region_id = region_id,
    run_label = run_label,
    execution_mode = execution_mode,
    primary_rule = unique(masks$mask_rule_primary)[[1L]],
    cells_deleted = 0L,
    stringsAsFactors = FALSE
  )
  gates <- data.frame(
    region_id = region_id,
    gate = c("required_files", "matrix_integrity", "fixed_primary_mask", "metadata", "overall"),
    status = c(
      if (nrow(inventory) && all(inventory$exists)) "PASS" else "HOLD",
      if (nrow(integrity) == 1L && isTRUE(integrity$dimension_match[[1L]])) "PASS" else "HOLD",
      if (!anyNA(masks$primary_include)) "PASS" else "HOLD",
      "PASS",
      "PENDING"
    ),
    details = c(
      "All minimum region inputs exist",
      "Sparse matrix/features/barcodes/cells align",
      "Fixed exclusive primary mask contains no missing values",
      "Four-row verified scWAT manifest loaded",
      "Derived from current technical gates"
    ),
    stringsAsFactors = FALSE
  )
  gates$status[gates$gate == "overall"] <- overall_readiness(gates$status[gates$gate != "overall"])
  alarms <- empty_alarm_table()
  if (!"region_id" %in% names(alarms)) alarms$region_id <- character(nrow(alarms))
  for (field in c("run_label", "execution_mode")) {
    cell_qc$summary[[field]] <- configuration[[field]][[1L]]
    cell_qc$thresholds[[field]] <- configuration[[field]][[1L]]
    gates[[field]] <- configuration[[field]][[1L]]
  }

  write_tsv(configuration, file.path(output_dir, "run_configuration.tsv"), project_root)
  write_tsv(manifest[manifest$region_id == region_id, , drop = FALSE], file.path(output_dir, "sample_metadata.tsv"), project_root)
  write_tsv(inventory, file.path(output_dir, "input_inventory.tsv"), project_root)
  write_tsv(integrity, file.path(output_dir, "input_integrity.tsv"), project_root)
  write_tsv(cell_qc$summary, file.path(output_dir, "qc_summary.tsv"), project_root)
  write_tsv(cell_qc$thresholds, file.path(output_dir, "qc_thresholds.tsv"), project_root)
  write_tsv(gates, file.path(output_dir, "section_readiness_gates.tsv"), project_root)
  write_tsv(alarms, file.path(output_dir, "analysis_alerts.tsv"), project_root)
  write_gz_tsv(masks, file.path(output_dir, "cell_qc_metadata.tsv.gz"), project_root)

  qc_object <- list(
    region_id = region_id,
    run_label = run_label,
    execution_mode = execution_mode,
    counts = bundle$counts,
    cells = masks,
    features = bundle$features,
    thresholds = cell_qc$thresholds,
    summary = cell_qc$summary,
    gates = gates,
    raw_counts_preserved = TRUE,
    cells_deleted = 0L
  )
  saveRDS(qc_object, file.path(output_dir, paste0(region_id, ".phase0_2_qc.rds")), compress = "xz")
  save_section_plots(plot_section_qc(masks, region_id), file.path(output_dir, "figures"), region_id, project_root)
  writeLines(capture.output(sessionInfo()), file.path(output_dir, "session_info.txt"), useBytes = TRUE)

  data.frame(
    region_id = region_id,
    output_dir = normalized_output,
    cells = nrow(masks),
    primary_include = sum(masks$primary_include),
    cells_deleted = 0L,
    stringsAsFactors = FALSE
  )
}

#' Validate and reload one scWAT region QC bundle.
#'
#' @param output_dir Region output directory written by
#'   [write_scwat_region_qc_bundle()].
#' @param region_id Expected region identifier.
#' @param stop_on_error Stop with a diagnostic message when validation fails.
#' @return `TRUE` when every artifact exists and the sparse matrix, cell table,
#'   and mask columns reload in identical cell order; otherwise `FALSE` when
#'   `stop_on_error = FALSE`.
#' @export
validate_scwat_region_qc_bundle <- function(output_dir, region_id, stop_on_error = TRUE) {
  required <- c(
    "run_configuration.tsv", "sample_metadata.tsv", "input_inventory.tsv", "input_integrity.tsv",
    "qc_summary.tsv", "qc_thresholds.tsv", "section_readiness_gates.tsv", "analysis_alerts.tsv",
    "cell_qc_metadata.tsv.gz", paste0(region_id, ".phase0_2_qc.rds"), "session_info.txt",
    file.path("figures", paste0(region_id, "_qc_overview.pdf")),
    file.path("figures", paste0(region_id, "_counts.png")),
    file.path("figures", paste0(region_id, "_features.png")),
    file.path("figures", paste0(region_id, "_area.png")),
    file.path("figures", paste0(region_id, "_spatial.png"))
  )
  paths <- file.path(output_dir, required)
  valid <- all(file.exists(paths))
  message <- if (valid) "PASS" else sprintf(
    "Missing artifacts: %s", paste(basename(paths[!file.exists(paths)]), collapse = ", ")
  )
  if (valid) {
    object <- tryCatch(readRDS(file.path(output_dir, paste0(region_id, ".phase0_2_qc.rds"))), error = identity)
    cells <- tryCatch(utils::read.delim(gzfile(file.path(output_dir, "cell_qc_metadata.tsv.gz")), check.names = FALSE), error = identity)
    valid <- !inherits(object, "error") && !inherits(cells, "error") &&
      inherits(object$counts, "sparseMatrix") && isTRUE(object$raw_counts_preserved) &&
      identical(object$cells_deleted, 0L) && ncol(object$counts) == nrow(cells) &&
      identical(colnames(object$counts), as.character(cells$cell_id)) &&
      !anyNA(cells[, c("primary_include", "strict_include", "hotspot_sensitivity_include"), drop = FALSE])
    if (!valid) message <- "Saved sparse object or cell mask failed reload/alignment validation."
  }
  if (!valid && isTRUE(stop_on_error)) stop(message, call. = FALSE)
  valid
}

#' Read four coherent region bundles for the scWAT slide summary.
#'
#' @param run_root Run directory containing `sections/Region_1` through
#'   `sections/Region_4`.
#' @param expected_regions Ordered region identifiers.
#' @param run_label Optional required run label.
#' @param execution_mode Optional required execution mode.
#' @return A list of coverage, QC summaries, thresholds, readiness gates,
#'   alarms, and combined per-cell metadata.
#' @export
read_scwat_slide_qc_outputs <- function(run_root, expected_regions = expected_scwat_regions(),
                                        run_label = NULL, execution_mode = NULL) {
  section_dirs <- file.path(run_root, "sections", expected_regions)
  index_rows <- lapply(seq_along(expected_regions), function(i) {
    region <- expected_regions[[i]]
    output_dir <- section_dirs[[i]]
    validate_scwat_region_qc_bundle(output_dir, region)
    config <- utils::read.delim(file.path(output_dir, "run_configuration.tsv"), check.names = FALSE)
    data.frame(
      region_id = region,
      run_label = config$run_label[[1L]],
      execution_mode = config$execution_mode[[1L]],
      section_output_dir = output_dir,
      stringsAsFactors = FALSE
    )
  })
  coverage <- validate_scwat_section_bundle_index(do.call(rbind, index_rows), expected_regions)
  if (!is.null(run_label) && !identical(unique(coverage$run_label), run_label)) {
    stop("Section bundle run label does not match requested run label.", call. = FALSE)
  }
  if (!is.null(execution_mode) && !identical(unique(coverage$execution_mode), execution_mode)) {
    stop("Section bundle execution mode does not match requested execution mode.", call. = FALSE)
  }
  read_one <- function(filename, gzipped = FALSE) {
    do.call(rbind, lapply(seq_len(nrow(coverage)), function(i) {
      path <- file.path(coverage$section_output_dir[[i]], filename)
      table <- if (gzipped) {
        utils::read.delim(gzfile(path), check.names = FALSE)
      } else {
        utils::read.delim(path, check.names = FALSE)
      }
      if (!"region_id" %in% names(table)) table$region_id <- coverage$region_id[[i]]
      table
    }))
  }
  list(
    coverage = coverage,
    qc_summary = read_one("qc_summary.tsv"),
    thresholds = read_one("qc_thresholds.tsv"),
    gates = read_one("section_readiness_gates.tsv"),
    alarms = read_one("analysis_alerts.tsv"),
    cell_metadata = read_one("cell_qc_metadata.tsv.gz", TRUE)
  )
}

#' Write the reloadable four-region scWAT slide QC summary bundle.
#'
#' @param project_root Absolute Xenium project root.
#' @param output_dir Slide-summary directory below `adipose_analysis_B2`.
#' @param run_label Run label shared by the four section bundles.
#' @param execution_mode Execution mode shared by the four section bundles.
#' @param slide_data Result from [read_scwat_slide_qc_outputs()].
#' @param slide_summary Result from [summarise_slide_qc()].
#' @param mouse_summary Result from [summarise_scwat_mouse_sections()].
#' @return A one-row data frame reporting run identity, number of regions, and
#'   overall technical status. Tables, figures, session information, and a
#'   reloadable `scwat_slide_qc.rds` object are written to `output_dir`.
#' @export
write_scwat_slide_qc_bundle <- function(project_root, output_dir, run_label, execution_mode,
                                        slide_data, slide_summary, mouse_summary) {
  assert_path_within(project_root, output_dir)
  if (!grepl("/adipose_analysis_B2/", canonical_path(output_dir), fixed = TRUE)) {
    stop("scWAT slide output must be below adipose_analysis_B2.", call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  write_tsv(slide_data$coverage, file.path(output_dir, "section_coverage.tsv"), project_root)
  write_tsv(slide_summary$section_summary, file.path(output_dir, "section_qc_summary.tsv"), project_root)
  write_tsv(slide_summary$readiness, file.path(output_dir, "region_readiness.tsv"), project_root)
  write_tsv(mouse_summary$within_mouse, file.path(output_dir, "within_mouse_section_summary.tsv"), project_root)
  write_tsv(mouse_summary$mouse_summary, file.path(output_dir, "mouse_summary.tsv"), project_root)
  object <- list(
    run_label = run_label,
    execution_mode = execution_mode,
    slide_data = slide_data,
    slide_summary = slide_summary,
    mouse_summary = mouse_summary
  )
  saveRDS(object, file.path(output_dir, "scwat_slide_qc.rds"), compress = "xz")
  slide_plots <- plot_slide_qc(slide_data, slide_summary)
  figure_dir <- file.path(output_dir, "figures")
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  pdf_path <- file.path(figure_dir, "scwat_slide_qc_figures.pdf")
  grDevices::pdf(pdf_path, width = 8, height = 5.5, onefile = TRUE)
  for (plot in slide_plots) print(plot)
  grDevices::dev.off()
  invisible(vapply(names(slide_plots), function(name) {
    path <- file.path(figure_dir, paste0("slide_", name, ".png"))
    ggplot2::ggsave(path, slide_plots[[name]], width = 8, height = 5.5, units = "in", dpi = 300, bg = "white")
    path
  }, character(1)))
  writeLines(capture.output(sessionInfo()), file.path(output_dir, "session_info.txt"), useBytes = TRUE)
  data.frame(
    run_label = run_label,
    execution_mode = execution_mode,
    regions = nrow(slide_data$coverage),
    overall_status = slide_summary$overall_status,
    stringsAsFactors = FALSE
  )
}

#' Validate and reload the scWAT slide QC summary bundle.
#'
#' @param output_dir Directory written by [write_scwat_slide_qc_bundle()].
#' @param expected_regions Ordered region identifiers.
#' @param stop_on_error Stop when any artifact is missing or inconsistent.
#' @return `TRUE` for a complete, reloadable bundle; otherwise `FALSE` when
#'   `stop_on_error = FALSE`.
#' @export
validate_scwat_slide_qc_bundle <- function(output_dir, expected_regions = expected_scwat_regions(),
                                           stop_on_error = TRUE) {
  required <- c(
    "section_coverage.tsv", "section_qc_summary.tsv", "region_readiness.tsv",
    "within_mouse_section_summary.tsv", "mouse_summary.tsv", "scwat_slide_qc.rds",
    "session_info.txt", file.path("figures", "scwat_slide_qc_figures.pdf"),
    file.path("figures", paste0("slide_", c("cell_yield", "counts", "features", "review", "flags", "thresholds"), ".png"))
  )
  paths <- file.path(output_dir, required)
  valid <- all(file.exists(paths))
  if (valid) {
    object <- tryCatch(readRDS(file.path(output_dir, "scwat_slide_qc.rds")), error = identity)
    valid <- !inherits(object, "error") &&
      identical(as.character(object$slide_data$coverage$region_id), expected_regions)
  }
  if (!valid && isTRUE(stop_on_error)) {
    stop("scWAT slide QC bundle failed artifact or reload validation.", call. = FALSE)
  }
  valid
}
#' Add evidence provenance.
#'
#' @param table Required `table` input; validated before computation.
#' @param run_label Required `run_label` input; validated before computation.
#' @param execution_mode Required `execution_mode` input; validated before computation.
#' @param provenance Required `provenance` input; validated before computation.
#' @param source_artifact Required `source_artifact` input; validated before computation.
#' @param generated_utc Required `generated_utc` input; validated before computation.
#' @return A deterministic derived value or annotated copy; the supplied raw data are not overwritten.
add_evidence_provenance <- function(table, run_label, execution_mode, provenance,
                                    source_artifact, generated_utc) {
  table$run_label <- run_label
  table$execution_mode <- toupper(execution_mode)
  table$generated_utc <- generated_utc
  table$source_artifact <- source_artifact
  table$provenance <- provenance
  table
}

#' Build one section downstream decision.
#'
#' @param masks Required `masks` input; validated before computation.
#' @param run_label Required `run_label` input; validated before computation.
#' @param execution_mode Required `execution_mode` input; validated before computation.
#' @param provenance Required `provenance` input; validated before computation.
#' @param generated_utc Required `generated_utc` input; validated before computation.
#' @return A deterministic derived value or annotated copy; the supplied raw data are not overwritten.
build_one_section_downstream_decision <- function(masks, run_label, execution_mode,
                                                  provenance, generated_utc) {
  regions <- unique(as.character(masks$region_id))
  if (length(regions) != 1L) stop("One-section downstream decision requires exactly one region.", call. = FALSE)
  region <- regions[[1]]
  x <- masks
  out <- data.frame(
    region_id = region, section_status = section_downstream_status(region),
    input_cells = nrow(x), primary_include_cells = sum(x$primary_include),
    primary_exclude_cells = sum(!x$primary_include), strict_include_cells = sum(x$strict_include),
    strict_exclude_cells = sum(!x$strict_include),
    hotspot_sensitivity_include_cells = sum(x$hotspot_sensitivity_include),
    hotspot_review_cells = sum(x$hotspot_review_cell),
    cluster_discovery_eligible = region %in% c("Region_1", "Region_2", "Region_3"),
    primary_gene_result_eligible = region %in% c("Region_1", "Region_2", "Region_3"),
    region4_mapping_rule = if (region == "Region_4") "MAP_TO_FINAL_REGION_1_3_REFERENCE; LOW_CONFIDENCE=Uncertain" else "NOT_APPLICABLE",
    primary_mask_rule = "nFeature_Xenium > 5 AND nFeature_Xenium < 200 AND nCount_Xenium > 10 AND nCount_Xenium < 1000",
    strict_mask_rule = "primary_include AND NOT qc_review_flag",
    hotspot_sensitivity_mask_rule = "primary_include AND NOT Region_3 morphology-review hotspot cell",
    decision_rule = "Fixed evidence-only section decision approved 2026-08-15; masks annotate without deleting raw cells",
    stringsAsFactors = FALSE
  )
  add_evidence_provenance(
    out, run_label, execution_mode, provenance,
    "section cell and spatial QC artifacts", generated_utc
  )
}
#' Summarise slide qc.
#'
#' @param slide_data Required `slide_data` input; validated before computation.
#' @return Computed evidence or annotations aligned to the supplied rows/cells; raw counts are unchanged.
summarise_slide_qc <- function(slide_data) {
  summary <- slide_data$qc_summary
  summary$core_pass_fraction <- ifelse(summary$input_cells > 0, summary$core_qc_pass / summary$input_cells, NA_real_)
  summary$review_fraction <- ifelse(summary$input_cells > 0, summary$review_flagged / summary$input_cells, NA_real_)
  overall_rows <- slide_data$gates[slide_data$gates$gate == "overall", , drop = FALSE]
  overall_status <- if (nrow(overall_rows)) overall_readiness(overall_rows$status) else "HOLD"
  readiness <- derive_scwat_region_readiness(slide_data$gates, expected_scwat_regions())
  list(section_summary = summary, readiness = readiness, overall_status = overall_status)
}

#' Plot slide qc.
#'
#' @param slide_data Required `slide_data` input; validated before computation.
#' @param slide_summary Required `slide_summary` input; validated before computation.
#' @return A `ggplot` object or named list of plots; input objects are not modified.
plot_slide_qc <- function(slide_data, slide_summary) {
  require_package("ggplot2")
  palette <- section_palette()
  cells <- slide_data$cell_metadata
  cells$region_id <- factor(cells$region_id, levels = names(palette))
  summary <- slide_summary$section_summary
  summary$region_id <- factor(summary$region_id, levels = names(palette))
  region_scale <- ggplot2::scale_fill_manual(values = palette, drop = FALSE)
  colour_scale <- ggplot2::scale_colour_manual(values = palette, drop = FALSE)
  distribution_plot <- function(metric, label, title) {
    ggplot2::ggplot(cells, ggplot2::aes(x = region_id, y = .data[[metric]], fill = region_id)) +
      ggplot2::geom_violin(scale = "width", trim = TRUE, linewidth = 0.25, alpha = 0.85) +
      ggplot2::geom_boxplot(width = 0.14, outlier.shape = NA, fill = "white", linewidth = 0.3) +
      region_scale + ggplot2::labs(title = title, x = NULL, y = label) + cell_style_theme() + ggplot2::theme(legend.position = "none")
  }
  review_long <- rbind(
    data.frame(region_id = summary$region_id, status = "Core pass", cells = summary$core_qc_pass),
    data.frame(region_id = summary$region_id, status = "Core fail", cells = summary$core_qc_fail)
  )
  flag_names <- c("nucleus_missing", "multiple_nuclei", "segmentation_multiplet", "area_outlier", "high_control")
  flag_long <- do.call(rbind, lapply(flag_names, function(name) data.frame(region_id = summary$region_id, flag = name, cells = summary[[name]])))
  thresholds <- slide_data$thresholds
  thresholds$region_id <- factor(thresholds$region_id, levels = names(palette))
  list(
    cell_yield = ggplot2::ggplot(summary, ggplot2::aes(region_id, input_cells, fill = region_id)) + ggplot2::geom_col(width = 0.7) + region_scale + ggplot2::labs(title = "Segmented cell yield", x = NULL, y = "Cells") + cell_style_theme() + ggplot2::theme(legend.position = "none"),
    counts = distribution_plot("nCount_Xenium", "Gene-expression transcripts per cell", "Transcript counts by section"),
    features = distribution_plot("nFeature_Xenium", "Genes detected per cell", "Detected features by section"),
    review = ggplot2::ggplot(review_long, ggplot2::aes(region_id, cells, fill = status)) + ggplot2::geom_col(position = "fill", width = 0.7) + ggplot2::scale_fill_manual(values = c("Core pass" = "#4DAF4A", "Core fail" = "#D73027")) + ggplot2::scale_y_continuous(labels = function(x) paste0(round(100*x), "%")) + ggplot2::labs(title = "Core QC outcome", x = NULL, y = "Cells", fill = NULL) + cell_style_theme(),
    flags = ggplot2::ggplot(flag_long, ggplot2::aes(region_id, cells, fill = region_id)) + ggplot2::geom_col(width = 0.7) + region_scale + ggplot2::facet_wrap(~flag, scales = "free_y", ncol = 3) + ggplot2::labs(title = "QC review flags", x = NULL, y = "Flagged cells") + cell_style_theme() + ggplot2::theme(legend.position = "none"),
    thresholds = ggplot2::ggplot(thresholds, ggplot2::aes(region_id, value, colour = region_id)) + ggplot2::geom_errorbar(ggplot2::aes(ymin = lower, ymax = upper), width = 0.15, linewidth = 0.45) + ggplot2::geom_point(size = 2) + colour_scale + ggplot2::facet_wrap(~metric, scales = "free_y", ncol = 2) + ggplot2::labs(title = "Section-specific QC intervals", x = NULL, y = "Median and QC interval") + cell_style_theme() + ggplot2::theme(legend.position = "none")
  )
}
#' Region bundle to spatial seurat.
#'
#' @param region_data Required `region_data` input; validated before computation.
#' @param xenium_dir Optional `xenium_dir` input with the default shown in the function signature.
#' @param mask Optional `mask` input with the default shown in the function signature.
#' @param genes Optional `genes` input with the default shown in the function signature.
#' @param project Optional `project` input with the default shown in the function signature.
#' @param assay Optional `assay` input with the default shown in the function signature.
#' @param fov Optional `fov` input with the default shown in the function signature.
#' @param include_cell_segmentation Optional `include_cell_segmentation` input with the default shown in the function signature.
#' @param include_nucleus_segmentation Optional `include_nucleus_segmentation` input with the default shown in the function signature.
#' @return A deterministic derived value or annotated copy; the supplied raw data are not overwritten.
region_bundle_to_spatial_seurat <- function(
    region_data,
    xenium_dir = NULL,
    mask = NULL,
    genes = NULL,
    project = NULL,
    assay = "Xenium",
    fov = "fov",
    include_cell_segmentation = TRUE,
    include_nucleus_segmentation = TRUE
) {

  # ============================================================
  # 1. Validate downstream bundle
  # ============================================================

  required <- c(
    "counts",
    "cell_metadata",
    "region_id",
    "gene_sets"
  )

  missing <- setdiff(
    required,
    names(region_data)
  )

  if (length(missing)) {
    stop(
      "region_data missing required elements: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  counts <- region_data$counts
  metadata <- region_data$cell_metadata

  if (!inherits(counts, "sparseMatrix")) {
    stop(
      "region_data$counts must be a sparse Matrix.",
      call. = FALSE
    )
  }

  if (
    is.null(rownames(counts)) ||
    is.null(colnames(counts))
  ) {
    stop(
      "region_data$counts must have gene rownames and cell colnames.",
      call. = FALSE
    )
  }

  if (anyDuplicated(rownames(counts))) {
    stop(
      "Count matrix contains duplicated gene names.",
      call. = FALSE
    )
  }

  if (anyDuplicated(colnames(counts))) {
    stop(
      "Count matrix contains duplicated cell IDs.",
      call. = FALSE
    )
  }


  # ------------------------------------------------------------
  # Required spatial metadata
  # ------------------------------------------------------------

  required_meta <- c(
    "cell_id",
    "x_centroid",
    "y_centroid"
  )

  missing_meta <- setdiff(
    required_meta,
    names(metadata)
  )

  if (length(missing_meta)) {
    stop(
      "Spatial metadata missing required columns: ",
      paste(missing_meta, collapse = ", "),
      call. = FALSE
    )
  }

  if (anyDuplicated(metadata$cell_id)) {
    stop(
      "region_data$cell_metadata contains duplicated cell_id values.",
      call. = FALSE
    )
  }

  if (
    any(!is.finite(metadata$x_centroid)) ||
    any(!is.finite(metadata$y_centroid))
  ) {
    stop(
      "Spatial centroid coordinates contain non-finite values.",
      call. = FALSE
    )
  }


  # ============================================================
  # 2. Cell selection
  # ============================================================

  if (!is.null(mask)) {

    if (
      length(mask) != 1L ||
      !mask %in% names(metadata)
    ) {
      stop(
        "Unknown or invalid cell mask: ",
        paste(mask, collapse = ", "),
        call. = FALSE
      )
    }

    keep <- as.logical(
      metadata[[mask]]
    )

    if (anyNA(keep)) {
      stop(
        "Cell mask contains NA values: ",
        mask,
        call. = FALSE
      )
    }

    metadata <- metadata[
      keep,
      ,
      drop = FALSE
    ]
  }


  # ------------------------------------------------------------
  # Cell IDs
  # ------------------------------------------------------------

  cell_ids <- as.character(
    metadata$cell_id
  )

  if (!length(cell_ids)) {
    stop(
      "No cells remain after cell selection.",
      call. = FALSE
    )
  }

  missing_cells <- setdiff(
    cell_ids,
    colnames(counts)
  )

  if (length(missing_cells)) {
    stop(
      length(missing_cells),
      " selected cells are absent from the count matrix.",
      call. = FALSE
    )
  }


  # ============================================================
  # 3. Gene selection
  # ============================================================

  if (is.null(genes)) {

    genes <- rownames(counts)

  } else if (
    length(genes) == 1L &&
    genes %in% names(region_data$gene_sets)
  ) {

    genes <- region_data$gene_sets[[genes]]

  }


  genes <- unique(
    as.character(genes)
  )

  genes <- genes[
    !is.na(genes) &
      nzchar(genes)
  ]

  if (!length(genes)) {
    stop(
      "No genes remain after gene selection.",
      call. = FALSE
    )
  }

  missing_genes <- setdiff(
    genes,
    rownames(counts)
  )

  if (length(missing_genes)) {
    stop(
      length(missing_genes),
      " requested genes are absent from the downstream count matrix. ",
      "Examples: ",
      paste(
        head(missing_genes, 10),
        collapse = ", "
      ),
      call. = FALSE
    )
  }


  # ============================================================
  # 4. Subset raw counts
  # ============================================================

  counts <- counts[
    genes,
    cell_ids,
    drop = FALSE
  ]


  # ============================================================
  # 5. Align metadata exactly to count matrix
  # ============================================================

  metadata <- metadata[
    match(
      colnames(counts),
      metadata$cell_id
    ),
    ,
    drop = FALSE
  ]

  rownames(metadata) <- as.character(
    metadata$cell_id
  )


  # ------------------------------------------------------------
  # Alignment validation
  # ------------------------------------------------------------

  if (!identical(
    colnames(counts),
    rownames(metadata)
  )) {
    stop(
      "Count matrix and metadata alignment failed.",
      call. = FALSE
    )
  }


  # ============================================================
  # 6. Create Seurat object
  # ============================================================

  if (is.null(project)) {
    project <- as.character(
      region_data$region_id
    )
  }

  object <- Seurat::CreateSeuratObject(
    counts = counts,
    meta.data = metadata,
    assay = assay,
    project = project,
    min.cells = 0L,
    min.features = 0L
  )


  # ============================================================
  # 7. Create spatial centroid object
  # ============================================================

  centroids <- data.frame(
    x = as.numeric(metadata$x_centroid),
    y = as.numeric(metadata$y_centroid),
    row.names = metadata$cell_id,
    check.names = FALSE
  )

  spatial_centroids <-
    SeuratObject::CreateCentroids(
      coords = centroids
    )


  # ============================================================
  # 8. Create initial FOV from centroids
  # ============================================================

  spatial_fov <-
    SeuratObject::CreateFOV(
      coords = spatial_centroids,
      type = "centroids",
      molecules = NULL,
      assay = assay,
      key = paste0(fov, "_")
    )

  object[[fov]] <- spatial_fov


  # ============================================================
  # 9. Optional native Xenium segmentation
  # ============================================================

  segmentation_loaded <- character()


  # ------------------------------------------------------------
  # Boundary reader helper
  # ------------------------------------------------------------

  read_xenium_boundary <- function(
      path,
      selected_cells
  ) {

    if (!file.exists(path)) {
      return(NULL)
    }

    boundary <- arrow::read_parquet(
      path,
      as_data_frame = TRUE
    )


    # ----------------------------------------------------------
    # Required cell ID
    # ----------------------------------------------------------

    if (!"cell_id" %in% names(boundary)) {
      stop(
        "Boundary table does not contain cell_id: ",
        path,
        call. = FALSE
      )
    }


    # ----------------------------------------------------------
    # Resolve x/y coordinate columns
    # ----------------------------------------------------------

    x_candidates <- c(
      "vertex_x",
      "x",
      "x_location",
      "x_centroid"
    )

    y_candidates <- c(
      "vertex_y",
      "y",
      "y_location",
      "y_centroid"
    )

    x_col <- intersect(
      x_candidates,
      names(boundary)
    )

    y_col <- intersect(
      y_candidates,
      names(boundary)
    )

    if (
      !length(x_col) ||
      !length(y_col)
    ) {
      stop(
        "Cannot identify polygon x/y columns in ",
        basename(path),
        ". Columns found: ",
        paste(
          names(boundary),
          collapse = ", "
        ),
        call. = FALSE
      )
    }

    x_col <- x_col[[1]]
    y_col <- y_col[[1]]


    # ----------------------------------------------------------
    # Restrict polygons to cells retained in object
    # ----------------------------------------------------------

    boundary <- boundary[
      as.character(boundary$cell_id) %in%
        selected_cells,
      ,
      drop = FALSE
    ]

    if (!nrow(boundary)) {
      stop(
        "No selected cells were found in ",
        basename(path),
        call. = FALSE
      )
    }


    # ----------------------------------------------------------
    # Convert to Seurat segmentation coordinate format
    # ----------------------------------------------------------

    polygon <- data.frame(
      x = as.numeric(
        boundary[[x_col]]
      ),
      y = as.numeric(
        boundary[[y_col]]
      ),
      cell = as.character(
        boundary$cell_id
      ),
      stringsAsFactors = FALSE
    )


    # ----------------------------------------------------------
    # Remove malformed vertices
    # ----------------------------------------------------------

    polygon <- polygon[
      is.finite(polygon$x) &
        is.finite(polygon$y) &
        !is.na(polygon$cell) &
        nzchar(polygon$cell),
      ,
      drop = FALSE
    ]

    if (!nrow(polygon)) {
      stop(
        "No valid polygon vertices remain in ",
        basename(path),
        call. = FALSE
      )
    }


    # ----------------------------------------------------------
    # Keep vertices grouped by cell
    # ----------------------------------------------------------

    polygon <- polygon[
      order(polygon$cell),
      ,
      drop = FALSE
    ]

    rownames(polygon) <- NULL


    # ----------------------------------------------------------
    # Ensure only selected cells remain
    # ----------------------------------------------------------

    unexpected_cells <- setdiff(
      unique(polygon$cell),
      selected_cells
    )

    if (length(unexpected_cells)) {
      stop(
        "Boundary reader retained unexpected cells.",
        call. = FALSE
      )
    }


    polygon
  }


  # ============================================================
  # 9A. Import native geometry
  # ============================================================

  if (!is.null(xenium_dir)) {

    if (!dir.exists(xenium_dir)) {
      stop(
        "xenium_dir does not exist: ",
        xenium_dir,
        call. = FALSE
      )
    }

    if (
      !requireNamespace(
        "arrow",
        quietly = TRUE
      )
    ) {
      stop(
        "Package 'arrow' is required to import Xenium polygons.",
        call. = FALSE
      )
    }


    # ==========================================================
    # 9A-1. Cell segmentation
    # ==========================================================

    if (isTRUE(
      include_cell_segmentation
    )) {

      cell_boundary_path <- file.path(
        xenium_dir,
        "cell_boundaries.parquet"
      )


      if (file.exists(
        cell_boundary_path
      )) {

        cell_polygon <-
          read_xenium_boundary(
            path = cell_boundary_path,
            selected_cells = cell_ids
          )


        message(
          "Creating compact cell segmentation from ",
          format(
            nrow(cell_polygon),
            big.mark = ","
          ),
          " vertices across ",
          format(
            length(
              unique(
                cell_polygon$cell
              )
            ),
            big.mark = ","
          ),
          " cells."
        )


        cell_segmentation <-
          SeuratObject::CreateSegmentation(
            coords = cell_polygon,
            compact = TRUE
          )


        object[[fov]][[
          "segmentation"
        ]] <- cell_segmentation


        segmentation_loaded <- c(
          segmentation_loaded,
          "segmentation"
        )


        rm(
          cell_polygon,
          cell_segmentation
        )

        invisible(gc())
      } else {

        warning(
          "Cell segmentation requested but ",
          "cell_boundaries.parquet was not found.",
          call. = FALSE
        )
      }
    }


    # ==========================================================
    # 9A-2. Nucleus segmentation
    # ==========================================================

    if (isTRUE(
      include_nucleus_segmentation
    )) {

      nucleus_boundary_path <- file.path(
        xenium_dir,
        "nucleus_boundaries.parquet"
      )


      if (file.exists(
        nucleus_boundary_path
      )) {

        nucleus_polygon <-
          read_xenium_boundary(
            path = nucleus_boundary_path,
            selected_cells = cell_ids
          )


        message(
          "Creating compact nucleus segmentation from ",
          format(
            nrow(nucleus_polygon),
            big.mark = ","
          ),
          " vertices across ",
          format(
            length(
              unique(
                nucleus_polygon$cell
              )
            ),
            big.mark = ","
          ),
          " cells."
        )


        nucleus_segmentation <-
          SeuratObject::CreateSegmentation(
            coords = nucleus_polygon,
            compact = TRUE
          )


        object[[fov]][[
          "nucleus_segmentation"
        ]] <- nucleus_segmentation


        segmentation_loaded <- c(
          segmentation_loaded,
          "nucleus_segmentation"
        )


        rm(
          nucleus_polygon,
          nucleus_segmentation
        )

        invisible(gc())
      } else {

        warning(
          "Nucleus segmentation requested but ",
          "nucleus_boundaries.parquet was not found.",
          call. = FALSE
        )
      }
    }

  }


  # ============================================================
  # 10. Set default FOV
  # ============================================================

  SeuratObject::DefaultFOV(
    object
  ) <- fov


  # ============================================================
  # 11. Determine available boundaries
  # ============================================================

  available_boundaries <-
    SeuratObject::Boundaries(
      object[[fov]]
    )


  # ------------------------------------------------------------
  # Keep centroids as default for whole-section visualization
  # ------------------------------------------------------------

  if (
    "centroids" %in%
      available_boundaries
  ) {

    SeuratObject::DefaultBoundary(
      object[[fov]]
    ) <- "centroids"
  }


  # ============================================================
  # 12. Provenance
  # ============================================================

  object@misc$downstream_bundle <- list(

    region_id =
      as.character(
        region_data$region_id
      ),

    section_status =
      region_data$section_status,

    cell_mask =
      if (is.null(mask)) {
        "ALL_CELLS"
      } else {
        mask
      },

    gene_count =
      length(genes),

    genes =
      genes,

    raw_counts_preserved =
      isTRUE(
        region_data$raw_counts_preserved
      ),

    spatial_source =
      if (is.null(xenium_dir)) {

        "DOWNSTREAM_CELL_METADATA_CENTROIDS"

      } else {

        "DOWNSTREAM_COUNTS_PLUS_NATIVE_XENIUM_GEOMETRY"

      },

    xenium_dir =
      if (is.null(xenium_dir)) {

        NA_character_

      } else {

        normalizePath(
          xenium_dir,
          winslash = "/",
          mustWork = TRUE
        )

      },

    spatial_boundaries =
      available_boundaries,

    segmentation_loaded =
      segmentation_loaded,

    expression_source =
      "QC_APPROVED_DOWNSTREAM_INPUT_RDS",

    segmentation_source =
      if (length(
        segmentation_loaded
      )) {
        "NATIVE_XENIUM_BOUNDARY_PARQUET"
      } else {
        "NOT_LOADED"
      }
  )


  # ============================================================
  # 13. Final integrity checks
  # ============================================================

  # ------------------------------------------------------------
  # Expression / metadata alignment
  # ------------------------------------------------------------

  if (!identical(
    colnames(object),
    rownames(object@meta.data)
  )) {

    stop(
      "Final Seurat cell/metadata alignment failed.",
      call. = FALSE
    )
  }


  # ------------------------------------------------------------
  # FOV cell coverage
  # ------------------------------------------------------------

  spatial_cells <- Cells(
    object[[fov]]
  )

  missing_spatial <- setdiff(
    colnames(object),
    spatial_cells
  )

  if (length(missing_spatial)) {

    stop(
      length(missing_spatial),
      " Seurat cells are missing from the FOV.",
      call. = FALSE
    )
  }


  # ------------------------------------------------------------
  # Verify segmentation contains only object cells
  # ------------------------------------------------------------

  if ("segmentation" %in% available_boundaries) {

    segmentation_cells <- Cells(
      object[[fov]][[
        "segmentation"
      ]]
    )

    unexpected <- setdiff(
      segmentation_cells,
      colnames(object)
    )

    if (length(unexpected)) {
      stop(
        "Cell segmentation contains cells absent from the Seurat object.",
        call. = FALSE
      )
    }
  }


  if (
    "nucleus_segmentation" %in%
      available_boundaries
  ) {

    nucleus_cells <- Cells(
      object[[fov]][[
        "nucleus_segmentation"
      ]]
    )

    unexpected <- setdiff(
      nucleus_cells,
      colnames(object)
    )

    if (length(unexpected)) {
      stop(
        "Nucleus segmentation contains cells absent from the Seurat object.",
        call. = FALSE
      )
    }
  }


  # ============================================================
  # 14. Report result
  # ============================================================

  message(
    "Created spatial Seurat object: ",
    project
  )

  message(
    "  Cells: ",
    format(
      ncol(object),
      big.mark = ","
    )
  )

  message(
    "  Genes: ",
    format(
      nrow(object),
      big.mark = ","
    )
  )

  message(
    "  FOV: ",
    fov
  )

  message(
    "  Boundaries: ",
    paste(
      available_boundaries,
      collapse = ", "
    )
  )

  if (length(
    segmentation_loaded
  )) {

    message(
      "  Native Xenium segmentations loaded: ",
      paste(
        segmentation_loaded,
        collapse = ", "
      )
    )
  }


  # ============================================================
  # 15. Return
  # ============================================================

  object
}

#' Calculate marker scores.
#'
#' @param object Required `object` input; validated before computation.
#' @param marker_df Required `marker_df` input; validated before computation.
#' @param group_col Required `group_col` input; validated before computation.
#' @param assay Optional `assay` input with the default shown in the function signature.
#' @param layer Optional `layer` input with the default shown in the function signature.
#' @param prefix Optional `prefix` input with the default shown in the function signature.
#' @param z_cap Optional `z_cap` input with the default shown in the function signature.
#' @return Computed evidence or annotations aligned to the supplied rows/cells; raw counts are unchanged.
calculate_marker_scores <- function(
  object,
  marker_df,
  group_col,
  assay = "Xenium",
  layer = "data",
  prefix = "Score",
  z_cap = 3
) {

  stopifnot(
    all(c("Gene_Symbol", group_col) %in% colnames(marker_df))
  )

  marker_df2 <- marker_df %>%
    filter(
      Gene_Symbol %in% rownames(object),
      !is.na(.data[[group_col]])
    ) %>%
    distinct(
      Gene_Symbol,
      .data[[group_col]]
    )

  genes <- unique(marker_df2$Gene_Symbol)

  # genes x cells
  mat <- LayerData(
    object = object,
    assay = assay,
    layer = layer
  )[genes, , drop = FALSE]

  # gene-wise mean / SD
  gene_mean <- Matrix::rowMeans(mat)

  gene_sq_mean <- Matrix::rowMeans(mat^2)

  gene_sd <- sqrt(
    pmax(
      gene_sq_mean - gene_mean^2,
      0
    )
  )

  # Drop genes with no variation
  keep <- is.finite(gene_sd) & gene_sd > 0

  mat <- mat[keep, , drop = FALSE]
  gene_mean <- gene_mean[keep]
  gene_sd <- gene_sd[keep]

  marker_df2 <- marker_df2 %>%
    filter(Gene_Symbol %in% rownames(mat))

  # Dense only for the small marker matrix, not the whole Xenium matrix
  z <- as.matrix(mat)

  z <- sweep(
    z,
    1,
    gene_mean,
    "-"
  )

  z <- sweep(
    z,
    1,
    gene_sd,
    "/"
  )

  # Prevent one extreme gene from dominating
  z[z >  z_cap] <-  z_cap
  z[z < -z_cap] <- -z_cap

  groups <- unique(marker_df2[[group_col]])

  score_df <- matrix(
    NA_real_,
    nrow = ncol(z),
    ncol = length(groups),
    dimnames = list(
      colnames(z),
      paste0(prefix, "_", make.names(groups))
    )
  )

  detection_df <- score_df

  for (grp in groups) {

    genes_grp <- marker_df2 %>%
      filter(.data[[group_col]] == grp) %>%
      pull(Gene_Symbol) %>%
      unique()

    genes_grp <- intersect(
      genes_grp,
      rownames(z)
    )

    # Number of groups each gene contributes to
    n_memberships <- marker_df2 %>%
      filter(Gene_Symbol %in% genes_grp) %>%
      count(Gene_Symbol, name = "n_group")

    weights <- 1 / n_memberships$n_group
    names(weights) <- n_memberships$Gene_Symbol

    weights <- weights[genes_grp]
    weights <- weights / sum(weights)

    score_df[
      ,
      paste0(prefix, "_", make.names(grp))
    ] <- as.numeric(
      crossprod(
        weights,
        z[genes_grp, , drop = FALSE]
      )
    )

    # Fraction of marker genes detected in each cell
    raw_grp <- mat[
      genes_grp,
      ,
      drop = FALSE
    ]

    detection_df[
      ,
      paste0(prefix, "_", make.names(grp))
    ] <- Matrix::colMeans(raw_grp > 0)
  }

  list(
    score = as.data.frame(score_df),
    detection = as.data.frame(detection_df)
  )
}

#' Score marker groups by cluster.
#'
#' @param object Required `object` input; validated before computation.
#' @param marker_df Required `marker_df` input; validated before computation.
#' @param group_col Required `group_col` input; validated before computation.
#' @param cluster_col Optional `cluster_col` input with the default shown in the function signature.
#' @param assay Optional `assay` input with the default shown in the function signature.
#' @param layer Optional `layer` input with the default shown in the function signature.
#' @param score_prefix Optional `score_prefix` input with the default shown in the function signature.
#' @param detect_prefix Optional `detect_prefix` input with the default shown in the function signature.
#' @param label_prefix Optional `label_prefix` input with the default shown in the function signature.
#' @param z_cap Optional `z_cap` input with the default shown in the function signature.
#' @return Computed evidence or annotations aligned to the supplied rows/cells; raw counts are unchanged.
score_marker_groups_by_cluster <- function(
  object,
  marker_df,
  group_col,
  cluster_col = "cluster_res_1_2",
  assay = "Xenium",
  layer = "data",
  score_prefix = "MainScore",
  detect_prefix = "MainDetect",
  label_prefix = "Main",
  z_cap = 3
) {

  stopifnot(
    group_col %in% colnames(marker_df),
    cluster_col %in% colnames(object@meta.data)
  )

  # ------------------------------------------------------------
  # 1. Calculate cell-level marker scores
  # ------------------------------------------------------------

  scores <- calculate_marker_scores(
    object = object,
    marker_df = marker_df,
    group_col = group_col,
    assay = assay,
    layer = layer,
    prefix = score_prefix,
    z_cap = z_cap
  )

  # ------------------------------------------------------------
  # 2. Add score metadata
  # ------------------------------------------------------------

  object <- Seurat::AddMetaData(
    object = object,
    metadata = scores$score
  )

  # ------------------------------------------------------------
  # 3. Add marker-detection metadata
  # ------------------------------------------------------------

  detect_df <- scores$detection

  colnames(detect_df) <- sub(
    paste0("^", score_prefix, "_"),
    paste0(detect_prefix, "_"),
    colnames(detect_df)
  )

  object <- Seurat::AddMetaData(
    object = object,
    metadata = detect_df
  )

  # ------------------------------------------------------------
  # 4. Find score columns
  # ------------------------------------------------------------

  score_cols <- grep(
    paste0("^", score_prefix, "_"),
    colnames(object@meta.data),
    value = TRUE
  )

  if (length(score_cols) < 2) {
    stop(
      "Fewer than two score groups were found for ",
      group_col,
      "."
    )
  }

  # ------------------------------------------------------------
  # 5. Summarise scores by cluster
  # ------------------------------------------------------------

  cluster_scores <- object@meta.data %>%
    dplyr::group_by(
      .data[[cluster_col]]
    ) %>%
    dplyr::summarise(
      n_cells = dplyr::n(),

      dplyr::across(
        dplyr::all_of(score_cols),
        list(
          mean = ~ mean(.x, na.rm = TRUE),
          median = ~ median(.x, na.rm = TRUE)
        )
      ),

      .groups = "drop"
    )

  # ------------------------------------------------------------
  # 6. Extract mean scores for ranking
  # ------------------------------------------------------------

  mean_cols <- grep(
    paste0(
      "^",
      score_prefix,
      "_.*_mean$"
    ),
    colnames(cluster_scores),
    value = TRUE
  )

  if (length(mean_cols) < 2) {
    stop(
      "Could not identify at least two mean score columns."
    )
  }

  score_mat <- as.matrix(
    cluster_scores[
      ,
      mean_cols,
      drop = FALSE
    ]
  )

  rownames(score_mat) <- as.character(
    cluster_scores[[cluster_col]]
  )

  # ------------------------------------------------------------
  # 7. Recover biological group names
  # ------------------------------------------------------------

  labels <- sub(
    "_mean$",
    "",
    sub(
      paste0("^", score_prefix, "_"),
      "",
      mean_cols
    )
  )

  # ------------------------------------------------------------
  # 8. Identify best and second-best groups
  # ------------------------------------------------------------

  annotation <- lapply(
    seq_len(nrow(score_mat)),
    function(i) {

      x <- score_mat[i, ]

      # Deal safely with NA/NaN
      x[
        !is.finite(x)
      ] <- NA_real_

      valid <- which(
        !is.na(x)
      )

      if (length(valid) == 0) {

        return(
          tibble::tibble(
            cluster = rownames(score_mat)[i],
            Best_label = NA_character_,
            Best_score = NA_real_,
            Second_label = NA_character_,
            Second_score = NA_real_,
            Score_margin = NA_real_
          )
        )
      }

      ord <- valid[
        order(
          x[valid],
          decreasing = TRUE
        )
      ]

      best_idx <- ord[1]

      if (length(ord) >= 2) {

        second_idx <- ord[2]

        second_label <- labels[second_idx]
        second_score <- x[second_idx]

        margin <-
          x[best_idx] -
          x[second_idx]

      } else {

        second_label <- NA_character_
        second_score <- NA_real_
        margin <- NA_real_
      }

      tibble::tibble(
        cluster =
          rownames(score_mat)[i],

        Best_label =
          labels[best_idx],

        Best_score =
          x[best_idx],

        Second_label =
          second_label,

        Second_score =
          second_score,

        Score_margin =
          margin
      )
    }
  ) %>%
    dplyr::bind_rows()

  # ------------------------------------------------------------
  # 9. Rename output columns appropriately
  # ------------------------------------------------------------

  colnames(annotation)[
    colnames(annotation) == "cluster"
  ] <- cluster_col

  colnames(annotation)[
    colnames(annotation) == "Best_label"
  ] <- paste0(
    label_prefix,
    "_label"
  )

  colnames(annotation)[
    colnames(annotation) == "Best_score"
  ] <- paste0(
    label_prefix,
    "_score"
  )

  colnames(annotation)[
    colnames(annotation) == "Second_label"
  ] <- paste0(
    label_prefix,
    "_second"
  )

  colnames(annotation)[
    colnames(annotation) == "Second_score"
  ] <- paste0(
    label_prefix,
    "_second_score"
  )

  colnames(annotation)[
    colnames(annotation) == "Score_margin"
  ] <- paste0(
    label_prefix,
    "_margin"
  )

  # ------------------------------------------------------------
  # 10. Sort clusters numerically where possible
  # ------------------------------------------------------------

  annotation <- annotation %>%
    dplyr::mutate(
      .cluster_numeric =
        suppressWarnings(
          as.numeric(
            as.character(
              .data[[cluster_col]]
            )
          )
        )
    ) %>%
    dplyr::arrange(
      .cluster_numeric,
      .data[[cluster_col]]
    ) %>%
    dplyr::select(
      -.cluster_numeric
    )

  # ------------------------------------------------------------
  # 11. Return everything
  # ------------------------------------------------------------

  return(
    list(
      object = object,

      cell_scores =
        scores$score,

      cell_detection =
        detect_df,

      cluster_scores =
        cluster_scores,

      score_matrix =
        score_mat,

      annotation =
        annotation
    )
  )
}

#' Refine xenium celltypes.
#'
#' @param object Required `object` input; validated before computation.
#' @param reduction Optional `reduction` input with the default shown in the function signature.
#' @param dims Optional `dims` input with the default shown in the function signature.
#' @param k Optional `k` input with the default shown in the function signature.
#' @param self_weight Optional `self_weight` input with the default shown in the function signature.
#' @param wang_prefix Optional `wang_prefix` input with the default shown in the function signature.
#' @param cluster_main_col Optional `cluster_main_col` input with the default shown in the function signature.
#' @param cluster_subtype_col Optional `cluster_subtype_col` input with the default shown in the function signature.
#' @param ontology_gap_max_immune_prob Optional `ontology_gap_max_immune_prob` input with the default shown in the function signature.
#' @param wang_main_review_score Optional `wang_main_review_score` input with the default shown in the function signature.
#' @param wang_main_review_margin Optional `wang_main_review_margin` input with the default shown in the function signature.
#' @param wang_subtype_review_margin Optional `wang_subtype_review_margin` input with the default shown in the function signature.
#' @param xenium_subtype_review_margin Optional `xenium_subtype_review_margin` input with the default shown in the function signature.
#' @param verbose Optional `verbose` input with the default shown in the function signature.
#' @return Computed evidence or annotations aligned to the supplied rows/cells; raw counts are unchanged.
refine_xenium_celltypes <- function(
  object,
  reduction = "pca",
  dims = 1:12,
  k = 15,
  self_weight = 0.70,
  wang_prefix = "RefAll",
  cluster_main_col = "Xenium_cluster_main",
  cluster_subtype_col = "Xenium_cluster_subtype",

  # Wang does not contain Schwann / lymphatic EC.
  # Only rescue these when the Wang neighborhood is not predominantly immune.
  ontology_gap_max_immune_prob = 0.50,

  # Optional empirical review thresholds.
  # Leave NULL initially; inspect distributions before choosing thresholds.
  wang_main_review_score = NULL,
  wang_main_review_margin = NULL,
  wang_subtype_review_margin = NULL,
  xenium_subtype_review_margin = NULL,

  verbose = TRUE
) {

  # ============================================================
  # 0. REQUIREMENTS
  # ============================================================

  if (!requireNamespace("FNN", quietly = TRUE)) {
    stop(
      "Package 'FNN' is required. Install with install.packages('FNN')."
    )
  }

  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package 'Matrix' is required.")
  }

  stopifnot(
    inherits(object, "Seurat"),
    reduction %in% Reductions(object),
    cluster_main_col %in% colnames(object@meta.data),
    cluster_subtype_col %in% colnames(object@meta.data),
    self_weight >= 0,
    self_weight <= 1,
    k >= 2
  )

  cells <- colnames(object)

  # Ensure metadata and embedding are aligned to exactly the same cell order.
  md <- object@meta.data[
    cells,
    ,
    drop = FALSE
  ]

  embedding <- Embeddings(
    object,
    reduction = reduction
  )

  if (!all(cells %in% rownames(embedding))) {
    stop("Not all Seurat cells are present in the requested reduction.")
  }

  if (max(dims) > ncol(embedding)) {
    stop(
      "Requested dims exceed the dimensions available in reduction '",
      reduction,
      "'."
    )
  }

  embedding <- embedding[
    cells,
    dims,
    drop = FALSE
  ]


  # ============================================================
  # 1. INTERNAL HELPER FUNCTIONS
  # ============================================================

  # ------------------------------------------------------------
  # Extract a matrix of metadata columns sharing a prefix.
  #
  # Example:
  # RefAll_main_prediction.score.Macrophage
  # becomes column name:
  # Macrophage
  # ------------------------------------------------------------

  extract_metadata_matrix <- function(
    metadata,
    prefix,
    exclude = character()
  ) {

    cols <- grep(
      paste0("^", prefix),
      colnames(metadata),
      value = TRUE
    )

    cols <- setdiff(
      cols,
      exclude
    )

    if (!length(cols)) {
      stop(
        "No metadata columns found with prefix: ",
        prefix
      )
    }

    x <- as.matrix(
      metadata[
        ,
        cols,
        drop = FALSE
      ]
    )

    storage.mode(x) <- "double"

    colnames(x) <- sub(
      paste0("^", prefix),
      "",
      colnames(x)
    )

    x
  }


  # ------------------------------------------------------------
  # Calculate top / second-best class and margin.
  #
  # Works for:
  #   Wang probabilities
  #   Xenium marker scores
  # ------------------------------------------------------------

  get_top_prediction <- function(x) {

    x2 <- x

    x2[
      !is.finite(x2)
    ] <- -Inf

    n <- nrow(x2)

    best_idx <- max.col(
      x2,
      ties.method = "first"
    )

    best_score <- x2[
      cbind(
        seq_len(n),
        best_idx
      )
    ]

    best_label <- colnames(x2)[
      best_idx
    ]

    if (ncol(x2) >= 2) {

      x_second <- x2

      x_second[
        cbind(
          seq_len(n),
          best_idx
        )
      ] <- -Inf

      second_idx <- max.col(
        x_second,
        ties.method = "first"
      )

      second_score <- x_second[
        cbind(
          seq_len(n),
          second_idx
        )
      ]

      second_label <- colnames(x2)[
        second_idx
      ]

      margin <- best_score -
        second_score

    } else {

      second_score <- rep(
        NA_real_,
        n
      )

      second_label <- rep(
        NA_character_,
        n
      )

      margin <- rep(
        NA_real_,
        n
      )
    }

    list(
      label = best_label,
      score = best_score,
      second = second_label,
      second_score = second_score,
      margin = margin
    )
  }


  # ------------------------------------------------------------
  # Restrict subtype classification to the selected broad lineage.
  #
  # Example:
  # if broad lineage = Macrophage,
  # only compare:
  #
  #   LYVE1_resident_Mac
  #   TREM2_LAM
  #   Inflammatory_Mac
  #
  # rather than allowing cDC2 / B / T etc. to compete.
  #
  # This hierarchical restriction is scientifically important.
  # ------------------------------------------------------------

  restricted_top_prediction <- function(
    score_matrix,
    broad_labels,
    subtype_to_main
  ) {

    n <- nrow(score_matrix)

    out <- list(
      label = rep(NA_character_, n),
      score = rep(NA_real_, n),
      second = rep(NA_character_, n),
      second_score = rep(NA_real_, n),
      margin = rep(NA_real_, n)
    )

    for (main_type in unique(
      broad_labels[
        !is.na(broad_labels)
      ]
    )) {

      candidate_subtypes <- names(
        subtype_to_main
      )[
        subtype_to_main == main_type
      ]

      candidate_subtypes <- intersect(
        candidate_subtypes,
        colnames(score_matrix)
      )

      idx <- which(
        broad_labels == main_type
      )

      if (
        !length(idx) ||
        !length(candidate_subtypes)
      ) {
        next
      }

      tmp <- get_top_prediction(
        score_matrix[
          idx,
          candidate_subtypes,
          drop = FALSE
        ]
      )

      out$label[idx] <- tmp$label
      out$score[idx] <- tmp$score
      out$second[idx] <- tmp$second
      out$second_score[idx] <- tmp$second_score
      out$margin[idx] <- tmp$margin
    }

    out
  }


  # ------------------------------------------------------------
  # Human-readable label formatting.
  # ------------------------------------------------------------

  display_main <- function(x) {

    dplyr::recode(
      x,
      "T.cell" = "T cell",
      "B.cell" = "B cell",
      "Mast.cell" = "Mast cell",
      .default = x
    )
  }


  display_subtype <- function(x) {

    dplyr::recode(
      x,
      "γδ.T" = "γδ T",
      "T.cell" = "T cell",
      "Mast.cell" = "Mast cell",
      "MC" = "Mast cell",
      .default = x
    )
  }


  # ============================================================
  # 2. BUILD TRANSCRIPTOMIC kNN
  # ============================================================
  #
  # IMPORTANT:
  # These are transcriptomic PCA neighbors, NOT spatial neighbors.
  #
  # Spatial proximity tells us which cells live next to one another.
  # Transcriptomic proximity is more appropriate for cell identity.
  # ============================================================

  if (verbose) {
    message(
      "Building transcriptomic kNN: k = ",
      k,
      ", reduction = ",
      reduction,
      ", dims = ",
      paste(range(dims), collapse = ":")
    )
  }

  knn <- FNN::get.knn(
    embedding,
    k = k
  )


  # ============================================================
  # 3. DISTANCE-WEIGHT THE kNN
  # ============================================================
  #
  # Close neighbors contribute more than distant neighbors.
  #
  # A Gaussian kernel is calculated independently for each cell,
  # with the median neighbor distance used as the local bandwidth.
  # ============================================================

  d <- knn$nn.dist

  sigma <- apply(
    d,
    1,
    median
  )

  sigma[
    !is.finite(sigma) |
      sigma <= 0
  ] <- 1

  weights <- exp(
    -(d^2) /
      (2 * sigma^2)
  )

  weights <- weights /
    rowSums(weights)


  # ------------------------------------------------------------
  # Construct sparse cell x cell kNN weight matrix.
  #
  # This is much faster than looping over ~70,000 cells.
  # ------------------------------------------------------------

  n_cells <- nrow(embedding)

  W <- Matrix::sparseMatrix(
    i = rep(
      seq_len(n_cells),
      each = k
    ),
    j = as.vector(
      t(knn$nn.index)
    ),
    x = as.vector(
      t(weights)
    ),
    dims = c(
      n_cells,
      n_cells
    )
  )


  # ------------------------------------------------------------
  # Smooth continuous evidence:
  #
  # final evidence =
  #     self_weight * own evidence
  #   + (1-self_weight) * neighborhood evidence
  #
  # This is deliberately NOT pure kNN majority voting.
  #
  # Retaining substantial self-weight protects rare populations
  # from being erased by abundant neighboring cell types.
  # ------------------------------------------------------------

  smooth_knn_matrix <- function(x) {

    neighbour_signal <- as.matrix(
      W %*% x
    )

    out <-
      self_weight * x +
      (1 - self_weight) *
      neighbour_signal

    rownames(out) <- cells

    out
  }


  # ============================================================
  # 4. WANG RefAll — BROAD CELL-TYPE PROBABILITIES
  # ============================================================

  wang_main_prefix <- paste0(
    wang_prefix,
    "_main_prediction.score."
  )

  wang_main_max_col <- paste0(
    wang_prefix,
    "_main_prediction.score.max"
  )

  wang_main_prob <- extract_metadata_matrix(
    metadata = md,
    prefix = wang_main_prefix,
    exclude = wang_main_max_col
  )

  rownames(wang_main_prob) <- cells

  wang_main_knn_prob <- smooth_knn_matrix(
    wang_main_prob
  )

  wang_main_top <- get_top_prediction(
    wang_main_knn_prob
  )


  # ============================================================
  # 5. WANG IMMUNE PROBABILITY
  # ============================================================
  #
  # Instead of:
  #
  #   Ptprc > 0 = immune
  #
  # sum all Wang probability mass assigned to immune lineages.
  #
  # Ptprc / Itgam remain independent validation variables.
  # ============================================================

  immune_main_raw <- c(
    "Macrophage",
    "Monocyte",
    "Neutrophil",
    "Mast.cell",
    "Eosinophil",
    "T.cell",
    "NK",
    "DC",
    "B.cell",
    "Plasma",
    "ILC"
  )

  immune_prob_cols <- intersect(
    immune_main_raw,
    colnames(wang_main_knn_prob)
  )

  wang_immune_probability <- rowSums(
    wang_main_knn_prob[
      ,
      immune_prob_cols,
      drop = FALSE
    ]
  )


  # ============================================================
  # 6. WANG RefAll — SUBTYPE PROBABILITIES
  # ============================================================

  wang_subtype_prefix <- paste0(
    wang_prefix,
    "_subtype_prediction.score."
  )

  wang_subtype_max_col <- paste0(
    wang_prefix,
    "_subtype_prediction.score.max"
  )

  wang_subtype_prob <- extract_metadata_matrix(
    metadata = md,
    prefix = wang_subtype_prefix,
    exclude = wang_subtype_max_col
  )

  rownames(wang_subtype_prob) <- cells

  wang_subtype_knn_prob <- smooth_knn_matrix(
    wang_subtype_prob
  )

  wang_subtype_top <- get_top_prediction(
    wang_subtype_knn_prob
  )


  # ============================================================
  # 7. XENIUM BROAD MARKER-SCORE VECTORS
  # ============================================================

  xenium_main_score <- extract_metadata_matrix(
    metadata = md,
    prefix = "MainScore_"
  )

  rownames(xenium_main_score) <- cells

  xenium_main_knn_score <- smooth_knn_matrix(
    xenium_main_score
  )

  xenium_main_top <- get_top_prediction(
    xenium_main_knn_score
  )


  # ============================================================
  # 8. XENIUM SUBTYPE MARKER-SCORE VECTORS
  # ============================================================

  xenium_subtype_score <- extract_metadata_matrix(
    metadata = md,
    prefix = "SubtypeScore_"
  )

  rownames(xenium_subtype_score) <- cells

  xenium_subtype_knn_score <- smooth_knn_matrix(
    xenium_subtype_score
  )

  xenium_subtype_top <- get_top_prediction(
    xenium_subtype_knn_score
  )


  # ============================================================
  # 9. DEFINE HIERARCHIES
  # ============================================================
  #
  # Wang subtype -> Wang broad lineage
  #
  # These machine-safe names correspond to the columns generated
  # by TransferData().
  # ============================================================

  wang_subtype_to_main <- c(

    # Stromal
    "ASC" = "Stromal_Fibroblast",
    "APC" = "Stromal_Fibroblast",

    # Vascular / structural
    "Mural" = "Mural",
    "Endothelial" = "Endothelial",
    "Adipocyte" = "Adipocyte",
    "Mesothelial" = "Mesothelial",
    "Epithelial" = "Epithelial",

    # Lymphoid
    "T" = "T.cell",
    "γδ.T" = "T.cell",
    "NK" = "NK",
    "ILC2" = "ILC",
    "B" = "B.cell",
    "Plasma" = "Plasma",

    # Monocyte
    "CCR2_inflammatory_Monocyte" = "Monocyte",
    "CX3CR1_Monocyte" = "Monocyte",

    # Macrophage
    "Inflammatory_Mac" = "Macrophage",
    "LYVE1_resident_Mac" = "Macrophage",
    "TREM2_LAM" = "Macrophage",

    # DC
    "cDC1" = "DC",
    "cDC2" = "DC",
    "CCR7_migratory_DC" = "DC",

    # Other immune
    "Neutrophil" = "Neutrophil",
    "Eosinophil" = "Eosinophil",
    "Mast.cell" = "Mast.cell"
  )


  # ------------------------------------------------------------
  # Xenium subtype -> Xenium broad lineage
  # ------------------------------------------------------------

  xenium_subtype_to_main <- c(

    # EC
    "Capillary_EC" = "Endothelial",
    "Arterial_EC" = "Endothelial",
    "Venous_EC" = "Endothelial",
    "Lymphatic_EC" = "Endothelial",

    # Adipocyte
    "Adipocyte" = "Adipocyte",

    # Stromal
    "ASC" = "Stromal_Fibroblast",
    "Fibroblast" = "Stromal_Fibroblast",

    # Mural
    "VSMC" = "Mural",
    "Pericyte" = "Mural",

    # Neural
    "Schwann" = "Neural",

    # Myeloid / immune
    "LYVE_macrophage" = "Macrophage",
    "Scavenging_macrophage" = "Macrophage",
    "Monocyte" = "Monocyte",
    "Neutrophil" = "Neutrophil",
    "MC" = "Mast.cell",
    "Eosinophil" = "Eosinophil",

    # Lymphoid
    "T.cell" = "T.cell",
    "NK" = "NK",

    # DC
    "DC" = "DC",
    "cDC" = "DC",

    # Mesothelial
    "Mesothelial" = "Mesothelial"
  )


  # ============================================================
  # 10. INITIAL BROAD LINEAGE = WANG RefAll
  # ============================================================
  #
  # Wang provides the cell-level scaffold.
  #
  # This is NOT yet the final subtype.
  # ============================================================

  final_main_raw <- wang_main_top$label

  final_main_score <- wang_main_top$score
  final_main_margin <- wang_main_top$margin

  final_source <- rep(
    "Wang_RefAll_broad",
    n_cells
  )

  final_reason <- rep(
    paste0(
      "Broad lineage assigned from pooled Wang Science 2025 ",
      "reference after conservative transcriptomic-kNN smoothing."
    ),
    n_cells
  )


  # ============================================================
  # 11. ONTOLOGY-GAP RESCUE
  # ============================================================
  #
  # Wang does NOT contain some biologically relevant populations.
  #
  # We should not interpret its nearest available class as evidence
  # against a class that was absent from the reference.
  #
  # Require concordance between:
  #
  #   1. Xenium cluster-level marker annotation
  #   2. cell/kNN Xenium broad score
  #   3. cell/kNN Xenium subtype score
  #
  # This makes the override intentionally conservative.
  # ============================================================

  cluster_subtype <- as.character(
    md[[cluster_subtype_col]]
  )


  # -----------------------------
  # Schwann
  # -----------------------------

  schwann_support <-
    cluster_subtype == "Schwann" &
    xenium_main_top$label == "Neural" &
    xenium_subtype_top$label == "Schwann" &
    wang_immune_probability <
      ontology_gap_max_immune_prob

  final_main_raw[
    schwann_support
  ] <- "Neural"

  final_main_score[
    schwann_support
  ] <- xenium_main_top$score[
    schwann_support
  ]

  final_main_margin[
    schwann_support
  ] <- xenium_main_top$margin[
    schwann_support
  ]

  final_source[
    schwann_support
  ] <- "Xenium_ontology_gap_Schwann"

  final_reason[
    schwann_support
  ] <- paste0(
    "Schwann is absent from the Wang ontology; ",
    "Xenium cluster, broad-score neighborhood and subtype-score ",
    "neighborhood independently support Schwann identity."
  )


  # -----------------------------
  # Lymphatic endothelial cells
  # -----------------------------
  #
  # Wang contains arterial/capillary-derived Endothelial but not
  # an explicit lymphatic EC population. Lymphatic EC may therefore
  # map to Endothelial, Mesothelial or stromal-like reference cells.
  # -----------------------------

  lymphatic_support <-
    cluster_subtype == "Lymphatic_EC" &
    xenium_main_top$label == "Endothelial" &
    xenium_subtype_top$label == "Lymphatic_EC" &
    wang_immune_probability <
      ontology_gap_max_immune_prob

  final_main_raw[
    lymphatic_support
  ] <- "Endothelial"

  final_main_score[
    lymphatic_support
  ] <- xenium_main_top$score[
    lymphatic_support
  ]

  final_main_margin[
    lymphatic_support
  ] <- xenium_main_top$margin[
    lymphatic_support
  ]

  final_source[
    lymphatic_support
  ] <- "Xenium_ontology_gap_Lymphatic_EC"

  final_reason[
    lymphatic_support
  ] <- paste0(
    "Lymphatic EC is absent as an explicit Wang subtype; ",
    "concordant Xenium cluster and neighborhood marker evidence ",
    "supports endothelial/lymphatic identity."
  )


  # ============================================================
  # 12. WANG SUBTYPE — RESTRICTED TO FINAL BROAD LINEAGE
  # ============================================================
  #
  # Prevent biologically impossible combinations such as:
  #
  #   Main = Macrophage
  #   Subtype = cDC2
  #
  # Only subtypes belonging to the selected broad lineage compete.
  # ============================================================

  wang_subtype_restricted <-
    restricted_top_prediction(
      score_matrix = wang_subtype_knn_prob,
      broad_labels = final_main_raw,
      subtype_to_main = wang_subtype_to_main
    )


  # ============================================================
  # 13. XENIUM SUBTYPE — RESTRICTED TO FINAL BROAD LINEAGE
  # ============================================================

  xenium_subtype_restricted <-
    restricted_top_prediction(
      score_matrix = xenium_subtype_knn_score,
      broad_labels = final_main_raw,
      subtype_to_main = xenium_subtype_to_main
    )


  # ============================================================
  # 14. FINAL SUBTYPE ARBITRATION
  # ============================================================

  final_subtype_raw <- rep(
    NA_character_,
    n_cells
  )

  final_subtype_score <- rep(
    NA_real_,
    n_cells
  )

  final_subtype_margin <- rep(
    NA_real_,
    n_cells
  )


  # ------------------------------------------------------------
  # A. IMMUNE
  #
  # Wang is primary because it contains a substantially richer
  # immune taxonomy than the targeted Xenium marker panel.
  # ------------------------------------------------------------

  immune_idx <- final_main_raw %in%
    immune_main_raw

  final_subtype_raw[
    immune_idx
  ] <- wang_subtype_restricted$label[
    immune_idx
  ]

  final_subtype_score[
    immune_idx
  ] <- wang_subtype_restricted$score[
    immune_idx
  ]

  final_subtype_margin[
    immune_idx
  ] <- wang_subtype_restricted$margin[
    immune_idx
  ]

  final_source[
    immune_idx
  ] <- "Wang_RefAll_immune"

  final_reason[
    immune_idx
  ] <- paste0(
    "Immune broad lineage and subtype are reference-driven; ",
    "Xenium marker scores and neighborhood evidence are retained ",
    "as orthogonal validation rather than competing equally."
  )


  # ------------------------------------------------------------
  # B. ENDOTHELIAL
  #
  # Wang establishes broad endothelial identity.
  # Xenium resolves:
  #
  #   Capillary
  #   Arterial (if present)
  #   Venous
  #   Lymphatic
  # ------------------------------------------------------------

  endothelial_idx <-
    final_main_raw == "Endothelial"

  final_subtype_raw[
    endothelial_idx
  ] <- xenium_subtype_restricted$label[
    endothelial_idx
  ]

  final_subtype_score[
    endothelial_idx
  ] <- xenium_subtype_restricted$score[
    endothelial_idx
  ]

  final_subtype_margin[
    endothelial_idx
  ] <- xenium_subtype_restricted$margin[
    endothelial_idx
  ]

  # Do not overwrite explicit lymphatic ontology-gap provenance.
  regular_endothelial <-
    endothelial_idx &
    !lymphatic_support

  final_source[
    regular_endothelial
  ] <- "Wang_broad_Xenium_EC_subtype"

  final_reason[
    regular_endothelial
  ] <- paste0(
    "Wang determines broad endothelial identity; ",
    "Xenium lineage-restricted marker scores resolve EC subtype."
  )


  # ------------------------------------------------------------
  # C. MURAL
  #
  # Wang establishes Mural.
  # Xenium compares only:
  #
  #   Pericyte
  #   VSMC
  #
  # This is particularly important for mixed vascular clusters.
  # ------------------------------------------------------------

  mural_idx <-
    final_main_raw == "Mural"

  final_subtype_raw[
    mural_idx
  ] <- xenium_subtype_restricted$label[
    mural_idx
  ]

  final_subtype_score[
    mural_idx
  ] <- xenium_subtype_restricted$score[
    mural_idx
  ]

  final_subtype_margin[
    mural_idx
  ] <- xenium_subtype_restricted$margin[
    mural_idx
  ]

  final_source[
    mural_idx
  ] <- "Wang_broad_Xenium_mural_subtype"

  final_reason[
    mural_idx
  ] <- paste0(
    "Wang determines broad mural identity; ",
    "Xenium lineage-restricted scores resolve Pericyte versus VSMC."
  )


    # ============================================================
    # D. STROMAL / FIBROBLAST
    #
    # Wang provides ASC vs APC.
    # Xenium independently provides ASC vs Fibroblast.
    #
    # APC is retained as its own biological/reference subtype and
    # is NOT automatically treated as synonymous with Fibroblast.
    # ============================================================
    
    stromal_idx <-
      final_main_raw == "Stromal_Fibroblast"
    
    
    # Wang stromal subtype restricted to ASC/APC
    wang_stromal <- wang_subtype_restricted$label
    
    
    # Xenium stromal subtype restricted to ASC/Fibroblast
    xenium_stromal <- xenium_subtype_restricted$label
    
    
    # Start unresolved
    final_subtype_raw[
      stromal_idx
    ] <- NA_character_
    
    
    # ------------------------------------------------------------
    # 1. ASC concordance
    # ------------------------------------------------------------
    
    stromal_ASC <-
      stromal_idx &
      wang_stromal == "ASC" &
      xenium_stromal == "ASC"
    
    final_subtype_raw[
      stromal_ASC
    ] <- "ASC"
    
    final_subtype_score[
      stromal_ASC
    ] <- wang_subtype_restricted$score[
      stromal_ASC
    ]
    
    final_subtype_margin[
      stromal_ASC
    ] <- wang_subtype_restricted$margin[
      stromal_ASC
    ]
    
    final_source[
      stromal_ASC
    ] <- "Wang_Xenium_stromal_consensus"
    
    final_reason[
      stromal_ASC
    ] <- paste0(
      "Wang and Xenium independently support ASC identity."
    )
    
    
    # ------------------------------------------------------------
    # 2. Wang APC + Xenium Fibroblast-like program
    #
    # This is biologically compatible, but APC is retained as the
    # final subtype because Wang explicitly contains that category.
    # Fibroblast-like Xenium expression is used as supporting evidence.
    # ------------------------------------------------------------
    
    stromal_APC <-
      stromal_idx &
      wang_stromal == "APC" &
      xenium_stromal == "Fibroblast"
    
    final_subtype_raw[
      stromal_APC
    ] <- "APC"
    
    final_subtype_score[
      stromal_APC
    ] <- wang_subtype_restricted$score[
      stromal_APC
    ]
    
    final_subtype_margin[
      stromal_APC
    ] <- wang_subtype_restricted$margin[
      stromal_APC
    ]
    
    final_source[
      stromal_APC
    ] <- "Wang_APC_Xenium_fibroblast_supported"
    
    final_reason[
      stromal_APC
    ] <- paste0(
      "Wang supports APC while Xenium shows a fibroblast-like stromal ",
      "program; APC is retained rather than equated with Fibroblast."
    )
    
    
    # ------------------------------------------------------------
    # 3. Discordant stromal cells
    #
    # Wang APC vs Xenium ASC
    # or
    # Wang ASC vs Xenium Fibroblast
    #
    # Keep Wang subtype provisionally but flag for review.
    # ------------------------------------------------------------
    
    stromal_discordant <-
      stromal_idx &
      !stromal_ASC &
      !stromal_APC
    
    final_subtype_raw[
      stromal_discordant
    ] <- wang_stromal[
      stromal_discordant
    ]
    
    final_subtype_score[
      stromal_discordant
    ] <- wang_subtype_restricted$score[
      stromal_discordant
    ]
    
    final_subtype_margin[
      stromal_discordant
    ] <- wang_subtype_restricted$margin[
      stromal_discordant
    ]
    
    final_source[
      stromal_discordant
    ] <- "Stromal_Wang_Xenium_discordant"
    
    final_reason[
      stromal_discordant
    ] <- paste0(
      "Wang ASC/APC prediction and Xenium ASC/Fibroblast marker ",
      "program disagree; subtype retained provisionally and flagged ",
      "for review."
    )


  # ------------------------------------------------------------
  # E. ADIPOCYTE
  # ------------------------------------------------------------

  adipocyte_idx <-
    final_main_raw == "Adipocyte"

  final_subtype_raw[
    adipocyte_idx
  ] <- "Adipocyte"

  final_source[
    adipocyte_idx
  ] <- "Wang_Xenium_Adipocyte"


  # ------------------------------------------------------------
  # F. NEURAL / SCHWANN
  # ------------------------------------------------------------

  neural_idx <-
    final_main_raw == "Neural"

  final_subtype_raw[
    neural_idx
  ] <- "Schwann"

  final_subtype_score[
    neural_idx
  ] <- xenium_subtype_restricted$score[
    neural_idx
  ]

  final_subtype_margin[
    neural_idx
  ] <- xenium_subtype_restricted$margin[
    neural_idx
  ]


  # ------------------------------------------------------------
  # G. MESOTHELIAL
  # ------------------------------------------------------------

  mesothelial_idx <-
    final_main_raw == "Mesothelial"

  final_subtype_raw[
    mesothelial_idx
  ] <- "Mesothelial"

  final_source[
    mesothelial_idx
  ] <- "Wang_Xenium_Mesothelial"


  # ------------------------------------------------------------
  # H. EPITHELIAL
  #
  # No dedicated Xenium epithelial subtype scoring currently.
  # Keep the Wang class.
  # ------------------------------------------------------------

  epithelial_idx <-
    final_main_raw == "Epithelial"

  final_subtype_raw[
    epithelial_idx
  ] <- "Epithelial"

  final_source[
    epithelial_idx
  ] <- "Wang_RefAll_Epithelial"


  # ============================================================
  # 15. BROAD-LINEAGE AGREEMENT
  # ============================================================

  broad_agreement <-
    final_main_raw ==
    xenium_main_top$label


  # ============================================================
  # 16. REVIEW FLAGS
  # ============================================================
  #
  # Review flags do NOT automatically invalidate a cell.
  # They identify cells deserving closer inspection.
  # ============================================================

  review_flag <- rep(
    FALSE,
    n_cells
  )


  # Missing subtype
  review_flag[
    is.na(final_subtype_raw)
  ] <- TRUE


  # ------------------------------------------------------------
  # Partial ontology-gap evidence.
  #
  # Example:
  # cluster says Schwann but neighborhood does not.
  # ------------------------------------------------------------

  schwann_partial <-
    (
      cluster_subtype == "Schwann" |
      xenium_subtype_top$label == "Schwann"
    ) &
    !schwann_support

  lymphatic_partial <-
    (
      cluster_subtype == "Lymphatic_EC" |
      xenium_subtype_top$label == "Lymphatic_EC"
    ) &
    !lymphatic_support

  review_flag[
    schwann_partial |
      lymphatic_partial
  ] <- TRUE


    review_flag[
      stromal_discordant
    ] <- TRUE
        
  # ------------------------------------------------------------
  # For nonimmune cells, disagreement between Wang broad lineage
  # and Xenium broad marker evidence is worth reviewing.
  #
  # We deliberately do NOT impose this rule on immune cells,
  # because B/T/NK/DC resolution is uneven in the targeted panel.
  # ------------------------------------------------------------

  nonimmune_idx <- !immune_idx

  review_flag[
    nonimmune_idx &
      !broad_agreement
  ] <- TRUE


  # ------------------------------------------------------------
  # Existing segmentation/QC flags
  # ------------------------------------------------------------

  if (
    "segmentation_multiplet_flag" %in%
    colnames(md)
  ) {

    review_flag[
      md$segmentation_multiplet_flag %in%
        TRUE
    ] <- TRUE
  }


  if (
    "do_not_interpret" %in%
    colnames(md)
  ) {

    review_flag[
      md$do_not_interpret %in%
        TRUE
    ] <- TRUE
  }


  # ============================================================
  # 17. OPTIONAL SCORE/MARGIN REVIEW THRESHOLDS
  # ============================================================
  #
  # These are intentionally optional.
  #
  # Do NOT choose thresholds just because they "look reasonable".
  # Inspect score/margin distributions first.
  # ============================================================

  if (!is.null(wang_main_review_score)) {

    review_flag[
      wang_main_top$score <
        wang_main_review_score
    ] <- TRUE
  }


  if (!is.null(wang_main_review_margin)) {

    review_flag[
      wang_main_top$margin <
        wang_main_review_margin
    ] <- TRUE
  }


  if (!is.null(wang_subtype_review_margin)) {

    review_flag[
      immune_idx &
        final_subtype_margin <
        wang_subtype_review_margin
    ] <- TRUE
  }


  if (!is.null(xenium_subtype_review_margin)) {

    review_flag[
      nonimmune_idx &
        is.finite(final_subtype_margin) &
        final_subtype_margin <
        xenium_subtype_review_margin
    ] <- TRUE
  }


  # ============================================================
  # 18. QUALITATIVE CONFIDENCE
  # ============================================================
  #
  # Continuous score and margin columns remain the primary
  # quantitative evidence.
  #
  # These qualitative labels are deliberately conservative.
  # ============================================================

  confidence <- dplyr::case_when(

    review_flag ~
      "REVIEW",

    schwann_support |
      lymphatic_support ~
      "ONTOLOGY_GAP_SUPPORTED",

    broad_agreement ~
      "HIGH_CONCORDANT",

    immune_idx ~
      "REFERENCE_PRIMARY",

    TRUE ~
      "HYBRID_SUPPORTED"
  )


  # ============================================================
  # 19. WRITE EVIDENCE + FINAL LABELS INTO metadata
  # ============================================================

  # Wang kNN evidence
  md$WangKNN_main <-
    display_main(
      wang_main_top$label
    )

  md$WangKNN_main_score <-
    wang_main_top$score

  md$WangKNN_main_second <-
    display_main(
      wang_main_top$second
    )

  md$WangKNN_main_margin <-
    wang_main_top$margin

  md$WangKNN_immune_probability <-
    wang_immune_probability


  md$WangKNN_subtype <-
    display_subtype(
      wang_subtype_top$label
    )

  md$WangKNN_subtype_score <-
    wang_subtype_top$score

  md$WangKNN_subtype_margin <-
    wang_subtype_top$margin


  # Xenium kNN evidence
  md$XeniumKNN_main <-
    display_main(
      xenium_main_top$label
    )

  md$XeniumKNN_main_score <-
    xenium_main_top$score

  md$XeniumKNN_main_margin <-
    xenium_main_top$margin


  md$XeniumKNN_subtype <-
    display_subtype(
      xenium_subtype_top$label
    )

  md$XeniumKNN_subtype_score <-
    xenium_subtype_top$score

  md$XeniumKNN_subtype_margin <-
    xenium_subtype_top$margin


  # Ontology-gap information
  md$Ontology_gap_support <- dplyr::case_when(

    schwann_support ~
      "Schwann",

    lymphatic_support ~
      "Lymphatic_EC",

    TRUE ~
      NA_character_
  )


  # Final annotation
  md$Final_CellType_main <-
    display_main(
      final_main_raw
    )

  md$Final_CellType_subtype <-
    display_subtype(
      final_subtype_raw
    )

  md$Final_main_score <-
    final_main_score

  md$Final_main_margin <-
    final_main_margin

  md$Final_subtype_score <-
    final_subtype_score

  md$Final_subtype_margin <-
    final_subtype_margin

  md$Final_annotation_source <-
    final_source

  md$Final_annotation_confidence <-
    confidence

  md$Final_review_flag <-
    review_flag

  md$Final_annotation_reason <-
    final_reason


  # Restore metadata
  object@meta.data <- md


  # ============================================================
  # 20. SUMMARY OUTPUT
  # ============================================================

  summary_table <- md %>%
    tibble::rownames_to_column(
      "cell"
    ) %>%
    dplyr::count(
      Final_CellType_main,
      Final_CellType_subtype,
      Final_annotation_source,
      Final_annotation_confidence,
      Final_review_flag,
      name = "n_cells"
    ) %>%
    dplyr::arrange(
      Final_CellType_main,
      dplyr::desc(n_cells)
    )


  if (verbose) {

    message(
      "Annotation complete."
    )

    message(
      "Review cells: ",
      sum(review_flag),
      " / ",
      n_cells,
      " (",
      round(
        100 * mean(review_flag),
        2
      ),
      "%)"
    )

    message(
      "Schwann ontology-gap rescue: ",
      sum(schwann_support)
    )

    message(
      "Lymphatic EC ontology-gap rescue: ",
      sum(lymphatic_support)
    )
  }


  # ============================================================
  # Return both annotated object and useful diagnostics
  # ============================================================

  list(
    object = object,

    summary = summary_table,

    parameters = list(
      reduction = reduction,
      dims = dims,
      k = k,
      self_weight = self_weight,
      wang_prefix = wang_prefix,
      ontology_gap_max_immune_prob =
        ontology_gap_max_immune_prob
    )
  )
}

# KNOWN INVALID NOTEBOOK-LOCAL NEAREST-CELL ANALYSIS
# The zero-distance Eos error is not produced by a function in this source file.
# It occurs in notebooks/B1_Region3_primary_479.ipynb after rescued Eos are
# assigned through Final_CellType_subtype_refined, while the later non-Eos
# reference pool is filtered with the older Final_CellType_subtype column.
# Rescued Eos can therefore occur in both the Eos query and non-Eos reference
# pools, and FNN::get.knnx() returns zero-distance self-matches. Those legacy
# nearest-cell and neighbourhood outputs are not presently valid. A corrected
# analysis must freeze one refined label for both pools and verify
# length(intersect(eos_cell_ids, reference_cell_ids)) == 0 before any kNN call.

#' Refine eosinophil identity.
#'
#' @param object Required `object` input; validated before computation.
#' @param assay Optional `assay` input with the default shown in the function signature.
#' @param reduction Optional `reduction` input with the default shown in the function signature.
#' @param dims Optional `dims` input with the default shown in the function signature.
#' @param k Optional `k` input with the default shown in the function signature.
#' @param self_weight Optional `self_weight` input with the default shown in the function signature.
#' @param wang_main_eos_col Optional `wang_main_eos_col` input with the default shown in the function signature.
#' @param wang_subtype_eos_col Optional `wang_subtype_eos_col` input with the default shown in the function signature.
#' @param wang_subtype_prefix Optional `wang_subtype_prefix` input with the default shown in the function signature.
#' @param xenium_eos_score_col Optional `xenium_eos_score_col` input with the default shown in the function signature.
#' @param eos_core_genes Optional `eos_core_genes` input with the default shown in the function signature.
#' @param eos_support_genes Optional `eos_support_genes` input with the default shown in the function signature.
#' @param competitor_score_cols Optional `competitor_score_cols` input with the default shown in the function signature.
#' @param wang_high Optional `wang_high` input with the default shown in the function signature.
#' @param wang_support Optional `wang_support` input with the default shown in the function signature.
#' @param min_core_high Optional `min_core_high` input with the default shown in the function signature.
#' @param min_core_probable Optional `min_core_probable` input with the default shown in the function signature.
#' @param xenium_margin_high Optional `xenium_margin_high` input with the default shown in the function signature.
#' @param verbose Optional `verbose` input with the default shown in the function signature.
#' @return Computed evidence or annotations aligned to the supplied rows/cells; raw counts are unchanged.
refine_eosinophil_identity <- function(
  object,

  # ----------------------------------------------------------
  # Seurat / kNN settings
  # ----------------------------------------------------------
  assay = "Xenium",
  reduction = "pca",
  dims = 1:12,
  k = 15,
  self_weight = 0.70,

  # ----------------------------------------------------------
  # Wang reference columns
  # Uses pooled/all-age Wang reference only.
  # ----------------------------------------------------------
  wang_main_eos_col =
    "RefAll_main_prediction.score.Eosinophil",

  wang_subtype_eos_col =
    "RefAll_subtype_prediction.score.Eosinophil",

  wang_subtype_prefix =
    "RefAll_subtype_prediction.score.",

  # ----------------------------------------------------------
  # Xenium Eos identity score
  # ----------------------------------------------------------
  xenium_eos_score_col =
    "SubtypeScore_Eosinophil",

  # ----------------------------------------------------------
  # Direct Eos identity genes
  #
  # IMPORTANT:
  # These should be identity genes only.
  # Do NOT include short-lived / long-lived Eos state genes
  # used later for biological state analysis.
  # ----------------------------------------------------------
  eos_core_genes = c(
    "Siglecf",
    "Ccr3",
    "Il5ra",
    "Prg2",
    "Epx"
  ),

  eos_support_genes = c(
    "Alox15",
    "Ear1",
    "Ear2",
    "Ltc4s"
  ),

  # ----------------------------------------------------------
  # Main competing Xenium populations
  #
  # These are populations most likely to generate a false Eos
  # call in a targeted immune panel.
  # ----------------------------------------------------------
  competitor_score_cols = c(
    "SubtypeScore_MC",
    "SubtypeScore_Neutrophil",
    "SubtypeScore_Monocyte",
    "SubtypeScore_LYVE_macrophage",
    "SubtypeScore_Scavenging_macrophage",
    "SubtypeScore_DC",
    "SubtypeScore_cDC",
    "SubtypeScore_T.cell",
    "SubtypeScore_NK"
  ),

  # ----------------------------------------------------------
  # Starting thresholds for Eos confidence.
  #
  # These are deliberately configurable.
  # They should later be checked against the distributions
  # observed in the actual Xenium dataset.
  # ----------------------------------------------------------

  # Strong Wang Eos support
  wang_high = 0.50,

  # Some Wang Eos support
  wang_support = 0.20,

  # Minimum directly detected canonical Eos genes
  # for a high-confidence molecular call
  min_core_high = 2,

  # Minimum direct Eos genes for a probable call
  min_core_probable = 1,

  # Minimum Xenium Eos advantage over closest competitor
  # for a strongly specific Eos marker program.
  xenium_margin_high = 0.20,

  verbose = TRUE
) {

  # ==========================================================
  # 0. CHECK INPUTS
  # ==========================================================

  if (!inherits(object, "Seurat")) {
    stop("'object' must be a Seurat object.")
  }

  if (!requireNamespace("FNN", quietly = TRUE)) {
    stop(
      "Package 'FNN' is required. ",
      "Install with install.packages('FNN')."
    )
  }

  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package 'Matrix' is required.")
  }

  if (!reduction %in% Reductions(object)) {
    stop(
      "Reduction '",
      reduction,
      "' is not present in the object."
    )
  }

  if (!assay %in% Assays(object)) {
    stop(
      "Assay '",
      assay,
      "' is not present in the object."
    )
  }

  cells <- colnames(object)

  md <- object@meta.data[
    cells,
    ,
    drop = FALSE
  ]


  # ----------------------------------------------------------
  # Check essential metadata columns.
  # ----------------------------------------------------------

  required_cols <- c(
    wang_main_eos_col,
    wang_subtype_eos_col,
    xenium_eos_score_col
  )

  missing_required <- setdiff(
    required_cols,
    colnames(md)
  )

  if (length(missing_required) > 0) {
    stop(
      "Missing required metadata columns: ",
      paste(
        missing_required,
        collapse = ", "
      )
    )
  }


  # ==========================================================
  # 1. BUILD TRANSCRIPTOMIC kNN
  # ==========================================================
  #
  # These are PCA-expression neighbors, NOT spatial neighbors.
  #
  # We use them as conservative molecular context:
  #
  #   70% = cell's own evidence
  #   30% = transcriptionally similar neighbors
  #
  # This avoids replacing rare Eos cells by the local majority.
  # ==========================================================

  embedding <- Embeddings(
    object,
    reduction = reduction
  )

  embedding <- embedding[
    cells,
    ,
    drop = FALSE
  ]

  if (max(dims) > ncol(embedding)) {
    stop(
      "Requested PCA dimensions exceed those available."
    )
  }

  embedding <- embedding[
    ,
    dims,
    drop = FALSE
  ]

  if (verbose) {
    message(
      "Building transcriptomic kNN: k = ",
      k,
      "; PCs = ",
      paste(range(dims), collapse = ":")
    )
  }

  knn <- FNN::get.knn(
    embedding,
    k = k
  )


  # ==========================================================
  # 2. DISTANCE-WEIGHT THE NEIGHBORS
  # ==========================================================
  #
  # Closer transcriptomic neighbors contribute more strongly.
  # ==========================================================

  d <- knn$nn.dist

  sigma <- apply(
    d,
    1,
    median
  )

  sigma[
    !is.finite(sigma) |
      sigma <= 0
  ] <- 1

  weights <- exp(
    -(d^2) /
      (2 * sigma^2)
  )

  weights <- weights /
    rowSums(weights)


  # Sparse neighbor-weight matrix.
  W <- Matrix::sparseMatrix(
    i = rep(
      seq_len(nrow(embedding)),
      each = k
    ),
    j = as.vector(
      t(knn$nn.index)
    ),
    x = as.vector(
      t(weights)
    ),
    dims = c(
      nrow(embedding),
      nrow(embedding)
    )
  )


  # ----------------------------------------------------------
  # Helper:
  # smooth any continuous evidence matrix over kNN.
  # ----------------------------------------------------------

  smooth_knn <- function(x) {

    if (is.vector(x)) {
      x <- matrix(
        x,
        ncol = 1
      )
    }

    neighbour_signal <-
      as.matrix(
        W %*% x
      )

    out <-
      self_weight * x +
      (1 - self_weight) *
      neighbour_signal

    rownames(out) <- cells

    out
  }


  # ==========================================================
  # 3. WANG EOSINOPHIL EVIDENCE
  # ==========================================================
  #
  # Keep BOTH:
  #
  #   main Eosinophil probability
  #   subtype Eosinophil probability
  #
  # We do not simply use:
  #
  #   predicted.id == "Eosinophil"
  #
  # because a winning probability of 0.30 is very different
  # from a winning probability of 0.90.
  # ==========================================================

  wang_main_eos_raw <-
    as.numeric(
      md[[wang_main_eos_col]]
    )

  wang_subtype_eos_raw <-
    as.numeric(
      md[[wang_subtype_eos_col]]
    )


  wang_main_eos_knn <-
    as.numeric(
      smooth_knn(
        wang_main_eos_raw
      )[, 1]
    )

  wang_subtype_eos_knn <-
    as.numeric(
      smooth_knn(
        wang_subtype_eos_raw
      )[, 1]
    )


  # ==========================================================
  # 4. DETERMINE WANG'S BEST IMMUNE SUBTYPE
  # ==========================================================
  #
  # This tells us whether Eosinophil is actually the leading
  # Wang immune identity rather than merely having some
  # non-zero probability.
  # ==========================================================

  wang_subtype_cols <- grep(
    paste0(
      "^",
      gsub(
        "\\.",
        "\\\\.",
        wang_subtype_prefix
      )
    ),
    colnames(md),
    value = TRUE
  )

  wang_subtype_cols <- setdiff(
    wang_subtype_cols,
    paste0(
      sub(
        "\\.$",
        "",
        wang_subtype_prefix
      ),
      ".max"
    )
  )


  # Restrict this comparison to immune populations.
  wang_immune_subtypes <- c(
    "CCR2_inflammatory_Monocyte",
    "CX3CR1_Monocyte",
    "Inflammatory_Mac",
    "LYVE1_resident_Mac",
    "TREM2_LAM",
    "cDC1",
    "cDC2",
    "CCR7_migratory_DC",
    "Neutrophil",
    "Eosinophil",
    "Mast.cell",
    "T",
    "γδ.T",
    "NK",
    "ILC2",
    "B",
    "Plasma"
  )

  wang_immune_cols <- paste0(
    wang_subtype_prefix,
    wang_immune_subtypes
  )

  wang_immune_cols <- intersect(
    wang_immune_cols,
    colnames(md)
  )

  wang_immune_matrix <- as.matrix(
    md[
      ,
      wang_immune_cols,
      drop = FALSE
    ]
  )

  storage.mode(
    wang_immune_matrix
  ) <- "double"

  colnames(
    wang_immune_matrix
  ) <- sub(
    paste0(
      "^",
      gsub(
        "\\.",
        "\\\\.",
        wang_subtype_prefix
      )
    ),
    "",
    colnames(
      wang_immune_matrix
    )
  )

  wang_immune_knn <-
    smooth_knn(
      wang_immune_matrix
    )

  wang_immune_best_idx <-
    max.col(
      wang_immune_knn,
      ties.method = "first"
    )

  wang_immune_best <-
    colnames(
      wang_immune_knn
    )[
      wang_immune_best_idx
    ]

  wang_immune_best_score <-
    wang_immune_knn[
      cbind(
        seq_len(nrow(wang_immune_knn)),
        wang_immune_best_idx
      )
    ]

  wang_eos_is_best <-
    wang_immune_best ==
      "Eosinophil"


  # ==========================================================
  # 5. XENIUM EOS SCORE
  # ==========================================================
  #
  # Xenium Eos marker score is kept independent of Wang.
  # ==========================================================

  xenium_eos_raw <-
    as.numeric(
      md[[xenium_eos_score_col]]
    )

  xenium_eos_knn <-
    as.numeric(
      smooth_knn(
        xenium_eos_raw
      )[, 1]
    )


  # ==========================================================
  # 6. XENIUM COMPETING IMMUNE PROGRAMS
  # ==========================================================
  #
  # The important question is NOT simply:
  #
  #   "Is the Eos score positive?"
  #
  # but:
  #
  #   "Does the Eos program beat plausible competing
  #    immune identities?"
  #
  # ==========================================================

  competitor_score_cols_use <-
    intersect(
      competitor_score_cols,
      colnames(md)
    )

  if (!length(competitor_score_cols_use)) {
    stop(
      "None of the requested competitor-score columns ",
      "are present."
    )
  }

  competitor_matrix <- as.matrix(
    md[
      ,
      competitor_score_cols_use,
      drop = FALSE
    ]
  )

  storage.mode(
    competitor_matrix
  ) <- "double"

  colnames(
    competitor_matrix
  ) <- sub(
    "^SubtypeScore_",
    "",
    colnames(
      competitor_matrix
    )
  )

  competitor_knn <-
    smooth_knn(
      competitor_matrix
    )

  competitor_best_idx <-
    max.col(
      competitor_knn,
      ties.method = "first"
    )

  eos_competitor <-
    colnames(
      competitor_knn
    )[
      competitor_best_idx
    ]

  eos_competitor_score <-
    competitor_knn[
      cbind(
        seq_len(nrow(competitor_knn)),
        competitor_best_idx
      )
    ]


  # ----------------------------------------------------------
  # Positive value:
  # Eos marker program beats every tested competitor.
  #
  # Negative value:
  # another immune program is stronger.
  # ----------------------------------------------------------

  eos_vs_competitor_margin <-
    xenium_eos_knn -
    eos_competitor_score

  xenium_eos_is_best <-
    eos_vs_competitor_margin > 0


  # ==========================================================
  # 7. DIRECT EOS IDENTITY-GENE DETECTION
  # ==========================================================
  #
  # This is deliberately separate from module scores.
  #
  # It creates an easily interpretable variable such as:
  #
  #   3 / 5 core Eos genes detected.
  #
  # Raw counts are used because the biological question is
  # simply whether a transcript was detected in that cell.
  # ==========================================================

  old_assay <- DefaultAssay(object)

  DefaultAssay(object) <- assay

  eos_core_use <- intersect(
    eos_core_genes,
    rownames(object[[assay]])
  )

  eos_support_use <- intersect(
    eos_support_genes,
    rownames(object[[assay]])
  )


  if (verbose) {

    message(
      "Core Eos identity genes available: ",
      paste(
        eos_core_use,
        collapse = ", "
      )
    )

    message(
      "Support Eos identity genes available: ",
      paste(
        eos_support_use,
        collapse = ", "
      )
    )
  }


  # Core genes
  if (length(eos_core_use) > 0) {

    core_expr <- FetchData(
      object,
      vars = eos_core_use,
      layer = "counts"
    )

    core_expr <- core_expr[
      cells,
      ,
      drop = FALSE
    ]

    eos_core_n_detected <-
      rowSums(
        core_expr > 0
      )

  } else {

    eos_core_n_detected <-
      rep(
        0L,
        length(cells)
      )
  }


  # Support genes
  if (length(eos_support_use) > 0) {

    support_expr <- FetchData(
      object,
      vars = eos_support_use,
      layer = "counts"
    )

    support_expr <- support_expr[
      cells,
      ,
      drop = FALSE
    ]

    eos_support_n_detected <-
      rowSums(
        support_expr > 0
      )

  } else {

    eos_support_n_detected <-
      rep(
        0L,
        length(cells)
      )
  }


  DefaultAssay(object) <- old_assay


  # ==========================================================
  # 8. DEFINE INDEPENDENT SUPPORT FLAGS
  # ==========================================================
  #
  # Do NOT numerically average these scores together.
  #
  # Wang probabilities and Xenium z-score marker scores are
  # different quantities on different scales.
  #
  # We instead combine them logically.
  # ==========================================================


  # Strong external-reference evidence
  wang_high_support <-
    wang_eos_is_best &
    wang_subtype_eos_knn >=
      wang_high


  # Moderate external-reference evidence
  wang_some_support <-
    wang_subtype_eos_knn >=
      wang_support


  # Strong Xenium-specific evidence
  xenium_high_support <-
    xenium_eos_is_best &
    eos_vs_competitor_margin >=
      xenium_margin_high


  # Any Xenium-specific support
  xenium_some_support <-
    xenium_eos_is_best


  # Direct canonical transcript evidence
  direct_high_support <-
    eos_core_n_detected >=
      min_core_high

  direct_some_support <-
    eos_core_n_detected >=
      min_core_probable


  # ==========================================================
  # 9. EOSINOPHIL CONFIDENCE CLASSIFICATION
  # ==========================================================
  #
  # The hierarchy is intentionally conservative because Eos
  # identity drives the downstream biological conclusions.
  #
  # HIGH_CONFIDENCE_EOS:
  #
  #   external reference agrees
  #   +
  #   Xenium-specific Eos program wins
  #   +
  #   multiple canonical Eos genes detected
  #
  #
  # PROBABLE_EOS:
  #
  #   strong evidence from two systems but one component is
  #   weakened by expected targeted-panel dropout.
  #
  #
  # AMBIGUOUS_EOS:
  #
  #   some Eos evidence exists, but independent evidence
  #   disagrees or is insufficient.
  #
  #
  # NON_EOS:
  #
  #   no compelling Eos evidence.
  # ==========================================================

  eos_call <- rep(
    "NON_EOS",
    length(cells)
  )


  # ----------------------------------------------------------
  # HIGH-CONFIDENCE EOS
  #
  # Require concordance across:
  #
  #   1. Wang
  #   2. Xenium marker program
  #   3. direct identity-gene detection
  # ----------------------------------------------------------

  high_confidence <-
    wang_high_support &
    xenium_high_support &
    direct_high_support

  eos_call[
    high_confidence
  ] <- "HIGH_CONFIDENCE_EOS"


  # ----------------------------------------------------------
  # PROBABLE EOS
  #
  # Route A:
  # Strong Wang + Eos marker program + >=1 core gene.
  #
  # This allows some transcript dropout.
  # ----------------------------------------------------------

  probable_A <-
    !high_confidence &
    wang_high_support &
    xenium_some_support &
    direct_some_support


  # ----------------------------------------------------------
  # Route B:
  # Strong Xenium/direct molecular evidence,
  # but Wang mapping is weaker.
  #
  # This avoids forcing external-reference false negatives.
  # ----------------------------------------------------------

  probable_B <-
    !high_confidence &
    xenium_high_support &
    direct_high_support &
    wang_some_support


  probable <-
    probable_A |
    probable_B

  eos_call[
    probable
  ] <- "PROBABLE_EOS"


  # ----------------------------------------------------------
  # AMBIGUOUS EOS CANDIDATE
  #
  # Any meaningful Eos signal that did not satisfy the
  # concordant criteria above is retained for review rather
  # than silently classified as non-Eos.
  # ----------------------------------------------------------

  ambiguous <-
    !high_confidence &
    !probable &
    (
      wang_eos_is_best |
      xenium_eos_is_best |
      direct_high_support
    )

  eos_call[
    ambiguous
  ] <- "AMBIGUOUS_EOS"


  # ==========================================================
  # 10. GENERATE HUMAN-READABLE REVIEW REASONS
  # ==========================================================

  eos_review_reason <- rep(
    NA_character_,
    length(cells)
  )


  eos_review_reason[
    high_confidence
  ] <- paste0(
    "Concordant Wang reference, Xenium Eos marker program, ",
    "and direct canonical Eos-gene detection."
  )


  eos_review_reason[
    probable_A
  ] <- paste0(
    "Strong Wang and Xenium Eos evidence with limited direct ",
    "canonical-gene detection, compatible with transcript dropout."
  )


  eos_review_reason[
    probable_B
  ] <- paste0(
    "Strong Xenium marker and direct canonical-gene evidence; ",
    "Wang Eos probability is supportive but weaker."
  )


  # More detailed ambiguous reasons

  idx <- ambiguous &
    wang_eos_is_best &
    !xenium_eos_is_best

  eos_review_reason[idx] <- paste0(
    "Wang favors Eosinophil, but a competing Xenium immune ",
    "marker program is stronger."
  )


  idx <- ambiguous &
    !wang_eos_is_best &
    xenium_eos_is_best

  eos_review_reason[idx] <- paste0(
    "Xenium favors Eosinophil, but Wang reference does not; ",
    "requires review."
  )


  idx <- ambiguous &
    direct_high_support &
    !wang_eos_is_best &
    !xenium_eos_is_best

  eos_review_reason[idx] <- paste0(
    "Multiple canonical Eos genes are detected, but neither ",
    "Wang nor the Xenium score classifier selects Eosinophil."
  )


  # ==========================================================
  # 11. ADD ALL EOS-SPECIFIC EVIDENCE TO METADATA
  # ==========================================================

  md$Eos_Wang_main_probability_raw <-
    wang_main_eos_raw

  md$Eos_Wang_main_probability_knn <-
    wang_main_eos_knn

  md$Eos_Wang_subtype_probability_raw <-
    wang_subtype_eos_raw

  md$Eos_Wang_subtype_probability_knn <-
    wang_subtype_eos_knn

  md$Eos_Wang_best_immune_subtype <-
    wang_immune_best

  md$Eos_Wang_best_immune_score <-
    wang_immune_best_score

  md$Eos_Wang_is_best <-
    wang_eos_is_best


  md$Eos_Xenium_score_raw <-
    xenium_eos_raw

  md$Eos_Xenium_score_knn <-
    xenium_eos_knn

  md$Eos_competitor <-
    eos_competitor

  md$Eos_competitor_score <-
    eos_competitor_score

  md$Eos_vs_competitor_margin <-
    eos_vs_competitor_margin

  md$Eos_Xenium_is_best <-
    xenium_eos_is_best


  md$Eos_core_n_detected <-
    eos_core_n_detected

  md$Eos_core_fraction_detected <-
    if (length(eos_core_use) > 0) {
      eos_core_n_detected /
        length(eos_core_use)
    } else {
      NA_real_
    }

  md$Eos_support_n_detected <-
    eos_support_n_detected


  md$Eos_call <-
    factor(
      eos_call,
      levels = c(
        "HIGH_CONFIDENCE_EOS",
        "PROBABLE_EOS",
        "AMBIGUOUS_EOS",
        "NON_EOS"
      )
    )

  md$Eos_primary_include <-
    eos_call ==
      "HIGH_CONFIDENCE_EOS"

  md$Eos_sensitivity_include <-
    eos_call %in%
      c(
        "HIGH_CONFIDENCE_EOS",
        "PROBABLE_EOS"
      )

  md$Eos_manual_review <-
    eos_call ==
      "AMBIGUOUS_EOS"

  md$Eos_review_reason <-
    eos_review_reason


  object@meta.data <- md


  # ==========================================================
  # 12. SUMMARY TABLE
  # ==========================================================

  summary <- md %>%
    tibble::rownames_to_column(
      "cell"
    ) %>%
    dplyr::count(
      Eos_call,
      name = "n_cells"
    ) %>%
    dplyr::mutate(
      fraction =
        n_cells /
        sum(n_cells)
    )


  # ==========================================================
  # 13. COMPETITOR SUMMARY FOR EOS CANDIDATES
  # ==========================================================

  candidate_summary <- md %>%
    tibble::rownames_to_column(
      "cell"
    ) %>%
    dplyr::filter(
      Eos_call !=
        "NON_EOS"
    ) %>%
    dplyr::count(
      Eos_call,
      Eos_competitor,
      name = "n_cells"
    ) %>%
    dplyr::group_by(
      Eos_call
    ) %>%
    dplyr::mutate(
      fraction =
        n_cells /
        sum(n_cells)
    ) %>%
    dplyr::ungroup()


  if (verbose) {

    message(
      "Eosinophil refinement complete."
    )

    message(
      "HIGH_CONFIDENCE_EOS: ",
      sum(
        eos_call ==
          "HIGH_CONFIDENCE_EOS"
      )
    )

    message(
      "PROBABLE_EOS: ",
      sum(
        eos_call ==
          "PROBABLE_EOS"
      )
    )

    message(
      "AMBIGUOUS_EOS: ",
      sum(
        eos_call ==
          "AMBIGUOUS_EOS"
      )
    )

    message(
      "Core identity genes used: ",
      paste(
        eos_core_use,
        collapse = ", "
      )
    )
  }


  # ==========================================================
  # 14. RETURN
  # ==========================================================

  list(
    object = object,

    summary = summary,

    candidate_summary =
      candidate_summary,

    genes = list(
      core = eos_core_use,
      support = eos_support_use
    ),

    parameters = list(
      reduction = reduction,
      dims = dims,
      k = k,
      self_weight = self_weight,
      wang_high = wang_high,
      wang_support = wang_support,
      min_core_high = min_core_high,
      min_core_probable =
        min_core_probable,
      xenium_margin_high =
        xenium_margin_high
    )
  )
}

#' Score eosinophil likeness.
#'
#' @param reference Required `reference` input; validated before computation.
#' @param query Required `query` input; validated before computation.
#' @param reference_group_col Optional `reference_group_col` input with the default shown in the function signature.
#' @param eos_label Optional `eos_label` input with the default shown in the function signature.
#' @param reference_sample_col Optional `reference_sample_col` input with the default shown in the function signature.
#' @param reference_assay Optional `reference_assay` input with the default shown in the function signature.
#' @param query_assay Optional `query_assay` input with the default shown in the function signature.
#' @param core_markers Optional `core_markers` input with the default shown in the function signature.
#' @param support_markers Optional `support_markers` input with the default shown in the function signature.
#' @param tier3_marker Optional `tier3_marker` input with the default shown in the function signature.
#' @param context_markers Optional `context_markers` input with the default shown in the function signature.
#' @param reference_immune_labels Optional `reference_immune_labels` input with the default shown in the function signature.
#' @param query_main_col Optional `query_main_col` input with the default shown in the function signature.
#' @param query_immune_labels Optional `query_immune_labels` input with the default shown in the function signature.
#' @param wang_predicted_col Optional `wang_predicted_col` input with the default shown in the function signature.
#' @param wang_eos_score_col Optional `wang_eos_score_col` input with the default shown in the function signature.
#' @param min_core_eos_cells Optional `min_core_eos_cells` input with the default shown in the function signature.
#' @param min_core_eos_pct Optional `min_core_eos_pct` input with the default shown in the function signature.
#' @param min_pair_eos_cells Optional `min_pair_eos_cells` input with the default shown in the function signature.
#' @param min_pair_eos_pct Optional `min_pair_eos_pct` input with the default shown in the function signature.
#' @param min_competitor_cells Optional `min_competitor_cells` input with the default shown in the function signature.
#' @param smoothing Optional `smoothing` input with the default shown in the function signature.
#' @param allow_nonimmune_rescue Optional `allow_nonimmune_rescue` input with the default shown in the function signature.
#' @param verbose Optional `verbose` input with the default shown in the function signature.
#' @return Computed evidence or annotations aligned to the supplied rows/cells; raw counts are unchanged.
score_eosinophil_likeness <- function(

  reference,
  query,

  # ==========================================================
  # REFERENCE ANNOTATION
  # ==========================================================

  reference_group_col = "CellType_subtype",
  eos_label = "Eosinophil",

  # Retained for compatibility / possible future
  # sample-wise validation.
  reference_sample_col = NULL,


  # ==========================================================
  # ASSAYS
  # ==========================================================

  reference_assay = "RNA",
  query_assay = "Xenium",


  # ==========================================================
  # TIER DEFINITIONS
  # ==========================================================

  # Tier 1:
  # Siglecf + Ccr3
  core_markers = c(
    "Siglecf",
    "Ccr3"
  ),

  # Tier 2:
  # one core marker + Il5ra/Alox15
  support_markers = c(
    "Il5ra",
    "Alox15"
  ),

  # Tier 3:
  # one core marker + Itgam
  tier3_marker = "Itgam",

  # Additional markers contributing to Tier 4.
  context_markers = c(
    "Ear1",
    "Ear2"
  ),


  # ==========================================================
  # REFERENCE IMMUNE POPULATIONS
  # ==========================================================

  reference_immune_labels = c(
    "Eosinophil",

    "LYVE1_resident_Mac",
    "Inflammatory_Mac",
    "TREM2_LAM",

    "CCR2_inflammatory_Monocyte",
    "CX3CR1_Monocyte",

    "Neutrophil",
    "Mast cell",

    "cDC1",
    "cDC2",
    "CCR7_migratory_DC",

    "NK",
    "T",
    "γδ T",
    "B",
    "Plasma",
    "ILC2"
  ),


  # ==========================================================
  # QUERY IMMUNE DEFINITION
  # ==========================================================

  query_main_col = "Final_CellType_main",

  query_immune_labels = c(
    "Macrophage",
    "Monocyte",
    "Neutrophil",
    "Mast cell",
    "Eosinophil",
    "T cell",
    "NK",
    "DC",
    "B cell",
    "Plasma",
    "ILC"
  ),


  # ==========================================================
  # WANG EVIDENCE
  #
  # Kept completely separate from Tier assignment.
  # ==========================================================

  wang_predicted_col =
    "RefAll_subtype_predicted.id",

  wang_eos_score_col =
    "RefAll_subtype_prediction.score.Eosinophil",


  # ==========================================================
  # REFERENCE VALIDATION
  # ==========================================================

  min_core_eos_cells = 3,
  min_core_eos_pct = 5,

  min_pair_eos_cells = 3,
  min_pair_eos_pct = 1,

  min_competitor_cells = 20,

  smoothing = 0.5,


  # Nonimmune Tier-1 + Wang-Eos cells are surfaced
  # for review, never automatically relabelled.
  allow_nonimmune_rescue = TRUE,

  verbose = TRUE
) {


  # ==========================================================
  # 0. BASIC CHECKS
  # ==========================================================

  if (!inherits(reference, "Seurat")) {
    stop("'reference' must be a Seurat object.")
  }

  if (!inherits(query, "Seurat")) {
    stop("'query' must be a Seurat object.")
  }

  if (!reference_assay %in% Assays(reference)) {
    stop(
      "Reference assay '",
      reference_assay,
      "' was not found."
    )
  }

  if (!query_assay %in% Assays(query)) {
    stop(
      "Query assay '",
      query_assay,
      "' was not found."
    )
  }

  if (length(core_markers) != 2) {
    stop(
      "Tier 1 requires exactly two core markers."
    )
  }

  if (length(tier3_marker) != 1) {
    stop(
      "'tier3_marker' must contain exactly one marker."
    )
  }


  # ==========================================================
  # 1. ALIGN METADATA
  # ==========================================================

  ref_cells <- colnames(reference)
  query_cells <- colnames(query)

  ref_md <- reference@meta.data[
    ref_cells,
    ,
    drop = FALSE
  ]

  query_md <- query@meta.data[
    query_cells,
    ,
    drop = FALSE
  ]

  if (!reference_group_col %in% colnames(ref_md)) {
    stop(
      "Reference metadata column '",
      reference_group_col,
      "' was not found."
    )
  }

  if (!query_main_col %in% colnames(query_md)) {
    stop(
      "Query metadata column '",
      query_main_col,
      "' was not found."
    )
  }

  if (
    !is.null(reference_sample_col) &&
    !reference_sample_col %in% colnames(ref_md)
  ) {
    stop(
      "Reference sample column '",
      reference_sample_col,
      "' was not found."
    )
  }


  # ==========================================================
  # 2. MARKERS
  # ==========================================================

  core_markers <- unique(core_markers)
  support_markers <- unique(support_markers)
  context_markers <- unique(context_markers)

  all_requested_markers <- unique(
    c(
      core_markers,
      support_markers,
      tier3_marker,
      context_markers
    )
  )

  ref_features <- rownames(
    reference[[reference_assay]]
  )

  query_features <- rownames(
    query[[query_assay]]
  )

  all_markers_use <- all_requested_markers[
    all_requested_markers %in% ref_features &
      all_requested_markers %in% query_features
  ]

  missing_markers <- setdiff(
    all_requested_markers,
    all_markers_use
  )


  # These are required because they explicitly define
  # Tier 1-3.
  required_tier_markers <- unique(
    c(
      core_markers,
      support_markers,
      tier3_marker
    )
  )

  missing_required <- setdiff(
    required_tier_markers,
    all_markers_use
  )

  if (length(missing_required) > 0) {
    stop(
      "Markers required for Tier 1-3 are missing: ",
      paste(
        missing_required,
        collapse = ", "
      )
    )
  }


  core_markers_use <- core_markers
  support_markers_use <- support_markers

  context_markers_use <- intersect(
    context_markers,
    all_markers_use
  )


  if (verbose) {

    message(
      "Tier 1 core markers: ",
      paste(
        core_markers_use,
        collapse = ", "
      )
    )

    message(
      "Tier 2 support markers: ",
      paste(
        support_markers_use,
        collapse = ", "
      )
    )

    message(
      "Tier 3 marker: ",
      tier3_marker
    )

    message(
      "Tier 4 additional markers: ",
      ifelse(
        length(context_markers_use) > 0,
        paste(
          context_markers_use,
          collapse = ", "
        ),
        "None"
      )
    )

    if (length(missing_markers) > 0) {
      message(
        "Unavailable optional markers: ",
        paste(
          missing_markers,
          collapse = ", "
        )
      )
    }
  }


  # ==========================================================
  # 3. RAW COUNTS -> DETECTION
  #
  # Cross-platform rule:
  #
  # count > 0
  # ==========================================================

  old_ref_assay <- DefaultAssay(reference)
  old_query_assay <- DefaultAssay(query)

  DefaultAssay(reference) <- reference_assay
  DefaultAssay(query) <- query_assay

  ref_counts <- FetchData(
    reference,
    vars = all_markers_use,
    layer = "counts"
  )

  query_counts <- FetchData(
    query,
    vars = all_markers_use,
    layer = "counts"
  )

  DefaultAssay(reference) <- old_ref_assay
  DefaultAssay(query) <- old_query_assay


  ref_counts <- ref_counts[
    ref_cells,
    all_markers_use,
    drop = FALSE
  ]

  query_counts <- query_counts[
    query_cells,
    all_markers_use,
    drop = FALSE
  ]


  ref_binary <- as.matrix(
    ref_counts > 0
  )

  query_binary <- as.matrix(
    query_counts > 0
  )


  # ==========================================================
  # 4. REFERENCE IMMUNE CELLS
  # ==========================================================

  ref_celltype <- as.character(
    ref_md[[reference_group_col]]
  )

  ref_is_immune <- ref_celltype %in%
    reference_immune_labels

  ref_is_eos <- ref_celltype ==
    eos_label

  keep_ref <- (
    ref_is_immune &
      !is.na(ref_celltype)
  )

  ref_binary_immune <- ref_binary[
    keep_ref,
    ,
    drop = FALSE
  ]

  ref_celltype_immune <- ref_celltype[
    keep_ref
  ]

  ref_is_eos_immune <- ref_is_eos[
    keep_ref
  ]

  n_ref_eos <- sum(
    ref_is_eos_immune
  )

  n_ref_other <- sum(
    !ref_is_eos_immune
  )

  if (n_ref_eos == 0) {
    stop(
      "No reference Eosinophils remain after immune filtering."
    )
  }

  if (n_ref_other == 0) {
    stop(
      "No non-Eosinophil immune reference cells remain."
    )
  }


  # ==========================================================
  # 5. COMPETING REFERENCE CELL TYPES
  # ==========================================================

  get_competing_types <- function(
    celltype,
    eos_status
  ) {

    competitor_table <- table(
      celltype[
        !eos_status &
          !is.na(celltype)
      ]
    )

    competing_types <- names(
      competitor_table[
        competitor_table >= min_competitor_cells
      ]
    )

    if (length(competing_types) == 0) {
      competing_types <- names(
        competitor_table
      )
    }

    competing_types
  }


  # ==========================================================
  # 6. REFERENCE SIGNAL STATISTICS
  # ==========================================================

  summarize_signal <- function(
    signal,
    eos_status,
    celltype
  ) {

    eos_signal <- signal[
      eos_status
    ]

    eos_n <- length(
      eos_signal
    )

    eos_positive <- sum(
      eos_signal,
      na.rm = TRUE
    )

    eos_pct <- 100 *
      eos_positive /
      eos_n


    competing_types <- get_competing_types(
      celltype,
      eos_status
    )


    competitor_stats <- lapply(
      competing_types,
      function(ct) {

        idx <- (
          !eos_status &
            celltype == ct
        )

        n_ct <- sum(
          idx
        )

        positive_ct <- sum(
          signal[idx],
          na.rm = TRUE
        )

        pct_ct <- 100 *
          positive_ct /
          n_ct

        p_ct <- (
          positive_ct +
            smoothing
        ) / (
          n_ct +
            2 * smoothing
        )

        data.frame(
          CellType = ct,
          n = n_ct,
          positive = positive_ct,
          pct = pct_ct,
          p_smoothed = p_ct,
          stringsAsFactors = FALSE
        )
      }
    )

    competitor_stats <- do.call(
      rbind,
      competitor_stats
    )


    other_macro_pct <- mean(
      competitor_stats$pct,
      na.rm = TRUE
    )

    p_other_macro <- mean(
      competitor_stats$p_smoothed,
      na.rm = TRUE
    )

    worst_idx <- which.max(
      competitor_stats$p_smoothed
    )

    max_other_pct <-
      competitor_stats$pct[
        worst_idx
      ]

    max_other_celltype <-
      competitor_stats$CellType[
        worst_idx
      ]


    p_eos <- (
      eos_positive +
        smoothing
    ) / (
      eos_n +
        2 * smoothing
    )


    eps <- 1e-8

    p_eos <- pmin(
      pmax(p_eos, eps),
      1 - eps
    )

    p_other_macro <- pmin(
      pmax(p_other_macro, eps),
      1 - eps
    )


    log2_or_macro <- log2(
      (
        p_eos /
          (1 - p_eos)
      ) /
        (
          p_other_macro /
            (1 - p_other_macro)
        )
    )


    list(
      n_Eos = eos_n,
      Eos_positive = eos_positive,
      Eos_pct = eos_pct,
      Other_macro_pct = other_macro_pct,
      Max_other_pct = max_other_pct,
      Max_other_celltype = max_other_celltype,
      log2_detection_OR_macro = log2_or_macro
    )
  }


  # ==========================================================
  # 7. INDIVIDUAL MARKER STATISTICS
  # ==========================================================

  marker_stats <- do.call(
    rbind,
    lapply(
      all_markers_use,
      function(g) {

        stats <- summarize_signal(
          ref_binary_immune[, g],
          ref_is_eos_immune,
          ref_celltype_immune
        )

        marker_role <- dplyr::case_when(
          g %in% core_markers_use ~
            "Tier1_core",

          g %in% support_markers_use ~
            "Tier2_support",

          g == tier3_marker ~
            "Tier3_support",

          TRUE ~
            "Tier4_context"
        )

        data.frame(
          Gene = g,
          Marker_role = marker_role,

          n_Eos =
            stats$n_Eos,

          Eos_positive =
            stats$Eos_positive,

          Eos_pct =
            stats$Eos_pct,

          Other_macro_pct =
            stats$Other_macro_pct,

          Max_other_pct =
            stats$Max_other_pct,

          Max_other_celltype =
            stats$Max_other_celltype,

          log2_detection_OR_macro =
            stats$log2_detection_OR_macro,

          stringsAsFactors = FALSE
        )
      }
    )
  )


  # ==========================================================
  # 8. CORE MARKER VALIDATION
  # ==========================================================

  marker_stats$core_anchor_valid <- (
    marker_stats$Marker_role ==
      "Tier1_core" &
      marker_stats$Eos_positive >=
      min_core_eos_cells &
      marker_stats$Eos_pct >=
      min_core_eos_pct &
      marker_stats$log2_detection_OR_macro >
      0 &
      marker_stats$Eos_pct >
      marker_stats$Max_other_pct
  )


  core_stats <- marker_stats[
    marker_stats$Marker_role ==
      "Tier1_core",
    ,
    drop = FALSE
  ]


  if (!all(core_stats$core_anchor_valid)) {
    warning(
      "One or more Tier-1 core markers do not satisfy ",
      "the requested reference-specificity criteria. ",
      "Tier calls will still be generated, but inspect marker_stats."
    )
  }


  # ==========================================================
  # 9. CORE SCORE
  #
  # Diagnostic / ranking only.
  #
  # It does NOT determine the Tier.
  # ==========================================================

  core_specificity <- pmax(
    core_stats$log2_detection_OR_macro,
    0
  )

  if (sum(core_specificity) > 0) {

    core_weights <- core_specificity /
      sum(core_specificity)

  } else {

    core_weights <- rep(
      1 / length(core_markers_use),
      length(core_markers_use)
    )
  }

  names(core_weights) <- core_stats$Gene


  marker_stats$core_weight <- 0

  marker_stats$core_weight[
    match(
      names(core_weights),
      marker_stats$Gene
    )
  ] <- core_weights


  # ==========================================================
  # 10. DIAGNOSTIC PAIR STATISTICS
  #
  # IMPORTANT:
  #
  # These statistics NO LONGER determine Tier membership.
  #
  # They only tell us how the fixed Tier2 / Tier3
  # combinations perform in the reference.
  # ==========================================================

  tier2_grid <- expand.grid(
    Core_marker = core_markers_use,
    Partner_marker = support_markers_use,
    Tier_rule = "Tier2",
    stringsAsFactors = FALSE
  )

  tier3_grid <- data.frame(
    Core_marker = core_markers_use,
    Partner_marker = tier3_marker,
    Tier_rule = "Tier3",
    stringsAsFactors = FALSE
  )

  pair_grid <- rbind(
    tier2_grid,
    tier3_grid
  )


  pair_stats <- do.call(
    rbind,
    lapply(
      seq_len(nrow(pair_grid)),
      function(i) {

        core_gene <-
          pair_grid$Core_marker[i]

        partner_gene <-
          pair_grid$Partner_marker[i]


        signal <- (
          ref_binary_immune[, core_gene] &
            ref_binary_immune[, partner_gene]
        )


        stats <- summarize_signal(
          signal,
          ref_is_eos_immune,
          ref_celltype_immune
        )


        reference_supported <- (
          stats$Eos_positive >=
            min_pair_eos_cells &
            stats$Eos_pct >=
            min_pair_eos_pct &
            stats$log2_detection_OR_macro >
            0 &
            stats$Eos_pct >
            stats$Max_other_pct
        )


        data.frame(
          Tier_rule =
            pair_grid$Tier_rule[i],

          Core_marker =
            core_gene,

          Partner_marker =
            partner_gene,

          Pair =
            paste0(
              core_gene,
              " + ",
              partner_gene
            ),

          n_Eos =
            stats$n_Eos,

          Eos_positive =
            stats$Eos_positive,

          Eos_pct =
            stats$Eos_pct,

          Other_macro_pct =
            stats$Other_macro_pct,

          Max_other_pct =
            stats$Max_other_pct,

          Max_other_celltype =
            stats$Max_other_celltype,

          log2_detection_OR_macro =
            stats$log2_detection_OR_macro,

          reference_supported =
            reference_supported,

          stringsAsFactors = FALSE
        )
      }
    )
  )


  # ==========================================================
  # 11. FIXED HIERARCHICAL TIER RULE
  # ==========================================================

  apply_tier_rule <- function(
    binary_matrix
  ) {

    n_cells <- nrow(
      binary_matrix
    )


    # --------------------------------------------------------
    # Marker detection components
    # --------------------------------------------------------

    core_binary <- binary_matrix[
      ,
      core_markers_use,
      drop = FALSE
    ]


    n_core_detected <- rowSums(
      core_binary
    )


    any_core <- (
      n_core_detected >= 1
    )


    both_core <- (
      n_core_detected == 2
    )


    n_tier2_detected <- rowSums(
      binary_matrix[
        ,
        support_markers_use,
        drop = FALSE
      ]
    )


    any_tier2_support <- (
      n_tier2_detected >= 1
    )


    tier3_detected <- as.logical(
      binary_matrix[
        ,
        tier3_marker
      ]
    )


    n_all_markers <- rowSums(
      binary_matrix[
        ,
        all_markers_use,
        drop = FALSE
      ]
    )


    # --------------------------------------------------------
    # TIER 1
    #
    # Siglecf + Ccr3
    # --------------------------------------------------------

    tier1 <- both_core


    # --------------------------------------------------------
    # TIER 2
    #
    # one core
    # +
    # Il5ra and/or Alox15
    #
    # Tier1 always takes precedence.
    # --------------------------------------------------------

    tier2 <- (
      !tier1 &
        any_core &
        any_tier2_support
    )


    # --------------------------------------------------------
    # TIER 3
    #
    # one core
    # +
    # Itgam
    #
    # Tier1 and Tier2 take precedence.
    # --------------------------------------------------------

    tier3 <- (
      !tier1 &
        !tier2 &
        any_core &
        tier3_detected
    )


    # --------------------------------------------------------
    # TIER 4
    #
    # Any remaining >=2-marker combination.
    # --------------------------------------------------------

    tier4 <- (
      !tier1 &
        !tier2 &
        !tier3 &
        n_all_markers >= 2
    )


    # --------------------------------------------------------
    # FINAL MUTUALLY EXCLUSIVE TIER
    # --------------------------------------------------------

    tier <- dplyr::case_when(

      tier1 ~
        "Conf Tier1",

      tier2 ~
        "Conf Tier2",

      tier3 ~
        "Conf Tier3",

      tier4 ~
        "Conf Tier4",

      TRUE ~
        "Rest"
    )


    tier_reason <- dplyr::case_when(

      tier1 ~
        "Siglecf + Ccr3",

      tier2 ~
        "Core + Il5ra/Alox15",

      tier3 ~
        "Core + Itgam",

      tier4 ~
        "Other >=2-marker combination",

      TRUE ~
        "No qualifying >=2-marker combination"
    )


    # --------------------------------------------------------
    # Exact detected-marker combination
    # --------------------------------------------------------

    marker_combination <- apply(
      binary_matrix[
        ,
        all_markers_use,
        drop = FALSE
      ],
      1,
      function(x) {

        genes <- all_markers_use[
          as.logical(x)
        ]

        if (length(genes) == 0) {

          "None"

        } else {

          paste(
            genes,
            collapse = " + "
          )
        }
      }
    )


    core_pattern <- apply(
      core_binary,
      1,
      function(x) {

        genes <- core_markers_use[
          as.logical(x)
        ]

        if (length(genes) == 0) {

          "None"

        } else {

          paste(
            genes,
            collapse = " + "
          )
        }
      }
    )


    data.frame(

      EosRef_tier =
        tier,

      EosRef_tier_reason =
        tier_reason,

      EosRef_n_core_detected =
        n_core_detected,

      EosRef_core_pattern =
        core_pattern,

      EosRef_n_tier2_support_detected =
        n_tier2_detected,

      EosRef_tier3_marker_detected =
        tier3_detected,

      EosRef_n_all_markers_detected =
        n_all_markers,

      EosRef_marker_combination =
        marker_combination,

      stringsAsFactors = FALSE
    )
  }


  # ==========================================================
  # 12. APPLY TIER RULE TO REFERENCE
  # ==========================================================

  reference_rule_all <- apply_tier_rule(
    ref_binary
  )


  reference_rule_immune <- reference_rule_all[
    keep_ref,
    ,
    drop = FALSE
  ]


  # ==========================================================
  # 13. REFERENCE RULE PERFORMANCE
  #
  # CUMULATIVE logic:
  #
  # Tier1
  # Tier1 + Tier2
  # Tier1 + Tier2 + Tier3
  # Tier1 + Tier2 + Tier3 + Tier4
  # ==========================================================

  evaluate_rule <- function(
    predicted,
    eos_status,
    celltype
  ) {

    sensitivity <- mean(
      predicted[eos_status]
    )

    pooled_fpr <- mean(
      predicted[!eos_status]
    )


    competing_types <- get_competing_types(
      celltype,
      eos_status
    )


    lineage_fpr <- sapply(
      competing_types,
      function(ct) {

        idx <- (
          !eos_status &
            celltype == ct
        )

        mean(
          predicted[idx]
        )
      }
    )


    macro_fpr <- mean(
      lineage_fpr
    )

    max_fpr <- max(
      lineage_fpr
    )

    worst_lineage <- names(
      lineage_fpr
    )[
      which.max(
        lineage_fpr
      )
    ]


    data.frame(

      sensitivity =
        sensitivity,

      pooled_false_positive_rate =
        pooled_fpr,

      macro_false_positive_rate =
        macro_fpr,

      macro_specificity =
        1 - macro_fpr,

      max_lineage_false_positive_rate =
        max_fpr,

      worst_false_positive_lineage =
        worst_lineage,

      stringsAsFactors = FALSE
    )
  }


  reference_tier <- reference_rule_immune$EosRef_tier


  cumulative_rules <- list(

    "TIER1_ONLY" =
      reference_tier ==
      "Conf Tier1",

    "TIER1_TO_TIER2" =
      reference_tier %in%
      c(
        "Conf Tier1",
        "Conf Tier2"
      ),

    "TIER1_TO_TIER3" =
      reference_tier %in%
      c(
        "Conf Tier1",
        "Conf Tier2",
        "Conf Tier3"
      ),

    "TIER1_TO_TIER4" =
      reference_tier %in%
      c(
        "Conf Tier1",
        "Conf Tier2",
        "Conf Tier3",
        "Conf Tier4"
      )
  )


  reference_rule_performance <- do.call(
    rbind,
    lapply(
      names(cumulative_rules),
      function(rule_name) {

        cbind(

          data.frame(
            Rule = rule_name,
            calibration_method =
              "FIXED_HIERARCHICAL_TIER_RULE",
            stringsAsFactors = FALSE
          ),

          evaluate_rule(
            cumulative_rules[[rule_name]],
            ref_is_eos_immune,
            ref_celltype_immune
          )
        )
      }
    )
  )


  # ==========================================================
  # 14. REFERENCE TIER DISTRIBUTION
  #
  # Useful directly for your stacked barplot.
  # ==========================================================

  reference_tier_distribution <- data.frame(

    CellType =
      ref_celltype,

    EosRef_tier =
      reference_rule_all$EosRef_tier,

    stringsAsFactors = FALSE
  ) %>%

    dplyr::count(
      CellType,
      EosRef_tier,
      name = "n"
    ) %>%

    dplyr::group_by(
      CellType
    ) %>%

    dplyr::mutate(

      n_total =
        sum(n),

      fraction =
        n / n_total
    ) %>%

    dplyr::ungroup()


  # ==========================================================
  # 15. APPLY TIERS TO QUERY
  # ==========================================================

  query_rule <- apply_tier_rule(
    query_binary
  )


  query_is_immune <- as.character(
    query_md[[query_main_col]]
  ) %in%
    query_immune_labels


  # ==========================================================
  # 16. FINAL IMMUNE-GATED EOS CALL
  #
  # Molecular Tier is retained for EVERY CELL.
  #
  # EosRef_call is immune-gated.
  # ==========================================================

  tier_to_call <- c(

    "Conf Tier1" =
      "REF_EOS_TIER1",

    "Conf Tier2" =
      "REF_EOS_TIER2",

    "Conf Tier3" =
      "REF_EOS_TIER3",

    "Conf Tier4" =
      "REF_EOS_TIER4",

    "Rest" =
      "REF_EOS_REST"
  )


  eos_call <- rep(
    "OUTSIDE_IMMUNE",
    length(query_cells)
  )


  eos_call[
    query_is_immune
  ] <- unname(
    tier_to_call[
      query_rule$EosRef_tier[
        query_is_immune
      ]
    ]
  )


  # ==========================================================
  # 17. WANG EOS EVIDENCE
  # ==========================================================

  wang_predicted_eos <- rep(
    NA,
    length(query_cells)
  )

  wang_eos_score <- rep(
    NA_real_,
    length(query_cells)
  )


  if (
    wang_predicted_col %in%
    colnames(query_md)
  ) {

    wang_predicted_eos <- (
      as.character(
        query_md[[wang_predicted_col]]
      ) ==
        eos_label
    )
  }


  if (
    wang_eos_score_col %in%
    colnames(query_md)
  ) {

    wang_eos_score <- as.numeric(
      query_md[[wang_eos_score_col]]
    )
  }


  # ==========================================================
  # 18. NONIMMUNE REVIEW RESCUE
  #
  # Only:
  #
  # Tier1
  # +
  # Wang predicts Eosinophil
  #
  # Never automatically relabel.
  # ==========================================================

  if (
    allow_nonimmune_rescue &&
      any(
        !is.na(
          wang_predicted_eos
        )
      )
  ) {

    rescue <- (
      !query_is_immune &
        query_rule$EosRef_tier ==
        "Conf Tier1" &
        wang_predicted_eos %in% TRUE
    )

    eos_call[
      rescue
    ] <-
      "REVIEW_NONIMMUNE_EOS_RESCUE"
  }


  # ==========================================================
  # 19. CORE SCORE
  #
  # Diagnostic ranking only.
  # ==========================================================

  query_core_score <- as.numeric(

    query_binary[
      ,
      names(core_weights),
      drop = FALSE
    ] %*%
      core_weights
  )


  ref_core_score <- as.numeric(

    ref_binary[
      ,
      names(core_weights),
      drop = FALSE
    ] %*%
      core_weights
  )


  # ==========================================================
  # 20. RANK QUERY IMMUNE CELLS
  # ==========================================================

  rank_priority <- c(

    "REF_EOS_TIER1" = 1,

    "REF_EOS_TIER2" = 2,

    "REF_EOS_TIER3" = 3,

    "REF_EOS_TIER4" = 4,

    "REF_EOS_REST" = 5
  )


  eos_rank <- rep(
    NA_integer_,
    length(query_cells)
  )


  immune_idx <- which(
    query_is_immune
  )


  if (length(immune_idx) > 0) {

    priority_value <- unname(
      rank_priority[
        eos_call[
          immune_idx
        ]
      ]
    )


    rank_order <- order(

      priority_value,

      -query_core_score[
        immune_idx
      ],

      -query_rule$EosRef_n_all_markers_detected[
        immune_idx
      ]
    )


    eos_rank[
      immune_idx[
        rank_order
      ]
    ] <- seq_along(
      rank_order
    )
  }


  # ==========================================================
  # 21. QUERY METADATA
  # ==========================================================

  new_md <- data.frame(

    row.names =
      query_cells,


    EosRef_is_immune =
      query_is_immune,


    # Molecular Tier assigned regardless of broad annotation.
    EosRef_tier =
      query_rule$EosRef_tier,


    EosRef_tier_reason =
      query_rule$EosRef_tier_reason,


    # Immune-gated final call.
    EosRef_call =
      eos_call,


    EosRef_core_score =
      query_core_score,


    EosRef_rank =
      eos_rank,


    EosRef_n_core_detected =
      query_rule$EosRef_n_core_detected,


    EosRef_core_pattern =
      query_rule$EosRef_core_pattern,


    EosRef_n_tier2_support_detected =
      query_rule$EosRef_n_tier2_support_detected,


    EosRef_tier3_marker_detected =
      query_rule$EosRef_tier3_marker_detected,


    EosRef_n_all_markers_detected =
      query_rule$EosRef_n_all_markers_detected,


    EosRef_marker_combination =
      query_rule$EosRef_marker_combination,


    EosRef_Wang_predicted_Eos =
      wang_predicted_eos,


    EosRef_Wang_Eos_score =
      wang_eos_score,


    stringsAsFactors = FALSE
  )


  # ==========================================================
  # 22. PER-GENE DETECTION FLAGS
  # ==========================================================

  for (g in all_markers_use) {

    new_md[[
      paste0(
        "EosRef_detect_",
        g
      )
    ]] <- query_binary[
      ,
      g
    ]
  }


  # ==========================================================
  # 23. ADD TO SEURAT OBJECT
  # ==========================================================

  query <- AddMetaData(
    query,
    metadata = new_md
  )


  # ==========================================================
  # 24. REFERENCE CELL TABLE
  # ==========================================================

  reference_score_table <- data.frame(

    cell =
      ref_cells,

    CellType =
      ref_celltype,

    IsEosinophil =
      ref_is_eos,

    IsImmuneReference =
      ref_is_immune,

    EosRef_tier =
      reference_rule_all$EosRef_tier,

    EosRef_tier_reason =
      reference_rule_all$EosRef_tier_reason,

    EosRef_core_score =
      ref_core_score,

    EosRef_n_all_markers_detected =
      reference_rule_all$EosRef_n_all_markers_detected,

    EosRef_marker_combination =
      reference_rule_all$EosRef_marker_combination,

    stringsAsFactors = FALSE
  )


  # ==========================================================
  # 25. QUERY CANDIDATE TABLE
  # ==========================================================

  candidate_table <- cbind(

    data.frame(

      cell =
        query_cells,

      Final_CellType_main =
        as.character(
          query_md[[query_main_col]]
        ),

      stringsAsFactors = FALSE
    ),

    new_md
  )


  candidate_table <- candidate_table[
    order(
      candidate_table$EosRef_rank,
      na.last = TRUE
    ),
    ,
    drop = FALSE
  ]


  # ==========================================================
  # 26. QUERY SUMMARY
  # ==========================================================

  summary_table <- as.data.frame(

    table(
      new_md$EosRef_call,
      useNA = "ifany"
    ),

    stringsAsFactors = FALSE
  )

  colnames(
    summary_table
  ) <- c(
    "EosRef_call",
    "n_cells"
  )


  summary_table$fraction <-
    summary_table$n_cells /
    sum(summary_table$n_cells)


  tier_summary <- as.data.frame(

    table(
      new_md$EosRef_tier,
      useNA = "ifany"
    ),

    stringsAsFactors = FALSE
  )

  colnames(
    tier_summary
  ) <- c(
    "EosRef_tier",
    "n_cells"
  )


  tier_summary$fraction <-
    tier_summary$n_cells /
    sum(tier_summary$n_cells)


  # ==========================================================
  # 27. VERBOSE REPORT
  # ==========================================================

  if (verbose) {

    message(
      "Reference Eosinophils: ",
      n_ref_eos
    )

    message(
      "Reference other immune cells: ",
      n_ref_other
    )

    message(
      "Tier hierarchy:"
    )

    message(
      "  Tier1 = Siglecf + Ccr3"
    )

    message(
      "  Tier2 = core + Il5ra/Alox15"
    )

    message(
      "  Tier3 = core + Itgam"
    )

    message(
      "  Tier4 = any remaining >=2-marker combination"
    )

    message(
      "  Rest  = all remaining cells"
    )

    message(
      "Core weights (ranking only): ",
      paste(
        paste0(
          names(core_weights),
          "=",
          round(
            core_weights,
            3
          )
        ),
        collapse = "; "
      )
    )

    message(
      "Tier assignment complete."
    )
  }


  # ==========================================================
  # 28. RETURN
  # ==========================================================

  list(

    object =
      query,


    marker_stats =
      marker_stats,


    # Tier2/Tier3 combinations evaluated in reference.
    # Diagnostic only; they do not determine Tier assignment.
    pair_stats =
      pair_stats,


    # Cumulative sensitivity / FPR:
    # T1, T1-2, T1-3, T1-4.
    rule_performance =
      reference_rule_performance,


    calibration =
      reference_rule_performance,


    # All reference cells.
    reference_scores =
      reference_score_table,


    # Ready for stacked-bar plotting.
    reference_tier_distribution =
      reference_tier_distribution,


    candidate_table =
      candidate_table,


    summary =
      summary_table,


    tier_summary =
      tier_summary,


    parameters = list(

      core_markers =
        core_markers_use,

      tier2_support_markers =
        support_markers_use,

      tier3_marker =
        tier3_marker,

      tier4_context_markers =
        context_markers_use,

      all_markers =
        all_markers_use,

      reference_group_col =
        reference_group_col,

      reference_sample_col =
        reference_sample_col,

      tier_rule =
        "FIXED_HIERARCHICAL",

      min_competitor_cells =
        min_competitor_cells
    )
  )
}

#' Plot eos with celltypes.
#'
#' @param object Required `object` input; validated before computation.
#' @param celltypes Required `celltypes` input; validated before computation.
#' @param subtype_col Optional `subtype_col` input with the default shown in the function signature.
#' @param eos_label Optional `eos_label` input with the default shown in the function signature.
#' @param fov Optional `fov` input with the default shown in the function signature.
#' @param eos_col Optional `eos_col` input with the default shown in the function signature.
#' @param other_cols Optional `other_cols` input with the default shown in the function signature.
#' @param background_col Optional `background_col` input with the default shown in the function signature.
#' @param background_size Optional `background_size` input with the default shown in the function signature.
#' @param other_size Optional `other_size` input with the default shown in the function signature.
#' @param eos_size Optional `eos_size` input with the default shown in the function signature.
#' @param background_alpha Optional `background_alpha` input with the default shown in the function signature.
#' @param highlight_alpha Optional `highlight_alpha` input with the default shown in the function signature.
#' @param flip_xy Optional `flip_xy` input with the default shown in the function signature.
#' @return A `ggplot` object or named list of plots; input objects are not modified.
plot_eos_with_celltypes <- function(
  object,
  celltypes,
  subtype_col = "Final_CellType_subtype_refined",
  eos_label = "Eosinophil",
  fov = "fov",

  # Colours
  eos_col = "#E89A8F",
  other_cols = NULL,
  background_col = "#D9D9D9",

  # Point sizes
  background_size = 0.10,
  other_size = 0.65,
  eos_size = 1.00,

  background_alpha = 0.6,
  highlight_alpha = 0.95,

  flip_xy = FALSE
) {

  stopifnot(
    inherits(object, "Seurat"),
    subtype_col %in% colnames(object@meta.data),
    fov %in% names(object@images)
  )

  celltypes <- unique(as.character(celltypes))

  # Don't duplicate Eosinophil if accidentally supplied
  celltypes <- setdiff(
    celltypes,
    eos_label
  )

  # ----------------------------------------------------------
  # 1. Get centroid coordinates
  # ----------------------------------------------------------

  coords <- Seurat::GetTissueCoordinates(
    object[[fov]],
    which = "centroids"
  ) |>
    as.data.frame()

  # Usually rownames are cell IDs
  if (!"cell" %in% colnames(coords)) {
    coords$cell <- rownames(coords)
  }

  # ----------------------------------------------------------
  # 2. Match annotations
  # ----------------------------------------------------------

  md <- object@meta.data

  coords$celltype <- as.character(
    md[
      match(
        coords$cell,
        rownames(md)
      ),
      subtype_col
    ]
  )

  if (anyNA(coords$celltype)) {
    warning(
      sum(is.na(coords$celltype)),
      " centroid(s) could not be matched to ",
      subtype_col
    )
  }

  # ----------------------------------------------------------
  # 3. Define display groups
  # ----------------------------------------------------------

  coords$plot_group <- dplyr::case_when(
    coords$celltype == eos_label ~ eos_label,
    coords$celltype %in% celltypes ~ coords$celltype,
    TRUE ~ "Other"
  )

  # ----------------------------------------------------------
  # 4. Default comparison colours
  # ----------------------------------------------------------

  default_cols <- c(
    "#7E9AD9",
    "#BFDCC8",
    "#B9B3D7",
    "#F4D6A0",
    "#9FC5C6",
    "#C8B6A6",
    "#A9BEDC",
    "#D6B5CF"
  )

  if (is.null(other_cols)) {

    if (length(celltypes) > length(default_cols)) {
      stop(
        "More comparison cell types requested than default colours. ",
        "Please supply named 'other_cols'."
      )
    }

    other_cols <- default_cols[
      seq_along(celltypes)
    ]

    names(other_cols) <- celltypes

  } else {

    # If unnamed colours supplied, assign in celltypes order
    if (is.null(names(other_cols))) {

      if (length(other_cols) < length(celltypes)) {
        stop(
          "'other_cols' must contain at least one colour ",
          "for each requested cell type."
        )
      }

      other_cols <- other_cols[
        seq_along(celltypes)
      ]

      names(other_cols) <- celltypes
    }
  }

  # ----------------------------------------------------------
  # 5. Base ImageDimPlot
  # ----------------------------------------------------------

  display_cols <- c(
    "Other" = background_col,
    other_cols,
    setNames(eos_col, eos_label)
  )

  object$.__Eos_overlay_group__ <- coords$plot_group[
    match(
      colnames(object),
      coords$cell
    )
  ]

  p <- Seurat::ImageDimPlot(
    object,
    fov = fov,
    group.by = ".__Eos_overlay_group__",
    cols = display_cols,
    size = background_size,
    border.color = NA,
    dark.background = FALSE,
    flip_xy = flip_xy
  )

  # ----------------------------------------------------------
  # 6. Overlay comparison populations with larger points
  # ----------------------------------------------------------

  for (ct in celltypes) {

    tmp <- coords |>
      dplyr::filter(
        plot_group == ct
      )

    p <- p +
      ggplot2::geom_point(
        data = tmp,
        ggplot2::aes(
          x = x,
          y = y
        ),
        inherit.aes = FALSE,
        colour = other_cols[[ct]],
        size = other_size,
        alpha = highlight_alpha
      )
  }

  # ----------------------------------------------------------
  # 7. Overlay Eos last so they remain visible
  # ----------------------------------------------------------

  eos_df <- coords |>
    dplyr::filter(
      plot_group == eos_label
    )

  p <- p +
    ggplot2::geom_point(
      data = eos_df,
      ggplot2::aes(
        x = x,
        y = y
      ),
      inherit.aes = FALSE,
      colour = eos_col,
      size = eos_size,
      alpha = highlight_alpha
    )

  p
}

# -----------------------------------------------------------------------------
# Downstream tissue-branch, stability, and spatial-safety helpers
# -----------------------------------------------------------------------------

#' Recreate the approved primary downstream cell mask from aligned metadata
#'
#' @param cell_metadata Data frame already aligned to the Seurat object by cell
#'   ID and containing core-QC, high-control, and segmentation-multiplet flags.
#' @param core_col,high_control_col,multiplet_col Metadata column names.
#'
#' @return A logical vector in `cell_metadata` row order. A cell passes only
#'   when core QC is TRUE and both exclusion flags are FALSE.
derive_primary_include_revised <- function(
    cell_metadata,
    core_col = "qc_core_pass",
    high_control_col = "high_control_flag",
    multiplet_col = "segmentation_multiplet_flag"
) {
  required <- c(core_col, high_control_col, multiplet_col)
  missing <- setdiff(required, colnames(cell_metadata))
  if (length(missing)) {
    stop("Missing primary-mask fields: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  values <- lapply(cell_metadata[required], function(x) {
    if (!is.logical(x) || anyNA(x)) stop("Primary-mask fields must be complete logical vectors.", call. = FALSE)
    x
  })
  values[[1L]] & !values[[2L]] & !values[[3L]]
}

#' Build mutually exclusive all-cell, adipose, and lymph-node branch masks
#'
#' @param cell_metadata Data frame with unique cell IDs and an approved primary
#'   inclusion field.
#' @param lymph_node_cell_ids Character vector of QC-passed cells belonging to
#'   the frozen lymph-node tissue domain. An empty vector records no LN domain.
#' @param cell_id_col Name of the explicit cell-ID column.
#' @param primary_col Name of the approved QC-passed inclusion column.
#' @param provenance Character description of the domain source and parameters.
#'
#' @return A data frame in input row order containing the original cell ID,
#'   primary mask, frozen lymph-node mask, three branch masks, tissue-domain
#'   label, and provenance. Adipose and LN masks exactly partition QC-passed
#'   cells; failed-QC cells belong to neither child branch.
build_tissue_branch_manifest <- function(
    cell_metadata,
    lymph_node_cell_ids = character(),
    cell_id_col = "cell_id",
    primary_col = "primary_include_revised",
    provenance = NA_character_
) {
  if (!is.data.frame(cell_metadata)) stop("cell_metadata must be a data frame.", call. = FALSE)
  required <- c(cell_id_col, primary_col)
  missing <- setdiff(required, colnames(cell_metadata))
  if (length(missing)) stop("Missing branch-manifest fields: ", paste(missing, collapse = ", "), call. = FALSE)
  ids <- as.character(cell_metadata[[cell_id_col]])
  if (anyNA(ids) || any(!nzchar(ids))) stop("Cell IDs must be complete and non-empty.", call. = FALSE)
  if (anyDuplicated(ids)) stop("Duplicate cell IDs are not allowed.", call. = FALSE)
  primary <- cell_metadata[[primary_col]]
  if (!is.logical(primary) || anyNA(primary)) stop("The primary inclusion field must be complete and logical.", call. = FALSE)
  lymph_node_cell_ids <- unique(as.character(lymph_node_cell_ids))
  unknown <- setdiff(lymph_node_cell_ids, ids)
  if (length(unknown)) stop("Lymph-node IDs are not present in cell metadata: ", paste(head(unknown, 10L), collapse = ", "), call. = FALSE)
  failed_qc_ln <- lymph_node_cell_ids[!primary[match(lymph_node_cell_ids, ids)]]
  if (length(failed_qc_ln)) stop("Lymph-node IDs must be restricted to primary QC-passed cells.", call. = FALSE)
  lymph_node <- primary & ids %in% lymph_node_cell_ids
  adipose <- primary & !lymph_node
  if (!all((adipose | lymph_node) == primary) || any(adipose & lymph_node)) {
    stop("Adipose and lymph-node masks do not partition the primary cohort.", call. = FALSE)
  }
  data.frame(
    cell_id = ids,
    primary_include_revised = primary,
    lymph_node_include = lymph_node,
    all_qcpass_include = primary,
    adipose_only_include = adipose,
    lymph_node_only_include = lymph_node,
    tissue_domain = ifelse(!primary, "QC_EXCLUDED", ifelse(lymph_node, "LYMPH_NODE", "ADIPOSE")),
    branch_manifest_provenance = rep(as.character(provenance)[1L], length(ids)),
    stringsAsFactors = FALSE
  )
}

#' Select the ordered cell IDs for one frozen tissue branch
#'
#' @param branch_manifest Output from `build_tissue_branch_manifest()`.
#' @param branch One of `all_qcpass`, `adipose_only`, or `lymph_node_only`.
#' @param min_cells Minimum required number of cells; the function stops below
#'   this gate instead of forcing an underpowered analysis.
#'
#' @return Character vector of selected cell IDs in manifest order.
select_tissue_branch_ids <- function(branch_manifest, branch, min_cells = 100L) {
  branch <- match.arg(branch, c("all_qcpass", "adipose_only", "lymph_node_only"))
  column <- paste0(branch, "_include")
  required <- c("cell_id", column)
  missing <- setdiff(required, colnames(branch_manifest))
  if (length(missing)) stop("Branch manifest is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  include <- branch_manifest[[column]]
  if (!is.logical(include) || anyNA(include)) stop("Branch inclusion field must be complete and logical.", call. = FALSE)
  ids <- as.character(branch_manifest$cell_id[include])
  min_cells <- as.integer(min_cells)
  if (!is.finite(min_cells) || min_cells < 1L) stop("min_cells must be a positive integer.", call. = FALSE)
  if (length(ids) < min_cells) {
    stop("INSUFFICIENT_", toupper(branch), "_CELLS: ", length(ids), " < ", min_cells, call. = FALSE)
  }
  ids
}

#' Validate that spatial query and reference pools are truly disjoint
#'
#' @param query_ids,reference_ids Character cell-ID vectors.
#' @param query_coordinates,reference_coordinates Optional data frames with
#'   `cell_id`, `x`, and `y`; when supplied, duplicate cross-pool coordinates
#'   are prohibited because they generate zero-distance matches.
#'
#' @return A one-row list with query/reference counts, ID-overlap count, and
#'   duplicate-coordinate-pair count. The function stops on any unsafe overlap.
validate_disjoint_spatial_pools <- function(
    query_ids,
    reference_ids,
    query_coordinates = NULL,
    reference_coordinates = NULL
) {
  query_ids <- as.character(query_ids)
  reference_ids <- as.character(reference_ids)
  if (!length(query_ids) || !length(reference_ids)) stop("Spatial pools must both be non-empty.", call. = FALSE)
  if (anyNA(query_ids) || anyNA(reference_ids) || any(!nzchar(query_ids)) || any(!nzchar(reference_ids))) {
    stop("Spatial pool IDs must be complete and non-empty.", call. = FALSE)
  }
  if (anyDuplicated(query_ids) || anyDuplicated(reference_ids)) stop("Duplicate IDs within a spatial pool are not allowed.", call. = FALSE)
  overlap <- intersect(query_ids, reference_ids)
  if (length(overlap)) stop("Spatial query/reference cell-ID overlap detected: ", paste(head(overlap, 10L), collapse = ", "), call. = FALSE)
  duplicate_coordinates <- 0L
  if (xor(is.null(query_coordinates), is.null(reference_coordinates))) {
    stop("Supply both coordinate tables or neither.", call. = FALSE)
  }
  if (!is.null(query_coordinates)) {
    validate_coordinates <- function(x, ids, label) {
      required <- c("cell_id", "x", "y")
      missing <- setdiff(required, colnames(x))
      if (length(missing)) stop(label, " coordinates are missing: ", paste(missing, collapse = ", "), call. = FALSE)
      if (anyDuplicated(x$cell_id)) stop(label, " coordinates contain duplicate cell IDs.", call. = FALSE)
      index <- match(ids, as.character(x$cell_id))
      if (anyNA(index)) stop(label, " coordinates do not cover every pool cell.", call. = FALSE)
      out <- x[index, required, drop = FALSE]
      if (any(!is.finite(out$x)) || any(!is.finite(out$y))) stop(label, " coordinates must be finite.", call. = FALSE)
      out
    }
    q <- validate_coordinates(query_coordinates, query_ids, "Query")
    r <- validate_coordinates(reference_coordinates, reference_ids, "Reference")
    q_key <- paste(format(q$x, digits = 17), format(q$y, digits = 17), sep = "\r")
    r_key <- paste(format(r$x, digits = 17), format(r$y, digits = 17), sep = "\r")
    duplicate_coordinates <- sum(q_key %in% r_key)
    if (duplicate_coordinates) stop("Cross-pool duplicate coordinate pairs would create zero-distance matches.", call. = FALSE)
  }
  list(
    n_query = length(query_ids),
    n_reference = length(reference_ids),
    n_overlap_ids = length(overlap),
    n_duplicate_coordinate_pairs = as.integer(duplicate_coordinates)
  )
}

#' Enforce the positive-distance gate after nearest-neighbour matching
#'
#' @param distances Numeric vector or matrix of spatial distances.
#' @param zero_tolerance Values less than or equal to this tolerance are unsafe.
#'
#' @return A one-row list containing count, minimum, median, maximum, and the
#'   number of zero/negative distances. Stops if an unsafe distance is present.
validate_neighbour_distances <- function(distances, zero_tolerance = 0) {
  distances <- as.numeric(distances)
  if (!length(distances) || any(!is.finite(distances))) stop("Neighbour distances must be non-empty and finite.", call. = FALSE)
  unsafe <- distances <= zero_tolerance
  if (any(unsafe)) stop("Zero or negative nearest-neighbour distance detected; spatial result is invalid.", call. = FALSE)
  list(
    n_distances = length(distances),
    minimum_distance = min(distances),
    median_distance = stats::median(distances),
    maximum_distance = max(distances),
    n_zero_or_negative = sum(unsafe)
  )
}

#' Build cell-ID-aligned Eosinophil and non-Eosinophil spatial pools
#'
#' @param cell_metadata Data frame with one row per cell and a unique `cell_id`.
#' @param coordinates Data frame with unique `cell_id`, `x`, and `y` columns.
#' @param eos_col Complete logical field defining the Eosinophil query pool.
#' @param cell_type_col Downstream cell-type field for reference labels.
#' @param state_col Continuous Eosinophil-state field.
#' @param extreme_col Optional descriptive Eosinophil-tail field.
#'
#' @return A typed list with aligned `all`, `query`, and `reference` tables plus
#'   the disjoint-pool validation result. An empty query or reference pool is a
#'   non-fatal `SKIPPED_EMPTY_EOS_OR_REFERENCE_POOL` result so small tissue
#'   branches can complete their audit outputs.
build_eos_spatial_pools <- function(
    cell_metadata,
    coordinates,
    eos_col = "Eos_inclusive",
    cell_type_col = "Final_CellType_subtype",
    state_col = "EosState_balance",
    extreme_col = "EosState_extreme"
) {
  metadata <- as.data.frame(cell_metadata, stringsAsFactors = FALSE)
  coords <- as.data.frame(coordinates, stringsAsFactors = FALSE)
  metadata_required <- c("cell_id", eos_col, cell_type_col)
  coordinate_required <- c("cell_id", "x", "y")
  missing_metadata <- setdiff(metadata_required, names(metadata))
  missing_coordinates <- setdiff(coordinate_required, names(coords))
  if (length(missing_metadata)) {
    stop("Cell metadata missing: ", paste(missing_metadata, collapse = ", "), call. = FALSE)
  }
  if (length(missing_coordinates)) {
    stop("Coordinates missing: ", paste(missing_coordinates, collapse = ", "), call. = FALSE)
  }
  metadata$cell_id <- as.character(metadata$cell_id)
  coords$cell_id <- as.character(coords$cell_id)
  if (anyDuplicated(metadata$cell_id) || anyDuplicated(coords$cell_id)) {
    stop("Cell metadata and coordinates require unique cell IDs.", call. = FALSE)
  }
  index <- match(metadata$cell_id, coords$cell_id)
  if (anyNA(index)) {
    stop("Coordinates do not cover every metadata cell ID.", call. = FALSE)
  }
  eos <- metadata[[eos_col]]
  if (!is.logical(eos) || anyNA(eos)) {
    stop(eos_col, " must be a complete logical field.", call. = FALSE)
  }
  cell_type <- as.character(metadata[[cell_type_col]])
  if (anyNA(cell_type) || any(!nzchar(cell_type))) {
    stop(cell_type_col, " must contain complete downstream labels.", call. = FALSE)
  }
  state <- if (state_col %in% names(metadata)) as.numeric(metadata[[state_col]]) else rep(NA_real_, nrow(metadata))
  extreme <- if (extreme_col %in% names(metadata)) as.character(metadata[[extreme_col]]) else rep(NA_character_, nrow(metadata))
  all_cells <- data.frame(
    cell_id = metadata$cell_id,
    x = as.numeric(coords$x[index]),
    y = as.numeric(coords$y[index]),
    Eos_inclusive = eos,
    cell_type = cell_type,
    EosState_balance = state,
    EosState_extreme = extreme,
    stringsAsFactors = FALSE
  )
  if (any(!is.finite(all_cells$x)) || any(!is.finite(all_cells$y))) {
    stop("All spatial coordinates must be finite.", call. = FALSE)
  }
  query <- all_cells[all_cells$Eos_inclusive, , drop = FALSE]
  reference <- all_cells[!all_cells$Eos_inclusive, , drop = FALSE]
  if (!nrow(query) || !nrow(reference)) {
    return(list(
      status = "SKIPPED_EMPTY_EOS_OR_REFERENCE_POOL",
      message = "Spatial query and reference pools must both be non-empty.",
      all = all_cells, query = query, reference = reference,
      gate = list(
        status = "SKIPPED_EMPTY_EOS_OR_REFERENCE_POOL",
        n_query = nrow(query), n_reference = nrow(reference),
        n_overlap_ids = 0L, n_duplicate_coordinate_pairs = 0L
      )
    ))
  }
  gate <- validate_disjoint_spatial_pools(
    query$cell_id, reference$cell_id,
    query[, c("cell_id", "x", "y"), drop = FALSE],
    reference[, c("cell_id", "x", "y"), drop = FALSE]
  )
  list(
    status = "PASS",
    message = "Spatial Eosinophil and reference pools are non-empty and disjoint.",
    all = all_cells, query = query, reference = reference, gate = gate
  )
}

#' Calculate Eosinophil-to-reference K-nearest-neighbour edge tables
#'
#' @param pools Output from `build_eos_spatial_pools()`.
#' @param k_values Positive requested neighbour counts.
#'
#' @return A named list (`k1`, `k15`, and so on) of edge data frames. Requested
#'   k is capped at the available reference-cell count and recorded as
#'   `k_actual`.
calculate_eos_knn_edges <- function(pools, k_values = c(1L, 15L)) {
  require_package("FNN")
  if (!is.list(pools) || !all(c("query", "reference", "gate") %in% names(pools))) {
    stop("pools must be returned by build_eos_spatial_pools().", call. = FALSE)
  }
  query <- pools$query
  reference <- pools$reference
  if (!nrow(query) || !nrow(reference)) {
    stop("Spatial query and reference pools must both be non-empty.", call. = FALSE)
  }
  k_values <- sort(unique(as.integer(k_values)))
  if (!length(k_values) || any(!is.finite(k_values)) || any(k_values < 1L)) {
    stop("k_values must contain positive integers.", call. = FALSE)
  }
  result <- lapply(k_values, function(k_requested) {
    k_actual <- min(k_requested, nrow(reference))
    fit <- FNN::get.knnx(
      data = as.matrix(reference[, c("x", "y"), drop = FALSE]),
      query = as.matrix(query[, c("x", "y"), drop = FALSE]),
      k = k_actual
    )
    distance <- as.vector(t(fit$nn.dist))
    validate_neighbour_distances(distance)
    reference_index <- as.vector(t(fit$nn.index))
    data.frame(
      eos_cell_id = rep(query$cell_id, each = k_actual),
      reference_cell_id = reference$cell_id[reference_index],
      neighbour_rank = rep(seq_len(k_actual), times = nrow(query)),
      k_requested = k_requested,
      k_actual = k_actual,
      distance = distance,
      reference_cell_type = reference$cell_type[reference_index],
      EosState_balance = rep(query$EosState_balance, each = k_actual),
      EosState_extreme = rep(query$EosState_extreme, each = k_actual),
      stringsAsFactors = FALSE
    )
  })
  names(result) <- paste0("k", k_values)
  result
}

#' Summarize Eosinophil KNN composition overall and by descriptive state
#'
#' @param edges One edge table returned by `calculate_eos_knn_edges()`.
#'
#' @return A list with `overall` and `by_state` count/fraction tables.
summarise_eos_knn_composition <- function(edges) {
  required <- c("reference_cell_type", "EosState_extreme")
  missing <- setdiff(required, names(edges))
  if (length(missing)) stop("KNN edges missing: ", paste(missing, collapse = ", "), call. = FALSE)
  overall <- as.data.frame(table(reference_cell_type = as.character(edges$reference_cell_type)), stringsAsFactors = FALSE)
  names(overall)[names(overall) == "Freq"] <- "n_edges"
  overall <- overall[overall$n_edges > 0L, , drop = FALSE]
  overall$fraction <- overall$n_edges / sum(overall$n_edges)
  state_keep <- !is.na(edges$EosState_extreme) & nzchar(as.character(edges$EosState_extreme))
  if (any(state_keep)) {
    by_state <- as.data.frame(table(
      EosState_extreme = as.character(edges$EosState_extreme[state_keep]),
      reference_cell_type = as.character(edges$reference_cell_type[state_keep])
    ), stringsAsFactors = FALSE)
    names(by_state)[names(by_state) == "Freq"] <- "n_edges"
    by_state <- by_state[by_state$n_edges > 0L, , drop = FALSE]
    totals <- ave(by_state$n_edges, by_state$EosState_extreme, FUN = sum)
    by_state$fraction <- by_state$n_edges / totals
  } else {
    by_state <- data.frame(
      EosState_extreme = character(), reference_cell_type = character(),
      n_edges = integer(), fraction = numeric(), stringsAsFactors = FALSE
    )
  }
  list(overall = overall, by_state = by_state)
}

#' Calculate each Eosinophil cell's distance to every eligible cell type
#'
#' @param pools Output from `build_eos_spatial_pools()`.
#' @param min_reference_cells Minimum reference cells required for a type.
#'
#' @return A list with cell-level distances and cell-type summaries containing
#'   median, quartiles and descriptive continuous-state Spearman correlation.
calculate_eos_distance_by_cell_type <- function(pools, min_reference_cells = 20L) {
  require_package("FNN")
  query <- pools$query
  reference <- pools$reference
  min_reference_cells <- as.integer(min_reference_cells)
  type_counts <- table(reference$cell_type)
  eligible <- names(type_counts[type_counts >= min_reference_cells])
  eligible <- scwat_cell_type_order(eligible)
  if (!length(eligible)) {
    return(list(cell_level = data.frame(), summary = data.frame()))
  }
  cell_level <- do.call(rbind, lapply(eligible, function(cell_type_name) {
    type_reference <- reference[reference$cell_type == cell_type_name, , drop = FALSE]
    fit <- FNN::get.knnx(
      as.matrix(type_reference[, c("x", "y"), drop = FALSE]),
      as.matrix(query[, c("x", "y"), drop = FALSE]),
      k = 1L
    )
    distance <- as.numeric(fit$nn.dist[, 1L])
    validate_neighbour_distances(distance)
    data.frame(
      eos_cell_id = query$cell_id,
      reference_cell_type = cell_type_name,
      distance = distance,
      EosState_balance = query$EosState_balance,
      stringsAsFactors = FALSE
    )
  }))
  rownames(cell_level) <- NULL
  split_distance <- split(cell_level, cell_level$reference_cell_type)
  summary <- do.call(rbind, lapply(names(split_distance), function(cell_type_name) {
    x <- split_distance[[cell_type_name]]
    complete <- is.finite(x$distance) & is.finite(x$EosState_balance)
    rho <- if (sum(complete) >= 3L && stats::sd(x$distance[complete]) > 0 && stats::sd(x$EosState_balance[complete]) > 0) {
      stats::cor(x$distance[complete], x$EosState_balance[complete], method = "spearman")
    } else {
      NA_real_
    }
    data.frame(
      reference_cell_type = cell_type_name,
      n_eos = nrow(x),
      n_reference = unname(type_counts[[cell_type_name]]),
      median_distance = stats::median(x$distance),
      q25_distance = unname(stats::quantile(x$distance, 0.25)),
      q75_distance = unname(stats::quantile(x$distance, 0.75)),
      spearman_rho_state = rho,
      interpretation = "DESCRIPTIVE_SPATIALLY_AUTOCORRELATED_NO_CELL_LEVEL_P_VALUE",
      stringsAsFactors = FALSE
    )
  }))
  rownames(summary) <- NULL
  list(cell_level = cell_level, summary = summary)
}

#' Rank cell types associated with the continuous Eosinophil KNN state
#'
#' @param edges_k15 KNN edge table containing Eosinophil ID, reference type and
#'   continuous state.
#' @param biological_order Preferred cell-type tie-breaking order.
#' @param min_eos Minimum Eosinophils with finite state.
#' @param top_n Maximum selected types in each state direction.
#'
#' @return A list with the full association table, per-Eosinophil composition,
#'   and deterministic short- and long-associated top-type vectors.
rank_eos_state_knn_associations <- function(
    edges_k15,
    biological_order = scwat_cell_type_order(),
    min_eos = 20L,
    top_n = 3L
) {
  required <- c("eos_cell_id", "reference_cell_type", "EosState_balance")
  missing <- setdiff(required, names(edges_k15))
  if (length(missing)) stop("KNN edges missing: ", paste(missing, collapse = ", "), call. = FALSE)
  edges <- as.data.frame(edges_k15, stringsAsFactors = FALSE)
  state_by_eos <- tapply(edges$EosState_balance, edges$eos_cell_id, function(x) unique(x[is.finite(x)]))
  valid_state <- vapply(state_by_eos, length, integer(1)) == 1L
  eos_ids <- names(state_by_eos)[valid_state]
  if (length(eos_ids) < as.integer(min_eos)) {
    return(list(
      status = "SKIPPED_INSUFFICIENT_EOS",
      full = data.frame(), per_eos = data.frame(),
      short_top = character(), long_top = character()
    ))
  }
  types <- unique(as.character(edges$reference_cell_type))
  grid <- expand.grid(
    eos_cell_id = eos_ids,
    reference_cell_type = types,
    stringsAsFactors = FALSE
  )
  counts <- aggregate(
    rep(1L, nrow(edges[edges$eos_cell_id %in% eos_ids, , drop = FALSE])),
    by = list(
      eos_cell_id = edges$eos_cell_id[edges$eos_cell_id %in% eos_ids],
      reference_cell_type = edges$reference_cell_type[edges$eos_cell_id %in% eos_ids]
    ),
    FUN = sum
  )
  names(counts)[[3L]] <- "n_edges"
  per_eos <- merge(grid, counts, by = c("eos_cell_id", "reference_cell_type"), all.x = TRUE, sort = FALSE)
  per_eos$n_edges[is.na(per_eos$n_edges)] <- 0L
  total_by_eos <- table(edges$eos_cell_id[edges$eos_cell_id %in% eos_ids])
  per_eos$k_actual <- as.integer(total_by_eos[per_eos$eos_cell_id])
  per_eos$neighbour_fraction <- per_eos$n_edges / per_eos$k_actual
  per_eos$EosState_balance <- as.numeric(vapply(
    state_by_eos[per_eos$eos_cell_id], `[[`, numeric(1), 1L
  ))
  split_type <- split(per_eos, per_eos$reference_cell_type)
  full <- do.call(rbind, lapply(names(split_type), function(cell_type_name) {
    x <- split_type[[cell_type_name]]
    rho <- if (stats::sd(x$neighbour_fraction) > 0 && stats::sd(x$EosState_balance) > 0) {
      stats::cor(x$neighbour_fraction, x$EosState_balance, method = "spearman")
    } else {
      NA_real_
    }
    data.frame(
      reference_cell_type = cell_type_name,
      n_eos = nrow(x),
      n_edges = sum(x$n_edges),
      mean_neighbour_fraction = mean(x$neighbour_fraction),
      spearman_rho = rho,
      direction = if (is.na(rho)) "UNINFORMATIVE" else if (rho > 0) "LONG_ASSOCIATED" else if (rho < 0) "SHORT_ASSOCIATED" else "NEUTRAL",
      interpretation = "DESCRIPTIVE_SPATIALLY_AUTOCORRELATED_NO_CELL_LEVEL_P_VALUE",
      stringsAsFactors = FALSE
    )
  }))
  rownames(full) <- NULL
  order_index <- match(full$reference_cell_type, biological_order)
  order_index[is.na(order_index)] <- length(biological_order) + seq_len(sum(is.na(order_index)))
  short_rows <- which(is.finite(full$spearman_rho) & full$spearman_rho < 0)
  long_rows <- which(is.finite(full$spearman_rho) & full$spearman_rho > 0)
  short_rows <- short_rows[order(full$spearman_rho[short_rows], -full$n_edges[short_rows], order_index[short_rows])]
  long_rows <- long_rows[order(-full$spearman_rho[long_rows], -full$n_edges[long_rows], order_index[long_rows])]
  list(
    status = "PASS",
    full = full,
    per_eos = per_eos,
    short_top = head(full$reference_cell_type[short_rows], as.integer(top_n)),
    long_top = head(full$reference_cell_type[long_rows], as.integer(top_n))
  )
}

#' Derive adequately sized Eosinophil groups for group-based CellChat
#'
#' @param state Numeric continuous Eosinophil-state values.
#' @param lower_probability Lower quantile defining the short-enriched group.
#' @param upper_probability Upper quantile defining the long-enriched group.
#' @param min_cells Minimum cells required in each retained state group.
#'
#' @return A list with typed status, per-input group vector, quantile thresholds
#'   and group counts. The middle observations remain unassigned.
derive_eos_cellchat_groups <- function(
    state,
    lower_probability = 0.30,
    upper_probability = 0.70,
    min_cells = 10L
) {
  state <- as.numeric(state)
  if (!is.finite(lower_probability) || !is.finite(upper_probability) ||
      lower_probability <= 0 || upper_probability >= 1 ||
      lower_probability >= upper_probability) {
    stop("CellChat state probabilities must satisfy 0 < lower < upper < 1.", call. = FALSE)
  }
  group <- rep(NA_character_, length(state))
  finite <- is.finite(state)
  if (!any(finite)) {
    return(list(
      status = "SKIPPED_INSUFFICIENT_STATE_GROUP_CELLS",
      message = "No finite Eosinophil-state values are available.",
      group = group, lower_threshold = NA_real_, upper_threshold = NA_real_,
      counts = integer()
    ))
  }
  lower <- unname(stats::quantile(state[finite], lower_probability, type = 7))
  upper <- unname(stats::quantile(state[finite], upper_probability, type = 7))
  group[finite & state <= lower] <- "Eos_short_enriched"
  group[finite & state >= upper] <- "Eos_long_enriched"
  counts <- table(factor(
    group,
    levels = c("Eos_short_enriched", "Eos_long_enriched")
  ), useNA = "no")
  status <- if (all(counts >= as.integer(min_cells))) {
    "PASS"
  } else {
    "SKIPPED_INSUFFICIENT_STATE_GROUP_CELLS"
  }
  list(
    status = status,
    message = if (status == "PASS") {
      "Prespecified continuous-state tails meet the CellChat group-size gate."
    } else {
      sprintf(
        "CellChat requires at least %d cells in each Eosinophil state group.",
        as.integer(min_cells)
      )
    },
    group = group,
    lower_threshold = lower,
    upper_threshold = upper,
    counts = counts
  )
}

#' Prepare spatial CellChat inputs for Eosinophils and selected neighbours
#'
#' @param object Seurat object containing normalized Xenium expression and
#'   required metadata.
#' @param coordinates Data frame with unique `cell_id`, `x`, and `y`.
#' @param top_short,top_long Cell types selected from continuous-state KNN
#'   association in the short and long directions.
#' @param eos_col Complete logical Eosinophil-selection field.
#' @param cell_type_col Downstream non-Eosinophil grouping field.
#' @param state_col Continuous Eosinophil-state field.
#' @param cell_area_col Cell-area field in square microns.
#' @param assay Seurat assay holding non-negative normalized expression.
#' @param min_cells Minimum cells required per CellChat group.
#'
#' @return A typed list containing sparse normalized expression, aligned group
#'   metadata, aligned coordinates, spatial scale factors and group counts.
prepare_eos_cellchat_inputs <- function(
    object,
    coordinates,
    top_short,
    top_long,
    eos_col = "Eos_inclusive",
    cell_type_col = "Final_CellType_subtype",
    state_col = "EosState_balance",
    cell_area_col = "cell_area",
    assay = "Xenium",
    min_cells = 10L
) {
  if (!inherits(object, "Seurat")) stop("object must be a Seurat object.", call. = FALSE)
  required_metadata <- c(eos_col, cell_type_col, state_col, cell_area_col)
  missing <- setdiff(required_metadata, colnames(object@meta.data))
  if (length(missing)) stop("Seurat metadata missing: ", paste(missing, collapse = ", "), call. = FALSE)
  top_types <- scwat_cell_type_order(unique(c(as.character(top_short), as.character(top_long))))
  top_types <- top_types[!is.na(top_types) & nzchar(top_types)]
  if (!length(top_types)) {
    return(list(status = "SKIPPED_NO_TOP_NEIGHBOUR_TYPES", message = "No KNN-selected neighbour types."))
  }
  metadata <- object@meta.data
  eos <- as.logical(metadata[[eos_col]])
  if (anyNA(eos)) stop(eos_col, " must be complete.", call. = FALSE)
  state_groups <- derive_eos_cellchat_groups(metadata[[state_col]][eos], min_cells = min_cells)
  if (state_groups$status != "PASS") {
    return(c(state_groups[c("status", "message")], list(state_groups = state_groups)))
  }
  eos_ids <- rownames(metadata)[eos]
  eos_group <- stats::setNames(state_groups$group, eos_ids)
  eos_selected <- eos_ids[!is.na(eos_group)]
  non_eos_selected <- rownames(metadata)[
    !eos & as.character(metadata[[cell_type_col]]) %in% top_types
  ]
  selected <- c(eos_selected, non_eos_selected)
  group <- c(
    unname(eos_group[eos_selected]),
    as.character(metadata[non_eos_selected, cell_type_col])
  )
  names(group) <- selected
  group_counts <- table(group)
  if (!length(non_eos_selected) || any(group_counts < as.integer(min_cells))) {
    return(list(
      status = "SKIPPED_INSUFFICIENT_CELLCHAT_GROUP_CELLS",
      message = sprintf("Every selected group requires at least %d cells.", as.integer(min_cells)),
      group_counts = group_counts,
      state_groups = state_groups
    ))
  }
  coords <- as.data.frame(coordinates, stringsAsFactors = FALSE)
  if (!all(c("cell_id", "x", "y") %in% names(coords)) || anyDuplicated(coords$cell_id)) {
    stop("coordinates require unique cell_id, x and y columns.", call. = FALSE)
  }
  coordinate_index <- match(selected, as.character(coords$cell_id))
  if (anyNA(coordinate_index)) stop("CellChat coordinates do not cover every selected cell.", call. = FALSE)
  aligned_coordinates <- data.frame(
    x = as.numeric(coords$x[coordinate_index]),
    y = as.numeric(coords$y[coordinate_index]),
    row.names = selected,
    stringsAsFactors = FALSE
  )
  if (any(!is.finite(as.matrix(aligned_coordinates)))) {
    stop("CellChat coordinates must be finite.", call. = FALSE)
  }
  expression <- SeuratObject::GetAssayData(object, assay = assay, layer = "data")
  expression <- expression[, selected, drop = FALSE]
  if (!identical(colnames(expression), selected) || any(expression@x < 0)) {
    stop("CellChat requires aligned non-negative normalized expression.", call. = FALSE)
  }
  cell_area <- as.numeric(metadata[selected, cell_area_col])
  valid_area <- is.finite(cell_area) & cell_area > 0
  if (!any(valid_area)) stop("A positive cell area is required for the spatial scale.", call. = FALSE)
  equivalent_diameter <- 2 * sqrt(cell_area[valid_area] / pi)
  meta <- data.frame(
    cellchat_group = factor(group, levels = c("Eos_short_enriched", "Eos_long_enriched", top_types)),
    samples = factor(rep("sample1", length(selected))),
    original_cell_type = as.character(metadata[selected, cell_type_col]),
    EosState_balance = as.numeric(metadata[selected, state_col]),
    row.names = selected,
    stringsAsFactors = FALSE
  )
  list(
    status = "PASS",
    message = "Spatial CellChat inputs passed alignment and group-size gates.",
    data = expression,
    meta = meta,
    coordinates = as.matrix(aligned_coordinates),
    scale_factors = list(
      # Xenium centroid coordinates are already in micrometres. In the legacy
      # CellChat API, spot.diameter / spot is the coordinate-to-micron ratio,
      # so equal values encode ratio = 1.
      spot = stats::median(equivalent_diameter),
      spot.diameter = stats::median(equivalent_diameter)
    ),
    spatial_factors = data.frame(
      ratio = 1,
      tol = stats::median(equivalent_diameter) / 2,
      row.names = "sample1",
      stringsAsFactors = FALSE
    ),
    group_counts = as.data.frame(group_counts, stringsAsFactors = FALSE),
    state_groups = state_groups,
    top_short = intersect(top_types, as.character(top_short)),
    top_long = intersect(top_types, as.character(top_long))
  )
}

#' Create a spatial CellChat object across legacy and current APIs
#'
#' CellChat up to the archived API accepts `scale.factors`; current jinworks
#' CellChat accepts `spatial.factors`. This adapter inspects the callable's
#' formal arguments and supplies exactly one representation.
#'
#' @param inputs Passed result from `prepare_eos_cellchat_inputs()`.
#' @param create_fun CellChat object-construction function. Defaults to the
#'   installed package export and can be injected for compatibility tests.
#'
#' @return The object returned by `create_fun`.
create_cellchat_object_compatible <- function(inputs, create_fun = NULL) {
  if (!identical(inputs$status, "PASS")) {
    stop("CellChat inputs must pass before object creation.", call. = FALSE)
  }
  if (is.null(create_fun)) {
    require_package("CellChat")
    create_fun <- getExportedValue("CellChat", "createCellChat")
  }
  formal_names <- names(formals(create_fun))
  arguments <- list(
    object = inputs$data,
    meta = inputs$meta,
    group.by = "cellchat_group",
    datatype = "spatial",
    coordinates = inputs$coordinates
  )
  if ("spatial.factors" %in% formal_names) {
    arguments$spatial.factors <- inputs$spatial_factors
  } else if ("scale.factors" %in% formal_names) {
    arguments$scale.factors <- inputs$scale_factors
  } else {
    stop(
      "Unsupported CellChat createCellChat() API: neither spatial.factors nor scale.factors is available.",
      call. = FALSE
    )
  }
  do.call(create_fun, arguments)
}

#' Run the optional spatial CellChat pipeline with typed failure status
#'
#' @param inputs Output from `prepare_eos_cellchat_inputs()`.
#' @param database Optional CellChat mouse database override.
#' @param seed Random seed for CellChat probability calculations.
#' @param min_cells Minimum cells used by `filterCommunication()`.
#'
#' @return A list with status, message, package version, failed/completed stage,
#'   optional CellChat object and extracted communication table.
run_eos_spatial_cellchat <- function(
    inputs,
    database = NULL,
    seed = 1234L,
    min_cells = 10L
) {
  input_status <- as.character(inputs$status %||% "INVALID_INPUT")
  if (!identical(input_status, "PASS")) {
    return(list(
      status = input_status,
      message = as.character(inputs$message %||% "CellChat input preparation did not pass."),
      package_version = NA_character_, stage = "INPUT_GATE",
      object = NULL, communication = data.frame()
    ))
  }
  if (!requireNamespace("CellChat", quietly = TRUE)) {
    return(list(
      status = "SKIPPED_PACKAGE_UNAVAILABLE",
      message = "Optional package 'CellChat' is unavailable.",
      package_version = NA_character_, stage = "PACKAGE_GATE",
      object = NULL, communication = data.frame()
    ))
  }
  package_version <- as.character(utils::packageVersion("CellChat"))
  stage <- "CREATE_OBJECT"
  result <- tryCatch({
    set.seed(as.integer(seed))
    cellchat <- create_cellchat_object_compatible(inputs)
    stage <- "SET_DATABASE"
    cellchat@DB <- database %||% getExportedValue("CellChat", "CellChatDB.mouse")
    stage <- "SUBSET_DATA"
    cellchat <- CellChat::subsetData(cellchat)
    stage <- "OVEREXPRESSED_GENES"
    cellchat <- CellChat::identifyOverExpressedGenes(cellchat)
    stage <- "OVEREXPRESSED_INTERACTIONS"
    cellchat <- CellChat::identifyOverExpressedInteractions(cellchat)
    stage <- "COMMUNICATION_PROBABILITY"
    cellchat <- CellChat::computeCommunProb(
      cellchat,
      type = "triMean",
      distance.use = TRUE,
      raw.use = TRUE
    )
    stage <- "FILTER_COMMUNICATION"
    cellchat <- CellChat::filterCommunication(cellchat, min.cells = as.integer(min_cells))
    stage <- "PATHWAY_PROBABILITY"
    cellchat <- CellChat::computeCommunProbPathway(cellchat)
    stage <- "AGGREGATE_NETWORK"
    cellchat <- CellChat::aggregateNet(cellchat)
    stage <- "EXTRACT_COMMUNICATION"
    communication <- CellChat::subsetCommunication(cellchat)
    list(object = cellchat, communication = as.data.frame(communication, stringsAsFactors = FALSE))
  }, error = identity)
  if (inherits(result, "error")) {
    return(list(
      status = "FAILED_CELLCHAT_RUNTIME",
      message = conditionMessage(result), package_version = package_version,
      stage = stage, object = NULL, communication = data.frame()
    ))
  }
  list(
    status = "PASS",
    message = "Spatial CellChat completed; interpret within this section only.",
    package_version = package_version,
    stage = "COMPLETE",
    object = result$object,
    communication = result$communication
  )
}

#' Filter and annotate Eosinophil-focused CellChat interactions
#'
#' @param table Communication table returned by CellChat.
#' @param raw_p_max Maximum raw CellChat permutation p-value.
#' @param adjusted_p_max Maximum BH-adjusted p-value across reported rows.
#' @param top_short,top_long Prespecified KNN-selected neighbour types for the
#'   short-enriched and long-enriched Eosinophil groups. When supplied, each
#'   state is restricted to its matching neighbour set.
#'
#' @return A list containing the complete annotated table and the significant
#'   outgoing Eosinophil-to-neighbour subset.
filter_eos_cellchat_interactions <- function(
    table,
    raw_p_max = 0.05,
    adjusted_p_max = 0.10,
    top_short = NULL,
    top_long = NULL
) {
  communication <- as.data.frame(table, stringsAsFactors = FALSE)
  required <- c("source", "target", "interaction_name", "prob", "pval")
  missing <- setdiff(required, names(communication))
  if (length(missing)) stop("CellChat table missing: ", paste(missing, collapse = ", "), call. = FALSE)
  if (!nrow(communication)) {
    communication$p_adjust_bh <- numeric()
    communication$direction <- character()
    communication$eos_state_group <- character()
    communication$neighbour_cell_type <- character()
    communication$matches_state_top_neighbour <- logical()
    communication$analysis_scope <- character()
    return(list(all = communication, significant = communication))
  }
  communication$p_adjust_bh <- stats::p.adjust(communication$pval, method = "BH")
  source_eos <- startsWith(as.character(communication$source), "Eos_")
  target_eos <- startsWith(as.character(communication$target), "Eos_")
  communication$direction <- ifelse(
    source_eos & !target_eos, "EOS_TO_NEIGHBOUR",
    ifelse(!source_eos & target_eos, "NEIGHBOUR_TO_EOS", "OUTSIDE_REQUESTED_DIRECTION")
  )
  communication$eos_state_group <- ifelse(
    source_eos, as.character(communication$source),
    ifelse(target_eos, as.character(communication$target), NA_character_)
  )
  communication$neighbour_cell_type <- ifelse(
    source_eos & !target_eos, as.character(communication$target),
    ifelse(!source_eos & target_eos, as.character(communication$source), NA_character_)
  )
  short_allowed <- if (is.null(top_short)) rep(TRUE, nrow(communication)) else
    communication$neighbour_cell_type %in% as.character(top_short)
  long_allowed <- if (is.null(top_long)) rep(TRUE, nrow(communication)) else
    communication$neighbour_cell_type %in% as.character(top_long)
  communication$matches_state_top_neighbour <- ifelse(
    communication$eos_state_group == "Eos_short_enriched", short_allowed,
    ifelse(communication$eos_state_group == "Eos_long_enriched", long_allowed, FALSE)
  )
  communication$analysis_scope <- "EXPLORATORY_WITHIN_SECTION_CELLCHAT"
  significant <- communication[
    is.finite(communication$pval) & communication$pval < raw_p_max &
      is.finite(communication$p_adjust_bh) & communication$p_adjust_bh < adjusted_p_max &
      communication$direction == "EOS_TO_NEIGHBOUR" &
      communication$matches_state_top_neighbour,
    , drop = FALSE
  ]
  significant <- significant[order(significant$p_adjust_bh, -significant$prob), , drop = FALSE]
  rownames(significant) <- NULL
  list(all = communication, significant = significant)
}

#' Deterministically sample cell IDs within annotation groups
#'
#' @param ids Unique cell identifiers.
#' @param groups Group label for every identifier.
#' @param max_per_group Maximum retained cells per group.
#' @param seed Random seed.
#'
#' @return Character vector of sampled IDs in original input order.
sample_ids_by_group <- function(ids, groups, max_per_group = 1000L, seed = 1234L) {
  ids <- as.character(ids)
  groups <- as.character(groups)
  if (length(ids) != length(groups) || !length(ids)) {
    stop("ids and groups must have equal positive length.", call. = FALSE)
  }
  if (anyNA(ids) || anyNA(groups) || anyDuplicated(ids)) {
    stop("Sampling requires unique IDs and complete group labels.", call. = FALSE)
  }
  max_per_group <- as.integer(max_per_group)
  if (!is.finite(max_per_group) || max_per_group < 1L) {
    stop("max_per_group must be a positive integer.", call. = FALSE)
  }
  set.seed(as.integer(seed))
  selected <- unlist(lapply(unique(groups), function(group_name) {
    candidates <- ids[groups == group_name]
    if (length(candidates) <= max_per_group) candidates else sample(candidates, max_per_group)
  }), use.names = FALSE)
  ids[ids %in% selected]
}

#' Balance a Wang Seurat reference by harmonized subtype
#'
#' @param reference Wang Seurat reference object.
#' @param subtype_col Metadata column defining reference subtypes.
#' @param max_per_subtype Maximum cells retained per subtype.
#' @param seed Random seed.
#'
#' @return A cell-subset Seurat object with deterministic subtype balancing.
sample_wang_reference <- function(
    reference,
    subtype_col = "Wang_subtype_harmonized",
    max_per_subtype = 1000L,
    seed = 1234L
) {
  if (!inherits(reference, "Seurat")) stop("reference must be a Seurat object.", call. = FALSE)
  if (!subtype_col %in% colnames(reference@meta.data)) {
    stop("Wang reference metadata missing: ", subtype_col, call. = FALSE)
  }
  selected <- sample_ids_by_group(
    colnames(reference), reference@meta.data[[subtype_col]],
    max_per_group = max_per_subtype, seed = seed
  )
  subset(reference, cells = selected)
}

#' Prefix all Seurat cell IDs before cross-dataset integration
#'
#' @param object Seurat object.
#' @param prefix Non-empty dataset prefix without a trailing separator.
#'
#' @return A renamed Seurat object whose cell IDs are `<prefix>_<old_id>`.
prefix_seurat_cell_ids <- function(object, prefix) {
  if (!inherits(object, "Seurat")) stop("object must be a Seurat object.", call. = FALSE)
  prefix <- as.character(prefix)[[1L]]
  if (is.na(prefix) || !nzchar(prefix)) stop("prefix must be non-empty.", call. = FALSE)
  new_ids <- paste0(prefix, "_", colnames(object))
  if (anyDuplicated(new_ids)) stop("Prefixed Seurat IDs are not unique.", call. = FALSE)
  SeuratObject::RenameCells(object, new.names = new_ids)
}

#' Validate the Wang–Xenium shared feature set against the fixed panel
#'
#' @param reference_genes Genes present in the Wang reference assay.
#' @param query_genes Genes present in the Xenium query assay.
#' @param panel_genes Ordered complete Xenium panel genes.
#' @param min_shared Minimum shared genes required for integration.
#'
#' @return A list with typed status, ordered shared features and a per-panel-gene
#'   presence manifest.
validate_shared_feature_set <- function(
    reference_genes,
    query_genes,
    panel_genes,
    min_shared = 100L
) {
  panel <- unique(as.character(panel_genes))
  reference <- unique(as.character(reference_genes))
  query <- unique(as.character(query_genes))
  manifest <- data.frame(
    gene = panel,
    reference_present = panel %in% reference,
    query_present = panel %in% query,
    stringsAsFactors = FALSE
  )
  manifest$shared <- manifest$reference_present & manifest$query_present
  features <- manifest$gene[manifest$shared]
  status <- if (length(features) >= as.integer(min_shared)) {
    "PASS"
  } else {
    "SKIPPED_INSUFFICIENT_SHARED_GENES"
  }
  list(
    status = status,
    message = sprintf("%d of %d panel genes are shared; minimum required is %d.", length(features), length(panel), as.integer(min_shared)),
    features = features,
    manifest = manifest
  )
}

#' Run an exploratory Seurat CCA integration of Wang and Xenium
#'
#' @param reference Wang Seurat object.
#' @param query Xenium Seurat object.
#' @param features Explicit ordered shared feature vector.
#' @param dims Integration/PCA dimensions; defaults to PCs 1–30.
#' @param seed Random seed.
#' @param reference_assay Wang expression assay.
#' @param query_assay Xenium expression assay.
#' @param min_shared Minimum shared features required before fitting.
#' @param resolution Exploratory joint-clustering resolution.
#'
#' @return A typed result with stage, message, optional integrated Seurat object,
#'   anchors and parameter table. Original objects are not modified.
run_wang_xenium_integration <- function(
    reference,
    query,
    features,
    dims = 1:30,
    seed = 1234L,
    reference_assay = "RNA",
    query_assay = "Xenium",
    min_shared = 100L,
    resolution = 0.8
) {
  if (!inherits(reference, "Seurat") || !inherits(query, "Seurat")) {
    stop("reference and query must be Seurat objects.", call. = FALSE)
  }
  features <- unique(as.character(features))
  features <- features[
    features %in% rownames(reference[[reference_assay]]) &
      features %in% rownames(query[[query_assay]])
  ]
  parameters <- data.frame(
    method = "SEURAT_CCA_LOGNORMALIZE",
    n_shared_features = length(features),
    dims = paste(range(as.integer(dims)), collapse = "-"),
    seed = as.integer(seed),
    resolution = as.numeric(resolution),
    analysis_scope = "EXPLORATORY_WANG_XENIUM_CONCORDANCE",
    stringsAsFactors = FALSE
  )
  if (length(features) < as.integer(min_shared)) {
    return(list(
      status = "SKIPPED_INSUFFICIENT_SHARED_GENES",
      message = sprintf("Need at least %d shared genes; found %d.", as.integer(min_shared), length(features)),
      stage = "SHARED_GENE_GATE", object = NULL, anchors = NULL,
      parameters = parameters
    ))
  }
  stage <- "PREFIX_IDS"
  output <- tryCatch({
    reference_use <- prefix_seurat_cell_ids(reference, "WANG")
    query_use <- prefix_seurat_cell_ids(query, "XENIUM")
    reference_use$integration_dataset <- "WANG"
    query_use$integration_dataset <- "XENIUM"
    SeuratObject::DefaultAssay(reference_use) <- reference_assay
    SeuratObject::DefaultAssay(query_use) <- query_assay
    stage <- "NORMALIZE"
    reference_use <- Seurat::NormalizeData(reference_use, assay = reference_assay, verbose = FALSE)
    query_use <- Seurat::NormalizeData(query_use, assay = query_assay, verbose = FALSE)
    stage <- "FIND_ANCHORS"
    anchors <- Seurat::FindIntegrationAnchors(
      object.list = list(reference_use, query_use),
      assay = c(reference_assay, query_assay),
      anchor.features = features,
      normalization.method = "LogNormalize",
      reduction = "cca",
      dims = as.integer(dims),
      verbose = FALSE
    )
    stage <- "INTEGRATE_DATA"
    integrated <- Seurat::IntegrateData(
      anchorset = anchors,
      normalization.method = "LogNormalize",
      dims = as.integer(dims),
      features.to.integrate = features,
      verbose = FALSE
    )
    SeuratObject::DefaultAssay(integrated) <- "integrated"
    stage <- "PCA_UMAP_CLUSTER"
    integrated <- Seurat::ScaleData(integrated, features = features, verbose = FALSE)
    integrated <- Seurat::RunPCA(
      integrated, features = features, npcs = max(as.integer(dims)),
      seed.use = as.integer(seed), verbose = FALSE
    )
    integrated <- Seurat::FindNeighbors(integrated, reduction = "pca", dims = as.integer(dims), verbose = FALSE)
    integrated <- Seurat::FindClusters(
      integrated, resolution = resolution, random.seed = as.integer(seed),
      cluster.name = "wang_xenium_cluster", verbose = FALSE
    )
    integrated <- Seurat::RunUMAP(
      integrated, reduction = "pca", dims = as.integer(dims),
      reduction.name = "wang_xenium_umap", seed.use = as.integer(seed),
      verbose = FALSE
    )
    list(object = integrated, anchors = anchors)
  }, error = identity)
  if (inherits(output, "error")) {
    return(list(
      status = "FAILED_WANG_XENIUM_INTEGRATION",
      message = conditionMessage(output), stage = stage,
      object = NULL, anchors = NULL, parameters = parameters
    ))
  }
  list(
    status = "PASS",
    message = "Exploratory Wang–Xenium joint embedding completed.",
    stage = "COMPLETE", object = output$object, anchors = output$anchors,
    parameters = parameters
  )
}

#' Summarize cross-dataset Eosinophil neighbourhood concordance
#'
#' @param embeddings Numeric cell-by-dimension matrix with unique row names.
#' @param metadata Cell metadata aligned by row name.
#' @param dataset_col Two-level dataset label.
#' @param eos_col Complete logical Eosinophil indicator.
#' @param cluster_col Optional integrated-cluster field.
#' @param k Cross-dataset neighbour count.
#'
#' @return A list with per-cell cross-dataset distances/Eosinophil fractions,
#'   dataset-level summaries and dataset-by-cluster composition/enrichment.
summarise_cross_dataset_eos_neighbours <- function(
    embeddings,
    metadata,
    dataset_col = "dataset",
    eos_col = "is_eosinophil",
    cluster_col = "integrated_cluster",
    k = 15L
) {
  require_package("FNN")
  embedding <- as.matrix(embeddings)
  metadata <- as.data.frame(metadata, stringsAsFactors = FALSE)
  if (is.null(rownames(embedding)) || anyDuplicated(rownames(embedding))) {
    stop("embeddings require unique cell row names.", call. = FALSE)
  }
  missing <- setdiff(c(dataset_col, eos_col), names(metadata))
  if (length(missing)) stop("Cross-dataset metadata missing: ", paste(missing, collapse = ", "), call. = FALSE)
  index <- match(rownames(embedding), rownames(metadata))
  if (anyNA(index)) stop("Metadata does not cover every embedding cell.", call. = FALSE)
  metadata <- metadata[index, , drop = FALSE]
  dataset <- as.character(metadata[[dataset_col]])
  datasets <- unique(dataset)
  if (length(datasets) != 2L) stop("Exactly two integration datasets are required.", call. = FALSE)
  is_eos <- as.logical(metadata[[eos_col]])
  if (anyNA(is_eos)) stop(eos_col, " must be complete and logical.", call. = FALSE)
  per_dataset <- lapply(datasets, function(query_dataset) {
    query_index <- which(dataset == query_dataset)
    reference_index <- which(dataset != query_dataset)
    k_actual <- min(as.integer(k), length(reference_index))
    fit <- FNN::get.knnx(embedding[reference_index, , drop = FALSE], embedding[query_index, , drop = FALSE], k = k_actual)
    validate_neighbour_distances(fit$nn.dist)
    neighbour_index <- matrix(reference_index[fit$nn.index], nrow = nrow(fit$nn.index), ncol = ncol(fit$nn.index))
    data.frame(
      cell_id = rownames(embedding)[query_index],
      dataset = query_dataset,
      is_eosinophil = is_eos[query_index],
      k_actual = k_actual,
      minimum_cross_dataset_distance = fit$nn.dist[, 1L],
      eos_neighbour_fraction = rowMeans(matrix(is_eos[neighbour_index], nrow = nrow(neighbour_index))),
      stringsAsFactors = FALSE
    )
  })
  per_cell <- do.call(rbind, per_dataset)
  rownames(per_cell) <- NULL
  dataset_summary <- do.call(rbind, lapply(split(per_cell, per_cell$dataset), function(x) {
    data.frame(
      dataset = x$dataset[[1L]], n_cells = nrow(x),
      n_eosinophil = sum(x$is_eosinophil),
      mean_eos_neighbour_fraction_all = mean(x$eos_neighbour_fraction),
      mean_eos_neighbour_fraction_eos = if (any(x$is_eosinophil)) mean(x$eos_neighbour_fraction[x$is_eosinophil]) else NA_real_,
      median_cross_dataset_distance = stats::median(x$minimum_cross_dataset_distance),
      stringsAsFactors = FALSE
    )
  }))
  if (cluster_col %in% names(metadata)) {
    cluster_enrichment <- as.data.frame(table(
      dataset = dataset,
      integrated_cluster = as.character(metadata[[cluster_col]])
    ), stringsAsFactors = FALSE)
    names(cluster_enrichment)[names(cluster_enrichment) == "Freq"] <- "n_cells"
    eos_counts <- aggregate(
      as.integer(is_eos),
      by = list(dataset = dataset, integrated_cluster = as.character(metadata[[cluster_col]])),
      FUN = sum
    )
    names(eos_counts)[[3L]] <- "n_eosinophil"
    cluster_enrichment <- merge(
      cluster_enrichment, eos_counts,
      by = c("dataset", "integrated_cluster"), all.x = TRUE, sort = FALSE
    )
    cluster_enrichment$n_eosinophil[is.na(cluster_enrichment$n_eosinophil)] <- 0L
    cluster_enrichment$eosinophil_fraction <- ifelse(
      cluster_enrichment$n_cells > 0,
      cluster_enrichment$n_eosinophil / cluster_enrichment$n_cells,
      NA_real_
    )
  } else {
    cluster_enrichment <- data.frame()
  }
  list(
    per_cell = per_cell,
    dataset_summary = dataset_summary,
    cluster_enrichment = cluster_enrichment
  )
}

#' Calculate the adjusted Rand index for two cluster assignments
#'
#' @param labels_a,labels_b Equal-length complete cluster-label vectors.
#'
#' @return A numeric scalar in the adjusted-Rand scale. Identical partitions,
#'   including label permutations, return exactly 1.
adjusted_rand_index <- function(labels_a, labels_b) {
  if (length(labels_a) != length(labels_b) || !length(labels_a)) stop("Cluster label vectors must have equal positive length.", call. = FALSE)
  if (anyNA(labels_a) || anyNA(labels_b)) stop("Cluster labels cannot contain missing values.", call. = FALSE)
  tab <- table(as.character(labels_a), as.character(labels_b))
  choose2 <- function(x) x * (x - 1) / 2
  sum_cells <- sum(choose2(tab))
  sum_rows <- sum(choose2(rowSums(tab)))
  sum_cols <- sum(choose2(colSums(tab)))
  total <- choose2(length(labels_a))
  expected <- if (total == 0) 0 else sum_rows * sum_cols / total
  maximum <- (sum_rows + sum_cols) / 2
  denominator <- maximum - expected
  if (denominator == 0) return(if (identical(as.integer(tab > 0), as.integer(diag(nrow(tab)) > 0))) 1 else 0)
  as.numeric((sum_cells - expected) / denominator)
}

#' Summarize pairwise stability across repeated cluster assignments
#'
#' @param cluster_assignments Data frame or matrix with cells in rows and one
#'   clustering run, seed, or resolution in each column.
#'
#' @return A list with a pairwise adjusted-Rand table and one-row summary of
#'   comparison count, minimum, median, mean, and maximum ARI.
summarise_cluster_stability <- function(cluster_assignments) {
  x <- as.data.frame(cluster_assignments, stringsAsFactors = FALSE)
  if (ncol(x) < 2L || nrow(x) < 2L) stop("At least two cells and two clustering runs are required.", call. = FALSE)
  if (anyNA(x)) stop("Cluster assignments cannot contain missing values.", call. = FALSE)
  pairs <- utils::combn(colnames(x), 2L, simplify = FALSE)
  pairwise <- do.call(rbind, lapply(pairs, function(pair) {
    data.frame(run_a = pair[[1L]], run_b = pair[[2L]], ari = adjusted_rand_index(x[[pair[[1L]]]], x[[pair[[2L]]]]), stringsAsFactors = FALSE)
  }))
  rownames(pairwise) <- NULL
  list(
    pairwise = pairwise,
    summary = data.frame(
      n_comparisons = nrow(pairwise),
      minimum_ari = min(pairwise$ari),
      median_ari = stats::median(pairwise$ari),
      mean_ari = mean(pairwise$ari),
      maximum_ari = max(pairwise$ari),
      stringsAsFactors = FALSE
    )
  )
}

#' Calculate tidy PC-to-QC Spearman correlations
#'
#' @param pca_embeddings Numeric cell-by-PC matrix.
#' @param cell_metadata Data frame in exactly the same cell order.
#' @param qc_metrics Candidate metadata columns to test.
#'
#' @return A data frame with one row per PC/available-QC-metric pair and fields
#'   `pc`, `qc_metric`, `rho`, `abs_rho`, and `technical_review_flag`.
summarise_pca_qc_correlations <- function(
    pca_embeddings,
    cell_metadata,
    qc_metrics = c("nCount_Xenium", "nFeature_Xenium", "cell_area"),
    review_abs_rho = 0.5
) {
  embeddings <- as.matrix(pca_embeddings)
  if (!is.numeric(embeddings) || nrow(embeddings) != nrow(cell_metadata)) stop("PCA embeddings and metadata must have the same cell rows.", call. = FALSE)
  if (is.null(colnames(embeddings))) colnames(embeddings) <- paste0("PC_", seq_len(ncol(embeddings)))
  metrics <- intersect(qc_metrics, colnames(cell_metadata))
  if (!length(metrics)) stop("No requested QC metrics are present.", call. = FALSE)
  rows <- lapply(colnames(embeddings), function(pc) lapply(metrics, function(metric) {
    rho <- suppressWarnings(stats::cor(embeddings[, pc], cell_metadata[[metric]], method = "spearman", use = "pairwise.complete.obs"))
    data.frame(pc = pc, qc_metric = metric, rho = as.numeric(rho), abs_rho = abs(as.numeric(rho)), technical_review_flag = is.finite(rho) && abs(rho) >= review_abs_rho, stringsAsFactors = FALSE)
  }))
  do.call(rbind, unlist(rows, recursive = FALSE))
}

#' Derive a spatial lymph-node candidate domain from local lymphoid enrichment
#'
#' @param cell_metadata Data frame containing unique IDs, coordinates, and a
#'   conservative cell-type label for QC-passed cells.
#' @param lymphoid_labels Labels considered direct lymphoid/DC evidence.
#' @param cell_id_col,x_col,y_col,label_col Input column names.
#' @param k Number of neighbouring cells used for local enrichment.
#' @param lymphoid_fraction_threshold Minimum local fraction defining core cells.
#' @param expansion_radius Maximum coordinate-space distance used to include
#'   stromal/vascular cells surrounding an LN core.
#' @param min_core_cells Minimum enriched core size required to call a candidate.
#'
#' @return A list with status, parameters, a cell-level diagnostic table, and
#'   counts. `LN_NOT_DETECTED` returns an all-FALSE domain instead of forcing LN.
derive_lymph_node_domain <- function(
    cell_metadata,
    lymphoid_labels,
    cell_id_col = "cell_id",
    x_col = "x",
    y_col = "y",
    label_col = "cell_type",
    k = 30L,
    lymphoid_fraction_threshold = 0.5,
    expansion_radius = 80,
    min_core_cells = 100L
) {
  require_package("FNN")
  required <- c(cell_id_col, x_col, y_col, label_col)
  missing <- setdiff(required, colnames(cell_metadata))
  if (length(missing)) stop("Missing lymph-node-domain fields: ", paste(missing, collapse = ", "), call. = FALSE)
  ids <- as.character(cell_metadata[[cell_id_col]])
  if (anyDuplicated(ids) || anyNA(ids)) stop("Lymph-node-domain cell IDs must be unique and complete.", call. = FALSE)
  coordinates <- as.matrix(cell_metadata[c(x_col, y_col)])
  storage.mode(coordinates) <- "double"
  if (any(!is.finite(coordinates))) stop("Lymph-node-domain coordinates must be finite.", call. = FALSE)
  n <- nrow(coordinates)
  k <- as.integer(k)
  if (n < 3L || k < 1L || k >= n) stop("k must be between 1 and the number of cells minus one.", call. = FALSE)
  if (!is.finite(lymphoid_fraction_threshold) || lymphoid_fraction_threshold < 0 || lymphoid_fraction_threshold > 1) stop("lymphoid_fraction_threshold must be between zero and one.", call. = FALSE)
  if (!is.finite(expansion_radius) || expansion_radius < 0) stop("expansion_radius must be non-negative.", call. = FALSE)
  neighbour_index <- FNN::get.knn(coordinates, k = k)$nn.index
  lymphoid <- as.character(cell_metadata[[label_col]]) %in% lymphoid_labels
  local_fraction <- rowMeans(matrix(lymphoid[neighbour_index], nrow = n, ncol = k))
  core <- local_fraction >= lymphoid_fraction_threshold
  minimum <- as.integer(min_core_cells)
  status <- if (sum(core) >= minimum) "LN_CANDIDATE" else "LN_NOT_DETECTED"
  distance_to_core <- rep(Inf, n)
  include <- rep(FALSE, n)
  if (status == "LN_CANDIDATE") {
    nearest <- FNN::get.knnx(coordinates[core, , drop = FALSE], coordinates, k = 1L)
    distance_to_core <- as.numeric(nearest$nn.dist[, 1L])
    include <- distance_to_core <= expansion_radius
  }
  cell_table <- data.frame(
    cell_id = ids,
    direct_lymphoid_evidence = lymphoid,
    local_lymphoid_fraction = local_fraction,
    lymph_node_core = core,
    distance_to_lymph_node_core = distance_to_core,
    lymph_node_include = include,
    stringsAsFactors = FALSE
  )
  list(
    status = status,
    parameters = data.frame(k = k, lymphoid_fraction_threshold = lymphoid_fraction_threshold, expansion_radius = expansion_radius, min_core_cells = minimum),
    cell_table = cell_table,
    counts = data.frame(n_cells = n, n_direct_lymphoid = sum(lymphoid), n_core = sum(core), n_domain = sum(include), stringsAsFactors = FALSE)
  )
}

#' Quantify lymph-node boundary agreement across parameter settings
#'
#' @param masks Named list of equal-length complete logical inclusion vectors.
#'
#' @return A list with pairwise Jaccard overlap and a one-row summary. Empty-set
#'   agreement is defined as one only when both compared masks are empty.
summarise_ln_boundary_sensitivity <- function(masks) {
  if (!is.list(masks) || length(masks) < 2L || is.null(names(masks)) || any(!nzchar(names(masks)))) stop("Supply at least two named LN masks.", call. = FALSE)
  lengths <- vapply(masks, length, integer(1))
  if (length(unique(lengths)) != 1L || any(vapply(masks, function(x) !is.logical(x) || anyNA(x), logical(1)))) stop("LN masks must be complete logical vectors of equal length.", call. = FALSE)
  pairs <- utils::combn(names(masks), 2L, simplify = FALSE)
  pairwise <- do.call(rbind, lapply(pairs, function(pair) {
    a <- masks[[pair[[1L]]]]
    b <- masks[[pair[[2L]]]]
    union_n <- sum(a | b)
    jaccard <- if (union_n == 0L) 1 else sum(a & b) / union_n
    data.frame(mask_a = pair[[1L]], mask_b = pair[[2L]], jaccard = jaccard, n_a = sum(a), n_b = sum(b), stringsAsFactors = FALSE)
  }))
  rownames(pairwise) <- NULL
  list(
    pairwise = pairwise,
    summary = data.frame(n_comparisons = nrow(pairwise), minimum_jaccard = min(pairwise$jaccard), median_jaccard = stats::median(pairwise$jaccard), mean_jaccard = mean(pairwise$jaccard), stringsAsFactors = FALSE)
  )
}

#' Save and reload-validate a Seurat checkpoint
#'
#' @param object Seurat object to checkpoint without modifying it.
#' @param path Destination `.rds` path.
#' @param project_root Approved project root that must contain `path`.
#' @param stage Stable stage label written to the checkpoint manifest.
#' @param compress Compression argument passed to `saveRDS`; `FALSE` is faster
#'   for large HPC checkpoints and is the default.
#'
#' @return A one-row data frame containing stage, absolute path, byte size, MD5,
#'   cell/feature counts, timestamp, and `PASS` validation status. Validation
#'   reloads the object and requires identical class, dimensions, feature IDs,
#'   cell IDs, assay names, and reduction names.
write_validated_seurat_checkpoint <- function(
    object,
    path,
    project_root,
    stage,
    compress = FALSE
) {
  if (!inherits(object, "Seurat")) stop("object must be a Seurat object.", call. = FALSE)
  if (length(path) != 1L || is.na(path) || !nzchar(path)) stop("path must be one non-empty value.", call. = FALSE)
  if (length(stage) != 1L || is.na(stage) || !nzchar(stage)) stop("stage must be one non-empty value.", call. = FALSE)
  assert_path_within(project_root, path)
  parent <- dirname(path)
  assert_path_within(project_root, parent)
  dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  saveRDS(object, path, compress = compress)
  if (!file.exists(path)) stop("Checkpoint was not created: ", path, call. = FALSE)
  restored <- readRDS(path)
  checks <- c(
    inherits(restored, "Seurat"),
    identical(class(restored), class(object)),
    identical(dim(restored), dim(object)),
    identical(rownames(restored), rownames(object)),
    identical(colnames(restored), colnames(object)),
    identical(names(restored@assays), names(object@assays)),
    identical(names(restored@reductions), names(object@reductions))
  )
  if (!all(checks)) stop("Reloaded checkpoint failed Seurat identity validation: ", path, call. = FALSE)
  info <- file.info(path)
  data.frame(
    stage = stage,
    path = normalizePath(path, winslash = "/", mustWork = TRUE),
    bytes = as.numeric(info$size),
    md5 = unname(tools::md5sum(path)),
    n_cells = ncol(object),
    n_features = nrow(object),
    validation_status = "PASS",
    timestamp_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    stringsAsFactors = FALSE
  )
}

#' Import a Xenium region as a spatial Seurat object
#'
#' @param xenium_dir Existing 10x Xenium output directory containing the
#'   cell-feature matrix, cell centroids, and requested segmentation files.
#' @param genes Optional character vector of panel genes to retain.
#' @param cells Optional character vector of cell IDs to retain.
#' @param project Optional Seurat project name; defaults to the directory name.
#' @param assay Name assigned to the count assay.
#' @param fov Name assigned to the Xenium field of view.
#' @param include_cell_segmentation Whether to attach cell polygons when present.
#' @param include_nucleus_segmentation Whether to attach nucleus polygons when
#'   present.
#'
#' @return A Seurat object containing sparse Xenium counts, cell metadata,
#'   centroid coordinates, a spatial FOV, and each requested available boundary.
create_spatial_seurat_from_xenium <- function(
    xenium_dir,
    genes = NULL,
    cells = NULL,
    project = NULL,
    assay = "Xenium",
    fov = "fov",
    include_cell_segmentation = TRUE,
    include_nucleus_segmentation = TRUE
) {

  # ============================================================
  # 0. Basic checks
  # ============================================================

  if (!dir.exists(xenium_dir)) {
    stop(
      "xenium_dir does not exist: ",
      xenium_dir,
      call. = FALSE
    )
  }

  if (is.null(project)) {
    project <- basename(
      normalizePath(
        xenium_dir,
        winslash = "/",
        mustWork = TRUE
      )
    )
  }

  message("Reading Xenium output directly from:")
  message("  ", xenium_dir)


  # ============================================================
  # 1. Read native Xenium count matrix
  # ============================================================

  h5_path <- file.path(
    xenium_dir,
    "cell_feature_matrix.h5"
  )

  matrix_dir <- file.path(
    xenium_dir,
    "cell_feature_matrix"
  )

  if (file.exists(h5_path)) {

    message("Reading cell_feature_matrix.h5 ...")

    counts_raw <- Seurat::Read10X_h5(
      filename = h5_path,
      use.names = TRUE,
      unique.features = TRUE
    )

  } else if (dir.exists(matrix_dir)) {

    message("Reading cell_feature_matrix/ directory ...")

    counts_raw <- Seurat::Read10X(
      data.dir = matrix_dir,
      gene.column = 2
    )

  } else {

    stop(
      paste0(
        "Could not find either:\n",
        "  cell_feature_matrix.h5\n",
        "or\n",
        "  cell_feature_matrix/\n",
        "inside:\n",
        xenium_dir
      ),
      call. = FALSE
    )
  }


  # ============================================================
  # 2. Handle feature-type list if Read10X returns one
  # ============================================================

  if (is.list(counts_raw)) {

    message(
      "Read10X returned feature types: ",
      paste(
        names(counts_raw),
        collapse = ", "
      )
    )

    if ("Gene Expression" %in% names(counts_raw)) {

      counts <- counts_raw[["Gene Expression"]]

    } else if ("GeneExpression" %in% names(counts_raw)) {

      counts <- counts_raw[["GeneExpression"]]

    } else {

      # Show available names rather than silently guessing
      stop(
        "Could not identify Gene Expression matrix. ",
        "Available entries: ",
        paste(
          names(counts_raw),
          collapse = ", "
        ),
        call. = FALSE
      )
    }

  } else {

    counts <- counts_raw
  }

  rm(counts_raw)
  invisible(gc())


  if (!inherits(counts, "sparseMatrix")) {
    counts <- methods::as(
      counts,
      "dgCMatrix"
    )
  }

  if (is.null(rownames(counts)) ||
      is.null(colnames(counts))) {

    stop(
      "Xenium count matrix lacks gene/cell names.",
      call. = FALSE
    )
  }


  message(
    "Native matrix: ",
    format(nrow(counts), big.mark = ","),
    " features × ",
    format(ncol(counts), big.mark = ","),
    " cells"
  )


  # ============================================================
  # 3. Read native Xenium cell metadata
  # ============================================================

  cells_parquet <- file.path(
    xenium_dir,
    "cells.parquet"
  )

  cells_csv <- file.path(
    xenium_dir,
    "cells.csv.gz"
  )

  if (file.exists(cells_parquet)) {

    if (!requireNamespace(
      "arrow",
      quietly = TRUE
    )) {
      stop(
        "Package 'arrow' is required to read cells.parquet.",
        call. = FALSE
      )
    }

    message("Reading cells.parquet ...")

    metadata <- arrow::read_parquet(
      cells_parquet,
      as_data_frame = TRUE
    )

  } else if (file.exists(cells_csv)) {

    message("Reading cells.csv.gz ...")

    metadata <- read.csv(
      cells_csv,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

  } else {

    stop(
      "Neither cells.parquet nor cells.csv.gz was found.",
      call. = FALSE
    )
  }


  # ============================================================
  # 4. Validate native metadata
  # ============================================================

  required_meta <- c(
    "cell_id",
    "x_centroid",
    "y_centroid"
  )

  missing_meta <- setdiff(
    required_meta,
    colnames(metadata)
  )

  if (length(missing_meta)) {
    stop(
      "Missing required Xenium cell columns: ",
      paste(
        missing_meta,
        collapse = ", "
      ),
      call. = FALSE
    )
  }

  metadata$cell_id <- as.character(
    metadata$cell_id
  )

  if (anyDuplicated(metadata$cell_id)) {
    stop(
      "Duplicated cell IDs in Xenium metadata.",
      call. = FALSE
    )
  }


  # ============================================================
  # 5. Align native matrix and metadata
  # ============================================================

  common_cells <- intersect(
    colnames(counts),
    metadata$cell_id
  )

  message(
    "Matrix/metadata overlapping cells: ",
    format(
      length(common_cells),
      big.mark = ","
    )
  )

  if (!length(common_cells)) {

    stop(
      paste0(
        "No overlapping cell IDs between ",
        "cell_feature_matrix and cells metadata.\n",
        "First matrix cell: ",
        head(colnames(counts), 1),
        "\nFirst metadata cell: ",
        head(metadata$cell_id, 1)
      ),
      call. = FALSE
    )
  }


  # Optional explicit cell subset
  if (!is.null(cells)) {

    cells <- unique(
      as.character(cells)
    )

    common_cells <- intersect(
      common_cells,
      cells
    )

    if (!length(common_cells)) {
      stop(
        "No requested cells remain.",
        call. = FALSE
      )
    }
  }


  counts <- counts[
    ,
    common_cells,
    drop = FALSE
  ]

  metadata <- metadata[
    match(
      common_cells,
      metadata$cell_id
    ),
    ,
    drop = FALSE
  ]

  rownames(metadata) <- metadata$cell_id


  stopifnot(
    identical(
      colnames(counts),
      rownames(metadata)
    )
  )


  # ============================================================
  # 6. Optional gene subset
  # ============================================================

  if (!is.null(genes)) {

    genes <- unique(
      as.character(genes)
    )

    genes_present <- intersect(
      genes,
      rownames(counts)
    )

    genes_missing <- setdiff(
      genes,
      rownames(counts)
    )

    message(
      "Requested genes present: ",
      length(genes_present),
      "/",
      length(genes)
    )

    if (length(genes_missing)) {

      message(
        "Missing genes: ",
        paste(
          head(genes_missing, 20),
          collapse = ", "
        )
      )
    }

    if (!length(genes_present)) {
      stop(
        "None of the requested genes are present.",
        call. = FALSE
      )
    }

    counts <- counts[
      genes_present,
      ,
      drop = FALSE
    ]
  }


  # ============================================================
  # 7. Create Seurat object
  # ============================================================

  object <- Seurat::CreateSeuratObject(
    counts = counts,
    meta.data = metadata,
    assay = assay,
    project = project,
    min.cells = 0,
    min.features = 0
  )


  # ============================================================
  # 8. Add Xenium centroids
  # ============================================================

  centroids <- data.frame(
    x = as.numeric(metadata$x_centroid),
    y = as.numeric(metadata$y_centroid),
    row.names = metadata$cell_id,
    check.names = FALSE
  )

  spatial_centroids <-
    SeuratObject::CreateCentroids(
      coords = centroids
    )

  spatial_fov <-
    SeuratObject::CreateFOV(
      coords = spatial_centroids,
      type = "centroids",
      molecules = NULL,
      assay = assay,
      key = paste0(fov, "_")
    )

  object[[fov]] <- spatial_fov


  # ============================================================
  # 9. Helper to read native Xenium boundaries
  # ============================================================

  read_xenium_boundary <- function(
      path,
      selected_cells
  ) {

    if (!file.exists(path)) {
      return(NULL)
    }

    if (!requireNamespace(
      "arrow",
      quietly = TRUE
    )) {
      stop(
        "Package 'arrow' required for Xenium polygons.",
        call. = FALSE
      )
    }

    boundary <- arrow::read_parquet(
      path,
      as_data_frame = TRUE
    )

    if (!"cell_id" %in% colnames(boundary)) {
      stop(
        "Boundary file lacks cell_id: ",
        basename(path),
        call. = FALSE
      )
    }


    # Xenium versions can differ slightly
    x_candidates <- c(
      "vertex_x",
      "x",
      "x_location",
      "x_centroid"
    )

    y_candidates <- c(
      "vertex_y",
      "y",
      "y_location",
      "y_centroid"
    )

    x_col <- intersect(
      x_candidates,
      colnames(boundary)
    )

    y_col <- intersect(
      y_candidates,
      colnames(boundary)
    )

    if (!length(x_col) ||
        !length(y_col)) {

      stop(
        "Cannot identify x/y polygon columns in ",
        basename(path),
        call. = FALSE
      )
    }

    x_col <- x_col[[1]]
    y_col <- y_col[[1]]


    boundary <- boundary[
      as.character(boundary$cell_id) %in%
        selected_cells,
      ,
      drop = FALSE
    ]


    polygon <- data.frame(
      x = as.numeric(
        boundary[[x_col]]
      ),
      y = as.numeric(
        boundary[[y_col]]
      ),
      cell = as.character(
        boundary$cell_id
      ),
      stringsAsFactors = FALSE
    )


    polygon <- polygon[
      is.finite(polygon$x) &
      is.finite(polygon$y) &
      !is.na(polygon$cell) &
      nzchar(polygon$cell),
      ,
      drop = FALSE
    ]

    polygon
  }


  # ============================================================
  # 10. Native cell segmentation
  # ============================================================

  segmentation_loaded <- character()

  if (isTRUE(
    include_cell_segmentation
  )) {

    cell_boundary_path <- file.path(
      xenium_dir,
      "cell_boundaries.parquet"
    )

    if (file.exists(cell_boundary_path)) {

      cell_polygon <- read_xenium_boundary(
        cell_boundary_path,
        selected_cells = colnames(object)
      )

      message(
        "Creating cell segmentation: ",
        format(
          nrow(cell_polygon),
          big.mark = ","
        ),
        " vertices"
      )

      cell_segmentation <-
        SeuratObject::CreateSegmentation(
          coords = cell_polygon,
          compact = TRUE
        )

      object[[fov]][["segmentation"]] <-
        cell_segmentation

      segmentation_loaded <- c(
        segmentation_loaded,
        "segmentation"
      )

      rm(
        cell_polygon,
        cell_segmentation
      )

      invisible(gc())

    } else {

      warning(
        "cell_boundaries.parquet not found.",
        call. = FALSE
      )
    }
  }


  # ============================================================
  # 11. Native nucleus segmentation
  # ============================================================

  if (isTRUE(
    include_nucleus_segmentation
  )) {

    nucleus_boundary_path <- file.path(
      xenium_dir,
      "nucleus_boundaries.parquet"
    )

    if (file.exists(
      nucleus_boundary_path
    )) {

      nucleus_polygon <- read_xenium_boundary(
        nucleus_boundary_path,
        selected_cells = colnames(object)
      )

      message(
        "Creating nucleus segmentation: ",
        format(
          nrow(nucleus_polygon),
          big.mark = ","
        ),
        " vertices"
      )

      nucleus_segmentation <-
        SeuratObject::CreateSegmentation(
          coords = nucleus_polygon,
          compact = TRUE
        )

      object[[fov]][["nucleus_segmentation"]] <-
        nucleus_segmentation

      segmentation_loaded <- c(
        segmentation_loaded,
        "nucleus_segmentation"
      )

      rm(
        nucleus_polygon,
        nucleus_segmentation
      )

      invisible(gc())

    } else {

      warning(
        "nucleus_boundaries.parquet not found.",
        call. = FALSE
      )
    }
  }


  # ============================================================
  # 12. FOV defaults
  # ============================================================

  SeuratObject::DefaultFOV(object) <- fov

  available_boundaries <-
    SeuratObject::Boundaries(
      object[[fov]]
    )

  if ("centroids" %in%
      available_boundaries) {

    SeuratObject::DefaultBoundary(
      object[[fov]]
    ) <- "centroids"
  }


  # ============================================================
  # 13. Provenance
  # ============================================================

  object@misc$xenium_import <- list(

    source =
      "NATIVE_XENIUM_OUTPUT",

    xenium_dir =
      normalizePath(
        xenium_dir,
        winslash = "/",
        mustWork = TRUE
      ),

    expression_source =
      if (file.exists(h5_path)) {
        "cell_feature_matrix.h5"
      } else {
        "cell_feature_matrix/"
      },

    metadata_source =
      if (file.exists(cells_parquet)) {
        "cells.parquet"
      } else {
        "cells.csv.gz"
      },

    genes =
      rownames(object),

    n_genes =
      nrow(object),

    n_cells =
      ncol(object),

    spatial_boundaries =
      available_boundaries,

    segmentation_loaded =
      segmentation_loaded
  )


  # ============================================================
  # 14. Final validation
  # ============================================================

  if (!identical(
    colnames(object),
    rownames(object@meta.data)
  )) {
    stop(
      "Final Seurat/metadata alignment failed.",
      call. = FALSE
    )
  }


  spatial_cells <- Cells(
    object[[fov]]
  )

  missing_spatial <- setdiff(
    colnames(object),
    spatial_cells
  )

  if (length(missing_spatial)) {
    stop(
      length(missing_spatial),
      " Seurat cells missing from FOV.",
      call. = FALSE
    )
  }


  # ============================================================
  # 15. Summary
  # ============================================================

  message("")
  message(
    "Created native Xenium spatial Seurat object: ",
    project
  )

  message(
    "  Cells: ",
    format(
      ncol(object),
      big.mark = ","
    )
  )

  message(
    "  Genes: ",
    format(
      nrow(object),
      big.mark = ","
    )
  )

  message(
    "  FOV: ",
    fov
  )

  message(
    "  Boundaries: ",
    paste(
      available_boundaries,
      collapse = ", "
    )
  )

  return(object)
}


#' Attach externally generated QC masks to a Seurat object by cell ID
#'
#' @param object Seurat object whose column names are Xenium cell IDs.
#' @param masks Data frame containing one unique cell-ID column and mask fields.
#' @param cell_id_col Name of the mask-table cell-ID column.
#' @param cols Optional mask columns to attach; `NULL` uses every non-ID column.
#' @param overwrite Whether existing metadata columns may be replaced.
#' @param require_complete_match Whether every Seurat cell must occur in `masks`.
#'
#' @return The input Seurat object with requested mask fields aligned and added
#'   to `object@meta.data` by matched cell ID, never by input row position.
add_masks_to_seurat <- function(
    object,
    masks,
    cell_id_col = "cell_id",
    cols = NULL,
    overwrite = FALSE,
    require_complete_match = FALSE
) {

  # ============================================================
  # 1. Validate inputs
  # ============================================================

  if (!inherits(object, "Seurat")) {
    stop("object must be a Seurat object.", call. = FALSE)
  }

  if (!cell_id_col %in% colnames(masks)) {
    stop(
      "masks does not contain cell ID column: ",
      cell_id_col,
      call. = FALSE
    )
  }

  masks <- as.data.frame(
    masks,
    stringsAsFactors = FALSE
  )

  masks[[cell_id_col]] <- as.character(
    masks[[cell_id_col]]
  )

  if (anyDuplicated(masks[[cell_id_col]])) {
    stop(
      "Duplicated cell IDs detected in masks.",
      call. = FALSE
    )
  }


  # ============================================================
  # 2. Match masks to Seurat cells
  # ============================================================

  object_cells <- colnames(object)

  idx <- match(
    object_cells,
    masks[[cell_id_col]]
  )

  n_matched <- sum(!is.na(idx))
  n_missing <- sum(is.na(idx))

  message(
    "Seurat cells: ",
    format(length(object_cells), big.mark = ",")
  )

  message(
    "Matched to masks: ",
    format(n_matched, big.mark = ",")
  )

  message(
    "Unmatched Seurat cells: ",
    format(n_missing, big.mark = ",")
  )


  if (require_complete_match && n_missing > 0) {

    missing_cells <- object_cells[
      is.na(idx)
    ]

    stop(
      n_missing,
      " Seurat cells were not found in masks. Examples: ",
      paste(
        head(missing_cells, 10),
        collapse = ", "
      ),
      call. = FALSE
    )
  }


  # ============================================================
  # 3. Decide which columns to transfer
  # ============================================================

  if (is.null(cols)) {

    cols <- setdiff(
      colnames(masks),
      cell_id_col
    )

  } else {

    missing_cols <- setdiff(
      cols,
      colnames(masks)
    )

    if (length(missing_cols)) {
      warning(
        "Columns not present in masks: ",
        paste(missing_cols, collapse = ", "),
        call. = FALSE
      )
    }

    cols <- intersect(
      cols,
      colnames(masks)
    )

    cols <- setdiff(
      cols,
      cell_id_col
    )
  }


  # ============================================================
  # 4. Handle columns already present in Seurat metadata
  # ============================================================

  existing_cols <- intersect(
    cols,
    colnames(object@meta.data)
  )

  if (length(existing_cols) && !overwrite) {

    message(
      "Skipping existing metadata columns: ",
      paste(existing_cols, collapse = ", ")
    )

    cols <- setdiff(
      cols,
      existing_cols
    )
  }

  if (!length(cols)) {
    message("No new columns to add.")
    return(object)
  }


  # ============================================================
  # 5. Add columns in EXACT Seurat cell order
  # ============================================================

  for (v in cols) {

    x <- masks[[v]][idx]

    # Protect against problematic list columns
    if (is.list(x)) {
      warning(
        "Column '", v,
        "' is a list column; converting to character.",
        call. = FALSE
      )

      x <- vapply(
        x,
        function(z) {
          if (length(z) == 0 || all(is.na(z))) {
            NA_character_
          } else {
            paste(z, collapse = ";")
          }
        },
        character(1)
      )
    }

    object@meta.data[[v]] <- x
  }


  # ============================================================
  # 6. Verify alignment
  # ============================================================

  if (!identical(
    rownames(object@meta.data),
    colnames(object)
  )) {
    stop(
      "Seurat metadata/cell alignment changed unexpectedly.",
      call. = FALSE
    )
  }


  # ============================================================
  # 7. Summary
  # ============================================================

  message(
    "Added ",
    length(cols),
    " metadata columns."
  )

  message(
    "Metadata dimensions: ",
    nrow(object@meta.data),
    " × ",
    ncol(object@meta.data)
  )

  return(object)
}

#' Transfer Wang-reference main and subtype labels to a Xenium query
#'
#' @param reference Seurat reference object with an RNA assay and Wang labels.
#' @param query Seurat Xenium query object.
#' @param annotation_genes Character vector of candidate shared features.
#' @param main_col Reference metadata column containing main cell types.
#' @param subtype_col Reference metadata column containing harmonized subtypes.
#' @param prefix Prefix used for transferred prediction-score columns.
#'
#' @return A list with `main` and `subtype` prediction data frames, the Seurat
#'   anchor set, and the exact shared feature vector used for transfer.
run_wang_transfer <- function(
  reference,
  query,
  annotation_genes,
  main_col = "Wang_main",
  subtype_col = "Wang_subtype_harmonized",
  prefix
) {

  DefaultAssay(reference) <- "RNA"

  # Restrict features to genes present in BOTH objects
  features_use <- Reduce(
    intersect,
    list(
      annotation_genes,
      rownames(reference),
      rownames(query)
    )
  )

  message(
    prefix,
    ": using ",
    length(features_use),
    " shared annotation genes"
  )

  reference <- NormalizeData(
    reference,
    assay = "RNA",
    verbose = FALSE
  )

  reference <- ScaleData(
    reference,
    assay = "RNA",
    features = features_use,
    verbose = FALSE
  )

  reference <- RunPCA(
    reference,
    assay = "RNA",
    features = features_use,
    npcs = 30,
    seed.use = 1234,
    verbose = FALSE
  )

  anchors <- FindTransferAnchors(
    reference = reference,
    query = query,
    reference.assay = "RNA",
    query.assay = "Xenium",
    normalization.method = "LogNormalize",
    reduction = "pcaproject",
    features = features_use,
    dims = 1:30,
    verbose = FALSE
  )

  pred_main <- TransferData(
    anchorset = anchors,
    refdata = reference[[main_col]][, 1],
    dims = 1:30,
    verbose = FALSE
  )

  pred_subtype <- TransferData(
    anchorset = anchors,
    refdata = reference[[subtype_col]][, 1],
    dims = 1:30,
    verbose = FALSE
  )

  # Rename output columns so 2.5m / 12m / all can coexist
  colnames(pred_main) <- paste0(
    prefix,
    "_main_",
    colnames(pred_main)
  )

  colnames(pred_subtype) <- paste0(
    prefix,
    "_subtype_",
    colnames(pred_subtype)
  )

  list(
    main = pred_main,
    subtype = pred_subtype,
    anchors = anchors,
    features = features_use
  )
}

#' Invoke an mclust-style fitter with its BIC function in the caller frame
#'
#' `mclust::Mclust()` rewrites its call to the unqualified symbol
#' `mclustBIC` and evaluates that call in its caller. This small adapter makes
#' that dependency explicit without attaching the package to the search path.
#'
#' @param data Finite numeric observations supplied to the mixture fitter.
#' @param G Integer candidate component counts.
#' @param mclust_fun Function with the public `Mclust()` calling contract.
#' @param mclust_bic_fun Function with the public `mclustBIC()` contract.
#' @param verbose Whether the fitter may print progress.
#'
#' @return The fitted object returned by `mclust_fun`.
run_mclust_with_binding <- function(
    data,
    G,
    mclust_fun,
    mclust_bic_fun,
    verbose = FALSE
) {
  mclustBIC <- mclust_bic_fun
  mclust_fun(data = data, G = G, verbose = verbose)
}

#' Run a namespace-safe Gaussian-mixture diagnostic for an Eosinophil score
#'
#' This optional diagnostic never creates biological state labels. Runtime and
#' sample-size failures are returned as typed statuses so they cannot terminate
#' an otherwise valid 479-gene notebook.
#'
#' @param x Numeric state-score vector. Non-finite entries are excluded and
#'   counted.
#' @param G Integer candidate component counts.
#' @param seed Random seed used by mclust.
#' @param min_n Minimum finite observations required to attempt fitting.
#' @param min_per_component Minimum observations required per candidate
#'   component; candidates exceeding this support are removed.
#'
#' @return A list containing status, message, package version, input/finite
#'   counts, selected model information, the optional fitted object and a BIC
#'   table.
run_mclust_diagnostic <- function(
    x,
    G = 1:3,
    seed = 1234L,
    min_n = 20L,
    min_per_component = 5L
) {
  x <- as.numeric(x)
  x_use <- x[is.finite(x)]
  available <- requireNamespace("mclust", quietly = TRUE)
  result <- list(
    status = NA_character_,
    message = NA_character_,
    package_version = if (available) {
      as.character(utils::packageVersion("mclust"))
    } else {
      NA_character_
    },
    n_input = length(x),
    n_finite = length(x_use),
    selected_G = NA_integer_,
    model_name = NA_character_,
    fit = NULL,
    bic_table = data.frame()
  )

  if (!available) {
    result$status <- "SKIPPED_PACKAGE_UNAVAILABLE"
    result$message <- "Optional package 'mclust' is unavailable."
    return(result)
  }

  if (length(x_use) < as.integer(min_n)) {
    result$status <- "SKIPPED_INSUFFICIENT_DATA"
    result$message <- sprintf(
      "Need at least %d finite observations; found %d.",
      as.integer(min_n), length(x_use)
    )
    return(result)
  }

  G_use <- sort(unique(as.integer(G)))
  G_use <- G_use[
    is.finite(G_use) & G_use >= 1L &
      G_use * as.integer(min_per_component) <= length(x_use)
  ]
  if (!length(G_use)) {
    result$status <- "SKIPPED_INSUFFICIENT_DATA"
    result$message <- "No requested component count has adequate observations."
    return(result)
  }

  set.seed(as.integer(seed))
  fit_or_error <- tryCatch(
    run_mclust_with_binding(
      data = x_use,
      G = G_use,
      mclust_fun = getExportedValue("mclust", "Mclust"),
      mclust_bic_fun = getExportedValue("mclust", "mclustBIC"),
      verbose = FALSE
    ),
    error = identity
  )
  if (inherits(fit_or_error, "error")) {
    result$status <- "FAILED_MCLUST_RUNTIME"
    result$message <- conditionMessage(fit_or_error)
    return(result)
  }

  result$status <- "PASS"
  result$message <- paste(
    "Gaussian-mixture diagnostic completed;",
    "the fitted components are not biological state assignments."
  )
  result$selected_G <- as.integer(fit_or_error$G)
  result$model_name <- as.character(fit_or_error$modelName)
  result$fit <- fit_or_error
  result$bic_table <- as.data.frame(fit_or_error$BIC)
  result
}

#' Plot Eosinophil evidence-call composition within annotated subtypes
#'
#' @param object Seurat object containing subtype and Eosinophil-call metadata.
#' @param subtype_col Metadata column holding final or provisional subtypes.
#' @param eos_call_col Metadata column holding ordered Eosinophil evidence calls.
#' @param eos_first Subtype displayed first in the plotted ordering.
#' @param title Plot title.
#' @param base_size Base ggplot text size.
#' @param show_n Whether subtype labels include cell counts.
#' @param return_data Whether to return plotting data with the plot.
#'
#' @return A ggplot object, or when `return_data = TRUE`, a list containing the
#'   plot, its summarized plotting data, subtype order, and displayed labels.
plot_eos_call_by_subtype <- function(
    object,
    subtype_col = "Final_CellType_subtype",
    eos_call_col = "EosRef_call",
    eos_first = "Eosinophil",
    title = "Eosinophil marker evidence across Xenium cell types",
    base_size = 11,
    show_n = TRUE,
    return_data = FALSE
) {

  # ============================================================
  # Packages
  # ============================================================

  requireNamespace("dplyr")
  requireNamespace("tidyr")
  requireNamespace("ggplot2")
  requireNamespace("scales")


  # ============================================================
  # 1. Check metadata columns
  # ============================================================

  md <- object@meta.data

  required_cols <- c(
    subtype_col,
    eos_call_col
  )

  missing_cols <- setdiff(
    required_cols,
    colnames(md)
  )

  if (length(missing_cols)) {
    stop(
      "Missing metadata columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }


  # ============================================================
  # 2. Eos call order / colours / labels
  # ============================================================

  call_order <- c(
    "REF_EOS_TIER1",
    "REF_EOS_TIER2",
    "REF_EOS_TIER3",
    "REF_EOS_TIER4",
    "REF_EOS_REST",
    "REVIEW_NONIMMUNE_EOS_RESCUE",
    "OUTSIDE_IMMUNE"
  )

  call_cols <- c(
    "REF_EOS_TIER1"               = "#E89A8F",
    "REF_EOS_TIER2"               = "#F4C6C3",
    "REF_EOS_TIER3"               = "#F4D6A0",
    "REF_EOS_TIER4"               = "#BFDCC8",
    "REF_EOS_REST"                = "#D9E5DF",
    "REVIEW_NONIMMUNE_EOS_RESCUE" = "#B9B3D7",
    "OUTSIDE_IMMUNE"              = "#E8E8E8"
  )

  call_labels <- c(
    "REF_EOS_TIER1"               = "Tier 1",
    "REF_EOS_TIER2"               = "Tier 2",
    "REF_EOS_TIER3"               = "Tier 3",
    "REF_EOS_TIER4"               = "Tier 4",
    "REF_EOS_REST"                = "Rest",
    "REVIEW_NONIMMUNE_EOS_RESCUE" = "Non-immune review",
    "OUTSIDE_IMMUNE"              = "Outside immune"
  )


  # ============================================================
  # 3. Prepare / summarise metadata
  # ============================================================

  plot_df <- md %>%
    dplyr::filter(
      !is.na(.data[[subtype_col]]),
      !is.na(.data[[eos_call_col]])
    ) %>%

    dplyr::transmute(
      subtype = as.character(
        .data[[subtype_col]]
      ),

      eos_call = factor(
        as.character(
          .data[[eos_call_col]]
        ),
        levels = call_order
      )
    ) %>%

    dplyr::count(
      subtype,
      eos_call,
      name = "n"
    ) %>%

    tidyr::complete(
      subtype,
      eos_call = factor(
        call_order,
        levels = call_order
      ),
      fill = list(
        n = 0
      )
    ) %>%

    dplyr::group_by(
      subtype
    ) %>%

    dplyr::mutate(
      n_total = sum(n),

      fraction = dplyr::if_else(
        n_total > 0,
        n / n_total,
        0
      ),

      percent = 100 * fraction
    ) %>%

    dplyr::ungroup()


  # ============================================================
  # 4. Order cell types
  #
  # Eosinophil first, then all other subtypes by cumulative
  # Tier 1-4 evidence.
  # ============================================================

  subtype_order <- plot_df %>%
    dplyr::group_by(
      subtype
    ) %>%

    dplyr::summarise(
      eos_evidence_pct = sum(
        percent[
          eos_call %in% c(
            "REF_EOS_TIER1",
            "REF_EOS_TIER2",
            "REF_EOS_TIER3",
            "REF_EOS_TIER4"
          )
        ],
        na.rm = TRUE
      ),

      .groups = "drop"
    ) %>%

    dplyr::mutate(
      eos_first_order = dplyr::if_else(
        subtype == eos_first,
        0L,
        1L
      )
    ) %>%

    dplyr::arrange(
      eos_first_order,
      dplyr::desc(
        eos_evidence_pct
      )
    ) %>%

    dplyr::pull(
      subtype
    )


  # coord_flip() places the last factor level at the top,
  # so reverse the desired top-to-bottom ordering.
  plot_df <- plot_df %>%
    dplyr::mutate(
      subtype = factor(
        subtype,
        levels = rev(
          subtype_order
        )
      )
    )


  # ============================================================
  # 5. Add cell numbers to subtype labels
  # ============================================================

  subtype_labels <- plot_df %>%
    dplyr::distinct(
      subtype,
      n_total
    ) %>%

    dplyr::mutate(
      subtype_chr =
        as.character(subtype),

      label = if (show_n) {

        paste0(
          subtype_chr,
          "  (n=",
          scales::comma(n_total),
          ")"
        )

      } else {

        subtype_chr
      }
    )


  subtype_label_vec <- stats::setNames(
    subtype_labels$label,
    subtype_labels$subtype_chr
  )


  # ============================================================
  # 6. Plot
  # ============================================================

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = subtype,
      y = fraction,
      fill = eos_call
    )
  ) +

    ggplot2::geom_col(
      position =
        ggplot2::position_stack(
          reverse = TRUE
        ),
      width = 0.78,
      colour = "black",
      linewidth = 0.25
    ) +

    ggplot2::coord_flip() +

    ggplot2::scale_fill_manual(
      values = call_cols,
      breaks = call_order,
      labels = call_labels,
      drop = FALSE
    ) +

    ggplot2::scale_x_discrete(
      labels = subtype_label_vec
    ) +

    ggplot2::scale_y_continuous(
      labels =
        scales::percent_format(
          accuracy = 1
        ),
      breaks =
        seq(
          0,
          1,
          by = 0.2
        ),
      limits = c(
        0,
        1
      ),
      expand = c(
        0,
        0
      )
    ) +

    ggplot2::labs(
      title = title,
      x = NULL,
      y = "Cells within subtype (%)",
      fill = "Eosinophil\nevidence"
    ) +

    cell_style_theme(
      base_size = base_size
    ) +

    ggplot2::theme(
      legend.position = "right",

      axis.text.y =
        ggplot2::element_text(
          size = 9,
          colour = "black"
        ),

      plot.title =
        ggplot2::element_text(
          face = "bold",
          size = 13
        )
    )


  # ============================================================
  # 7. Return
  # ============================================================

  if (return_data) {

    return(
      list(
        plot = p,
        plot_data = plot_df,
        subtype_order = subtype_order,
        subtype_labels = subtype_labels
      )
    )
  }

  p
}

#' Draw a continuously ordered Eosinophil-state expression heatmap
#'
#' @param eos_obj Seurat object restricted to the Eosinophil analysis cohort.
#' @param eos_gene_sets Data frame with `gene` and `gene_set` columns.
#' @param state_col Metadata column containing descriptive state categories.
#' @param balance_col Numeric continuous state-balance metadata column.
#' @param assay Seurat assay containing normalized expression.
#' @param layer Assay layer used for the heatmap.
#' @param remove_short_ribosomal Whether to remove ribosomal genes from the
#'   short-lived signature before plotting.
#' @param z_cap Absolute cap applied to row-wise expression z-scores.
#' @param cluster_rows Whether genes are clustered within signature blocks.
#' @param clustering_distance_rows ComplexHeatmap row-distance setting.
#' @param clustering_method_rows Hierarchical clustering method for rows.
#' @param column_title Heatmap title.
#' @param show_row_names Whether gene names are drawn.
#' @param row_name_size Gene-label font size.
#' @param draw_heatmap Whether to draw immediately.
#' @param return_data Whether to return matrices, orders, and annotations.
#'
#' @return If `return_data = TRUE`, a list containing the heatmap, optional
#'   drawn object, row-z-scored matrix, cell order, gene groups, and annotation
#'   vectors. Otherwise returns the heatmap object, invisibly returning the
#'   drawn object when `draw_heatmap = TRUE`.
plot_eos_state_heatmap <- function(
    eos_obj,
    eos_gene_sets,
    state_col = "EosState_extreme",
    balance_col = "EosState_balance",
    assay = "Xenium",
    layer = "data",
    remove_short_ribosomal = TRUE,
    z_cap = 2.5,
    cluster_rows = TRUE,
    clustering_distance_rows = "pearson",
    clustering_method_rows = "ward.D2",
    column_title = "Eosinophils ordered from short-lived-like to long-lived-like",
    show_row_names = TRUE,
    row_name_size = 7,
    draw_heatmap = TRUE,
    return_data = FALSE
) {

  # ============================================================
  # Packages
  # ============================================================

  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    stop("Package 'ComplexHeatmap' is required.", call. = FALSE)
  }

  if (!requireNamespace("circlize", quietly = TRUE)) {
    stop("Package 'circlize' is required.", call. = FALSE)
  }

  if (!requireNamespace("grid", quietly = TRUE)) {
    stop("Package 'grid' is required.", call. = FALSE)
  }


  # ============================================================
  # 1. Validate inputs
  # ============================================================

  if (!inherits(eos_obj, "Seurat")) {
    stop("eos_obj must be a Seurat object.", call. = FALSE)
  }

  required_gene_cols <- c(
    "gene",
    "gene_set"
  )

  missing_gene_cols <- setdiff(
    required_gene_cols,
    colnames(eos_gene_sets)
  )

  if (length(missing_gene_cols)) {
    stop(
      "eos_gene_sets missing columns: ",
      paste(missing_gene_cols, collapse = ", "),
      call. = FALSE
    )
  }

  required_meta <- c(
    state_col,
    balance_col
  )

  missing_meta <- setdiff(
    required_meta,
    colnames(eos_obj@meta.data)
  )

  if (length(missing_meta)) {
    stop(
      "eos_obj metadata missing columns: ",
      paste(missing_meta, collapse = ", "),
      call. = FALSE
    )
  }

  if (!assay %in% names(eos_obj@assays)) {
    stop(
      "Assay not found: ",
      assay,
      call. = FALSE
    )
  }


  # ============================================================
  # 2. Define short- and long-lived signatures
  # ============================================================

  short_genes <- unique(
    as.character(
      eos_gene_sets$gene[
        eos_gene_sets$gene_set == "short_lived"
      ]
    )
  )

  long_genes <- unique(
    as.character(
      eos_gene_sets$gene[
        eos_gene_sets$gene_set == "long_lived"
      ]
    )
  )

  short_use <- intersect(
    short_genes,
    rownames(eos_obj[[assay]])
  )

  long_use <- intersect(
    long_genes,
    rownames(eos_obj[[assay]])
  )


  if (remove_short_ribosomal) {

    short_primary <- short_use[
      !grepl(
        "^Rp[ls]",
        short_use
      )
    ]

  } else {

    short_primary <- short_use
  }


  state_genes_primary <- unique(
    c(
      short_primary,
      long_use
    )
  )


  message(
    "Short-lived genes available: ",
    length(short_use),
    " / ",
    length(short_genes)
  )

  message(
    "Short-lived genes used: ",
    length(short_primary)
  )

  message(
    "Long-lived genes available: ",
    length(long_use),
    " / ",
    length(long_genes)
  )

  message(
    "Total heatmap genes: ",
    length(state_genes_primary)
  )


  if (!length(state_genes_primary)) {
    stop(
      "No state genes are present in the object.",
      call. = FALSE
    )
  }


  # ============================================================
  # 3. Extract normalized expression
  # ============================================================

  expr <- SeuratObject::GetAssayData(
    eos_obj,
    assay = assay,
    layer = layer
  )

  expr <- as.matrix(
    expr[
      state_genes_primary,
      ,
      drop = FALSE
    ]
  )


  # ============================================================
  # 4. Row z-score
  # ============================================================

  expr_z <- t(
    scale(
      t(expr)
    )
  )


  # Remove zero-variance / non-finite genes
  keep_gene <- apply(
    expr_z,
    1,
    function(x) {
      all(is.finite(x))
    }
  )

  removed_zero_variance <- rownames(expr_z)[
    !keep_gene
  ]

  if (length(removed_zero_variance)) {
    message(
      "Removed ",
      length(removed_zero_variance),
      " zero-variance/non-finite genes."
    )
  }

  expr_z <- expr_z[
    keep_gene,
    ,
    drop = FALSE
  ]


  # Cap extreme z-scores
  expr_z[
    expr_z > z_cap
  ] <- z_cap

  expr_z[
    expr_z < -z_cap
  ] <- -z_cap


  # ============================================================
  # 5. Order cells by continuous Eos state balance
  # ============================================================

  md <- eos_obj@meta.data


  # Only retain cells with a valid balance score
  valid_cells <- rownames(md)[
    is.finite(
      md[[balance_col]]
    )
  ]


  cell_order_state <- valid_cells[
    order(
      md[
        valid_cells,
        balance_col
      ],
      decreasing = FALSE
    )
  ]


  # Ensure cells actually exist in expression matrix
  cell_order_state <- intersect(
    cell_order_state,
    colnames(expr_z)
  )


  expr_z <- expr_z[
    ,
    cell_order_state,
    drop = FALSE
  ]


  message(
    "Cells shown in heatmap: ",
    ncol(expr_z)
  )


  # ============================================================
  # 6. Column annotations
  # ============================================================

  state_annotation <- as.character(
    md[
      cell_order_state,
      state_col
    ]
  )

  state_annotation <- factor(
    state_annotation,
    levels = c(
      "Short-lived-like",
      "Intermediate",
      "Long-lived-like"
    )
  )


  balance_annotation <- as.numeric(
    md[
      cell_order_state,
      balance_col
    ]
  )


  state_cols <- c(
    "Short-lived-like" = "#7E9AD9",
    "Intermediate"     = "#D9D9D9",
    "Long-lived-like"  = "#E89A8F"
  )


  balance_lim <- max(
    abs(balance_annotation),
    na.rm = TRUE
  )


  if (!is.finite(balance_lim) ||
      balance_lim == 0) {

    balance_lim <- 1
  }


  balance_fun <- circlize::colorRamp2(
    c(
      -balance_lim,
      0,
      balance_lim
    ),
    c(
      "#7E9AD9",
      "#F7F7F7",
      "#E89A8F"
    )
  )


  ha_top <- ComplexHeatmap::HeatmapAnnotation(

    State = state_annotation,

    Balance = ComplexHeatmap::anno_simple(
      balance_annotation,
      col = balance_fun,
      border = FALSE
    ),

    col = list(
      State = state_cols
    ),

    annotation_name_gp = grid::gpar(
      fontsize = 9,
      fontface = "bold"
    ),

    annotation_legend_param = list(
      State = list(
        title = "Eosinophil state"
      )
    )
  )


  # ============================================================
  # 7. Row / gene-set annotation
  # ============================================================

  gene_state <- ifelse(
    rownames(expr_z) %in%
      short_primary,
    "Short-lived signature",
    "Long-lived signature"
  )


  gene_state <- factor(
    gene_state,
    levels = c(
      "Short-lived signature",
      "Long-lived signature"
    )
  )


  gene_set_cols <- c(
    "Short-lived signature" = "#7E9AD9",
    "Long-lived signature"  = "#E89A8F"
  )


  ha_row <- ComplexHeatmap::rowAnnotation(

    Signature = gene_state,

    col = list(
      Signature = gene_set_cols
    ),

    show_annotation_name = FALSE
  )


  # ============================================================
  # 8. Expression colour scale
  # ============================================================

  expr_col_fun <- circlize::colorRamp2(
    c(
      -z_cap,
      0,
      z_cap
    ),
    c(
      "#5B7DB1",
      "#F7F7F7",
      "#D77A72"
    )
  )


  # ============================================================
  # 9. Heatmap
  # ============================================================

  ht <- ha_row +

    ComplexHeatmap::Heatmap(

      expr_z,

      name = "Row z-score",

      col = expr_col_fun,

      top_annotation = ha_top,


      # --------------------------------------------------------
      # Cells remain ordered continuously by balance
      # --------------------------------------------------------

      cluster_columns = FALSE,


      # --------------------------------------------------------
      # Genes split into Short / Long signatures
      # --------------------------------------------------------

      row_split = gene_state,

      cluster_rows = cluster_rows,

      clustering_distance_rows =
        clustering_distance_rows,

      clustering_method_rows =
        clustering_method_rows,


      # --------------------------------------------------------
      # Appearance
      # --------------------------------------------------------

      show_column_names = FALSE,

      show_row_names = show_row_names,

      row_names_gp = grid::gpar(
        fontsize = row_name_size
      ),

      row_names_side = "left",

      column_title = column_title,

      column_title_gp = grid::gpar(
        fontsize = 12,
        fontface = "bold"
      ),

      row_title_gp = grid::gpar(
        fontsize = 10,
        fontface = "bold"
      ),

      row_gap = grid::unit(
        2,
        "mm"
      ),

      border = FALSE,

      rect_gp = grid::gpar(
        col = NA
      ),

      heatmap_legend_param = list(

        title = "Expression\n(z-score)",

        at = c(
          -2,
          0,
          2
        )
      )
    )


  # ============================================================
  # 10. Draw
  # ============================================================

  ht_drawn <- NULL

  if (draw_heatmap) {

    ht_drawn <- ComplexHeatmap::draw(

      ht,

      heatmap_legend_side = "right",

      annotation_legend_side = "right",

      merge_legends = TRUE
    )
  }


  # ============================================================
  # 11. Return
  # ============================================================

  if (return_data) {

    return(
      list(

        heatmap = ht,

        heatmap_drawn = ht_drawn,

        expression_z = expr_z,

        cell_order = cell_order_state,

        gene_state = gene_state,

        short_genes = short_primary,

        long_genes = long_use,

        removed_genes = removed_zero_variance,

        state_annotation = state_annotation,

        balance_annotation = balance_annotation
      )
    )
  }


  if (draw_heatmap) {
    return(
      invisible(ht_drawn)
    )
  }

  ht
}
