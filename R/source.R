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

`%||%` <- function(x, y) {
  if (is.null(x) || !length(x)) return(y)
  first <- x[[1]]
  missing_scalar <- is.atomic(first) && length(first) == 1L && (is.na(first) || !nzchar(as.character(first)))
  if (missing_scalar) y else x
}

# Purpose: Canonical path.
# Inputs: required: path.
# Output: Returns a deterministic scalar, vector, path, status, or empty-schema object used by downstream functions.
canonical_path <- function(path) {
  value <- normalizePath(path, winslash = "/", mustWork = FALSE)
  value <- sub("/+$", "", value)
  if (.Platform$OS.type == "windows") value <- tolower(value)
  value
}

# Purpose: Assert path within.
# Inputs: required: project_root, candidate.
# Output: Returns validation evidence/TRUE (or the validated value) and stops with an informative error when the contract fails.
assert_path_within <- function(project_root, candidate) {
  root <- canonical_path(project_root)
  path <- canonical_path(candidate)
  valid <- identical(path, root) || startsWith(path, paste0(root, "/"))
  if (!valid) stop(sprintf("Unsafe path outside project root: %s", path), call. = FALSE)
  invisible(TRUE)
}

# Purpose: Validate runtime paths.
# Inputs: required: project_root, input_root, output_root; optional/defaulted: temp_root.
# Output: Returns validation evidence/TRUE (or the validated value) and stops with an informative error when the contract fails.
validate_runtime_paths <- function(project_root, input_root, output_root, temp_root = tempdir()) {
  if (!dir.exists(project_root)) stop(sprintf("Project root does not exist: %s", project_root), call. = FALSE)
  if (!dir.exists(input_root)) stop(sprintf("Input root does not exist: %s", input_root), call. = FALSE)
  for (path in c(input_root, output_root, temp_root)) assert_path_within(project_root, path)
  invisible(TRUE)
}

# Purpose: Discover xenium sections.
# Inputs: required: input_root; optional/defaulted: expected_section_count.
# Output: Returns the derived R object described by the function name; no files are written unless an explicit output path is an input.
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

# Purpose: Discover one section.
# Inputs: required: input_root, region_id.
# Output: Returns the derived R object described by the function name; no files are written unless an explicit output path is an input.
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

# Purpose: Create synthetic manifest.
# Inputs: required: region_ids; optional/defaulted: seed.
# Output: Returns a newly constructed or annotated R object while preserving the supplied raw object/data rows.
create_synthetic_manifest <- function(region_ids, seed = 20260814L) {
  region_ids <- sort(unique(as.character(region_ids)))
  if (length(region_ids) != 4L) stop("Synthetic scWAT metadata requires exactly four regions.", call. = FALSE)
  design <- expand.grid(mouse_id = c("SyntheticMouse_1", "SyntheticMouse_2"), side = c("Left", "Right"), stringsAsFactors = FALSE)
  set.seed(as.integer(seed))
  design <- design[sample(seq_len(nrow(design))), , drop = FALSE]
  data.frame(
    tissue = "scWAT", region_id = region_ids, mouse_id = design$mouse_id, side = design$side,
    section_id = paste0("SyntheticSection_", seq_along(region_ids)),
    biological_replicate_id = design$mouse_id, technical_replicate_id = region_ids,
    sex = "UNKNOWN_SYNTHETIC", age = "UNKNOWN_SYNTHETIC", genotype = "WT_PLACEHOLDER",
    diet_or_condition = "UNKNOWN_SYNTHETIC", metadata_status = "SYNTHETIC_PLACEHOLDER",
    synthetic_seed = as.integer(seed), do_not_interpret = TRUE,
    notes = "Temporary randomized mapping for pipeline testing only; replace before biological analysis.",
    stringsAsFactors = FALSE
  )
}

# Purpose: Validate sample manifest.
# Inputs: required: manifest, expected_regions.
# Output: Returns validation evidence/TRUE (or the validated value) and stops with an informative error when the contract fails.
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

# Purpose: Xenium required files.
# Inputs: none.
# Output: Returns a deterministic scalar, vector, path, status, or empty-schema object used by downstream functions.
xenium_required_files <- function() {
  c(
    "experiment.xenium", "metrics_summary.csv", "analysis_summary.html", "gene_panel.json",
    "cells.csv.gz", "cell_feature_matrix/features.tsv.gz",
    "cell_feature_matrix/barcodes.tsv.gz", "cell_feature_matrix/matrix.mtx.gz"
  )
}

# Purpose: Inventory section files.
# Inputs: required: region_dir, region_id; optional/defaulted: calculate_md5.
# Output: Returns the derived R object described by the function name; no files are written unless an explicit output path is an input.
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

# Purpose: Read mtx dimensions.
# Inputs: required: path.
# Output: Returns parsed, validated R data (vector, data frame, sparse matrix bundle, or named list according to the input format).
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

# Purpose: Read gz rows.
# Inputs: required: path; optional/defaulted: header.
# Output: Returns parsed, validated R data (vector, data frame, sparse matrix bundle, or named list according to the input format).
read_gz_rows <- function(path, header = FALSE) {
  utils::read.delim(gzfile(path), header = header, quote = "", check.names = FALSE, stringsAsFactors = FALSE)
}

# Purpose: Validate section integrity.
# Inputs: required: region_dir, region_id.
# Output: Returns validation evidence/TRUE (or the validated value) and stops with an informative error when the contract fails.
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

# Purpose: Write tsv.
# Inputs: required: x, path, project_root.
# Output: Writes validated artifact file(s) and returns their path(s) invisibly or as a named path list.
write_tsv <- function(x, path, project_root) {
  assert_path_within(project_root, path)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  invisible(path)
}

# Purpose: Require package.
# Inputs: required: package.
# Output: Returns invisibly after validation and stops if the required dependency or argument is unavailable.
require_package <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) stop(sprintf("Required R package '%s' is unavailable.", package), call. = FALSE)
}

# Purpose: Read extended qc config.
# Inputs: required: path.
# Output: Returns parsed, validated R data (vector, data frame, sparse matrix bundle, or named list according to the input format).
read_extended_qc_config <- function(path) {
  if (!file.exists(path)) stop(sprintf("Extended QC config not found: %s", path), call. = FALSE)
  config <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  required <- c("key", "value", "type")
  missing <- setdiff(required, names(config))
  if (length(missing)) stop(sprintf("Extended QC config missing columns: %s", paste(missing, collapse = ", ")), call. = FALSE)
  if (anyDuplicated(config$key)) stop("Extended QC config contains duplicate keys.", call. = FALSE)
  if (anyNA(config[, required, drop = FALSE]) || any(trimws(as.matrix(config[, required, drop = FALSE])) == "")) {
    stop("Extended QC config contains missing or blank required values.", call. = FALSE)
  }
  convert <- function(value, type, key) {
    type <- tolower(type)
    result <- switch(type,
      integer = suppressWarnings(as.integer(value)),
      numeric = suppressWarnings(as.numeric(value)),
      logical = {
        normalized <- toupper(value)
        if (!normalized %in% c("TRUE", "FALSE")) NA else identical(normalized, "TRUE")
      },
      character = as.character(value),
      stop(sprintf("Unknown extended QC config type '%s' for key '%s'.", type, key), call. = FALSE)
    )
    if (length(result) != 1L || is.na(result)) stop(sprintf("Cannot coerce extended QC config key '%s' as %s.", key, type), call. = FALSE)
    result
  }
  values <- Map(convert, config$value, config$type, config$key)
  stats::setNames(values, config$key)
}

# Purpose: Resolve extended qc mode.
# Inputs: required: region_dir; optional/defaulted: requested_mode.
# Output: Returns a deterministic scalar, vector, path, status, or empty-schema object used by downstream functions.
resolve_extended_qc_mode <- function(requested_mode = "AUTO", region_dir) {
  mode <- toupper(trimws(as.character(requested_mode)))
  allowed <- c("AUTO", "LOCAL_SUBSET", "FULL_HPC")
  if (length(mode) != 1L || !mode %in% allowed) {
    stop(sprintf("Extended QC mode must be one of: %s.", paste(allowed, collapse = ", ")), call. = FALSE)
  }
  if (mode != "AUTO") return(mode)
  if (file.exists(file.path(region_dir, "transcripts.parquet"))) "FULL_HPC" else "LOCAL_SUBSET"
}

# Purpose: Extended qc preflight.
# Inputs: required: mode, region_dir, config.
# Output: Returns the derived R object described by the function name; no files are written unless an explicit output path is an input.
extended_qc_preflight <- function(mode, region_dir, config) {
  mode <- resolve_extended_qc_mode(mode, region_dir)
  if (!is.list(config) || !length(config)) stop("Extended QC config must be a non-empty named list.", call. = FALSE)
  checks <- data.frame(
    check = c("transcripts_parquet", "arrow", "dplyr", "RANN"),
    available = c(
      file.exists(file.path(region_dir, "transcripts.parquet")),
      requireNamespace("arrow", quietly = TRUE),
      requireNamespace("dplyr", quietly = TRUE),
      requireNamespace("RANN", quietly = TRUE)
    ),
    stringsAsFactors = FALSE
  )
  checks$required <- mode == "FULL_HPC"
  checks$status <- ifelse(checks$available, "PASS", ifelse(checks$required, "FAIL", "SKIP_ALLOWED"))
  checks$details <- c(
    file.path(region_dir, "transcripts.parquet"),
    "R package for projected Parquet aggregation",
    "R package for lazy Arrow grouping and aggregation",
    "R package for scalable nearest-neighbour calculations"
  )
  checks$mode <- mode
  checks[, c("check", "required", "available", "status", "details", "mode")]
}

# Purpose: Build cycle alarm evidence.
# Inputs: required: alarms, region_id.
# Output: Returns a newly constructed or annotated R object while preserving the supplied raw object/data rows.
build_cycle_alarm_evidence <- function(alarms, region_id) {
  if (!is.data.frame(alarms)) stop("alarms must be a data.frame.", call. = FALSE)
  required <- c("raised", "title", "message", "level", "id")
  missing <- setdiff(required, names(alarms))
  if (length(missing)) stop(sprintf("Alarm table missing columns: %s", paste(missing, collapse = ", ")), call. = FALSE)
  if (!nrow(alarms)) {
    return(data.frame(
      region_id = region_id, raised = FALSE, title = NA_character_, message = NA_character_,
      level = NA_character_, id = NA_character_, evidence_status = "NO_ALARM_REPORTED",
      cycle_identity_status = "NO_POOR_CYCLE_ALARM_REPORTED",
      gene_effect_status = "NO_POOR_CYCLE_ALARM_REPORTED", stringsAsFactors = FALSE
    ))
  }
  evidence <- alarms[, required, drop = FALSE]
  evidence$region_id <- region_id
  evidence$evidence_status <- ifelse(evidence$raised, "DIRECT_EVIDENCE", "ALARM_NOT_RAISED")
  poor_cycle <- evidence$id == "poor_quality_cycles_detected" & evidence$raised
  evidence$cycle_identity_status <- ifelse(poor_cycle, "CYCLE_IDENTITY_UNRESOLVED_REQUIRES_10X", "NOT_APPLICABLE")
  evidence$gene_effect_status <- ifelse(poor_cycle, "GENE_EFFECT_UNCONFIRMED", "NOT_APPLICABLE")
  evidence[, c("region_id", required, "evidence_status", "cycle_identity_status", "gene_effect_status")]
}

# Purpose: Resolve transcript schema.
# Inputs: required: columns.
# Output: Returns a deterministic scalar, vector, path, status, or empty-schema object used by downstream functions.
resolve_transcript_schema <- function(columns) {
  columns <- as.character(columns)
  select_alias <- function(aliases, label, required = TRUE) {
    matched <- intersect(aliases, columns)
    if (length(matched) > 1L) stop(sprintf("Ambiguous transcript %s fields: %s", label, paste(matched, collapse = ", ")), call. = FALSE)
    if (!length(matched)) {
      if (required) stop(sprintf("Transcript schema is missing a %s field.", label), call. = FALSE)
      return(NA_character_)
    }
    matched[[1]]
  }
  list(
    gene = select_alias(c("feature_name", "gene", "target_name"), "gene"),
    qv = select_alias(c("qv", "quality_value"), "QV"),
    codeword = select_alias(c("codeword_index", "codeword"), "codeword", required = FALSE)
  )
}

# Purpose: Summarise transcript quality table.
# Inputs: required: transcripts, region_id; optional/defaulted: qv_threshold.
# Output: Returns computed QC evidence as a vector, data frame, or named summary list; it does not modify raw input files.
summarise_transcript_quality_table <- function(transcripts, region_id, qv_threshold = 20) {
  if (!is.data.frame(transcripts)) stop("transcripts must be a data.frame.", call. = FALSE)
  schema <- resolve_transcript_schema(names(transcripts))
  gene <- as.character(transcripts[[schema$gene]])
  qv <- suppressWarnings(as.numeric(transcripts[[schema$qv]]))
  keep <- !is.na(gene) & nzchar(gene) & is.finite(qv)
  gene <- gene[keep]; qv <- qv[keep]
  codeword <- if (!is.na(schema$codeword)) transcripts[[schema$codeword]][keep] else rep(NA, sum(keep))
  if (!length(gene)) return(data.frame(
    region_id = character(), gene = character(), transcript_rows = integer(), mean_qv = numeric(),
    fraction_q20 = numeric(), represented_codewords = integer(), transcript_status = character(), stringsAsFactors = FALSE
  ))
  groups <- split(seq_along(gene), gene)
  out <- do.call(rbind, lapply(names(groups), function(name) {
    index <- groups[[name]]
    represented <- if (all(is.na(codeword[index]))) NA_integer_ else length(unique(codeword[index][!is.na(codeword[index])]))
    data.frame(
      region_id = region_id, gene = name, transcript_rows = length(index), mean_qv = mean(qv[index]),
      fraction_q20 = mean(qv[index] >= qv_threshold), represented_codewords = represented,
      transcript_status = "MEASURED", stringsAsFactors = FALSE
    )
  }))
  rownames(out) <- NULL
  out[order(out$gene), , drop = FALSE]
}

# Purpose: Build transcript quality queries.
# Inputs: required: projected; optional/defaulted: qv_threshold, has_codeword.
# Output: Returns a newly constructed or annotated R object while preserving the supplied raw object/data rows.
build_transcript_quality_queries <- function(projected, qv_threshold = 20, has_codeword = TRUE) {
  require_package("dplyr")
  required <- c("gene", "qv", if (isTRUE(has_codeword)) "codeword")
  missing <- setdiff(required, names(projected))
  if (length(missing)) {
    stop(sprintf("Projected transcript table missing columns: %s", paste(missing, collapse = ", ")), call. = FALSE)
  }
  filtered <- projected |>
    dplyr::filter(!is.na(.data$gene), .data$gene != "", !is.na(.data$qv))
  summary <- filtered |>
    dplyr::group_by(.data$gene) |>
    dplyr::summarise(
      transcript_rows = dplyr::n(),
      mean_qv = mean(.data$qv),
      fraction_q20 = mean(.data$qv >= qv_threshold),
      .groups = "drop"
    )
  codewords <- if (isTRUE(has_codeword)) {
    filtered |>
      dplyr::filter(!is.na(.data$codeword)) |>
      dplyr::distinct(.data$gene, .data$codeword)
  } else {
    NULL
  }
  list(summary = summary, codewords = codewords)
}

# Purpose: Summarise transcript quality arrow.
# Inputs: required: path, region_id; optional/defaulted: qv_threshold.
# Output: Returns computed QC evidence as a vector, data frame, or named summary list; it does not modify raw input files.
summarise_transcript_quality_arrow <- function(path, region_id, qv_threshold = 20) {
  require_package("arrow")
  require_package("dplyr")
  if (!file.exists(path)) stop(sprintf("Transcript Parquet not found: %s", path), call. = FALSE)
  dataset <- arrow::open_dataset(path, format = "parquet")
  schema <- resolve_transcript_schema(names(dataset$schema))
  if (!is.na(schema$codeword)) {
    projected <- dplyr::select(
      dataset,
      gene = dplyr::all_of(schema$gene),
      qv = dplyr::all_of(schema$qv),
      codeword = dplyr::all_of(schema$codeword)
    )
  } else {
    projected <- dplyr::select(
      dataset,
      gene = dplyr::all_of(schema$gene),
      qv = dplyr::all_of(schema$qv)
    )
  }
  queries <- build_transcript_quality_queries(
    projected,
    qv_threshold = qv_threshold,
    has_codeword = !is.na(schema$codeword)
  )
  out <- dplyr::collect(queries$summary)
  if (!is.na(schema$codeword)) {
    codewords <- dplyr::collect(queries$codewords)
    represented <- if (nrow(codewords)) {
      counts <- table(codewords$gene)
      data.frame(
        gene = names(counts),
        represented_codewords = as.integer(counts),
        stringsAsFactors = FALSE
      )
    } else {
      data.frame(gene = character(), represented_codewords = integer(), stringsAsFactors = FALSE)
    }
    out <- merge(out, represented, by = "gene", all.x = TRUE, sort = FALSE)
    out$represented_codewords[is.na(out$represented_codewords)] <- 0L
  } else {
    out$represented_codewords <- NA_integer_
  }
  out$mean_qv[is.nan(out$mean_qv)] <- NA_real_
  out$fraction_q20[is.nan(out$fraction_q20)] <- NA_real_
  out$region_id <- region_id
  out$transcript_status <- "MEASURED_ARROW_PROJECTED_AGGREGATE"
  out <- out[, c("region_id", "gene", "transcript_rows", "mean_qv", "fraction_q20", "represented_codewords", "transcript_status")]
  out[order(out$gene), , drop = FALSE]
}

# Purpose: Summarise gene matrix qc.
# Inputs: required: counts, region_id; optional/defaulted: gene_sets.
# Output: Returns computed QC evidence as a vector, data frame, or named summary list; it does not modify raw input files.
summarise_gene_matrix_qc <- function(counts, region_id, gene_sets = NULL) {
  require_package("Matrix")
  if (is.null(rownames(counts)) || is.null(colnames(counts))) stop("Gene count matrix requires row and column names.", call. = FALSE)
  raw_counts <- as.numeric(Matrix::rowSums(counts))
  detected_cells <- as.numeric(Matrix::rowSums(counts > 0))
  total <- sum(raw_counts)
  genes <- rownames(counts)
  data.frame(
    region_id = region_id, gene = genes,
    gene_set = if (is.null(gene_sets)) NA_character_ else unname(gene_sets[match(genes, names(gene_sets))]),
    raw_counts = raw_counts,
    counts_per_10000 = if (total > 0) raw_counts / total * 10000 else 0,
    detected_cells = detected_cells,
    detection_fraction = if (ncol(counts) > 0) detected_cells / ncol(counts) else NA_real_,
    matrix_cells = ncol(counts), stringsAsFactors = FALSE
  )
}

# Purpose: Combine gene quality.
# Inputs: required: matrix_qc, transcript_qc.
# Output: Returns computed QC evidence as a vector, data frame, or named summary list; it does not modify raw input files.
combine_gene_quality <- function(matrix_qc, transcript_qc) {
  required_matrix <- c("region_id", "gene")
  if (!all(required_matrix %in% names(matrix_qc))) stop("matrix_qc requires region_id and gene.", call. = FALSE)
  if (!all(required_matrix %in% names(transcript_qc))) stop("transcript_qc requires region_id and gene.", call. = FALSE)
  out <- merge(matrix_qc, transcript_qc, by = c("region_id", "gene"), all.x = TRUE, sort = FALSE)
  out <- out[match(paste(matrix_qc$region_id, matrix_qc$gene), paste(out$region_id, out$gene)), , drop = FALSE]
  out$transcript_rows[is.na(out$transcript_rows)] <- 0L
  out$transcript_status[is.na(out$transcript_status)] <- "NO_TRANSCRIPTS"
  rownames(out) <- NULL
  out
}

# Purpose: Validate spatial cells.
# Inputs: required: cells.
# Output: Returns validation evidence/TRUE (or the validated value) and stops with an informative error when the contract fails.
validate_spatial_cells <- function(cells) {
  required <- c("cell_id", "x_centroid", "y_centroid", "qc_review_flag")
  missing <- setdiff(required, names(cells))
  if (length(missing)) stop(sprintf("Spatial cell table missing columns: %s", paste(missing, collapse = ", ")), call. = FALSE)
  if (anyDuplicated(cells$cell_id)) stop("Spatial cell IDs must be unique.", call. = FALSE)
  if (any(!is.finite(cells$x_centroid)) || any(!is.finite(cells$y_centroid))) stop("Spatial centroids must be finite.", call. = FALSE)
  invisible(TRUE)
}

# Purpose: Knn index distance.
# Inputs: required: cells; optional/defaulted: k, mode.
# Output: Returns the derived R object described by the function name; no files are written unless an explicit output path is an input.
knn_index_distance <- function(cells, k = 15L, mode = "LOCAL_SUBSET") {
  validate_spatial_cells(cells)
  n <- nrow(cells); k <- as.integer(k); mode <- toupper(mode)
  if (n < 2L || k < 1L || k >= n) stop("k must be at least 1 and smaller than the number of cells.", call. = FALSE)
  coordinates <- as.matrix(cells[, c("x_centroid", "y_centroid")])
  if (mode == "FULL_HPC") {
    require_package("RANN")
    result <- RANN::nn2(coordinates, k = k + 1L)
    return(list(index = result$nn.idx[, -1L, drop = FALSE], distance = result$nn.dists[, -1L, drop = FALSE]))
  }
  if (n > 2000L) stop("Base-R nearest-neighbour fallback is limited to 2,000 cells; use FULL_HPC with RANN.", call. = FALSE)
  index <- matrix(NA_integer_, nrow = n, ncol = k)
  distance <- matrix(NA_real_, nrow = n, ncol = k)
  for (i in seq_len(n)) {
    squared <- rowSums((coordinates - matrix(coordinates[i, ], nrow = n, ncol = 2L, byrow = TRUE))^2)
    squared[i] <- Inf
    selected <- order(squared)[seq_len(k)]
    index[i, ] <- selected
    distance[i, ] <- sqrt(squared[selected])
  }
  list(index = index, distance = distance)
}

# Purpose: Calculate knn density.
# Inputs: required: cells; optional/defaulted: k, mode.
# Output: Returns computed QC evidence as a vector, data frame, or named summary list; it does not modify raw input files.
calculate_knn_density <- function(cells, k = 15L, mode = "LOCAL_SUBSET") {
  neighbors <- knn_index_distance(cells, k, mode)
  radius <- neighbors$distance[, ncol(neighbors$distance)]
  radius[radius <= 0] <- min(radius[radius > 0], na.rm = TRUE)
  as.numeric(k) / (pi * radius^2)
}

# Purpose: Assign spatial grid.
# Inputs: required: cells; optional/defaulted: grid_size_um.
# Output: Returns a newly constructed or annotated R object while preserving the supplied raw object/data rows.
assign_spatial_grid <- function(cells, grid_size_um = 100) {
  validate_spatial_cells(cells)
  if (length(grid_size_um) != 1L || !is.finite(grid_size_um) || grid_size_um <= 0) stop("grid_size_um must be positive.", call. = FALSE)
  out <- cells
  origin_x <- min(out$x_centroid); origin_y <- min(out$y_centroid)
  out$grid_x <- as.integer(floor((out$x_centroid - origin_x) / grid_size_um))
  out$grid_y <- as.integer(floor((out$y_centroid - origin_y) / grid_size_um))
  out$grid_id <- paste(out$grid_x, out$grid_y, sep = ":")
  occupied <- unique(out$grid_id)
  unique_bins <- unique(out[, c("grid_x", "grid_y", "grid_id")])
  unique_bins$edge_proxy <- vapply(seq_len(nrow(unique_bins)), function(i) {
    offsets <- expand.grid(dx = -1:1, dy = -1:1)
    offsets <- offsets[!(offsets$dx == 0 & offsets$dy == 0), , drop = FALSE]
    neighbors <- paste(unique_bins$grid_x[i] + offsets$dx, unique_bins$grid_y[i] + offsets$dy, sep = ":")
    any(!neighbors %in% occupied)
  }, logical(1))
  out$edge_proxy <- unique_bins$edge_proxy[match(out$grid_id, unique_bins$grid_id)]
  out$grid_size_um <- grid_size_um
  out
}

# Purpose: Summarise spatial enrichment.
# Inputs: required: annotated_cells.
# Output: Returns computed QC evidence as a vector, data frame, or named summary list; it does not modify raw input files.
summarise_spatial_enrichment <- function(annotated_cells) {
  if (!all(c("qc_review_flag", "edge_proxy") %in% names(annotated_cells))) stop("Spatial enrichment requires qc_review_flag and edge_proxy.", call. = FALSE)
  make_rows <- function(class_type, positive, positive_label, negative_label) {
    rows <- do.call(rbind, lapply(list(positive, !positive), function(index) {
      label <- if (identical(index, positive)) positive_label else negative_label
      cells <- sum(index); flagged <- sum(annotated_cells$qc_review_flag[index])
      data.frame(class_type = class_type, class = label, cells = cells, flagged = flagged,
                 review_rate = if (cells) flagged / cells else NA_real_, stringsAsFactors = FALSE)
    }))
    if (any(rows$cells == 0L)) {
      rows$risk_ratio <- NA_real_
      rows$absolute_rate_difference <- NA_real_
    } else {
      corrected_rate <- (rows$flagged + 0.5) / (rows$cells + 1)
      rows$risk_ratio <- corrected_rate[[1]] / corrected_rate[[2]]
      rows$absolute_rate_difference <- rows$review_rate[[1]] - rows$review_rate[[2]]
    }
    rows
  }
  out <- make_rows("edge_proxy", annotated_cells$edge_proxy, "edge", "interior")
  if ("dense_aggregate" %in% names(annotated_cells)) {
    out <- rbind(out, make_rows("local_density", annotated_cells$dense_aggregate, "dense", "non_dense"))
  }
  rownames(out) <- NULL
  out
}

# Purpose: Test spatial flag clustering.
# Inputs: required: cells; optional/defaulted: k, permutations, seed, mode.
# Output: Returns computed QC evidence as a vector, data frame, or named summary list; it does not modify raw input files.
test_spatial_flag_clustering <- function(cells, k = 15L, permutations = 999L, seed = 20260814L, mode = "LOCAL_SUBSET") {
  validate_spatial_cells(cells)
  flags <- as.numeric(as.logical(cells$qc_review_flag))
  if (length(unique(flags)) < 2L) return(data.frame(
    status = "NOT_ESTIMABLE", statistic = NA_real_, empirical_p = NA_real_, permutations = as.integer(permutations),
    k = as.integer(k), seed = as.integer(seed), stringsAsFactors = FALSE
  ))
  neighbors <- knn_index_distance(cells, k, mode)$index
  statistic <- function(values) {
    centered <- values - mean(values)
    denominator <- sum(centered^2)
    if (denominator == 0) return(NA_real_)
    nrow(neighbors) / length(neighbors) * sum(centered[row(neighbors)] * centered[neighbors]) / denominator
  }
  observed <- statistic(flags)
  set.seed(as.integer(seed))
  permuted <- replicate(as.integer(permutations), statistic(sample(flags, replace = FALSE)))
  data.frame(
    status = "ESTIMATED", statistic = observed,
    empirical_p = (1 + sum(permuted >= observed, na.rm = TRUE)) / (1 + as.integer(permutations)),
    permutations = as.integer(permutations), k = as.integer(k), seed = as.integer(seed), stringsAsFactors = FALSE
  )
}

# Purpose: Find spatial qc hotspots.
# Inputs: required: annotated_cells; optional/defaulted: permutations, min_bin_cells, fdr, seed.
# Output: Returns computed QC evidence as a vector, data frame, or named summary list; it does not modify raw input files.
find_spatial_qc_hotspots <- function(annotated_cells, permutations = 999L, min_bin_cells = 20L, fdr = 0.05, seed = 20260814L) {
  required <- c("grid_id", "grid_x", "grid_y", "x_centroid", "y_centroid", "qc_review_flag", "grid_size_um")
  missing <- setdiff(required, names(annotated_cells))
  if (length(missing)) stop(sprintf("Hotspot table missing columns: %s", paste(missing, collapse = ", ")), call. = FALSE)
  bin_factor <- factor(annotated_cells$grid_id, levels = unique(annotated_cells$grid_id))
  cells_per_bin <- as.integer(table(bin_factor))
  flags <- as.integer(as.logical(annotated_cells$qc_review_flag))
  flagged_per_bin <- as.integer(rowsum(flags, bin_factor, reorder = FALSE))
  eligible <- cells_per_bin >= as.integer(min_bin_cells)
  bins <- levels(bin_factor)[eligible]
  if (!length(bins)) return(data.frame(
    grid_id = character(), cells = integer(), flagged = integer(), review_rate = numeric(), global_rate = numeric(),
    rate_difference = numeric(), empirical_p = numeric(), adjusted_p = numeric(), hotspot_status = character(),
    x_min = numeric(), x_max = numeric(), y_min = numeric(), y_max = numeric(),
    fdr_threshold = numeric(), min_bin_cells_threshold = integer(), permutations = integer(), seed = integer(),
    stringsAsFactors = FALSE
  ))
  observed <- flagged_per_bin[eligible]
  denominators <- cells_per_bin[eligible]
  global_rate <- mean(flags)
  set.seed(as.integer(seed))
  exceedances <- integer(length(bins))
  for (iteration in seq_len(as.integer(permutations))) {
    permuted_counts <- as.integer(rowsum(sample(flags, replace = FALSE), bin_factor, reorder = FALSE))[eligible]
    exceedances <- exceedances + as.integer(permuted_counts >= observed)
  }
  empirical_p <- (1 + exceedances) / (1 + as.integer(permutations))
  adjusted_p <- stats::p.adjust(empirical_p, method = "BH")
  bounds <- do.call(rbind, lapply(bins, function(id) {
    index <- annotated_cells$grid_id == id
    data.frame(x_min = min(annotated_cells$x_centroid[index]), x_max = max(annotated_cells$x_centroid[index]),
               y_min = min(annotated_cells$y_centroid[index]), y_max = max(annotated_cells$y_centroid[index]))
  }))
  data.frame(
    grid_id = bins, cells = denominators, flagged = observed, review_rate = observed / denominators,
    global_rate = global_rate, rate_difference = observed / denominators - global_rate,
    empirical_p = empirical_p, adjusted_p = adjusted_p,
    hotspot_status = ifelse(adjusted_p <= fdr & observed / denominators > global_rate, "MORPHOLOGY_REVIEW_REQUIRED", "NO_HOTSPOT_EVIDENCE"),
    bounds, fdr_threshold = as.numeric(fdr), min_bin_cells_threshold = as.integer(min_bin_cells),
    permutations = as.integer(permutations), seed = as.integer(seed), stringsAsFactors = FALSE
  )
}

# Purpose: Rank candidate cycle genes.
# Inputs: required: gene_quality, config; optional/defaulted: reference_region.
# Output: Returns computed QC evidence as a vector, data frame, or named summary list; it does not modify raw input files.
rank_candidate_cycle_genes <- function(gene_quality, config, reference_region = "Region_3") {
  required <- c("region_id", "gene", "counts_per_10000", "detection_fraction", "fraction_q20")
  missing <- setdiff(required, names(gene_quality))
  if (length(missing)) stop(sprintf("Gene-quality table missing columns: %s", paste(missing, collapse = ", ")), call. = FALSE)
  if (anyDuplicated(gene_quality[, c("region_id", "gene")])) stop("Gene-quality table contains duplicate region/gene rows.", call. = FALSE)
  genes <- unique(gene_quality$gene)
  reference <- gene_quality[gene_quality$region_id == reference_region, , drop = FALSE]
  if (!setequal(reference$gene, genes)) stop("Reference region must contain exactly one row for every gene.", call. = FALSE)
  comparison_regions <- setdiff(unique(gene_quality$region_id), reference_region)
  out <- do.call(rbind, lapply(genes, function(gene) {
    ref <- reference[reference$gene == gene, , drop = FALSE]
    observed <- gene_quality[gene_quality$gene == gene & gene_quality$region_id %in% comparison_regions, , drop = FALSE]
    data.frame(
      gene = gene, region_id = observed$region_id, reference_region = reference_region,
      comparison_type = ifelse(observed$region_id == "Region_4", "WITHIN_MOUSE_TECHNICAL_PAIR", "CROSS_MOUSE_DIAGNOSTIC_REFERENCE"),
      observed_counts_per_10000 = observed$counts_per_10000,
      reference_counts_per_10000 = ref$counts_per_10000,
      log2_count_ratio = log2((observed$counts_per_10000 + 0.5) / (ref$counts_per_10000 + 0.5)),
      detection_fraction_difference = observed$detection_fraction - ref$detection_fraction,
      observed_fraction_q20 = observed$fraction_q20,
      reference_fraction_q20 = ref$fraction_q20,
      q20_difference = observed$fraction_q20 - ref$fraction_q20,
      stringsAsFactors = FALSE
    )
  }))
  out$abundance_depletion <- is.finite(out$log2_count_ratio) & out$log2_count_ratio <= -as.numeric(config$candidate_abs_log2_ratio)
  out$quality_degradation <- is.finite(out$q20_difference) & out$q20_difference <= -as.numeric(config$candidate_abs_q20_delta)
  out$insufficient_quality <- !is.finite(out$q20_difference)
  out$section_candidate_flag <- out$abundance_depletion | out$quality_degradation
  out$section_evidence_status <- ifelse(
    out$abundance_depletion & out$quality_degradation,
    "DEPLETION_AND_Q20_LOSS",
    ifelse(
      out$abundance_depletion,
      "DEPLETION_ONLY",
      ifelse(
        out$quality_degradation,
        "Q20_LOSS_ONLY",
        ifelse(out$insufficient_quality, "INSUFFICIENT_Q20_EVIDENCE", "NO_SECTION_LEVEL_SIGNAL")
      )
    )
  )
  tiers <- vapply(split(seq_len(nrow(out)), out$gene), function(index) {
    paired <- any(out$region_id[index] == "Region_4" & out$abundance_depletion[index] & out$quality_degradation[index])
    recurring <- sum(out$abundance_depletion[index] & out$quality_degradation[index]) >= 2L
    partial <- any(xor(out$abundance_depletion[index], out$quality_degradation[index]))
    if (paired) "Tier_A" else if (recurring) "Tier_B" else if (partial) "Tier_C" else "Unranked"
  }, character(1))
  out$evidence_tier <- unname(tiers[match(out$gene, names(tiers))])
  out$candidate_status <- "CANDIDATE_NOT_CONFIRMED"
  out$exact_cycle_status <- "REQUIRES_10X_DIAGNOSTICS"
  out$threshold_abs_log2_ratio <- as.numeric(config$candidate_abs_log2_ratio)
  out$threshold_abs_q20_delta <- as.numeric(config$candidate_abs_q20_delta)
  rownames(out) <- NULL
  out[order(out$gene, out$region_id), , drop = FALSE]
}

# Purpose: Build gene downstream decision.
# Inputs: required: candidates, panel_genes, run_label, execution_mode, provenance.
# Output: Returns a newly constructed or annotated R object while preserving the supplied raw object/data rows.
build_gene_downstream_decision <- function(candidates, panel_genes, run_label,
                                           execution_mode, provenance) {
  required <- c(
    "gene", "region_id", "section_candidate_flag", "section_evidence_status",
    "threshold_abs_log2_ratio", "threshold_abs_q20_delta"
  )
  missing <- setdiff(required, names(candidates))
  if (length(missing)) {
    stop(sprintf("Candidate evidence lacks gene-decision fields: %s", paste(missing, collapse = ", ")), call. = FALSE)
  }
  if (is.data.frame(panel_genes)) {
    gene_column <- intersect(c("gene", "feature_name", "name"), names(panel_genes))
    if (!length(gene_column)) stop("Panel features require a gene column.", call. = FALSE)
    panel_genes <- panel_genes[[gene_column[[1]]]]
  }
  panel_genes <- unique(as.character(panel_genes))
  panel_genes <- panel_genes[!is.na(panel_genes) & nzchar(panel_genes)]
  if (!length(panel_genes)) stop("Panel genes cannot be empty.", call. = FALSE)
  if (!all(panel_genes %in% candidates$gene)) {
    stop("Candidate evidence must contain every panel gene.", call. = FALSE)
  }
  alarm_regions <- c("Region_1", "Region_2", "Region_4")
  candidate_alarm <- candidates[candidates$region_id %in% alarm_regions, , drop = FALSE]
  recurrence <- vapply(panel_genes, function(gene) {
    rows <- candidate_alarm[candidate_alarm$gene == gene, , drop = FALSE]
    length(unique(rows$region_id[as.logical(rows$section_candidate_flag)]))
  }, integer(1))
  flag_for <- function(gene, region) {
    rows <- candidate_alarm[candidate_alarm$gene == gene & candidate_alarm$region_id == region, , drop = FALSE]
    if (!nrow(rows)) return(NA)
    any(as.logical(rows$section_candidate_flag))
  }
  evidence_for <- function(gene, region) {
    rows <- candidate_alarm[candidate_alarm$gene == gene & candidate_alarm$region_id == region, , drop = FALSE]
    if (!nrow(rows)) return(NA_character_)
    paste(unique(as.character(rows$section_evidence_status)), collapse = ";")
  }
  out <- data.frame(
    gene = panel_genes,
    raw_panel_status = "RAW_COMPLETE_PANEL",
    alarm_positive_section_count = recurrence,
    conservative_evidence_status = ifelse(
      recurrence == 0L, "CONSERVATIVE_NO_SIGNAL_DETECTED", "NOT_IN_CONSERVATIVE_ZERO_ALARM_SET"
    ),
    primary_feature_status = ifelse(
      recurrence <= 1L, "PROVISIONAL_PRIMARY_FEATURES", "EXCLUDED_FROM_PRIMARY_FEATURES"
    ),
    technical_risk_status = ifelse(
      recurrence >= 2L, "TECHNICAL_RISK_SENSITIVITY_ONLY", "NOT_IN_TECHNICAL_RISK_SET"
    ),
    region_1_flag = vapply(panel_genes, flag_for, logical(1), region = "Region_1"),
    region_2_flag = vapply(panel_genes, flag_for, logical(1), region = "Region_2"),
    region_4_flag = vapply(panel_genes, flag_for, logical(1), region = "Region_4"),
    region_1_evidence = vapply(panel_genes, evidence_for, character(1), region = "Region_1"),
    region_2_evidence = vapply(panel_genes, evidence_for, character(1), region = "Region_2"),
    region_4_evidence = vapply(panel_genes, evidence_for, character(1), region = "Region_4"),
    threshold_abs_log2_ratio = unique(candidates$threshold_abs_log2_ratio)[[1]],
    threshold_abs_q20_delta = unique(candidates$threshold_abs_q20_delta)[[1]],
    expected_full_raw_count = 479L,
    expected_full_conservative_count = 67L,
    expected_full_provisional_count = 245L,
    expected_full_technical_risk_count = 234L,
    decision_rule = "0 alarm sections=conservative subset; 0-1=provisional primary; >=2=technical-risk sensitivity-only",
    cycle_mapping_status = "UNAVAILABLE_EVIDENCE_ONLY_DECISION",
    run_label = run_label, execution_mode = toupper(execution_mode),
    provenance = provenance, stringsAsFactors = FALSE
  )
  out[order(out$gene), , drop = FALSE]
}

# Purpose: Build eos gene decision.
# Inputs: required: gene_decision, eos_gene_sets, run_label, execution_mode, provenance.
# Output: Returns a newly constructed or annotated R object while preserving the supplied raw object/data rows.
build_eos_gene_decision <- function(gene_decision, eos_gene_sets, run_label,
                                    execution_mode, provenance) {
  required_gene <- c("gene", "primary_feature_status", "conservative_evidence_status", "technical_risk_status")
  required_eos <- c("gene", "gene_set")
  if (length(setdiff(required_gene, names(gene_decision)))) stop("Gene decision lacks Eos decision fields.", call. = FALSE)
  if (length(setdiff(required_eos, names(eos_gene_sets)))) stop("Eos gene set requires gene and gene_set.", call. = FALSE)
  if (anyDuplicated(eos_gene_sets$gene)) stop("Eos gene set contains duplicated genes.", call. = FALSE)
  matched <- match(as.character(eos_gene_sets$gene), gene_decision$gene)
  out <- eos_gene_sets[, required_eos, drop = FALSE]
  out$panel_membership <- !is.na(matched)
  out$retained_provisional <- !is.na(matched) &
    gene_decision$primary_feature_status[matched] == "PROVISIONAL_PRIMARY_FEATURES"
  out$conservative_sensitivity <- !is.na(matched) &
    gene_decision$conservative_evidence_status[matched] == "CONSERVATIVE_NO_SIGNAL_DETECTED"
  out$technical_risk <- !is.na(matched) &
    gene_decision$technical_risk_status[matched] == "TECHNICAL_RISK_SENSITIVITY_ONLY"
  out$provisional_signature_status <- ifelse(
    out$retained_provisional, "EOS_PROVISIONAL_PRIMARY_53", "NOT_IN_EOS_PROVISIONAL_PRIMARY"
  )
  out$complete_signature_status <- "RAW_COMPLETE_EOS_100"
  expected_by_set <- c(common = 4L, short_lived = 27L, long_lived = 22L)
  out$expected_full_retained_total <- 53L
  out$expected_full_gene_set_count <- unname(expected_by_set[as.character(out$gene_set)])
  out$selection_threshold <- "alarm_positive_section_count <= 1"
  out$decision_rule <- "Retain Eos genes eligible as PROVISIONAL_PRIMARY_FEATURES; compare with complete Eos signature in Region_3"
  out$run_label <- run_label
  out$execution_mode <- toupper(execution_mode)
  out$provenance <- provenance
  out
}

# Purpose: Compare subset full qc.
# Inputs: required: full_summary, subset_reference.
# Output: Returns computed QC evidence as a vector, data frame, or named summary list; it does not modify raw input files.
compare_subset_full_qc <- function(full_summary, subset_reference) {
  full_required <- c("region_id", "input_cells", "review_flagged")
  subset_required <- c("region_id", "review_fraction", "subset_rank")
  if (length(setdiff(full_required, names(full_summary)))) stop("Full QC summary lacks required review columns.", call. = FALSE)
  if (length(setdiff(subset_required, names(subset_reference)))) stop("Subset reference lacks required rank columns.", call. = FALSE)
  if (anyDuplicated(full_summary$region_id) || anyDuplicated(subset_reference$region_id)) stop("QC ranking inputs require unique regions.", call. = FALSE)
  if (!setequal(full_summary$region_id, subset_reference$region_id)) stop("Full and subset QC region sets differ.", call. = FALSE)
  full <- full_summary[, full_required, drop = FALSE]
  full$full_review_fraction <- ifelse(full$input_cells > 0, full$review_flagged / full$input_cells, NA_real_)
  order_index <- order(-full$full_review_fraction, full$region_id)
  full$full_rank <- NA_integer_; full$full_rank[order_index] <- seq_along(order_index)
  ranking <- merge(
    subset_reference[, subset_required, drop = FALSE], full,
    by = "region_id", all = FALSE, sort = FALSE
  )
  ranking <- ranking[match(full_summary$region_id, ranking$region_id), , drop = FALSE]
  ranking$rank_change <- ranking$full_rank - ranking$subset_rank
  ranking$review_fraction_difference <- ranking$full_review_fraction - ranking$review_fraction
  agreement <- data.frame(
    sections = nrow(ranking),
    spearman_rho = suppressWarnings(stats::cor(ranking$subset_rank, ranking$full_rank, method = "spearman")),
    kendall_tau = suppressWarnings(stats::cor(ranking$subset_rank, ranking$full_rank, method = "kendall")),
    interpretation = "DESCRIPTIVE_FOUR_SECTIONS", stringsAsFactors = FALSE
  )
  list(ranking = ranking, agreement = agreement)
}

# Purpose: Build section pairs.
# Inputs: required: manifest.
# Output: Returns a newly constructed or annotated R object while preserving the supplied raw object/data rows.
build_section_pairs <- function(manifest) {
  required <- c("region_id", "mouse_id", "section_id")
  missing <- setdiff(required, names(manifest))
  if (length(missing)) stop(sprintf("Manifest missing pair columns: %s", paste(missing, collapse = ", ")), call. = FALSE)
  if (anyDuplicated(manifest$region_id) || anyDuplicated(manifest$section_id)) stop("Manifest section and region IDs must be unique.", call. = FALSE)
  rows <- lapply(unique(as.character(manifest$mouse_id)), function(mouse) {
    sample_rows <- manifest[manifest$mouse_id == mouse, , drop = FALSE]
    if (nrow(sample_rows) != 2L) return(data.frame(
      mouse_id = mouse, section_a = NA_character_, section_b = NA_character_,
      region_a = NA_character_, region_b = NA_character_, pair_status = "NOT_ESTIMABLE",
      details = sprintf("Expected two sections but found %d", nrow(sample_rows)), stringsAsFactors = FALSE
    ))
    region_number <- suppressWarnings(as.integer(sub("^Region_", "", sample_rows$region_id)))
    sample_rows <- sample_rows[order(region_number), , drop = FALSE]
    data.frame(
      mouse_id = mouse, section_a = as.character(sample_rows$section_id[[1]]), section_b = as.character(sample_rows$section_id[[2]]),
      region_a = as.character(sample_rows$region_id[[1]]), region_b = as.character(sample_rows$region_id[[2]]),
      pair_status = "ESTIMABLE", details = "Two verified technical sections", stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows); rownames(out) <- NULL; out
}

# Purpose: Quantile distribution distance.
# Inputs: required: a, b; optional/defaulted: probabilities.
# Output: Returns the derived R object described by the function name; no files are written unless an explicit output path is an input.
quantile_distribution_distance <- function(a, b, probabilities = seq(0.01, 0.99, 0.01)) {
  a <- a[is.finite(a)]; b <- b[is.finite(b)]
  if (!length(a) || !length(b)) return(NA_real_)
  qa <- as.numeric(stats::quantile(a, probabilities, names = FALSE, type = 7))
  qb <- as.numeric(stats::quantile(b, probabilities, names = FALSE, type = 7))
  scale <- stats::median(c(a, b))
  if (!is.finite(scale) || scale == 0) scale <- 1
  mean(abs(qa - qb)) / abs(scale)
}

# Purpose: Calculate within mouse concordance.
# Inputs: required: manifest, section_summary, cell_metadata, gene_quality, config.
# Output: Returns computed QC evidence as a vector, data frame, or named summary list; it does not modify raw input files.
calculate_within_mouse_concordance <- function(manifest, section_summary, cell_metadata, gene_quality, config) {
  pairs <- build_section_pairs(manifest)
  summary_required <- c("region_id", "input_cells", "review_flagged")
  cell_required <- c("region_id", "nCount_Xenium", "nFeature_Xenium", "cell_area", "control_fraction")
  gene_required <- c("region_id", "gene", "counts_per_10000", "detection_fraction")
  if (length(setdiff(summary_required, names(section_summary)))) stop("Section summary lacks concordance columns.", call. = FALSE)
  if (length(setdiff(cell_required, names(cell_metadata)))) stop("Cell metadata lacks concordance columns.", call. = FALSE)
  if (length(setdiff(gene_required, names(gene_quality)))) stop("Gene-quality table lacks concordance columns.", call. = FALSE)
  summary_rows <- list(); gene_rows <- list()
  for (index in seq_len(nrow(pairs))) {
    pair <- pairs[index, , drop = FALSE]
    if (pair$pair_status != "ESTIMABLE") {
      summary_rows[[index]] <- data.frame(
        mouse_id = pair$mouse_id, section_a = pair$section_a, section_b = pair$section_b,
        region_a = pair$region_a, region_b = pair$region_b, pair_status = pair$pair_status,
        review_rate_a = NA_real_, review_rate_b = NA_real_, review_rate_difference = NA_real_,
        median_count_ratio = NA_real_, median_feature_ratio = NA_real_, median_area_ratio = NA_real_,
        median_control_fraction_ratio = NA_real_, count_quantile_distance = NA_real_,
        feature_quantile_distance = NA_real_, area_quantile_distance = NA_real_, control_quantile_distance = NA_real_,
        gene_count_spearman = NA_real_, gene_detection_spearman = NA_real_,
        review_criterion = NA, count_criterion = NA, feature_criterion = NA,
        gene_count_criterion = NA, gene_detection_criterion = NA,
        concordance_status = "NOT_ESTIMABLE",
        interpretation = "ADVISORY_TECHNICAL_CONCORDANCE_NOT_BIOLOGICAL_TEST", stringsAsFactors = FALSE
      )
      next
    }
    qa <- section_summary[section_summary$region_id == pair$region_a, , drop = FALSE]
    qb <- section_summary[section_summary$region_id == pair$region_b, , drop = FALSE]
    ca <- cell_metadata[cell_metadata$region_id == pair$region_a, , drop = FALSE]
    cb <- cell_metadata[cell_metadata$region_id == pair$region_b, , drop = FALSE]
    if (nrow(qa) != 1L || nrow(qb) != 1L || !nrow(ca) || !nrow(cb)) stop(sprintf("Incomplete concordance inputs for %s.", pair$mouse_id), call. = FALSE)
    ga <- gene_quality[gene_quality$region_id == pair$region_a, gene_required[-1], drop = FALSE]
    gb <- gene_quality[gene_quality$region_id == pair$region_b, gene_required[-1], drop = FALSE]
    genes <- merge(ga, gb, by = "gene", suffixes = c("_a", "_b"), all = TRUE, sort = TRUE)
    genes$mouse_id <- pair$mouse_id; genes$region_a <- pair$region_a; genes$region_b <- pair$region_b
    complete_counts <- is.finite(genes$counts_per_10000_a) & is.finite(genes$counts_per_10000_b)
    complete_detection <- is.finite(genes$detection_fraction_a) & is.finite(genes$detection_fraction_b)
    count_correlation <- if (sum(complete_counts) >= 20L) stats::cor(log1p(genes$counts_per_10000_a[complete_counts]), log1p(genes$counts_per_10000_b[complete_counts]), method = "spearman") else NA_real_
    detection_correlation <- if (sum(complete_detection) >= 20L) stats::cor(genes$detection_fraction_a[complete_detection], genes$detection_fraction_b[complete_detection], method = "spearman") else NA_real_
    review_a <- qa$review_flagged / qa$input_cells; review_b <- qb$review_flagged / qb$input_cells
    median_ratio <- function(name) stats::median(cb[[name]], na.rm = TRUE) / stats::median(ca[[name]], na.rm = TRUE)
    count_ratio <- median_ratio("nCount_Xenium"); feature_ratio <- median_ratio("nFeature_Xenium")
    criteria <- c(
      review = abs(review_b - review_a) <= as.numeric(config$concordance_review_rate_difference),
      counts = count_ratio >= as.numeric(config$concordance_count_ratio_lower) & count_ratio <= as.numeric(config$concordance_count_ratio_upper),
      features = feature_ratio >= as.numeric(config$concordance_feature_ratio_lower) & feature_ratio <= as.numeric(config$concordance_feature_ratio_upper),
      gene_counts = is.finite(count_correlation) && count_correlation >= as.numeric(config$concordance_gene_spearman),
      gene_detection = is.finite(detection_correlation) && detection_correlation >= as.numeric(config$concordance_gene_spearman)
    )
    summary_rows[[index]] <- data.frame(
      mouse_id = pair$mouse_id, section_a = pair$section_a, section_b = pair$section_b,
      region_a = pair$region_a, region_b = pair$region_b, pair_status = pair$pair_status,
      review_rate_a = review_a, review_rate_b = review_b, review_rate_difference = abs(review_b - review_a),
      median_count_ratio = count_ratio, median_feature_ratio = feature_ratio,
      median_area_ratio = median_ratio("cell_area"), median_control_fraction_ratio = median_ratio("control_fraction"),
      count_quantile_distance = quantile_distribution_distance(ca$nCount_Xenium, cb$nCount_Xenium),
      feature_quantile_distance = quantile_distribution_distance(ca$nFeature_Xenium, cb$nFeature_Xenium),
      area_quantile_distance = quantile_distribution_distance(ca$cell_area, cb$cell_area),
      control_quantile_distance = quantile_distribution_distance(ca$control_fraction, cb$control_fraction),
      gene_count_spearman = count_correlation, gene_detection_spearman = detection_correlation,
      review_criterion = criteria[["review"]], count_criterion = criteria[["counts"]],
      feature_criterion = criteria[["features"]], gene_count_criterion = criteria[["gene_counts"]],
      gene_detection_criterion = criteria[["gene_detection"]],
      concordance_status = if (all(criteria)) "CONCORDANT" else "REVIEW",
      interpretation = "ADVISORY_TECHNICAL_CONCORDANCE_NOT_BIOLOGICAL_TEST", stringsAsFactors = FALSE
    )
    gene_rows[[index]] <- genes
  }
  summary <- do.call(rbind, summary_rows); rownames(summary) <- NULL
  genes <- if (length(gene_rows)) do.call(rbind, gene_rows) else data.frame()
  if (nrow(genes)) rownames(genes) <- NULL
  list(summary = summary, genes = genes)
}

# Purpose: Empty alarm table.
# Inputs: none.
# Output: Returns a deterministic scalar, vector, path, status, or empty-schema object used by downstream functions.
empty_alarm_table <- function() {
  data.frame(
    raw_value = logical(), formatted_value = character(), raised = logical(),
    title = character(), message = character(), level = character(), id = character(),
    stringsAsFactors = FALSE
  )
}

# Purpose: Extract analysis alarms.
# Inputs: required: path.
# Output: Returns computed QC evidence as a vector, data frame, or named summary list; it does not modify raw input files.
extract_analysis_alarms <- function(path) {
  require_package("jsonlite")
  if (!file.exists(path)) stop(sprintf("Analysis summary not found: %s", path), call. = FALSE)
  html <- readChar(path, nchars = file.info(path)$size, useBytes = TRUE)
  prefix <- '"alarms":{"alarms":'
  start <- regexpr(prefix, html, fixed = TRUE)[1]
  if (start < 0L) return(empty_alarm_table())
  remainder <- substr(html, start + nchar(prefix, type = "bytes"), nchar(html))
  finish <- regexpr('},"sample":', remainder, fixed = TRUE)[1]
  if (finish < 0L) stop(sprintf("Could not parse alarms block in %s", path), call. = FALSE)
  parsed <- jsonlite::fromJSON(substr(remainder, 1L, finish - 1L), simplifyDataFrame = TRUE)
  if (!length(parsed)) return(empty_alarm_table())
  if (!is.data.frame(parsed)) parsed <- as.data.frame(parsed, stringsAsFactors = FALSE)
  wanted <- names(empty_alarm_table())
  for (name in setdiff(wanted, names(parsed))) parsed[[name]] <- NA
  parsed[, wanted, drop = FALSE]
}

# Purpose: Read custom panel genes.
# Inputs: required: path.
# Output: Returns parsed, validated R data (vector, data frame, sparse matrix bundle, or named list according to the input format).
read_custom_panel_genes <- function(path) {
  require_package("jsonlite")
  if (!file.exists(path)) stop(sprintf("Panel JSON not found: %s", path), call. = FALSE)
  panel <- jsonlite::fromJSON(path, simplifyVector = FALSE)
  targets <- panel$payload$targets %||% list()
  genes <- vapply(targets, function(target) {
    if (identical(target$source$category, "current") && identical(target$type$descriptor, "gene")) target$type$data$name else NA_character_
  }, character(1))
  sort(unique(stats::na.omit(genes)))
}

# Purpose: Reconcile panel.
# Inputs: required: expected_genes, installed_genes; optional/defaulted: gene_sets.
# Output: Returns computed QC evidence as a vector, data frame, or named summary list; it does not modify raw input files.
reconcile_panel <- function(expected_genes, installed_genes, gene_sets = NULL) {
  genes <- sort(unique(c(as.character(expected_genes), as.character(installed_genes))))
  out <- data.frame(
    gene = genes, expected = genes %in% expected_genes, installed = genes %in% installed_genes,
    stringsAsFactors = FALSE
  )
  out$status <- ifelse(out$expected & out$installed, "MATCH", ifelse(out$expected, "MISSING", "EXTRA"))
  if (!is.null(gene_sets)) out$gene_set <- gene_sets[match(out$gene, names(gene_sets))]
  out
}

# Purpose: Read xenium features.
# Inputs: required: path.
# Output: Returns parsed, validated R data (vector, data frame, sparse matrix bundle, or named list according to the input format).
read_xenium_features <- function(path) {
  features <- read_gz_rows(path, header = FALSE)
  if (ncol(features) != 3L) stop(sprintf("Expected three columns in %s", path), call. = FALSE)
  names(features) <- c("feature_id", "feature_name", "feature_type")
  features
}

# Purpose: Read xenium barcodes.
# Inputs: required: path.
# Output: Returns parsed, validated R data (vector, data frame, sparse matrix bundle, or named list according to the input format).
read_xenium_barcodes <- function(path) {
  scan(gzfile(path), what = character(), quiet = TRUE)
}

# Purpose: Import xenium mex.
# Inputs: required: region_dir.
# Output: Returns parsed, validated R data (vector, data frame, sparse matrix bundle, or named list according to the input format).
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

# Purpose: Safe quantile.
# Inputs: required: x, probability.
# Output: Returns the derived R object described by the function name; no files are written unless an explicit output path is an input.
safe_quantile <- function(x, probability) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  unname(stats::quantile(x, probability, names = FALSE, na.rm = TRUE, type = 7))
}

# Purpose: Robust interval.
# Inputs: required: x; optional/defaulted: lower_mads, upper_mads, floor_value.
# Output: Returns the derived R object described by the function name; no files are written unless an explicit output path is an input.
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

# Purpose: Overall readiness.
# Inputs: required: status.
# Output: Returns a deterministic scalar, vector, path, status, or empty-schema object used by downstream functions.
overall_readiness <- function(status) {
  status <- toupper(as.character(status))
  if (any(status %in% c("FAIL", "HOLD", "BLOCKED"))) return("HOLD")
  if (any(status %in% c("PENDING", "WARN"))) return("PENDING")
  "PASS"
}

# Purpose: Calculate readiness gates.
# Inputs: required: inventory, integrity, panel_reconciliation, alarms, manifest.
# Output: Returns computed QC evidence as a vector, data frame, or named summary list; it does not modify raw input files.
calculate_readiness_gates <- function(inventory, integrity, panel_reconciliation, alarms, manifest) {
  alarm_levels <- if (nrow(alarms) && "level" %in% names(alarms)) toupper(alarms$level) else character()
  panel_bad <- if (nrow(panel_reconciliation)) sum(panel_reconciliation$status %in% c("MISSING", "EXTRA")) else 0L
  gates <- rbind(
    data.frame(gate = "required_files", status = if (nrow(inventory) && all(inventory$exists)) "PASS" else "HOLD", details = sprintf("%d missing required files", sum(!inventory$exists))),
    data.frame(gate = "matrix_integrity", status = if (nrow(integrity) == 1L && isTRUE(integrity$dimension_match[[1]])) "PASS" else "HOLD", details = "Matrix, features, barcodes, and cells must align"),
    data.frame(gate = "panel_reconciliation", status = if (panel_bad == 0L) "PASS" else "HOLD", details = sprintf("%d missing/extra expected-panel entries", panel_bad)),
    data.frame(gate = "metadata", status = if (any(manifest$metadata_status == "SYNTHETIC_PLACEHOLDER")) "PENDING" else "PASS", details = "Synthetic metadata blocks biological interpretation"),
    data.frame(gate = "xenium_analysis_alerts", status = if (any(alarm_levels == "ERROR")) "HOLD" else if (length(alarm_levels)) "WARN" else "PASS", details = sprintf("%d alarms; %d errors", nrow(alarms), sum(alarm_levels == "ERROR"))),
    stringsAsFactors = FALSE
  )
  gates <- rbind(gates, data.frame(gate = "overall", status = overall_readiness(gates$status), details = "Worst-case readiness across gates", stringsAsFactors = FALSE))
  rownames(gates) <- NULL
  gates
}

# Purpose: Section palette.
# Inputs: none.
# Output: Returns a deterministic scalar, vector, path, status, or empty-schema object used by downstream functions.
section_palette <- function() {
  c(Region_1 = "#3C5488", Region_2 = "#00A087", Region_3 = "#E64B35", Region_4 = "#F39B7F")
}

# Purpose: Section downstream status.
# Inputs: required: region_id.
# Output: Returns a deterministic scalar, vector, path, status, or empty-schema object used by downstream functions.
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

# Purpose: Cell style theme.
# Inputs: optional/defaulted: base_size.
# Output: Returns the derived R object described by the function name; no files are written unless an explicit output path is an input.
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

# Purpose: Plot section qc.
# Inputs: required: cell_metadata, region_id.
# Output: Returns a ggplot object or named list of plots; plotting does not mutate the input object.
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
    spatial = ggplot2::ggplot(spatial_data, ggplot2::aes(x = x_centroid, y = y_centroid, colour = QC_status)) +
      ggplot2::geom_point(size = 0.35, alpha = 0.75) +
      ggplot2::scale_colour_manual(values = c(Pass = "#BDBDBD", Review = "#D73027"), drop = FALSE) +
      ggplot2::coord_fixed() + ggplot2::labs(title = "Spatial QC review map", subtitle = region_id, x = "X centroid", y = "Y centroid", colour = "QC") + cell_style_theme()
  )
}

# Purpose: Save section plots.
# Inputs: required: plots, figure_dir, region_id, project_root.
# Output: Writes validated artifact file(s) and returns their path(s) invisibly or as a named path list.
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

# Purpose: Write gz tsv.
# Inputs: required: x, path, project_root.
# Output: Writes validated artifact file(s) and returns their path(s) invisibly or as a named path list.
write_gz_tsv <- function(x, path, project_root) {
  assert_path_within(project_root, path)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  con <- gzfile(path, "wt"); on.exit(close(con), add = TRUE)
  utils::write.table(x, con, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  invisible(path)
}

# Purpose: Plot extended spatial qc.
# Inputs: required: spatial_cells, spatial_edge_density, spatial_hotspots, region_id.
# Output: Returns a ggplot object or named list of plots; plotting does not mutate the input object.
plot_extended_spatial_qc <- function(spatial_cells, spatial_edge_density, spatial_hotspots, region_id) {
  require_package("ggplot2")
  validate_spatial_cells(spatial_cells)
  required <- c("edge_proxy", "qc_review_flag")
  missing <- setdiff(required, names(spatial_cells))
  if (length(missing)) stop(sprintf("Extended spatial plot data missing columns: %s", paste(missing, collapse = ", ")), call. = FALSE)
  cells <- spatial_cells
  cells$review_status <- factor(ifelse(cells$qc_review_flag, "Review", "Pass"), levels = c("Pass", "Review"))
  cells$location_class <- ifelse(cells$edge_proxy, "Edge proxy", "Interior")
  if ("dense_aggregate" %in% names(cells)) {
    cells$density_class <- ifelse(cells$dense_aggregate, "Dense aggregate proxy", "Other")
  } else {
    cells$density_class <- "Density not calculated"
  }
  base_map <- ggplot2::ggplot(cells, ggplot2::aes(x = x_centroid, y = y_centroid)) +
    ggplot2::coord_fixed() + cell_style_theme() +
    ggplot2::labs(subtitle = region_id, x = "X centroid (microns)", y = "Y centroid (microns)")
  hotspot_plot <- base_map +
    ggplot2::geom_point(ggplot2::aes(colour = review_status), size = 0.45, alpha = 0.75) +
    ggplot2::scale_colour_manual(values = c(Pass = "#BDBDBD", Review = "#D73027"), drop = FALSE) +
    ggplot2::labs(title = "QC review flags and candidate spatial hotspots", colour = "QC")
  review_hotspots <- spatial_hotspots[spatial_hotspots$hotspot_status == "MORPHOLOGY_REVIEW_REQUIRED", , drop = FALSE]
  if (nrow(review_hotspots)) {
    hotspot_plot <- hotspot_plot + ggplot2::geom_rect(
      data = review_hotspots,
      ggplot2::aes(xmin = x_min, xmax = x_max, ymin = y_min, ymax = y_max),
      inherit.aes = FALSE, fill = NA, colour = "#6A3D9A", linewidth = 0.7
    )
  }
  list(
    review_map = base_map +
      ggplot2::geom_point(ggplot2::aes(colour = review_status, shape = location_class), size = 0.5, alpha = 0.75) +
      ggplot2::scale_colour_manual(values = c(Pass = "#BDBDBD", Review = "#D73027"), drop = FALSE) +
      ggplot2::labs(title = "Spatial QC review map with tissue-edge proxy", colour = "QC", shape = "Location"),
    edge_density = ggplot2::ggplot(spatial_edge_density[is.finite(spatial_edge_density$review_rate), , drop = FALSE], ggplot2::aes(x = class, y = review_rate, fill = class_type)) +
      ggplot2::geom_col(width = 0.7, alpha = 0.9) +
      ggplot2::facet_wrap(~class_type, scales = "free_x") +
      ggplot2::scale_y_continuous(labels = function(x) paste0(round(100 * x, 1), "%")) +
      ggplot2::labs(title = "QC review enrichment at edge and density proxies", subtitle = region_id,
                    x = NULL, y = "Review-flag rate", fill = "Proxy") + cell_style_theme(),
    hotspots = hotspot_plot
  )
}

# Purpose: Extended section required artifacts.
# Inputs: required: region_id; optional/defaulted: mode.
# Output: Returns the derived R object described by the function name; no files are written unless an explicit output path is an input.
extended_section_required_artifacts <- function(region_id, mode = "LOCAL_SUBSET") {
  mode <- toupper(as.character(mode))
  if (length(mode) != 1L || !mode %in% c("LOCAL_SUBSET", "FULL_HPC")) stop("Extended section artifact mode must be LOCAL_SUBSET or FULL_HPC.", call. = FALSE)
  c(
    "extended_qc_preflight.tsv", "cycle_alarm_evidence.tsv", "gene_transcript_quality.tsv",
    "spatial_qc_global.tsv", "spatial_qc_edge_density.tsv", "spatial_qc_hotspots.tsv",
    "spatial_qc_cell_annotations.tsv.gz", "spatial_manual_review_manifest.tsv",
    "extended_qc_status.tsv", file.path("figures", paste0(region_id, "_extended_spatial_qc.pdf"))
  )
}

# Purpose: Save extended spatial plots.
# Inputs: required: plots, figure_dir, region_id, project_root.
# Output: Writes validated artifact file(s) and returns their path(s) invisibly or as a named path list.
save_extended_spatial_plots <- function(plots, figure_dir, region_id, project_root) {
  require_package("ggplot2")
  if (!length(plots) || any(!vapply(plots, inherits, logical(1), what = "ggplot"))) stop("Extended spatial plots must be a non-empty named list of ggplot objects.", call. = FALSE)
  assert_path_within(project_root, figure_dir)
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  pdf_path <- file.path(figure_dir, paste0(region_id, "_extended_spatial_qc.pdf"))
  grDevices::pdf(pdf_path, width = 8, height = 5.5, onefile = TRUE)
  on.exit(if (grDevices::dev.cur() > 1L) grDevices::dev.off(), add = TRUE)
  for (plot in plots) print(plot)
  grDevices::dev.off()
  pdf_path
}

# Purpose: Write extended section artifacts.
# Inputs: required: project_root, output_dir, region_id, mode, preflight, cycle_alarm_evidence, gene_quality, spatial_global, spatial_edge_density, spatial_hotspots, spatial_cells, manual_review_manifest, plots.
# Output: Writes validated artifact file(s) and returns their path(s) invisibly or as a named path list.
write_extended_section_artifacts <- function(project_root, output_dir, region_id, mode, preflight,
                                             cycle_alarm_evidence, gene_quality, spatial_global,
                                             spatial_edge_density, spatial_hotspots, spatial_cells,
                                             manual_review_manifest, plots) {
  mode <- toupper(as.character(mode))
  extended_section_required_artifacts(region_id, mode)
  assert_path_within(project_root, output_dir)
  if (!"transcript_status" %in% names(gene_quality)) stop("Gene-quality table requires transcript_status.", call. = FALSE)
  if (mode == "FULL_HPC" && any(gene_quality$transcript_status == "NOT_RUN_LOCAL_SUBSET", na.rm = TRUE)) {
    stop("FULL_HPC artifacts cannot contain NOT_RUN_LOCAL_SUBSET transcript status.", call. = FALSE)
  }
  if (mode == "FULL_HPC" && any(preflight$status == "FAIL")) stop("FULL_HPC preflight contains FAIL checks.", call. = FALSE)
  tables <- list(
    extended_qc_preflight.tsv = preflight,
    cycle_alarm_evidence.tsv = cycle_alarm_evidence,
    gene_transcript_quality.tsv = gene_quality,
    spatial_qc_global.tsv = spatial_global,
    spatial_qc_edge_density.tsv = spatial_edge_density,
    spatial_qc_hotspots.tsv = spatial_hotspots,
    spatial_manual_review_manifest.tsv = manual_review_manifest
  )
  table_paths <- vapply(names(tables), function(name) write_tsv(tables[[name]], file.path(output_dir, name), project_root), character(1))
  cell_path <- write_gz_tsv(spatial_cells, file.path(output_dir, "spatial_qc_cell_annotations.tsv.gz"), project_root)
  status <- data.frame(
    region_id = region_id, mode = mode,
    transcript_status = if (all(gene_quality$transcript_status == "NOT_RUN_LOCAL_SUBSET")) "NOT_RUN_LOCAL_SUBSET" else "COMPUTED",
    diagnostic_status = if (mode == "LOCAL_SUBSET") "LOCAL_SUBSET_COMPLETE_WITH_EXPECTED_SKIPS" else "FULL_HPC_COMPLETE",
    cells = nrow(spatial_cells), cells_deleted = 0L,
    cycle_identity_status = "CYCLE_IDENTITY_UNRESOLVED_REQUIRES_10X",
    morphology_status = "COORDINATE_DIAGNOSTICS_REQUIRE_IMAGE_REVIEW",
    stringsAsFactors = FALSE
  )
  status_path <- write_tsv(status, file.path(output_dir, "extended_qc_status.tsv"), project_root)
  plot_path <- save_extended_spatial_plots(plots, file.path(output_dir, "figures"), region_id, project_root)
  paths <- c(unname(table_paths), cell_path, status_path, plot_path)
  if (!validate_extended_section_artifacts(output_dir, region_id, mode)) stop("Extended section artifact reload validation failed.", call. = FALSE)
  paths
}

# Purpose: Read extended section artifacts.
# Inputs: required: output_dir, region_id; optional/defaulted: mode.
# Output: Returns parsed, validated R data (vector, data frame, sparse matrix bundle, or named list according to the input format).
read_extended_section_artifacts <- function(output_dir, region_id, mode = "LOCAL_SUBSET") {
  if (!validate_extended_section_artifacts(output_dir, region_id, mode, stop_on_error = TRUE)) stop("Extended section artifact validation failed.", call. = FALSE)
  read_table <- function(name) utils::read.delim(file.path(output_dir, name), check.names = FALSE, stringsAsFactors = FALSE)
  list(
    preflight = read_table("extended_qc_preflight.tsv"),
    cycle_alarm_evidence = read_table("cycle_alarm_evidence.tsv"),
    gene_quality = read_table("gene_transcript_quality.tsv"),
    spatial_global = read_table("spatial_qc_global.tsv"),
    spatial_edge_density = read_table("spatial_qc_edge_density.tsv"),
    spatial_hotspots = read_table("spatial_qc_hotspots.tsv"),
    spatial_cells = utils::read.delim(gzfile(file.path(output_dir, "spatial_qc_cell_annotations.tsv.gz")), check.names = FALSE, stringsAsFactors = FALSE),
    manual_review_manifest = read_table("spatial_manual_review_manifest.tsv"),
    status = read_table("extended_qc_status.tsv")
  )
}

# Purpose: Validate extended section artifacts.
# Inputs: required: output_dir, region_id; optional/defaulted: mode, stop_on_error.
# Output: Returns validation evidence/TRUE (or the validated value) and stops with an informative error when the contract fails.
validate_extended_section_artifacts <- function(output_dir, region_id, mode = "LOCAL_SUBSET", stop_on_error = FALSE) {
  fail <- function(message) {
    if (isTRUE(stop_on_error)) stop(message, call. = FALSE)
    FALSE
  }
  mode <- toupper(as.character(mode))
  required <- tryCatch(extended_section_required_artifacts(region_id, mode), error = function(error) return(NULL))
  if (is.null(required)) return(fail("Invalid extended section artifact mode."))
  absent <- required[!file.exists(file.path(output_dir, required))]
  if (length(absent)) return(fail(sprintf("Missing extended section artifacts: %s", paste(absent, collapse = ", "))))
  gene_quality <- tryCatch(utils::read.delim(file.path(output_dir, "gene_transcript_quality.tsv"), check.names = FALSE, stringsAsFactors = FALSE), error = identity)
  if (inherits(gene_quality, "error") || !"transcript_status" %in% names(gene_quality)) return(fail("Invalid gene_transcript_quality.tsv."))
  if (mode == "FULL_HPC" && any(gene_quality$transcript_status == "NOT_RUN_LOCAL_SUBSET", na.rm = TRUE)) {
    return(fail("FULL_HPC validation rejects NOT_RUN_LOCAL_SUBSET transcript status."))
  }
  status <- tryCatch(utils::read.delim(file.path(output_dir, "extended_qc_status.tsv"), check.names = FALSE, stringsAsFactors = FALSE), error = identity)
  if (inherits(status, "error") || nrow(status) != 1L || !identical(as.character(status$region_id), region_id) || !identical(as.character(status$mode), mode)) {
    return(fail("Extended QC status does not match section or mode."))
  }
  TRUE
}

# Purpose: Evidence only required artifacts.
# Inputs: none.
# Output: Returns the derived R object described by the function name; no files are written unless an explicit output path is an input.
evidence_only_required_artifacts <- function() {
  c(
    "cell_downstream_masks.tsv.gz",
    "section_downstream_decision.tsv",
    "gene_downstream_decision.tsv",
    "eos_gene_decision_summary.tsv",
    "hotspot_sensitivity_decision.tsv",
    "evidence_only_qc_release.tsv"
  )
}

# Purpose: Add evidence provenance.
# Inputs: required: table, run_label, execution_mode, provenance, source_artifact, generated_utc.
# Output: Returns a newly constructed or annotated R object while preserving the supplied raw object/data rows.
add_evidence_provenance <- function(table, run_label, execution_mode, provenance,
                                    source_artifact, generated_utc) {
  table$run_label <- run_label
  table$execution_mode <- toupper(execution_mode)
  table$generated_utc <- generated_utc
  table$source_artifact <- source_artifact
  table$provenance <- provenance
  table
}

# Purpose: Build one section downstream decision.
# Inputs: required: masks, run_label, execution_mode, provenance, generated_utc.
# Output: Returns a newly constructed or annotated R object while preserving the supplied raw object/data rows.
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

# Purpose: Build section downstream decision.
# Inputs: required: masks, run_label, execution_mode, provenance, generated_utc.
# Output: Returns a newly constructed or annotated R object while preserving the supplied raw object/data rows.
build_section_downstream_decision <- function(masks, run_label, execution_mode,
                                              provenance, generated_utc) {
  regions <- paste0("Region_", 1:4)
  if (!identical(sort(unique(as.character(masks$region_id))), regions)) {
    stop("Slide section decisions require Region_1 through Region_4.", call. = FALSE)
  }
  rows <- lapply(regions, function(region) {
    build_one_section_downstream_decision(
      masks[masks$region_id == region, , drop = FALSE], run_label,
      execution_mode, provenance, generated_utc
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# Purpose: Build hotspot sensitivity decision.
# Inputs: required: spatial_hotspots, masks, run_label, execution_mode, provenance, generated_utc.
# Output: Returns a newly constructed or annotated R object while preserving the supplied raw object/data rows.
build_hotspot_sensitivity_decision <- function(spatial_hotspots, masks, run_label,
                                               execution_mode, provenance, generated_utc) {
  positive <- spatial_hotspots[
    spatial_hotspots$region_id == "Region_3" &
      spatial_hotspots$hotspot_status == "MORPHOLOGY_REVIEW_REQUIRED", , drop = FALSE
  ]
  if (nrow(positive)) {
    positive$hotspot_cells_from_mask <- vapply(as.character(positive$grid_id), function(id) {
      sum(masks$region_id == "Region_3" & masks$grid_id == id, na.rm = TRUE)
    }, integer(1))
    positive$primary_handling <- "KEEP_IN_PRIMARY"
    positive$sensitivity_handling <- "EXCLUDE_IN_HOTSPOT_SENSITIVITY"
    positive$morphology_status <- "NOT_CONFIRMED_AS_ARTIFACT"
    positive$decision_rule <- "FDR-positive Region_3 bins are retained in primary and excluded only in hotspot sensitivity"
    out <- positive
  } else {
    out <- data.frame(
      region_id = "Region_3", grid_id = NA_character_, cells = 0L, flagged = 0L,
      review_rate = NA_real_, global_rate = NA_real_, rate_difference = NA_real_,
      empirical_p = NA_real_, adjusted_p = NA_real_, hotspot_status = "NO_FDR_POSITIVE_BIN_IN_THIS_RUN",
      x_min = NA_real_, x_max = NA_real_, y_min = NA_real_, y_max = NA_real_,
      fdr_threshold = NA_real_, min_bin_cells_threshold = NA_integer_,
      permutations = NA_integer_, seed = NA_integer_,
      hotspot_cells_from_mask = 0L, primary_handling = "KEEP_IN_PRIMARY",
      sensitivity_handling = "NO_CELLS_TO_EXCLUDE",
      morphology_status = "NOT_ESTIMABLE_LOCAL_SUBSET",
      decision_rule = "Full-data hotspot reconciliation required; no primary exclusion is made",
      stringsAsFactors = FALSE
    )
  }
  add_evidence_provenance(
    out, run_label, execution_mode, provenance,
    "slide_summary/combined_spatial_hotspots.tsv", generated_utc
  )
}

# Purpose: Build evidence only release.
# Inputs: required: section_decision, masks, gene_decision, eos_decision, hotspot_decision, run_label, execution_mode, provenance, generated_utc.
# Output: Returns a newly constructed or annotated R object while preserving the supplied raw object/data rows.
build_evidence_only_release <- function(section_decision, masks, gene_decision,
                                        eos_decision, hotspot_decision,
                                        run_label, execution_mode, provenance,
                                        generated_utc) {
  mode <- toupper(execution_mode)
  counts <- c(
    raw = nrow(gene_decision),
    conservative = sum(gene_decision$conservative_evidence_status == "CONSERVATIVE_NO_SIGNAL_DETECTED"),
    provisional = sum(gene_decision$primary_feature_status == "PROVISIONAL_PRIMARY_FEATURES"),
    risk = sum(gene_decision$technical_risk_status == "TECHNICAL_RISK_SENSITIVITY_ONLY"),
    eos = sum(eos_decision$retained_provisional)
  )
  full_expected <- identical(mode, "FULL_HPC")
  gene_reconciled <- identical(
    as.integer(unname(counts[c("raw", "conservative", "provisional", "risk")])),
    c(479L, 67L, 245L, 234L)
  )
  eos_partition <- table(factor(eos_decision$gene_set[eos_decision$retained_provisional], levels = c("common", "short_lived", "long_lived")))
  eos_reconciled <- sum(eos_partition) == 53L && identical(as.integer(eos_partition), c(4L, 27L, 22L))
  qc_rows <- data.frame(
    gate_id = c("fixed_section_decisions", "cell_mask_reconciliation", "gene_tier_reconciliation",
                "eos_gene_reconciliation", "region4_excluded_from_reference_definition"),
    scope = c("sections", "cells", "genes", "Eos genes", "Region_4"),
    gate_status = c(
      if (identical(section_decision$section_status, c("PRIMARY_CONDITIONAL", "PRIMARY_CONDITIONAL", "PRIMARY", "SENSITIVITY_ONLY"))) "PASS" else "STOP",
      if (!anyNA(masks[, c("primary_include", "strict_include", "hotspot_sensitivity_include")]) && all(!masks$strict_include | masks$primary_include) && all(!masks$hotspot_sensitivity_include | masks$primary_include)) "PASS" else "STOP",
      if (full_expected) if (gene_reconciled) "PASS" else "STOP" else "PENDING_DOWNSTREAM_ANALYSIS",
      if (full_expected) if (eos_reconciled) "PASS" else "STOP" else "PENDING_DOWNSTREAM_ANALYSIS",
      if (!section_decision$cluster_discovery_eligible[section_decision$region_id == "Region_4"]) "PASS" else "STOP"
    ),
    measured_value = c(
      paste(section_decision$region_id, section_decision$section_status, sep = "=", collapse = ";"),
      sprintf("primary=%d;strict=%d;hotspot_sensitivity=%d;total=%d", sum(masks$primary_include), sum(masks$strict_include), sum(masks$hotspot_sensitivity_include), nrow(masks)),
      sprintf("raw=%d;conservative=%d;provisional=%d;risk=%d", counts[["raw"]], counts[["conservative"]], counts[["provisional"]], counts[["risk"]]),
      sprintf("retained=%d;common=%d;short_lived=%d;long_lived=%d", counts[["eos"]], eos_partition[[1]], eos_partition[[2]], eos_partition[[3]]),
      "cluster_discovery_eligible=FALSE"
    ),
    threshold = c("R1/R2 conditional; R3 primary; R4 sensitivity-only", "strict and hotspot masks must be subsets of primary", "FULL_HPC: 479/67/245/234", "FULL_HPC: 53=4/27/22", "FALSE"),
    evidence_source = c("section_downstream_decision.tsv", "cell_downstream_masks.tsv.gz", "gene_downstream_decision.tsv", "eos_gene_decision_summary.tsv", "section_downstream_decision.tsv"),
    interpretation = c("Fixed approved decisions", "Non-destructive masks reconcile", "Local subset cannot freeze full-data counts", "Local subset cannot freeze full-data Eos retention", "Region_4 is mapping-only"),
    next_required_artifact = c("none", "none", if (full_expected) "none" else "FULL_HPC summary", if (full_expected) "none" else "FULL_HPC summary", "Region_1-3 reference"),
    stringsAsFactors = FALSE
  )
  pending <- data.frame(
    gate_id = c("section_identity_dominance", "primary_strict_celltype_stability", "feature245_vs_gene67_celltype_agreement",
                "eos_state_stability", "region3_hotspot_conclusion_stability", "region4_mapping_quality", "technical_risk_gene_dependence"),
    scope = c("PCA/clustering", "major cell types", "major cell types", "Eos state", "Region_3", "Region_4", "primary conclusions"),
    gate_status = "PENDING_DOWNSTREAM_ANALYSIS", measured_value = NA_character_,
    threshold = c("section identity must not dominate", "stable assignments", "no contradictory major cell types", "no reversal or unstable assignments", "main conclusion unchanged", "confidence acceptable and not strongly cell-type dependent", "conclusion not primarily dependent on 234 risk genes"),
    evidence_source = c("downstream PCA and cluster diagnostics", "primary/strict label comparison", "245/67 sensitivity comparison", "Region_3 Eos robustness grid", "Region_3 hotspot sensitivity", "held-out calibrated label transfer", "risk-gene dependence sensitivity"),
    interpretation = "Primary release cannot pass until downstream evidence is supplied",
    next_required_artifact = c("PCA section-mixing diagnostics", "cell-type stability table", "245-vs-67 comparison", "Eos stability table", "hotspot sensitivity comparison", "Region_4 mapping metrics", "technical-risk dependence analysis"),
    stringsAsFactors = FALSE
  )
  out <- rbind(qc_rows, pending)
  overall <- if (any(out$gate_status == "STOP")) "STOP" else if (any(out$gate_status == "PENDING_DOWNSTREAM_ANALYSIS")) "PENDING_DOWNSTREAM_ANALYSIS" else "PASS"
  out <- rbind(out, data.frame(
    gate_id = "overall_primary_release", scope = "scWAT primary release", gate_status = overall,
    measured_value = overall, threshold = "PASS only after every required gate passes",
    evidence_source = "all evidence-only and downstream gates",
    interpretation = "Automatic worst-case release status", next_required_artifact = "all pending downstream artifacts",
    stringsAsFactors = FALSE
  ))
  add_evidence_provenance(out, run_label, execution_mode, provenance, "evidence-only QC decision engine", generated_utc)
}

# Purpose: Summarise evidence only qc.
# Inputs: required: extended_slide_data, candidates, eos_gene_sets, run_label, execution_mode, provenance.
# Output: Returns computed QC evidence as a vector, data frame, or named summary list; it does not modify raw input files.
summarise_evidence_only_qc <- function(extended_slide_data, candidates, eos_gene_sets,
                                       run_label, execution_mode, provenance) {
  generated_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  masks <- build_cell_downstream_masks(
    extended_slide_data$spatial_cells, extended_slide_data$spatial_hotspots, provenance
  )
  masks <- add_evidence_provenance(
    masks, run_label, execution_mode, provenance,
    "sections/*/spatial_qc_cell_annotations.tsv.gz", generated_utc
  )
  masks$section_input_cells <- ave(rep(1L, nrow(masks)), masks$region_id, FUN = sum)
  masks$section_primary_include_cells <- ave(as.integer(masks$primary_include), masks$region_id, FUN = sum)
  masks$section_strict_include_cells <- ave(as.integer(masks$strict_include), masks$region_id, FUN = sum)
  masks$section_hotspot_sensitivity_include_cells <- ave(as.integer(masks$hotspot_sensitivity_include), masks$region_id, FUN = sum)
  masks$qc_threshold_source <- file.path("sections", masks$region_id, "qc_thresholds.tsv")
  sections <- build_section_downstream_decision(
    masks, run_label, execution_mode, provenance, generated_utc
  )
  genes <- build_gene_downstream_decision(
    candidates, unique(as.character(extended_slide_data$gene_quality$gene)),
    run_label, execution_mode, provenance
  )
  genes$generated_utc <- generated_utc
  genes$source_artifact <- "slide_summary/candidate_cycle_affected_genes.tsv"
  eos <- build_eos_gene_decision(genes, eos_gene_sets, run_label, execution_mode, provenance)
  eos$generated_utc <- generated_utc
  eos$source_artifact <- "config/eos_gene_sets.tsv + gene_downstream_decision.tsv"
  eos$retained_total <- sum(eos$retained_provisional)
  eos$retained_gene_set_count <- ave(as.integer(eos$retained_provisional), eos$gene_set, FUN = sum)
  hotspots <- build_hotspot_sensitivity_decision(
    extended_slide_data$spatial_hotspots, masks, run_label, execution_mode,
    provenance, generated_utc
  )
  release <- build_evidence_only_release(
    sections, masks, genes, eos, hotspots, run_label, execution_mode,
    provenance, generated_utc
  )
  list(cell_masks = masks, sections = sections, genes = genes, eos = eos,
       hotspots = hotspots, release = release)
}

# Purpose: Write evidence only qc artifacts.
# Inputs: required: project_root, run_root, evidence_summary.
# Output: Writes validated artifact file(s) and returns their path(s) invisibly or as a named path list.
write_evidence_only_qc_artifacts <- function(project_root, run_root, evidence_summary) {
  require_package("Matrix")
  assert_path_within(project_root, run_root)
  slide_dir <- file.path(run_root, "slide_summary")
  downstream_dir <- file.path(run_root, "downstream_inputs")
  dir.create(slide_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(downstream_dir, recursive = TRUE, showWarnings = FALSE)
  tables <- list(
    section_downstream_decision.tsv = evidence_summary$sections,
    gene_downstream_decision.tsv = evidence_summary$genes,
    eos_gene_decision_summary.tsv = evidence_summary$eos,
    hotspot_sensitivity_decision.tsv = evidence_summary$hotspots,
    evidence_only_qc_release.tsv = evidence_summary$release
  )
  paths <- vapply(names(tables), function(name) {
    write_tsv(tables[[name]], file.path(slide_dir, name), project_root)
  }, character(1))
  mask_path <- write_gz_tsv(
    evidence_summary$cell_masks, file.path(slide_dir, "cell_downstream_masks.tsv.gz"), project_root
  )
  region_paths <- character(4)
  manifest_rows <- vector("list", 4)
  for (index in seq_len(4L)) {
    region <- paste0("Region_", index)
    section_dir <- file.path(run_root, "sections", region)
    section_rds <- file.path(section_dir, paste0(region, ".phase0_2_qc.rds"))
    if (!file.exists(section_rds)) stop(sprintf("Missing raw section object: %s", section_rds), call. = FALSE)
    section_object <- readRDS(section_rds)
    if (!isTRUE(section_object$raw_counts_preserved) || !inherits(section_object$counts, "sparseMatrix")) {
      stop(sprintf("Section object does not preserve sparse raw counts for %s.", region), call. = FALSE)
    }
    region_masks <- evidence_summary$cell_masks[evidence_summary$cell_masks$region_id == region, , drop = FALSE]
    order_index <- match(colnames(section_object$counts), region_masks$cell_id)
    if (anyNA(order_index)) stop(sprintf("Mask alignment failed for %s.", region), call. = FALSE)
    region_masks <- region_masks[order_index, , drop = FALSE]
    gene_sets <- list(
      provisional_primary_features = evidence_summary$genes$gene[evidence_summary$genes$primary_feature_status == "PROVISIONAL_PRIMARY_FEATURES"],
      conservative_no_signal_detected = evidence_summary$genes$gene[evidence_summary$genes$conservative_evidence_status == "CONSERVATIVE_NO_SIGNAL_DETECTED"],
      technical_risk_sensitivity_only = evidence_summary$genes$gene[evidence_summary$genes$technical_risk_status == "TECHNICAL_RISK_SENSITIVITY_ONLY"],
      raw_complete_panel = evidence_summary$genes$gene,
      eos_provisional_primary = evidence_summary$eos$gene[evidence_summary$eos$retained_provisional]
    )
    bundle <- list(
      schema_version = "evidence_only_qc_v1", region_id = region,
      section_status = section_downstream_status(region),
      counts = section_object$counts, features = section_object$features,
      cell_metadata = region_masks, gene_sets = gene_sets,
      raw_counts_preserved = TRUE,
      downstream_contract = if (region == "Region_4") "MAP_TO_REGION_1_3_REFERENCE_WITH_UNCERTAIN" else "REFERENCE_ELIGIBILITY_FROM_CELL_MASKS",
      mapping_requirements = if (region == "Region_4") c("mapping_confidence", "reference_distance", "second_best_label", "confidence_margin", "Uncertain") else character(),
      provenance = unique(region_masks$provenance)
    )
    region_path <- file.path(downstream_dir, paste0(region, ".downstream_input.rds"))
    saveRDS(bundle, region_path, compress = FALSE)
    region_paths[[index]] <- region_path
    write_gz_tsv(region_masks, file.path(section_dir, "cell_downstream_masks.tsv.gz"), project_root)
    write_tsv(evidence_summary$sections[evidence_summary$sections$region_id == region, , drop = FALSE],
              file.path(section_dir, "section_downstream_decision.tsv"), project_root)
    manifest_rows[[index]] <- data.frame(
      region_id = region, section_status = bundle$section_status,
      cells = ncol(bundle$counts), genes = nrow(bundle$counts),
      primary_include_cells = sum(bundle$cell_metadata$primary_include),
      strict_include_cells = sum(bundle$cell_metadata$strict_include),
      hotspot_sensitivity_include_cells = sum(bundle$cell_metadata$hotspot_sensitivity_include),
      path = normalizePath(region_path, winslash = "/", mustWork = TRUE),
      schema_version = bundle$schema_version, raw_counts_preserved = TRUE,
      downstream_contract = bundle$downstream_contract, stringsAsFactors = FALSE
    )
  }
  manifest <- do.call(rbind, manifest_rows)
  manifest_path <- write_tsv(manifest, file.path(downstream_dir, "downstream_input_manifest.tsv"), project_root)
  all_paths <- c(mask_path, unname(paths), region_paths, manifest_path)
  if (!validate_evidence_only_qc_artifacts(run_root, stop_on_error = TRUE)) stop("Evidence-only QC reload validation failed.", call. = FALSE)
  all_paths
}

# Purpose: Validate evidence only qc artifacts.
# Inputs: required: run_root; optional/defaulted: stop_on_error.
# Output: Returns validation evidence/TRUE (or the validated value) and stops with an informative error when the contract fails.
validate_evidence_only_qc_artifacts <- function(run_root, stop_on_error = FALSE) {
  require_package("Matrix")
  fail <- function(message) {
    if (isTRUE(stop_on_error)) stop(message, call. = FALSE)
    FALSE
  }
  slide_dir <- file.path(run_root, "slide_summary")
  required <- file.path(slide_dir, evidence_only_required_artifacts())
  absent <- required[!file.exists(required)]
  if (length(absent)) return(fail(sprintf("Missing evidence-only artifacts: %s", paste(basename(absent), collapse = ", "))))
  sections <- tryCatch(utils::read.delim(file.path(slide_dir, "section_downstream_decision.tsv"), check.names = FALSE), error = identity)
  genes <- tryCatch(utils::read.delim(file.path(slide_dir, "gene_downstream_decision.tsv"), check.names = FALSE), error = identity)
  release <- tryCatch(utils::read.delim(file.path(slide_dir, "evidence_only_qc_release.tsv"), check.names = FALSE), error = identity)
  masks <- tryCatch(utils::read.delim(gzfile(file.path(slide_dir, "cell_downstream_masks.tsv.gz")), check.names = FALSE), error = identity)
  if (any(vapply(list(sections, genes, release, masks), inherits, logical(1), what = "error"))) return(fail("An evidence-only table cannot be reloaded."))
  if (!identical(as.character(sections$region_id), paste0("Region_", 1:4))) return(fail("Section decisions must contain Region_1 through Region_4 in order."))
  if (any(grepl("CONFIRMED_(AFFECTED|UNAFFECTED)", unlist(genes)))) return(fail("Prohibited confirmed gene terminology detected."))
  if (!all(release$gate_status %in% c("PASS", "STOP", "PENDING_DOWNSTREAM_ANALYSIS"))) return(fail("Invalid release-gate status."))
  for (region in paste0("Region_", 1:4)) {
    path <- file.path(run_root, "downstream_inputs", paste0(region, ".downstream_input.rds"))
    if (!file.exists(path)) return(fail(sprintf("Missing downstream bundle for %s.", region)))
    object <- tryCatch(readRDS(path), error = identity)
    if (inherits(object, "error") || !isTRUE(object$raw_counts_preserved) ||
        ncol(object$counts) != nrow(object$cell_metadata) ||
        !identical(colnames(object$counts), as.character(object$cell_metadata$cell_id))) {
      return(fail(sprintf("Invalid downstream bundle for %s.", region)))
    }
  }
  TRUE
}

# Purpose: Extended slide required artifacts.
# Inputs: none.
# Output: Returns the derived R object described by the function name; no files are written unless an explicit output path is an input.
extended_slide_required_artifacts <- function() {
  c(
    "combined_cycle_alarm_evidence.tsv", "combined_gene_transcript_quality.tsv",
    "candidate_cycle_affected_genes.tsv", "candidate_cycle_affected_genes_affected_only.tsv",
    "subset_full_qc_ranking.tsv",
    "subset_full_qc_rank_agreement.tsv", "combined_spatial_qc.tsv",
    "combined_spatial_hotspots.tsv", "within_mouse_section_concordance.tsv",
    "within_mouse_gene_concordance.tsv", "extended_slide_qc_status.tsv",
    file.path("figures", "scwat_extended_qc_diagnostics.pdf")
  )
}

# Purpose: Rbind fill.
# Inputs: required: tables.
# Output: Returns the derived R object described by the function name; no files are written unless an explicit output path is an input.
rbind_fill <- function(tables) {
  tables <- tables[vapply(tables, is.data.frame, logical(1))]
  if (!length(tables)) return(data.frame())
  columns <- unique(unlist(lapply(tables, names), use.names = FALSE))
  normalized <- lapply(tables, function(table) {
    for (name in setdiff(columns, names(table))) table[[name]] <- NA
    table[, columns, drop = FALSE]
  })
  out <- do.call(rbind, normalized)
  rownames(out) <- NULL
  out
}

# Purpose: Validate four extended section outputs.
# Inputs: required: run_root; optional/defaulted: expected_regions.
# Output: Returns validation evidence/TRUE (or the validated value) and stops with an informative error when the contract fails.
validate_four_extended_section_outputs <- function(run_root, expected_regions = paste0("Region_", 1:4)) {
  sections_root <- file.path(run_root, "sections")
  dirs <- list.dirs(sections_root, recursive = FALSE, full.names = TRUE)
  region_ids <- basename(dirs)
  if (length(dirs) != 4L || !setequal(region_ids, expected_regions) || anyDuplicated(region_ids)) {
    stop(sprintf("Expected exactly four unique extended section outputs (%s).", paste(expected_regions, collapse = ", ")), call. = FALSE)
  }
  dirs <- dirs[match(expected_regions, region_ids)]
  statuses <- lapply(seq_along(dirs), function(index) {
    path <- file.path(dirs[[index]], "extended_qc_status.tsv")
    if (!file.exists(path)) stop(sprintf("Missing extended section status: %s", path), call. = FALSE)
    value <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
    if (nrow(value) != 1L || !all(c("region_id", "mode") %in% names(value))) stop(sprintf("Invalid extended section status: %s", path), call. = FALSE)
    if (!identical(as.character(value$region_id), expected_regions[[index]])) stop(sprintf("Extended section status region mismatch for %s.", expected_regions[[index]]), call. = FALSE)
    value
  })
  modes <- vapply(statuses, function(value) as.character(value$mode[[1]]), character(1))
  if (length(unique(modes)) != 1L) stop("Extended slide aggregation rejects mixed execution modes.", call. = FALSE)
  for (index in seq_along(dirs)) validate_extended_section_artifacts(dirs[[index]], expected_regions[[index]], modes[[index]], stop_on_error = TRUE)
  data.frame(
    region_id = expected_regions, section_output_dir = normalizePath(dirs, winslash = "/", mustWork = TRUE),
    mode = modes, stringsAsFactors = FALSE
  )
}

# Purpose: Read extended slide qc outputs.
# Inputs: required: run_root; optional/defaulted: expected_regions.
# Output: Returns parsed, validated R data (vector, data frame, sparse matrix bundle, or named list according to the input format).
read_extended_slide_qc_outputs <- function(run_root, expected_regions = paste0("Region_", 1:4)) {
  coverage <- validate_four_extended_section_outputs(run_root, expected_regions)
  bundles <- lapply(seq_len(nrow(coverage)), function(index) {
    read_extended_section_artifacts(
      coverage$section_output_dir[[index]], coverage$region_id[[index]], coverage$mode[[index]]
    )
  })
  combine <- function(name) rbind_fill(lapply(bundles, `[[`, name))
  spatial_global <- combine("spatial_global"); spatial_global$diagnostic_type <- "global_clustering"
  spatial_enrichment <- combine("spatial_edge_density"); spatial_enrichment$diagnostic_type <- "edge_density_enrichment"
  list(
    coverage = coverage, mode = unique(coverage$mode),
    cycle_alarm_evidence = combine("cycle_alarm_evidence"),
    gene_quality = combine("gene_quality"),
    spatial_global = spatial_global, spatial_edge_density = spatial_enrichment,
    spatial_qc = rbind_fill(list(spatial_global, spatial_enrichment)),
    spatial_hotspots = combine("spatial_hotspots"),
    spatial_cells = combine("spatial_cells"),
    manual_review_manifest = combine("manual_review_manifest"),
    section_status = combine("status")
  )
}

# Purpose: Summarise extended slide qc.
# Inputs: required: extended_slide_data, section_summary, manifest, config, subset_reference.
# Output: Returns computed QC evidence as a vector, data frame, or named summary list; it does not modify raw input files.
summarise_extended_slide_qc <- function(extended_slide_data, section_summary, manifest, config, subset_reference) {
  if (!identical(sort(unique(extended_slide_data$gene_quality$region_id)), paste0("Region_", 1:4))) stop("Extended slide gene-quality input must contain Region_1 through Region_4.", call. = FALSE)
  candidates <- rank_candidate_cycle_genes(extended_slide_data$gene_quality, config)
  mode <- unique(as.character(extended_slide_data$mode))
  if (length(mode) != 1L) stop("Extended slide summary requires one execution mode.", call. = FALSE)
  if (mode == "FULL_HPC") {
    ranking <- compare_subset_full_qc(section_summary, subset_reference)
    ranking$ranking$comparison_status <- "FULL_DATA_COMPARISON"
    ranking$agreement$comparison_status <- "FULL_DATA_COMPARISON"
  } else {
    observed <- section_summary[, c("region_id", "input_cells", "review_flagged"), drop = FALSE]
    observed$observed_review_fraction <- ifelse(observed$input_cells > 0, observed$review_flagged / observed$input_cells, NA_real_)
    ranking_table <- merge(subset_reference, observed, by = "region_id", all = TRUE, sort = FALSE)
    ranking_table <- ranking_table[match(paste0("Region_", 1:4), ranking_table$region_id), , drop = FALSE]
    ranking_table$full_rank <- NA_integer_
    ranking_table$comparison_status <- "NOT_RUN_LOCAL_SUBSET"
    ranking <- list(
      ranking = ranking_table,
      agreement = data.frame(spearman_rho = NA_real_, kendall_tau = NA_real_, interpretation = "FULL_DATA_REQUIRED", comparison_status = "NOT_RUN_LOCAL_SUBSET", stringsAsFactors = FALSE)
    )
  }
  cells <- extended_slide_data$spatial_cells
  if (!"control_fraction" %in% names(cells) && "control_fraction_cell" %in% names(cells)) cells$control_fraction <- cells$control_fraction_cell
  concordance <- calculate_within_mouse_concordance(manifest, section_summary, cells, extended_slide_data$gene_quality, config)
  alarm_regions <- extended_slide_data$cycle_alarm_evidence$region_id[extended_slide_data$cycle_alarm_evidence$evidence_status == "DIRECT_EVIDENCE"]
  status <- data.frame(
    scope = "scWAT_extended_slide", mode = mode, sections = 4L,
    alarm_positive_regions = paste(alarm_regions, collapse = ","),
    cycle_identity_status = "CYCLE_IDENTITY_UNRESOLVED_REQUIRES_10X",
    candidate_gene_status = "CANDIDATE_NOT_CONFIRMED",
    subset_full_status = unique(ranking$ranking$comparison_status)[[1]],
    spatial_interpretation = "COORDINATE_EVIDENCE_REQUIRES_MORPHOLOGY_REVIEW",
    concordance_interpretation = "ADVISORY_TECHNICAL_CONCORDANCE_NOT_BIOLOGICAL_TEST",
    cells_deleted = 0L, stringsAsFactors = FALSE
  )
  list(
    candidates = candidates, ranking = ranking$ranking, rank_agreement = ranking$agreement,
    spatial_qc = extended_slide_data$spatial_qc, spatial_hotspots = extended_slide_data$spatial_hotspots,
    concordance = concordance, status = status
  )
}

# Purpose: Plot extended slide qc.
# Inputs: required: extended_slide_data, extended_slide_summary.
# Output: Returns a ggplot object or named list of plots; plotting does not mutate the input object.
plot_extended_slide_qc <- function(extended_slide_data, extended_slide_summary) {
  require_package("ggplot2")
  palette <- section_palette()
  alarm <- extended_slide_data$cycle_alarm_evidence
  alarm$alarm_state <- ifelse(alarm$evidence_status == "DIRECT_EVIDENCE", "Alarm reported", "No alarm reported")
  candidate <- extended_slide_summary$candidates
  candidate <- candidate[candidate$section_candidate_flag, , drop = FALSE]
  candidate_counts <- as.data.frame(table(candidate$region_id, candidate$section_evidence_status), stringsAsFactors = FALSE)
  names(candidate_counts) <- c("region_id", "section_evidence_status", "genes")
  ranking <- extended_slide_summary$ranking
  ranking_long <- rbind(
    data.frame(region_id = ranking$region_id, source = "Validated subset", review_fraction = ranking$review_fraction),
    data.frame(region_id = ranking$region_id, source = ifelse(ranking$comparison_status == "FULL_DATA_COMPARISON", "Full data", "Current local subset"), review_fraction = ranking$observed_review_fraction %||% ranking$full_review_fraction)
  )
  global <- extended_slide_data$spatial_global
  concordance <- extended_slide_summary$concordance$summary
  list(
    alarm_evidence = ggplot2::ggplot(alarm, ggplot2::aes(region_id, 1, fill = alarm_state)) +
      ggplot2::geom_col(width = 0.7) + ggplot2::scale_fill_manual(values = c("Alarm reported" = "#D73027", "No alarm reported" = "#4DAF4A")) +
      ggplot2::labs(title = "Direct poor-cycle alarm evidence", x = NULL, y = NULL, fill = NULL) + cell_style_theme() +
      ggplot2::theme(axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank()),
    candidate_genes = ggplot2::ggplot(candidate_counts, ggplot2::aes(region_id, genes, fill = section_evidence_status)) +
      ggplot2::geom_col(width = 0.7) + ggplot2::labs(title = "Section-level candidate affected genes",x = NULL, y = "Candidate genes", fill = "Evidence") + cell_style_theme(),
    ranking = ggplot2::ggplot(ranking_long, ggplot2::aes(region_id, review_fraction, colour = source, group = source)) +
      ggplot2::geom_line(linewidth = 0.7) + ggplot2::geom_point(size = 2) +
      ggplot2::scale_y_continuous(labels = function(x) paste0(round(100 * x, 1), "%")) +
      ggplot2::labs(title = "Subset versus full-data review burden", x = NULL, y = "Review-flag rate", colour = NULL) + cell_style_theme(),
    spatial = ggplot2::ggplot(global, ggplot2::aes(region_id, statistic, colour = region_id)) +
      ggplot2::geom_hline(yintercept = 0, colour = "#BDBDBD") + ggplot2::geom_point(size = 2) +
      ggplot2::scale_colour_manual(values = palette, drop = FALSE) +
      ggplot2::labs(title = "Global spatial clustering of QC review flags", x = NULL, y = "kNN clustering statistic") + cell_style_theme() + ggplot2::theme(legend.position = "none"),
    concordance = ggplot2::ggplot(concordance, ggplot2::aes(mouse_id, review_rate_difference, fill = concordance_status)) +
      ggplot2::geom_col(width = 0.65) + ggplot2::scale_fill_manual(values = c(CONCORDANT = "#4DAF4A", REVIEW = "#D73027", NOT_ESTIMABLE = "#BDBDBD"), drop = FALSE) +
      ggplot2::labs(title = "Within-mouse technical concordance", x = NULL, y = "Absolute review-rate difference", fill = NULL) + cell_style_theme()
  )
}

# Purpose: Write extended slide qc artifacts.
# Inputs: required: project_root, run_root, extended_slide_data, extended_slide_summary, plots.
# Output: Writes validated artifact file(s) and returns their path(s) invisibly or as a named path list.
write_extended_slide_qc_artifacts <- function(project_root, run_root, extended_slide_data, extended_slide_summary, plots) {
  assert_path_within(project_root, run_root)
  output_dir <- file.path(run_root, "slide_summary")
  figure_dir <- file.path(output_dir, "figures")
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  tables <- list(
    combined_cycle_alarm_evidence.tsv = extended_slide_data$cycle_alarm_evidence,
    combined_gene_transcript_quality.tsv = extended_slide_data$gene_quality,
    candidate_cycle_affected_genes.tsv = extended_slide_summary$candidates,
    candidate_cycle_affected_genes_affected_only.tsv = extended_slide_summary$candidates[extended_slide_summary$candidates$section_candidate_flag, , drop = FALSE],
    subset_full_qc_ranking.tsv = extended_slide_summary$ranking,
    subset_full_qc_rank_agreement.tsv = extended_slide_summary$rank_agreement,
    combined_spatial_qc.tsv = extended_slide_summary$spatial_qc,
    combined_spatial_hotspots.tsv = extended_slide_summary$spatial_hotspots,
    within_mouse_section_concordance.tsv = extended_slide_summary$concordance$summary,
    within_mouse_gene_concordance.tsv = extended_slide_summary$concordance$genes,
    extended_slide_qc_status.tsv = extended_slide_summary$status
  )
  paths <- vapply(names(tables), function(name) write_tsv(tables[[name]], file.path(output_dir, name), project_root), character(1))
  pdf_path <- file.path(figure_dir, "scwat_extended_qc_diagnostics.pdf")
  grDevices::pdf(pdf_path, width = 8, height = 5.5, onefile = TRUE)
  on.exit(if (grDevices::dev.cur() > 1L) grDevices::dev.off(), add = TRUE)
  for (plot in plots) print(plot)
  grDevices::dev.off()
  all_paths <- c(unname(paths), pdf_path)
  if (!validate_extended_slide_qc_artifacts(run_root)) stop("Extended slide artifact validation failed.", call. = FALSE)
  all_paths
}

# Purpose: Validate extended slide qc artifacts.
# Inputs: required: run_root; optional/defaulted: stop_on_error.
# Output: Returns validation evidence/TRUE (or the validated value) and stops with an informative error when the contract fails.
validate_extended_slide_qc_artifacts <- function(run_root, stop_on_error = FALSE) {
  fail <- function(message) {
    if (isTRUE(stop_on_error)) stop(message, call. = FALSE)
    FALSE
  }
  output_dir <- file.path(run_root, "slide_summary")
  required <- extended_slide_required_artifacts()
  absent <- required[!file.exists(file.path(output_dir, required))]
  if (length(absent)) return(fail(sprintf("Missing extended slide artifacts: %s", paste(absent, collapse = ", "))))
  status <- tryCatch(utils::read.delim(file.path(output_dir, "extended_slide_qc_status.tsv"), check.names = FALSE, stringsAsFactors = FALSE), error = identity)
  if (inherits(status, "error") || nrow(status) != 1L || status$sections[[1]] != 4L || status$cells_deleted[[1]] != 0L) return(fail("Invalid extended slide QC status."))
  concordance <- tryCatch(utils::read.delim(file.path(output_dir, "within_mouse_section_concordance.tsv"), check.names = FALSE, stringsAsFactors = FALSE), error = identity)
  if (inherits(concordance, "error") || nrow(concordance) != 2L) return(fail("Extended slide concordance must contain two mouse pairs."))
  TRUE
}

# Purpose: Section required artifacts.
# Inputs: required: region_id.
# Output: Returns a deterministic scalar, vector, path, status, or empty-schema object used by downstream functions.
section_required_artifacts <- function(region_id) {
  c(
    "configuration.tsv", "section_manifest.tsv", "environment_preflight.tsv", "file_inventory.tsv",
    "integrity_summary.tsv", "feature_type_summary.tsv", "panel_reconciliation.tsv", "analysis_alerts.tsv",
    "qc_thresholds.tsv", "qc_summary.tsv", "cell_qc_metadata.tsv.gz",
    paste0(region_id, ".phase0_2_qc.rds"), "section_readiness_gates.tsv", "sessionInfo.txt",
    file.path("figures", paste0(region_id, "_qc_overview.pdf")),
    file.path("figures", paste0(region_id, "_counts.png")),
    file.path("figures", paste0(region_id, "_features.png")),
    file.path("figures", paste0(region_id, "_area.png")),
    file.path("figures", paste0(region_id, "_spatial.png"))
  )
}

# Purpose: Write section artifacts.
# Inputs: required: project_root, output_dir, region_id, configuration, manifest, environment, inventory, integrity, feature_type_summary, panel_reconciliation, alarms, qc, counts, features; optional/defaulted: strict_mode.
# Output: Writes validated artifact file(s) and returns their path(s) invisibly or as a named path list.
write_section_artifacts <- function(project_root, output_dir, region_id, configuration, manifest, environment,
                                    inventory, integrity, feature_type_summary, panel_reconciliation, alarms,
                                    qc, counts, features, strict_mode = FALSE) {
  assert_path_within(project_root, output_dir)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  gates <- calculate_readiness_gates(inventory, integrity, panel_reconciliation, alarms, manifest)
  if (isTRUE(strict_mode) && gates$status[gates$gate == "overall"] == "HOLD") stop(sprintf("Strict mode stopped %s because readiness is HOLD.", region_id), call. = FALSE)
  tables <- list(
    configuration.tsv = configuration, section_manifest.tsv = manifest, environment_preflight.tsv = environment,
    file_inventory.tsv = inventory, integrity_summary.tsv = integrity, feature_type_summary.tsv = feature_type_summary,
    panel_reconciliation.tsv = panel_reconciliation, analysis_alerts.tsv = alarms,
    qc_thresholds.tsv = qc$thresholds, qc_summary.tsv = qc$summary, section_readiness_gates.tsv = gates
  )
  table_paths <- vapply(names(tables), function(name) write_tsv(tables[[name]], file.path(output_dir, name), project_root), character(1))
  cell_path <- write_gz_tsv(qc$cell_metadata, file.path(output_dir, "cell_qc_metadata.tsv.gz"), project_root)
  rds_path <- file.path(output_dir, paste0(region_id, ".phase0_2_qc.rds"))
  saveRDS(list(counts = counts, cells = qc$cell_metadata, features = features, region_id = region_id, raw_counts_preserved = TRUE), rds_path, compress = FALSE)
  plot_paths <- save_section_plots(plot_section_qc(qc$cell_metadata, region_id), file.path(output_dir, "figures"), region_id, project_root)
  session_path <- file.path(output_dir, "sessionInfo.txt")
  capture.output(sessionInfo(), file = session_path)
  paths <- c(unname(table_paths), cell_path, rds_path, plot_paths, session_path)
  if (!validate_section_artifacts(output_dir, region_id)) stop("Section artifact reload validation failed.", call. = FALSE)
  paths
}

# Purpose: Validate section artifacts.
# Inputs: required: output_dir, region_id.
# Output: Returns validation evidence/TRUE (or the validated value) and stops with an informative error when the contract fails.
validate_section_artifacts <- function(output_dir, region_id) {
  require_package("Matrix")
  paths <- file.path(output_dir, section_required_artifacts(region_id))
  if (!all(file.exists(paths))) return(FALSE)
  object <- readRDS(file.path(output_dir, paste0(region_id, ".phase0_2_qc.rds")))
  inherits(object$counts, "sparseMatrix") && ncol(object$counts) == nrow(object$cells) && identical(colnames(object$counts), object$cells$cell_id)
}

# Purpose: Slide section required files.
# Inputs: none.
# Output: Returns the derived R object described by the function name; no files are written unless an explicit output path is an input.
slide_section_required_files <- function() {
  c("qc_summary.tsv", "qc_thresholds.tsv", "section_readiness_gates.tsv", "analysis_alerts.tsv", "cell_qc_metadata.tsv.gz")
}

# Purpose: Validate four section outputs.
# Inputs: required: run_root; optional/defaulted: expected_regions.
# Output: Returns validation evidence/TRUE (or the validated value) and stops with an informative error when the contract fails.
validate_four_section_outputs <- function(run_root, expected_regions = paste0("Region_", 1:4)) {
  sections_root <- file.path(run_root, "sections")
  dirs <- list.dirs(sections_root, recursive = FALSE, full.names = TRUE)
  region_ids <- basename(dirs)
  valid_dirs <- region_ids %in% expected_regions
  dirs <- dirs[valid_dirs]; region_ids <- region_ids[valid_dirs]
  if (length(dirs) != 4L || !setequal(region_ids, expected_regions) || anyDuplicated(region_ids)) {
    stop(sprintf("Expected exactly four unique section outputs (%s).", paste(expected_regions, collapse = ", ")), call. = FALSE)
  }
  order_index <- match(expected_regions, region_ids)
  dirs <- dirs[order_index]; region_ids <- region_ids[order_index]
  missing <- unlist(lapply(seq_along(dirs), function(index) {
    paths <- file.path(dirs[[index]], slide_section_required_files())
    absent <- paths[!file.exists(paths)]
    if (!length(absent)) character() else paste(region_ids[[index]], basename(absent), sep = "/")
  }))
  if (length(missing)) stop(sprintf("Missing slide-summary inputs: %s", paste(missing, collapse = ", ")), call. = FALSE)
  data.frame(region_id = region_ids, section_output_dir = normalizePath(dirs, winslash = "/", mustWork = TRUE), stringsAsFactors = FALSE)
}

# Purpose: Read slide qc outputs.
# Inputs: required: run_root; optional/defaulted: expected_regions.
# Output: Returns parsed, validated R data (vector, data frame, sparse matrix bundle, or named list according to the input format).
read_slide_qc_outputs <- function(run_root, expected_regions = paste0("Region_", 1:4)) {
  coverage <- validate_four_section_outputs(run_root, expected_regions)
  read_one <- function(filename, gzipped = FALSE) {
    tables <- lapply(seq_len(nrow(coverage)), function(index) {
      path <- file.path(coverage$section_output_dir[[index]], filename)
      table <- if (gzipped) utils::read.delim(gzfile(path), check.names = FALSE) else utils::read.delim(path, check.names = FALSE)
      if (!"region_id" %in% names(table)) table$region_id <- rep(coverage$region_id[[index]], nrow(table))
      table
    })
    do.call(rbind, tables)
  }
  list(
    coverage = coverage,
    qc_summary = read_one("qc_summary.tsv"), thresholds = read_one("qc_thresholds.tsv"),
    gates = read_one("section_readiness_gates.tsv"), alarms = read_one("analysis_alerts.tsv"),
    cell_metadata = read_one("cell_qc_metadata.tsv.gz", gzipped = TRUE)
  )
}

# Purpose: Summarise slide qc.
# Inputs: required: slide_data.
# Output: Returns computed QC evidence as a vector, data frame, or named summary list; it does not modify raw input files.
summarise_slide_qc <- function(slide_data) {
  summary <- slide_data$qc_summary
  summary$core_pass_fraction <- ifelse(summary$input_cells > 0, summary$core_qc_pass / summary$input_cells, NA_real_)
  summary$review_fraction <- ifelse(summary$input_cells > 0, summary$review_flagged / summary$input_cells, NA_real_)
  overall_rows <- slide_data$gates[slide_data$gates$gate == "overall", , drop = FALSE]
  overall_status <- if (nrow(overall_rows)) overall_readiness(overall_rows$status) else "HOLD"
  readiness <- data.frame(
    region_id = overall_rows$region_id, status = overall_rows$status,
    biological_interpretation_allowed = overall_rows$status == "PASS", stringsAsFactors = FALSE
  )
  list(section_summary = summary, readiness = readiness, overall_status = overall_status)
}

# Purpose: Plot slide qc.
# Inputs: required: slide_data, slide_summary.
# Output: Returns a ggplot object or named list of plots; plotting does not mutate the input object.
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

# Purpose: Write slide qc artifacts.
# Inputs: required: project_root, run_root, slide_data, slide_summary, slide_plots.
# Output: Writes validated artifact file(s) and returns their path(s) invisibly or as a named path list.
write_slide_qc_artifacts <- function(project_root, run_root, slide_data, slide_summary, slide_plots) {
  assert_path_within(project_root, run_root)
  output_dir <- file.path(run_root, "slide_summary")
  figure_dir <- file.path(output_dir, "figures")
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  table_map <- list(
    combined_qc_summary.tsv = slide_summary$section_summary,
    combined_qc_thresholds.tsv = slide_data$thresholds,
    combined_readiness.tsv = slide_summary$readiness,
    combined_analysis_alerts.tsv = slide_data$alarms
  )
  paths <- vapply(names(table_map), function(name) write_tsv(table_map[[name]], file.path(output_dir, name), project_root), character(1))
  paths <- c(paths, write_gz_tsv(slide_data$cell_metadata, file.path(output_dir, "combined_cell_qc_metadata.tsv.gz"), project_root))
  rds_path <- file.path(output_dir, "slide_qc_summary.rds")
  saveRDS(list(data = slide_data, summary = slide_summary, raw_counts_in_section_objects = TRUE), rds_path, compress = FALSE)
  pdf_path <- file.path(figure_dir, "scwat_slide_qc_figures.pdf")
  grDevices::pdf(pdf_path, width = 8, height = 5.5, onefile = TRUE)
  for (plot in slide_plots) print(plot)
  grDevices::dev.off()
  png_paths <- vapply(names(slide_plots), function(name) {
    path <- file.path(figure_dir, paste0("slide_", name, ".png"))
    ggplot2::ggsave(path, slide_plots[[name]], width = 8, height = 5.5, units = "in", dpi = 300, bg = "white")
    path
  }, character(1))
  session_path <- file.path(output_dir, "sessionInfo.txt"); capture.output(sessionInfo(), file = session_path)
  status_path <- write_tsv(data.frame(scope = "scWAT_slide", status = slide_summary$overall_status, sections = 4L, cells = nrow(slide_data$cell_metadata), cells_deleted = 0L), file.path(output_dir, "slide_qc_status.tsv"), project_root)
  all_paths <- c(unname(paths), rds_path, pdf_path, unname(png_paths), session_path, status_path)
  if (!all(file.exists(all_paths))) stop("Slide QC artifact validation failed.", call. = FALSE)
  all_paths
}

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

# Purpose: Region bundle to spatial seurat.
# Inputs: required: region_data; optional/defaulted: xenium_dir, mask, genes, project, assay, fov, include_cell_segmentation, include_nucleus_segmentation.
# Output: Returns the derived R object described by the function name; no files are written unless an explicit output path is an input.
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

# Purpose: Calculate marker scores.
# Inputs: required: object, marker_df, group_col; optional/defaulted: assay, layer, prefix, z_cap.
# Output: Returns computed QC evidence as a vector, data frame, or named summary list; it does not modify raw input files.
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

# Purpose: Score marker groups by cluster.
# Inputs: required: object, marker_df, group_col; optional/defaulted: cluster_col, assay, layer, score_prefix, detect_prefix, label_prefix, z_cap.
# Output: Returns scores or refined annotations aligned to the supplied cells/features; raw counts are unchanged.
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

# Purpose: Refine xenium celltypes.
# Inputs: required: object; optional/defaulted: reduction, dims, k, self_weight, wang_prefix, cluster_main_col, cluster_subtype_col, ontology_gap_max_immune_prob, wang_main_review_score, wang_main_review_margin, wang_subtype_review_margin, xenium_subtype_review_margin, verbose.
# Output: Returns scores or refined annotations aligned to the supplied cells/features; raw counts are unchanged.
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

# Purpose: Refine eosinophil identity.
# Inputs: required: object; optional/defaulted: assay, reduction, dims, k, self_weight, wang_main_eos_col, wang_subtype_eos_col, wang_subtype_prefix, xenium_eos_score_col, eos_core_genes, eos_support_genes, competitor_score_cols, wang_high, wang_support, min_core_high, min_core_probable, xenium_margin_high, verbose.
# Output: Returns scores or refined annotations aligned to the supplied cells/features; raw counts are unchanged.
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

# Purpose: Score eosinophil likeness.
# Inputs: required: reference, query; optional/defaulted: reference_group_col, eos_label, reference_sample_col, reference_assay, query_assay, core_markers, support_markers, tier3_marker, context_markers, reference_immune_labels, query_main_col, query_immune_labels, wang_predicted_col, wang_eos_score_col, min_core_eos_cells, min_core_eos_pct, min_pair_eos_cells, min_pair_eos_pct, min_competitor_cells, smoothing, allow_nonimmune_rescue, verbose.
# Output: Returns scores or refined annotations aligned to the supplied cells/features; raw counts are unchanged.
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

# Purpose: Plot eos with celltypes.
# Inputs: required: object, celltypes; optional/defaulted: subtype_col, eos_label, fov, eos_col, other_cols, background_col, background_size, other_size, eos_size, background_alpha, highlight_alpha, flip_xy.
# Output: Returns a ggplot object or named list of plots; plotting does not mutate the input object.
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
