options(stringsAsFactors = FALSE)

`%||%` <- function(x, y) {
  if (is.null(x) || !length(x)) return(y)
  first <- x[[1]]
  missing_scalar <- is.atomic(first) && length(first) == 1L && (is.na(first) || !nzchar(as.character(first)))
  if (missing_scalar) y else x
}

canonical_path <- function(path) {
  value <- normalizePath(path, winslash = "/", mustWork = FALSE)
  value <- sub("/+$", "", value)
  if (.Platform$OS.type == "windows") value <- tolower(value)
  value
}

assert_path_within <- function(project_root, candidate) {
  root <- canonical_path(project_root)
  path <- canonical_path(candidate)
  valid <- identical(path, root) || startsWith(path, paste0(root, "/"))
  if (!valid) stop(sprintf("Unsafe path outside project root: %s", path), call. = FALSE)
  invisible(TRUE)
}

validate_runtime_paths <- function(project_root, input_root, output_root, temp_root = tempdir()) {
  if (!dir.exists(project_root)) stop(sprintf("Project root does not exist: %s", project_root), call. = FALSE)
  if (!dir.exists(input_root)) stop(sprintf("Input root does not exist: %s", input_root), call. = FALSE)
  for (path in c(input_root, output_root, temp_root)) assert_path_within(project_root, path)
  invisible(TRUE)
}

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

xenium_required_files <- function() {
  c(
    "experiment.xenium", "metrics_summary.csv", "analysis_summary.html", "gene_panel.json",
    "cells.csv.gz", "cell_feature_matrix/features.tsv.gz",
    "cell_feature_matrix/barcodes.tsv.gz", "cell_feature_matrix/matrix.mtx.gz"
  )
}

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

read_gz_rows <- function(path, header = FALSE) {
  utils::read.delim(gzfile(path), header = header, quote = "", check.names = FALSE, stringsAsFactors = FALSE)
}

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

write_tsv <- function(x, path, project_root) {
  assert_path_within(project_root, path)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  invisible(path)
}

require_package <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) stop(sprintf("Required R package '%s' is unavailable.", package), call. = FALSE)
}

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

resolve_extended_qc_mode <- function(requested_mode = "AUTO", region_dir) {
  mode <- toupper(trimws(as.character(requested_mode)))
  allowed <- c("AUTO", "LOCAL_SUBSET", "FULL_HPC")
  if (length(mode) != 1L || !mode %in% allowed) {
    stop(sprintf("Extended QC mode must be one of: %s.", paste(allowed, collapse = ", ")), call. = FALSE)
  }
  if (mode != "AUTO") return(mode)
  if (file.exists(file.path(region_dir, "transcripts.parquet"))) "FULL_HPC" else "LOCAL_SUBSET"
}

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

validate_spatial_cells <- function(cells) {
  required <- c("cell_id", "x_centroid", "y_centroid", "qc_review_flag")
  missing <- setdiff(required, names(cells))
  if (length(missing)) stop(sprintf("Spatial cell table missing columns: %s", paste(missing, collapse = ", ")), call. = FALSE)
  if (anyDuplicated(cells$cell_id)) stop("Spatial cell IDs must be unique.", call. = FALSE)
  if (any(!is.finite(cells$x_centroid)) || any(!is.finite(cells$y_centroid))) stop("Spatial centroids must be finite.", call. = FALSE)
  invisible(TRUE)
}

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

calculate_knn_density <- function(cells, k = 15L, mode = "LOCAL_SUBSET") {
  neighbors <- knn_index_distance(cells, k, mode)
  radius <- neighbors$distance[, ncol(neighbors$distance)]
  radius[radius <= 0] <- min(radius[radius > 0], na.rm = TRUE)
  as.numeric(k) / (pi * radius^2)
}

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

quantile_distribution_distance <- function(a, b, probabilities = seq(0.01, 0.99, 0.01)) {
  a <- a[is.finite(a)]; b <- b[is.finite(b)]
  if (!length(a) || !length(b)) return(NA_real_)
  qa <- as.numeric(stats::quantile(a, probabilities, names = FALSE, type = 7))
  qb <- as.numeric(stats::quantile(b, probabilities, names = FALSE, type = 7))
  scale <- stats::median(c(a, b))
  if (!is.finite(scale) || scale == 0) scale <- 1
  mean(abs(qa - qb)) / abs(scale)
}

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

empty_alarm_table <- function() {
  data.frame(
    raw_value = logical(), formatted_value = character(), raised = logical(),
    title = character(), message = character(), level = character(), id = character(),
    stringsAsFactors = FALSE
  )
}

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

read_xenium_features <- function(path) {
  features <- read_gz_rows(path, header = FALSE)
  if (ncol(features) != 3L) stop(sprintf("Expected three columns in %s", path), call. = FALSE)
  names(features) <- c("feature_id", "feature_name", "feature_type")
  features
}

read_xenium_barcodes <- function(path) {
  scan(gzfile(path), what = character(), quiet = TRUE)
}

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

safe_quantile <- function(x, probability) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  unname(stats::quantile(x, probability, names = FALSE, na.rm = TRUE, type = 7))
}

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

calculate_xenium_cell_qc <- function(counts, cells, region_id) {
  require_package("Matrix")
  required <- c("cell_id", "total_counts", "control_probe_counts", "genomic_control_counts", "control_codeword_counts", "cell_area", "nucleus_count")
  missing <- setdiff(required, names(cells))
  if (length(missing)) stop(sprintf("Cell metadata missing QC columns: %s", paste(missing, collapse = ",")), call. = FALSE)
  if (!identical(colnames(counts), cells$cell_id)) stop("Count columns and cell metadata are not aligned.", call. = FALSE)
  n_count <- as.numeric(Matrix::colSums(counts))
  n_feature <- as.numeric(Matrix::colSums(counts > 0))
  count_bounds <- robust_interval(n_count, 3, 5, 1)
  feature_bounds <- robust_interval(n_feature, 3, 5, 1)
  area_bounds <- robust_interval(cells$cell_area, 5, 5, 0)
  control_count <- cells$control_probe_counts + cells$genomic_control_counts + cells$control_codeword_counts
  control_fraction <- ifelse(cells$total_counts > 0, control_count / cells$total_counts, 0)
  control_upper <- max(0.05, safe_quantile(control_fraction, 0.995))
  qc_core_pass <- n_count >= as.numeric(count_bounds["lower"]) & n_count <= as.numeric(count_bounds["upper"]) &
    n_feature >= as.numeric(feature_bounds["lower"]) & n_feature <= as.numeric(feature_bounds["upper"])
  nucleus_missing <- cells$nucleus_count == 0
  multiple_nuclei <- cells$nucleus_count > 1
  area_outlier <- cells$cell_area < as.numeric(area_bounds["lower"]) | cells$cell_area > as.numeric(area_bounds["upper"])
  high_control <- control_fraction > control_upper
  high_complexity <- n_count > as.numeric(count_bounds["upper"]) | n_feature > as.numeric(feature_bounds["upper"])
  segmentation_multiplet <- multiple_nuclei | (high_complexity & cells$cell_area > as.numeric(area_bounds["upper"]))
  out <- cells
  out$region_id <- region_id; out$nCount_Xenium <- n_count; out$nFeature_Xenium <- n_feature
  out$control_fraction_cell <- control_fraction; out$nucleus_missing_flag <- nucleus_missing
  out$multiple_nuclei_flag <- multiple_nuclei; out$cell_area_outlier_flag <- area_outlier
  out$high_control_flag <- high_control; out$segmentation_multiplet_flag <- segmentation_multiplet
  out$qc_core_pass <- qc_core_pass
  out$qc_review_flag <- nucleus_missing | segmentation_multiplet | area_outlier | high_control | !qc_core_pass
  thresholds <- rbind(
    data.frame(metric = "nCount_Xenium", lower = as.numeric(count_bounds["lower"]), upper = as.numeric(count_bounds["upper"]), value = as.numeric(count_bounds["median"]), method = count_bounds["method"]),
    data.frame(metric = "nFeature_Xenium", lower = as.numeric(feature_bounds["lower"]), upper = as.numeric(feature_bounds["upper"]), value = as.numeric(feature_bounds["median"]), method = feature_bounds["method"]),
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

overall_readiness <- function(status) {
  status <- toupper(as.character(status))
  if (any(status %in% c("FAIL", "HOLD", "BLOCKED"))) return("HOLD")
  if (any(status %in% c("PENDING", "WARN"))) return("PENDING")
  "PASS"
}

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

section_palette <- function() {
  c(Region_1 = "#3C5488", Region_2 = "#00A087", Region_3 = "#E64B35", Region_4 = "#F39B7F")
}

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
  out$primary_include <- as.logical(out$qc_core_pass) &
    !as.logical(out$segmentation_multiplet_flag) & !as.logical(out$high_control_flag)
  out$strict_include <- !as.logical(out$qc_review_flag)
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
  out$mask_rule_primary <- "qc_core_pass AND NOT segmentation_multiplet_flag AND NOT high_control_flag"
  out$mask_rule_strict <- "NOT qc_review_flag"
  out$mask_rule_hotspot_sensitivity <- "primary_include AND NOT Region_3 morphology-review hotspot cell"
  out$provenance <- provenance
  out
}

cell_style_theme <- function(base_size = 10) {
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

write_gz_tsv <- function(x, path, project_root) {
  assert_path_within(project_root, path)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  con <- gzfile(path, "wt"); on.exit(close(con), add = TRUE)
  utils::write.table(x, con, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  invisible(path)
}

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

add_evidence_provenance <- function(table, run_label, execution_mode, provenance,
                                    source_artifact, generated_utc) {
  table$run_label <- run_label
  table$execution_mode <- toupper(execution_mode)
  table$generated_utc <- generated_utc
  table$source_artifact <- source_artifact
  table$provenance <- provenance
  table
}

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
    primary_mask_rule = "qc_core_pass AND NOT segmentation_multiplet_flag AND NOT high_control_flag",
    strict_mask_rule = "NOT qc_review_flag",
    hotspot_sensitivity_mask_rule = "primary_include AND NOT Region_3 morphology-review hotspot cell",
    decision_rule = "Fixed evidence-only section decision approved 2026-08-15; masks annotate without deleting raw cells",
    stringsAsFactors = FALSE
  )
  add_evidence_provenance(
    out, run_label, execution_mode, provenance,
    "section cell and spatial QC artifacts", generated_utc
  )
}

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
  gene_reconciled <- identical(unname(counts[c("raw", "conservative", "provisional", "risk")]), c(479, 67, 245, 234))
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
      ggplot2::geom_col(width = 0.7) + ggplot2::labs(title = "Section-level candidate affected genes", subtitle = "Candidates are not confirmed; exact cycle identity requires 10x diagnostics", x = NULL, y = "Candidate genes", fill = "Evidence") + cell_style_theme(),
    ranking = ggplot2::ggplot(ranking_long, ggplot2::aes(region_id, review_fraction, colour = source, group = source)) +
      ggplot2::geom_line(linewidth = 0.7) + ggplot2::geom_point(size = 2) +
      ggplot2::scale_y_continuous(labels = function(x) paste0(round(100 * x, 1), "%")) +
      ggplot2::labs(title = "Subset versus full-data review burden", x = NULL, y = "Review-flag rate", colour = NULL) + cell_style_theme(),
    spatial = ggplot2::ggplot(global, ggplot2::aes(region_id, statistic, colour = region_id)) +
      ggplot2::geom_hline(yintercept = 0, colour = "#BDBDBD") + ggplot2::geom_point(size = 2) +
      ggplot2::scale_colour_manual(values = palette, drop = FALSE) +
      ggplot2::labs(title = "Global spatial clustering of QC review flags", subtitle = "Coordinate statistic; morphology labels require image review", x = NULL, y = "kNN clustering statistic") + cell_style_theme() + ggplot2::theme(legend.position = "none"),
    concordance = ggplot2::ggplot(concordance, ggplot2::aes(mouse_id, review_rate_difference, fill = concordance_status)) +
      ggplot2::geom_col(width = 0.65) + ggplot2::scale_fill_manual(values = c(CONCORDANT = "#4DAF4A", REVIEW = "#D73027", NOT_ESTIMABLE = "#BDBDBD"), drop = FALSE) +
      ggplot2::labs(title = "Within-mouse technical concordance", subtitle = "Advisory section-pair comparison", x = NULL, y = "Absolute review-rate difference", fill = NULL) + cell_style_theme()
  )
}

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

validate_section_artifacts <- function(output_dir, region_id) {
  require_package("Matrix")
  paths <- file.path(output_dir, section_required_artifacts(region_id))
  if (!all(file.exists(paths))) return(FALSE)
  object <- readRDS(file.path(output_dir, paste0(region_id, ".phase0_2_qc.rds")))
  inherits(object$counts, "sparseMatrix") && ncol(object$counts) == nrow(object$cells) && identical(colnames(object$counts), object$cells$cell_id)
}

slide_section_required_files <- function() {
  c("qc_summary.tsv", "qc_thresholds.tsv", "section_readiness_gates.tsv", "analysis_alerts.tsv", "cell_qc_metadata.tsv.gz")
}

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
