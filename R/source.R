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
  required <- c("tissue", "region_id", "mouse_id", "side", "section_id", "biological_replicate_id", "technical_replicate_id", "metadata_status")
  issues <- character()
  missing <- setdiff(required, names(manifest))
  if (length(missing)) issues <- c(issues, paste0("missing_columns:", paste(missing, collapse = ",")))
  if (!length(missing)) {
    if (anyDuplicated(manifest$region_id)) issues <- c(issues, "duplicate_region_id")
    if (anyDuplicated(manifest$section_id)) issues <- c(issues, "duplicate_section_id")
    if (!setequal(manifest$region_id, expected_regions)) issues <- c(issues, "region_set_mismatch")
    if (any(!manifest$side %in% c("Left", "Right"))) issues <- c(issues, "invalid_side")
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
