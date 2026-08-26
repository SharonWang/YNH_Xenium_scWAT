# ARCHIVAL SOURCE — NOT LOADED BY ACTIVE QC NOTEBOOKS
#
# Purpose: preserve functions with no active notebook, test, or support-script
# consumer as of 2026-08-26. This file is retained for provenance and possible
# later recovery; it must not be sourced by the Phase 0-2 QC notebooks.
#
# Inputs: source R/source.R first if an archived function is intentionally
# reactivated, because archived workflows may call shared active helpers.
# Output: function definitions only; sourcing this file writes no artifacts.

# Archived function: annotate_seurat_clusters
annotate_seurat_clusters <- function(object, marker_config, config) {
  assert_downstream_model_environment(config)
  scores <- cluster_marker_evidence(object, marker_config, config)
  clusters <- sort(unique(scores$cluster))
  support <- do.call(rbind, lapply(clusters, function(cluster) {
    candidate <- scores[
      scores$cluster == cluster & scores$supported_primary_marker_count >= 2L &
        is.finite(scores$primary_score), , drop = FALSE
    ]
    if (!nrow(candidate)) {
      return(data.frame(
        cluster = cluster, reference_label = paste0("Unresolved_", cluster),
        support_status = "INSUFFICIENT_ABSOLUTE_MARKER_EVIDENCE",
        primary_marker_count = 0L, marker_detection_fraction = 0,
        score_margin = NA_real_, exclusion_contradiction = NA_real_,
        risk_only_support = any(scores$risk_marker_count[scores$cluster == cluster] >= 2L),
        stringsAsFactors = FALSE
      ))
    }
    candidate <- candidate[order(candidate$primary_score, decreasing = TRUE), , drop = FALSE]
    best <- candidate[1, , drop = FALSE]
    second_score <- if (nrow(candidate) > 1L) candidate$primary_score[[2]] else 0
    margin <- best$primary_score[[1]] - second_score
    competing <- scores[
      scores$cluster == cluster & scores$reference_label != best$reference_label[[1]] &
        scores$supported_primary_marker_count >= 2L & is.finite(scores$primary_score), , drop = FALSE
    ]
    competing_score <- if (nrow(competing)) max(competing$primary_score) else 0
    contradiction <- if (best$primary_score[[1]] > 0) {
      max(0, competing_score) / best$primary_score[[1]]
    } else {
      Inf
    }
    supported <- margin >= as.numeric(config$marker_score_margin) &&
      contradiction <= as.numeric(config$exclusion_contradiction_max)
    data.frame(
      cluster = cluster,
      reference_label = if (supported) best$reference_label[[1]] else paste0("Unresolved_", cluster),
      support_status = if (supported) "SUPPORTED" else "AMBIGUOUS_OR_CONTRADICTORY_MARKER_EVIDENCE",
      primary_marker_count = best$supported_primary_marker_count[[1]],
      marker_detection_fraction = best$marker_detection_fraction[[1]],
      score_margin = margin, exclusion_contradiction = contradiction,
      risk_only_support = FALSE, stringsAsFactors = FALSE
    )
  }))
  label_map <- setNames(support$reference_label, support$cluster)
  object$reference_label <- unname(label_map[as.character(object$seurat_clusters)])
  list(object = object, marker_support = support, marker_scores = scores)
}

# Archived function: apply_manual_admission_review
apply_manual_admission_review <- function(admissions, review_path = NULL) {
  if (is.null(review_path) || !nzchar(review_path)) return(admissions)
  if (!file.exists(review_path)) stop(sprintf("Missing manual admission review: %s", review_path), call. = FALSE)
  review <- utils::read.delim(review_path, check.names = FALSE, stringsAsFactors = FALSE)
  required <- c("region_id", "review_decision", "reviewer", "reviewed_utc", "rationale")
  if (length(setdiff(required, names(review))) || anyDuplicated(review$region_id)) {
    stop("Manual admission review requires unique region_id, review_decision, reviewer, reviewed_utc, and rationale.", call. = FALSE)
  }
  allowed <- c("ADMIT_AFTER_MARKER_MORPHOLOGY_REVIEW", "KEEP_SENSITIVITY_ONLY")
  if (any(!review$review_decision %in% allowed)) stop("Manual admission review contains an unsupported decision.", call. = FALSE)
  if (any(!nzchar(review$reviewer)) || any(!nzchar(review$reviewed_utc)) || any(!nzchar(review$rationale))) {
    stop("Manual admission review requires non-empty reviewer, reviewed_utc, and rationale.", call. = FALSE)
  }
  for (region in intersect(names(admissions), review$region_id)) {
    row <- review[review$region_id == region, , drop = FALSE]
    original <- admissions[[region]]$admission$decision
    if (row$review_decision == "ADMIT_AFTER_MARKER_MORPHOLOGY_REVIEW") {
      if (!identical(original, "REVIEW_POTENTIALLY_UNREPRESENTED")) {
        stop(sprintf("%s cannot be manually admitted from status %s.", region, original), call. = FALSE)
      }
      admissions[[region]]$admission$decision <- "ADMITTED_TO_CONSENSUS"
      admissions[[region]]$admission$summary$decision <- "ADMITTED_TO_CONSENSUS"
      admissions[[region]]$admission$summary$manual_review <- paste(row$reviewer, row$reviewed_utc, row$rationale, sep = "|")
    }
  }
  attr(admissions, "manual_review") <- review
  admissions
}


# Archived function: assert_downstream_model_environment
assert_downstream_model_environment <- function(config) {
  preflight <- downstream_model_preflight(config)
  if (any(preflight$status != "PASS")) {
    failed <- paste(preflight$package[preflight$status != "PASS"], preflight$status[preflight$status != "PASS"], sep = "=")
    stop(sprintf("Seurat/Harmony HPC model environment is not ready: %s", paste(failed, collapse = ", ")), call. = FALSE)
  }
  invisible(preflight)
}


# Archived function: assess_mapped_marker_coherence
assess_mapped_marker_coherence <- function(query, mapping, marker_config, config) {
  data <- get_seurat_normalized_data(query, "RNA")
  mapped <- mapping[mapping$mapping_status == "Mapped", , drop = FALSE]
  if (!nrow(mapped)) return(data.frame(label = character(), cells = integer(), coherent = logical()))
  rows <- lapply(sort(unique(mapped$predicted_label)), function(label) {
    ids <- mapped$cell_id[mapped$predicted_label == label]
    genes <- marker_config$gene[
      marker_config$cell_type == label & marker_config$use_policy == "PRIMARY_SUPPORT"
    ]
    genes <- intersect(genes, rownames(data))
    detected <- if (length(genes)) Matrix::rowMeans(data[genes, ids, drop = FALSE] > 0) else numeric()
    supported <- sum(detected >= config$marker_detection_min)
    data.frame(
      label = label, cells = length(ids), available_markers = length(genes),
      supported_markers = supported, coherent = supported >= 2L,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}


# Archived function: assess_region3_morphology_review
assess_region3_morphology_review <- function(hotspot_decision, review_path = NULL) {
  required_hotspots <- hotspot_decision[
    hotspot_decision$region_id == "Region_3" &
      hotspot_decision$hotspot_status == "MORPHOLOGY_REVIEW_REQUIRED", , drop = FALSE
  ]
  if (!nrow(required_hotspots)) {
    return(data.frame(
      gate_id = "region3_morphology_review", reviewed = 0L, required = 0L,
      gate_status = "PASS_NO_HOTSPOTS", details = "No FDR-positive Region 3 hotspots",
      stringsAsFactors = FALSE
    ))
  }
  if (is.null(review_path) || !nzchar(review_path) || !file.exists(review_path)) {
    return(data.frame(
      gate_id = "region3_morphology_review", reviewed = 0L, required = nrow(required_hotspots),
      gate_status = "REVIEW_MORPHOLOGY_PENDING",
      details = "DAPI/morphology/cell-boundary review is required before final exploratory release",
      stringsAsFactors = FALSE
    ))
  }
  review <- utils::read.delim(review_path, check.names = FALSE, stringsAsFactors = FALSE)
  required <- c("grid_id", "review_decision", "reviewer", "reviewed_utc")
  if (length(setdiff(required, names(review))) || anyDuplicated(review$grid_id)) {
    stop("Morphology review requires unique grid_id, review_decision, reviewer, and reviewed_utc columns.", call. = FALSE)
  }
  allowed <- c("VALID_ANATOMY", "ARTIFACT_EXCLUDE_IN_SENSITIVITY", "UNCERTAIN_RETAIN_PRIMARY")
  if (any(!review$review_decision %in% allowed)) stop("Morphology review contains an unsupported decision.", call. = FALSE)
  matched <- match(required_hotspots$grid_id, review$grid_id)
  complete <- !anyNA(matched)
  data.frame(
    gate_id = "region3_morphology_review", reviewed = sum(!is.na(matched)), required = nrow(required_hotspots),
    gate_status = if (complete) "PASS_REVIEW_COMPLETE" else "REVIEW_MORPHOLOGY_PENDING",
    details = if (complete) "All Region 3 hotspots reviewed; primary retention and sensitivity exclusions remain explicit" else "Some required hotspots lack review",
    stringsAsFactors = FALSE
  )
}


# Archived function: assign_spatial_mapping_folds
assign_spatial_mapping_folds <- function(cells, folds, seed) {
  if (!is.data.frame(cells) || !"cell_id" %in% names(cells)) stop("Calibration cells require cell_id.", call. = FALSE)
  folds <- as.integer(folds)
  if (all(c("x_centroid", "y_centroid") %in% names(cells))) {
    x_rank <- rank(cells$x_centroid, ties.method = "first")
    y_rank <- rank(cells$y_centroid, ties.method = "first")
    x_bin <- pmin(folds - 1L, floor((x_rank - 1L) / nrow(cells) * folds))
    y_bin <- pmin(folds - 1L, floor((y_rank - 1L) / nrow(cells) * folds))
    fold <- (x_bin + 2L * y_bin) %% folds + 1L
    method <- "SPATIALLY_BLOCKED_INTERNAL_CALIBRATION"
  } else {
    set.seed(as.integer(seed))
    fold <- sample(rep(seq_len(folds), length.out = nrow(cells)))
    method <- "STRATIFIED_RANDOM_FALLBACK_NO_COORDINATES"
  }
  data.frame(
    cell_id = as.character(cells$cell_id), calibration_fold = as.integer(fold),
    calibration_method = method,
    interpretation = "INTERNAL_REPRODUCIBILITY_NOT_BIOLOGICAL_ACCURACY",
    stringsAsFactors = FALSE
  )
}


# Archived function: build_eligible_consensus
build_eligible_consensus <- function(handoff, anchor, admissions, marker_config, config) {
  admitted <- names(admissions)[vapply(admissions, function(x) identical(x$admission$decision, "ADMITTED_TO_CONSENSUS"), logical(1))]
  eligible <- c("Region_3", admitted)
  primary_genes <- handoff$regions$Region_3$gene_sets$provisional_primary_features
  region_parts <- lapply(eligible, function(region) {
    bundle <- handoff$regions[[region]]
    cells <- bundle$cell_metadata[bundle$cell_metadata$primary_include, , drop = FALSE]
    list(counts = bundle$counts[, cells$cell_id, drop = FALSE], cells = cells)
  })
  counts <- do.call(cbind, lapply(region_parts, `[[`, "counts"))
  cells <- do.call(rbind, lapply(region_parts, `[[`, "cells"))
  consensus <- build_seurat_reference(counts, cells, primary_genes, config, "ELIGIBLE_UNCORRECTED_CONSENSUS", config$seed)
  annotation <- annotate_seurat_clusters(consensus, marker_config, config)
  consensus <- annotation$object
  consensus@misc$primary_reduction_policy <- "UNCORRECTED_PCA_PRIMARY"
  harmony_status <- "NOT_RUN_REGION3_ONLY"
  if (length(eligible) > 1L) {
    dims <- seq_len(min(config$n_pcs, ncol(Seurat::Embeddings(consensus, "pca"))))
    consensus <- harmony::RunHarmony(
      consensus, group.by.vars = "region_id", reduction.use = "pca", dims.use = dims,
      theta = config$harmony_theta, lambda = config$harmony_lambda, sigma = config$harmony_sigma,
      max.iter.harmony = config$harmony_max_iter, reference_values = config$harmony_reference_region,
      reduction.save = "harmony_sensitivity", verbose = FALSE
    )
    consensus@misc$harmony_policy <- "HARMONY_SENSITIVITY_ONLY"
    harmony_status <- "HARMONY_SENSITIVITY_ONLY"
  }
  expected_labels <- anchor$branches$primary_245$cell_labels
  observed <- data.frame(cell_id = colnames(consensus), consensus_label = as.character(consensus$reference_label))
  anchor_compare <- merge(expected_labels[, c("cell_id", "reference_label")], observed, by = "cell_id")
  anchor_agreement <- mean(anchor_compare$reference_label == anchor_compare$consensus_label)
  cluster_region <- table(consensus$seurat_clusters, consensus$region_id)
  dominant_fraction <- apply(cluster_region, 1L, function(x) max(x) / sum(x))
  cluster_size <- rowSums(cluster_region)
  dominant_review <- dominant_fraction >= config$section_cluster_dominance & cluster_size >= config$section_cluster_min_cells
  gates <- data.frame(
    gate_id = c("region3_label_preservation", "section_dominant_cluster_review"),
    observed = c(anchor_agreement, sum(dominant_review)),
    threshold = c(config$consensus_label_stability_min, 0),
    gate_status = c(
      if (anchor_agreement >= config$consensus_label_stability_min) "PASS" else "STOP",
      if (any(dominant_review)) "REVIEW_MORPHOLOGY_AND_MARKERS" else "PASS"
    ), stringsAsFactors = FALSE
  )
  hard_stop <- any(gates$gate_status == "STOP")
  final_reference <- if (hard_stop) anchor$branches$primary_245$object else consensus
  decision <- if (hard_stop) "FALLBACK_TO_REGION3" else if (any(grepl("REVIEW", gates$gate_status))) "CONSENSUS_REVIEW" else "CONSENSUS_PASS"
  list(
    eligible_regions = eligible, admitted_regions = admitted,
    consensus_object = consensus, final_reference = final_reference,
    decision = decision, gates = gates, marker_support = annotation$marker_support,
    harmony_status = harmony_status,
    interpretation = "Section dominance is reviewed, not assumed technical; Harmony is sensitivity-only"
  )
}


# Archived function: build_seurat_reference
build_seurat_reference <- function(counts, cells, genes, config, role, seed = config$seed) {
  assert_downstream_model_environment(config)
  if (!inherits(counts, "sparseMatrix")) stop("Reference counts must be a sparse Matrix.", call. = FALSE)
  if (!is.data.frame(cells) || !"cell_id" %in% names(cells) || anyDuplicated(cells$cell_id)) {
    stop("Reference cells must contain unique cell_id values.", call. = FALSE)
  }
  if (!length(genes) || anyDuplicated(genes) || length(setdiff(genes, rownames(counts)))) {
    stop("Reference genes must be unique and present in the count matrix.", call. = FALSE)
  }
  cell_ids <- as.character(cells$cell_id)
  if (length(setdiff(cell_ids, colnames(counts)))) stop("Reference cell IDs are absent from the count matrix.", call. = FALSE)
  if (length(cell_ids) < 3L || length(genes) < 3L) stop("Reference fitting requires at least three cells and genes.", call. = FALSE)
  model_counts <- counts[, cell_ids, drop = FALSE]
  metadata <- cells[match(cell_ids, cells$cell_id), , drop = FALSE]
  rownames(metadata) <- cell_ids
  set.seed(as.integer(seed))
  object <- Seurat::CreateSeuratObject(counts = model_counts, meta.data = metadata, project = role, min.cells = 0L, min.features = 0L)
  object <- Seurat::NormalizeData(
    object, normalization.method = "LogNormalize",
    scale.factor = as.numeric(config$normalization_scale_factor), verbose = FALSE
  )
  object <- Seurat::ScaleData(object, features = genes, verbose = FALSE)
  npcs <- min(as.integer(config$n_pcs), length(genes) - 1L, ncol(object) - 1L)
  object <- Seurat::RunPCA(object, features = genes, npcs = npcs, seed.use = as.integer(seed), verbose = FALSE)
  dims <- seq_len(npcs)
  object <- Seurat::FindNeighbors(object, reduction = "pca", dims = dims, k.param = min(as.integer(config$knn_k), ncol(object) - 1L), verbose = FALSE)
  object <- Seurat::FindClusters(
    object, resolution = as.numeric(config$leiden_resolution), algorithm = 4L,
    random.seed = as.integer(seed), verbose = FALSE
  )
  object <- Seurat::RunUMAP(
    object, reduction = "pca", dims = dims, seed.use = as.integer(seed),
    return.model = TRUE, reduction.name = "umap_pca", verbose = FALSE
  )
  object@misc$downstream_reference_contract <- list(
    role = role, feature_policy = primary_model_feature_policy(), feature_names = genes,
    raw_feature_names = rownames(model_counts), raw_gene_count = nrow(model_counts),
    n_pcs = npcs, seed = as.integer(seed), normalization_method = "LogNormalize",
    normalization_scale_factor = as.numeric(config$normalization_scale_factor),
    package_versions = setNames(
      vapply(c("Seurat", "SeuratObject"), function(package) as.character(utils::packageVersion(package)), character(1)),
      c("Seurat", "SeuratObject")
    )
  )
  object
}


# Archived function: calibrate_mapping_thresholds
calibrate_mapping_thresholds <- function(mapping, config) {
  required <- c("predicted_label", "prediction_confidence", "confidence_margin", "reference_distance")
  if (length(setdiff(required, names(mapping)))) stop("Calibration mapping table is incomplete.", call. = FALSE)
  labels <- c("__GLOBAL__", sort(unique(mapping$predicted_label)))
  rows <- lapply(labels, function(label) {
    x <- if (label == "__GLOBAL__") mapping else mapping[mapping$predicted_label == label, , drop = FALSE]
    enough <- label == "__GLOBAL__" || nrow(x) >= as.integer(config$mapping_min_label_cells)
    if (!enough) return(NULL)
    data.frame(
      reference_label = label, cells = nrow(x),
      confidence_min = as.numeric(stats::quantile(x$prediction_confidence, config$mapping_confidence_quantile, na.rm = TRUE)),
      margin_min = as.numeric(stats::quantile(x$confidence_margin, config$mapping_confidence_quantile, na.rm = TRUE)),
      distance_max = as.numeric(stats::quantile(x$reference_distance, config$mapping_distance_quantile, na.rm = TRUE)),
      calibration_scope = "INTERNAL_REPRODUCIBILITY_NOT_BIOLOGICAL_ACCURACY",
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
}


# Archived function: calibrate_region3_mapping
calibrate_region3_mapping <- function(bundle, anchor, config) {
  primary <- anchor$branches$primary_245
  cells <- bundle$cell_metadata[bundle$cell_metadata$primary_include, , drop = FALSE]
  folds <- assign_spatial_mapping_folds(cells, config$mapping_folds, config$seed)
  truth <- primary$cell_labels[, c("cell_id", "reference_label"), drop = FALSE]
  mappings <- list()
  for (fold in seq_len(as.integer(config$mapping_folds))) {
    test_ids <- folds$cell_id[folds$calibration_fold == fold]
    train_ids <- setdiff(primary$cell_labels$cell_id, test_ids)
    training_cells <- bundle$cell_metadata[match(train_ids, bundle$cell_metadata$cell_id), , drop = FALSE]
    reference <- build_seurat_reference(
      bundle$counts, training_cells, primary$feature_names, config,
      role = paste0("REGION3_CALIBRATION_FOLD_", fold), seed = config$seed + fold
    )
    reference$reference_label <- truth$reference_label[match(colnames(reference), truth$cell_id)]
    test_bundle <- bundle
    test_bundle$cell_metadata$calibration_test <- test_bundle$cell_metadata$cell_id %in% test_ids
    query <- prepare_seurat_query(test_bundle, "calibration_test", primary$feature_names, config, paste0("REGION3_QUERY_FOLD_", fold))
    mapped <- map_query_to_frozen_reference(reference, query, primary$feature_names, config)
    mapped$calibration_fold <- fold
    mapped$anchor_label <- truth$reference_label[match(mapped$cell_id, truth$cell_id)]
    mapped$label_reproducible <- mapped$predicted_label == mapped$anchor_label
    mappings[[fold]] <- mapped
  }
  mapping <- do.call(rbind, mappings)
  thresholds <- calibrate_mapping_thresholds(mapping, config)
  mapping <- classify_mapping_uncertainty(mapping, thresholds)
  list(
    folds = folds, mapping = mapping, thresholds = thresholds,
    summary = data.frame(
      cells = nrow(mapping), label_reproducibility = mean(mapping$label_reproducible),
      mapped_fraction = mean(mapping$mapping_status == "Mapped"),
      interpretation = "INTERNAL_REPRODUCIBILITY_NOT_BIOLOGICAL_ACCURACY",
      stringsAsFactors = FALSE
    )
  )
}


# Archived function: classify_mapping_uncertainty
classify_mapping_uncertainty <- function(mapping, thresholds) {
  global <- thresholds[thresholds$reference_label == "__GLOBAL__", , drop = FALSE]
  if (nrow(global) != 1L) stop("Mapping thresholds require exactly one global row.", call. = FALSE)
  matched <- match(mapping$predicted_label, thresholds$reference_label)
  confidence_min <- ifelse(is.na(matched), global$confidence_min, thresholds$confidence_min[matched])
  margin_min <- ifelse(is.na(matched), global$margin_min, thresholds$margin_min[matched])
  distance_max <- ifelse(is.na(matched), global$distance_max, thresholds$distance_max[matched])
  confidence_fail <- mapping$prediction_confidence < confidence_min
  margin_fail <- mapping$confidence_margin < margin_min
  distance_fail <- mapping$reference_distance > distance_max
  mapping$mapping_status <- ifelse(
    distance_fail & confidence_fail, "Potentially_unrepresented",
    ifelse(confidence_fail | margin_fail | distance_fail, "Uncertain", "Mapped")
  )
  mapping$final_label <- ifelse(mapping$mapping_status == "Mapped", mapping$predicted_label, mapping$mapping_status)
  mapping$confidence_threshold <- confidence_min
  mapping$margin_threshold <- margin_min
  mapping$distance_threshold <- distance_max
  mapping
}


# Archived function: cluster_marker_evidence
cluster_marker_evidence <- function(object, marker_config, config) {
  required <- c("cell_type", "gene", "use_policy")
  if (!inherits(object, "Seurat")) stop("Cluster marker evidence requires a Seurat object.", call. = FALSE)
  if (!is.data.frame(marker_config) || length(setdiff(required, names(marker_config)))) {
    stop("Canonical marker configuration is incomplete.", call. = FALSE)
  }
  cluster_id <- as.character(object$seurat_clusters)
  data <- get_seurat_normalized_data(object, "RNA")
  markers <- marker_config[marker_config$gene %in% rownames(data), , drop = FALSE]
  clusters <- sort(unique(cluster_id))
  types <- sort(unique(markers$cell_type))
  marker_genes <- unique(markers$gene)
  mean_matrix <- vapply(clusters, function(cluster) {
    Matrix::rowMeans(data[marker_genes, cluster_id == cluster, drop = FALSE])
  }, numeric(length(marker_genes)))
  detection_matrix <- vapply(clusters, function(cluster) {
    Matrix::rowMeans(data[marker_genes, cluster_id == cluster, drop = FALSE] > 0)
  }, numeric(length(marker_genes)))
  if (is.null(dim(mean_matrix))) mean_matrix <- matrix(mean_matrix, ncol = 1L)
  if (is.null(dim(detection_matrix))) detection_matrix <- matrix(detection_matrix, ncol = 1L)
  rownames(mean_matrix) <- rownames(detection_matrix) <- marker_genes
  colnames(mean_matrix) <- colnames(detection_matrix) <- clusters
  relative_matrix <- matrix(
    0, nrow = length(marker_genes), ncol = length(clusters),
    dimnames = list(marker_genes, clusters)
  )
  for (gene in marker_genes) {
    values <- mean_matrix[gene, ]
    deviation <- stats::sd(values)
    if (!is.na(deviation) && deviation > 0) relative_matrix[gene, ] <- (values - mean(values)) / deviation
  }
  rows <- list()
  for (cluster in clusters) {
    cluster_mean <- mean_matrix[, cluster]
    cluster_detection <- detection_matrix[, cluster]
    for (cell_type in types) {
      primary <- markers$gene[
        markers$cell_type == cell_type & markers$use_policy == "PRIMARY_SUPPORT"
      ]
      risk <- markers$gene[
        markers$cell_type == cell_type & markers$use_policy == "VALIDATION_ONLY_TECHNICAL_RISK"
      ]
      primary_detected <- primary[
        cluster_detection[primary] >= as.numeric(config$marker_detection_min) &
          cluster_mean[primary] >= as.numeric(config$marker_average_log_expression_min)
      ]
      rows[[length(rows) + 1L]] <- data.frame(
        cluster = cluster, reference_label = cell_type,
        primary_marker_count = length(primary),
        supported_primary_marker_count = length(primary_detected),
        marker_detection_fraction = if (length(primary)) length(primary_detected) / length(primary) else 0,
        primary_score = if (length(primary_detected)) mean(relative_matrix[primary_detected, cluster]) else NA_real_,
        risk_marker_count = length(risk),
        risk_validation_score = if (length(risk)) mean(cluster_mean[risk]) else NA_real_,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}


# Archived function: downstream_model_dependencies
downstream_model_dependencies <- function() {
  data.frame(
    package = c("Seurat", "SeuratObject", "harmony", "leidenbase"),
    minimum_version = c("4.3.0", "4.1.0", "1.2.0", "0.1.0"),
    role = c("reference_clustering_and_mapping", "sparse_object_container", "eligible_consensus_integration", "leiden_clustering"),
    stringsAsFactors = FALSE
  )
}


# Archived function: downstream_model_preflight
downstream_model_preflight <- function(config) {
  dependencies <- downstream_model_dependencies()
  dependencies$minimum_version[dependencies$package == "Seurat"] <- as.character(config$seurat_min_version)
  dependencies$minimum_version[dependencies$package == "harmony"] <- as.character(config$harmony_min_version)
  dependencies$available <- vapply(dependencies$package, requireNamespace, logical(1), quietly = TRUE)
  dependencies$installed_version <- vapply(seq_len(nrow(dependencies)), function(index) {
    if (!dependencies$available[[index]]) return(NA_character_)
    as.character(utils::packageVersion(dependencies$package[[index]]))
  }, character(1))
  dependencies$version_ok <- vapply(seq_len(nrow(dependencies)), function(index) {
    dependencies$available[[index]] &&
      utils::compareVersion(dependencies$installed_version[[index]], dependencies$minimum_version[[index]]) >= 0L
  }, logical(1))
  dependencies$status <- ifelse(
    !dependencies$available,
    "SKIP_LOCAL_MODEL_TEST_HPC_REQUIRED",
    ifelse(dependencies$version_ok, "PASS", "FAIL_HPC_PACKAGE_VERSION")
  )
  dependencies
}


# Archived function: evaluate_section_admission
evaluate_section_admission <- function(region_id, primary_mapping, strict_mapping,
                                       conservative_mapping, marker_coherence, config) {
  shared_strict <- merge(
    primary_mapping[, c("cell_id", "predicted_label")],
    strict_mapping[, c("cell_id", "predicted_label")], by = "cell_id", suffixes = c("_primary", "_strict")
  )
  shared_conservative <- merge(
    primary_mapping[, c("cell_id", "predicted_label")],
    conservative_mapping[, c("cell_id", "predicted_label")], by = "cell_id", suffixes = c("_245", "_67")
  )
  accepted <- mean(primary_mapping$mapping_status == "Mapped")
  unrepresented <- mean(primary_mapping$mapping_status == "Potentially_unrepresented")
  strict_agreement <- if (nrow(shared_strict)) mean(shared_strict$predicted_label_primary == shared_strict$predicted_label_strict) else NA_real_
  conservative_agreement <- if (nrow(shared_conservative)) mean(shared_conservative$predicted_label_245 == shared_conservative$predicted_label_67) else NA_real_
  marker_fraction <- if (nrow(marker_coherence)) {
    sum(marker_coherence$cells * marker_coherence$coherent) / sum(marker_coherence$cells)
  } else 0
  gates <- data.frame(
    gate_id = c("mapping_coverage", "canonical_marker_coherence", "primary_strict_mapping_stability", "gene245_vs_gene67_mapping"),
    observed = c(accepted, marker_fraction, strict_agreement, conservative_agreement),
    threshold = c(config$mapping_accepted_fraction_min, config$mapping_marker_coherence_min,
                  config$mapping_label_stability_min, config$gene_sensitivity_concordance_min),
    stringsAsFactors = FALSE
  )
  gates$gate_status <- ifelse(!is.na(gates$observed) & gates$observed >= gates$threshold, "PASS", "STOP")
  review_novel <- unrepresented >= config$mapping_unrepresented_fraction_review
  decision <- if (any(gates$gate_status == "STOP")) {
    "SENSITIVITY_ONLY"
  } else if (review_novel) {
    "REVIEW_POTENTIALLY_UNREPRESENTED"
  } else {
    "ADMITTED_TO_CONSENSUS"
  }
  list(
    region_id = region_id, decision = decision, gates = gates,
    summary = data.frame(
      region_id = region_id, decision = decision, mapped_fraction = accepted,
      potentially_unrepresented_fraction = unrepresented,
      marker_coherence_fraction = marker_fraction,
      strict_label_agreement = strict_agreement,
      conservative_label_agreement = conservative_agreement,
      interpretation = "Poor mapping may be technical or biologically unrepresented; review coherent out-of-reference cells",
      stringsAsFactors = FALSE
    )
  )
}


# Archived function: finalize_downstream_release
finalize_downstream_release <- function(anchor_validation, consensus, eos_result,
                                        region4_result, admissions, morphology_gate = NULL) {
  primary_hard_pass <- anchor_validation$decision != "STOP_ANCHOR" &&
    consensus$decision != "FALLBACK_TO_REGION3"
  primary_status <- if (primary_hard_pass) {
    if (anchor_validation$decision == "REVIEW_ANCHOR" || consensus$decision == "CONSENSUS_REVIEW") "REVIEW" else "PASS"
  } else if (anchor_validation$decision != "STOP_ANCHOR" && consensus$decision == "FALLBACK_TO_REGION3") {
    "PASS_REGION3_FALLBACK"
  } else {
    "STOP"
  }
  if (!is.null(morphology_gate) && any(grepl("PENDING", morphology_gate$gate_status)) &&
      primary_status != "STOP") primary_status <- "REVIEW"
  eos_status <- if (nrow(eos_result$stability) && all(eos_result$stability$gate_status == "PASS_EOS")) "PASS" else "STOP"
  region4_status <- if (all(region4_result$gate$gate_status == "PASS_SENSITIVITY")) "PASS" else "STOP"
  data.frame(
    release_domain = c(
      "PRIMARY_EXPLORATORY_REFERENCE_RELEASE",
      "EOS_DESCRIPTIVE_ANALYSIS_RELEASE",
      "REGION4_SENSITIVITY_RELEASE"
    ),
    release_status = c(primary_status, eos_status, region4_status),
    evidence = c(
      paste("anchor=", anchor_validation$decision, ";consensus=", consensus$decision, sep = ""),
      paste(eos_result$stability$gate_id, eos_result$stability$gate_status, collapse = ";"),
      paste(region4_result$gate$gate_id, region4_result$gate$gate_status, collapse = ";")
    ),
    inference_scope = c(
      "Exploratory broad-cell reference; one clean anchor section and two biological mice",
      "Descriptive Eosinophil identity/state only; no population-level inference",
      "Sensitivity-only mapped labels; Region 4 never contributes to reference fitting"
    ),
    generated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    stringsAsFactors = FALSE
  )
}


# Archived function: find_primary_markers
find_primary_markers <- function(object, genes, config) {
  assert_downstream_model_environment(config)
  if (!inherits(object, "Seurat")) stop("Primary marker input must be a Seurat object.", call. = FALSE)
  genes <- intersect(as.character(genes), rownames(object))
  if (!length(genes)) stop("No primary marker genes are present in the Seurat object.", call. = FALSE)
  object <- join_seurat_layers_if_needed(object, "RNA")
  Seurat::Idents(object) <- "seurat_clusters"
  markers <- Seurat::FindAllMarkers(
    object, assay = "RNA", features = genes, only.pos = TRUE,
    min.pct = 0.05, logfc.threshold = 0.1, verbose = FALSE
  )
  if (nrow(markers) && any(!markers$gene %in% genes)) stop("Primary marker result escaped the frozen feature set.", call. = FALSE)
  markers
}



# Archived function: fit_region3_anchor_branches
fit_region3_anchor_branches <- function(bundle, marker_config, config) {
  if (!identical(bundle$region_id, "Region_3") || !identical(bundle$section_status, "PRIMARY")) {
    stop("Region 3 PRIMARY bundle is required for anchor fitting.", call. = FALSE)
  }
  definitions <- region3_anchor_branch_definitions()
  branches <- setNames(vector("list", nrow(definitions)), definitions$branch_id)
  for (index in seq_len(nrow(definitions))) {
    definition <- definitions[index, , drop = FALSE]
    mask <- definition$mask_name[[1]]
    gene_set <- definition$gene_set_name[[1]]
    cells <- bundle$cell_metadata[bundle$cell_metadata[[mask]], , drop = FALSE]
    genes <- bundle$gene_sets[[gene_set]]
    model <- build_seurat_reference(
      bundle$counts, cells, genes, config,
      role = paste0("REGION3_", toupper(definition$branch_id[[1]])), seed = config$seed
    )
    markers <- find_primary_markers(model, genes, config)
    annotation <- annotate_seurat_clusters(model, marker_config, config)
    labels <- data.frame(
      cell_id = colnames(annotation$object),
      cluster = as.character(annotation$object$seurat_clusters),
      reference_label = as.character(annotation$object$reference_label),
      stringsAsFactors = FALSE
    )
    branches[[definition$branch_id[[1]]]] <- list(
      object = annotation$object, cell_labels = labels, primary_markers = markers,
      marker_support = annotation$marker_support, marker_scores = annotation$marker_scores,
      mask_name = mask, gene_set_name = gene_set, feature_names = genes,
      scientific_role = definition$scientific_role[[1]]
    )
  }
  validation <- validate_anchor_branches(
    branches$primary_245, branches$strict_245, branches$hotspot_245,
    branches$conservative_67, branches$complete_479, marker_config, config
  )
  list(branch_definitions = definitions, branches = branches, validation = validation)
}


# Archived function: get_seurat_normalized_data
get_seurat_normalized_data <- function(object, assay = "RNA") {
  if (utils::packageVersion("SeuratObject") >= "5.0.0" && exists("LayerData", envir = asNamespace("SeuratObject"), inherits = FALSE)) {
    return(SeuratObject::LayerData(object, assay = assay, layer = "data"))
  }
  Seurat::GetAssayData(object, assay = assay, slot = "data")
}




# Archived function: identify_eosinophils_independently
identify_eosinophils_independently <- function(reference, marker_config, config) {
  data <- get_seurat_normalized_data(reference, "RNA")
  eos_genes <- marker_config$gene[
    marker_config$cell_type == "Eosinophil" & marker_config$use_policy == "PRIMARY_SUPPORT"
  ]
  eos_genes <- intersect(unique(eos_genes), rownames(data))
  exclusion_types <- c("Macrophage", "Neutrophil", "Mast_cell")
  exclusion <- marker_config$gene[
    marker_config$cell_type %in% exclusion_types & marker_config$use_policy == "PRIMARY_SUPPORT"
  ]
  exclusion <- intersect(unique(exclusion), rownames(data))
  if (!length(eos_genes)) stop("No eosinophil identity markers are available in the frozen panel.", call. = FALSE)
  eos_count <- Matrix::colSums(data[eos_genes, , drop = FALSE] > 0)
  exclusion_count <- if (length(exclusion)) Matrix::colSums(data[exclusion, , drop = FALSE] > 0) else rep(0, ncol(data))
  eos_fraction <- as.numeric(eos_count) / length(eos_genes)
  exclusion_fraction <- if (length(exclusion)) as.numeric(exclusion_count) / length(exclusion) else rep(0, ncol(data))
  label <- as.character(reference$reference_label)
  validated <- label == "Eosinophil" &
    eos_count >= as.integer(config$eos_identity_min_markers) &
    eos_fraction >= as.numeric(config$eos_identity_detection_min) &
    exclusion_fraction < eos_fraction
  data.frame(
    cell_id = colnames(reference), reference_label = label,
    eos_identity_marker_count = as.integer(eos_count),
    eos_identity_marker_fraction = eos_fraction,
    exclusion_marker_count = as.integer(exclusion_count),
    exclusion_marker_fraction = exclusion_fraction,
    validated_eosinophil = validated,
    eos_identity_rule = paste0(
      "Reference-label confirmation using canonical eosinophil identity markers; requires >=",
      as.integer(config$eos_identity_min_markers), " detected markers, eosinophil marker fraction >=",
      as.numeric(config$eos_identity_detection_min),
      ", and eosinophil marker fraction greater than pooled macrophage/neutrophil/mast exclusion fraction; state genes not used for selection"
    ),
    stringsAsFactors = FALSE
  )
}


# Archived function: join_seurat_layers_if_needed
join_seurat_layers_if_needed <- function(object, assay = "RNA") {
  if (utils::packageVersion("SeuratObject") >= "5.0.0" && exists("JoinLayers", envir = asNamespace("SeuratObject"), inherits = FALSE)) {
    object <- SeuratObject::JoinLayers(object, assay = assay)
  }
  object
}


# Archived function: label_agreement_on_major
label_agreement_on_major <- function(primary_labels, branch_labels, major_labels) {
  joined <- merge(
    primary_labels[, c("cell_id", "reference_label")],
    branch_labels[, c("cell_id", "reference_label")], by = "cell_id",
    suffixes = c("_primary", "_branch"), sort = FALSE
  )
  joined <- joined[joined$reference_label_primary %in% major_labels, , drop = FALSE]
  if (!nrow(joined)) return(NA_real_)
  mean(joined$reference_label_primary == joined$reference_label_branch)
}


# Archived function: map_and_evaluate_conditional_region
map_and_evaluate_conditional_region <- function(bundle, anchor, calibration,
                                                marker_config, config) {
  if (!bundle$region_id %in% c("Region_1", "Region_2")) stop("Only Region 1 or 2 can enter conditional admission.", call. = FALSE)
  primary_reference <- anchor$branches$primary_245$object
  conservative_reference <- anchor$branches$conservative_67$object
  primary_query <- prepare_seurat_query(bundle, "primary_include", anchor$branches$primary_245$feature_names, config, paste0(bundle$region_id, "_PRIMARY_QUERY"))
  strict_query <- prepare_seurat_query(bundle, "strict_include", anchor$branches$primary_245$feature_names, config, paste0(bundle$region_id, "_STRICT_QUERY"))
  conservative_query <- prepare_seurat_query(bundle, "primary_include", anchor$branches$conservative_67$feature_names, config, paste0(bundle$region_id, "_CONSERVATIVE_QUERY"))
  primary_mapping <- map_query_to_frozen_reference(primary_reference, primary_query, anchor$branches$primary_245$feature_names, config, calibration$thresholds)
  strict_mapping <- map_query_to_frozen_reference(primary_reference, strict_query, anchor$branches$primary_245$feature_names, config, calibration$thresholds)
  conservative_mapping <- map_query_to_frozen_reference(conservative_reference, conservative_query, anchor$branches$conservative_67$feature_names, config)
  coherence <- assess_mapped_marker_coherence(primary_query, primary_mapping, marker_config, config)
  admission <- evaluate_section_admission(bundle$region_id, primary_mapping, strict_mapping, conservative_mapping, coherence, config)
  list(
    region_id = bundle$region_id, primary_query = primary_query,
    primary_mapping = primary_mapping, strict_mapping = strict_mapping,
    conservative_mapping = conservative_mapping, marker_coherence = coherence,
    admission = admission
  )
}


# Archived function: map_query_to_frozen_reference
map_query_to_frozen_reference <- function(reference, query, features, config,
                                          thresholds = NULL) {
  assert_downstream_model_environment(config)
  if (!inherits(reference, "Seurat") || !inherits(query, "Seurat")) stop("Mapping requires Seurat reference and query objects.", call. = FALSE)
  features <- intersect(features, intersect(rownames(reference), rownames(query)))
  dims <- seq_len(min(config$n_pcs, ncol(Seurat::Embeddings(reference, "pca"))))
  anchors <- Seurat::FindTransferAnchors(
    reference = reference, query = query, normalization.method = "LogNormalize",
    reference.reduction = "pca", reduction = "pcaproject", features = features,
    dims = dims, k.score = as.integer(config$mapping_k), verbose = FALSE
  )
  predictions <- Seurat::TransferData(
    anchorset = anchors, refdata = as.character(reference$reference_label),
    dims = dims, k.weight = min(as.integer(config$mapping_k), ncol(reference) - 1L), verbose = FALSE
  )
  score_columns <- prediction_score_columns(predictions)
  score_matrix <- as.matrix(predictions[, score_columns, drop = FALSE])
  ordered <- t(apply(score_matrix, 1L, sort, decreasing = TRUE))
  margin <- if (ncol(ordered) >= 2L) ordered[, 1] - ordered[, 2] else ordered[, 1]
  second_best <- if (ncol(score_matrix) >= 2L) {
    sub("^prediction.score\\.", "", apply(score_matrix, 1L, function(values) names(sort(values, decreasing = TRUE))[[2]]))
  } else {
    rep(NA_character_, nrow(score_matrix))
  }
  mapped_query <- Seurat::MapQuery(
    anchorset = anchors, query = query, reference = reference,
    refdata = list(reference_label = "reference_label"),
    new.reduction.name = "ref.pca", reference.reduction = "pca", reduction.model = "umap_pca",
    transferdata.args = list(k.weight = min(as.integer(config$mapping_k), ncol(reference) - 1L)),
    verbose = FALSE
  )
  query_embedding <- Seurat::Embeddings(mapped_query, reduction = "ref.pca")
  query_embedding <- query_embedding[rownames(predictions), , drop = FALSE]
  centroids <- reference_label_centroids(reference)
  predicted <- as.character(predictions$predicted.id)
  distance <- vapply(seq_len(nrow(query_embedding)), function(index) {
    label <- predicted[[index]]
    if (!label %in% rownames(centroids)) return(Inf)
    sqrt(sum((query_embedding[index, ] - centroids[label, ])^2))
  }, numeric(1))
  out <- data.frame(
    cell_id = rownames(predictions), predicted_label = predicted,
    prediction_confidence = as.numeric(predictions$prediction.score.max),
    confidence_margin = as.numeric(margin), second_best_label = second_best,
    reference_distance = distance, stringsAsFactors = FALSE
  )
  if (!is.null(thresholds)) out <- classify_mapping_uncertainty(out, thresholds)
  attr(out, "anchors") <- anchors
  attr(out, "mapped_query") <- mapped_query
  out
}


# Archived function: map_region4_sensitivity
map_region4_sensitivity <- function(bundle, frozen_reference, features, thresholds,
                                    marker_config, config) {
  if (!identical(bundle$region_id, "Region_4") ||
      !identical(bundle$downstream_contract, "MAP_TO_REGION_1_3_REFERENCE_WITH_UNCERTAIN")) {
    stop("Region 4 mapping-only bundle is required.", call. = FALSE)
  }
  query <- prepare_seurat_query(bundle, "primary_include", features, config, "REGION4_MAPPING_ONLY_QUERY")
  mapping <- map_query_to_frozen_reference(frozen_reference, query, features, config, thresholds)
  coherence <- assess_mapped_marker_coherence(query, mapping, marker_config, config)
  mapped_fraction <- mean(mapping$mapping_status == "Mapped")
  gate <- data.frame(
    gate_id = "region4_mapping_quality", observed = mapped_fraction,
    threshold = config$mapping_accepted_fraction_min,
    gate_status = if (mapped_fraction >= config$mapping_accepted_fraction_min) "PASS_SENSITIVITY" else "STOP_SENSITIVITY",
    interpretation = "Region 4 never trained or altered the frozen reference",
    stringsAsFactors = FALSE
  )
  list(query = query, mapping = mapping, marker_coherence = coherence, gate = gate)
}


# Archived function: match_cluster_labels
match_cluster_labels <- function(reference_cluster, candidate_cluster) {
  if (is.null(names(reference_cluster)) || is.null(names(candidate_cluster))) {
    stop("Cluster label vectors must be named by cell ID.", call. = FALSE)
  }
  shared <- intersect(names(reference_cluster), names(candidate_cluster))
  if (!length(shared)) stop("No shared cells are available for cluster matching.", call. = FALSE)
  reference <- as.character(reference_cluster[shared])
  candidate <- as.character(candidate_cluster[shared])
  reference_levels <- sort(unique(reference))
  candidate_levels <- sort(unique(candidate))
  score <- matrix(0, nrow = length(reference_levels), ncol = length(candidate_levels),
                  dimnames = list(reference_levels, candidate_levels))
  for (reference_id in reference_levels) {
    reference_cells <- shared[reference == reference_id]
    for (candidate_id in candidate_levels) {
      candidate_cells <- shared[candidate == candidate_id]
      score[reference_id, candidate_id] <- length(intersect(reference_cells, candidate_cells)) /
        length(union(reference_cells, candidate_cells))
    }
  }
  size <- max(nrow(score), ncol(score))
  padded <- matrix(0, nrow = size, ncol = size)
  padded[seq_len(nrow(score)), seq_len(ncol(score))] <- score
  assignment <- clue::solve_LSAP(padded, maximum = TRUE)
  rows <- seq_len(nrow(score))
  columns <- as.integer(assignment[rows])
  keep <- columns <= ncol(score)
  rows <- rows[keep]
  columns <- columns[keep]
  data.frame(
    reference_cluster = rownames(score)[rows],
    candidate_cluster = colnames(score)[columns],
    jaccard = score[cbind(rows, columns)],
    shared_cells = vapply(seq_along(rows), function(index) {
      sum(reference == rownames(score)[rows[[index]]] & candidate == colnames(score)[columns[[index]]])
    }, integer(1)),
    stringsAsFactors = FALSE
  )
}



# Archived function: parse_downstream_cli
parse_downstream_cli <- function(args = commandArgs(trailingOnly = TRUE)) {
  values <- list()
  for (argument in args) {
    if (!startsWith(argument, "--") || !grepl("=", argument, fixed = TRUE)) next
    parts <- strsplit(sub("^--", "", argument), "=", fixed = TRUE)[[1]]
    values[[parts[[1]]]] <- paste(parts[-1], collapse = "=")
  }
  values
}


# Archived function: plot_annotation_overlap_heatmap
plot_annotation_overlap_heatmap <- function(
    object,
    row_var,
    col_var,
    normalize = c("row", "column", "none"),
    min_label = 5,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_order = NULL,
    col_order = NULL,
    legend_title = "Cells (%)",
    row_title = "Current annotation",
    column_title = "Reference annotation",
    show_values = TRUE,
    digits = 0
) {

  requireNamespace("ComplexHeatmap")
  requireNamespace("circlize")
  requireNamespace("grid")

  normalize <- match.arg(normalize)

  ## ----------------------------
  ## 1. Extract metadata
  ## ----------------------------
  meta <- object[[]]

  if (!row_var %in% colnames(meta)) {
    stop("row_var not found in object metadata: ", row_var)
  }

  if (!col_var %in% colnames(meta)) {
    stop("col_var not found in object metadata: ", col_var)
  }

  df <- meta[, c(row_var, col_var), drop = FALSE]

  # Remove NA annotations
  df <- df[
    !is.na(df[[row_var]]) &
      !is.na(df[[col_var]]),
    ,
    drop = FALSE
  ]

  ## ----------------------------
  ## 2. Build contingency table
  ## ----------------------------
  count_mat <- table(
    df[[row_var]],
    df[[col_var]]
  )

  ## ----------------------------
  ## 3. Convert to percentages
  ## ----------------------------
  if (normalize == "row") {

    prop_mat <- prop.table(
      count_mat,
      margin = 1
    ) * 100

  } else if (normalize == "column") {

    prop_mat <- prop.table(
      count_mat,
      margin = 2
    ) * 100

  } else {

    prop_mat <- count_mat / sum(count_mat) * 100
  }

  prop_mat <- as.matrix(prop_mat)

  # Remove empty rows/columns
  prop_mat <- prop_mat[
    rowSums(prop_mat) > 0,
    colSums(prop_mat) > 0,
    drop = FALSE
  ]

  ## ----------------------------
  ## 4. Optional manual ordering
  ## ----------------------------
  if (!is.null(row_order)) {

    row_order <- intersect(
      row_order,
      rownames(prop_mat)
    )

    remaining_rows <- setdiff(
      rownames(prop_mat),
      row_order
    )

    prop_mat <- prop_mat[
      c(row_order, remaining_rows),
      ,
      drop = FALSE
    ]
  }

  if (!is.null(col_order)) {

    col_order <- intersect(
      col_order,
      colnames(prop_mat)
    )

    remaining_cols <- setdiff(
      colnames(prop_mat),
      col_order
    )

    prop_mat <- prop_mat[
      ,
      c(col_order, remaining_cols),
      drop = FALSE
    ]
  }

  ## ----------------------------
  ## 5. Cell-style color scale
  ## ----------------------------
  col_fun <- circlize::colorRamp2(
    c(0, 25, 50, 75, 100),
    c(
      "#FFFFFF",
      "#FDE0DD",
      "#FCAE91",
      "#FB6A4A",
      "#CB181D"
    )
  )

  ## ----------------------------
  ## 6. Heatmap
  ## ----------------------------
  ht <- ComplexHeatmap::Heatmap(
    prop_mat,

    name = "cell_percentage",
    col = col_fun,

    cluster_rows = cluster_rows,
    cluster_columns = cluster_columns,

    row_title = row_title,
    column_title = column_title,

    row_names_side = "left",

    row_names_gp = grid::gpar(
      fontsize = 10
    ),

    column_names_gp = grid::gpar(
      fontsize = 10
    ),

    column_names_rot = 45,

    rect_gp = grid::gpar(
      col = "white",
      lwd = 1
    ),

    border = TRUE,

    cell_fun = if (show_values) {

      function(j, i, x, y, width, height, fill) {

        val <- prop_mat[i, j]

        if (!is.na(val) && val >= min_label) {

          grid::grid.text(
            sprintf(
              paste0("%.", digits, "f"),
              val
            ),
            x,
            y,
            gp = grid::gpar(
              fontsize = 8,
              fontface = ifelse(
                val >= 50,
                "bold",
                "plain"
              ),
              col = ifelse(
                val >= 60,
                "white",
                "black"
              )
            )
          )
        }
      }

    } else {
      NULL
    },

    heatmap_legend_param = list(
      title = legend_title,
      at = c(0, 25, 50, 75, 100),
      labels = c("0", "25", "50", "75", "100"),
      legend_height = grid::unit(3.5, "cm"),
      title_gp = grid::gpar(
        fontsize = 10,
        fontface = "bold"
      ),
      labels_gp = grid::gpar(
        fontsize = 9
      )
    )
  )

  ## ----------------------------
  ## 7. Draw
  ## ----------------------------
  ComplexHeatmap::draw(
    ht,
    heatmap_legend_side = "right"
  )

  invisible(
    list(
      heatmap = ht,
      percentage_matrix = prop_mat,
      count_matrix = count_mat
    )
  )
}


# Archived function: plot_spatial_discrete_overlay
plot_spatial_discrete_overlay <- function(
    object,
    group.by,
    fov = NULL,
    highlight = NULL,
    background_col = "#D9D9D9",
    highlight_cols = NULL,
    base_size = 0.25,
    highlight_size = 0.9,
    base_alpha = 0.7,
    highlight_alpha = 1,
    flip_xy = FALSE,
    dark.background = FALSE,
    axes = FALSE,
    title = NULL,
    subtitle = NULL
) {

  # ============================================================
  # 1. Validate object and metadata
  # ============================================================

  if (!inherits(object, "Seurat")) {
    stop("object must be a Seurat object.", call. = FALSE)
  }

  if (!group.by %in% colnames(object@meta.data)) {
    stop(
      "Metadata variable not found: ",
      group.by,
      call. = FALSE
    )
  }

  if (is.null(fov)) {
    fov <- SeuratObject::DefaultFOV(object)
  }

  if (!fov %in% Seurat::Images(object)) {
    stop(
      "FOV not found: ",
      fov,
      call. = FALSE
    )
  }


  # ============================================================
  # 2. Prepare metadata
  # ============================================================

  meta <- object@meta.data
  meta$cell <- rownames(meta)

  values <- meta[[group.by]]

  if (is.logical(values)) {
    values <- factor(
      values,
      levels = c(TRUE, FALSE)
    )
  } else {
    values <- factor(values)
  }

  meta[[group.by]] <- values
  levels_use <- levels(values)


  # ============================================================
  # 3. Validate highlight groups
  # ============================================================

  if (is.null(highlight)) {
    highlight <- character()
  }

  highlight <- as.character(highlight)

  unknown_highlight <- setdiff(
    highlight,
    levels_use
  )

  if (length(unknown_highlight)) {
    stop(
      "Highlight values not present in ",
      group.by,
      ": ",
      paste(unknown_highlight, collapse = ", "),
      call. = FALSE
    )
  }


  # ============================================================
  # 4. Base ImageDimPlot colors
  #
  # Draw every category grey first. Selected categories are
  # redrawn later as separate layers.
  # ============================================================

  base_cols <- stats::setNames(
    rep(
      background_col,
      length(levels_use)
    ),
    levels_use
  )


  # ============================================================
  # 5. Highlight colors
  # ============================================================

  if (length(highlight)) {

    if (is.null(highlight_cols)) {

      default_cols <- c(
        "#D73027",
        "#0072B2",
        "#009E73",
        "#CC79A7",
        "#E69F00",
        "#56B4E9"
      )

      highlight_cols <- stats::setNames(
        rep(
          default_cols,
          length.out = length(highlight)
        ),
        highlight
      )

    } else if (is.null(names(highlight_cols))) {

      if (length(highlight_cols) != length(highlight)) {
        stop(
          "Unnamed highlight_cols must match highlight length.",
          call. = FALSE
        )
      }

      highlight_cols <- stats::setNames(
        highlight_cols,
        highlight
      )

    } else {

      missing_cols <- setdiff(
        highlight,
        names(highlight_cols)
      )

      if (length(missing_cols)) {
        stop(
          "Missing highlight colors for: ",
          paste(missing_cols, collapse = ", "),
          call. = FALSE
        )
      }

      highlight_cols <- highlight_cols[
        highlight
      ]
    }
  }


  # ============================================================
  # 6. Base ImageDimPlot
  # ============================================================

  p <- Seurat::ImageDimPlot(
    object,
    fov = fov,
    group.by = group.by,
    cols = unname(base_cols),
    size = base_size,
    alpha = base_alpha,
    flip_xy = flip_xy,
    dark.background = dark.background,
    axes = axes
  )


  # ============================================================
  # 7. Get centroid coordinates from the FOV
  # ============================================================

  coords <- SeuratObject::GetTissueCoordinates(
    object[[fov]],
    which = "centroids"
  )

  if (!"cell" %in% names(coords)) {
    coords$cell <- rownames(coords)
  }

  if (!all(c("x", "y") %in% names(coords))) {
    stop(
      "Centroid coordinates do not contain x/y columns.",
      call. = FALSE
    )
  }


  # ============================================================
  # 8. Align metadata without merge reordering
  # ============================================================

  coords[[group.by]] <- meta[
    match(coords$cell, meta$cell),
    group.by
  ]


  # ============================================================
  # 9. Convert tissue coordinates to ImageDimPlot plot coordinates
  #
  # Important:
  # Seurat's SingleImagePlot internally maps:
  #
  #       plot x <- tissue y
  #       plot y <- tissue x
  #
  # and then flip_xy controls whether coord_flip() is applied.
  #
  # Therefore:
  #
  # flip_xy = TRUE:
  #     ggplot layer coordinates = (y, x)
  #
  # flip_xy = FALSE:
  #     coord_flip() effectively displays them as (x, y)
  #
  # For added ggplot layers we need to supply coordinates in the
  # coordinate system expected BEFORE coord_flip().
  # ============================================================

  if (isTRUE(flip_xy)) {

    coords$plot_x <- coords$y
    coords$plot_y <- coords$x

  } else {

    # ImageDimPlot will apply coord_flip(), so the custom layer
    # must also enter as x=y, y=x.
    coords$plot_x <- coords$y
    coords$plot_y <- coords$x
  }


  # ============================================================
  # 10. Overlay highlighted groups
  # ============================================================

  if (length(highlight)) {

    for (group_value in highlight) {

      tmp <- coords[
        !is.na(coords[[group.by]]) &
          as.character(coords[[group.by]]) == group_value,
        ,
        drop = FALSE
      ]

      if (!nrow(tmp)) {
        next
      }

      p <- p +
        ggplot2::geom_point(
          data = tmp,
          ggplot2::aes(
            x = plot_x,
            y = plot_y
          ),
          inherit.aes = FALSE,
          colour = unname(
            highlight_cols[group_value]
          ),
          size = highlight_size,
          alpha = highlight_alpha
        )
    }
  }


  # ============================================================
  # 11. Labels
  # ============================================================

  p <- p +
    ggplot2::labs(
      title = title,
      subtitle = subtitle
    )


  # ============================================================
  # 12. Return
  # ============================================================

  p
}


# Archived function: prediction_score_columns
prediction_score_columns <- function(predictions) {
  setdiff(grep("^prediction.score\\.", names(predictions), value = TRUE), "prediction.score.max")
}


# Archived function: prepare_seurat_query
prepare_seurat_query <- function(bundle, mask_name, genes, config, role) {
  assert_downstream_model_environment(config)
  cells <- bundle$cell_metadata[bundle$cell_metadata[[mask_name]], , drop = FALSE]
  ids <- as.character(cells$cell_id)
  metadata <- cells[match(ids, cells$cell_id), , drop = FALSE]
  rownames(metadata) <- ids
  object <- Seurat::CreateSeuratObject(
    counts = bundle$counts[, ids, drop = FALSE], meta.data = metadata,
    project = role, min.cells = 0L, min.features = 0L
  )
  object <- Seurat::NormalizeData(
    object, normalization.method = "LogNormalize",
    scale.factor = config$normalization_scale_factor, verbose = FALSE
  )
  object <- Seurat::ScaleData(object, features = genes, verbose = FALSE)
  object
}


# Archived function: primary_model_feature_policy
primary_model_feature_policy <- function() "FIXED_245_NO_VARIABLE_FEATURE_SELECTION"


# Archived function: read_canonical_marker_config
read_canonical_marker_config <- function(path, panel_genes, gene_decision) {
  if (!file.exists(path)) stop(sprintf("Missing canonical marker config: %s", path), call. = FALSE)
  markers <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  required <- c("cell_type", "gene", "direction", "marker_group", "use_policy")
  if (length(setdiff(required, names(markers)))) {
    stop("Canonical marker config lacks required columns.", call. = FALSE)
  }
  if (!length(panel_genes) || anyDuplicated(panel_genes)) stop("Panel genes must be non-empty and unique.", call. = FALSE)
  gene_required <- c("gene", "technical_risk_status")
  if (length(setdiff(gene_required, names(gene_decision))) || anyDuplicated(gene_decision$gene)) {
    stop("Gene decision must contain unique gene and technical_risk_status columns.", call. = FALSE)
  }
  markers <- markers[markers$gene %in% panel_genes, required, drop = FALSE]
  if (!nrow(markers)) stop("No configured canonical markers are present in the panel.", call. = FALSE)
  risk_genes <- gene_decision$gene[
    gene_decision$technical_risk_status == "TECHNICAL_RISK_SENSITIVITY_ONLY"
  ]
  markers$use_policy <- ifelse(
    markers$gene %in% risk_genes,
    "VALIDATION_ONLY_TECHNICAL_RISK",
    "PRIMARY_SUPPORT"
  )
  marker_counts <- table(markers$cell_type)
  markers$cell_type_panel_marker_count <- as.integer(marker_counts[markers$cell_type])
  markers$cell_type_support_status <- ifelse(
    markers$cell_type_panel_marker_count >= 2L,
    "SUPPORTED",
    "INSUFFICIENT_PANEL_SUPPORT"
  )
  rownames(markers) <- NULL
  markers
}


# Archived function: read_downstream_handoff
read_downstream_handoff <- function(qc_run_root, release_path = NULL) {
  if (!dir.exists(qc_run_root)) stop(sprintf("Missing Phase 0-2 QC run root: %s", qc_run_root), call. = FALSE)
  region_ids <- paste0("Region_", seq_len(4L))
  bundle_paths <- file.path(qc_run_root, "downstream_inputs", paste0(region_ids, ".downstream_input.rds"))
  slide_root <- file.path(qc_run_root, "slide_summary")
  table_paths <- c(
    gene_decision = file.path(slide_root, "gene_downstream_decision.tsv"),
    eos_decision = file.path(slide_root, "eos_gene_decision_summary.tsv"),
    section_decision = file.path(slide_root, "section_downstream_decision.tsv"),
    masks = file.path(slide_root, "cell_downstream_masks.tsv.gz"),
    release = release_path %||% file.path(slide_root, "evidence_only_qc_release.tsv")
  )
  missing <- c(bundle_paths[!file.exists(bundle_paths)], table_paths[!file.exists(table_paths)])
  if (length(missing)) stop(sprintf("Missing downstream handoff artifacts: %s", paste(missing, collapse = ", ")), call. = FALSE)
  regions <- setNames(lapply(bundle_paths, readRDS), region_ids)
  list(
    qc_run_root = normalizePath(qc_run_root, winslash = "/", mustWork = TRUE),
    release_path = normalizePath(table_paths[["release"]], winslash = "/", mustWork = TRUE),
    regions = regions,
    gene_decision = utils::read.delim(table_paths[["gene_decision"]], check.names = FALSE),
    eos_decision = utils::read.delim(table_paths[["eos_decision"]], check.names = FALSE),
    section_decision = utils::read.delim(table_paths[["section_decision"]], check.names = FALSE),
    masks = utils::read.delim(gzfile(table_paths[["masks"]]), check.names = FALSE),
    release = utils::read.delim(table_paths[["release"]], check.names = FALSE)
  )
}


# Archived function: read_downstream_reference_config
read_downstream_reference_config <- function(path) {
  if (!file.exists(path)) stop(sprintf("Missing downstream reference config: %s", path), call. = FALSE)
  table <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  required_columns <- c("key", "value", "description")
  if (length(setdiff(required_columns, names(table)))) {
    stop("Downstream reference config must contain key, value, and description columns.", call. = FALSE)
  }
  if (!nrow(table) || any(!nzchar(table$key)) || anyDuplicated(table$key)) {
    stop("Downstream reference config keys must be non-empty and unique.", call. = FALSE)
  }
  required_keys <- c(
    "seed", "normalization_scale_factor", "primary_gene_count", "conservative_gene_count",
    "technical_risk_gene_count", "raw_gene_count", "n_pcs", "knn_k", "leiden_resolution", "marker_score_margin",
    "seurat_min_version", "harmony_min_version", "harmony_theta", "harmony_lambda",
    "harmony_sigma", "harmony_max_iter", "harmony_reference_region", "mapping_folds",
    "mapping_k", "mapping_confidence_quantile", "mapping_distance_quantile",
    "mapping_min_label_cells", "major_label_fraction", "major_label_cells",
    "section_cluster_dominance", "section_cluster_min_cells", "label_stability_min",
    "gene_sensitivity_concordance_min", "section_predictability_permutations",
    "section_predictability_margin", "eos_assignment_jaccard_min", "eos_score_spearman_min",
    "eos_module_nbin", "eos_module_ctrl", "marker_detection_min",
    "marker_average_log_expression_min", "exclusion_contradiction_max",
    "mapping_accepted_fraction_min", "mapping_unrepresented_fraction_review",
    "mapping_marker_coherence_min", "mapping_label_stability_min",
    "consensus_label_stability_min", "eos_identity_min_markers",
    "eos_identity_detection_min", "eos_state_min_genes_detected"
  )
  missing_keys <- setdiff(required_keys, table$key)
  if (length(missing_keys)) {
    stop(sprintf("Downstream reference config is missing required keys: %s", paste(missing_keys, collapse = ", ")), call. = FALSE)
  }
  integer_keys <- c(
    "seed", "primary_gene_count", "conservative_gene_count", "technical_risk_gene_count",
    "raw_gene_count", "n_pcs", "knn_k", "harmony_max_iter", "mapping_folds",
    "mapping_k", "mapping_min_label_cells", "major_label_cells",
    "section_cluster_min_cells", "section_predictability_permutations",
    "eos_module_nbin", "eos_module_ctrl", "eos_identity_min_markers",
    "eos_state_min_genes_detected"
  )
  numeric_keys <- c(
    "normalization_scale_factor", "leiden_resolution", "marker_score_margin", "harmony_theta", "harmony_lambda",
    "harmony_sigma", "mapping_confidence_quantile", "mapping_distance_quantile",
    "major_label_fraction", "section_cluster_dominance", "label_stability_min",
    "gene_sensitivity_concordance_min", "section_predictability_margin",
    "eos_assignment_jaccard_min", "eos_score_spearman_min", "marker_score_margin",
    "marker_detection_min", "marker_average_log_expression_min",
    "exclusion_contradiction_max", "mapping_accepted_fraction_min",
    "mapping_unrepresented_fraction_review", "mapping_marker_coherence_min",
    "mapping_label_stability_min", "consensus_label_stability_min",
    "eos_identity_detection_min"
  )
  values <- setNames(as.list(table$value), table$key)
  for (key in intersect(integer_keys, names(values))) values[[key]] <- as.integer(values[[key]])
  for (key in intersect(numeric_keys, names(values))) values[[key]] <- as.numeric(values[[key]])
  invalid_integer <- integer_keys[!vapply(values[integer_keys], function(value) {
    length(value) == 1L && !is.na(value) && is.finite(value) && value >= 1L
  }, logical(1))]
  if (length(invalid_integer)) {
    stop(sprintf("Downstream reference config requires a positive integer for: %s", paste(invalid_integer, collapse = ", ")), call. = FALSE)
  }
  invalid_numeric <- numeric_keys[!vapply(values[numeric_keys], function(value) {
    length(value) == 1L && !is.na(value) && is.finite(value)
  }, logical(1))]
  if (length(invalid_numeric)) {
    stop(sprintf("Downstream reference config requires a finite numeric value for: %s", paste(invalid_numeric, collapse = ", ")), call. = FALSE)
  }
  open_unit_keys <- c("mapping_confidence_quantile", "mapping_distance_quantile")
  invalid_open_unit <- open_unit_keys[!vapply(values[open_unit_keys], function(value) value > 0 && value < 1, logical(1))]
  if (length(invalid_open_unit)) {
    stop(sprintf("Downstream reference config values must be between 0 and 1 (exclusive) for: %s", paste(invalid_open_unit, collapse = ", ")), call. = FALSE)
  }
  closed_unit_keys <- c(
    "major_label_fraction", "section_cluster_dominance", "label_stability_min",
    "gene_sensitivity_concordance_min", "section_predictability_margin",
    "eos_assignment_jaccard_min", "eos_score_spearman_min",
    "marker_detection_min", "exclusion_contradiction_max",
    "mapping_accepted_fraction_min", "mapping_unrepresented_fraction_review",
    "mapping_marker_coherence_min", "mapping_label_stability_min",
    "consensus_label_stability_min", "eos_identity_detection_min"
  )
  invalid_closed_unit <- closed_unit_keys[!vapply(values[closed_unit_keys], function(value) value >= 0 && value <= 1, logical(1))]
  if (length(invalid_closed_unit)) {
    stop(sprintf("Downstream reference config values must be between 0 and 1 for: %s", paste(invalid_closed_unit, collapse = ", ")), call. = FALSE)
  }
  values
}


# Archived function: read_region3_anchor_artifacts
read_region3_anchor_artifacts <- function(output_root) {
  definitions <- utils::read.delim(file.path(output_root, "region3_anchor_branches.tsv"), check.names = FALSE)
  validation <- readRDS(file.path(output_root, "region3_anchor_validation.rds"))
  branches <- setNames(lapply(definitions$branch_id, function(branch) {
    root <- file.path(output_root, branch)
    object <- readRDS(file.path(root, "reference_object.rds"))
    list(
      object = object,
      cell_labels = utils::read.delim(file.path(root, "cell_labels.tsv"), check.names = FALSE),
      marker_support = utils::read.delim(file.path(root, "marker_support.tsv"), check.names = FALSE),
      marker_scores = utils::read.delim(file.path(root, "marker_scores.tsv"), check.names = FALSE),
      feature_names = object@misc$downstream_reference_contract$feature_names,
      mask_name = definitions$mask_name[definitions$branch_id == branch],
      gene_set_name = definitions$gene_set_name[definitions$branch_id == branch],
      scientific_role = definitions$scientific_role[definitions$branch_id == branch]
    )
  }), definitions$branch_id)
  list(branch_definitions = definitions, branches = branches, validation = validation)
}


# Archived function: reference_label_centroids
reference_label_centroids <- function(reference, label_column = "reference_label") {
  embedding <- Seurat::Embeddings(reference, reduction = "pca")
  labels <- as.character(reference[[label_column, drop = TRUE]])
  split_rows <- split(seq_len(nrow(embedding)), labels)
  do.call(rbind, lapply(names(split_rows), function(label) {
    values <- matrix(colMeans(embedding[split_rows[[label]], , drop = FALSE]), nrow = 1L)
    rownames(values) <- label
    values
  }))
}


# Archived function: region3_anchor_branch_definitions
region3_anchor_branch_definitions <- function() {
  data.frame(
    branch_id = c("primary_245", "strict_245", "hotspot_245", "conservative_67", "complete_479"),
    mask_name = c(
      "primary_include", "strict_include", "hotspot_sensitivity_include",
      "primary_include", "primary_include"
    ),
    gene_set_name = c(
      "provisional_primary_features", "provisional_primary_features",
      "provisional_primary_features", "conservative_no_signal_detected",
      "raw_complete_panel"
    ),
    scientific_role = c(
      "PRIMARY", "CELL_QC_SENSITIVITY", "SPATIAL_HOTSPOT_SENSITIVITY",
      "CONSERVATIVE_LOW_INFORMATION_STRESS_TEST", "CLEAN_SECTION_COMPLETE_PANEL_SENSITIVITY"
    ),
    stringsAsFactors = FALSE
  )
}


# Archived function: require_downstream_argument
require_downstream_argument <- function(arguments, name) {
  value <- arguments[[name]] %||% ""
  if (length(value) != 1L || is.na(value) || !nzchar(value)) {
    stop(sprintf("Missing required --%s= argument.", name), call. = FALSE)
  }
  value
}


# Archived function: run_eosinophil_robustness
run_eosinophil_robustness <- function(reference, eos_decision, marker_config, config) {
  identity <- identify_eosinophils_independently(reference, marker_config, config)
  eos_ids <- identity$cell_id[identity$validated_eosinophil]
  if (length(eos_ids) < 10L) {
    return(list(
      identity = identity, scores = data.frame(), stability = data.frame(
        gate_id = "eos_cell_count", observed = length(eos_ids), threshold = 10L,
        gate_status = "STOP_EOS", interpretation = "Too few independently validated Eosinophils", stringsAsFactors = FALSE
      ), gene_coherence = data.frame(), leave_one_gene_out = data.frame(),
      mask_stability = data.frame(), sensitivity_method_status = data.frame(),
      complexity_diagnostics = data.frame()
    ))
  }
  data <- get_seurat_normalized_data(reference, "RNA")[, eos_ids, drop = FALSE]
  set_for <- function(set, retained = NULL) {
    rows <- eos_decision$gene_set == set
    if (!is.null(retained)) rows <- rows & eos_decision[[retained]]
    as.character(eos_decision$gene[rows])
  }
  short53 <- set_for("short_lived", "retained_provisional")
  long53 <- set_for("long_lived", "retained_provisional")
  short100 <- set_for("short_lived")
  long100 <- set_for("long_lived")
  short_nonrib <- short100[!grepl("^Rp[sl]", short100)]
  scores <- data.frame(
    cell_id = eos_ids,
    short_primary_53 = score_standardized_gene_set(data, short53, config$eos_state_min_genes_detected),
    long_primary_53 = score_standardized_gene_set(data, long53, config$eos_state_min_genes_detected),
    short_complete_100 = score_standardized_gene_set(data, short100, config$eos_state_min_genes_detected),
    long_complete_100 = score_standardized_gene_set(data, long100, config$eos_state_min_genes_detected),
    short_nonribosomal = score_standardized_gene_set(data, short_nonrib, config$eos_state_min_genes_detected),
    stringsAsFactors = FALSE
  )
  scores$state_primary <- ifelse(
    is.na(scores$short_primary_53) | is.na(scores$long_primary_53), "Unresolved",
    ifelse(scores$short_primary_53 > scores$long_primary_53, "AT_short_like", "AT_long_like")
  )
  eos_metadata <- reference@meta.data[match(eos_ids, rownames(reference@meta.data)), , drop = FALSE]
  for (column in intersect(c("region_id", "mouse_id", "strict_include", "hotspot_sensitivity_include", "nCount_RNA", "nFeature_RNA"), names(eos_metadata))) {
    scores[[column]] <- eos_metadata[[column]]
  }
  gene_coherence_for <- function(genes, composite, state_name) {
    genes <- intersect(genes, rownames(data))
    if (!length(genes)) return(data.frame())
    data.frame(
      state = state_name, gene = genes,
      detected_fraction = Matrix::rowMeans(data[genes, , drop = FALSE] > 0),
      spearman_with_composite = vapply(genes, function(gene) {
        suppressWarnings(stats::cor(as.numeric(data[gene, ]), composite, method = "spearman", use = "pairwise.complete.obs"))
      }, numeric(1)), stringsAsFactors = FALSE
    )
  }
  gene_coherence <- rbind(
    gene_coherence_for(short53, scores$short_primary_53, "short_primary_53"),
    gene_coherence_for(long53, scores$long_primary_53, "long_primary_53")
  )
  leave_one_out_for <- function(genes, full_score, state_name) {
    genes <- intersect(genes, rownames(data))
    if (length(genes) < 2L) return(data.frame())
    do.call(rbind, lapply(genes, function(omitted) {
      loo <- score_standardized_gene_set(data, setdiff(genes, omitted), max(1L, config$eos_state_min_genes_detected - 1L))
      data.frame(
        state = state_name, omitted_gene = omitted,
        spearman_with_full = suppressWarnings(stats::cor(loo, full_score, method = "spearman", use = "pairwise.complete.obs")),
        stringsAsFactors = FALSE
      )
    }))
  }
  leave_one_gene_out <- rbind(
    leave_one_out_for(short53, scores$short_primary_53, "short_primary_53"),
    leave_one_out_for(long53, scores$long_primary_53, "long_primary_53")
  )
  identity_set <- identity$cell_id[identity$validated_eosinophil]
  jaccard_subset <- function(mask_name) {
    if (!mask_name %in% names(reference@meta.data)) return(NA_real_)
    kept <- rownames(reference@meta.data)[as.logical(reference@meta.data[[mask_name]])]
    length(intersect(identity_set, kept)) / length(union(identity_set, intersect(identity_set, kept)))
  }
  mask_stability <- data.frame(
    comparison = c("primary_vs_strict_eos_assignment", "primary_vs_hotspot_sensitivity_eos_assignment"),
    jaccard = c(jaccard_subset("strict_include"), jaccard_subset("hotspot_sensitivity_include")),
    threshold = config$eos_assignment_jaccard_min, stringsAsFactors = FALSE
  )
  mask_stability$gate_status <- ifelse(
    is.na(mask_stability$jaccard), "NOT_ESTIMABLE",
    ifelse(mask_stability$jaccard >= mask_stability$threshold, "PASS_EOS", "STOP_EOS")
  )
  sensitivity_method_status <- data.frame(
    method = c("GENEWISE_STANDARDIZED_MEAN", "ADDMODULESCORE_TARGETED_PANEL", "UCELL_TARGETED_PANEL"),
    scientific_role = c("PRIMARY", "SENSITIVITY_ONLY", "SENSITIVITY_ONLY"),
    status = c("COMPUTED", "NOT_COMPUTED", "NOT_COMPUTED"), stringsAsFactors = FALSE
  )
  module_object <- reference[, eos_ids]
  module_result <- tryCatch({
    x <- Seurat::AddModuleScore(
      module_object, features = list(intersect(short53, rownames(module_object))),
      nbin = as.integer(config$eos_module_nbin), ctrl = as.integer(config$eos_module_ctrl),
      name = "short_targeted_module", seed = as.integer(config$seed)
    )
    x <- Seurat::AddModuleScore(
      x, features = list(intersect(long53, rownames(x))),
      nbin = as.integer(config$eos_module_nbin), ctrl = as.integer(config$eos_module_ctrl),
      name = "long_targeted_module", seed = as.integer(config$seed)
    )
    x
  }, error = identity)
  if (!inherits(module_result, "error")) {
    module_object <- module_result
    scores$short_addmodule_sensitivity <- module_object$short_targeted_module1[match(scores$cell_id, colnames(module_object))]
    scores$long_addmodule_sensitivity <- module_object$long_targeted_module1[match(scores$cell_id, colnames(module_object))]
    sensitivity_method_status$status[sensitivity_method_status$method == "ADDMODULESCORE_TARGETED_PANEL"] <- "COMPUTED_SENSITIVITY_ONLY"
  } else {
    scores$short_addmodule_sensitivity <- NA_real_
    scores$long_addmodule_sensitivity <- NA_real_
    sensitivity_method_status$status[sensitivity_method_status$method == "ADDMODULESCORE_TARGETED_PANEL"] <- "FAILED_NONBLOCKING_TARGETED_CONTROLS"
  }
  if (requireNamespace("UCell", quietly = TRUE)) {
    ucell_result <- tryCatch(
      UCell::AddModuleScore_UCell(
        reference[, eos_ids], features = list(short_ucell = short53, long_ucell = long53),
        assay = "RNA", name = NULL
      ), error = identity
    )
    if (!inherits(ucell_result, "error")) {
      short_column <- grep("^short_ucell.*UCell$", names(ucell_result@meta.data), value = TRUE)[1]
      long_column <- grep("^long_ucell.*UCell$", names(ucell_result@meta.data), value = TRUE)[1]
      if (!is.na(short_column) && !is.na(long_column)) {
        scores$short_ucell_sensitivity <- ucell_result@meta.data[[short_column]][match(scores$cell_id, rownames(ucell_result@meta.data))]
        scores$long_ucell_sensitivity <- ucell_result@meta.data[[long_column]][match(scores$cell_id, rownames(ucell_result@meta.data))]
        sensitivity_method_status$status[sensitivity_method_status$method == "UCELL_TARGETED_PANEL"] <- "COMPUTED_SENSITIVITY_ONLY"
      } else {
        scores$short_ucell_sensitivity <- NA_real_
        scores$long_ucell_sensitivity <- NA_real_
        sensitivity_method_status$status[sensitivity_method_status$method == "UCELL_TARGETED_PANEL"] <- "FAILED_NONBLOCKING_OUTPUT_SCHEMA"
      }
    } else {
      scores$short_ucell_sensitivity <- NA_real_
      scores$long_ucell_sensitivity <- NA_real_
      sensitivity_method_status$status[sensitivity_method_status$method == "UCELL_TARGETED_PANEL"] <- "FAILED_NONBLOCKING"
    }
  } else {
    scores$short_ucell_sensitivity <- NA_real_
    scores$long_ucell_sensitivity <- NA_real_
    sensitivity_method_status$status[sensitivity_method_status$method == "UCELL_TARGETED_PANEL"] <- "PACKAGE_UNAVAILABLE_NONBLOCKING"
  }
  complexity_diagnostics <- data.frame()
  if (all(c("nCount_RNA", "nFeature_RNA") %in% names(scores))) {
    complexity_diagnostics <- data.frame(
      score = c("short_primary_53", "long_primary_53"),
      spearman_nCount = c(
        suppressWarnings(stats::cor(scores$short_primary_53, scores$nCount_RNA, method = "spearman", use = "pairwise.complete.obs")),
        suppressWarnings(stats::cor(scores$long_primary_53, scores$nCount_RNA, method = "spearman", use = "pairwise.complete.obs"))
      ),
      spearman_nFeature = c(
        suppressWarnings(stats::cor(scores$short_primary_53, scores$nFeature_RNA, method = "spearman", use = "pairwise.complete.obs")),
        suppressWarnings(stats::cor(scores$long_primary_53, scores$nFeature_RNA, method = "spearman", use = "pairwise.complete.obs"))
      ),
      interpretation = "Diagnostic only; region and mouse are not regressed from biological scores",
      stringsAsFactors = FALSE
    )
  }
  correlation <- function(x, y) suppressWarnings(stats::cor(x, y, method = "spearman", use = "pairwise.complete.obs"))
  stability <- data.frame(
    gate_id = c("short53_vs_short100", "long53_vs_long100", "short100_vs_nonribosomal"),
    observed = c(
      correlation(scores$short_primary_53, scores$short_complete_100),
      correlation(scores$long_primary_53, scores$long_complete_100),
      correlation(scores$short_complete_100, scores$short_nonribosomal)
    ),
    threshold = config$eos_score_spearman_min, stringsAsFactors = FALSE
  )
  stability$gate_status <- ifelse(
    !is.na(stability$observed) & stability$observed >= stability$threshold, "PASS_EOS", "STOP_EOS"
  )
  stability$interpretation <- "Gene-wise standardized mean is primary; targeted-panel module scores are sensitivity only"
  stability <- rbind(
    stability,
    data.frame(
      gate_id = paste0("mask_", mask_stability$comparison),
      observed = mask_stability$jaccard, threshold = mask_stability$threshold,
      gate_status = mask_stability$gate_status,
      interpretation = "Eosinophil identity must remain stable across QC masks",
      stringsAsFactors = FALSE
    )
  )
  list(
    identity = identity, scores = scores, stability = stability,
    gene_coherence = gene_coherence, leave_one_gene_out = leave_one_gene_out,
    mask_stability = mask_stability, sensitivity_method_status = sensitivity_method_status,
    complexity_diagnostics = complexity_diagnostics
  )
}


# Archived function: score_standardized_gene_set
score_standardized_gene_set <- function(data, genes, minimum_detected) {
  genes <- intersect(genes, rownames(data))
  if (!length(genes)) return(rep(NA_real_, ncol(data)))
  values <- as.matrix(data[genes, , drop = FALSE])
  standardized <- t(scale(t(values)))
  standardized[!is.finite(standardized)] <- 0
  detected <- colSums(values > 0)
  score <- colMeans(standardized)
  score[detected < as.integer(minimum_detected)] <- NA_real_
  score
}


# Archived function: validate_anchor_branches
validate_anchor_branches <- function(primary, strict, hotspot, conservative, complete,
                                     markers, config) {
  branches <- list(primary = primary, strict = strict, hotspot = hotspot,
                   conservative = conservative, complete = complete)
  for (name in names(branches)) {
    labels <- branches[[name]]$cell_labels
    if (!is.data.frame(labels) || anyDuplicated(labels$cell_id) ||
        length(setdiff(c("cell_id", "reference_label"), names(labels)))) {
      stop(sprintf("Anchor branch %s has an invalid label table.", name), call. = FALSE)
    }
  }
  primary_labels <- primary$cell_labels
  counts <- table(primary_labels$reference_label)
  major <- names(counts)[
    counts >= as.integer(config$major_label_cells) |
      counts / nrow(primary_labels) >= as.numeric(config$major_label_fraction)
  ]
  primary_support <- primary$marker_support
  primary_supported <- major[major %in% primary_support$reference_label[
    primary_support$support_status == "SUPPORTED" & !primary_support$risk_only_support
  ]]
  marker_ok <- length(primary_supported) == length(major) && length(major) > 0L
  strict_agreement <- label_agreement_on_major(primary_labels, strict$cell_labels, major)
  hotspot_agreement <- label_agreement_on_major(primary_labels, hotspot$cell_labels, major)
  complete_agreement <- label_agreement_on_major(primary_labels, complete$cell_labels, major)
  conservative_supported <- unique(conservative$marker_support$reference_label[
    conservative$marker_support$support_status == "SUPPORTED"
  ])
  conservative_estimable <- length(major) > 0L && all(major %in% conservative_supported)
  conservative_agreement <- if (conservative_estimable) {
    label_agreement_on_major(primary_labels, conservative$cell_labels, major)
  } else {
    NA_real_
  }
  gates <- data.frame(
    gate_id = c(
      "major_label_marker_support", "primary_strict_stability", "hotspot_stability",
      "complete479_vs_primary245_stability", "gene245_vs_gene67_stability"
    ),
    observed = c(if (marker_ok) 1 else 0, strict_agreement, hotspot_agreement,
                 complete_agreement, conservative_agreement),
    threshold = c(1, config$label_stability_min, config$label_stability_min,
                  config$label_stability_min, config$gene_sensitivity_concordance_min),
    gate_status = NA_character_, stringsAsFactors = FALSE
  )
  gates$gate_status[1:4] <- ifelse(
    !is.na(gates$observed[1:4]) & gates$observed[1:4] >= gates$threshold[1:4], "PASS", "STOP"
  )
  gates$gate_status[5] <- if (!conservative_estimable) {
    "NOT_ESTIMABLE_GENE67"
  } else if (conservative_agreement >= config$gene_sensitivity_concordance_min) {
    "PASS"
  } else {
    "REVIEW_GENE67"
  }
  decision <- if (any(gates$gate_status == "STOP")) {
    "STOP_ANCHOR"
  } else if (any(grepl("REVIEW|NOT_ESTIMABLE", gates$gate_status))) {
    "REVIEW_ANCHOR"
  } else {
    "PASS_ANCHOR"
  }
  list(
    decision = decision, gate_table = gates, major_labels = major,
    major_label_counts = as.data.frame(counts, stringsAsFactors = FALSE),
    scientific_interpretation = "Region_3 is a technical anchor, not biological ground truth"
  )
}


# Archived function: validate_downstream_handoff
validate_downstream_handoff <- function(qc_run_root, release_path = NULL, stop_on_error = TRUE) {
  handoff <- read_downstream_handoff(qc_run_root, release_path)
  region_ids <- paste0("Region_", seq_len(4L))
  expected_status <- c("PRIMARY_CONDITIONAL", "PRIMARY_CONDITIONAL", "PRIMARY", "SENSITIVITY_ONLY")
  expected_contract <- c(rep("REFERENCE_ELIGIBILITY_FROM_CELL_MASKS", 3L), "MAP_TO_REGION_1_3_REFERENCE_WITH_UNCERTAIN")
  checks <- list()
  add_check <- function(check, passed, details) {
    checks[[length(checks) + 1L]] <<- data.frame(
      check = check, status = if (isTRUE(passed)) "PASS" else "FAIL",
      details = as.character(details), stringsAsFactors = FALSE
    )
  }

  add_check("four_region_bundles", identical(names(handoff$regions), region_ids), paste(names(handoff$regions), collapse = ","))
  bundle_regions <- vapply(handoff$regions, function(bundle) as.character(bundle$region_id %||% ""), character(1))
  bundle_status <- vapply(handoff$regions, function(bundle) as.character(bundle$section_status %||% ""), character(1))
  bundle_contract <- vapply(handoff$regions, function(bundle) as.character(bundle$downstream_contract %||% ""), character(1))
  add_check("bundle_region_identity", identical(unname(bundle_regions), region_ids), paste(bundle_regions, collapse = ","))
  add_check("section_status_contract", identical(unname(bundle_status), expected_status), paste(bundle_status, collapse = ","))
  add_check("region4_mapping_only_contract", identical(unname(bundle_contract), expected_contract), paste(bundle_contract, collapse = ","))

  sparse_ok <- vapply(handoff$regions, function(bundle) inherits(bundle$counts, "sparseMatrix"), logical(1))
  raw_ok <- vapply(handoff$regions, function(bundle) isTRUE(bundle$raw_counts_preserved), logical(1))
  add_check("sparse_raw_counts", all(sparse_ok & raw_ok), sprintf("sparse=%d/4; preserved=%d/4", sum(sparse_ok), sum(raw_ok)))

  alignment_ok <- vapply(handoff$regions, function(bundle) {
    counts <- bundle$counts
    cells <- bundle$cell_metadata
    genes <- bundle$gene_sets$raw_complete_panel
    !is.null(counts) && is.data.frame(cells) && length(genes) == nrow(counts) &&
      identical(colnames(counts), as.character(cells$cell_id)) &&
      setequal(rownames(counts), as.character(genes)) &&
      !anyDuplicated(cells$cell_id) && !anyDuplicated(rownames(counts)) && !anyDuplicated(genes)
  }, logical(1))
  add_check("matrix_cell_gene_alignment", all(alignment_ok), paste(names(alignment_ok)[!alignment_ok], collapse = ","))

  all_cells <- do.call(rbind, lapply(handoff$regions, function(bundle) bundle$cell_metadata[, c("region_id", "cell_id"), drop = FALSE]))
  add_check("global_cell_identity", !anyDuplicated(paste(all_cells$region_id, all_cells$cell_id, sep = "|")), sprintf("cells=%d", nrow(all_cells)))

  required_mask_columns <- c("region_id", "cell_id", "primary_include", "strict_include", "hotspot_sensitivity_include")
  bundle_cell_keys <- sort(paste(all_cells$region_id, all_cells$cell_id, sep = "|"))
  mask_cell_keys <- if (all(c("region_id", "cell_id") %in% names(handoff$masks))) {
    sort(paste(handoff$masks$region_id, handoff$masks$cell_id, sep = "|"))
  } else {
    character()
  }
  masks_ok <- !length(setdiff(required_mask_columns, names(handoff$masks))) &&
    !anyDuplicated(paste(handoff$masks$region_id, handoff$masks$cell_id, sep = "|")) &&
    nrow(handoff$masks) == nrow(all_cells) &&
    identical(mask_cell_keys, bundle_cell_keys) &&
    all(!handoff$masks$strict_include | handoff$masks$primary_include) &&
    all(!handoff$masks$hotspot_sensitivity_include | handoff$masks$primary_include)
  add_check("cell_mask_contract", masks_ok, sprintf("mask_rows=%d; bundle_cells=%d", nrow(handoff$masks), nrow(all_cells)))

  gene <- handoff$gene_decision
  gene_ok <- all(c("gene", "primary_feature_status", "conservative_evidence_status", "technical_risk_status") %in% names(gene)) &&
    nrow(gene) == 479L && !anyDuplicated(gene$gene) &&
    sum(gene$primary_feature_status == "PROVISIONAL_PRIMARY_FEATURES") == 245L &&
    sum(gene$conservative_evidence_status == "CONSERVATIVE_NO_SIGNAL_DETECTED") == 67L &&
    sum(gene$technical_risk_status == "TECHNICAL_RISK_SENSITIVITY_ONLY") == 234L
  add_check("gene_tier_contract", gene_ok, sprintf("rows=%d", nrow(gene)))

  eos <- handoff$eos_decision
  eos_ok <- all(c("gene", "retained_provisional") %in% names(eos)) &&
    sum(as.logical(eos$retained_provisional), na.rm = TRUE) == 53L
  add_check("eosinophil_gene_contract", eos_ok, sprintf("retained=%d", sum(as.logical(eos$retained_provisional), na.rm = TRUE)))

  sections <- handoff$section_decision
  section_ok <- all(c("region_id", "section_status") %in% names(sections)) &&
    identical(as.character(sections$region_id), region_ids) &&
    identical(as.character(sections$section_status), expected_status)
  add_check("slide_section_decisions", section_ok, paste(sections$section_status, collapse = ","))

  required_qc_gates <- c(
    "fixed_section_decisions", "cell_mask_reconciliation", "gene_tier_reconciliation",
    "eos_gene_reconciliation", "region4_excluded_from_reference_definition"
  )
  release <- handoff$release
  gate_ok <- all(c("gate_id", "gate_status") %in% names(release)) &&
    all(required_qc_gates %in% release$gate_id) &&
    all(release$gate_status[match(required_qc_gates, release$gate_id)] == "PASS")
  add_check("completed_qc_release_gates", gate_ok, paste(release$gate_status[match(required_qc_gates, release$gate_id)], collapse = ","))

  result <- do.call(rbind, checks)
  if (isTRUE(stop_on_error) && any(result$status == "FAIL")) {
    failed <- result$check[result$status == "FAIL"]
    stop(sprintf("Downstream handoff validation failed: %s", paste(failed, collapse = ", ")), call. = FALSE)
  }
  result
}


# Archived function: validate_stage_files
validate_stage_files <- function(root, required, rds = character(), stop_on_error = TRUE) {
  paths <- file.path(root, required)
  missing <- paths[!file.exists(paths)]
  valid_rds <- vapply(file.path(root, rds), function(path) {
    file.exists(path) && !inherits(tryCatch(readRDS(path), error = identity), "error")
  }, logical(1))
  passed <- !length(missing) && all(valid_rds)
  result <- data.frame(
    check = c("required_files", "reloadable_rds"),
    status = c(if (!length(missing)) "PASS" else "FAIL", if (all(valid_rds)) "PASS" else "FAIL"),
    details = c(
      if (!length(missing)) paste("files=", length(required), sep = "") else paste(basename(missing), collapse = ";"),
      if (all(valid_rds)) paste("rds=", length(rds), sep = "") else paste(basename(file.path(root, rds)[!valid_rds]), collapse = ";")
    ), stringsAsFactors = FALSE
  )
  if (isTRUE(stop_on_error) && !passed) stop(sprintf("Stage artifact validation failed under %s.", root), call. = FALSE)
  result
}


# Archived function: write_consensus_artifacts
write_consensus_artifacts <- function(project_root, output_root, consensus) {
  assert_path_within(project_root, output_root); dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
  saveRDS(consensus$consensus_object, file.path(output_root, "eligible_consensus_object.rds"), compress = FALSE)
  saveRDS(consensus$final_reference, file.path(output_root, "final_frozen_reference.rds"), compress = FALSE)
  write_tsv(consensus$gates, file.path(output_root, "consensus_gates.tsv"), project_root)
  write_tsv(consensus$marker_support, file.path(output_root, "consensus_marker_support.tsv"), project_root)
  write_tsv(data.frame(
    decision = consensus$decision, eligible_regions = paste(consensus$eligible_regions, collapse = ";"),
    admitted_regions = paste(consensus$admitted_regions, collapse = ";"),
    primary_reduction = "UNCORRECTED_PCA_PRIMARY", harmony = consensus$harmony_status,
    interpretation = consensus$interpretation, stringsAsFactors = FALSE
  ), file.path(output_root, "consensus_summary.tsv"), project_root)
  validate_stage_files(
    output_root,
    c("eligible_consensus_object.rds", "final_frozen_reference.rds", "consensus_gates.tsv",
      "consensus_marker_support.tsv", "consensus_summary.tsv"),
    c("eligible_consensus_object.rds", "final_frozen_reference.rds"), TRUE
  )
}


# Archived function: write_mapping_calibration_artifacts
write_mapping_calibration_artifacts <- function(project_root, output_root, calibration) {
  assert_path_within(project_root, output_root); dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
  write_tsv(calibration$folds, file.path(output_root, "region3_mapping_folds.tsv"), project_root)
  write_tsv(calibration$mapping, file.path(output_root, "region3_mapping_calibration.tsv"), project_root)
  write_tsv(calibration$thresholds, file.path(output_root, "mapping_thresholds.tsv"), project_root)
  write_tsv(calibration$summary, file.path(output_root, "mapping_calibration_summary.tsv"), project_root)
  saveRDS(calibration, file.path(output_root, "mapping_calibration.rds"), compress = FALSE)
  validate_stage_files(
    output_root,
    c("region3_mapping_folds.tsv", "region3_mapping_calibration.tsv", "mapping_thresholds.tsv",
      "mapping_calibration_summary.tsv", "mapping_calibration.rds"),
    "mapping_calibration.rds", TRUE
  )
}


# Archived function: write_region_admission_artifacts
write_region_admission_artifacts <- function(project_root, output_root, result) {
  assert_path_within(project_root, output_root); dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
  write_tsv(result$primary_mapping, file.path(output_root, "primary_mapping.tsv"), project_root)
  write_tsv(result$strict_mapping, file.path(output_root, "strict_mapping.tsv"), project_root)
  write_tsv(result$conservative_mapping, file.path(output_root, "conservative_mapping.tsv"), project_root)
  write_tsv(result$marker_coherence, file.path(output_root, "marker_coherence.tsv"), project_root)
  write_tsv(result$admission$gates, file.path(output_root, "admission_gates.tsv"), project_root)
  write_tsv(result$admission$summary, file.path(output_root, "admission_summary.tsv"), project_root)
  lightweight <- result
  lightweight$primary_query <- NULL
  for (name in c("primary_mapping", "strict_mapping", "conservative_mapping")) {
    attr(lightweight[[name]], "anchors") <- NULL
    attr(lightweight[[name]], "mapped_query") <- NULL
  }
  saveRDS(lightweight, file.path(output_root, "admission_result.rds"), compress = FALSE)
  validate_stage_files(
    output_root,
    c("primary_mapping.tsv", "strict_mapping.tsv", "conservative_mapping.tsv",
      "marker_coherence.tsv", "admission_gates.tsv", "admission_summary.tsv", "admission_result.rds"),
    "admission_result.rds", TRUE
  )
}


# Archived function: write_region3_anchor_artifacts
write_region3_anchor_artifacts <- function(project_root, output_root, anchor) {
  assert_path_within(project_root, output_root)
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
  paths <- c(
    write_tsv(anchor$branch_definitions, file.path(output_root, "region3_anchor_branches.tsv"), project_root),
    write_tsv(anchor$validation$gate_table, file.path(output_root, "region3_anchor_gates.tsv"), project_root),
    write_tsv(anchor$validation$major_label_counts, file.path(output_root, "region3_major_label_counts.tsv"), project_root)
  )
  for (branch in names(anchor$branches)) {
    branch_root <- file.path(output_root, branch)
    dir.create(branch_root, recursive = TRUE, showWarnings = FALSE)
    saveRDS(anchor$branches[[branch]]$object, file.path(branch_root, "reference_object.rds"), compress = FALSE)
    paths <- c(
      paths,
      write_tsv(anchor$branches[[branch]]$cell_labels, file.path(branch_root, "cell_labels.tsv"), project_root),
      write_tsv(anchor$branches[[branch]]$marker_support, file.path(branch_root, "marker_support.tsv"), project_root),
      write_tsv(anchor$branches[[branch]]$marker_scores, file.path(branch_root, "marker_scores.tsv"), project_root)
    )
  }
  summary <- list(
    schema_version = "region3_anchor_v2", validation = anchor$validation,
    branch_definitions = anchor$branch_definitions,
    primary_object_path = file.path(output_root, "primary_245", "reference_object.rds"),
    generated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
  )
  saveRDS(anchor$validation, file.path(output_root, "region3_anchor_validation.rds"))
  saveRDS(summary, file.path(output_root, "region3_anchor_summary.rds"))
  validate_stage_files(
    output_root,
    c("region3_anchor_branches.tsv", "region3_anchor_gates.tsv", "region3_major_label_counts.tsv",
      "region3_anchor_validation.rds", "region3_anchor_summary.rds",
      file.path("primary_245", "reference_object.rds"),
      file.path("complete_479", "reference_object.rds")),
    c("region3_anchor_validation.rds", "region3_anchor_summary.rds",
      file.path("primary_245", "reference_object.rds"),
      file.path("complete_479", "reference_object.rds")), TRUE
  )
  invisible(paths)
}

# End of archived functions.
