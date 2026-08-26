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

# -----------------------------------------------------------------------------
# Archived 2026-08-27 after Colon-parity initial-QC notebook replacement
# -----------------------------------------------------------------------------

# Archived function: assign_spatial_grid

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

# Archived function: build_cycle_alarm_evidence

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

# Archived function: build_eos_gene_decision

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

# Archived function: build_evidence_only_release

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

# Archived function: build_gene_downstream_decision

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

# Archived function: build_hotspot_sensitivity_decision

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

# Archived function: build_section_downstream_decision

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

# Archived function: build_section_pairs

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

# Archived function: build_transcript_quality_queries

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

# Archived function: calculate_knn_density

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

# Archived function: calculate_readiness_gates

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

# Archived function: calculate_within_mouse_concordance

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

# Archived function: combine_gene_quality

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

# Archived function: compare_subset_full_qc

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

# Archived function: create_synthetic_manifest

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

# Archived function: evidence_only_required_artifacts

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

# Archived function: extended_qc_preflight

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

# Archived function: extended_section_required_artifacts

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

# Archived function: extended_slide_required_artifacts

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

# Archived function: extract_analysis_alarms

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

# Archived function: find_spatial_qc_hotspots

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

# Archived function: knn_index_distance

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

# Archived function: plot_extended_slide_qc

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

# Archived function: plot_extended_spatial_qc

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

# Archived function: quantile_distribution_distance

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

# Archived function: rank_candidate_cycle_genes

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

# Archived function: rbind_fill

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

# Archived function: read_custom_panel_genes

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

# Archived function: read_extended_qc_config

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

# Archived function: read_extended_section_artifacts

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

# Archived function: read_extended_slide_qc_outputs

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

# Archived function: read_slide_qc_outputs

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

#' Summarise four-region scWAT slide QC evidence.
#'
#' @param slide_data Combined output from [read_scwat_slide_qc_outputs()].
#' @return A list containing section-level pass/review fractions, validated
#'   technical readiness for each region, and the worst-case overall status.
#'   The function does not redefine any cell mask.
#' @export

# Archived function: reconcile_panel

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

# Archived function: resolve_extended_qc_mode

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

# Archived function: resolve_transcript_schema

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

# Archived function: save_extended_spatial_plots

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

# Archived function: section_required_artifacts

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

# Archived function: slide_section_required_files

# Purpose: Slide section required files.
# Inputs: none.
# Output: Returns the derived R object described by the function name; no files are written unless an explicit output path is an input.
slide_section_required_files <- function() {
  c("qc_summary.tsv", "qc_thresholds.tsv", "section_readiness_gates.tsv", "analysis_alerts.tsv", "cell_qc_metadata.tsv.gz")
}

# Purpose: Validate four section outputs.
# Inputs: required: run_root; optional/defaulted: expected_regions.
# Output: Returns validation evidence/TRUE (or the validated value) and stops with an informative error when the contract fails.

# Archived function: summarise_evidence_only_qc

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

# Archived function: summarise_extended_slide_qc

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

# Archived function: summarise_gene_matrix_qc

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

# Archived function: summarise_spatial_enrichment

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

# Archived function: summarise_transcript_quality_arrow

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

# Archived function: summarise_transcript_quality_table

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

# Archived function: test_spatial_flag_clustering

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

# Archived function: validate_evidence_only_qc_artifacts

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

# Archived function: validate_extended_section_artifacts

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

# Archived function: validate_extended_slide_qc_artifacts

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

# Archived function: validate_four_extended_section_outputs

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

# Archived function: validate_four_section_outputs

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

# Archived function: validate_section_artifacts

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

# Archived function: validate_spatial_cells

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

# Archived function: write_evidence_only_qc_artifacts

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

# Archived function: write_extended_section_artifacts

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

# Archived function: write_extended_slide_qc_artifacts

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

# Archived function: write_section_artifacts

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

# Archived function: write_slide_qc_artifacts

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
