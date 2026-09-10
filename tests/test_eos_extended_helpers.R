#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
test_path <- normalizePath(sub("^--file=", "", script_arg[[1L]]), winslash = "/", mustWork = TRUE)
repo_root <- dirname(dirname(test_path))
source(file.path(repo_root, "R", "source.R"))

# Break caught: mclust::Mclust() evaluates an unqualified mclustBIC call in
# its caller, so a namespace-only invocation fails unless the exported BIC
# function is deliberately bound in that caller frame.
fake_mclust_bic <- function(data, G = NULL, verbose = FALSE, ...) {
  structure(matrix(c(-10, -8), nrow = 2L), G = G, modelNames = "V")
}
fake_mclust <- function(data, G = NULL, verbose = FALSE, ...) {
  mc <- match.call(expand.dots = TRUE)
  mc[[1L]] <- as.name("mclustBIC")
  mc[[2L]] <- data
  bic <- eval(mc, parent.frame())
  list(G = 2L, modelName = "V", BIC = bic)
}
core_fit <- run_mclust_with_binding(
  data = seq(-1, 1, length.out = 40L),
  G = 1:2,
  mclust_fun = fake_mclust,
  mclust_bic_fun = fake_mclust_bic
)
stopifnot(core_fit$G == 2L, core_fit$modelName == "V")

diagnostic <- run_mclust_diagnostic(seq(-1, 1, length.out = 40L), G = 1:3)
if (requireNamespace("mclust", quietly = TRUE)) {
  stopifnot(
    diagnostic$status == "PASS",
    diagnostic$n_finite == 40L,
    diagnostic$selected_G %in% 1:3,
    is.data.frame(diagnostic$bic_table)
  )
} else {
  stopifnot(diagnostic$status == "SKIPPED_PACKAGE_UNAVAILABLE")
}

small <- run_mclust_diagnostic(1:10, min_n = 20L)
if (requireNamespace("mclust", quietly = TRUE)) {
  stopifnot(small$status == "SKIPPED_INSUFFICIENT_DATA")
} else {
  stopifnot(small$status == "SKIPPED_PACKAGE_UNAVAILABLE")
}

# Break caught: plots silently fall back to alphabetical ordering or assign
# unstable colours when a section contains a previously unseen subtype.
labels <- c("T", "Adipocyte", "Eosinophil", "Capillary_EC", "new_type")
ordered_labels <- apply_scwat_cell_type_order(labels)
stopifnot(
  identical(
    levels(ordered_labels),
    c("Adipocyte", "Capillary_EC", "Eosinophil", "T", "new_type")
  ),
  identical(as.character(ordered_labels), labels)
)
palette <- cell_macaron_palette(labels)
stopifnot(
  setequal(names(palette), unique(labels)),
  all(grepl("^#[0-9A-Fa-f]{6}$", unname(palette)))
)
stopifnot(identical(palette, cell_macaron_palette(labels)))

overlay_fixture <- data.frame(
  x = 1:5, y = 5:1,
  group = c("Other", "Short-lived-like", "Other", "Long-lived-like", "Other"),
  stringsAsFactors = FALSE
)
overlay_plot <- plot_target_overlay(
  overlay_fixture, x_col = "x", y_col = "y", group_col = "group",
  target_labels = c("Short-lived-like", "Long-lived-like"),
  palette = c("Other" = "#D9D9D9", "Short-lived-like" = "#3979A8", "Long-lived-like" = "#B84E4B"),
  background_size = 0.5, target_size = 2
)
stopifnot(
  inherits(overlay_plot, "ggplot"),
  length(overlay_plot$layers) == 2L,
  nrow(overlay_plot$layers[[1L]]$data) == 3L,
  nrow(overlay_plot$layers[[2L]]$data) == 2L,
  overlay_plot$layers[[1L]]$aes_params$size == 0.5,
  overlay_plot$layers[[2L]]$aes_params$size == 2,
  is.null(overlay_plot$labels$subtitle)
)

markers <- data.frame(
  Gene_Symbol = c("Cd3d", "Pck1", "Siglecf", "Kdr"),
  CellType_subtype = c("T", "Adipocyte", "Eosinophil", "Capillary_EC"),
  stringsAsFactors = FALSE
)
marker_order <- order_marker_features(markers, available_genes = markers$Gene_Symbol)
stopifnot(
  identical(marker_order$marker_table$Gene_Symbol, c("Pck1", "Kdr", "Siglecf", "Cd3d")),
  identical(names(marker_order$feature_groups), c("Adipocyte", "Capillary_EC", "Eosinophil", "T"))
)

plain_plot <- ggplot2::ggplot(data.frame(x = 1:3, y = 1:3), ggplot2::aes(x, y)) +
  ggplot2::geom_point()
styled_plot <- style_cell_plot(plain_plot)
stopifnot(inherits(styled_plot, "ggplot"))
plot_output_dir <- file.path(Sys.getenv("TMPDIR"), "cell_plot_contract")
plot_files <- save_cell_plot(
  styled_plot, stem = "test_plot", output_dir = plot_output_dir,
  project_root = Sys.getenv("TMPDIR"), width = 4, height = 3, dpi = 72
)
stopifnot(
  setequal(names(plot_files), c("png", "pdf")),
  all(file.exists(unname(plot_files)))
)

# Break caught: write.table() fails on list/matrix-valued columns returned by
# optional packages, preventing the final audit-output cell from completing.
nested_table <- data.frame(id = c("a", "b"), stringsAsFactors = FALSE)
nested_table$matrix_metric <- I(matrix(1:4, nrow = 2, dimnames = list(NULL, c("x", "y"))))
nested_table$list_metric <- I(list(c("u", "v"), "w"))
nested_path <- file.path(Sys.getenv("TMPDIR"), "nested_table.tsv")
write_tsv(nested_table, nested_path, Sys.getenv("TMPDIR"))
nested_roundtrip <- read.delim(nested_path, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(
  nrow(nested_roundtrip) == 2L,
  all(c("matrix_metric.x", "matrix_metric.y", "list_metric") %in% colnames(nested_roundtrip)),
  nested_roundtrip$list_metric[[1L]] == "u | v"
)

# Break caught on HPC: some optional-package outputs contain a nested data-frame
# column whose internal row count differs from the enclosing table. The writer
# must retain the outer rows and preserve that nested payload once, rather than
# aborting the complete final-save checkpoint.
mismatched_nested_table <- structure(
  list(
    id = c("a", "b"),
    x = I(data.frame(component = 1:3, score = c(0.1, 0.2, 0.3)))
  ),
  class = "data.frame",
  row.names = .set_row_names(2L)
)
mismatched_nested_path <- file.path(Sys.getenv("TMPDIR"), "mismatched_nested_table.tsv")
write_tsv(mismatched_nested_table, mismatched_nested_path, Sys.getenv("TMPDIR"))
mismatched_nested_roundtrip <- read.delim(
  mismatched_nested_path, check.names = FALSE, stringsAsFactors = FALSE
)
stopifnot(
  nrow(mismatched_nested_roundtrip) == 2L,
  all(c(
    "x.__nested_serialized__", "x.__nested_rows__", "x.__nested_cols__"
  ) %in% colnames(mismatched_nested_roundtrip)),
  grepl("component", mismatched_nested_roundtrip$x.__nested_serialized__[[1L]], fixed = TRUE),
  is.na(mismatched_nested_roundtrip$x.__nested_serialized__[[2L]]),
  mismatched_nested_roundtrip$x.__nested_rows__[[1L]] == 3L,
  mismatched_nested_roundtrip$x.__nested_cols__[[1L]] == 2L
)

# Break caught: spatial pools are joined by row position, Eosinophils leak into
# the reference pool, or k-neighbour edge counts/distances are incorrect.
spatial_fixture <- data.frame(
  cell_id = c("e1", "e2", "a1", "a2", "m1", "m2", "t1", "t2"),
  Eos_inclusive = c(TRUE, TRUE, rep(FALSE, 6)),
  Final_CellType_subtype = c(
    "Eosinophil", "Eosinophil", "ASC", "ASC",
    "Macrophage", "Macrophage", "T", "T"
  ),
  EosState_balance = c(-1, 1, rep(NA_real_, 6)),
  EosState_extreme = c("Short-lived-like", "Long-lived-like", rep(NA_character_, 6)),
  stringsAsFactors = FALSE
)
coordinate_fixture <- data.frame(
  cell_id = rev(spatial_fixture$cell_id),
  x = c(13, 3, 12, 2, 11, 1, 10, 0),
  y = 0,
  stringsAsFactors = FALSE
)
pools <- build_eos_spatial_pools(spatial_fixture, coordinate_fixture)
stopifnot(
  identical(pools$query$cell_id, c("e1", "e2")),
  pools$gate$n_overlap_ids == 0L,
  pools$gate$n_duplicate_coordinate_pairs == 0L
)
empty_eos_fixture <- spatial_fixture
empty_eos_fixture$Eos_inclusive <- FALSE
empty_pools <- build_eos_spatial_pools(empty_eos_fixture, coordinate_fixture)
stopifnot(
  empty_pools$status == "SKIPPED_EMPTY_EOS_OR_REFERENCE_POOL",
  nrow(empty_pools$query) == 0L,
  nrow(empty_pools$reference) == nrow(empty_eos_fixture)
)
edges <- calculate_eos_knn_edges(pools, k_values = c(1L, 3L))
stopifnot(
  nrow(edges$k1) == 2L,
  nrow(edges$k3) == 6L,
  all(edges$k1$distance > 0),
  identical(edges$k1$reference_cell_type, c("ASC", "ASC"))
)
composition <- summarise_eos_knn_composition(edges$k3)
stopifnot(sum(composition$overall$n_edges) == 6L)

distance_by_type <- calculate_eos_distance_by_cell_type(pools, min_reference_cells = 2L)
stopifnot(
  nrow(distance_by_type$cell_level) == 6L,
  setequal(distance_by_type$summary$reference_cell_type, c("ASC", "Macrophage", "T")),
  "EosState_extreme" %in% colnames(distance_by_type$cell_level),
  setequal(unique(distance_by_type$cell_level$EosState_extreme), c("Short-lived-like", "Long-lived-like"))
)

association_edges <- data.frame(
  eos_cell_id = rep(c("e1", "e2", "e3"), each = 3L),
  reference_cell_type = c("ASC", "ASC", "T", "ASC", "Macrophage", "T", "Macrophage", "Macrophage", "T"),
  EosState_balance = rep(c(-1, 0, 1), each = 3L),
  EosState_extreme = rep(c("Short-lived-like", "Intermediate", "Long-lived-like"), each = 3L),
  stringsAsFactors = FALSE
)
ranking <- rank_eos_state_knn_associations(
  association_edges,
  biological_order = c("ASC", "Macrophage", "T"),
  min_eos = 3L,
  top_n = 2L
)
stopifnot(
  ranking$full$spearman_rho[ranking$full$reference_cell_type == "ASC"] == -1,
  ranking$full$spearman_rho[ranking$full$reference_cell_type == "Macrophage"] == 1,
  identical(ranking$short_top[[1L]], "ASC"),
  identical(ranking$long_top[[1L]], "Macrophage"),
  "EosState_extreme" %in% colnames(ranking$per_eos),
  all(ranking$per_eos$EosState_extreme[ranking$per_eos$eos_cell_id == "e1"] == "Short-lived-like")
)

# Break caught: continuous Eosinophil scores are passed directly to CellChat,
# underpowered tails are forced, or the reported significant table ignores
# multiplicity and sender/receiver direction.
cellchat_groups <- derive_eos_cellchat_groups(1:100, min_cells = 10L)
stopifnot(
  cellchat_groups$status == "PASS",
  sum(cellchat_groups$group == "Eos_short_enriched", na.rm = TRUE) == 30L,
  sum(cellchat_groups$group == "Eos_long_enriched", na.rm = TRUE) == 30L
)
too_small_groups <- derive_eos_cellchat_groups(1:20, min_cells = 10L)
stopifnot(too_small_groups$status == "SKIPPED_INSUFFICIENT_STATE_GROUP_CELLS")

cellchat_counts <- Matrix::Matrix(
  matrix(rep(c(1, 2, 0, 3), 60L), nrow = 4L), sparse = TRUE
)
rownames(cellchat_counts) <- c("A", "B", "C", "D")
colnames(cellchat_counts) <- paste0("cc", seq_len(60L))
cellchat_object <- Seurat::CreateSeuratObject(cellchat_counts, assay = "Xenium")
cellchat_object <- Seurat::NormalizeData(cellchat_object, assay = "Xenium", verbose = FALSE)
cellchat_object$Eos_inclusive <- c(rep(TRUE, 40L), rep(FALSE, 20L))
cellchat_object$Final_CellType_subtype <- c(rep("Eosinophil", 40L), rep("ASC", 10L), rep("Macrophage", 10L))
cellchat_object$EosState_balance <- c(seq(-2, 2, length.out = 40L), rep(NA_real_, 20L))
cellchat_object$cell_area <- rep(100, 60L)
cellchat_coordinates <- data.frame(
  cell_id = rev(colnames(cellchat_object)), x = rev(seq_len(60L)), y = 0,
  stringsAsFactors = FALSE
)
cellchat_inputs <- prepare_eos_cellchat_inputs(
  cellchat_object, cellchat_coordinates,
  top_short = c("ASC"), top_long = c("Macrophage"), min_cells = 5L
)
stopifnot(
  cellchat_inputs$status == "PASS",
  ncol(cellchat_inputs$data) == 44L,
  identical(rownames(cellchat_inputs$meta), colnames(cellchat_inputs$data)),
  cellchat_inputs$scale_factors$spot == 2 * sqrt(100 / pi),
  isTRUE(all.equal(cellchat_inputs$scale_factors$spot.diameter, 2 * sqrt(100 / pi))),
  cellchat_inputs$spatial_factors$ratio == 1,
  isTRUE(all.equal(cellchat_inputs$spatial_factors$tol, sqrt(100 / pi))),
  identical(as.character(cellchat_inputs$meta$samples), rep("sample1", nrow(cellchat_inputs$meta)))
)

old_create <- function(object, meta, group.by, datatype, coordinates, scale.factors) {
  list(api = "scale.factors", scale = scale.factors, coordinates = coordinates)
}
new_create <- function(object, meta, group.by, datatype, coordinates, spatial.factors) {
  list(api = "spatial.factors", scale = spatial.factors, coordinates = coordinates)
}
old_created <- create_cellchat_object_compatible(cellchat_inputs, create_fun = old_create)
new_created <- create_cellchat_object_compatible(cellchat_inputs, create_fun = new_create)
stopifnot(
  old_created$api == "scale.factors",
  new_created$api == "spatial.factors",
  old_created$scale$spot == 2 * sqrt(100 / pi),
  new_created$scale$ratio == 1
)

lr_fixture <- data.frame(
  source = c("Eos_short_enriched", "Eos_long_enriched", "ASC"),
  target = c("ASC", "Macrophage", "Eos_short_enriched"),
  interaction_name = c("A_B", "C_D", "E_F"),
  pathway_name = c("P1", "P2", "P3"),
  ligand = c("A", "C", "E"),
  receptor = c("B", "D", "F"),
  prob = c(0.4, 0.2, 0.3),
  pval = c(0.001, 0.002, 0.01),
  stringsAsFactors = FALSE
)
filtered_lr <- filter_eos_cellchat_interactions(
  lr_fixture, top_short = "ASC", top_long = "Macrophage"
)
stopifnot(
  nrow(filtered_lr$significant) == 2L,
  all(filtered_lr$significant$direction == "EOS_TO_NEIGHBOUR"),
  all(filtered_lr$significant$p_adjust_bh < 0.10),
  all(filtered_lr$significant$matches_state_top_neighbour)
)

cellchat_skip <- run_eos_spatial_cellchat(list(status = "SKIPPED_TEST_INPUT"))
stopifnot(cellchat_skip$status == "SKIPPED_TEST_INPUT")

# Break caught: Wang sampling is unbalanced/non-deterministic, cell IDs collide
# after merge, or shared-feature and cross-dataset Eosinophil concordance gates
# are bypassed.
sampled_ids <- sample_ids_by_group(
  ids = paste0("w", seq_len(12L)),
  groups = rep(c("Eosinophil", "ASC"), each = 6L),
  max_per_group = 3L,
  seed = 5L
)
stopifnot(
  length(sampled_ids) == 6L,
  identical(
    sampled_ids,
    sample_ids_by_group(
      paste0("w", seq_len(12L)),
      rep(c("Eosinophil", "ASC"), each = 6L),
      max_per_group = 3L,
      seed = 5L
    )
  )
)

shared <- validate_shared_feature_set(
  reference_genes = c("A", "B", "C"),
  query_genes = c("B", "C", "D"),
  panel_genes = c("A", "B", "C", "D"),
  min_shared = 2L
)
stopifnot(shared$status == "PASS", identical(shared$features, c("B", "C")))
shared_fail <- validate_shared_feature_set(c("A"), c("A"), c("A"), min_shared = 2L)
stopifnot(shared_fail$status == "SKIPPED_INSUFFICIENT_SHARED_GENES")

prefixed_object <- prefix_seurat_cell_ids(cellchat_object[, 1:3], "WANG")
stopifnot(identical(colnames(prefixed_object), paste0("WANG_", colnames(cellchat_object)[1:3])))

cross_embeddings <- rbind(
  WANG_eos = c(0, 0), WANG_asc = c(10, 0),
  XENIUM_eos = c(0.2, 0), XENIUM_asc = c(10.2, 0)
)
colnames(cross_embeddings) <- c("PC_1", "PC_2")
cross_metadata <- data.frame(
  dataset = c("WANG", "WANG", "XENIUM", "XENIUM"),
  is_eosinophil = c(TRUE, FALSE, TRUE, FALSE),
  integrated_cluster = c("0", "1", "0", "1"),
  row.names = rownames(cross_embeddings),
  stringsAsFactors = FALSE
)
cross_summary <- summarise_cross_dataset_eos_neighbours(
  cross_embeddings, cross_metadata, k = 1L
)
stopifnot(
  identical(cross_summary$per_cell$eos_neighbour_fraction, c(1, 0, 1, 0)),
  all(cross_summary$per_cell$minimum_cross_dataset_distance > 0),
  nrow(cross_summary$cluster_enrichment) == 4L
)

integration_gate <- run_wang_xenium_integration(
  reference = cellchat_object[, 1:10],
  query = cellchat_object[, 11:20],
  features = c("A", "B"),
  reference_assay = "Xenium",
  query_assay = "Xenium",
  min_shared = 10L
)
stopifnot(integration_gate$status == "SKIPPED_INSUFFICIENT_SHARED_GENES")

cat("Extended Eosinophil helper tests passed.\n")
