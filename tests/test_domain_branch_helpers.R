#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
test_path <- normalizePath(sub("^--file=", "", script_arg[[1L]]), winslash = "/", mustWork = TRUE)
repo_root <- dirname(dirname(test_path))
source(file.path(repo_root, "R", "source.R"))

expect_error <- function(expr, pattern) {
  error <- tryCatch({ force(expr); NULL }, error = identity)
  stopifnot(inherits(error, "error"), grepl(pattern, conditionMessage(error), ignore.case = TRUE))
}

# Break caught: a row-order join or incomplete partition silently puts a cell
# into the wrong tissue branch.
metadata <- data.frame(
  cell_id = paste0("cell", 1:6),
  primary_include_revised = c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE),
  stringsAsFactors = FALSE
)
manifest <- build_tissue_branch_manifest(
  cell_metadata = metadata,
  lymph_node_cell_ids = c("cell3", "cell2")
)
stopifnot(
  identical(manifest$cell_id, metadata$cell_id),
  identical(manifest$lymph_node_include, c(FALSE, TRUE, TRUE, FALSE, FALSE, FALSE)),
  identical(manifest$all_qcpass_include, c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE)),
  identical(manifest$adipose_only_include, c(TRUE, FALSE, FALSE, TRUE, TRUE, FALSE)),
  identical(manifest$lymph_node_only_include, c(FALSE, TRUE, TRUE, FALSE, FALSE, FALSE)),
  all((manifest$adipose_only_include | manifest$lymph_node_only_include) == manifest$all_qcpass_include)
)
stopifnot(identical(select_tissue_branch_ids(manifest, "all_qcpass", min_cells = 1L), paste0("cell", 1:5)))
stopifnot(identical(select_tissue_branch_ids(manifest, "adipose_only", min_cells = 1L), c("cell1", "cell4", "cell5")))
stopifnot(identical(select_tissue_branch_ids(manifest, "lymph_node_only", min_cells = 1L), c("cell2", "cell3")))
expect_error(select_tissue_branch_ids(manifest, "lymph_node_only", min_cells = 3L), "insufficient")
expect_error(build_tissue_branch_manifest(rbind(metadata, metadata[1, ]), c("cell2")), "duplicate")
expect_error(build_tissue_branch_manifest(metadata, c("not-a-cell")), "not present")

# Break caught: rescued Eos cells leak into the non-Eos reference pool, or
# duplicate coordinates create zero-distance matches despite disjoint IDs.
pool_gate <- validate_disjoint_spatial_pools(
  query_ids = c("e1", "e2"),
  reference_ids = c("r1", "r2"),
  query_coordinates = data.frame(cell_id = c("e1", "e2"), x = c(0, 1), y = c(0, 1)),
  reference_coordinates = data.frame(cell_id = c("r1", "r2"), x = c(2, 3), y = c(2, 3))
)
stopifnot(pool_gate$n_overlap_ids == 0L, pool_gate$n_duplicate_coordinate_pairs == 0L)
expect_error(validate_disjoint_spatial_pools(c("e1", "e2"), c("e2", "r1")), "overlap")
expect_error(validate_disjoint_spatial_pools(
  c("e1"), c("r1"),
  data.frame(cell_id = "e1", x = 0, y = 0),
  data.frame(cell_id = "r1", x = 0, y = 0)
), "coordinate")

distance_gate <- validate_neighbour_distances(c(0.2, 1.5, 3), zero_tolerance = 0)
stopifnot(distance_gate$n_zero_or_negative == 0L, distance_gate$minimum_distance == 0.2)
expect_error(validate_neighbour_distances(c(0, 1), zero_tolerance = 0), "zero")

# Break caught: clustering stability is asserted from visual inspection rather
# than a label-invariant statistic.
stopifnot(
  identical(adjusted_rand_index(c(1, 1, 2, 2), c("A", "A", "B", "B")), 1),
  identical(adjusted_rand_index(c(1, 1, 2, 2), c(2, 2, 1, 1)), 1),
  isTRUE(all.equal(adjusted_rand_index(c(1, 1, 2, 2), c(1, 2, 1, 2)), -0.5))
)
assignments <- data.frame(
  seed_1 = c(1, 1, 2, 2), seed_2 = c(2, 2, 1, 1), seed_3 = c(1, 2, 1, 2)
)
stability <- summarise_cluster_stability(assignments)
stopifnot(nrow(stability$pairwise) == 3L, stability$summary$n_comparisons == 3L)

# Break caught: a technically dominated PC passes without a visible numerical
# diagnostic. Expected correlations are derived directly from the tiny fixture.
embeddings <- cbind(PC_1 = 1:4, PC_2 = c(1, -1, 1, -1))
qc_metadata <- data.frame(nCount_Xenium = 1:4, nFeature_Xenium = 4:1)
pc_qc <- summarise_pca_qc_correlations(embeddings, qc_metadata)
stopifnot(pc_qc$rho[pc_qc$pc == "PC_1" & pc_qc$qc_metric == "nCount_Xenium"] == 1)
stopifnot(pc_qc$rho[pc_qc$pc == "PC_1" & pc_qc$qc_metric == "nFeature_Xenium"] == -1)

# Break caught: the reusable Eosinophil composition plot depends on a
# notebook-local `theme_cell()` helper instead of functions sourced from
# R/source.R. The split notebooks attach dplyr during setup, so reproduce that
# execution context while deliberately leaving `theme_cell()` undefined.
suppressPackageStartupMessages(library(dplyr))
plot_counts <- Matrix::Matrix(
  matrix(c(2, 0, 1, 0, 3, 1), nrow = 2L),
  sparse = TRUE
)
rownames(plot_counts) <- c("Siglecf", "Adgre1")
colnames(plot_counts) <- c("eos1", "mac1", "other1")
plot_object <- Seurat::CreateSeuratObject(counts = plot_counts, assay = "Xenium")
plot_object$Final_CellType_subtype <- c("Eosinophil", "Macrophage", "Stromal")
plot_object$EosRef_call <- c("REF_EOS_TIER1", "REF_EOS_REST", "OUTSIDE_IMMUNE")
eos_composition <- plot_eos_call_by_subtype(plot_object, return_data = TRUE)
stopifnot(
  inherits(eos_composition$plot, "ggplot"),
  nrow(eos_composition$plot_data) == 21L,
  identical(eos_composition$subtype_order[[1L]], "Eosinophil")
)

# Break caught: a lymph-node domain is forced from isolated lymphoid calls or
# a tissue boundary is silently unstable across reasonable parameters.
ln_fixture <- data.frame(
  cell_id = paste0("c", 1:20),
  x = c(seq(0, 0.9, length.out = 10), seq(100, 100.9, length.out = 10)),
  y = rep(seq(0, 0.9, length.out = 10), 2),
  cell_type = c(rep(c("B", "T"), 5), rep("Adipocyte", 10)),
  stringsAsFactors = FALSE
)
ln_domain <- derive_lymph_node_domain(
  cell_metadata = ln_fixture,
  lymphoid_labels = c("B", "T"),
  k = 3L,
  lymphoid_fraction_threshold = 2 / 3,
  expansion_radius = 2,
  min_core_cells = 3L
)
stopifnot(
  identical(ln_domain$status, "LN_CANDIDATE"),
  all(ln_domain$cell_table$lymph_node_include[1:10]),
  !any(ln_domain$cell_table$lymph_node_include[11:20])
)
ln_sensitivity <- summarise_ln_boundary_sensitivity(list(
  primary = ln_domain$cell_table$lymph_node_include,
  identical = ln_domain$cell_table$lymph_node_include,
  narrower = c(rep(TRUE, 8), rep(FALSE, 12))
))
stopifnot(nrow(ln_sensitivity$pairwise) == 3L, ln_sensitivity$summary$minimum_jaccard == 0.8)

cat("Tissue-branch and spatial-safety helper tests passed.\n")
