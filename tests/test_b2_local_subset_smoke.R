#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
test_path <- normalizePath(sub("^--file=", "", script_arg[[1L]]), winslash = "/", mustWork = TRUE)
repo_root <- dirname(dirname(test_path))
source(file.path(repo_root, "R", "source.R"))

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

input_root <- Sys.getenv(
  "SCWAT_LOCAL_SUBSET_INPUT_ROOT",
  "D:/Xiaonan/CODEX_projects/Yanan_Xenium/adipose_analysis/subset_input/adipose_data"
)
qc_root <- Sys.getenv(
  "SCWAT_LOCAL_SUBSET_QC_ROOT",
  "D:/Xiaonan/CODEX_projects/Yanan_Xenium/adipose_analysis/scwat_qc_outputs/local_colon_parity_all/sections/Region_1"
)
if (!dir.exists(input_root) || !file.exists(file.path(qc_root, "cell_qc_metadata.tsv.gz"))) {
  cat("SKIP: local D-drive subset or matching Region 1 QC mask is unavailable.\n")
  quit(status = 0L)
}

set.seed(1234L)
region_dir <- discover_one_section(input_root, "Region_1")$region_dir[[1L]]
object <- create_spatial_seurat_from_xenium(
  region_dir, project = "Region_1", assay = "Xenium", fov = "fov",
  include_cell_segmentation = TRUE, include_nucleus_segmentation = TRUE
)
masks <- read.delim(file.path(qc_root, "cell_qc_metadata.tsv.gz"), check.names = FALSE)
object <- add_masks_to_seurat(object, masks, require_complete_match = TRUE)
object$primary_include_revised <- derive_primary_include_revised(object@meta.data)
object <- subset(object, cells = colnames(object)[object$primary_include_revised])
stopifnot(nrow(object) == 479L)

object <- NormalizeData(
  object, assay = "Xenium", normalization.method = "LogNormalize", verbose = FALSE
)
object <- ScaleData(object, assay = "Xenium", features = rownames(object), verbose = FALSE)
object <- RunPCA(
  object, assay = "Xenium", features = rownames(object), npcs = 30,
  seed.use = 1234L, verbose = FALSE
)
pc_qc <- summarise_pca_qc_correlations(Embeddings(object, "pca"), object@meta.data)
object <- FindNeighbors(
  object, reduction = "pca", dims = 1:30, k.param = 20,
  graph.name = c("subset_nn", "subset_snn"), verbose = FALSE
)
for (seed in c(1234L, 2026L, 31415L)) {
  object <- FindClusters(
    object, graph.name = "subset_snn", algorithm = 1, resolution = 0.8,
    cluster.name = paste0("seed_", seed), random.seed = seed, verbose = FALSE
  )
}
stability <- summarise_cluster_stability(
  object@meta.data[, c("seed_1234", "seed_2026", "seed_31415"), drop = FALSE]
)
checkpoint_dir <- file.path(tempdir(), "b2_checkpoint_contract")
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
on.exit(unlink(checkpoint_dir, recursive = TRUE, force = TRUE), add = TRUE)
checkpoint <- write_validated_seurat_checkpoint(
  object = object,
  path = file.path(checkpoint_dir, "subset_pca.rds"),
  project_root = tempdir(),
  stage = "LOCAL_SUBSET_PCA"
)

# Exercise the complete-panel notebook's ordering and corrected spatial safety
# contracts without pretending that this subset has biological annotations.
synthetic_labels <- rep(c("T", "Adipocyte", "Eosinophil", "Macrophage"), length.out = ncol(object))
ordered_labels <- apply_scwat_cell_type_order(synthetic_labels)
stopifnot(
  identical(levels(ordered_labels), c("Adipocyte", "Macrophage", "Eosinophil", "T")),
  identical(as.character(ordered_labels), synthetic_labels)
)
coordinates <- Seurat::GetTissueCoordinates(object[["fov"]], which = "centroids") |> as.data.frame()
if (!"cell" %in% colnames(coordinates)) coordinates$cell <- rownames(coordinates)
spatial_metadata <- data.frame(
  cell_id = colnames(object),
  Eos_inclusive = colnames(object) %in% head(colnames(object), 20L),
  Final_CellType_subtype = ifelse(
    colnames(object) %in% head(colnames(object), 20L),
    "Eosinophil", rep(c("Adipocyte", "Macrophage"), length.out = ncol(object))
  ),
  EosState_balance = seq(-1, 1, length.out = ncol(object)),
  EosState_extreme = NA_character_,
  stringsAsFactors = FALSE
)
spatial_pools <- build_eos_spatial_pools(
  spatial_metadata,
  data.frame(cell_id = coordinates$cell, x = coordinates$x, y = coordinates$y),
  cell_type_col = "Final_CellType_subtype"
)
spatial_edges <- calculate_eos_knn_edges(spatial_pools, k_values = c(1L, 15L))
optional_mclust <- run_mclust_diagnostic(1:10, min_n = 20L)
optional_cellchat <- run_eos_spatial_cellchat(list(status = "SKIPPED_LOCAL_SUBSET"))
stopifnot(
  checkpoint$validation_status == "PASS",
  checkpoint$n_cells == ncol(object),
  checkpoint$n_features == nrow(object),
  nzchar(checkpoint$md5),
  spatial_pools$status == "PASS",
  nrow(spatial_edges$k1) == 20L,
  nrow(spatial_edges$k15) == 300L,
  all(spatial_edges$k15$distance > 0),
  optional_mclust$status %in% c("SKIPPED_PACKAGE_UNAVAILABLE", "SKIPPED_INSUFFICIENT_DATA"),
  optional_cellchat$status == "SKIPPED_LOCAL_SUBSET"
)
cat(
  "LOCAL_SMOKE_PASS cells=", ncol(object),
  " genes=", nrow(object),
  " pcs=", ncol(Embeddings(object, "pca")),
  " median_ari=", stability$summary$median_ari,
  " max_abs_pc_qc=", max(pc_qc$abs_rho, na.rm = TRUE),
  "\n", sep = ""
)
