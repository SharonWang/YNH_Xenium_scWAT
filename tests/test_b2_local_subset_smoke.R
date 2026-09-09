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
stopifnot(
  checkpoint$validation_status == "PASS",
  checkpoint$n_cells == ncol(object),
  checkpoint$n_features == nrow(object),
  nzchar(checkpoint$md5)
)
cat(
  "LOCAL_SMOKE_PASS cells=", ncol(object),
  " genes=", nrow(object),
  " pcs=", ncol(Embeddings(object, "pca")),
  " median_ari=", stability$summary$median_ari,
  " max_abs_pc_qc=", max(pc_qc$abs_rho, na.rm = TRUE),
  "\n", sep = ""
)
