#!/usr/bin/env python3
"""Build the twelve output-free, full-panel scWAT branch notebooks.

The notebooks are reader-facing HPC workflows. Region and branch are locked in
each file; all other paths may be overridden through SCWAT_* environment
variables. This builder uses only the Python standard library.
"""

from __future__ import annotations

import json
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
NOTEBOOK_DIR = REPO / "notebooks"
BRANCHES = ("all_QCpass", "adipose_only", "lymph_node_only")


def split_source(text: str) -> list[str]:
    lines = text.strip("\n").splitlines(keepends=True)
    if lines and not lines[-1].endswith("\n"):
        lines[-1] += "\n"
    return lines


def make_cell(cell_type: str, source: str, index: int) -> dict:
    payload = {
        "cell_type": cell_type,
        "id": f"cell-{index:03d}",
        "metadata": {},
        "source": split_source(source),
    }
    if cell_type == "code":
        payload.update(execution_count=None, outputs=[])
    return payload


def make_notebook(cells: list[tuple[str, str]]) -> dict:
    return {
        "cells": [make_cell(kind, source, index) for index, (kind, source) in enumerate(cells, 1)],
        "metadata": {
            "kernelspec": {"display_name": "R", "language": "R", "name": "ir"},
            "language_info": {"name": "R", "mimetype": "text/x-r-source", "file_extension": ".r"},
        },
        "nbformat": 4,
        "nbformat_minor": 5,
    }


def md(text: str) -> tuple[str, str]:
    return ("markdown", text)


def code(text: str) -> tuple[str, str]:
    return ("code", text)


def notebook_cells(region_number: int, branch: str) -> list[tuple[str, str]]:
    region = f"Region_{region_number}"
    section_role = "SENSITIVITY_ONLY" if region_number == 4 else ("PRIMARY" if region_number == 3 else "PRIMARY_CONDITIONAL")
    branch_label = {
        "all_QCpass": "all cells passing the revised primary mask",
        "adipose_only": "QC-passed cells outside the frozen lymph-node domain",
        "lymph_node_only": "QC-passed cells inside the frozen lymph-node domain",
    }[branch]
    child_branch = branch != "all_QCpass"

    cells: list[tuple[str, str]] = [
        md(f"""# {region}: {branch_label} — complete 479-gene panel

## Goal

Run a clean-kernel, stepwise Xenium analysis of **{region}** using **{branch_label}**. Every numerical decision is displayed before it is used. Raw counts and the original QC object are not modified.

**Section role:** `{section_role}`. Region 3 is the anchor. Regions 1–2 require admission against that anchor. Region 4 is exploratory/mapping-only and must not define the final reference, PCA, integration, clustering, markers, or labels.

### Interpretation boundary

The 479-gene workflow is a full-panel exploratory section analysis. It must later be reconciled with the 245-gene provisional and 67-gene conservative sensitivity analyses. Cells, not mice, are sampled here; cell-level tests cannot establish mouse-level differential expression."""),
        md("""### Why Harmony is not used in this notebook

Each notebook contains one physical section and therefore has no defensible within-notebook batch factor for Harmony. Applying Harmony to clusters, cell types or tissue domains would risk removing biology. Harmony belongs only in the later frozen eligible-reference workflow, after Region 1/2 anchor-admission decisions, with section as the technical batch and mouse retained as the biological replicate."""),
        md("""## 1. Reproducible setup

HPC paths are defaults and can be overridden with `SCWAT_*` environment variables. No package is installed from this notebook. All outputs, checkpoints and temporary files must remain under the project root. Run from a fresh R kernel in order."""),
        code(f'''REGION_ID <- "{region}"
ANALYSIS_BRANCH <- "{branch}"
SECTION_ROLE <- "{section_role}"
RANDOM_SEED <- 1234L
dims_use <- 1:30
CLUSTER_ALGORITHM <- as.integer(Sys.getenv("SCWAT_CLUSTER_ALGORITHM", "1"))
if (!CLUSTER_ALGORITHM %in% 1:4) stop("SCWAT_CLUSTER_ALGORITHM must be one of Seurat algorithms 1-4.")
WRITE_CHECKPOINTS <- identical(toupper(Sys.getenv("SCWAT_WRITE_CHECKPOINTS", "TRUE")), "TRUE")

PROJECT_ROOT <- Sys.getenv(
  "SCWAT_PROJECT_ROOT",
  "/dssg/home/acct-svetoslav_chakarov/svetoslav_chakarov/Lab_members/Yanan_Hu/YNH_Xenium"
)
PIPELINE_REPO <- Sys.getenv(
  "SCWAT_PIPELINE_REPO",
  file.path(PROJECT_ROOT, "adipose_analysis_B2", "YNH_Xenium_scWAT")
)
INPUT_ROOT <- Sys.getenv("SCWAT_INPUT_ROOT", file.path(PROJECT_ROOT, "adipose_data_B2"))
QC_RUN_ROOT <- Sys.getenv(
  "SCWAT_QC_RUN_ROOT",
  file.path(PROJECT_ROOT, "adipose_analysis_B2", "scwat_qc_outputs", "scwat_qc_hpc", "sections", REGION_ID)
)
RUN_LABEL <- Sys.getenv("SCWAT_B2_RUN_LABEL", "full_panel_branch_v1")
DOWNSTREAM_ROOT <- Sys.getenv(
  "SCWAT_B2_OUTPUT_ROOT",
  file.path(PROJECT_ROOT, "adipose_analysis_B2", "scwat_downstream_outputs", RUN_LABEL)
)
REGION_SHARED_ROOT <- file.path(DOWNSTREAM_ROOT, REGION_ID, "00_shared_domain")
BRANCH_ROOT <- file.path(DOWNSTREAM_ROOT, REGION_ID, ANALYSIS_BRANCH)
FIGURE_ROOT <- file.path(BRANCH_ROOT, "figures")
TEMP_ROOT <- file.path(DOWNSTREAM_ROOT, "tmp", REGION_ID, ANALYSIS_BRANCH)
set.seed(RANDOM_SEED)'''),
        code('''source(file.path(PIPELINE_REPO, "R", "source.R"))
validate_runtime_paths(PROJECT_ROOT, INPUT_ROOT, DOWNSTREAM_ROOT, TEMP_ROOT)
assert_path_within(PROJECT_ROOT, REGION_SHARED_ROOT)
assert_path_within(PROJECT_ROOT, BRANCH_ROOT)
assert_path_within(PROJECT_ROOT, FIGURE_ROOT)
dir.create(REGION_SHARED_ROOT, recursive = TRUE, showWarnings = FALSE)
dir.create(BRANCH_ROOT, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURE_ROOT, recursive = TRUE, showWarnings = FALSE)
dir.create(TEMP_ROOT, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(TMPDIR = TEMP_ROOT, TMP = TEMP_ROOT, TEMP = TEMP_ROOT)

required_packages <- c(
  "Seurat", "SeuratObject", "Matrix", "dplyr", "tidyr", "tibble",
  "ggplot2", "FNN", "patchwork"
)
if (ANALYSIS_BRANCH == "all_QCpass") required_packages <- c(required_packages, "dbscan")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Missing HPC packages: ", paste(missing_packages, collapse = ", "),
       ". Install them in the controlled HPC environment; this notebook installs nothing.")
}
optional_packages <- c("diptest", "mclust", "CellChat")
optional_available <- setNames(
  vapply(optional_packages, requireNamespace, logical(1), quietly = TRUE),
  optional_packages
)
optional_versions <- data.frame(
  package = optional_packages,
  available = unname(optional_available),
  version = vapply(optional_packages, function(package) {
    if (optional_available[[package]]) as.character(packageVersion(package)) else NA_character_
  }, character(1)),
  stringsAsFactors = FALSE
)
RUN_DIPTEST <- unname(optional_available[["diptest"]])
RUN_MCLUST <- unname(optional_available[["mclust"]])
RUN_CELLCHAT <- unname(optional_available[["CellChat"]]) &&
  identical(toupper(Sys.getenv("SCWAT_RUN_CELLCHAT", "TRUE")), "TRUE")
RUN_WANG_INTEGRATION <- identical(toupper(Sys.getenv("SCWAT_RUN_WANG_INTEGRATION", "TRUE")), "TRUE")
CELLCHAT_STATE_PROPORTION <- 0.30
print(optional_versions)
if (CLUSTER_ALGORITHM == 4L && !requireNamespace("leidenbase", quietly = TRUE)) {
  stop("SCWAT_CLUSTER_ALGORITHM=4 requires package 'leidenbase'. Use algorithm 1 or install it in the controlled HPC environment.")
}
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
})
theme_set(cell_style_theme(base_size = 12))
sessionInfo()'''),
        md("""### 1.1 Inputs, versions and hashes

The file inventory is recorded before import. Hashing can be disabled only when the HPC filesystem makes full input hashing impractical; that decision is saved explicitly."""),
        code('''region_info <- discover_one_section(INPUT_ROOT, REGION_ID)
REGION_XENIUM_DIR <- region_info$region_dir[[1L]]
HASH_INPUTS <- identical(toupper(Sys.getenv("SCWAT_HASH_INPUTS", "TRUE")), "TRUE")
input_inventory <- inventory_section_files(
  REGION_XENIUM_DIR, REGION_ID, calculate_md5 = HASH_INPUTS
)
source_hash <- unname(tools::md5sum(file.path(PIPELINE_REPO, "R", "source.R")))
mask_file <- file.path(QC_RUN_ROOT, "cell_qc_metadata.tsv.gz")
eos_gene_set_path <- file.path(PIPELINE_REPO, "config", "eos_gene_sets.tsv")
stopifnot(file.exists(mask_file), file.exists(eos_gene_set_path))
provenance <- data.frame(
  region_id = REGION_ID, analysis_branch = ANALYSIS_BRANCH,
  section_role = SECTION_ROLE, run_label = RUN_LABEL, random_seed = RANDOM_SEED,
  cluster_algorithm = CLUSTER_ALGORITHM, write_checkpoints = WRITE_CHECKPOINTS,
  source_md5 = source_hash, input_hashing_enabled = HASH_INPUTS,
  timestamp_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  stringsAsFactors = FALSE
)
input_inventory
provenance'''),
        md("""## 2. Import and attach QC masks by cell ID

The complete Xenium object is imported with all 479 genes. QC fields are joined by explicit cell ID. `primary_include_revised` is reconstructed from the already aligned object metadata: core QC must pass, and high-control and segmentation-multiplet flags must both be absent."""),
        code('''region_spatial_all <- create_spatial_seurat_from_xenium(
  xenium_dir = REGION_XENIUM_DIR, project = REGION_ID, assay = "Xenium", fov = "fov",
  include_cell_segmentation = TRUE, include_nucleus_segmentation = TRUE
)
masks <- read.delim(mask_file, stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(anyDuplicated(masks$cell_id) == 0L)
region_spatial_all <- add_masks_to_seurat(
  region_spatial_all, masks, overwrite = FALSE, require_complete_match = TRUE
)
region_spatial_all$primary_include_revised <- derive_primary_include_revised(
  region_spatial_all@meta.data
)
stopifnot(identical(rownames(region_spatial_all@meta.data), colnames(region_spatial_all)))
table(primary_include_revised = region_spatial_all$primary_include_revised)'''),
        md("""## 3. QC counts and spatial masks before subsetting

These plots show what is included and excluded before any branch object is created. Review flags remain metadata and the raw object is retained."""),
        code('''qc_columns <- c(
  "qc_core_pass", "primary_include_revised", "high_control_flag",
  "segmentation_multiplet_flag", "nucleus_missing_flag",
  "multiple_nuclei_flag", "cell_area_outlier_flag"
)
qc_columns <- intersect(qc_columns, colnames(region_spatial_all@meta.data))
qc_counts <- do.call(rbind, lapply(qc_columns, function(field) data.frame(
  region_id = REGION_ID, field = field,
  n_flagged_or_included = sum(region_spatial_all@meta.data[[field]] %in% TRUE, na.rm = TRUE),
  pct = 100 * mean(region_spatial_all@meta.data[[field]] %in% TRUE, na.rm = TRUE)
)))
qc_counts'''),
        code('''options(repr.plot.width = 14, repr.plot.height = 6)
qc_coordinates <- GetTissueCoordinates(region_spatial_all[["fov"]], which = "centroids") %>% as.data.frame()
if (!"cell" %in% colnames(qc_coordinates)) qc_coordinates$cell <- rownames(qc_coordinates)
qc_coordinates <- qc_coordinates %>%
  rename(cell_id = cell) %>%
  mutate(
    qc_core_pass = as.character(region_spatial_all$qc_core_pass[match(cell_id, colnames(region_spatial_all))]),
    primary_include_revised = as.character(region_spatial_all$primary_include_revised[match(cell_id, colnames(region_spatial_all))])
  )
qc_core_plot <- plot_target_overlay(
  qc_coordinates, x_col = "x", y_col = "y", group_col = "qc_core_pass", target_labels = "FALSE",
  palette = c("TRUE" = "#D9D9D9", "FALSE" = "#963E49"),
  background_size = 0.10, target_size = 0.70, title = "Core-QC failures", fixed_coordinates = TRUE
)
primary_mask_plot <- plot_target_overlay(
  qc_coordinates, x_col = "x", y_col = "y", group_col = "primary_include_revised", target_labels = "FALSE",
  palette = c("TRUE" = "#D9D9D9", "FALSE" = "#963E49"),
  background_size = 0.10, target_size = 0.70, title = "Cells excluded from primary analysis", fixed_coordinates = TRUE
)
print(qc_core_plot)
print(primary_mask_plot)'''),
    ]

    if child_branch:
        cells.extend([
            md("""## 4. Load the frozen tissue-domain partition

The lymph-node boundary is not recalculated here. It must come from the successfully completed `all_QCpass` notebook for this region. The manifest is matched to object cells and its partition is validated before subsetting."""),
            code('''domain_manifest_path <- file.path(REGION_SHARED_ROOT, "lymph_node_domain_manifest.tsv.gz")
stopifnot(file.exists(domain_manifest_path))
domain_manifest <- read.delim(domain_manifest_path, stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(
  anyDuplicated(domain_manifest$cell_id) == 0L,
  setequal(domain_manifest$cell_id, colnames(region_spatial_all))
)
manifest_index <- match(colnames(region_spatial_all), domain_manifest$cell_id)
stopifnot(!anyNA(manifest_index))
domain_manifest <- domain_manifest[manifest_index, , drop = FALSE]
stopifnot(identical(domain_manifest$cell_id, colnames(region_spatial_all)))
stopifnot(identical(
  domain_manifest$primary_include_revised,
  as.logical(region_spatial_all$primary_include_revised)
))
table(domain_manifest$tissue_domain, useNA = "ifany")'''),
            code('''branch_key <- if (ANALYSIS_BRANCH == "adipose_only") "adipose_only" else "lymph_node_only"
minimum_branch_cells <- if (ANALYSIS_BRANCH == "lymph_node_only") 100L else 500L
branch_cell_ids <- select_tissue_branch_ids(
  domain_manifest, branch = branch_key, min_cells = minimum_branch_cells
)
region_spatial <- subset(region_spatial_all, cells = branch_cell_ids)
stopifnot(identical(colnames(region_spatial), branch_cell_ids))
region_spatial'''),
        ])
    else:
        cells.extend([
            md("""## 4. Create the all-QC-passed branch

Only `primary_include_revised` is used here. The lymph-node domain is derived later, after conservative cell annotation, and then frozen for the two child notebooks."""),
            code('''branch_cell_ids <- rownames(region_spatial_all@meta.data)[
  region_spatial_all$primary_include_revised %in% TRUE
]
stopifnot(length(branch_cell_ids) > 0L)
region_spatial <- subset(region_spatial_all, cells = branch_cell_ids)
stopifnot(all(region_spatial$primary_include_revised %in% TRUE))
region_spatial'''),
        ])

    cells.append(code('''checkpoint_records <- list()
if (WRITE_CHECKPOINTS) {
  checkpoint_records$branch_input <- write_validated_seurat_checkpoint(
    region_spatial,
    file.path(BRANCH_ROOT, "01_branch_input_479.rds"),
    PROJECT_ROOT, stage = "01_BRANCH_INPUT"
  )
  checkpoint_records$branch_input
}'''))

    cells.extend([
        md("""## 5. Complete-panel normalization and PCA

All 479 panel genes are retained. Counts are normalized with Seurat `LogNormalize`, every panel gene is scaled, and exactly PCs 1–30 are used downstream. Scaling may densify the 479-by-cell matrix, but the raw count layer remains sparse and unchanged."""),
        code('''raw_panel_genes <- rownames(region_spatial)
FIX_GENESET <- raw_panel_genes
stopifnot(length(FIX_GENESET) == 479L)
DefaultAssay(region_spatial) <- "Xenium"
region_spatial <- NormalizeData(
  region_spatial, assay = "Xenium", normalization.method = "LogNormalize",
  scale.factor = 10000, verbose = FALSE
)
region_spatial <- ScaleData(
  region_spatial, assay = "Xenium", features = FIX_GENESET, verbose = FALSE
)
region_spatial <- RunPCA(
  region_spatial, assay = "Xenium", features = FIX_GENESET,
  npcs = 30, seed.use = RANDOM_SEED, verbose = FALSE
)
stopifnot(ncol(Embeddings(region_spatial, "pca")) >= 30L)
region_spatial'''),
        code('''if (WRITE_CHECKPOINTS) {
  checkpoint_records$pca <- write_validated_seurat_checkpoint(
    region_spatial,
    file.path(BRANCH_ROOT, "02_log_normalized_scaled_pca_479.rds"),
    PROJECT_ROOT, stage = "02_LOGNORMALIZE_SCALE_PCA30"
  )
  checkpoint_records$pca
}'''),
        code('''options(repr.plot.width = 10, repr.plot.height = 5)
ElbowPlot(region_spatial, ndims = 30)
options(repr.plot.width = 16, repr.plot.height = 10)
VizDimLoadings(region_spatial, dims = 1:12, reduction = "pca", nfeatures = 15)
VizDimLoadings(region_spatial, dims = 13:30, reduction = "pca", nfeatures = 10)
options(repr.plot.width = 16, repr.plot.height = 12)
DimHeatmap(region_spatial, dims = 1:12, cells = min(1000L, ncol(region_spatial)), balanced = TRUE)
DimHeatmap(region_spatial, dims = 13:30, cells = min(1000L, ncol(region_spatial)), balanced = TRUE)'''),
        md("""### 5.1 Technical-dominance gate

Strong PC correlations with counts, detected genes or cell area are review evidence. They do not automatically justify regressing the covariate, because complexity can be biologically coupled to cell type. Instead, the correlations and cluster-level QC distributions are exported for later release gating."""),
        code('''pca_qc_correlations <- summarise_pca_qc_correlations(
  Embeddings(region_spatial, "pca"), region_spatial@meta.data,
  qc_metrics = c("nCount_Xenium", "nFeature_Xenium", "cell_area"),
  review_abs_rho = 0.5
)
pca_qc_correlations %>% arrange(desc(abs_rho)) %>% head(20)'''),
        md("""## 6. Neighbours, UMAP and quantitative clustering stability

The same PC1–30 graph specification is used for every branch. Four resolutions are run across three random seeds. Adjusted Rand index is calculated within each resolution before a working resolution is selected. Marker support and QC association remain required validation; ARI alone does not select biology."""),
        code('''region_spatial <- FindNeighbors(
  region_spatial, reduction = "pca", dims = dims_use, k.param = 20,
  graph.name = c("Xenium_nn_PC30", "Xenium_snn_PC30"), verbose = FALSE
)
region_spatial <- RunUMAP(
  region_spatial, reduction = "pca", dims = dims_use,
  reduction.name = "umap_PC30", seed.use = RANDOM_SEED, verbose = FALSE
)
resolutions <- c(0.7, 0.8, 1.0, 1.2)
cluster_seeds <- c(1234L, 2026L, 31415L)
cluster_column <- function(resolution, seed) paste0(
  "cluster_res_", gsub("\\\\.", "_", format(resolution, trim = TRUE)), "_seed_", seed
)
for (resolution in resolutions) {
  for (seed in cluster_seeds) {
    region_spatial <- FindClusters(
      region_spatial, graph.name = "Xenium_snn_PC30", algorithm = CLUSTER_ALGORITHM,
      resolution = resolution, cluster.name = cluster_column(resolution, seed),
      random.seed = seed, verbose = FALSE
    )
  }
}'''),
        code('''cluster_stability <- bind_rows(lapply(resolutions, function(resolution) {
  columns <- vapply(cluster_seeds, function(seed) cluster_column(resolution, seed), character(1))
  result <- summarise_cluster_stability(region_spatial@meta.data[, columns, drop = FALSE])
  transform(result$pairwise, resolution = resolution)
}))
cluster_stability_summary <- cluster_stability %>%
  group_by(resolution) %>%
  summarise(minimum_ari = min(ari), median_ari = median(ari), mean_ari = mean(ari), .groups = "drop")
cluster_stability_summary
style_cell_plot(ggplot(cluster_stability, aes(factor(resolution), ari, fill = factor(resolution))) +
  geom_boxplot(width = 0.55, outlier.shape = NA, show.legend = FALSE) + geom_jitter(width = 0.08, size = 2, colour = "#81506E") +
  scale_fill_manual(values = cell_macaron_palette(as.character(resolutions))) +
  coord_cartesian(ylim = c(-0.05, 1.02)) +
  labs(title = "Clustering stability across random seeds", x = "Resolution", y = "Pairwise adjusted Rand index"))'''),
        code('''WORKING_RESOLUTION <- as.numeric(Sys.getenv("SCWAT_WORKING_RESOLUTION", "0.8"))
stopifnot(WORKING_RESOLUTION %in% resolutions)
working_cluster_col <- cluster_column(WORKING_RESOLUTION, RANDOM_SEED)
region_spatial$working_cluster <- factor(region_spatial@meta.data[[working_cluster_col]])
region_spatial$cluster_res_0_8 <- region_spatial$working_cluster
cluster_sizes <- as.data.frame(table(cluster = region_spatial$working_cluster))
cluster_sizes
options(repr.plot.width = 12, repr.plot.height = 8)
style_cell_plot(DimPlot(
  region_spatial, reduction = "umap_PC30", group.by = "working_cluster",
  label = TRUE, repel = TRUE, cols = cell_macaron_palette(levels(region_spatial$working_cluster))
))'''),
        code('''cluster_qc_summary <- region_spatial@meta.data %>%
  mutate(working_cluster = as.character(working_cluster)) %>%
  group_by(working_cluster) %>%
  summarise(
    n_cells = n(), median_nCount = median(nCount_Xenium),
    median_nFeature = median(nFeature_Xenium),
    median_cell_area = if ("cell_area" %in% colnames(cur_data())) median(cell_area) else NA_real_,
    .groups = "drop"
  )
cluster_qc_summary'''),
        md("""## 7. Canonical-marker annotation

Cluster-level scores provide one evidence stream. Marker coverage is shown first because a targeted panel cannot support every fine label equally. Low margins are not silently converted into confident labels."""),
        code('''marker_df <- tibble::tribble(
  ~Gene_Symbol, ~CellType_main, ~CellType_subtype,
  "Gpihbp1", "Endothelial", "Capillary_EC", "Kdr", "Endothelial", "Capillary_EC",
  "Vwf", "Endothelial", "Venous_EC", "Plvap", "Endothelial", "Venous_EC",
  "Mmrn1", "Endothelial", "Lymphatic_EC", "Prox1", "Endothelial", "Lymphatic_EC", "Lyve1", "Endothelial", "Lymphatic_EC",
  "Pck1", "Adipocyte", "Adipocyte", "Retn", "Adipocyte", "Adipocyte", "Aqp7", "Adipocyte", "Adipocyte", "Car3", "Adipocyte", "Adipocyte",
  "Pi16", "Stromal_Fibroblast", "ASC", "Pcolce2", "Stromal_Fibroblast", "ASC", "Mfap5", "Stromal_Fibroblast", "ASC",
  "Lum", "Stromal_Fibroblast", "Fibroblast", "Serpinf1", "Stromal_Fibroblast", "Fibroblast", "Dpt", "Stromal_Fibroblast", "Fibroblast",
  "Tagln", "Mural", "VSMC", "Myh11", "Mural", "VSMC", "Cnn1", "Mural", "VSMC",
  "Rgs5", "Mural", "Pericyte", "Higd1b", "Mural", "Pericyte", "Ndufa4l2", "Mural", "Pericyte",
  "Plp1", "Neural", "Schwann", "Prx", "Neural", "Schwann", "Pmp22", "Neural", "Schwann",
  "Pf4", "Macrophage", "LYVE1_resident_Mac", "Folr2", "Macrophage", "LYVE1_resident_Mac", "F13a1", "Macrophage", "LYVE1_resident_Mac",
  "Marco", "Macrophage", "Scavenging_macrophage", "Mpeg1", "Macrophage", "Scavenging_macrophage", "Cd5l", "Macrophage", "Scavenging_macrophage",
  "Cebpb", "Monocyte", "Monocyte", "Ctss", "Monocyte", "Monocyte", "Cd300a", "Monocyte", "Monocyte",
  "Csf3r", "Neutrophil", "Neutrophil", "Cxcr2", "Neutrophil", "Neutrophil", "Mmp9", "Neutrophil", "Neutrophil",
  "Mrgpra2a", "Mast cell", "Mast cell", "Il1rl1", "Mast cell", "Mast cell", "Cpa3", "Mast cell", "Mast cell",
  "Alox15", "Eosinophil", "Eosinophil", "Siglecf", "Eosinophil", "Eosinophil", "Ccr3", "Eosinophil", "Eosinophil", "Il5ra", "Eosinophil", "Eosinophil", "Prg2", "Eosinophil", "Eosinophil",
  "Cd3d", "T cell", "T", "Cd8a", "T cell", "T", "Epsti1", "T cell", "T",
  "Itgal", "NK", "NK", "Klra8", "NK", "NK",
  "Basp1", "DC", "DC", "Cytip", "DC", "DC", "Cd24a", "DC", "cDC",
  "Krt19", "Mesothelial", "Mesothelial", "Cav1", "Mesothelial", "Mesothelial", "Igfbp6", "Mesothelial", "Mesothelial"
)
marker_df_use <- marker_df %>% filter(Gene_Symbol %in% rownames(region_spatial)) %>% distinct()
marker_coverage <- marker_df %>% mutate(present = Gene_Symbol %in% rownames(region_spatial)) %>%
  group_by(CellType_main, CellType_subtype) %>%
  summarise(n_markers = n(), n_present = sum(present), coverage = mean(present), .groups = "drop")
marker_coverage'''),
        code('''main_result <- score_marker_groups_by_cluster(
  region_spatial, marker_df_use, group_col = "CellType_main",
  cluster_col = "cluster_res_0_8", assay = "Xenium", layer = "data",
  score_prefix = "MainScore", detect_prefix = "MainDetect", label_prefix = "Main"
)
region_spatial <- main_result$object
subtype_result <- score_marker_groups_by_cluster(
  region_spatial, marker_df_use, group_col = "CellType_subtype",
  cluster_col = "cluster_res_0_8", assay = "Xenium", layer = "data",
  score_prefix = "SubtypeScore", detect_prefix = "SubtypeDetect", label_prefix = "Subtype"
)
region_spatial <- subtype_result$object

subtype_to_main <- marker_df_use %>%
  transmute(CellType_subtype = as.character(CellType_subtype), CellType_main = as.character(CellType_main), subtype_key = make.names(CellType_subtype)) %>% distinct()
xenium_cluster_annotation <- subtype_result$annotation %>%
  mutate(subtype_key = make.names(as.character(Subtype_label))) %>%
  rename(
    Xenium_cluster_subtype_score = Subtype_score,
    Xenium_cluster_subtype_second = Subtype_second,
    Xenium_cluster_subtype_second_score = Subtype_second_score,
    Xenium_cluster_subtype_margin = Subtype_margin
  ) %>%
  left_join(subtype_to_main, by = "subtype_key") %>%
  mutate(Xenium_cluster_subtype = CellType_subtype) %>%
  select(-subtype_key, -CellType_subtype, -Subtype_label) %>%
  rename(Xenium_cluster_main = CellType_main) %>%
  left_join(main_result$annotation %>% transmute(
    cluster_res_0_8, Xenium_main_score_label = as.character(Main_label),
    Xenium_main_score = Main_score, Xenium_main_second = as.character(Main_second),
    Xenium_main_second_score = Main_second_score, Xenium_main_margin = Main_margin
  ), by = "cluster_res_0_8")

xenium_columns <- setdiff(colnames(xenium_cluster_annotation), "cluster_res_0_8")
metadata_joined <- region_spatial@meta.data %>% rownames_to_column("cell_id") %>%
  select(-any_of(xenium_columns)) %>% mutate(cluster_res_0_8 = as.character(cluster_res_0_8)) %>%
  left_join(xenium_cluster_annotation %>% mutate(cluster_res_0_8 = as.character(cluster_res_0_8)), by = "cluster_res_0_8") %>%
  column_to_rownames("cell_id")
stopifnot(identical(rownames(metadata_joined), colnames(region_spatial)))
region_spatial@meta.data <- metadata_joined
subtype_result$annotation'''),
        md("""## 8. Wang-reference transfer and conservative reconciliation

The complete all-age Wang reference is primary to maximize reference Eosinophil representation for transfer and joint-embedding diagnostics. The 2.5-month transfer is retained as an age-matched sensitivity comparison. Both use all shared genes from the 479-gene panel. Fine labels with review-level evidence receive an explicit `Uncertain` analysis label rather than a forced assignment."""),
        code('''WANG_REFERENCE_PATH <- Sys.getenv(
  "SCWAT_WANG_REFERENCE",
  "/dssg/home/acct-svetoslav_chakarov/svetoslav_chakarov/Data/External_Data/Wang_Science2025_mouse_human/analysis_mouse/Wang_Science2025_cleaned_seurat.rds"
)
stopifnot(file.exists(WANG_REFERENCE_PATH))
scwat_ref <- readRDS(WANG_REFERENCE_PATH)
DefaultAssay(scwat_ref) <- "RNA"
scwat_ref$Wang_subtype_original <- as.character(scwat_ref$CellType_subtype)
wang_map <- c(
  "ASC"="ASC", "IAP"="APC", "CP-1"="APC", "CP-A"="APC", "CP-2"="APC", "RegFibro"="APC",
  "MatureAdip"="Adipocyte", "Pericyte"="Mural", "Arterial_EC"="Endothelial", "Capillary_EC"="Endothelial",
  "Epithelial"="Epithelial", "Mesothelial"="Mesothelial", "T"="T", "γδ T"="γδ T", "NK"="NK", "ILC2"="ILC2",
  "B"="B", "Plasma"="Plasma", "CCR2_inflammatory_Monocyte"="CCR2_inflammatory_Monocyte",
  "CX3CR1_Monocyte"="CX3CR1_Monocyte", "Inflammatory_Mac"="Inflammatory_Mac",
  "LYVE1_resident_Mac"="LYVE1_resident_Mac", "TREM2_LAM"="TREM2_LAM", "IFN_response_myeloid"="IFN_response_myeloid",
  "cDC1"="cDC1", "cDC2"="cDC2", "CCR7_migratory_DC"="CCR7_migratory_DC", "Neutrophil"="Neutrophil",
  "Eosinophil"="Eosinophil", "Mast cell"="Mast cell"
)
unexpected <- setdiff(unique(scwat_ref$Wang_subtype_original), names(wang_map))
if (length(unexpected)) stop("Unexpected Wang labels: ", paste(unexpected, collapse = ", "))
scwat_ref$Wang_subtype_harmonized <- unname(wang_map[scwat_ref$Wang_subtype_original])
scwat_ref$Wang_main <- case_when(
  scwat_ref$Wang_subtype_harmonized %in% c("ASC", "APC") ~ "Stromal_Fibroblast",
  scwat_ref$Wang_subtype_harmonized == "Mural" ~ "Mural",
  scwat_ref$Wang_subtype_harmonized == "Adipocyte" ~ "Adipocyte",
  scwat_ref$Wang_subtype_harmonized == "Endothelial" ~ "Endothelial",
  scwat_ref$Wang_subtype_harmonized %in% c("Epithelial", "Mesothelial") ~ scwat_ref$Wang_subtype_harmonized,
  scwat_ref$Wang_subtype_harmonized %in% c("T", "γδ T") ~ "T cell",
  scwat_ref$Wang_subtype_harmonized %in% c("NK", "ILC2", "B", "Plasma") ~ recode(scwat_ref$Wang_subtype_harmonized, ILC2="ILC", B="B cell"),
  scwat_ref$Wang_subtype_harmonized %in% c("CCR2_inflammatory_Monocyte", "CX3CR1_Monocyte") ~ "Monocyte",
  scwat_ref$Wang_subtype_harmonized %in% c("Inflammatory_Mac", "LYVE1_resident_Mac", "TREM2_LAM") ~ "Macrophage",
  scwat_ref$Wang_subtype_harmonized == "IFN_response_myeloid" ~ "Myeloid",
  scwat_ref$Wang_subtype_harmonized %in% c("cDC1", "cDC2", "CCR7_migratory_DC") ~ "DC",
  TRUE ~ scwat_ref$Wang_subtype_harmonized
)
stopifnot(!anyNA(scwat_ref$Wang_main))
ref_2p5 <- subset(scwat_ref, subset = Condition == "2.5 months")
ref_all <- scwat_ref
list(
  all_age_primary = table(ref_all$Wang_subtype_harmonized),
  age_matched_sensitivity = table(ref_2p5$Wang_subtype_harmonized)
)'''),
        code('''annotation_genes <- Reduce(intersect, list(FIX_GENESET, rownames(ref_all), rownames(region_spatial)))
stopifnot(length(annotation_genes) >= 100L)
transfer_2p5 <- run_wang_transfer(ref_2p5, region_spatial, annotation_genes, prefix = "Ref2p5")
transfer_all <- run_wang_transfer(ref_all, region_spatial, annotation_genes, prefix = "RefAll")
region_spatial <- AddMetaData(region_spatial, transfer_2p5$main)
region_spatial <- AddMetaData(region_spatial, transfer_2p5$subtype)
region_spatial <- AddMetaData(region_spatial, transfer_all$main)
region_spatial <- AddMetaData(region_spatial, transfer_all$subtype)
annotation_result <- refine_xenium_celltypes(
  region_spatial, reduction = "pca", dims = dims_use, k = 15,
  self_weight = 0.70, wang_prefix = "RefAll"
)
region_spatial <- annotation_result$object
region_spatial$Final_CellType_subtype_with_uncertain <- ifelse(
  region_spatial$Final_annotation_confidence == "REVIEW" |
    region_spatial$Final_review_flag %in% TRUE,
  "Uncertain", as.character(region_spatial$Final_CellType_subtype)
)
# Preserve the uncertainty-aware field for review, but use the evidence-based
# definitive subtype field for downstream summaries and biological ordering.
region_spatial$Final_CellType_subtype <- apply_scwat_cell_type_order(
  region_spatial$Final_CellType_subtype
)
region_spatial$Final_CellType_subtype_with_uncertain <- apply_scwat_cell_type_order(
  region_spatial$Final_CellType_subtype_with_uncertain
)
table(region_spatial$Final_annotation_confidence, useNA = "ifany")
table(region_spatial$Final_CellType_subtype_with_uncertain, useNA = "ifany")'''),
        code('''if (WRITE_CHECKPOINTS) {
  checkpoint_records$annotation <- write_validated_seurat_checkpoint(
    region_spatial,
    file.path(BRANCH_ROOT, "03_annotated_479.rds"),
    PROJECT_ROOT, stage = "03_MARKER_WANG_ANNOTATION"
  )
  checkpoint_records$annotation
}'''),
        code('''reference_concordance <- data.frame(
  age_matched = as.character(region_spatial$Ref2p5_subtype_predicted.id),
  all_age = as.character(region_spatial$RefAll_subtype_predicted.id)
)
reference_concordance_summary <- reference_concordance %>%
  count(age_matched, all_age, name = "n_cells") %>% arrange(desc(n_cells))
reference_concordance_summary %>% head(30)
options(repr.plot.width = 15, repr.plot.height = 9)
annotation_palette <- cell_macaron_palette(levels(region_spatial$Final_CellType_subtype))
annotation_palette["Eosinophil"] <- "#8F3340"
annotation_umap <- as.data.frame(Embeddings(region_spatial, "umap_PC30")) %>%
  rownames_to_column("cell_id") %>%
  mutate(cell_type = as.character(region_spatial$Final_CellType_subtype[match(cell_id, colnames(region_spatial))]))
annotation_spatial <- GetTissueCoordinates(region_spatial[["fov"]], which = "centroids") %>% as.data.frame()
if (!"cell" %in% colnames(annotation_spatial)) annotation_spatial$cell <- rownames(annotation_spatial)
annotation_spatial <- annotation_spatial %>% rename(cell_id = cell) %>%
  mutate(cell_type = as.character(region_spatial$Final_CellType_subtype[match(cell_id, colnames(region_spatial))]))
annotation_umap_centres <- annotation_umap %>%
  group_by(cell_type) %>%
  summarise(label_x = median(.data[[colnames(annotation_umap)[2]]]),
            label_y = median(.data[[colnames(annotation_umap)[3]]]), .groups = "drop")
annotation_umap_plot <- plot_target_overlay(
  annotation_umap, x_col = colnames(annotation_umap)[2], y_col = colnames(annotation_umap)[3],
  group_col = "cell_type", target_labels = "Eosinophil", palette = annotation_palette,
  background_size = 0.22, target_size = 1.15, title = "Definitive subtype labels used downstream"
) + geom_text(data = annotation_umap_centres, aes(label_x, label_y, label = cell_type),
              inherit.aes = FALSE, size = 3, colour = "#3F3F3F", check_overlap = TRUE)
print(annotation_umap_plot)
print(plot_target_overlay(
  annotation_spatial, x_col = "x", y_col = "y", group_col = "cell_type", target_labels = "Eosinophil",
  palette = annotation_palette, background_size = 0.10, target_size = 0.85,
  title = "Spatial definitive subtype labels", fixed_coordinates = TRUE
))
marker_plot_features <- order_marker_features(
  marker_df_use, available_genes = rownames(region_spatial)
)$feature_groups
print(style_cell_plot(DotPlot(
  region_spatial, features = marker_plot_features,
  group.by = "Final_CellType_subtype", assay = "Xenium",
  dot.scale = 5, scale = TRUE, cols = c("#D9D9D9", "#5A2F5E")
)) + RotatedAxis() + labs(
  title = "Canonical-marker consistency (internal, not independent validation)",
  x = NULL, y = NULL
))'''),
        md("""## 9. Eosinophil identity: annotation plus Tier 1/2 sensitivity set

The inclusive set is the union of existing Eosinophil annotation and Tier 1/2 evidence, matching the prior analysis. Crucially, rescued cells do not overwrite the principal cell-type label. Origin, evidence tier and confidence remain separate fields."""),
        code('''eos_result <- score_eosinophil_likeness(
  reference = ref_all, query = region_spatial,
  reference_group_col = "Wang_subtype_harmonized",
  wang_predicted_col = "RefAll_subtype_predicted.id",
  wang_eos_score_col = "RefAll_subtype_prediction.score.Eosinophil"
)
region_spatial <- eos_result$object
existing_eos <- as.character(region_spatial$Final_CellType_subtype) == "Eosinophil"
tier12 <- region_spatial$EosRef_call %in% c("REF_EOS_TIER1", "REF_EOS_TIER2")
region_spatial$Eos_inclusive <- existing_eos | tier12
region_spatial$Eos_primary_confident <- existing_eos & tier12 &
  region_spatial$Final_annotation_confidence != "REVIEW"
region_spatial$Eos_origin <- case_when(
  existing_eos & tier12 ~ "ANNOTATION_AND_TIER12",
  existing_eos & !tier12 ~ "ANNOTATION_ONLY_REVIEW",
  !existing_eos & tier12 ~ "TIER12_RESCUE_REVIEW",
  TRUE ~ "NOT_EOS"
)
region_spatial$Eos_manual_review <- region_spatial$Eos_origin %in% c("ANNOTATION_ONLY_REVIEW", "TIER12_RESCUE_REVIEW")
eos_identity_summary <- region_spatial@meta.data %>% count(Eos_origin, EosRef_tier, EosRef_call, name = "n_cells")
eos_identity_summary'''),
        code('''eos_plot <- plot_eos_call_by_subtype(
  region_spatial, subtype_col = "Final_CellType_subtype", return_data = TRUE
)
options(repr.plot.width = 12, repr.plot.height = 7)
eos_plot$plot
options(repr.plot.width = 14, repr.plot.height = 7)
eos_origin_coordinates <- GetTissueCoordinates(region_spatial[["fov"]], which = "centroids") %>% as.data.frame()
if (!"cell" %in% colnames(eos_origin_coordinates)) eos_origin_coordinates$cell <- rownames(eos_origin_coordinates)
eos_origin_coordinates <- eos_origin_coordinates %>% rename(cell_id = cell) %>%
  mutate(Eos_origin = region_spatial$Eos_origin[match(cell_id, colnames(region_spatial))])
plot_target_overlay(
  eos_origin_coordinates, x_col = "x", y_col = "y", group_col = "Eos_origin",
  target_labels = c("ANNOTATION_AND_TIER12", "ANNOTATION_ONLY_REVIEW", "TIER12_RESCUE_REVIEW"),
  palette = c("ANNOTATION_AND_TIER12" = "#9B3C48", "ANNOTATION_ONLY_REVIEW" = "#B66A25",
              "TIER12_RESCUE_REVIEW" = "#694F85", "NOT_EOS" = "#D9D9D9"),
  background_size = 0.10, target_size = 0.95,
  title = "Spatial Eosinophil evidence origin", fixed_coordinates = TRUE
)'''),
        md("""### 9.1 Wang–Xenium joint integration diagnostic

This is deliberately separate from label transfer. A deterministic, subtype-balanced sample from the complete all-age Wang reference is jointly integrated with the Xenium section by Seurat CCA using only genes shared with the complete 479-gene panel. The larger reference Eosinophil pool improves the diagnostic's ability to reveal an Eosinophil-enriched joint cluster. The integrated UMAP and cross-dataset Eosinophil-neighbour summaries assess alignment; they do not replace the Xenium PCA, clusters or final labels."""),
        code('''wang_reference_balanced <- sample_wang_reference(
  ref_all, subtype_col = "Wang_subtype_harmonized",
  max_per_subtype = as.integer(Sys.getenv("SCWAT_WANG_MAX_PER_SUBTYPE", "1000")),
  seed = RANDOM_SEED
)
wang_reference_balanced$integration_cell_type <- as.character(wang_reference_balanced$Wang_subtype_harmonized)
wang_reference_balanced$integration_is_eosinophil <- wang_reference_balanced$Wang_subtype_harmonized == "Eosinophil"
region_spatial$integration_cell_type <- as.character(region_spatial$Final_CellType_subtype)
region_spatial$integration_is_eosinophil <- region_spatial$Eos_inclusive %in% TRUE
wang_shared_gene_gate <- validate_shared_feature_set(
  rownames(wang_reference_balanced[["RNA"]]), rownames(region_spatial[["Xenium"]]),
  FIX_GENESET, min_shared = 100L
)
wang_shared_gene_gate[c("status", "message")]
table(wang_reference_balanced$integration_cell_type)'''),
        code('''wang_xenium_integration <- if (RUN_WANG_INTEGRATION && wang_shared_gene_gate$status == "PASS") {
  run_wang_xenium_integration(
    wang_reference_balanced, region_spatial,
    features = wang_shared_gene_gate$features, dims = dims_use,
    seed = RANDOM_SEED, reference_assay = "RNA", query_assay = "Xenium",
    min_shared = 100L, resolution = 0.8
  )
} else {
  list(
    status = if (RUN_WANG_INTEGRATION) wang_shared_gene_gate$status else "SKIPPED_BY_CONFIGURATION",
    message = if (RUN_WANG_INTEGRATION) wang_shared_gene_gate$message else "SCWAT_RUN_WANG_INTEGRATION is FALSE.",
    stage = "CONFIGURATION_GATE", object = NULL, anchors = NULL,
    parameters = data.frame()
  )
}
wang_xenium_integration_status <- data.frame(
  status = wang_xenium_integration$status,
  stage = wang_xenium_integration$stage,
  message = wang_xenium_integration$message,
  stringsAsFactors = FALSE
)
wang_xenium_neighbour_summary <- list(per_cell = data.frame(), dataset_summary = data.frame(), cluster_enrichment = data.frame())
wang_xenium_integration_status'''),
        code('''if (identical(wang_xenium_integration$status, "PASS")) {
  wang_xenium_object <- wang_xenium_integration$object
  integration_metadata <- wang_xenium_object@meta.data
  integration_metadata$dataset <- integration_metadata$integration_dataset
  integration_metadata$is_eosinophil <- integration_metadata$integration_is_eosinophil %in% TRUE
  integration_metadata$integrated_cluster <- as.character(integration_metadata$wang_xenium_cluster)
  integration_embeddings <- Embeddings(wang_xenium_object, "pca")[, dims_use, drop = FALSE]
  wang_xenium_neighbour_summary <- summarise_cross_dataset_eos_neighbours(
    integration_embeddings, integration_metadata,
    dataset_col = "dataset", eos_col = "is_eosinophil",
    cluster_col = "integrated_cluster", k = 15L
  )
  wang_xenium_object$integration_eos_display <- ifelse(
    wang_xenium_object$integration_is_eosinophil %in% TRUE,
    paste(wang_xenium_object$integration_dataset, "Eosinophil"), "Other"
  )
  wang_xenium_object$integration_cell_type <- apply_scwat_cell_type_order(
    wang_xenium_object$integration_cell_type
  )
  wang_umap_data <- as.data.frame(Embeddings(wang_xenium_object, "wang_xenium_umap")) %>%
    rownames_to_column("cell_id") %>%
    mutate(
      integration_dataset = wang_xenium_object$integration_dataset[match(cell_id, colnames(wang_xenium_object))],
      integration_cell_type = as.character(wang_xenium_object$integration_cell_type[match(cell_id, colnames(wang_xenium_object))]),
      integration_eos_display = wang_xenium_object$integration_eos_display[match(cell_id, colnames(wang_xenium_object))]
    )
  options(repr.plot.width = 14, repr.plot.height = 6)
  wang_dataset_plot <- plot_target_overlay(
    wang_umap_data, x_col = colnames(wang_umap_data)[2], y_col = colnames(wang_umap_data)[3],
    group_col = "integration_dataset", target_labels = "XENIUM",
    palette = c(WANG = "#91C9B6", XENIUM = "#B84E4B"),
    background_size = 0.20, target_size = 0.65,
    title = "Wang–Xenium CCA integration: dataset mixing"
  )
  wang_subtype_palette <- cell_macaron_palette(unique(wang_umap_data$integration_cell_type))
  wang_subtype_palette["Eosinophil"] <- "#8F3340"
  wang_subtype_plot <- plot_target_overlay(
    wang_umap_data, x_col = colnames(wang_umap_data)[2], y_col = colnames(wang_umap_data)[3],
    group_col = "integration_cell_type", target_labels = "Eosinophil",
    palette = wang_subtype_palette, background_size = 0.20, target_size = 0.90,
    title = "Wang and Xenium subtype structure in joint space"
  )
  wang_eos_plot <- plot_target_overlay(
    wang_umap_data, x_col = colnames(wang_umap_data)[2], y_col = colnames(wang_umap_data)[3],
    group_col = "integration_eos_display", target_labels = c("WANG Eosinophil", "XENIUM Eosinophil"),
    palette = c("Other" = "#D9D9D9", "WANG Eosinophil" = "#3979A8", "XENIUM Eosinophil" = "#B84E4B"),
    background_size = 0.18, target_size = 1.00,
    title = "Wang and Xenium Eosinophils in joint integrated space"
  )
  wang_cross_distance_plot <- style_cell_plot(ggplot(
    wang_xenium_neighbour_summary$per_cell,
    aes(dataset, minimum_cross_dataset_distance, fill = dataset)
  ) + geom_boxplot(outlier.shape = NA) +
    scale_fill_manual(values = c(WANG = "#91C9B6", XENIUM = "#F29B8F")) +
    facet_wrap(~is_eosinophil, labeller = label_both) +
    labs(title = "Cross-dataset neighbour distance", x = NULL, y = "Minimum PCA distance", fill = NULL))
  print(wang_dataset_plot)
  print(wang_subtype_plot)
  print(wang_eos_plot)
  print(wang_cross_distance_plot)
  save_cell_plot(wang_dataset_plot, "wang_xenium_dataset_umap", FIGURE_ROOT, PROJECT_ROOT, width = 10, height = 7)
  save_cell_plot(wang_subtype_plot, "wang_xenium_subtype_umap", FIGURE_ROOT, PROJECT_ROOT, width = 13, height = 8)
  save_cell_plot(wang_eos_plot, "wang_xenium_eosinophil_umap", FIGURE_ROOT, PROJECT_ROOT, width = 10, height = 7)
  save_cell_plot(wang_cross_distance_plot, "wang_xenium_cross_dataset_distance", FIGURE_ROOT, PROJECT_ROOT, width = 9, height = 6)
  if (WRITE_CHECKPOINTS) {
    checkpoint_records$wang_integration <- write_validated_seurat_checkpoint(
      wang_xenium_object,
      file.path(BRANCH_ROOT, "03b_wang_xenium_integration.rds"),
      PROJECT_ROOT, stage = "03B_WANG_XENIUM_CCA_DIAGNOSTIC"
    )
  }
  print(wang_xenium_neighbour_summary$dataset_summary)
  print(wang_xenium_neighbour_summary$cluster_enrichment %>%
          arrange(desc(eosinophil_fraction)) %>% head(20))
} else {
  message(wang_xenium_integration$message)
}'''),
        md("""## 10. Eosinophil-state continuum

State is calculated only when at least 20 inclusive Eosinophil candidates are available. The continuous balance is primary. Dip and mixture-model outputs are diagnostics only: a mixture fit does not establish discrete biological subtypes. Quantile tails are descriptive and are not independent discovery groups."""),
        code('''eos_gene_sets <- read.delim(eos_gene_set_path, stringsAsFactors = FALSE)
eos_cells <- colnames(region_spatial)[region_spatial$Eos_inclusive %in% TRUE]
RUN_EOS_STATE <- length(eos_cells) >= 20L
cat("Inclusive Eos cells:", length(eos_cells), "\n")
if (!RUN_EOS_STATE) message("Eos state analysis skipped: fewer than 20 inclusive Eos candidates.")
region_spatial$EosState_balance <- NA_real_
region_spatial$EosState_extreme <- NA_character_
eos_obj <- if (RUN_EOS_STATE) subset(region_spatial, cells = eos_cells) else NULL'''),
        code('''mixture_diagnostic <- list(
  status = "SKIPPED_NO_EOS_STATE", message = "Eosinophil-state analysis did not meet its cell-count gate.",
  package_version = if (RUN_MCLUST) as.character(packageVersion("mclust")) else NA_character_,
  selected_G = NA_integer_, model_name = NA_character_, fit = NULL, bic_table = data.frame()
)
if (RUN_EOS_STATE) {
  short_genes <- eos_gene_sets$gene[eos_gene_sets$gene_set == "short_lived"]
  long_genes <- eos_gene_sets$gene[eos_gene_sets$gene_set == "long_lived"]
  short_use <- intersect(short_genes, rownames(eos_obj))
  long_use <- intersect(long_genes, rownames(eos_obj))
  short_nonrib <- setdiff(short_use, grep("^Rp[ls]", short_use, value = TRUE))
  state_genes <- unique(c(short_use, long_use))
  stopifnot(length(short_nonrib) > 0L, length(long_use) > 0L)
  eos_obj <- ScaleData(eos_obj, assay = "Xenium", features = state_genes, verbose = FALSE)
  state_scaled <- GetAssayData(eos_obj, assay = "Xenium", layer = "scale.data")[state_genes, , drop = FALSE]
  eos_obj$EosShort_score_full <- colMeans(state_scaled[short_use, , drop = FALSE])
  eos_obj$EosShort_score_nonrib <- colMeans(state_scaled[short_nonrib, , drop = FALSE])
  eos_obj$EosLong_score <- colMeans(state_scaled[long_use, , drop = FALSE])
  eos_obj$EosState_score <- eos_obj$EosLong_score - eos_obj$EosShort_score_nonrib
  eos_obj$EosState_score_full <- eos_obj$EosLong_score - eos_obj$EosShort_score_full
  eos_obj$EosShort_z <- as.numeric(scale(eos_obj$EosShort_score_nonrib))
  eos_obj$EosLong_z <- as.numeric(scale(eos_obj$EosLong_score))
  eos_obj$EosState_balance <- (eos_obj$EosLong_z - eos_obj$EosShort_z) / sqrt(2)
  eos_obj$EosState_activity <- (eos_obj$EosLong_z + eos_obj$EosShort_z) / sqrt(2)
  state_counts <- FetchData(eos_obj, vars = state_genes, layer = "counts")
  eos_obj$EosState_n_detected <- rowSums(state_counts > 0)
  q_low <- quantile(eos_obj$EosState_balance, 0.10, na.rm = TRUE)
  q_high <- quantile(eos_obj$EosState_balance, 0.90, na.rm = TRUE)
  eos_obj$EosState_extreme <- case_when(
    eos_obj$EosState_balance >= q_high & eos_obj$EosLong_z > 0 ~ "Long-lived-like",
    eos_obj$EosState_balance <= q_low & eos_obj$EosShort_z > 0 ~ "Short-lived-like",
    TRUE ~ "Intermediate"
  )
  eos_obj$EosState_extreme <- factor(eos_obj$EosState_extreme, levels = c("Short-lived-like", "Intermediate", "Long-lived-like"))
}
if (RUN_EOS_STATE) table(eos_obj$EosState_extreme)'''),
        code('''if (RUN_EOS_STATE) {
  eos_state_qc <- data.frame(
    metric = c("balance_vs_nCount", "balance_vs_nFeature", "balance_vs_detected", "activity_vs_detected", "full_vs_nonrib_balance"),
    spearman_rho = c(
      cor(eos_obj$EosState_balance, eos_obj$nCount_Xenium, method = "spearman"),
      cor(eos_obj$EosState_balance, eos_obj$nFeature_Xenium, method = "spearman"),
      cor(eos_obj$EosState_balance, eos_obj$EosState_n_detected, method = "spearman"),
      cor(eos_obj$EosState_activity, eos_obj$EosState_n_detected, method = "spearman"),
      cor(eos_obj$EosState_score_full, eos_obj$EosState_score, method = "spearman")
    )
  )
  print(eos_state_qc)
  if (RUN_DIPTEST) {
    dip_result <- diptest::dip.test(eos_obj$EosState_balance)
    print(dip_result)
  } else {
    message("Dip-test diagnostic skipped: optional package 'diptest' is unavailable.")
  }
  mixture_diagnostic <- run_mclust_diagnostic(
    eos_obj$EosState_balance, G = 1:3, seed = RANDOM_SEED
  )
  print(mixture_diagnostic[c("status", "message", "selected_G", "model_name", "bic_table")])
  if (identical(mixture_diagnostic$status, "PASS")) {
    mixture_plot_data <- data.frame(
      EosState_balance = eos_obj$EosState_balance[is.finite(eos_obj$EosState_balance)],
      component = factor(mixture_diagnostic$fit$classification)
    )
    mixture_plot <- ggplot(mixture_plot_data, aes(EosState_balance, fill = component)) +
      geom_density(alpha = 0.45, linewidth = 0.4) +
      scale_fill_manual(values = cell_macaron_palette(levels(mixture_plot_data$component))) +
      labs(title = "Gaussian-mixture diagnostic on continuous Eosinophil state",
           x = "EosState balance", y = "Density", fill = "Component")
    mixture_plot <- style_cell_plot(mixture_plot)
    print(mixture_plot)
    save_cell_plot(mixture_plot, "eos_state_mclust_diagnostic", FIGURE_ROOT, PROJECT_ROOT, width = 9, height = 5)
  }
}
mclust_status <- data.frame(
  status = mixture_diagnostic$status, message = mixture_diagnostic$message,
  package_version = mixture_diagnostic$package_version,
  selected_G = mixture_diagnostic$selected_G, model_name = mixture_diagnostic$model_name,
  interpretation = "DIAGNOSTIC_NOT_BIOLOGICAL_SUBTYPE_ASSIGNMENT",
  stringsAsFactors = FALSE
)
mclust_status'''),
        code('''if (RUN_EOS_STATE) {
  options(repr.plot.width = 13, repr.plot.height = 5)
  eos_state_scatter_data <- eos_obj@meta.data %>% rownames_to_column("cell_id")
  eos_state_intermediate <- eos_state_scatter_data %>% filter(EosState_extreme == "Intermediate")
  eos_state_tails <- eos_state_scatter_data %>%
    filter(EosState_extreme %in% c("Short-lived-like", "Long-lived-like"))
  state_range <- range(eos_state_scatter_data$EosState_balance, na.rm = TRUE)
  state_binwidth <- diff(state_range) / 30
  if (!is.finite(state_binwidth) || state_binwidth <= 0) state_binwidth <- 1
  eos_state_histogram <- ggplot() +
    geom_histogram(
      data = eos_state_intermediate, aes(EosState_balance, fill = EosState_extreme),
      binwidth = state_binwidth, boundary = state_range[[1L]], colour = "white", linewidth = 0.2
    ) +
    geom_histogram(
      data = eos_state_tails, aes(EosState_balance, fill = EosState_extreme),
      binwidth = state_binwidth, boundary = state_range[[1L]], colour = "white", linewidth = 0.2
    ) +
    geom_vline(xintercept = 0, linetype = 2) +
    scale_fill_manual(values = c("Short-lived-like" = "#3979A8", "Intermediate" = "#CFCFCF", "Long-lived-like" = "#B84E4B")) +
    labs(title = "Continuous Eosinophil-state balance", x = "Long-lived-like − short-lived-like", y = "Cells", fill = NULL)
  eos_state_scatter <- ggplot() +
    geom_point(data = eos_state_intermediate, aes(EosShort_z, EosLong_z),
               colour = "#CFCFCF", size = 0.75, alpha = 0.55) +
    geom_point(data = eos_state_tails, aes(EosShort_z, EosLong_z, colour = EosState_extreme),
               size = 1.65, alpha = 0.95) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey55") +
    scale_colour_manual(values = c("Short-lived-like" = "#3979A8", "Long-lived-like" = "#B84E4B")) +
    labs(title = "Eosinophil state is displayed as a continuum", x = "Short-lived signature z-score",
         y = "Long-lived signature z-score", colour = NULL)
  print(style_cell_plot(eos_state_histogram))
  print(style_cell_plot(eos_state_scatter))
  plot_eos_state_heatmap(eos_obj, eos_gene_sets)
}'''),
        md("""### 10.A Exploratory within-section state-associated markers

This is not mouse-level differential expression. The groups are cell-level descriptive tails defined by the state signatures. State genes are excluded from testing to reduce circularity, and results are labelled exploratory."""),
        code('''eos_state_association <- data.frame()
if (RUN_EOS_STATE) {
  tail_counts <- table(eos_obj$EosState_extreme)
  if (all(tail_counts[c("Short-lived-like", "Long-lived-like")] >= 20L)) {
    Idents(eos_obj) <- "EosState_extreme"
    association_features <- setdiff(FIX_GENESET, state_genes)
    eos_state_association <- FindMarkers(
      eos_obj, ident.1 = "Long-lived-like", ident.2 = "Short-lived-like",
      features = association_features, assay = "Xenium", test.use = "wilcox",
      logfc.threshold = 0, min.pct = 0
    ) %>% rownames_to_column("gene") %>%
      mutate(analysis_scope = "EXPLORATORY_WITHIN_SECTION_ASSOCIATION",
             biological_replicate_inference = FALSE)
  } else {
    message("Exploratory tail association skipped: fewer than 20 cells in one or both tails.")
  }
}
head(eos_state_association, 30)'''),
        code('''if (RUN_EOS_STATE) {
  state_fields <- c(
    "EosShort_score_full", "EosShort_score_nonrib", "EosLong_score",
    "EosState_score", "EosState_score_full", "EosShort_z", "EosLong_z",
    "EosState_balance", "EosState_activity", "EosState_n_detected", "EosState_extreme"
  )
  eos_metadata <- eos_obj@meta.data
  for (field in state_fields) {
    region_spatial@meta.data[[field]] <- eos_metadata[match(colnames(region_spatial), rownames(eos_metadata)), field]
  }
}'''),
        code('''if (RUN_EOS_STATE) {
  state_coordinates <- GetTissueCoordinates(region_spatial[["fov"]], which = "centroids") %>% as.data.frame()
  if (!"cell" %in% colnames(state_coordinates)) state_coordinates$cell <- rownames(state_coordinates)
  state_coordinates$EosState_balance <- region_spatial$EosState_balance[match(state_coordinates$cell, colnames(region_spatial))]
  state_coordinates$EosState_extreme <- region_spatial$EosState_extreme[match(state_coordinates$cell, colnames(region_spatial))]
  eos_state_coordinates <- state_coordinates %>% filter(!is.na(EosState_balance))
  balance_limit <- max(abs(eos_state_coordinates$EosState_balance), na.rm = TRUE)
  options(repr.plot.width = 14, repr.plot.height = 6)
  print(
    ggplot() +
      geom_point(data = state_coordinates, aes(x, y), colour = "#E8E8E8", size = 0.10, alpha = 0.45) +
      geom_point(data = eos_state_coordinates %>% filter(EosState_extreme == "Intermediate"),
                 aes(x, y), colour = "#CFCFCF", size = 0.45, alpha = 0.60) +
      geom_point(data = eos_state_coordinates %>% filter(EosState_extreme != "Intermediate"),
                 aes(x, y, colour = EosState_balance), size = 1.5, alpha = 0.95) +
      scale_colour_gradient2(low = "#3979A8", mid = "#CFCFCF", high = "#B84E4B", midpoint = 0,
                             limits = c(-balance_limit, balance_limit)) +
      coord_fixed() + theme_void() +
      labs(title = "Spatial Eosinophil-state continuum", colour = "State balance")
  )
  print(
    ggplot() +
      geom_point(data = state_coordinates, aes(x, y), colour = "#E8E8E8", size = 0.10, alpha = 0.45) +
      geom_point(data = eos_state_coordinates %>% filter(EosState_extreme == "Intermediate"),
        aes(x, y), colour = "#CFCFCF", size = 0.45, alpha = 0.60) +
      geom_point(
        data = eos_state_coordinates %>% filter(EosState_extreme != "Intermediate"),
        aes(x, y, colour = EosState_extreme), size = 1.5, alpha = 0.95
      ) +
      scale_colour_manual(values = c("Short-lived-like" = "#3979A8", "Long-lived-like" = "#B84E4B")) +
      coord_fixed() + theme_void() + labs(title = "Descriptive Eosinophil-state tails", colour = "Tail")
  )
}'''),
        md("""## 10.1 Build disjoint Eosinophil and non-Eosinophil spatial pools

Eosinophil query IDs and non-Eosinophil reference IDs are constructed from the same explicit `Eos_inclusive` field. The downstream subtype is `Final_CellType_subtype`; the uncertainty-aware field remains available only for review. Duplicate IDs, pool overlap and duplicate cross-pool coordinates are hard stops."""),
        code('''coordinates <- GetTissueCoordinates(region_spatial[["fov"]], which = "centroids") %>% as.data.frame()
if (!"cell" %in% colnames(coordinates)) coordinates$cell <- rownames(coordinates)
stopifnot(all(c("x", "y", "cell") %in% colnames(coordinates)))
coordinates <- coordinates %>%
  mutate(
    Eos_inclusive = region_spatial$Eos_inclusive[match(cell, colnames(region_spatial))],
    cell_type = as.character(region_spatial$Final_CellType_subtype[match(cell, colnames(region_spatial))]),
    EosState_balance = region_spatial$EosState_balance[match(cell, colnames(region_spatial))]
  )
spatial_metadata <- coordinates %>%
  transmute(cell_id = cell, Eos_inclusive, cell_type, EosState_balance,
            EosState_extreme = region_spatial$EosState_extreme[match(cell, colnames(region_spatial))])
spatial_pools <- build_eos_spatial_pools(
  spatial_metadata, coordinates %>% transmute(cell_id = cell, x, y),
  eos_col = "Eos_inclusive", cell_type_col = "cell_type",
  state_col = "EosState_balance", extreme_col = "EosState_extreme"
)
RUN_EOS_SPATIAL <- identical(spatial_pools$status, "PASS")
spatial_pool_gate <- spatial_pools$gate
print(spatial_pool_gate)
if (RUN_EOS_SPATIAL) {
  spatial_pool_plot_data <- spatial_pools$all %>% mutate(
    display_group = case_when(
      EosState_extreme == "Short-lived-like" ~ "Short-lived-like",
      EosState_extreme == "Long-lived-like" ~ "Long-lived-like",
      TRUE ~ "Reference / intermediate"
    )
  )
  spatial_pool_plot <- plot_target_overlay(
    spatial_pool_plot_data, x_col = "x", y_col = "y", group_col = "display_group",
    target_labels = c("Short-lived-like", "Long-lived-like"),
    palette = c("Reference / intermediate" = "#CFCFCF", "Short-lived-like" = "#3979A8", "Long-lived-like" = "#B84E4B"),
    background_size = 0.12, target_size = 1.45,
    title = "Step 10.1: disjoint pools with Eosinophil-state tails", fixed_coordinates = TRUE
  )
  print(spatial_pool_plot)
  save_cell_plot(spatial_pool_plot, "step10_1_eos_reference_pools", FIGURE_ROOT, PROJECT_ROOT, width = 10, height = 7)
}'''),
        md("""## 10.2 Quantify nearest-cell composition at k = 1 and k = 15

The k = 1 table answers the literal nearest-cell question. The k = 15 table is less sensitive to a single boundary cell and summarizes the local neighbourhood. These are descriptive spatial summaries, not independent-cell hypothesis tests."""),
        code('''spatial_knn_edges <- list(k1 = data.frame(), k15 = data.frame())
knn_composition_k1 <- list(overall = data.frame(), by_state = data.frame())
knn_composition_k15 <- list(overall = data.frame(), by_state = data.frame())
if (RUN_EOS_SPATIAL) {
  spatial_knn_edges <- calculate_eos_knn_edges(spatial_pools, k_values = c(1L, 15L))
  knn_composition_k1 <- summarise_eos_knn_composition(spatial_knn_edges$k1)
  knn_composition_k15 <- summarise_eos_knn_composition(spatial_knn_edges$k15)
}
print(knn_composition_k1$overall %>% arrange(desc(fraction)))
print(knn_composition_k15$overall %>% arrange(desc(fraction)))
print(knn_composition_k15$by_state %>% arrange(EosState_extreme, desc(fraction)))'''),
        code('''if (nrow(knn_composition_k15$overall)) {
  make_signed_composition <- function(by_state) {
    by_state %>%
      filter(EosState_extreme %in% c("Short-lived-like", "Long-lived-like")) %>%
      mutate(signed_fraction = ifelse(EosState_extreme == "Short-lived-like", -fraction, fraction)) %>%
      group_by(reference_cell_type) %>%
      mutate(order_score = sum(signed_fraction, na.rm = TRUE)) %>% ungroup() %>%
      arrange(desc(order_score), desc(signed_fraction)) %>%
      mutate(reference_cell_type = factor(reference_cell_type, levels = rev(unique(reference_cell_type))))
  }
  plot_k1_data <- make_signed_composition(knn_composition_k1$by_state)
  plot_k15_data <- make_signed_composition(knn_composition_k15$by_state)
  plot_k1 <- ggplot(plot_k1_data, aes(signed_fraction, reference_cell_type, fill = EosState_extreme)) +
    geom_vline(xintercept = 0, colour = "grey70") + geom_col() +
    scale_fill_manual(values = c("Short-lived-like" = "#3979A8", "Long-lived-like" = "#B84E4B")) +
    labs(title = "Step 10.2: nearest-cell composition (k = 1)", x = "Fraction (short-like left; long-like right)", y = NULL, fill = NULL)
  plot_k15 <- ggplot(plot_k15_data, aes(signed_fraction, reference_cell_type, fill = EosState_extreme)) +
    geom_vline(xintercept = 0, colour = "grey70") + geom_col() +
    scale_fill_manual(values = c("Short-lived-like" = "#3979A8", "Long-lived-like" = "#B84E4B")) +
    labs(title = "Step 10.2: neighbourhood composition (k = 15)", x = "Fraction of edges (short-like left; long-like right)", y = NULL, fill = NULL)
  print(style_cell_plot(plot_k1))
  print(style_cell_plot(plot_k15))
  save_cell_plot(style_cell_plot(plot_k1), "step10_2_k1_neighbour_composition", FIGURE_ROOT, PROJECT_ROOT, width = 9, height = 7)
  save_cell_plot(style_cell_plot(plot_k15), "step10_2_k15_neighbour_by_state", FIGURE_ROOT, PROJECT_ROOT, width = 10, height = 7)
}'''),
        md("""## 10.3 Estimate Eosinophil distance to every adequately represented cell type

For every Eosinophil, the minimum Euclidean centroid distance to each cell type containing at least 20 reference cells is calculated. The boxplot summarizes proximity; the state–distance panel is descriptive because spatial autocorrelation invalidates ordinary cell-level p-values."""),
        code('''eos_distance_result <- if (RUN_EOS_SPATIAL) {
  calculate_eos_distance_by_cell_type(spatial_pools, min_reference_cells = 20L)
} else {
  list(cell_level = data.frame(), summary = data.frame())
}
distance_to_each_cell_type <- eos_distance_result$cell_level
distance_state_correlations <- eos_distance_result$summary
distance_tail_summary <- data.frame()
print(distance_state_correlations %>% arrange(median_distance))'''),
        code('''if (nrow(distance_to_each_cell_type)) {
  distance_tail_summary <- distance_to_each_cell_type %>%
    filter(EosState_extreme %in% c("Short-lived-like", "Long-lived-like")) %>%
    group_by(reference_cell_type, EosState_extreme) %>%
    summarise(median_distance = median(distance, na.rm = TRUE), .groups = "drop") %>%
    tidyr::complete(reference_cell_type,
                    EosState_extreme = c("Short-lived-like", "Long-lived-like"),
                    fill = list(median_distance = NA_real_)) %>%
    tidyr::pivot_wider(names_from = EosState_extreme, values_from = median_distance) %>%
    mutate(proximity_difference = `Short-lived-like` - `Long-lived-like`) %>%
    arrange(desc(proximity_difference))
  distance_order <- rev(distance_tail_summary$reference_cell_type)
  distance_boxplot <- distance_to_each_cell_type %>%
    filter(EosState_extreme %in% c("Short-lived-like", "Long-lived-like")) %>%
    mutate(reference_cell_type = factor(reference_cell_type, levels = distance_order)) %>%
    ggplot(aes(distance, reference_cell_type, fill = EosState_extreme)) +
    geom_boxplot(outlier.shape = NA, linewidth = 0.25) +
    scale_fill_manual(values = c("Short-lived-like" = "#3979A8", "Long-lived-like" = "#B84E4B")) +
    labs(title = "Step 10.3: tail-specific distance to each reference cell type", x = "Minimum centroid distance", y = NULL, fill = NULL)
  distance_state_plot <- distance_tail_summary %>%
    mutate(reference_cell_type = factor(reference_cell_type, levels = distance_order),
           direction = ifelse(proximity_difference >= 0, "Long-like closer", "Short-like closer")) %>%
    ggplot(aes(proximity_difference, reference_cell_type, colour = direction)) +
    geom_vline(xintercept = 0, colour = "grey70") +
    geom_segment(aes(x = 0, xend = proximity_difference, yend = reference_cell_type), linewidth = 0.75) +
    geom_point(size = 3) +
    scale_colour_manual(values = c("Short-like closer" = "#3979A8", "Long-like closer" = "#B84E4B")) +
    labs(title = "Step 10.3: relative proximity of Eosinophil-state tails",
         x = "Median distance difference (short - long)", y = NULL, colour = NULL)
  print(style_cell_plot(distance_boxplot))
  print(style_cell_plot(distance_state_plot))
  save_cell_plot(style_cell_plot(distance_boxplot), "step10_3_distance_by_cell_type", FIGURE_ROOT, PROJECT_ROOT, width = 10, height = 8)
  save_cell_plot(style_cell_plot(distance_state_plot), "step10_3_distance_by_continuous_state", FIGURE_ROOT, PROJECT_ROOT, width = 14, height = 10)
}'''),
        md("""## 10.4 Rank neighbour cell types along the continuous Eosinophil-state axis

Neighbour fractions are correlated with the continuous state score. Negative rho denotes short-like association and positive rho long-like association. The five strongest informative types per direction define the prespecified recipient/sender types for the next CellChat diagnostic."""),
        code('''eos_knn_state_association <- if (RUN_EOS_SPATIAL) {
  rank_eos_state_knn_associations(
    spatial_knn_edges$k15, biological_order = scwat_cell_type_order(),
    min_eos = 20L, top_n = 5L
  )
} else {
  list(status = "SKIPPED_NO_SPATIAL_EDGES", full = data.frame(), per_eos = data.frame(),
       short_top = character(), long_top = character())
}
print(eos_knn_state_association$full %>% arrange(desc(abs(spearman_rho))))
cat("Short-like top neighbours:", paste(eos_knn_state_association$short_top, collapse = ", "), "\n")
cat("Long-like top neighbours:", paste(eos_knn_state_association$long_top, collapse = ", "), "\n")'''),
        code('''if (identical(eos_knn_state_association$status, "PASS")) {
  association_plot_data <- eos_knn_state_association$full %>%
    filter(is.finite(spearman_rho)) %>%
    arrange(desc(spearman_rho)) %>%
    mutate(reference_cell_type = factor(reference_cell_type, levels = rev(unique(reference_cell_type))))
  association_plot <- ggplot(association_plot_data, aes(spearman_rho, reference_cell_type, colour = direction)) +
    geom_vline(xintercept = 0, colour = "grey75") + geom_segment(aes(x = 0, xend = spearman_rho, yend = reference_cell_type), linewidth = 0.7) +
    geom_point(size = 3) +
    scale_colour_manual(values = c("SHORT_ASSOCIATED" = "#3979A8", "LONG_ASSOCIATED" = "#B84E4B", "NEUTRAL" = "#C8C8C8")) +
    labs(title = "Step 10.4: KNN composition association with continuous Eosinophil state",
         x = "Spearman rho", y = NULL, colour = NULL)
  top_types <- unique(c(eos_knn_state_association$short_top, eos_knn_state_association$long_top))
  per_eos_plot <- eos_knn_state_association$per_eos %>%
    filter(reference_cell_type %in% top_types,
           EosState_extreme %in% c("Short-lived-like", "Long-lived-like")) %>%
    mutate(reference_cell_type = factor(reference_cell_type, levels = scwat_cell_type_order(top_types))) %>%
    ggplot(aes(EosState_balance, neighbour_fraction, colour = EosState_extreme)) +
    geom_point(alpha = 0.70, size = 1.15) + geom_smooth(method = "lm", se = FALSE, linewidth = 0.7) +
    facet_wrap(~reference_cell_type) +
    scale_colour_manual(values = c("Short-lived-like" = "#3979A8", "Long-lived-like" = "#B84E4B")) +
    labs(title = "Step 10.4: per-Eosinophil neighbour fraction", x = "EosState balance", y = "Fraction among k = 15")
  print(style_cell_plot(association_plot))
  print(style_cell_plot(per_eos_plot))
  save_cell_plot(style_cell_plot(association_plot), "step10_4_knn_state_association", FIGURE_ROOT, PROJECT_ROOT, width = 10, height = 8)
  save_cell_plot(style_cell_plot(per_eos_plot), "step10_4_per_eos_neighbour_fraction", FIGURE_ROOT, PROJECT_ROOT, width = 13, height = 8)
}'''),
        md("""## 11. Spatial CellChat for the top five short-like and long-like neighbour types

CellChat is run on normalized Xenium expression, aligned centroid coordinates and a physical scale derived from median cell area. Because CellChat requires categorical groups, the lower and upper 30% of the continuous Eosinophil-state score define `Eos_short_enriched` and `Eos_long_enriched`; the continuous KNN ranking remains the selection analysis. The significant report is restricted to outgoing Eosinophil-state signals toward that state's top-five neighbours; the unfiltered table is retained for audit. Results are exploratory within-section communication probabilities, not mouse-level inference or proof of signalling."""),
        code('''cellchat_inputs <- prepare_eos_cellchat_inputs(
  region_spatial,
  coordinates %>% transmute(cell_id = cell, x, y),
  top_short = eos_knn_state_association$short_top,
  top_long = eos_knn_state_association$long_top,
  eos_col = "Eos_inclusive", cell_type_col = "Final_CellType_subtype",
  state_col = "EosState_balance", cell_area_col = "cell_area",
  assay = "Xenium", min_cells = 10L
)
cellchat_result <- if (RUN_CELLCHAT) {
  run_eos_spatial_cellchat(cellchat_inputs, seed = RANDOM_SEED, min_cells = 10L)
} else {
  list(status = "SKIPPED_PACKAGE_OR_CONFIGURATION", message = "CellChat is unavailable or disabled.",
       package_version = if (optional_available[["CellChat"]]) as.character(packageVersion("CellChat")) else NA_character_,
       stage = "CONFIGURATION_GATE", object = NULL, communication = data.frame())
}
cellchat_status <- data.frame(
  status = cellchat_result$status, stage = cellchat_result$stage,
  package_version = cellchat_result$package_version, message = cellchat_result$message,
  state_tail_proportion = CELLCHAT_STATE_PROPORTION,
  analysis_scope = "EXPLORATORY_WITHIN_SECTION_CELLCHAT",
  biological_replicate_inference = FALSE, stringsAsFactors = FALSE
)
cellchat_filtered <- if (all(c("source", "target", "interaction_name", "prob", "pval") %in%
                             colnames(cellchat_result$communication))) {
  filter_eos_cellchat_interactions(
    cellchat_result$communication, raw_p_max = 0.05, adjusted_p_max = 0.10,
    top_short = eos_knn_state_association$short_top,
    top_long = eos_knn_state_association$long_top
  )
} else {
  list(all = data.frame(), significant = data.frame())
}
cellchat_all_interactions <- cellchat_filtered$all
cellchat_significant_interactions <- cellchat_filtered$significant
print(cellchat_status)
if (nrow(cellchat_significant_interactions)) {
  print(cellchat_significant_interactions %>% arrange(p_adjust_bh) %>% head(30))
}'''),
        code('''if (nrow(cellchat_significant_interactions)) {
  lr_plot_data <- cellchat_significant_interactions %>%
    mutate(pair_label = paste0(interaction_name, " | ", source, " → ", target)) %>%
    slice_min(order_by = p_adjust_bh, n = 30, with_ties = FALSE) %>%
    mutate(pair_label = factor(pair_label, levels = rev(unique(pair_label))))
  cellchat_dotplot <- ggplot(lr_plot_data, aes(prob, pair_label, colour = pathway_name, size = -log10(pmax(p_adjust_bh, 1e-300)))) +
    geom_point(alpha = 0.85) + scale_colour_manual(values = cell_macaron_palette(unique(lr_plot_data$pathway_name))) +
    labs(title = "Exploratory spatial CellChat: significant Eosinophil-focused pairs",
         x = "CellChat communication probability", y = NULL, colour = "Pathway", size = "-log10 FDR")
  cellchat_heatmap <- ggplot(lr_plot_data, aes(source, pair_label, fill = prob)) +
    geom_tile(colour = "white") + facet_wrap(~target, scales = "free_x") +
    scale_fill_gradient(low = "#F8ECEF", high = "#9D5C73") +
    labs(title = "Exploratory spatial CellChat probability by source and target", x = NULL, y = NULL, fill = "Probability")
  print(style_cell_plot(cellchat_dotplot))
  print(style_cell_plot(cellchat_heatmap))
  save_cell_plot(style_cell_plot(cellchat_dotplot), "eos_spatial_cellchat_significant_dotplot", FIGURE_ROOT, PROJECT_ROOT, width = 13, height = 10)
  save_cell_plot(style_cell_plot(cellchat_heatmap), "eos_spatial_cellchat_significant_heatmap", FIGURE_ROOT, PROJECT_ROOT, width = 13, height = 10)
} else {
  message("No CellChat interaction passed the prespecified raw-p and BH-FDR gates, or CellChat was skipped.")
}'''),
    ])

    if not child_branch:
        cells.extend([
            md("""## 12. Derive, visualize and freeze the lymph-node domain

The preliminary core is based on local enrichment of confidently annotated lymphoid/DC cells. DBSCAN is then run only on those core coordinates with `eps = 80` and `minPts = 10`; cluster 0 is noise, and only the largest non-noise cluster is retained as the lymph-node domain. There is no spatial expansion beyond that cluster. Multiple local-enrichment settings are compared; poor boundary agreement is a review gate. If the retained largest cluster is too small, the status is `LN_NOT_DETECTED` and no LN-only analysis is forced."""),
            code('''ln_cell_types <- c("B", "T", "γδ T", "NK", "Plasma", "cDC1", "cDC2", "CCR7_migratory_DC")
ln_input <- coordinates %>% transmute(
  cell_id = cell, x = x, y = y,
  cell_type = ifelse(
    region_spatial$Final_annotation_confidence[match(cell, colnames(region_spatial))] == "REVIEW",
    NA_character_,
    as.character(region_spatial$Final_CellType_subtype[match(cell, colnames(region_spatial))])
  )
)
ln_parameter_grid <- tribble(
  ~setting, ~k, ~fraction,
  "primary", 30L, 0.50,
  "lower_fraction", 30L, 0.40,
  "higher_fraction", 30L, 0.60,
  "smaller_neighbourhood", 20L, 0.50,
  "larger_neighbourhood", 40L, 0.50
)
ln_results <- setNames(lapply(seq_len(nrow(ln_parameter_grid)), function(i) {
  derive_lymph_node_domain(
    ln_input, lymphoid_labels = ln_cell_types,
    k = ln_parameter_grid$k[[i]],
    lymphoid_fraction_threshold = ln_parameter_grid$fraction[[i]],
    min_core_cells = 100L, dbscan_eps = 80, dbscan_min_pts = 10L
  )
}), ln_parameter_grid$setting)
ln_counts <- bind_rows(lapply(names(ln_results), function(setting) cbind(
  setting = setting, ln_results[[setting]]$parameters, ln_results[[setting]]$counts,
  largest_cluster_id = ln_results[[setting]]$largest_cluster_id,
  status = ln_results[[setting]]$status
)))
ln_counts
primary_dbscan_cluster_sizes <- ln_results$primary$dbscan_cluster_sizes
primary_dbscan_cluster_sizes'''),
            code('''ln_masks <- lapply(ln_results, function(result) result$cell_table$lymph_node_include)
ln_boundary_sensitivity <- summarise_ln_boundary_sensitivity(ln_masks)
ln_boundary_sensitivity$pairwise
ln_boundary_sensitivity$summary
primary_ln <- ln_results$primary
ln_plot_data <- cbind(
  ln_input,
  primary_ln$cell_table[, c("direct_lymphoid_evidence", "local_lymphoid_fraction",
                            "lymph_node_preliminary_core", "lymph_node_dbscan_cluster",
                            "lymph_node_core", "lymph_node_include")]
)
options(repr.plot.width = 15, repr.plot.height = 6)
ln_plot_data$lymph_node_include_display <- as.character(ln_plot_data$lymph_node_include)
print(plot_target_overlay(
  ln_plot_data, x_col = "x", y_col = "y", group_col = "lymph_node_include_display",
  target_labels = "TRUE", palette = c("FALSE" = "#D9D9D9", "TRUE" = "#694F85"),
  background_size = 0.10, target_size = 0.85,
  title = paste(REGION_ID, primary_ln$status, "largest DBSCAN lymph-node cluster"),
  legend_title = "LN domain", fixed_coordinates = TRUE
))
ggplot(ln_plot_data, aes(x, y, colour = local_lymphoid_fraction)) +
  geom_point(size = 0.15) + coord_fixed() + scale_colour_gradient(low = "#F8EEF2", high = "#9FC8D8") +
  labs(title = "Local confident lymphoid/DC fraction", colour = "Local fraction") + theme_void()'''),
            code('''ln_cell_ids <- primary_ln$cell_table$cell_id[primary_ln$cell_table$lymph_node_include]
masks$primary_include_revised <- derive_primary_include_revised(masks)
domain_manifest <- build_tissue_branch_manifest(
  masks, lymph_node_cell_ids = ln_cell_ids,
  provenance = paste(
    RUN_LABEL, REGION_ID, "k=30", "lymphoid_fraction>=0.50", "dbscan_eps=80",
    "dbscan_minPts=10", paste0("largest_cluster=", primary_ln$largest_cluster_id),
    primary_ln$status, sep = "::"
  )
)
stopifnot(
  identical(domain_manifest$cell_id, masks$cell_id),
  all((domain_manifest$adipose_only_include | domain_manifest$lymph_node_only_include) == domain_manifest$all_qcpass_include)
)
write_gz_tsv(domain_manifest, file.path(REGION_SHARED_ROOT, "lymph_node_domain_manifest.tsv.gz"), PROJECT_ROOT)
write_tsv(ln_boundary_sensitivity$pairwise, file.path(REGION_SHARED_ROOT, "lymph_node_boundary_sensitivity.tsv"), PROJECT_ROOT)
write_tsv(ln_counts, file.path(REGION_SHARED_ROOT, "lymph_node_boundary_counts.tsv"), PROJECT_ROOT)
write_tsv(primary_dbscan_cluster_sizes, file.path(REGION_SHARED_ROOT, "lymph_node_dbscan_cluster_sizes.tsv"), PROJECT_ROOT)
write_gz_tsv(primary_ln$cell_table, file.path(REGION_SHARED_ROOT, "lymph_node_domain_diagnostics.tsv.gz"), PROJECT_ROOT)
table(domain_manifest$tissue_domain)'''),
        ])

    cells.extend([
        md(f"""## {13 if not child_branch else 12}. Save branch-specific audit outputs

Objects and tables use unique stage/branch paths; no checkpoint is silently overwritten. The raw object remains unchanged. These results are section-level diagnostics and must pass the Region 3 anchor/admission workflow before primary release."""),
        code('''write_tsv(provenance, file.path(BRANCH_ROOT, "analysis_provenance.tsv"), PROJECT_ROOT)
write_tsv(optional_versions, file.path(BRANCH_ROOT, "optional_package_versions.tsv"), PROJECT_ROOT)
write_tsv(input_inventory, file.path(BRANCH_ROOT, "input_inventory.tsv"), PROJECT_ROOT)
write_tsv(qc_counts, file.path(BRANCH_ROOT, "qc_counts.tsv"), PROJECT_ROOT)
write_tsv(pca_qc_correlations, file.path(BRANCH_ROOT, "pca_qc_correlations.tsv"), PROJECT_ROOT)
write_tsv(cluster_stability, file.path(BRANCH_ROOT, "cluster_stability_pairwise.tsv"), PROJECT_ROOT)
write_tsv(cluster_stability_summary, file.path(BRANCH_ROOT, "cluster_stability_summary.tsv"), PROJECT_ROOT)
write_tsv(cluster_sizes, file.path(BRANCH_ROOT, "cluster_sizes.tsv"), PROJECT_ROOT)
write_tsv(cluster_qc_summary, file.path(BRANCH_ROOT, "cluster_qc_summary.tsv"), PROJECT_ROOT)
write_tsv(marker_coverage, file.path(BRANCH_ROOT, "canonical_marker_coverage.tsv"), PROJECT_ROOT)
write_tsv(reference_concordance_summary, file.path(BRANCH_ROOT, "wang_reference_concordance.tsv"), PROJECT_ROOT)
write_tsv(wang_shared_gene_gate$manifest, file.path(BRANCH_ROOT, "wang_xenium_shared_gene_manifest.tsv"), PROJECT_ROOT)
write_tsv(wang_xenium_integration_status, file.path(BRANCH_ROOT, "wang_xenium_integration_status.tsv"), PROJECT_ROOT)
write_tsv(wang_xenium_neighbour_summary$per_cell, file.path(BRANCH_ROOT, "wang_xenium_cross_dataset_neighbours.tsv"), PROJECT_ROOT)
write_tsv(wang_xenium_neighbour_summary$dataset_summary, file.path(BRANCH_ROOT, "wang_xenium_eos_dataset_summary.tsv"), PROJECT_ROOT)
write_tsv(wang_xenium_neighbour_summary$cluster_enrichment, file.path(BRANCH_ROOT, "wang_xenium_eos_cluster_enrichment.tsv"), PROJECT_ROOT)
write_tsv(eos_identity_summary, file.path(BRANCH_ROOT, "eos_identity_summary.tsv"), PROJECT_ROOT)
write_tsv(mclust_status, file.path(BRANCH_ROOT, "eos_state_mclust_status.tsv"), PROJECT_ROOT)
write_tsv(mixture_diagnostic$bic_table, file.path(BRANCH_ROOT, "eos_state_mclust_bic.tsv"), PROJECT_ROOT)
write_tsv(eos_state_association, file.path(BRANCH_ROOT, "eos_state_association_exploratory.tsv"), PROJECT_ROOT)
write_tsv(spatial_knn_edges$k1, file.path(BRANCH_ROOT, "eos_knn_k1_edges.tsv"), PROJECT_ROOT)
write_gz_tsv(spatial_knn_edges$k15, file.path(BRANCH_ROOT, "eos_knn_k15_edges.tsv.gz"), PROJECT_ROOT)
write_tsv(knn_composition_k1$overall, file.path(BRANCH_ROOT, "eos_k1_neighbour_composition.tsv"), PROJECT_ROOT)
write_tsv(knn_composition_k15$overall, file.path(BRANCH_ROOT, "eos_k15_neighbour_composition.tsv"), PROJECT_ROOT)
write_tsv(knn_composition_k15$by_state, file.path(BRANCH_ROOT, "eos_k15_neighbour_composition_by_state.tsv"), PROJECT_ROOT)
write_gz_tsv(distance_to_each_cell_type, file.path(BRANCH_ROOT, "eos_distance_to_each_cell_type.tsv.gz"), PROJECT_ROOT)
write_tsv(distance_state_correlations, file.path(BRANCH_ROOT, "eos_distance_state_correlations_descriptive.tsv"), PROJECT_ROOT)
write_tsv(distance_tail_summary, file.path(BRANCH_ROOT, "eos_distance_tail_comparison.tsv"), PROJECT_ROOT)
write_tsv(eos_knn_state_association$full, file.path(BRANCH_ROOT, "eos_knn_state_association.tsv"), PROJECT_ROOT)
write_gz_tsv(eos_knn_state_association$per_eos, file.path(BRANCH_ROOT, "eos_knn_state_per_cell.tsv.gz"), PROJECT_ROOT)
write_tsv(cellchat_status, file.path(BRANCH_ROOT, "eos_spatial_cellchat_status.tsv"), PROJECT_ROOT)
write_tsv(cellchat_all_interactions, file.path(BRANCH_ROOT, "eos_spatial_cellchat_all.tsv"), PROJECT_ROOT)
write_tsv(cellchat_significant_interactions, file.path(BRANCH_ROOT, "eos_spatial_cellchat_significant.tsv"), PROJECT_ROOT)
if (WRITE_CHECKPOINTS && identical(cellchat_result$status, "PASS")) {
  cellchat_checkpoint <- file.path(BRANCH_ROOT, "eos_spatial_cellchat_object.rds")
  assert_path_within(PROJECT_ROOT, cellchat_checkpoint)
  saveRDS(cellchat_result$object, cellchat_checkpoint)
}
checkpoint_records$final <- write_validated_seurat_checkpoint(
  region_spatial,
  file.path(BRANCH_ROOT, paste0(REGION_ID, "_", ANALYSIS_BRANCH, "_479_final.rds")),
  PROJECT_ROOT, stage = "04_FINAL_BRANCH_OBJECT"
)
checkpoint_manifest <- bind_rows(checkpoint_records)
write_tsv(checkpoint_manifest, file.path(BRANCH_ROOT, "checkpoint_manifest.tsv"), PROJECT_ROOT)
if (RUN_EOS_STATE) saveRDS(eos_obj, file.path(BRANCH_ROOT, paste0(REGION_ID, "_", ANALYSIS_BRANCH, "_eos.rds")))
writeLines(capture.output(sessionInfo()), file.path(BRANCH_ROOT, "sessionInfo.txt"))
list.files(BRANCH_ROOT)'''),
        md(f"""## {14 if not child_branch else 13}. Review gate and next step

Before accepting this branch, review:

- PC–QC correlations and cluster-level complexity;
- clustering stability and marker support;
- 2.5-month versus all-age Wang concordance;
- the fraction labelled `Uncertain`;
- Eosinophil rescue origin and state–complexity correlations;
- the mandatory disjoint-ID and positive-distance gates;
- lymph-node boundary sensitivity for the all-QC-passed notebook;
- dependence on technical-risk genes in the later 245/67-gene sensitivity runs.

`{section_role}` remains the governing section role. Region 4 outputs must remain separately labelled sensitivity-only."""),
    ])
    return cells


def write_notebook(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")


def main() -> None:
    NOTEBOOK_DIR.mkdir(parents=True, exist_ok=True)
    for region_number in range(1, 5):
        for branch in BRANCHES:
            name = f"B2_Region{region_number}_{branch}_479.ipynb"
            write_notebook(NOTEBOOK_DIR / name, make_notebook(notebook_cells(region_number, branch)))


if __name__ == "__main__":
    main()
