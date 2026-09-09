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
TEMP_ROOT <- file.path(DOWNSTREAM_ROOT, "tmp", REGION_ID, ANALYSIS_BRANCH)
set.seed(RANDOM_SEED)'''),
        code('''source(file.path(PIPELINE_REPO, "R", "source.R"))
validate_runtime_paths(PROJECT_ROOT, INPUT_ROOT, DOWNSTREAM_ROOT, TEMP_ROOT)
assert_path_within(PROJECT_ROOT, REGION_SHARED_ROOT)
assert_path_within(PROJECT_ROOT, BRANCH_ROOT)
dir.create(REGION_SHARED_ROOT, recursive = TRUE, showWarnings = FALSE)
dir.create(BRANCH_ROOT, recursive = TRUE, showWarnings = FALSE)
dir.create(TEMP_ROOT, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(TMPDIR = TEMP_ROOT, TMP = TEMP_ROOT, TEMP = TEMP_ROOT)

required_packages <- c(
  "Seurat", "SeuratObject", "Matrix", "dplyr", "tidyr", "tibble",
  "ggplot2", "FNN", "patchwork"
)
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
RUN_DIPTEST <- unname(optional_available[["diptest"]])
RUN_MCLUST <- unname(optional_available[["mclust"]])
RUN_LR_COEXPRESSION <- unname(optional_available[["CellChat"]])
print(optional_available)
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
ImageDimPlot(
  region_spatial_all, fov = "fov",
  group.by = c("qc_core_pass", "primary_include_revised"),
  flip_xy = FALSE, dark.background = FALSE,
  cols = c("FALSE" = "#D73027", "TRUE" = "#D9D9D9")
)'''),
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
ggplot(cluster_stability, aes(factor(resolution), ari)) +
  geom_boxplot(width = 0.55, outlier.shape = NA) + geom_jitter(width = 0.08, size = 2) +
  coord_cartesian(ylim = c(-0.05, 1.02)) +
  labs(title = "Clustering stability across random seeds", x = "Resolution", y = "Pairwise adjusted Rand index")'''),
        code('''WORKING_RESOLUTION <- as.numeric(Sys.getenv("SCWAT_WORKING_RESOLUTION", "0.8"))
stopifnot(WORKING_RESOLUTION %in% resolutions)
working_cluster_col <- cluster_column(WORKING_RESOLUTION, RANDOM_SEED)
region_spatial$working_cluster <- factor(region_spatial@meta.data[[working_cluster_col]])
region_spatial$cluster_res_0_8 <- region_spatial$working_cluster
cluster_sizes <- as.data.frame(table(cluster = region_spatial$working_cluster))
cluster_sizes
options(repr.plot.width = 12, repr.plot.height = 8)
DimPlot(region_spatial, reduction = "umap_PC30", group.by = "working_cluster", label = TRUE, repel = TRUE)'''),
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

The 2.5-month reference is primary because the study mice are eight weeks old. The all-age transfer is retained as a sensitivity comparison. Both use all shared genes from the 479-gene panel. Fine labels with review-level evidence receive an explicit `Uncertain` analysis label rather than a forced assignment."""),
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
table(ref_2p5$Wang_subtype_harmonized)'''),
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
  self_weight = 0.70, wang_prefix = "Ref2p5"
)
region_spatial <- annotation_result$object
region_spatial$Final_CellType_subtype_with_uncertain <- ifelse(
  region_spatial$Final_annotation_confidence == "REVIEW" |
    region_spatial$Final_review_flag %in% TRUE,
  "Uncertain", as.character(region_spatial$Final_CellType_subtype)
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
DimPlot(region_spatial, reduction = "umap_PC30", group.by = "Final_CellType_subtype_with_uncertain", label = TRUE, repel = TRUE)
ImageDimPlot(region_spatial, fov = "fov", group.by = "Final_CellType_subtype_with_uncertain", flip_xy = FALSE, dark.background = FALSE)
marker_plot_genes <- unique(marker_df_use$Gene_Symbol)
DotPlot(
  region_spatial, features = marker_plot_genes,
  group.by = "Final_CellType_subtype_with_uncertain", assay = "Xenium",
  dot.scale = 5, scale = TRUE, cols = c("#D9D9D9", "#5A2F5E")
) + RotatedAxis() + labs(
  title = "Canonical-marker consistency (internal, not independent validation)",
  x = NULL, y = NULL
)'''),
        md("""## 9. Eosinophil identity: annotation plus Tier 1/2 sensitivity set

The inclusive set is the union of existing Eosinophil annotation and Tier 1/2 evidence, matching the prior analysis. Crucially, rescued cells do not overwrite the principal cell-type label. Origin, evidence tier and confidence remain separate fields."""),
        code('''eos_result <- score_eosinophil_likeness(
  reference = ref_2p5, query = region_spatial,
  reference_group_col = "Wang_subtype_harmonized",
  wang_predicted_col = "Ref2p5_subtype_predicted.id",
  wang_eos_score_col = "Ref2p5_subtype_prediction.score.Eosinophil"
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
ImageDimPlot(region_spatial, fov = "fov", group.by = "Eos_origin", flip_xy = FALSE, dark.background = FALSE)'''),
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
        code('''if (RUN_EOS_STATE) {
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
  if (RUN_MCLUST) {
    mixture_diagnostic <- mclust::Mclust(eos_obj$EosState_balance, G = 1:3)
    print(summary(mixture_diagnostic))
  } else {
    message("Mixture diagnostic skipped: optional package 'mclust' is unavailable.")
  }
}'''),
        code('''if (RUN_EOS_STATE) {
  options(repr.plot.width = 13, repr.plot.height = 5)
  print(
    ggplot(eos_obj@meta.data, aes(EosState_balance, fill = EosState_extreme)) +
      geom_histogram(bins = 30, colour = "black") + geom_vline(xintercept = 0, linetype = 2) +
      labs(title = "Continuous Eosinophil-state balance", x = "Long-lived-like − short-lived-like", y = "Cells")
  )
  print(
    ggplot(eos_obj@meta.data, aes(EosShort_z, EosLong_z, colour = EosState_extreme)) +
      geom_point(alpha = 0.75) + geom_abline(slope = 1, intercept = 0, linetype = 2) +
      labs(title = "Eosinophil state is displayed as a continuum")
  )
  plot_eos_state_heatmap(eos_obj, eos_gene_sets)
}'''),
        md("""### 10.1 Exploratory within-section state-associated markers

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
      geom_point(data = eos_state_coordinates, aes(x, y, colour = EosState_balance), size = 1.2, alpha = 0.9) +
      scale_colour_gradient2(low = "#5B8DB8", mid = "lightyellow", high = "#D97A7A", midpoint = 0,
                             limits = c(-balance_limit, balance_limit)) +
      coord_fixed() + theme_void() +
      labs(title = "Spatial Eosinophil-state continuum", colour = "State balance")
  )
  print(
    ggplot() +
      geom_point(data = state_coordinates, aes(x, y), colour = "#E8E8E8", size = 0.10, alpha = 0.45) +
      geom_point(
        data = eos_state_coordinates %>% filter(EosState_extreme != "Intermediate"),
        aes(x, y, colour = EosState_extreme), size = 1.2, alpha = 0.9
      ) +
      scale_colour_manual(values = c("Short-lived-like" = "#7E9AD9", "Long-lived-like" = "#E89A8F")) +
      coord_fixed() + theme_void() + labs(title = "Descriptive Eosinophil-state tails", colour = "Tail")
  )
}'''),
        md("""## 11. Corrected spatial-neighbour analysis

Eosinophil query IDs and non-Eosinophil reference IDs are constructed from the same explicit `Eos_inclusive` field and must be disjoint. Cross-pool duplicate coordinates and zero/negative nearest-neighbour distances stop the analysis. Tables are descriptive; ordinary cell-level p-values are intentionally omitted because cells are spatially autocorrelated."""),
        code('''coordinates <- GetTissueCoordinates(region_spatial[["fov"]], which = "centroids") %>% as.data.frame()
if (!"cell" %in% colnames(coordinates)) coordinates$cell <- rownames(coordinates)
stopifnot(all(c("x", "y", "cell") %in% colnames(coordinates)))
coordinates <- coordinates %>%
  mutate(
    Eos_inclusive = region_spatial$Eos_inclusive[match(cell, colnames(region_spatial))],
    cell_type = region_spatial$Final_CellType_subtype_with_uncertain[match(cell, colnames(region_spatial))],
    EosState_balance = region_spatial$EosState_balance[match(cell, colnames(region_spatial))]
  )
query_coordinates <- coordinates %>% filter(Eos_inclusive %in% TRUE) %>% transmute(cell_id = cell, x, y)
reference_coordinates <- coordinates %>% filter(!(Eos_inclusive %in% TRUE)) %>% transmute(cell_id = cell, x, y)
RUN_EOS_SPATIAL <- nrow(query_coordinates) > 0L && nrow(reference_coordinates) > 0L
spatial_pool_gate <- if (RUN_EOS_SPATIAL) {
  validate_disjoint_spatial_pools(
    query_coordinates$cell_id, reference_coordinates$cell_id,
    query_coordinates, reference_coordinates
  )
} else {
  list(
    status = "SKIPPED_EMPTY_EOS_OR_REFERENCE_POOL",
    n_query = nrow(query_coordinates), n_reference = nrow(reference_coordinates)
  )
}
spatial_pool_gate'''),
        code('''nearest_eos <- data.frame(
  eos_cell_id = character(), reference_cell_id = character(),
  nearest_distance = numeric(), nearest_cell_type = character(),
  EosState_balance = numeric(), stringsAsFactors = FALSE
)
neighbour_composition <- data.frame(
  cell_type = character(), fraction_of_k15_edges = numeric(), stringsAsFactors = FALSE
)
neighbour_composition_by_state <- data.frame(
  EosState_extreme = character(), neighbour_cell_type = character(),
  n_edges = integer(), fraction = numeric(), stringsAsFactors = FALSE
)
neighbour_edge_table <- data.frame(
  eos_cell_id = character(), neighbour_rank = integer(),
  neighbour_cell_type = character(), EosState_extreme = character(),
  stringsAsFactors = FALSE
)
if (RUN_EOS_SPATIAL) {
  nearest_fit <- FNN::get.knnx(
    as.matrix(reference_coordinates[, c("x", "y")]),
    as.matrix(query_coordinates[, c("x", "y")]), k = min(15L, nrow(reference_coordinates))
  )
  distance_gate <- validate_neighbour_distances(nearest_fit$nn.dist[, 1L])
  nearest_reference_ids <- reference_coordinates$cell_id[nearest_fit$nn.index[, 1L]]
  nearest_eos <- data.frame(
    eos_cell_id = query_coordinates$cell_id,
    reference_cell_id = nearest_reference_ids,
    nearest_distance = nearest_fit$nn.dist[, 1L],
    nearest_cell_type = coordinates$cell_type[match(nearest_reference_ids, coordinates$cell)],
    EosState_balance = coordinates$EosState_balance[match(query_coordinates$cell_id, coordinates$cell)],
    stringsAsFactors = FALSE
  )
  neighbour_types <- matrix(
    coordinates$cell_type[match(reference_coordinates$cell_id[nearest_fit$nn.index], coordinates$cell)],
    nrow = nrow(nearest_fit$nn.index), ncol = ncol(nearest_fit$nn.index)
  )
  neighbour_edge_table <- data.frame(
    eos_cell_id = rep(query_coordinates$cell_id, each = ncol(neighbour_types)),
    neighbour_rank = rep(seq_len(ncol(neighbour_types)), times = nrow(neighbour_types)),
    neighbour_cell_type = as.vector(t(neighbour_types)),
    stringsAsFactors = FALSE
  ) %>% mutate(
    EosState_extreme = region_spatial$EosState_extreme[match(eos_cell_id, colnames(region_spatial))]
  )
  neighbour_composition <- as.data.frame(prop.table(table(as.vector(neighbour_types))))
  colnames(neighbour_composition) <- c("cell_type", "fraction_of_k15_edges")
  neighbour_composition_by_state <- neighbour_edge_table %>%
    filter(!is.na(EosState_extreme)) %>%
    count(EosState_extreme, neighbour_cell_type, name = "n_edges") %>%
    group_by(EosState_extreme) %>% mutate(fraction = n_edges / sum(n_edges)) %>% ungroup()
  print(distance_gate)
}
nearest_eos %>% count(nearest_cell_type, name = "n_eos") %>% arrange(desc(n_eos))
neighbour_composition %>% arrange(desc(fraction_of_k15_edges))
neighbour_composition_by_state %>% arrange(EosState_extreme, desc(fraction))'''),
        code('''if (nrow(nearest_eos)) {
  options(repr.plot.width = 12, repr.plot.height = 6)
  ggplot(nearest_eos, aes(reorder(nearest_cell_type, nearest_distance, median), nearest_distance)) +
    geom_boxplot(outlier.shape = NA) + coord_flip() +
    labs(title = "Nearest non-Eosinophil distance by cell type", x = NULL, y = "Distance in Xenium coordinate units")
}'''),
        code('''distance_to_each_cell_type <- data.frame(
  eos_cell_id = character(), reference_cell_type = character(),
  distance = numeric(), EosState_balance = numeric(), stringsAsFactors = FALSE
)
distance_state_correlations <- data.frame(
  reference_cell_type = character(), n_eos = integer(), spearman_rho = numeric(),
  interpretation = character(), stringsAsFactors = FALSE
)
if (RUN_EOS_SPATIAL) {
  eligible_reference_types <- names(which(table(coordinates$cell_type[!(coordinates$Eos_inclusive %in% TRUE)]) >= 20L))
  distance_to_each_cell_type <- bind_rows(lapply(eligible_reference_types, function(cell_type_name) {
    type_coordinates <- coordinates %>%
      filter(!(Eos_inclusive %in% TRUE), cell_type == cell_type_name) %>%
      select(x, y)
    fit <- FNN::get.knnx(
      as.matrix(type_coordinates), as.matrix(query_coordinates[, c("x", "y")]), k = 1L
    )
    validate_neighbour_distances(fit$nn.dist[, 1L])
    data.frame(
      eos_cell_id = query_coordinates$cell_id,
      reference_cell_type = cell_type_name,
      distance = fit$nn.dist[, 1L],
      EosState_balance = region_spatial$EosState_balance[match(query_coordinates$cell_id, colnames(region_spatial))],
      stringsAsFactors = FALSE
    )
  }))
  distance_state_correlations <- distance_to_each_cell_type %>%
    group_by(reference_cell_type) %>%
    summarise(
      n_eos = sum(complete.cases(distance, EosState_balance)),
      spearman_rho = if (n_eos >= 20L) cor(distance, EosState_balance, method = "spearman", use = "complete.obs") else NA_real_,
      interpretation = "DESCRIPTIVE_SPATIALLY_AUTOCORRELATED_NO_CELL_LEVEL_P_VALUE",
      .groups = "drop"
    ) %>% arrange(desc(abs(spearman_rho)))
}
distance_state_correlations'''),
        md("""## 12. Exploratory panel-limited ligand–receptor spatial co-expression

This is not CellChat inference and does not establish signalling. It reports expression products on the validated nearest Eosinophil–non-Eosinophil edges for simple ligand/receptor pairs present in the 479-gene panel. Complex-subunit interactions are excluded."""),
        code('''lr_scores <- data.frame(
  interaction_name = character(), pathway_name = character(),
  ligand = character(), receptor = character(), n_edges = integer(),
  pct_nonzero_edges = numeric(), mean_expression_product = numeric(),
  median_expression_product = numeric(), interpretation = character(),
  stringsAsFactors = FALSE
)
if (nrow(nearest_eos) && RUN_LR_COEXPRESSION) {
  lr_db <- CellChat::CellChatDB.mouse$interaction
  lr_pairs <- lr_db %>%
    filter(!grepl("_", ligand), !grepl("_", receptor), ligand %in% FIX_GENESET, receptor %in% FIX_GENESET) %>%
    distinct(interaction_name, pathway_name, ligand, receptor)
  normalized <- GetAssayData(region_spatial, assay = "Xenium", layer = "data")
  lr_scores <- bind_rows(lapply(seq_len(nrow(lr_pairs)), function(index) {
    pair <- lr_pairs[index, ]
    ligand_expression <- as.numeric(normalized[pair$ligand, nearest_eos$eos_cell_id])
    receptor_expression <- as.numeric(normalized[pair$receptor, nearest_eos$reference_cell_id])
    product <- ligand_expression * receptor_expression
    data.frame(
      interaction_name = pair$interaction_name, pathway_name = pair$pathway_name,
      ligand = pair$ligand, receptor = pair$receptor,
      n_edges = length(product), pct_nonzero_edges = 100 * mean(product > 0),
      mean_expression_product = mean(product), median_expression_product = median(product),
      interpretation = "EXPLORATORY_SPATIAL_COEXPRESSION_NOT_CELLCHAT_INFERENCE",
      stringsAsFactors = FALSE
    )
  }))
} else if (!RUN_LR_COEXPRESSION) {
  message("Ligand-receptor co-expression skipped: optional package 'CellChat' is unavailable.")
}
lr_scores %>% arrange(desc(mean_expression_product)) %>% head(30)'''),
    ])

    if not child_branch:
        cells.extend([
            md("""## 13. Derive, visualize and freeze the lymph-node domain

The candidate domain is based on local enrichment of confidently annotated lymphoid/DC cells and a limited spatial expansion that includes nearby stromal and vascular cells. It does not use a convex hull. Multiple parameter settings are compared; poor boundary agreement is a review gate. If insufficient core cells are found, the status is `LN_NOT_DETECTED` and no LN-only analysis should be forced."""),
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
  ~setting, ~k, ~fraction, ~radius,
  "primary", 30L, 0.50, 80,
  "lower_fraction", 30L, 0.40, 80,
  "higher_fraction", 30L, 0.60, 80,
  "smaller_neighbourhood", 20L, 0.50, 60,
  "larger_neighbourhood", 40L, 0.50, 100
)
ln_results <- setNames(lapply(seq_len(nrow(ln_parameter_grid)), function(i) {
  derive_lymph_node_domain(
    ln_input, lymphoid_labels = ln_cell_types,
    k = ln_parameter_grid$k[[i]],
    lymphoid_fraction_threshold = ln_parameter_grid$fraction[[i]],
    expansion_radius = ln_parameter_grid$radius[[i]], min_core_cells = 100L
  )
}), ln_parameter_grid$setting)
ln_counts <- bind_rows(lapply(names(ln_results), function(setting) cbind(
  setting = setting, ln_results[[setting]]$parameters, ln_results[[setting]]$counts,
  status = ln_results[[setting]]$status
)))
ln_counts'''),
            code('''ln_masks <- lapply(ln_results, function(result) result$cell_table$lymph_node_include)
ln_boundary_sensitivity <- summarise_ln_boundary_sensitivity(ln_masks)
ln_boundary_sensitivity$pairwise
ln_boundary_sensitivity$summary
primary_ln <- ln_results$primary
ln_plot_data <- cbind(
  ln_input,
  primary_ln$cell_table[, c("direct_lymphoid_evidence", "local_lymphoid_fraction", "lymph_node_core", "lymph_node_include")]
)
options(repr.plot.width = 15, repr.plot.height = 6)
ggplot(ln_plot_data, aes(x, y, colour = lymph_node_include)) +
  geom_point(size = 0.15, alpha = 0.7) + coord_fixed() +
  scale_colour_manual(values = c("FALSE" = "#D9D9D9", "TRUE" = "#8E5AA9")) +
  labs(title = paste(REGION_ID, primary_ln$status, "primary lymph-node domain"), colour = "LN domain") + theme_void()
ggplot(ln_plot_data, aes(x, y, colour = local_lymphoid_fraction)) +
  geom_point(size = 0.15) + coord_fixed() + scale_colour_viridis_c() +
  labs(title = "Local confident lymphoid/DC fraction", colour = "Local fraction") + theme_void()'''),
            code('''ln_cell_ids <- primary_ln$cell_table$cell_id[primary_ln$cell_table$lymph_node_include]
masks$primary_include_revised <- derive_primary_include_revised(masks)
domain_manifest <- build_tissue_branch_manifest(
  masks, lymph_node_cell_ids = ln_cell_ids,
  provenance = paste(
    RUN_LABEL, REGION_ID, "k=30", "lymphoid_fraction>=0.50", "expansion_radius=80",
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
write_gz_tsv(primary_ln$cell_table, file.path(REGION_SHARED_ROOT, "lymph_node_domain_diagnostics.tsv.gz"), PROJECT_ROOT)
table(domain_manifest$tissue_domain)'''),
        ])

    cells.extend([
        md("""## 14. Save branch-specific audit outputs

Objects and tables use unique stage/branch paths; no checkpoint is silently overwritten. The raw object remains unchanged. These results are section-level diagnostics and must pass the Region 3 anchor/admission workflow before primary release."""),
        code('''write_tsv(provenance, file.path(BRANCH_ROOT, "analysis_provenance.tsv"), PROJECT_ROOT)
write_tsv(input_inventory, file.path(BRANCH_ROOT, "input_inventory.tsv"), PROJECT_ROOT)
write_tsv(qc_counts, file.path(BRANCH_ROOT, "qc_counts.tsv"), PROJECT_ROOT)
write_tsv(pca_qc_correlations, file.path(BRANCH_ROOT, "pca_qc_correlations.tsv"), PROJECT_ROOT)
write_tsv(cluster_stability, file.path(BRANCH_ROOT, "cluster_stability_pairwise.tsv"), PROJECT_ROOT)
write_tsv(cluster_stability_summary, file.path(BRANCH_ROOT, "cluster_stability_summary.tsv"), PROJECT_ROOT)
write_tsv(cluster_sizes, file.path(BRANCH_ROOT, "cluster_sizes.tsv"), PROJECT_ROOT)
write_tsv(cluster_qc_summary, file.path(BRANCH_ROOT, "cluster_qc_summary.tsv"), PROJECT_ROOT)
write_tsv(marker_coverage, file.path(BRANCH_ROOT, "canonical_marker_coverage.tsv"), PROJECT_ROOT)
write_tsv(reference_concordance_summary, file.path(BRANCH_ROOT, "wang_reference_concordance.tsv"), PROJECT_ROOT)
write_tsv(eos_identity_summary, file.path(BRANCH_ROOT, "eos_identity_summary.tsv"), PROJECT_ROOT)
write_tsv(eos_state_association, file.path(BRANCH_ROOT, "eos_state_association_exploratory.tsv"), PROJECT_ROOT)
write_tsv(nearest_eos, file.path(BRANCH_ROOT, "eos_nearest_non_eos.tsv"), PROJECT_ROOT)
write_tsv(neighbour_composition, file.path(BRANCH_ROOT, "eos_k15_neighbour_composition.tsv"), PROJECT_ROOT)
write_tsv(neighbour_composition_by_state, file.path(BRANCH_ROOT, "eos_k15_neighbour_composition_by_state.tsv"), PROJECT_ROOT)
write_gz_tsv(distance_to_each_cell_type, file.path(BRANCH_ROOT, "eos_distance_to_each_cell_type.tsv.gz"), PROJECT_ROOT)
write_tsv(distance_state_correlations, file.path(BRANCH_ROOT, "eos_distance_state_correlations_descriptive.tsv"), PROJECT_ROOT)
write_tsv(lr_scores, file.path(BRANCH_ROOT, "eos_lr_spatial_coexpression_exploratory.tsv"), PROJECT_ROOT)
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
        md(f"""## 15. Review gate and next step

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
