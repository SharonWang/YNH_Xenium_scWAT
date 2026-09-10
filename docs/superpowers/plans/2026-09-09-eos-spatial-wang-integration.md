# Eosinophil Spatial and Wang Integration Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Extend every generated scWAT B2 full-panel notebook with robust Eosinophil diagnostics, biologically ordered cell/macaron figures, four explicit spatial-neighbour stages, gated spatial CellChat inference, and a supporting Wang–Xenium joint embedding.

**Architecture:** Reusable computations and plotting contracts live in `R/source.R`; the deterministic Python builder composes those functions into 12 output-free notebooks. Small behavior tests establish each R contract before notebook generation, while full Wang and CellChat runs remain explicit HPC checkpoints.

**Tech Stack:** R 4.6.1, Seurat/SeuratObject, Matrix, dplyr/tidyr/ggplot2, FNN, optional mclust/CellChat, Python notebook builder, Jupyter IRkernel.

**Spec:** `docs/superpowers/specs/2026-09-09-eos-spatial-wang-integration-design.md`

## Global Constraints

- Use all 479 Xenium genes for per-region normalization, scaling, PCA and clustering.
- Use PCs 1–30 and preserve the current `LogNormalize` workflow.
- Retain `Final_CellType_subtype_with_uncertain`; use `Final_CellType_subtype` for downstream grouping.
- Preserve raw objects and attach results by matched cell ID, never row order.
- Region 4 remains sensitivity-only and cannot define primary labels or clusters.
- CellChat output is `EXPLORATORY_WITHIN_SECTION_CELLCHAT`, not mouse-level inference.
- Keep all created files, caches, tests and temporary data under the D: project.
- Do not install packages locally; unavailable optional modules return typed skip statuses.

---

### Task 1: Safe mclust diagnostic

**Files:**
- Modify: `R/source.R`
- Create: `tests/test_eos_extended_helpers.R`

**Interfaces:**
- Consumes: numeric Eosinophil-state vector, candidate `G`, seed.
- Produces: `run_mclust_diagnostic(x, G = 1:3, seed = 1234L, min_n = 20L, min_per_component = 5L)` returning `status`, `message`, `package_version`, `n_input`, `n_finite`, `selected_G`, `model_name`, `fit`, and `bic_table`.

- [ ] **Step 1: Write the failing behavior test**

```r
set.seed(10)
x <- c(rnorm(30, -1, 0.2), rnorm(30, 1, 0.2))
diagnostic <- run_mclust_diagnostic(x, G = 1:3, seed = 10)
if (requireNamespace("mclust", quietly = TRUE)) {
  stopifnot(
    diagnostic$status == "PASS",
    diagnostic$n_finite == 60L,
    diagnostic$selected_G %in% 1:3,
    is.data.frame(diagnostic$bic_table)
  )
}
small <- run_mclust_diagnostic(1:10, min_n = 20L)
stopifnot(small$status == "SKIPPED_INSUFFICIENT_DATA")
```

- [ ] **Step 2: Run the new test and confirm RED**

Run with D:-redirected `R_LIBS_USER`, `TMPDIR`, `TMP` and `TEMP`:

```powershell
& 'D:\Programs\R-4.6.1\bin\Rscript.exe' tests\test_eos_extended_helpers.R
```

Expected: failure because `run_mclust_diagnostic` is undefined.

- [ ] **Step 3: Implement the minimal wrapper with the caller-frame binding**

```r
run_mclust_diagnostic <- function(x, G = 1:3, seed = 1234L,
                                  min_n = 20L, min_per_component = 5L) {
  x <- as.numeric(x)
  finite <- is.finite(x)
  x_use <- x[finite]
  base <- list(
    status = NA_character_, message = NA_character_,
    package_version = if (requireNamespace("mclust", quietly = TRUE))
      as.character(utils::packageVersion("mclust")) else NA_character_,
    n_input = length(x), n_finite = length(x_use), selected_G = NA_integer_,
    model_name = NA_character_, fit = NULL, bic_table = data.frame()
  )
  if (!requireNamespace("mclust", quietly = TRUE)) {
    base$status <- "SKIPPED_PACKAGE_UNAVAILABLE"
    return(base)
  }
  if (length(x_use) < min_n) {
    base$status <- "SKIPPED_INSUFFICIENT_DATA"
    return(base)
  }
  G_use <- sort(unique(as.integer(G)))
  G_use <- G_use[G_use >= 1L & G_use * min_per_component <= length(x_use)]
  if (!length(G_use)) {
    base$status <- "SKIPPED_INSUFFICIENT_DATA"
    return(base)
  }
  set.seed(seed)
  result <- tryCatch({
    mclustBIC <- getExportedValue("mclust", "mclustBIC")
    fit <- mclust::Mclust(x_use, G = G_use, verbose = FALSE)
    list(fit = fit, bic = as.data.frame(fit$BIC))
  }, error = identity)
  if (inherits(result, "error")) {
    base$status <- "FAILED_MCLUST_RUNTIME"
    base$message <- conditionMessage(result)
    return(base)
  }
  base$status <- "PASS"
  base$message <- "Gaussian-mixture diagnostic completed; biological states are not inferred."
  base$selected_G <- as.integer(result$fit$G)
  base$model_name <- as.character(result$fit$modelName)
  base$fit <- result$fit
  base$bic_table <- result$bic
  base
}
```

- [ ] **Step 4: Run GREEN and existing source tests**

Run `test_eos_extended_helpers.R`, `test_source_contract.R` and `test_source.R`; require exit 0.

- [ ] **Step 5: Commit**

```powershell
git add R/source.R tests/test_eos_extended_helpers.R
git commit -m "fix: make mclust eos diagnostic namespace-safe"
```

---

### Task 2: Cell/macaron palette and biological ordering

**Files:**
- Modify: `R/source.R`
- Modify: `tests/test_eos_extended_helpers.R`

**Interfaces:**
- Produces: `scwat_cell_type_order(labels = NULL)`, `apply_scwat_cell_type_order(labels)`, `cell_macaron_palette(labels = NULL)`, `style_cell_plot(plot, base_size = 12, legend_position = "right")`, and `order_marker_features(marker_df, available_genes)`.

- [ ] **Step 1: Add RED tests for order, unknown-label retention and palette completeness**

```r
labels <- c("T", "Adipocyte", "Eosinophil", "Capillary_EC", "new_type")
ordered <- apply_scwat_cell_type_order(labels)
stopifnot(
  identical(levels(ordered), c("Adipocyte", "Capillary_EC", "Eosinophil", "T", "new_type")),
  identical(as.character(ordered), labels)
)
palette <- cell_macaron_palette(labels)
stopifnot(setequal(names(palette), unique(labels)), all(grepl("^#[0-9A-Fa-f]{6}$", palette)))
markers <- data.frame(
  Gene_Symbol = c("Cd3d", "Pck1", "Siglecf", "Kdr"),
  CellType_subtype = c("T", "Adipocyte", "Eosinophil", "Capillary_EC")
)
marker_order <- order_marker_features(markers, rownames_to_keep = markers$Gene_Symbol)
stopifnot(identical(marker_order$Gene_Symbol, c("Pck1", "Kdr", "Siglecf", "Cd3d")))
```

- [ ] **Step 2: Verify RED**

Expected: missing ordering/palette functions.

- [ ] **Step 3: Implement documented helpers**

Use a fixed compartment order, append unknown labels in first-observed order, and map known types to stable pastel blue, gold, orange, olive, pink and purple roots. Use deterministic HCL fallback colours only for unknown labels. `style_cell_plot` must accept ggplot and patchwork objects without mutating data.

- [ ] **Step 4: Verify GREEN and mutation cases**

Confirm tests fail if unknown labels are dropped, if Eosinophil sorts after T, or if palette names do not cover every observed label.

- [ ] **Step 5: Commit**

```powershell
git add R/source.R tests/test_eos_extended_helpers.R
git commit -m "feat: add scwat biological ordering and macaron plots"
```

---

### Task 3: Four-stage Eosinophil KNN analysis

**Files:**
- Modify: `R/source.R`
- Modify: `tests/test_eos_extended_helpers.R`

**Interfaces:**
- Produces:
  - `build_eos_spatial_pools(cell_metadata, coordinates, eos_col, cell_type_col, state_col)`;
  - `calculate_eos_knn_edges(pools, k_values = c(1L, 15L))`;
  - `summarise_eos_knn_composition(edges)`;
  - `calculate_eos_distance_by_cell_type(pools, min_reference_cells = 20L)`;
  - `rank_eos_state_knn_associations(edges_k15, biological_order, min_eos = 20L, top_n = 3L)`.

- [ ] **Step 1: Add a hand-checked RED fixture**

```r
spatial_fixture <- data.frame(
  cell_id = c("e1", "e2", "a1", "a2", "m1", "m2", "t1", "t2"),
  x = c(0, 10, 1, 11, 2, 12, 3, 13), y = 0,
  Eos_inclusive = c(TRUE, TRUE, rep(FALSE, 6)),
  Final_CellType_subtype = c("Eosinophil", "Eosinophil", "ASC", "ASC", "Macrophage", "Macrophage", "T", "T"),
  EosState_balance = c(-1, 1, rep(NA_real_, 6))
)
pools <- build_eos_spatial_pools(spatial_fixture, spatial_fixture[, c("cell_id", "x", "y")])
edges <- calculate_eos_knn_edges(pools, k_values = c(1L, 3L))
stopifnot(nrow(edges$k1) == 2L, nrow(edges$k3) == 6L, all(edges$k1$distance > 0))
ranking <- rank_eos_state_knn_associations(edges$k3, min_eos = 2L, top_n = 2L)
stopifnot(nrow(ranking$full) == 3L, length(ranking$short_top) <= 2L, length(ranking$long_top) <= 2L)
```

- [ ] **Step 2: Verify RED**

Expected: missing spatial-pipeline functions.

- [ ] **Step 3: Implement pool and edge contracts**

All joins use `cell_id`; call existing `validate_disjoint_spatial_pools()` and `validate_neighbour_distances()`. Edge rows contain `eos_cell_id`, `reference_cell_id`, `neighbour_rank`, `k`, `distance`, `reference_cell_type`, `EosState_balance` and `EosState_extreme`.

- [ ] **Step 4: Implement summaries, per-type distances and deterministic ranking**

For each Eosinophil and cell type, calculate `neighbour_fraction`; summarize Spearman rho without an inferential p-value; rank positive/negative directions with absolute rho, edge count and biological order as tie breakers.

- [ ] **Step 5: Verify GREEN plus safety failures**

Add failures for overlapping IDs, duplicated coordinates, zero distances and row-order perturbation. Run the focused and existing domain-helper suites.

- [ ] **Step 6: Commit**

```powershell
git add R/source.R tests/test_eos_extended_helpers.R
git commit -m "feat: add staged eos spatial neighbour analysis"
```

---

### Task 4: CellChat state groups, inputs and result filtering

**Files:**
- Modify: `R/source.R`
- Modify: `tests/test_eos_extended_helpers.R`

**Interfaces:**
- Produces:
  - `derive_eos_cellchat_groups(state, lower_probability = 0.30, upper_probability = 0.70, min_cells = 10L)`;
  - `prepare_eos_cellchat_inputs(object, coordinates, top_short, top_long, ...)`;
  - `run_eos_spatial_cellchat(inputs, database = NULL, seed = 1234L, min_cells = 10L)`;
  - `filter_eos_cellchat_interactions(table, raw_p_max = 0.05, adjusted_p_max = 0.10)`.

- [ ] **Step 1: Add RED tests for deterministic tails and skip gates**

```r
groups <- derive_eos_cellchat_groups(1:100, min_cells = 10L)
stopifnot(
  sum(groups$group == "Eos_short_enriched", na.rm = TRUE) == 30L,
  sum(groups$group == "Eos_long_enriched", na.rm = TRUE) == 30L,
  groups$status == "PASS"
)
too_small <- derive_eos_cellchat_groups(1:20, min_cells = 10L)
stopifnot(too_small$status == "SKIPPED_INSUFFICIENT_STATE_GROUP_CELLS")
```

- [ ] **Step 2: Verify RED**

Expected: missing CellChat helper functions.

- [ ] **Step 3: Implement data-only group and input preparation**

Require disjoint Eosinophil/non-Eosinophil IDs, use the union of top-three types, ensure every group has at least 10 cells, calculate `spot = 1` and `spot.diameter = median(2 * sqrt(cell_area / pi))`, and preserve a sparse normalized expression matrix.

- [ ] **Step 4: Implement an optional runtime adapter**

The adapter checks the installed CellChat version/API, creates a spatial object, runs the supported standard pipeline, and always returns a typed status plus raw extracted communication table. Runtime errors are captured with package version and stage name.

- [ ] **Step 5: Implement filtering and test it independently of CellChat**

```r
fixture <- data.frame(
  source = c("Eos_short_enriched", "Eos_long_enriched"),
  target = c("ASC", "Macrophage"), interaction_name = c("A_B", "C_D"),
  prob = c(0.4, 0.2), pval = c(0.001, 0.20)
)
filtered <- filter_eos_cellchat_interactions(fixture)
stopifnot(nrow(filtered$significant) == 1L, filtered$significant$interaction_name == "A_B")
```

- [ ] **Step 6: Run GREEN and commit**

```powershell
git add R/source.R tests/test_eos_extended_helpers.R
git commit -m "feat: add gated spatial CellChat adapter"
```

---

### Task 5: Wang–Xenium integration preparation and concordance

**Files:**
- Modify: `R/source.R`
- Modify: `tests/test_eos_extended_helpers.R`

**Interfaces:**
- Produces:
  - `sample_wang_reference(reference, subtype_col, max_per_subtype = 1000L, seed = 1234L)`;
  - `prefix_seurat_cell_ids(object, prefix)`;
  - `validate_wang_xenium_shared_genes(reference, query, panel_genes, min_shared = 100L)`;
  - `run_wang_xenium_integration(reference, query, features, dims = 1:30, seed = 1234L)`;
  - `summarise_cross_dataset_eos_neighbours(embeddings, metadata, k = 15L)`.

- [ ] **Step 1: Add RED tests for balancing, ID safety and shared genes**

```r
sampled_ids <- sample_ids_by_group(
  ids = paste0("c", 1:12), groups = rep(c("Eosinophil", "ASC"), each = 6),
  max_per_group = 3L, seed = 5L
)
stopifnot(length(sampled_ids) == 6L, identical(sampled_ids, sample_ids_by_group(
  paste0("c", 1:12), rep(c("Eosinophil", "ASC"), each = 6), 3L, 5L
)))
shared <- validate_shared_feature_set(c("A", "B", "C"), c("B", "C", "D"), c("A", "B", "C", "D"), min_shared = 2L)
stopifnot(identical(shared$features, c("B", "C")))
```

- [ ] **Step 2: Verify RED**

Expected: helper functions undefined.

- [ ] **Step 3: Implement data-only preparation helpers**

Sampling is deterministic within subtype, IDs receive `WANG_` or `XENIUM_` prefixes before merge, and the shared-gene manifest records panel genes present/absent in each modality.

- [ ] **Step 4: Implement Seurat CCA integration wrapper**

Normalize separately, use the explicit shared feature vector in `FindIntegrationAnchors(reduction = "cca", dims = 1:30)`, call `IntegrateData`, scale, PCA and UMAP. Return status, object, anchors, parameters and messages without altering the original objects.

- [ ] **Step 5: Implement cross-dataset Eosinophil concordance**

Use FNN on integrated PCA/UMAP embeddings with reference pools restricted to the opposite dataset; report per-cell Eosinophil-neighbour fraction, modality mixing and integrated-cluster enrichment. Validate positive distances and unique prefixed IDs.

- [ ] **Step 6: Run GREEN and commit**

```powershell
git add R/source.R tests/test_eos_extended_helpers.R
git commit -m "feat: add Wang Xenium integration diagnostics"
```

---

### Task 6: Rebuild all 12 stepwise notebooks

**Files:**
- Modify: `scripts/build_b2_split_notebooks.py`
- Regenerate: `notebooks/B2_Region{1..4}_{all_QCpass,adipose_only,lymph_node_only}_479.ipynb`
- Modify: `tests/test_b2_split_notebooks.py`
- Modify: `tests/test_b2_notebook_r_syntax.R`

**Interfaces:**
- Consumes: Tasks 1–5 helpers.
- Produces: deterministic notebooks with visible stepwise outputs and branch-local TSV/PNG/PDF artifacts.

- [ ] **Step 1: Add RED notebook-contract tests**

For every generated notebook assert:

```python
assert "run_mclust_diagnostic(" in code
assert "Final_CellType_subtype_with_uncertain" in code
assert 'cell_type_col = "Final_CellType_subtype"' in code
assert "### 10.1" in markdown and "### 10.4" in markdown
assert "run_eos_spatial_cellchat(" in code
assert "run_wang_xenium_integration(" in code
assert "cell_macaron_palette(" in code
```

Also assert prohibited downstream uses of `Final_CellType_subtype_with_uncertain` are absent from KNN and CellChat blocks.

- [ ] **Step 2: Verify RED against current notebooks**

Expected: missing staged helpers and joint integration calls.

- [ ] **Step 3: Update setup and annotation plots**

Record optional package versions, keep uncertainty diagnostics, order `Final_CellType_subtype`, use matching ordered marker blocks, and apply explicit named macaron palettes to DimPlot/ImageDimPlot/DotPlot.

- [ ] **Step 4: Replace the inline mixture block**

Call `run_mclust_diagnostic()`, print a bounded status table and BIC table, and plot mixture density only when status is `PASS`.

- [ ] **Step 5: Replace current spatial analysis with Steps 10.1–10.4**

Each subsection gets a methods markdown cell, one preparation/statistics cell, one plot cell and bounded previews. Write all tables named in the specification.

- [ ] **Step 6: Add spatial CellChat cells**

Build state groups, prepare inputs, run the adapter, print status, filter results, plot only non-empty significant tables and write all status/result/parameter files.

- [ ] **Step 7: Add Wang–Xenium integration cells beside label transfer**

Run only when `SCWAT_RUN_WANG_INTEGRATION` defaults to `TRUE`; write a status even when skipped/failed. Display modality, subtype, Eosinophil and concordance plots and save the optional checkpoint.

- [ ] **Step 8: Regenerate and verify GREEN**

Run the builder twice and assert the second run has no Git diff. Run Python notebook contracts and parse every R cell with R 4.6.1.

- [ ] **Step 9: Commit**

```powershell
git add scripts/build_b2_split_notebooks.py notebooks tests/test_b2_split_notebooks.py tests/test_b2_notebook_r_syntax.R
git commit -m "feat: expand stepwise eos and Wang notebooks"
```

---

### Task 7: Bounded local execution and HPC handoff

**Files:**
- Modify: `tests/test_b2_local_subset_smoke.R`
- Modify: `SCWAT_QUERY_AND_RESEARCH_PLAN.md`

**Interfaces:**
- Consumes: generated Region 3 all-QC-pass notebook and subset input.
- Produces: local test evidence and exact HPC commands; no full-data claim.

- [ ] **Step 1: Add smoke assertions**

Require the subset workflow to retain 479 genes, 30 PCs, ordered final labels, disjoint spatial pools, positive distances, non-empty Step 10 tables when sufficient Eosinophils exist, and typed optional-module statuses otherwise.

- [ ] **Step 2: Run bounded Region 3 subset smoke test**

Use only D:-scoped input, output and temp roots. Do not run full Wang or CellChat if local data/packages are insufficient; validate their preparation/status contracts.

- [ ] **Step 3: Review generated plots when local data reach the plotting gates**

Inspect PNG outputs for readable labels, biological ordering, explicit denominators, consistent palettes and non-clipped legends. If the subset is insufficient, record visual QA as an HPC checkpoint.

- [ ] **Step 4: Update progress and HPC commands**

Record exact pass/fail/skip results and provide environment exports plus one notebook-at-a-time Jupyter commands. State which output files must be downloaded for review.

- [ ] **Step 5: Run final verification**

Run all feature R tests, all generated-notebook tests, `git diff --check`, JSON validation and `git status`. Run the legacy suite separately and report inherited failures without conflating them with feature regressions.

- [ ] **Step 6: Commit**

```powershell
git add tests/test_b2_local_subset_smoke.R SCWAT_QUERY_AND_RESEARCH_PLAN.md
git commit -m "test: verify extended eos workflow on local subset"
```

