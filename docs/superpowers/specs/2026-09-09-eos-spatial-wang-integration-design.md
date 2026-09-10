# scWAT Eosinophil spatial, CellChat and Wang-integration design

Date: 2026-09-09  
Status: approved in chat; written specification awaiting final review  
Scope: the 12 generated `B2_Region{1..4}_{all_QCpass,adipose_only,lymph_node_only}_479.ipynb` notebooks and their reusable R functions  
Execution: bounded local tests on D:; full data, Wang integration and CellChat on HPC

## 1. Objectives

The generated B2 notebooks will retain the complete 479-gene Xenium workflow while improving seven linked areas:

1. make the optional `mclust` diagnostic safe and reproducible;
2. apply one documented cell-style theme and explicit macaron palettes;
3. retain the uncertainty-aware annotation while using the non-forced final subtype for downstream analysis;
4. impose a biological cell-type order and align marker blocks with that order;
5. expand Eosinophil spatial-neighbour analysis into four independently reviewable steps;
6. add spatial CellChat inference for state-enriched Eosinophils and their top three KNN-associated cell types; and
7. add an exploratory Wang–Xenium joint embedding alongside the existing Wang label transfer.

No raw Xenium object will be modified. Existing QC masks, the 479-gene normalized assay and all annotation evidence fields will be retained.

## 2. Files and generation model

The source of truth will remain:

- `R/source.R` for reusable, documented functions;
- `scripts/build_b2_split_notebooks.py` for notebook content;
- generated B2 notebooks under `notebooks/`;
- behavior and notebook-contract tests under `tests/`; and
- `SCWAT_QUERY_AND_RESEARCH_PLAN.md` for progress and validation history.

Generated notebooks will not be edited independently. The builder will regenerate all 12 notebooks deterministically with empty outputs and stable cell ordering.

## 3. Optional mixture diagnostic

### Root cause

`mclust::Mclust()` rewrites its call to the unqualified symbol `mclustBIC` and evaluates that call in the caller environment. A namespace-qualified call can therefore report `could not find function "mclustBIC"` even though `requireNamespace("mclust")` succeeded.

### Design

Add `run_mclust_diagnostic()` with the following contract:

- input: a finite numeric state vector, candidate component counts and a seed;
- validation: remove named/dimensional attributes, record non-finite removals, require at least 20 observations and at least five observations per requested component;
- execution: bind the exported `mclustBIC` function in the local caller frame before invoking `mclust::Mclust()`;
- error handling: return a typed `SKIPPED_PACKAGE_UNAVAILABLE`, `SKIPPED_INSUFFICIENT_DATA` or `FAILED_MCLUST_RUNTIME` result rather than terminating the notebook;
- output: status, package version, observation count, fitted object, selected component count/model, BIC table and message.

The diagnostic remains exploratory. A selected multi-component model will not create or rename biological Eosinophil states.

## 4. Annotation fields and cell-type order

Both annotation fields remain in the object:

- `Final_CellType_subtype`: downstream analysis label;
- `Final_CellType_subtype_with_uncertain`: review display that replaces review-level assignments with `Uncertain`.

All clustering summaries, canonical-marker plots, Eosinophil reference pools, KNN analyses, cell-type distance analyses and CellChat groups will use `Final_CellType_subtype`. The uncertainty-aware field will be shown in a dedicated concordance/count diagnostic but will not define downstream groups.

Add `scwat_cell_type_order()` and `apply_scwat_cell_type_order()`. Known labels will be ordered by compartment:

1. adipocyte;
2. ASC/APC/fibroblast;
3. VSMC/pericyte and other mural cells;
4. capillary/venous/lymphatic and other endothelial cells;
5. Schwann, mesothelial and epithelial cells;
6. resident/scavenging macrophage, monocyte, DC and neutrophil;
7. Eosinophil and mast cell;
8. ILC/NK/T/γδ T/B/plasma cells; and
9. other or unmatched labels, retained in first-observed order at the end.

The canonical-marker data frame will be sorted by the same subtype order. The DotPlot will use a named feature list so marker blocks remain labelled and adjacent, and the identity factor will use the matching biological order.

## 5. Plot system

Add reusable palette and plotting helpers rather than notebook-local themes:

- `cell_macaron_palette()` returns stable named colours for cell types, Eosinophil origins and state groups;
- `scale_*_cell_macaron()` returns explicit discrete scales;
- `style_seurat_plot()` adds `cell_style_theme()`, consistent legends and plot margins to ggplot/patchwork-compatible Seurat plots;
- spatial helpers use a white background, fixed coordinates, pale-grey context cells and explicit focal colours;
- heatmaps use a muted blue–cream–coral diverging scale centered at zero.

Colour will not be the only distinction where state direction matters: titles, labels, ordering and facets will also encode the comparison. Every plot will carry a neutral descriptive title and the relevant denominator, k value or cell count in its subtitle/caption.

## 6. Eosinophil state

The existing inclusive definition is unchanged: final Eosinophil annotation OR Tier 1/2 Eosinophil evidence. This field selects Eosinophil query cells but does not overwrite `Final_CellType_subtype`.

The continuous primary state remains `EosState_balance`. Existing 10% descriptive extremes remain `EosState_extreme` and are not treated as discovered subtypes. For CellChat only, add `EosState_CCC_group` using the lower and upper 30% of finite balance values:

- lower 30%: `Eos_short_enriched`;
- upper 30%: `Eos_long_enriched`;
- middle 40%: excluded from state-stratified CellChat.

Each retained state group must contain at least 10 cells. Otherwise CellChat is skipped with a typed status and no forced result.

## 7. Four spatial-neighbour steps

### Step 10.1 — Spatial pools and safety gates

- Construct Eosinophil queries and non-Eosinophil references from the same `Eos_inclusive` field.
- Use `Final_CellType_subtype` for reference labels.
- Stop on overlapping IDs, duplicated cross-pool coordinates, zero/negative nearest distances or coordinate/cell-ID mismatch.
- Plot all tissue cells, Eosinophil queries and eligible non-Eosinophil references.
- Export pool counts and validation status.

### Step 10.2 — Nearest-neighbour composition

- Calculate k=1 and k=15 neighbours among non-Eosinophil reference cells.
- Export one row per Eosinophil-neighbour edge with cell IDs, rank, distance, cell type and state.
- Summarize counts and fractions overall and by Eosinophil state.
- Plot ranked k=1 counts, k=15 composition and state-stratified composition.

### Step 10.3 — Distance to every eligible cell type

- Eligible reference types require at least 20 cells.
- Calculate each Eosinophil cell's minimum distance to each eligible reference type.
- Summarize median and interquartile distance, cell counts and Spearman correlation with continuous `EosState_balance`.
- Plot ordered distributions, a cell-by-type distance heatmap and continuous state-versus-distance panels.
- Correlations are descriptive because cells within a section are spatially autocorrelated.

### Step 10.4 — Continuous-state KNN association and top-three selection

- For each Eosinophil cell and reference type, calculate the fraction among its k=15 neighbours.
- Calculate the Spearman association between the neighbour-type fraction and continuous `EosState_balance`.
- Rank long-associated types by positive correlation and short-associated types by negative correlation, requiring at least 20 Eosinophil cells with finite state.
- Select the first three unique types in each direction; ties are resolved by absolute correlation, overall edge count and biological order.
- Export the full ranking and the two top-three lists.
- Plot signed association ranks, state-versus-neighbour-fraction panels and spatial locations of Eosinophils linked to each selected type.

No ordinary cell-level association p-value will be presented as biological-replicate inference. The ranking is a spatially descriptive candidate-selection step for CellChat.

## 8. Spatial CellChat analysis

CellChat uses normalized, non-integrated Xenium expression, spatial coordinates and `Final_CellType_subtype`-derived groups. It will not use the integrated Wang–Xenium assay.

### Included cells and groups

- `Eos_short_enriched` and `Eos_long_enriched` Eosinophils from the prespecified 30% tails;
- the union of Step 10.4's short-associated and long-associated top-three non-Eosinophil cell types;
- at least 10 cells per group;
- no Eosinophil cell in any non-Eosinophil reference group.

### Spatial inputs

Xenium centroid coordinates are treated as micron-scale coordinates. `spot = 1`; `spot.diameter` is the median equivalent cell diameter, `2 * sqrt(cell_area / pi)`, among included cells. Both values and their derivation are exported.

### Inference and reporting

- Use the installed CellChat mouse database, normalized Xenium assay and spatial mode.
- Run the standard overexpressed-gene, overexpressed-interaction, communication-probability, filtering and pathway aggregation stages supported by the installed CellChat version.
- Preserve complex/cofactor requirements; do not reduce the database to simple one-gene ligand/receptor strings.
- Extract Eosinophil-to-neighbour and neighbour-to-Eosinophil directions separately.
- Retain CellChat probability and permutation p-value and add BH adjustment over the reported section-level candidate table.
- A reported interaction requires raw CellChat p < 0.05, BH-adjusted p < 0.10, at least 10 cells in both groups and measured required genes in the panel.
- Plot significant interaction dot plots, sender–receiver heatmaps and pathway summaries using the cell/macaron style.

These results are labelled `EXPLORATORY_WITHIN_SECTION_CELLCHAT`. They are not mouse-level differential communication because a single section does not provide biological replication. Cross-region recurrence will be summarized later without treating technical sections as independent mice.

## 9. Wang–Xenium joint integration

Existing age-matched and all-age Wang label transfers remain unchanged. A separate exploratory 2.5-month Wang–Xenium joint embedding will be added.

### Input balancing

- Use all branch Xenium cells.
- Deterministically sample at most 1,000 Wang cells per harmonized subtype, retaining all Wang Eosinophils up to that cap.
- Prefix cell IDs by dataset before merging.
- Use every gene shared by the Wang RNA assay and the 479-gene Xenium panel; record absent genes and require at least 100 shared genes.

### Integration

- Log-normalize each object independently.
- Find Seurat CCA integration anchors with the explicit shared-gene vector and dimensions 1–30.
- Integrate data, scale shared features, run PCA 1–30 and UMAP with the project seed.
- Do not replace raw RNA/Xenium assays or primary Xenium clustering.

### Validation and plots

- Plot integrated UMAP by dataset, harmonized Wang subtype, Xenium final subtype and Eosinophil evidence.
- Quantify local dataset mixing.
- For each Xenium Eosinophil candidate, calculate the fraction of Wang Eosinophils among its Wang-only integrated neighbours.
- For each Wang Eosinophil, calculate the fraction of Xenium Eosinophil candidates among its Xenium-only integrated neighbours.
- Report integrated-cluster Eosinophil enrichment, modality composition and neighbour concordance.
- Treat a shared Eosinophil neighbourhood as supporting evidence only; absence of one does not independently disprove Eosinophil identity because the targeted panel and assay modalities differ.

Region 4 remains sensitivity-only. Its integration cannot change any primary reference label or cluster.

## 10. Outputs

Each eligible notebook branch will write bounded tables under its existing `BRANCH_ROOT`:

- `mclust_eos_state_diagnostic.tsv`;
- `cell_type_display_order.tsv`;
- `step10_1_spatial_pool_gate.tsv`;
- `step10_2_knn_edges.tsv.gz` and composition summaries;
- `step10_3_distance_by_cell_type.tsv.gz` and summaries;
- `step10_4_state_knn_association.tsv` and top-three tables;
- `cellchat_run_status.tsv`;
- `cellchat_significant_interactions.tsv`;
- `cellchat_spatial_parameters.tsv`;
- `wang_xenium_integration_status.tsv`;
- `wang_xenium_eos_concordance.tsv`; and
- reload-validated RDS checkpoints for CellChat and the joint integration when checkpoint writing is enabled.

Plots remain visible step by step in the notebooks and are also saved as PDF and PNG under branch-specific `plots/` directories so HPC results can be downloaded and reviewed.

## 11. Failure and skip policy

Required core analysis still stops on unsafe paths, ID mismatch, overlapping Eosinophil/reference pools, invalid coordinates or zero distances. Optional modules return explicit status tables rather than silently failing:

- mclust unavailable, insufficient or broken runtime;
- fewer than 20 inclusive Eosinophils for continuous state analysis;
- fewer than 10 cells in either CellChat state group;
- fewer than three eligible neighbour cell types;
- CellChat unavailable or API-incompatible;
- Wang reference missing, insufficient shared genes or failed integration.

Optional-module failure does not invalidate the preceding 479-gene clustering and annotation results, but it prevents the corresponding optional conclusion.

## 12. Testing and execution

Implementation will be test-first. Local tests will cover:

- the exact namespace-qualified `mclust::Mclust()` failure and safe-wrapper result;
- stable biological ordering with unknown labels appended;
- palette completeness and deterministic mapping;
- use of `Final_CellType_subtype` rather than the uncertainty-display field;
- disjoint KNN edges, correct k=1/k=15 counts and deterministic top-three ranking;
- CellChat group construction and skip gates without requiring a full local CellChat run;
- Wang sampling, ID prefixing, shared-gene validation and concordance metrics;
- deterministic regeneration, JSON validity and parsing of every R notebook cell.

A bounded Region 3 subset will exercise all locally available modules. Full CellChat and Wang integration are HPC execution checkpoints when local package or memory constraints prevent end-to-end execution. All local temporary files, tests and generated artifacts remain under the D: project.

