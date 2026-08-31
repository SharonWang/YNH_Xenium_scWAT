# B1 Region 3 complete-panel method audit

**Artifact:** `notebooks/B1_Region3_primary_479.ipynb`  
**SHA-256:** `739BD6648CF48057E14699694556D33AD24E09110D0FDF31A6DD6109214D1037`  
**Executed scope:** Region 3; 75,655 imported cells; 479 Gene Expression features  
**Purpose:** Map every unique B1 analysis to its scientifically reviewed B2 treatment. B1 remains unchanged as execution provenance.

## Decision definitions

- `RETAIN`: method and interpretation are suitable; recompute cleanly in B2.
- `RETAIN_WITH_CAVEAT`: useful descriptive evidence, but B2 must add provenance, uncertainty or interpretation limits.
- `CORRECT`: calculation, execution dependency or claim is not reliable enough to carry forward unchanged; B2 must recompute it with explicit checks.

## Audit table

| B1 component | Evidence observed in B1 | Decision | B2 treatment |
|---|---|---|---|
| Package loading, paths and seed | R packages loaded and seed 1234 set; several paths/checkpoints are introduced interactively | CORRECT | One setup block, path containment, package/version table, hashes and one checkpoint contract |
| Native Xenium import | 479 genes and 75,655 cells; cell/nucleus polygons loaded | RETAIN_WITH_CAVEAT | Retain import and object summaries; surface the metadata-overlap warning and validate final alignment |
| QC mask attachment | 75,655/75,655 cells matched to mask table | RETAIN | Match by `cell_id` with complete-match and duplicate-ID hard gates |
| `primary_include_revised` | 73,537 included; 2,118 excluded | CORRECT | Preserve the rule but compute it from ID-aligned Seurat metadata, never directly from external-table row order |
| Pre-subset QC plots | Core-pass map, revised-mask map and individual flag counts | RETAIN | Retain counts, cross-tabs and spatial maps before subsetting |
| Complete-panel definition | All 479 genes selected; 93 short/long state genes identified | RETAIN | Use all 479 for normalization, PCA, graph, UMAP and clustering as user specified |
| Log normalization and scaling | `LogNormalize`, scale factor 10,000; all 479 genes scaled | RETAIN | Preserve with explicit assay/layer checks |
| PCA diagnostics | 30 PCs, elbow/loadings/heatmap and QC correlations; PC2 correlated 0.795 with nFeature and 0.602 with nCount | RETAIN_WITH_CAVEAT | Retain plots/tables; state QC correlations without assuming biological relevance |
| PCA naming | Downstream uses dimensions 1–30 but graph/UMAP names contain `PC12` | CORRECT | Use dimensions 1–30 and names containing `PC1_30` consistently |
| Multi-resolution clustering | Leiden resolutions 0.7, 0.8, 1.0 and 1.2; 0.8 used downstream | RETAIN_WITH_CAVEAT | Retain UMAP panels and add cluster-size, transition and adjusted-Rand stability tables |
| Curated marker coverage | All listed markers available; cluster-level top/second scores and margins displayed | RETAIN | Retain marker table, coverage, score margins, feature plots and spatial views |
| Wang reference preparation | 2.5-month and all-age reference counts shown; ontology harmonization explicit | RETAIN | Preserve label harmonization and reference counts |
| Wang transfer | 467 shared genes; separate 2.5-month and all-reference transfers | RETAIN_WITH_CAVEAT | Preserve both; report shared features and make age-nearest transfer principal, all-age transfer a concordance check |
| Consolidated annotation | Marker, reference and kNN evidence combined; 24,924/73,537 cells (33.89%) flagged for review | RETAIN_WITH_CAVEAT | Preserve evidence columns and review counts; label insufficient/conflicting evidence `Uncertain` |
| Initial annotated Eos | Consolidated annotation contains 70 Eos cells (66 non-review and 4 review) | RETAIN | Preserve original annotation and confidence as immutable evidence columns |
| Reference-supported Eos tiers | Tier 1–4 marker/reference performance tables and cross-tabs shown | RETAIN_WITH_CAVEAT | Preserve calibration and plots; explain low single-marker detection and non-Eos tier calls |
| Inclusive Eos union | Annotation Eos plus Tier 1/2 yields 481 cells | RETAIN_WITH_CAVEAT | Apply the user-specified union and add `Eos_origin` for annotation, Tier 1 rescue and Tier 2 rescue |
| Short/long score construction | 47 short, 46 long and 30 non-ribosomal short genes; scaled mean scores, z-scores, balance and activity | RETAIN_WITH_CAVEAT | Preserve formulas and gene lists; show availability and complexity correlations |
| Eos score distributions | Histograms/scatters, detection summaries and complexity diagnostics | RETAIN | Recompute and retain all unique diagnostic plots and tables |
| Modality assessment | Dip test p=0.5227; Gaussian mixture selected one component for 481 Eos | RETAIN | Retain and state that the data do not support discrete state modes |
| Quantile-defined Eos extremes | 40 short-like, 401 intermediate and 40 long-like cells | RETAIN_WITH_CAVEAT | Retain q10/q90 descriptive labels; do not call them discovered subtypes |
| Eos heatmap | State genes row-scaled and cells ordered by score | RETAIN_WITH_CAVEAT | Retain as a visualization of the defining signatures, not independent validation |
| Cell-level Eos expression comparison | 40 versus 40 cells; 379 non-state genes tested with Wilcoxon | CORRECT | Retain effect sizes, detection fractions, FDR and plots as exploratory association only; no mouse-level DE claim |
| External DE-list comparisons | Directional overlap and log-fold-change correlation plots | RETAIN_WITH_CAVEAT | Preserve as external concordance checks with source paths and gene-overlap denominators |
| Spatial state transfer | State metadata are transferred through interactive state; a loop iterates colour values rather than state-column names | CORRECT | Use one tested metadata join by cell ID and assert exact non-missing Eos count |
| Eos/reference spatial pools | Current code selects `Final_CellType_subtype_refined`, but lacks an explicit overlap/self-match gate and retains a stale invalid-legacy warning | CORRECT | Build the reference as the exact cell-ID complement of inclusive Eos and stop on overlap, self-match or distance <= 0 |
| Closest-cell and per-type distances | Closest type, distance-by-type and continuous balance correlations displayed | RETAIN_WITH_CAVEAT | Recompute after hard gates; retain effects and plots with spatially constrained null tests |
| Fifteen-cell neighbourhood composition | State-group differences and continuous correlations displayed | RETAIN_WITH_CAVEAT | Recompute from validated edges; retain descriptive summaries and add constrained permutations |
| Suspicious-gene proximity check | Vwf/Eng proximity to endothelial cells inspected | RETAIN_WITH_CAVEAT | Retain as contamination/proximity diagnostic, not proof of transcript origin |
| Ligand–receptor catalogue filtering | 21 simple interactions retained from the mouse CellChat catalogue | RETAIN | Retain panel-coverage inventory |
| Ligand–receptor product scores | Products of log-normalized ligand and receptor expression across edges | CORRECT | Rename spatial co-detection/co-expression, report zero inflation and use constrained nulls; do not call CellChat inference |
| Final Eos spatial maps | Continuous balance and q10/q90 extreme maps displayed | RETAIN | Recompute from ID-aligned metadata and retain both maps |
| Notebook execution order | Several cells have missing/nonmonotonic execution counts; `_spatial.rds` is referenced although `_spatial_passQC.rds` is saved | CORRECT | B2 must execute top-to-bottom with no hidden state, stale reload names or empty cells |

## Release interpretation

The B1 notebook is valuable as an exploratory execution record, but it is not a reproducible analysis endpoint. B2 must reproduce all retained evidence and must replace every `CORRECT` row before Region 3 results are used for downstream decisions.
