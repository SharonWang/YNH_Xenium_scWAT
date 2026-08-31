# Region 3 complete-479 stepwise review notebook design

**Date:** 2026-09-01  
**Status:** Awaiting user review  
**Scope:** Region 3 only; complete 479-gene Xenium panel

## Objective

Create a new, reader-facing Region 3 notebook that preserves the step-by-step style of `B1_Region3_primary_479.ipynb` while removing hidden interactive state, correcting invalid spatial logic, making statistical limitations explicit, and moving genuinely reusable operations into `R/source.R`.

The existing B1 notebook and its executed outputs are retained unchanged as provenance. The proposed notebook is:

`notebooks/B2_Region3_complete479_stepwise_reviewed.ipynb`

## Fixed analytical decisions

1. Use all 479 Gene Expression features for normalization, scaling, PCA, graph construction, UMAP and clustering.
2. Use PCs 1–30 for all downstream graph, UMAP and transcriptomic-neighbour operations. Object names and captions must say `PC1_30` rather than `PC12`.
3. Use only `primary_include_revised` for the analysed cell cohort. It is defined after cell-ID alignment as:

   `qc_core_pass & !high_control_flag & !segmentation_multiplet_flag`

4. Attach masks by exact `cell_id` matching. Row-order assignment is prohibited. The notebook must stop on duplicated IDs, unmatched cells, missing mask values or failed object/metadata alignment.
5. The inclusive Eosinophil set is the union of:

   - cells already annotated as Eosinophil; and
   - cells called `REF_EOS_TIER1` or `REF_EOS_TIER2` by the existing reference-supported tier logic.

6. Preserve the original annotation and add `Eos_origin` with mutually exclusive provenance categories such as `ANNOTATION_EOS`, `TIER1_RESCUE`, `TIER2_RESCUE`, and `NOT_EOS`. The inclusive relabel must never erase the original label or tier evidence.
7. Region 3 is one section from one mouse. Cell-level tests are exploratory and cannot establish mouse-level biological differential expression.

## Notebook architecture

### 0. Reader guide and reproducibility contract

State the scientific purpose, complete-panel scope, interpretation limits, required inputs and generated outputs. Record project paths, run label, R/kernel version, package versions, random seed, source-file hash, notebook hash when available, input-file hashes and checkpoint paths. Assert that every write target is below the D:-mapped project root locally or the declared HPC project root remotely.

### 1. Import and input validation

Discover Region 3, validate the required Xenium files, import the Gene Expression matrix and spatial geometry, and assert exactly 479 genes. Load `cell_qc_metadata.tsv.gz`, match it to Seurat cells by `cell_id`, and report match coverage before adding metadata.

### 2. QC review and analysed cohort

Display input-cell counts, core-QC counts, individual review flags and spatial flag maps before subsetting. Compute `primary_include_revised` only from aligned Seurat metadata. Show the exact inclusion cross-tabulation and save an unfiltered checkpoint. Subset once to `primary_include_revised == TRUE` and save the analysed-cohort checkpoint without modifying the raw object.

### 3. Complete-panel normalization and PCA

Use Seurat `LogNormalize` with scale factor 10,000, `ScaleData` on all 479 genes and `RunPCA(..., npcs = 30)`. Display explained variance, elbow plot, PC loadings and Spearman correlations between PCs and QC metrics. These are diagnostics only: PCs 1–30 remain fixed for downstream analysis as requested.

### 4. Graph, UMAP and clustering stability

Build neighbours and UMAP from PCA dimensions 1–30. Run Leiden clustering at the requested resolution grid (0.7, 0.8, 1.0 and 1.2) using fixed seeds. Report cluster counts, minimum/median cluster sizes, adjacent-resolution adjusted Rand indices and an alluvial/clustree-style transition view when the optional package is available. Resolution 0.8 is the prespecified working candidate, but the notebook must visibly flag unsupported, tiny or marker-incoherent clusters.

### 5. Marker-based annotation

Score the curated scWAT marker groups at cluster level and show marker coverage, top/second scores and score margins. Preserve marker evidence separately from final labels. Generate UMAP, dot/heatmap and spatial views for the proposed labels.

### 6. Wang-reference annotation

Use the age-nearest 2.5-month Wang reference as the principal transfer and the complete Wang reference as a concordance check. Transfer can only use genes shared by query and reference; the notebook must report the shared count rather than implying that all 479 are present in the reference. Retain prediction scores, score margins and competing labels. Assign `Uncertain` when reference confidence is insufficient or marker and reference evidence conflict.

### 7. Consolidated annotation

Combine cluster-marker and Wang-reference evidence without silently forcing uncertain cells. Report label counts, confidence categories, disagreement tables and spatial distributions. The consolidation function must preserve all contributing evidence columns.

### 8. Inclusive Eosinophil definition

Run the existing tier-calibration logic, show reference marker performance and cross-tabulate Tier 1/2 calls against original cell types. Construct the inclusive Eos set exactly as requested: original annotation Eos OR Tier 1/2. Store `Eos_origin`, the original subtype, tier, tier rule and reference scores. Explicitly show how many cells enter through each route before creating the Eos object.

### 9. Eosinophil state analysis

Follow the existing short-/long-lived scoring logic, including the complete 47 short-lived and 46 long-lived panels and the non-ribosomal short-lived score. Show gene availability, detection rates, score distributions, complexity correlations, dip test and one- to three-component mixture-model comparison. Quantile-defined extreme groups remain descriptive labels. State prominently that unimodality or a one-component mixture does not support discrete Eos subtypes.

Any cell-level long-like versus short-like expression comparison must be labelled exploratory. Report effect sizes, detection fractions and adjusted p-values, but do not call it biological DE. Signature genes used to define the state cannot be presented as independent validation of that state; comparisons to external lists are concordance checks only.

### 10. Valid spatial-neighbour analysis

Define query cells from the inclusive Eos set and reference cells from its exact complement using cell IDs, not labels alone. Assert:

- no query/reference cell-ID overlap;
- no duplicated cell IDs;
- finite coordinates;
- positive query-to-reference nearest distances; and
- no self-matches.

Recalculate nearest cell type, nearest distance by cell type and 15-neighbour composition. Replace naive cell-level significance claims with spatially constrained null tests that permute Eos state labels or continuous balance values within prespecified spatial blocks and local-density strata using a fixed seed. Report observed effect, null interval, empirical p-value, multiplicity correction and the number of valid permutations. Keep conventional Wilcoxon or Spearman results only as descriptive diagnostics if retained.

### 11. Exploratory ligand–receptor co-expression

Do not call the calculation CellChat inference. Filter the mouse interaction database to panel-observed ligand and receptor genes, then calculate explicitly named spatial co-detection/co-expression scores across validated Eos-neighbour edges. Use the same spatially constrained null framework for state contrasts. Report panel coverage and zero-inflation, and label all outputs hypothesis-generating.

### 12. Outputs and handoff

Save checkpoints and compact TSV summaries beneath a new Region 3 output root, including provenance, analysed-cell IDs, clustering stability, annotation evidence, Eos provenance, Eos-state summaries, validated neighbour results and exploratory ligand–receptor scores. Finish with a markdown interpretation section that separates verified computation, exploratory evidence, limitations and HPC rerun requirements.

## Reusable functions for `R/source.R`

The notebook will call documented, testable helpers for:

- provenance and file-hash capture;
- exact mask-to-Seurat alignment and revised-mask construction;
- complete-panel PCA diagnostics;
- clustering stability summaries;
- Wang-reference preparation and transfer;
- annotation evidence consolidation with `Uncertain` handling;
- inclusive Eos provenance construction;
- Eos-state scoring and diagnostics;
- disjoint spatial query/reference construction;
- nearest-distance and neighbourhood composition calculations;
- spatially constrained permutations; and
- exploratory ligand–receptor co-expression scoring.

Notebook-specific plotting choices and biological interpretation remain in the notebook. Existing public function interfaces will not be changed unless required for correctness; new helpers receive roxygen-style purpose, parameter, return-value and failure-condition documentation.

## Validation strategy

1. Preserve all current uncommitted HPC-returned files.
2. Add contract tests first for cell-ID alignment, the exact revised mask, 479-gene enforcement, PC1–30 usage, Eos union/provenance and disjoint spatial pools.
3. Run R parse/source tests with `D:\Programs\R-4.6.1` and `D:\Programs\R_library`, with all temporary paths forced under the D: project.
4. Execute a bounded local subset when dependencies and memory permit.
5. Validate notebook JSON, R-cell syntax, top-to-bottom dependency order and absence of stored runtime errors.
6. Provide exact HPC execution commands for the full-data notebook. Full-data conclusions remain unvalidated until the executed notebook and outputs are returned.

## Explicit exclusions

- Do not modify or delete `B1_Region3_primary_479.ipynb`.
- Do not use Harmony for the single Region 3 analysis.
- Do not infer mouse-level biology from one section.
- Do not describe product-of-expression scores as CellChat probabilities.
- Do not write project, temporary, cache or test artifacts to C:.
