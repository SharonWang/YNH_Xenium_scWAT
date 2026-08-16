# scWAT Region 3 Anchor Reference Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a sequential downstream Xenium pipeline that freezes a Region 3 anchor, independently gates Regions 1-2 for admission, falls back to Region 3 when an expanded consensus fails, maps Region 4 without reference influence, and releases Eosinophil conclusions only when Region 3-centered sensitivity checks are stable.

**Architecture:** Phase 0-2 downstream bundles remain immutable inputs. Seurat constructs the Region 3 anchor and performs reference mapping with `FindTransferAnchors`, `TransferData`, and `MapQuery`; Harmony builds an expanded eligible-section embedding while holding Region 3 as the reference level. Jupyter notebooks orchestrate one gate at a time, retain unintegrated PCA diagnostics, and fall back to the frozen Region 3 anchor when the expanded Harmony consensus is unstable or section-driven.

**Tech Stack:** R 4.3+, Seurat 4.3+ or 5.x, SeuratObject, harmony 1.2+, leidenbase, Matrix, clue, base R TSV/gzip/RDS I/O, ggplot2, dplyr, jsonlite, Python 3 notebook generation/contract validation, Jupyter IRkernel, PowerShell local contract runner, Bash/Slurm HPC runner.

## Global Constraints

- Write every implementation file, build artifact, test output, cache, and temporary file below `D:/Xiaonan/CODEX_projects/Yanan_Xenium` locally; write nothing to C:.
- Use `/dssg/home/acct-svetoslav_chakarov/svetoslav_chakarov/Lab_members/Yanan_Hu/YNH_Xenium` as the HPC project root.
- Run only bounded cell subsets locally; run full matrices and final gates on HPC.
- Preserve raw sparse counts, all cells, all masks, and all 479 genes. Add derived objects and decisions without modifying raw objects.
- Use Region 3 as the immutable initial anchor. Regions 1 and 2 are evaluated independently against that same anchor.
- Build an expanded consensus from Region 3 plus every admitted Region 1-2 section. If the expanded consensus fails, release the frozen Region 3 anchor instead.
- Region 4 is always mapping-only and sensitivity-only. It must not affect normalization, feature selection, PCA, integration, graph construction, clustering, markers, thresholds, centroids, or reference labels.
- Primary PCA and clustering use the 245 `PROVISIONAL_PRIMARY_FEATURES`; the 67-gene set is a strict sensitivity; the 234 technical-risk genes never define primary clusters; all 479 genes remain available for exploratory inspection.
- Eosinophil primary evidence comes from Region 3. The frozen eligible reference provides validation; Region 4 is a separately labelled sensitivity result that cannot establish, reverse, or rescue the primary conclusion.
- Mouse is the biological unit and section is the technical unit. Do not use cells or sections as biological replicates for inference.
- Use seed `20260814L` for every stochastic operation and record it in every model and output.
- Use failing-first red/green tests for every reusable function, artifact schema, notebook contract, and runner branch.
- Use Seurat `LogNormalize` with a fixed scale factor of 10,000. Do not run data-driven variable-feature selection for primary analysis; pass the frozen 245-gene vector explicitly to `ScaleData`, `RunPCA`, transfer-anchor finding, and primary marker testing.
- Use Harmony only after Region 1-2 admission, only on eligible reference cells, with `region_id` as the integration variable and Region 3 protected through `reference_values="Region_3"`. Retain the unintegrated PCA as a required diagnostic.
- Do not harmonize, regress, or integrate on `mouse_id`; with two mice and section/mouse confounding, mouse effects must remain visible for validation rather than be forcibly removed.
- Local absence of Seurat/Harmony is an expected environment limitation. Do not replace these methods locally; run pure-R contracts and notebook structure locally and run every model-fitting red/green checkpoint on HPC.

## File structure

**Create:**

- `config/downstream_reference_defaults.tsv` — frozen model and gate parameters.
- `config/scwat_canonical_markers.tsv` — panel-aware positive marker sets and marker-use policy.
- `tests/test_downstream_reference.R` — pure-R contracts plus conditional Seurat/Harmony model, mapping, admission, fallback, Region 4 isolation, Eosinophil, and artifact tests.
- `scripts/build_downstream_subset_fixture.R` — bounded four-region fixture builder that preserves the full 479-gene contract.
- `scripts/execute_downstream_subset.ps1` — local, one-gate-at-a-time executor with D:-only temporary paths.
- `notebooks/03_region3_anchor_reference.ipynb` — anchor construction and internal sensitivities.
- `notebooks/04_region1_anchor_validation.ipynb` — independent Region 1 admission.
- `notebooks/05_region2_anchor_validation.ipynb` — independent Region 2 admission.
- `notebooks/06_region1_3_consensus_reference.ipynb` — eligible consensus and Region 3 fallback.
- `notebooks/07_region4_consensus_mapping.ipynb` — mapping-only Region 4 sensitivity.
- `notebooks/08_region3_eosinophil_robustness.ipynb` — Region 3-primary Eosinophil stability.
- `notebooks/09_downstream_release_summary.ipynb` — final reader-facing release report.
- `shell/run_downstream_reference_hpc.sh` — gated HPC launcher.
- `slurm/scwat_downstream_reference.sbatch` — one-target Slurm wrapper.

**Modify:**

- `R/source.R` — reusable downstream functions, appended after existing Phase 0-2 functions.
- `scripts/render_notebooks.py` — downstream notebook builders, parameter injection, and structural validation.
- `tests/test_notebook_contracts.py` — exact notebook names, order, parameters, inputs, outputs, and prohibited Region 4 behavior.
- `notebooks/02_slide_QC_summary.ipynb` — explicit handoff to Notebook 03 after completed evidence-only gates pass.
- `README.md` — local/HPC prerequisites, gate order, inputs, outputs, and commands.
- `D:/Xiaonan/CODEX_projects/Yanan_Xenium/PROJECT_RESEARCH_PLAN_AND_PROGRESS.md` — checkpoint results and changes after each verified task.

---

### Task 1: Freeze downstream configuration and canonical-marker contracts

**Files:**
- Create: `config/downstream_reference_defaults.tsv`
- Create: `config/scwat_canonical_markers.tsv`
- Create: `tests/test_downstream_reference.R`
- Modify: `R/source.R`

**Interfaces:**
- Consumes: versioned TSV configuration files.
- Produces: `read_downstream_reference_config(path) -> named list`; `read_canonical_marker_config(path, panel_genes, gene_decision) -> data.frame`.

- [ ] **Step 1: Write failing parser and policy tests**

```r
source(file.path(repo_root, "R", "source.R"))
cfg <- read_downstream_reference_config(file.path(repo_root, "config", "downstream_reference_defaults.tsv"))
stopifnot(cfg$seed == 20260814L, cfg$primary_gene_count == 245L)
stopifnot(cfg$conservative_gene_count == 67L, cfg$raw_gene_count == 479L)
stopifnot(cfg$n_pcs == 30L, cfg$knn_k == 30L, cfg$mapping_folds == 5L)
stopifnot(cfg$harmony_reference_region == "Region_3", cfg$harmony_theta == 2)
markers <- read_canonical_marker_config(marker_path, panel_genes, gene_decision)
stopifnot(all(markers$gene %in% panel_genes))
stopifnot(all(markers$use_policy %in% c("PRIMARY_SUPPORT", "VALIDATION_ONLY_TECHNICAL_RISK")))
stopifnot(all(markers$use_policy[markers$gene %in% risk_genes] == "VALIDATION_ONLY_TECHNICAL_RISK"))
```

- [ ] **Step 2: Run the new test and verify RED**

Run locally with all temporary variables set to `D:/Xiaonan/CODEX_projects/Yanan_Xenium/adipose_analysis/tmp`:

```powershell
& 'C:\Program Files\R\R-4.3.3\bin\Rscript.exe' tests/test_downstream_reference.R
```

Expected: FAIL because the configuration files and readers do not exist.

- [ ] **Step 3: Add the exact default configuration**

Write `key`, `value`, and `description` columns with at least:

```text
seed=20260814
normalization_scale_factor=10000
primary_gene_count=245
conservative_gene_count=67
technical_risk_gene_count=234
raw_gene_count=479
n_pcs=30
knn_k=30
leiden_resolution=0.6
seurat_min_version=4.3.0
harmony_min_version=1.2.0
harmony_theta=2
harmony_lambda=1
harmony_sigma=0.1
harmony_max_iter=20
harmony_reference_region=Region_3
mapping_folds=5
mapping_k=30
mapping_confidence_quantile=0.05
mapping_distance_quantile=0.95
mapping_min_label_cells=200
major_label_fraction=0.01
major_label_cells=500
section_cluster_dominance=0.80
section_cluster_min_cells=200
label_stability_min=0.85
gene_sensitivity_concordance_min=0.75
section_predictability_permutations=199
section_predictability_margin=0.05
eos_assignment_jaccard_min=0.80
eos_score_spearman_min=0.70
eos_module_nbin=12
eos_module_ctrl=5
```

The parser must reject duplicate keys, missing required keys, non-finite numeric values, counts below one, quantiles outside `(0,1)`, and proportions outside `[0,1]`.

- [ ] **Step 4: Add panel-aware marker configuration**

Use columns `cell_type`, `gene`, `direction`, `marker_group`, and `use_policy`. Include only genes observed in the 479-gene panel. Initial positive marker groups must cover, where available: adipocyte (`Aqp7`, `Pck1`, `Retn`, `Cox8b`), endothelial (`Kdr`, `Vwf`, `Gng11`, `Plvap`, `Eng`, `Sox17`, `Tie1`, `Clec14a`, `Gpihbp1`, `Aqp1`), lymphatic endothelial (`Lyve1`, `Prox1`, `Mmrn1`), stromal/fibroblast (`Pi16`, `Dpt`, `Lum`, `Mfap4`, `Mfap5`, `Pcolce2`, `Sfrp1`, `Bgn`, `Aspn`, `Serpinf1`, `Ccdc80`), mural (`Rgs5`, `Cnn1`, `Myh11`, `Tagln`, `Rbp1`, `Itga8`, `Higd1b`, `Ndufa4l2`), macrophage (`Adgre1`, `C5ar1`, `F13a1`, `Folr2`, `Marco`, `Mpeg1`, `Ms4a7`, `Pf4`, `Lst1`, `Ccl6`, `Cd5l`), neutrophil (`Csf3r`, `Cxcr2`, `Mmp8`, `Mmp9`, `Mpo`, `Pglyrp1`, `Ctsg`), Eosinophil (`Siglecf`, `Prg2`, `Ccr3`, `Il5ra`, `Alox15`, `Ear1`, `Ear2`, `Ltc4s`), T cell (`Cd3d`, `Cd8a`, `Ctla4`), mast cell (`Cpa3`, `Mrgpra2a`, `Il1rl1`), and neural/Schwann (`Plp1`, `Pmp22`, `Prx`, `Gfap`).

At read time, assign `VALIDATION_ONLY_TECHNICAL_RISK` to every marker in the 234-gene risk tier; otherwise assign `PRIMARY_SUPPORT`. A cell type with fewer than two available positive markers is `INSUFFICIENT_PANEL_SUPPORT` and cannot be auto-released.

- [ ] **Step 5: Implement readers, rerun the test, and verify GREEN**

- [ ] **Step 6: Commit**

```powershell
git add config/downstream_reference_defaults.tsv config/scwat_canonical_markers.tsv R/source.R tests/test_downstream_reference.R
git commit -m "feat: freeze downstream reference configuration"
```

### Task 2: Validate immutable Phase 0-2 handoff and create a bounded fixture

**Files:**
- Create: `scripts/build_downstream_subset_fixture.R`
- Modify: `R/source.R`
- Modify: `tests/test_downstream_reference.R`

**Interfaces:**
- Consumes: four Phase 0-2 downstream RDS bundles and the six slide-level decision tables.
- Produces: `validate_downstream_handoff(qc_run_root, stop_on_error=TRUE) -> data.frame`; `read_downstream_handoff(qc_run_root) -> list`; bounded fixture below `adipose_analysis/subset_input/downstream_reference/`.

- [ ] **Step 1: Write failing handoff tests**

```r
handoff <- read_downstream_handoff(fixture_root)
stopifnot(identical(names(handoff$regions), paste0("Region_", 1:4)))
stopifnot(nrow(handoff$gene_decision) == 479L)
stopifnot(sum(handoff$gene_decision$primary_feature_status == "PROVISIONAL_PRIMARY_FEATURES") == 245L)
stopifnot(sum(handoff$gene_decision$conservative_evidence_status == "CONSERVATIVE_NO_SIGNAL_DETECTED") == 67L)
stopifnot(sum(handoff$gene_decision$technical_risk_status == "TECHNICAL_RISK_SENSITIVITY_ONLY") == 234L)
stopifnot(all(vapply(handoff$regions, function(x) inherits(x$counts, "sparseMatrix"), logical(1))))
stopifnot(validate_downstream_handoff(fixture_root, stop_on_error=TRUE)$status == "PASS")
```

Add negative fixtures for duplicated cell IDs, altered count dimensions, missing masks, a non-PASS completed QC gate, and Region 4 incorrectly marked reference-eligible.

- [ ] **Step 2: Run and verify RED**

Expected: FAIL because the handoff functions do not exist.

- [ ] **Step 3: Implement strict readers and validators**

`validate_downstream_handoff` must confirm unique `region_id + cell_id`, exact matrix/cell alignment, raw sparse counts, 479/245/67/234 and 53 Eosinophil counts, section statuses, nested masks, PASS for all completed Phase 0-2 gates, and the Region 4 mapping-only contract. It must return a table of checks and stop before model fitting on any failure.

- [ ] **Step 4: Implement the bounded fixture builder**

The script takes `--source-qc-run`, `--output-root`, `--cells-per-region`, and `--seed`. It samples at most 500 cells per region, preserves all 479 genes and the audited slide-level gene/Eosinophil decisions, subsets every matrix and metadata table in identical cell order, writes hashes/provenance, and refuses an output path outside the project root.

Expected local target:

```text
D:/Xiaonan/CODEX_projects/Yanan_Xenium/adipose_analysis/subset_input/downstream_reference/
```

- [ ] **Step 5: Build the fixture, validate it, and verify GREEN**

- [ ] **Step 6: Commit code and tests; do not commit generated fixture data**

```powershell
git add R/source.R tests/test_downstream_reference.R scripts/build_downstream_subset_fixture.R
git commit -m "feat: validate downstream QC handoff"
```

### Task 3: Build the Seurat Region 3 anchor and sensitivity branches

**Files:**
- Modify: `R/source.R`
- Modify: `tests/test_downstream_reference.R`

**Interfaces:**
- Consumes: one sparse counts matrix, aligned cell metadata, fixed gene vector, mask, model configuration, marker configuration.
- Produces: `build_seurat_reference(counts, cells, genes, config, role, seed) -> Seurat`; `find_primary_markers(object, genes, config) -> data.frame`; `match_cluster_labels(reference_labels, candidate_labels, shared_cell_ids) -> data.frame`; `validate_anchor_branches(primary, strict, hotspot, conservative, markers, config) -> list`.

- [ ] **Step 1: Write failing environment and static contract tests locally**

```r
deps <- downstream_model_dependencies()
stopifnot(identical(deps$package, c("Seurat", "SeuratObject", "harmony")))
stopifnot(deps$minimum_version[deps$package == "Seurat"] == "4.3.0")
stopifnot(deps$minimum_version[deps$package == "harmony"] == "1.2.0")
stopifnot(identical(primary_model_feature_policy(), "FIXED_245_NO_VARIABLE_FEATURE_SELECTION"))
```

If Seurat/Harmony are absent locally, model execution tests must print `SKIP_LOCAL_MODEL_TEST_HPC_REQUIRED` and exit successfully after all pure-R contract tests. They must not substitute another clustering implementation.

- [ ] **Step 2: Run locally and verify RED for missing functions**

- [ ] **Step 3: Implement the Seurat anchor wrapper**

The wrapper must execute this fixed workflow:

```r
object <- Seurat::CreateSeuratObject(counts = counts[, cells$cell_id, drop=FALSE], meta.data = cells)
object <- Seurat::NormalizeData(object, normalization.method="LogNormalize", scale.factor=10000, verbose=FALSE)
object <- Seurat::ScaleData(object, features=genes, verbose=FALSE)
npcs <- min(config$n_pcs, length(genes)-1L, ncol(object)-1L)
object <- Seurat::RunPCA(object, features=genes, npcs=npcs, seed.use=seed, verbose=FALSE)
object <- Seurat::FindNeighbors(object, reduction="pca", dims=seq_len(npcs), k.param=config$knn_k, verbose=FALSE)
object <- Seurat::FindClusters(object, resolution=config$leiden_resolution, algorithm=4, random.seed=seed, verbose=FALSE)
object <- Seurat::RunUMAP(object, reduction="pca", dims=seq_len(npcs), seed.use=seed, return.model=TRUE, verbose=FALSE)
```

Do not call `FindVariableFeatures`. Store the exact 245/67 feature vector, mask, seed, package versions, and raw-count hash. Primary cluster-defining markers use `FindAllMarkers(features=primary_genes, assay="RNA", slot="data")`; technical-risk markers are summarized separately and cannot define a primary label.

- [ ] **Step 4: Implement four Region 3 branches**

Fit primary+245, strict+245, hotspot-sensitivity+245, and primary+67 Seurat objects. Match clusters with `clue::solve_LSAP` on shared-cell Jaccard overlap. Calculate shared-cell label agreement, canonical-marker support, hotspot dependence, 245-versus-67 compatibility, and technical-risk validation dependence.

- [ ] **Step 5: Write the failing HPC model test before implementation is accepted**

Run the bounded downstream fixture on HPC and assert deterministic cell order, PCA dimension, cluster labels under the fixed seed, absence of risk-only genes from PCA/primary markers, and `PASS_ANCHOR`/`REVIEW_ANCHOR`/`STOP_ANCHOR` behavior on synthetic perturbations.

- [ ] **Step 6: Run twice on HPC and verify GREEN**

- [ ] **Step 7: Commit**

```powershell
git add R/source.R tests/test_downstream_reference.R
git commit -m "feat: build and validate Seurat Region 3 anchor"
```

### Task 4: Map and independently gate Regions 1-2 with Seurat anchors

**Files:**
- Modify: `R/source.R`
- Modify: `tests/test_downstream_reference.R`

**Interfaces:**
- Consumes: frozen Region 3 Seurat anchor/labels, one candidate bundle, masks, 245/67 genes, markers, configuration.
- Produces: `calibrate_seurat_mapping(anchor, labels, config) -> list`; `map_seurat_query(reference, query_bundle, genes, calibration, config) -> list`; `evaluate_candidate_admission(anchor, calibration, candidate_bundle, gene_sets, markers, config) -> list`.

- [ ] **Step 1: Write failing pure-R threshold and state tests locally**

Test calibration quantiles, per-label/global fallback, `Uncertain` assignment, admission-table schemas, and proof that Region 1 output is independent of Region 2 state.

- [ ] **Step 2: Write failing Seurat mapping tests for HPC**

```r
cal <- calibrate_seurat_mapping(anchor, anchor$reference_label, cfg)
mapped <- map_seurat_query(anchor, candidate_bundle, primary_genes, cal, cfg)
required <- c("cell_id","predicted_id","prediction_score_max","reference_distance",
              "second_label","confidence_margin","final_label","threshold_source")
stopifnot(all(required %in% names(mapped$cell_mapping)))
stopifnot(all(mapped$cell_mapping$final_label[
  mapped$cell_mapping$prediction_score_max < mapped$cell_mapping$confidence_cutoff |
  mapped$cell_mapping$reference_distance > mapped$cell_mapping$distance_cutoff] == "Uncertain"))
```

- [ ] **Step 3: Implement held-out Region 3 calibration with Seurat label transfer**

Use five stratified folds. For each held-out fold, build transfer anchors from training Region 3 cells with:

```r
anchors <- Seurat::FindTransferAnchors(
  reference=reference, query=query, normalization.method="LogNormalize",
  reduction="pcaproject", reference.reduction="pca", features=genes,
  dims=seq_len(npcs), k.anchor=5, k.score=30, verbose=FALSE
)
prediction <- Seurat::TransferData(anchorset=anchors, refdata=reference$reference_label,
                                   dims=seq_len(npcs), verbose=FALSE)
```

Use `prediction.score.max` as confidence. Compute reference distance from the query's projected PCA coordinates to the assigned-label centroid. Use the 5th percentile of correct confidence and 95th percentile of correct distance globally and per label when at least 200 calibration cells exist.

- [ ] **Step 4: Implement candidate mapping with `MapQuery`**

Run `FindTransferAnchors` with the frozen 245 features and anchor PCA, then `MapQuery` using the anchor PCA and stored UMAP model. Region 1/2 never alter the anchor. Repeat mapping with strict cells and the 67-gene anchor branch for sensitivity.

- [ ] **Step 5: Implement admission diagnostics**

Calculate uncertainty, primary-versus-strict agreement, 245-versus-67 agreement, within-label pseudobulk Spearman correlation, canonical-marker support, and candidate-enriched clusters in anchor-projected coordinates. A diagnostic cluster is section-dominated at 80% candidate and at least 200 candidate cells; it fails admission only when dominance accompanies unsupported markers, mapping instability, conservative contradiction, or technical-risk dependence.

- [ ] **Step 6: Run the bounded HPC mapping tests and verify GREEN**

- [ ] **Step 7: Commit**

```powershell
git add R/source.R tests/test_downstream_reference.R
git commit -m "feat: gate candidate sections with Seurat mapping"
```

### Task 5: Build the eligible Seurat-Harmony consensus and Region 3 fallback

**Files:**
- Modify: `R/source.R`
- Modify: `tests/test_downstream_reference.R`

**Interfaces:**
- Consumes: frozen anchor, Region 1/2 admission tables, admitted bundles, configuration.
- Produces: `select_eligible_reference_regions(admission_tables) -> character`; `build_harmony_consensus(anchor, admitted_bundles, genes, config) -> Seurat`; `validate_expanded_consensus(anchor, consensus, branches, markers, config) -> list`; `freeze_final_reference(anchor, consensus, consensus_gate) -> list`.

- [ ] **Step 1: Write failing pure-R state-machine tests locally**

```r
stopifnot(identical(select_eligible_reference_regions(pass_pass), c("Region_1","Region_2","Region_3")))
stopifnot(identical(select_eligible_reference_regions(pass_fail), c("Region_1","Region_3")))
stopifnot(identical(select_eligible_reference_regions(fail_pass), c("Region_2","Region_3")))
stopifnot(identical(select_eligible_reference_regions(fail_fail), "Region_3"))
fallback <- freeze_final_reference(anchor, failed_consensus, list(decision="FAIL_CONSENSUS"))
stopifnot(fallback$decision == "REGION3_ANCHOR_FALLBACK")
stopifnot(identical(fallback$reference_hash, anchor_hash))
```

- [ ] **Step 2: Write failing Harmony isolation tests for HPC**

Assert that only eligible regions are present, `Region_3` is the Harmony `reference_values`, the 245 features define PCA, both `pca` and `harmony` reductions are retained, and no Region 4 barcode appears.

- [ ] **Step 3: Implement section-balanced discovery input**

Use Region 3 plus admitted regions only. Set each region's discovery cap to the smallest eligible primary-cell count and sample with the fixed seed, stratified by preliminary anchor label. Merge raw-count Seurat objects, normalize with `LogNormalize`, scale the 245 genes, and run unintegrated PCA. Preserve this PCA for section-dominance diagnostics.

- [ ] **Step 4: Run anchored Harmony and cluster in Harmony space**

```r
consensus <- harmony::RunHarmony(
  object=consensus, group.by.vars="region_id", reduction="pca",
  dims.use=seq_len(npcs), theta=config$harmony_theta,
  lambda=config$harmony_lambda, sigma=config$harmony_sigma,
  max.iter.harmony=config$harmony_max_iter,
  reference_values="Region_3", reduction.save="harmony",
  plot_convergence=TRUE, verbose=TRUE
)
consensus <- Seurat::FindNeighbors(consensus, reduction="harmony", dims=seq_len(npcs), k.param=config$knn_k)
consensus <- Seurat::FindClusters(consensus, resolution=config$leiden_resolution, algorithm=4, random.seed=config$seed)
consensus <- Seurat::RunUMAP(consensus, reduction="harmony", dims=seq_len(npcs), seed.use=config$seed, return.model=TRUE)
```

Retain both a PCA UMAP model (`umap_pca`) and a Harmony UMAP (`umap_harmony`). Use `FindAllMarkers(features=primary_genes)` for primary consensus markers, calling a version-aware `JoinLayers` helper first when Seurat v5 layers require it. Map non-discovery eligible cells back with Seurat transfer anchors; do not refit discovery clusters on those cells.

- [ ] **Step 5: Implement consensus gates and fallback**

Compare unintegrated PCA and Harmony section association, Region 3 representation/markers, primary-versus-strict stability, 245-versus-67 agreement, hotspot sensitivity, and technical-risk dependence. Section predictability uses cross-validation plus 199 label permutations. Expanded failure returns the original Region 3 Seurat object and hash unchanged as `REGION3_ANCHOR_FALLBACK`; the failed Harmony object is audit-only.

- [ ] **Step 6: Run bounded HPC tests for all four admission combinations and verify GREEN**

- [ ] **Step 7: Commit**

```powershell
git add R/source.R tests/test_downstream_reference.R
git commit -m "feat: add Seurat Harmony consensus fallback"
```

### Task 6: Map Region 4 with frozen Seurat label transfer only

**Files:**
- Modify: `R/source.R`
- Modify: `tests/test_downstream_reference.R`

**Interfaces:**
- Consumes: final frozen Seurat reference, frozen calibration, Region 4 bundle.
- Produces: `map_region4_sensitivity(final_reference, calibration, region4_bundle, config) -> list`.

- [ ] **Step 1: Write failing isolation contracts locally and model tests on HPC**

```r
before <- hash_seurat_reference(final_reference)
r4 <- map_region4_sensitivity(final_reference, calibration, region4_bundle, cfg)
stopifnot(hash_seurat_reference(final_reference) == before)
stopifnot(r4$reference_hash == before)
stopifnot(all(r4$cell_mapping$analysis_role == "SENSITIVITY_ONLY"))
stopifnot(!any(names(r4) %in% c("refit_pca","refit_harmony","refit_clusters","updated_labels")))
```

- [ ] **Step 2: Implement Seurat mapping-only behavior**

Use `FindTransferAnchors(reduction="pcaproject", reference.reduction="pca")` and `TransferData` against the frozen final reference and its Harmony-derived labels. Use `MapQuery` only with the stored PCA-based UMAP model (`umap_pca`); do not project Region 4 directly into the Harmony UMAP because Harmony was fit only on eligible reference sections. Do not call `NormalizeData`, `ScaleData`, `RunPCA`, `RunHarmony`, `FindNeighbors`, `FindClusters`, `FindAllMarkers`, or reference-mutating assignment on the reference object after its hash is recorded. Query normalization occurs only inside the new Region 4 query object.

- [ ] **Step 3: Apply frozen confidence/distance thresholds**

Report predicted label, prediction score, second label, margin, projected reference distance, threshold source, and `Uncertain`. Return `REGION4_MAPPING_NOT_RELEASABLE` when overall or lineage-specific mapping fails; never change primary-reference status.

- [ ] **Step 4: Run bounded HPC mapping twice and verify the reference hash is unchanged**

- [ ] **Step 5: Commit**

```powershell
git add R/source.R tests/test_downstream_reference.R
git commit -m "feat: isolate Region 4 Seurat mapping"
```

### Task 7: Implement Region 3-centered Seurat Eosinophil robustness and release aggregation

**Files:**
- Modify: `R/source.R`
- Modify: `tests/test_downstream_reference.R`

**Interfaces:**
- Consumes: Region 3 Seurat object/masks/hotspots, 53- and 100-gene Eosinophil sets, Region 3 anchor labels, final frozen labels, admitted-region mappings, optional Region 4 mapping.
- Produces: `add_eosinophil_scores(object, gene_sets, seed) -> Seurat`; `evaluate_eosinophil_robustness(...) -> list`; `aggregate_downstream_release(...) -> data.frame`.

- [ ] **Step 1: Write failing evidence-hierarchy tests locally**

Test that Region 3, eligible-reference, and Region 4 metrics enter separate reducers and that Eosinophil or Region 4 failure does not stop a valid reference.

- [ ] **Step 2: Implement Seurat program scoring on HPC**

Use `Seurat::AddModuleScore` on the `RNA` assay's normalized data for the 53-gene provisional set, 100-gene complete set, and common/short-lived/long-lived components, with `nbin=12`, `ctrl=5`, and the fixed seed because the panel contains only 479 genes. Record available-gene denominators and selected controls. Because the panel is targeted, also calculate the uncorrected mean normalized expression of each set as a sensitivity; conclusions must agree in direction between `AddModuleScore` and the mean-expression check.

- [ ] **Step 3: Implement Region 3 primary comparisons**

Compare 53 versus 100 genes, primary versus strict masks, hotspot retained versus excluded, and anchor versus frozen-reference assignments. Use Spearman score correlation and Jaccard overlap of Eosinophil assignments with minima 0.70 and 0.80; any direction reversal produces `EOS_PRIMARY_NOT_STABLE`.

- [ ] **Step 4: Add eligible-reference validation and separate Region 4 sensitivity**

Admitted regions validate direction and mapped assignments but cannot override Region 3. Region 4 is scored only among confidently mapped cells, written separately, and excluded from the primary reducer.

- [ ] **Step 5: Implement orthogonal release aggregation**

Only Region 3 anchor failure produces `STOP_REGION3_ANCHOR_FAILED`. Expanded Harmony failure maps to Region 3 fallback. Eosinophil instability and Region 4 mapping failure remain branch-specific.

- [ ] **Step 6: Run pure-R tests locally and Seurat score tests on the bounded HPC fixture; verify GREEN**

- [ ] **Step 7: Commit**

```powershell
git add R/source.R tests/test_downstream_reference.R
git commit -m "feat: gate Seurat Eosinophil analysis"
```

### Task 8: Add provenance-aware downstream artifacts and reload validators

**Files:**
- Modify: `R/source.R`
- Modify: `tests/test_downstream_reference.R`

**Interfaces:**
- Consumes: outputs from Tasks 3-7.
- Produces: `write_downstream_stage_artifacts(stage, result, output_root, project_root) -> character`; `validate_downstream_stage_artifacts(stage, output_root, stop_on_error=TRUE) -> logical`.

- [ ] **Step 1: Write failing artifact tests**

For each stage, write a fixture, reload it in a fresh R session, compare row counts, hashes, cell order, feature order, decisions, thresholds, and upstream reference hash. Add corrupt-file, wrong-reference-hash, missing-provenance, and Region 4 reference-mutation failures.

- [ ] **Step 2: Run and verify RED**

- [ ] **Step 3: Implement stage schemas and writers**

Required files include all specification outputs:

```text
03_region3_anchor/region3_anchor_reference.rds
03_region3_anchor/region3_anchor_cell_labels.tsv.gz
03_region3_anchor/region3_anchor_cluster_markers.tsv.gz
03_region3_anchor/region3_anchor_cluster_stability.tsv
03_region3_anchor/region3_anchor_qc_gate.tsv
04_region1_admission/region1_anchor_mapping.tsv.gz
04_region1_admission/region1_celltype_concordance.tsv
04_region1_admission/region1_section_cluster_diagnostics.tsv
04_region1_admission/region1_admission_gate.tsv
05_region2_admission/region2_anchor_mapping.tsv.gz
05_region2_admission/region2_celltype_concordance.tsv
05_region2_admission/region2_section_cluster_diagnostics.tsv
05_region2_admission/region2_admission_gate.tsv
06_consensus/final_frozen_reference.rds
06_consensus/final_frozen_reference_labels.tsv.gz
06_consensus/final_reference_decision.tsv
07_region4_mapping/region4_consensus_mapping.tsv.gz
07_region4_mapping/region4_mapping_by_celltype.tsv
07_region4_mapping/region4_mapping_qc_gate.tsv
08_eosinophil/eos_robustness_decision.tsv
09_release/downstream_release.tsv
```

Every TSV contains `schema_version`, `qc_run_label`, `downstream_run_label`, `generated_utc`, `seed`, `cell_mask`, `gene_tier`, `threshold_source`, `source_artifact`, and `reference_hash` where applicable.

- [ ] **Step 4: Implement reload validation and verify GREEN**

- [ ] **Step 5: Commit**

```powershell
git add R/source.R tests/test_downstream_reference.R
git commit -m "feat: validate downstream reference artifacts"
```

### Task 9: Generate Notebooks 03-09 and enforce notebook contracts

**Files:**
- Modify: `scripts/render_notebooks.py`
- Modify: `tests/test_notebook_contracts.py`
- Create: `notebooks/03_region3_anchor_reference.ipynb`
- Create: `notebooks/04_region1_anchor_validation.ipynb`
- Create: `notebooks/05_region2_anchor_validation.ipynb`
- Create: `notebooks/06_region1_3_consensus_reference.ipynb`
- Create: `notebooks/07_region4_consensus_mapping.ipynb`
- Create: `notebooks/08_region3_eosinophil_robustness.ipynb`
- Create: `notebooks/09_downstream_release_summary.ipynb`

**Interfaces:**
- Consumes: functions and stage schemas from Tasks 1-8.
- Produces: seven committed IRkernel notebooks with one tagged parameter cell and explicit stage dependencies.

- [ ] **Step 1: Write failing notebook contract tests**

```python
expected = [
    "03_region3_anchor_reference.ipynb",
    "04_region1_anchor_validation.ipynb",
    "05_region2_anchor_validation.ipynb",
    "06_region1_3_consensus_reference.ipynb",
    "07_region4_consensus_mapping.ipynb",
    "08_region3_eosinophil_robustness.ipynb",
    "09_downstream_release_summary.ipynb",
]
for name in expected:
    self.assertTrue((NOTEBOOKS / name).is_file())
    self.assertEqual(read_notebook(NOTEBOOKS / name)["metadata"]["kernelspec"]["name"], "ir")
```

Require every notebook to state Inputs, Outputs, Gate, Stop/fallback behavior, and interpretation limits. Require Notebooks 03-08 to name the Seurat/Harmony method and frozen gene set they use. Assert that Notebook 07 contains `FindTransferAnchors`, `TransferData`, and `MapQuery` through the reusable mapping wrapper but no `RunPCA`, `RunHarmony`, `FindNeighbors`, `FindClusters`, `FindAllMarkers`, or reference-label-definition call. Assert that Notebook 08 separates `REGION3`, `FROZEN_ELIGIBLE_REFERENCE`, and `SENSITIVITY_ONLY` panels.

- [ ] **Step 2: Run `python tests/test_notebook_contracts.py` and verify RED**

- [ ] **Step 3: Add downstream builders and validators**

Add a `DOWNSTREAM_NOTEBOOKS` mapping and one cell-builder per stage. Each parameter cell contains:

```r
PROJECT_ROOT <- "/dssg/home/acct-svetoslav_chakarov/svetoslav_chakarov/Lab_members/Yanan_Hu/YNH_Xenium"
PIPELINE_REPO <- file.path(PROJECT_ROOT, "adipose_analysis", "YNH_Xenium_scWAT")
QC_RUN_LABEL <- "full_notebook_qc_v2"
DOWNSTREAM_RUN_LABEL <- "region3_anchor_reference_v1"
DOWNSTREAM_CONFIG_PATH <- file.path(PIPELINE_REPO, "config", "downstream_reference_defaults.tsv")
MARKER_CONFIG_PATH <- file.path(PIPELINE_REPO, "config", "scwat_canonical_markers.tsv")
SEED <- 20260814L
EXECUTION_MODE <- "FULL_HPC"
```

Every notebook refuses to run when the required previous gate/artifact is absent. Notebook 06 chooses eligible regions from admission tables and implements fallback. Notebook 07 reads only the final frozen reference. Notebook 09 reads stage artifacts and does not refit models.

- [ ] **Step 4: Build notebooks with the renderer and verify structural contracts GREEN**

- [ ] **Step 5: Commit**

```powershell
git add scripts/render_notebooks.py tests/test_notebook_contracts.py notebooks/03_region3_anchor_reference.ipynb notebooks/04_region1_anchor_validation.ipynb notebooks/05_region2_anchor_validation.ipynb notebooks/06_region1_3_consensus_reference.ipynb notebooks/07_region4_consensus_mapping.ipynb notebooks/08_region3_eosinophil_robustness.ipynb notebooks/09_downstream_release_summary.ipynb
git commit -m "feat: add gate-first downstream notebooks"
```

### Task 10: Add local and HPC one-gate-at-a-time runners

**Files:**
- Create: `scripts/execute_downstream_subset.ps1`
- Create: `shell/run_downstream_reference_hpc.sh`
- Create: `slurm/scwat_downstream_reference.sbatch`
- Modify: `tests/test_notebook_contracts.py`

**Interfaces:**
- Consumes: `TARGET`, QC run label/root, downstream run label, config paths.
- Produces: executed notebook plus validated stage artifacts for exactly one target.

- [ ] **Step 1: Write failing runner-contract tests**

Require ordered targets:

```text
ANCHOR -> REGION1 -> REGION2 -> CONSENSUS -> REGION4 -> EOS -> RELEASE
```

Assert that each target resolves one notebook, validates required upstream gates, redirects `TMPDIR`, `TMP`, `TEMP`, `R_LIBS_USER`, and Python bytecode below the project, and rejects paths outside the project root. `ALL` may exist for automated HPC use but must call the same stage functions in order and stop at each failed prerequisite.

- [ ] **Step 2: Run and verify RED**

- [ ] **Step 3: Implement the bounded local contract runner**

Default to the downstream fixture and run the pure-R/state-machine tests, notebook parameter injection, JSON validation, path guards, and artifact-schema fixtures. If compatible Seurat/Harmony packages exist, permit execution of one bounded model target with `scripts/execute_r_notebook.R`. If they are absent, write `local_model_preflight.tsv` with `SKIP_LOCAL_MODEL_TEST_HPC_REQUIRED` and do not replace the model with another method. Never invoke full inputs locally.

- [ ] **Step 4: Implement the HPC runner and Slurm wrapper**

Use the fixed HPC root and require `Seurat >= 4.3.0`, `SeuratObject`, `harmony >= 1.2.0`, `leidenbase`, `Matrix`, `clue`, `ggplot2`, `dplyr`, `jsonlite`, and `IRkernel`. Print exact versions, filesystem capacity, input dimensions, and upstream artifact checks before execution. Use unlimited notebook timeout and post-run R validation. Default Slurm resources are 8 CPUs, 128 GB RAM, and 24 hours; record actual scheduler resources in provenance.

- [ ] **Step 5: Run all local contracts and record the model checkpoint honestly**

Verify the fixture, pure-R functions, state transitions, notebook structures, parameter injection, runner target resolution, and D:-only path guards. If Seurat/Harmony remain absent, verify `SKIP_LOCAL_MODEL_TEST_HPC_REQUIRED`; the first actual `ANCHOR` model run is then an HPC checkpoint. Stop and report any discrepancy rather than weakening thresholds.

- [ ] **Step 6: Commit**

```powershell
git add scripts/execute_downstream_subset.ps1 shell/run_downstream_reference_hpc.sh slurm/scwat_downstream_reference.sbatch tests/test_notebook_contracts.py
git commit -m "feat: add gated downstream execution runners"
```

### Task 11: Update the Phase 0-2 handoff, documentation, and progress record

**Files:**
- Modify: `scripts/render_notebooks.py`
- Modify: `notebooks/02_slide_QC_summary.ipynb`
- Modify: `tests/test_notebook_contracts.py`
- Modify: `README.md`
- Modify: `D:/Xiaonan/CODEX_projects/Yanan_Xenium/PROJECT_RESEARCH_PLAN_AND_PROGRESS.md`

**Interfaces:**
- Consumes: completed downstream pipeline and verified local checkpoint results.
- Produces: an explicit Phase 0-2-to-Notebook-03 handoff and copy/paste HPC runbook.

- [ ] **Step 1: Write failing handoff-documentation tests**

Require Notebook 02 and README to contain `Region 3 anchor`, `03_region3_anchor_reference.ipynb`, `REGION3_ANCHOR_FALLBACK`, `Region 4 mapping-only`, all seven targets, input roots, output roots, and the prohibition on interpreting sections as biological replicates.

- [ ] **Step 2: Run and verify RED**

- [ ] **Step 3: Update Notebook 02 and README**

Notebook 02 ends with the completed Phase 0-2 PASS evidence and exact Notebook 03 input paths; it does not perform downstream fitting. README provides one copy/paste block per HPC target, expected gate file, expected output directory, and the instruction to inspect each gate before submitting the next stage.

- [ ] **Step 4: Append verified progress**

Record changed files, local fixture source, cell/gene dimensions, each target's gate decision, test commands, limitations, and required HPC reruns. Do not claim full-data scientific results from the local subset.

- [ ] **Step 5: Run the complete bounded verification suite**

```powershell
& 'C:\Program Files\R\R-4.3.3\bin\Rscript.exe' tests/test_source.R
& 'C:\Program Files\R\R-4.3.3\bin\Rscript.exe' tests/test_extended_qc.R
& 'C:\Program Files\R\R-4.3.3\bin\Rscript.exe' tests/test_downstream_reference.R
& '<bundled-python-on-C:-read-only>' tests/test_notebook_contracts.py
```

Then validate all seven notebook JSON contracts and injected parameter cells. If Seurat/Harmony are absent, verify the explicit HPC-required skip and do not claim notebook model execution. If they are available, execute only the bounded fixture targets. In all cases, run `git diff --check` and verify that every generated local path begins with `D:/Xiaonan/CODEX_projects/Yanan_Xenium`.

- [ ] **Step 6: Commit**

```powershell
git add scripts/render_notebooks.py notebooks/02_slide_QC_summary.ipynb tests/test_notebook_contracts.py README.md
git commit -m "docs: document gate-first downstream workflow"
```

`PROJECT_RESEARCH_PLAN_AND_PROGRESS.md` is outside the nested Git repository. Update it on D: as required, verify its new checkpoint text directly, and do not attempt to stage it in this repository.

### Task 12: Provide and verify HPC execution checkpoints

**Files:**
- Modify if verification exposes issues: `shell/run_downstream_reference_hpc.sh`
- Modify if verification exposes issues: `slurm/scwat_downstream_reference.sbatch`
- Modify: `README.md`

**Interfaces:**
- Consumes: full Phase 0-2 HPC output with corrected PASS evidence.
- Produces: a reviewed, stepwise HPC execution handoff; no full-data result is claimed locally.

- [ ] **Step 1: Document preflight commands**

Include `bash -n` for both shell scripts; explicit `Seurat`, `SeuratObject`, and `harmony` version/API checks; input hashes/dimensions; project free space; and verification that `TMPDIR`, `TMP`, `TEMP`, `R_LIBS_USER`, and executed notebooks are below the HPC project root.

- [ ] **Step 2: Document one submission per gate**

Provide commands using the same `QC_RUN_LABEL` and `DOWNSTREAM_RUN_LABEL` for `ANCHOR`, `REGION1`, `REGION2`, `CONSENSUS`, `REGION4`, `EOS`, and `RELEASE`. After every job, show the exact TSV to inspect and the statuses that permit the next submission.

- [ ] **Step 3: Document download requirements**

Require the complete downstream output directory, executed notebooks, Slurm stdout/stderr, configuration files, and session/provenance manifests. Returned data go below a new versioned folder under `D:/Xiaonan/CODEX_projects/Yanan_Xenium/adipose_analysis/hpc_return/`.

- [ ] **Step 4: Run all locally available syntax/contract checks**

If Bash is unavailable locally, state `NOT_RUN_LOCAL_BASH_UNAVAILABLE` and make `bash -n` the first HPC checkpoint. Do not claim HPC execution.

- [ ] **Step 5: Commit any verified corrections**

```powershell
git add shell/run_downstream_reference_hpc.sh slurm/scwat_downstream_reference.sbatch README.md
git commit -m "docs: finalize downstream HPC checkpoints"
```

## Final acceptance checklist

- [ ] Region 3 anchor is fit only from Region 3 primary cells and 245 provisional genes.
- [ ] Anchor sensitivity uses strict cells, hotspot exclusion, and 67 conservative genes.
- [ ] Regions 1 and 2 are mapped and gated independently against the identical frozen anchor.
- [ ] Every admission combination selects the correct eligible reference sections.
- [ ] Expanded-consensus failure returns the byte-identical frozen Region 3 anchor.
- [ ] Region 4 cannot mutate or refit any reference component and has a separate release state.
- [ ] Eosinophil primary evidence is Region 3-centered; eligible-reference and Region 4 evidence remain subordinate and separated.
- [ ] Primary reference, Eosinophil, and Region 4 release states are orthogonal.
- [ ] All required artifacts reload with matching hashes, dimensions, cell order, feature order, thresholds, and provenance.
- [ ] All locally available pure-R, notebook, path, runner, and artifact tests pass without writing project artifacts to C:; any Seurat/Harmony model test not runnable locally is explicitly `SKIP_LOCAL_MODEL_TEST_HPC_REQUIRED`.
- [ ] All seven bounded Seurat/Harmony notebook targets pass on HPC before the full-data run proceeds.
- [ ] Full-data scientific decisions remain pending until the ordered HPC outputs are returned and audited.
