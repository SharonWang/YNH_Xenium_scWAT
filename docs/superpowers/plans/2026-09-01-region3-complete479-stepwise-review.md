# Region 3 complete-479 stepwise review implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a new top-to-bottom Region 3 notebook that retains every valid B1 analysis and plot, corrects invalid or overstated methods, uses all 479 genes and PCs 1–30, and produces auditable full-HPC outputs.

**Architecture:** Preserve the executed B1 notebook as immutable provenance. Add tested reusable functions to the existing active `R/source.R`, generate the new notebook reproducibly from a focused Python builder, and keep notebook-specific narrative and plots in the notebook. Every corrected analysis must display its checks and interpretation immediately after execution.

**Tech Stack:** R 4.6.1, Seurat/SeuratObject, Matrix, FNN, dplyr/tidyr/tibble, ggplot2, mclust, ComplexHeatmap, CellChat interaction database as a catalogue only, IRkernel/Jupyter, Python notebook JSON tooling, Git.

**Spec:** `docs/superpowers/specs/2026-09-01-region3-complete479-stepwise-review-design.md`

## Global Constraints

- All project files, temporary files, caches, tests and execution artifacts remain below `D:\Xiaonan\CODEX_projects\Yanan_Xenium` locally.
- Local R executable: `D:\Programs\R-4.6.1\bin\Rscript.exe`; local library: `D:\Programs\R_library`.
- HPC project root: `/dssg/home/acct-svetoslav_chakarov/svetoslav_chakarov/Lab_members/Yanan_Hu/YNH_Xenium`.
- Do not modify `notebooks/B1_Region3_primary_479.ipynb` or its executed outputs.
- Use exactly all 479 Gene Expression genes for normalization, scaling, PCA, graph construction, UMAP and clustering.
- Use exactly PCA dimensions 1–30 downstream and use names containing `PC1_30`.
- Analyse only cells with `primary_include_revised == TRUE`, where the mask is computed after cell-ID alignment as `qc_core_pass & !high_control_flag & !segmentation_multiplet_flag`.
- Inclusive Eos means original annotation Eos OR `REF_EOS_TIER1` OR `REF_EOS_TIER2`; preserve the route in `Eos_origin`.
- Region 3 is one section from one mouse; cell-level results are exploratory, not mouse-level biological inference.
- Full-data execution is performed on HPC. Local execution is bounded to a small subset and must not densify the full matrix.
- Full-data outputs use `scwat_downstream_outputs/full_notebook_qc_v2/04_region3_complete479_reviewed`; bounded local outputs use `adipose_analysis/scwat_downstream_outputs/local_region3_complete479_reviewed`.

---

### Task 1: Preserve and inventory the HPC-returned baseline

**Files:**
- Inspect: `R/source.R`
- Inspect: `notebooks/B1_Region3_primary_479.ipynb`
- Inspect: `notebooks/01_QC_Region1.ipynb` through `notebooks/02_slide_QC_summary.ipynb`
- Create: `docs/validation/2026-09-01-b1-region3-method-audit.md`
- Modify: `SCWAT_QUERY_AND_RESEARCH_PLAN.md`

**Interfaces:**
- Consumes: current dirty working tree and September 1 HPC-returned artifacts.
- Produces: an immutable baseline commit and a cell-by-cell audit mapping every unique B1 statistic/plot to `RETAIN`, `CORRECT`, or `RETAIN_WITH_CAVEAT`.

- [ ] **Step 1: Record exact baseline provenance**

Run:

```powershell
git -c safe.directory='D:/Xiaonan/CODEX_projects/Yanan_Xenium/adipose_analysis/YNH_Xenium_scWAT' status --short
git -c safe.directory='D:/Xiaonan/CODEX_projects/Yanan_Xenium/adipose_analysis/YNH_Xenium_scWAT' diff --stat
Get-FileHash -Algorithm SHA256 -LiteralPath 'D:\Xiaonan\CODEX_projects\Yanan_Xenium\adipose_analysis\YNH_Xenium_scWAT\notebooks\B1_Region3_primary_479.ipynb'
Get-FileHash -Algorithm SHA256 -LiteralPath 'D:\Xiaonan\CODEX_projects\Yanan_Xenium\adipose_analysis\YNH_Xenium_scWAT\R\source.R'
```

Expected: hashes are captured; no file is changed.

- [ ] **Step 2: Write the B1 method-audit table**

Create a markdown table with these mandatory rows:

```markdown
| B1 component | Decision | Reason | B2 treatment |
| QC counts and spatial masks | RETAIN | Direct descriptive evidence | Recompute after ID-safe attachment |
| LogNormalize/ScaleData/PCA | RETAIN_WITH_CAVEAT | Method is suitable; feature and reduction labels require correction | All 479 genes; PC1_30 names |
| Marker and Wang transfer | RETAIN_WITH_CAVEAT | Useful evidence; uncertainty was underrepresented | Retain scores and add Uncertain |
| Inclusive Eos union | RETAIN_WITH_CAVEAT | User-specified union; rescue route must remain visible | Add Eos_origin |
| Eos score distributions and heatmap | RETAIN_WITH_CAVEAT | Descriptive; extreme groups are score-derived | Retain with circularity warning |
| Cell-level Eos DE | CORRECT | One section and score-derived groups do not support biological DE | Exploratory association table |
| Nearest-cell analysis | CORRECT | Legacy pools could overlap | Disjoint IDs and positive-distance gate |
| Ligand-receptor products | CORRECT | Not CellChat inference | Spatial co-expression terminology and nulls |
```

- [ ] **Step 3: Validate the HPC-returned baseline before checkpointing**

Run notebook JSON parsing, R-cell parsing, `source.R` parsing and `git diff --check`. Record any failures as baseline evidence; do not repair them in this step.

- [ ] **Step 4: Commit the reviewed HPC baseline separately**

Stage only the reviewed current HPC-returned files and audit record. Confirm the staged list with `git diff --cached --name-status`, then commit with:

```powershell
git commit -m 'Checkpoint September HPC Region 3 analysis'
```

Expected: the later refactor diff can be reviewed independently from the imported HPC state.

---

### Task 2: Add failing contracts for the reviewed Region 3 workflow

**Files:**
- Create: `tests/test_region3_complete479.R`
- Modify: `tests/test_notebook_contracts.py`

**Interfaces:**
- Consumes: public functions sourced from `R/source.R` and the planned B2 notebook path.
- Produces: executable red gates for cell-ID alignment, revised mask semantics, 479 genes, PC1–30, Eos provenance and disjoint spatial pools.

- [ ] **Step 1: Write the R contract fixtures**

Add a four-cell metadata fixture with deliberately shuffled rows and assertions equivalent to:

```r
aligned <- align_seurat_metadata_by_cell_id(
  object_cells = c("c1", "c2", "c3", "c4"),
  metadata = data.frame(
    cell_id = c("c3", "c1", "c4", "c2"),
    qc_core_pass = c(TRUE, TRUE, FALSE, TRUE),
    high_control_flag = c(FALSE, FALSE, FALSE, TRUE),
    segmentation_multiplet_flag = c(FALSE, TRUE, FALSE, FALSE)
  )
)
stopifnot(identical(aligned$cell_id, c("c1", "c2", "c3", "c4")))
mask <- derive_primary_include_revised(aligned)
stopifnot(identical(mask, c(FALSE, FALSE, TRUE, FALSE)))
```

Also assert that duplicates, missing cells and missing flag values stop with informative messages.

- [ ] **Step 2: Write Eos-union and spatial-disjointness red tests**

Use explicit expected provenance:

```r
eos <- build_inclusive_eos_provenance(
  data.frame(
    original_subtype = c("Eosinophil", "B", "Macrophage", "ASC"),
    EosRef_call = c("REF_EOS_REST", "REF_EOS_TIER1", "REF_EOS_TIER2", "OUTSIDE_IMMUNE")
  ),
  annotation_col = "original_subtype",
  tier_call_col = "EosRef_call"
)
stopifnot(identical(eos$Eos_origin, c("ANNOTATION_EOS", "TIER1_RESCUE", "TIER2_RESCUE", "NOT_EOS")))
```

Create overlapping query/reference coordinates and assert `build_disjoint_spatial_pools()` stops; create disjoint coordinates and assert every nearest distance from `compute_spatial_neighbours()` is strictly positive.

- [ ] **Step 3: Write the B2 notebook red contract**

Require:

```python
path = NOTEBOOKS / "B2_Region3_complete479_stepwise_reviewed.ipynb"
text = notebook_text(read_notebook(path))
self.assertIn("primary_include_revised", text)
self.assertIn("FIX_GENESET <- rownames(region_spatial)", text)
self.assertIn("dims_use <- 1:30", text)
self.assertIn("Eos_origin", text)
self.assertIn("assert_disjoint_spatial_pools", text)
self.assertNotIn("install.packages(", text)
self.assertNotIn("umap_PC12", text)
```

- [ ] **Step 4: Run focused tests and verify red failures**

Run:

```powershell
& 'D:\Programs\R-4.6.1\bin\Rscript.exe' --vanilla tests/test_region3_complete479.R
& 'C:\Users\Xiaonan_Wang\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' -B -m unittest tests.test_notebook_contracts
```

Expected: missing functions and missing B2 notebook cause the intended failures.

- [ ] **Step 5: Commit the red contracts**

```powershell
git add tests/test_region3_complete479.R tests/test_notebook_contracts.py
git commit -m 'test: define Region 3 complete-panel contracts'
```

---

### Task 3: Implement alignment, PCA and clustering-stability helpers

**Files:**
- Modify: `R/source.R`
- Test: `tests/test_region3_complete479.R`

**Interfaces:**
- Produces:
  - `align_seurat_metadata_by_cell_id(object_cells, metadata, cell_id_col = "cell_id", require_complete = TRUE) -> data.frame`
  - `derive_primary_include_revised(metadata) -> logical`
  - `validate_complete_xenium_panel(object, assay = "Xenium", expected_genes = 479L) -> character`
  - `summarise_pca_qc(object, reduction = "pca", dims = 1:30, qc_vars) -> data.frame`
  - `adjusted_rand_index(labels_a, labels_b) -> numeric(1)`
  - `summarise_cluster_stability(metadata, cluster_cols) -> list(pairwise, sizes)`

- [ ] **Step 1: Implement exact ID alignment and revised-mask validation**

The alignment function must use `match(object_cells, metadata[[cell_id_col]])`, restore object-cell order, and assert `identical(aligned[[cell_id_col]], object_cells)`. The mask function must require logical/non-missing `qc_core_pass`, `high_control_flag` and `segmentation_multiplet_flag`, then return exactly:

```r
(metadata$qc_core_pass %in% TRUE) &
  !(metadata$high_control_flag %in% TRUE) &
  !(metadata$segmentation_multiplet_flag %in% TRUE)
```

- [ ] **Step 2: Implement panel and PCA diagnostics**

`validate_complete_xenium_panel()` must stop unless `nrow(object[[assay]]) == 479L`, gene names are non-empty and unique, and return the assay row names in their stored order. `summarise_pca_qc()` must return one row per PC/QC variable with Spearman rho and complete observation count.

- [ ] **Step 3: Implement dependency-free ARI and stability summaries**

Compute ARI from the contingency table so the pipeline does not require `mclust` merely for cluster comparison. Return `NA_real_` for undefined comparisons and include cluster sizes for every resolution.

- [ ] **Step 4: Run focused and existing source tests**

Expected: new contracts pass; `tests/test_source.R`, `tests/test_source_contract.R`, `tests/test_fixed_cell_qc.R` and `tests/test_scwat_initial_qc_contract.R` remain green.

- [ ] **Step 5: Commit**

```powershell
git add R/source.R tests/test_region3_complete479.R
git commit -m 'feat: add Region 3 cohort and stability helpers'
```

---

### Task 4: Implement reference-transfer, Eos-provenance and state helpers

**Files:**
- Modify: `R/source.R`
- Test: `tests/test_region3_complete479.R`

**Interfaces:**
- Produces:
  - `run_scwat_wang_transfer(reference, query, features, main_col, subtype_col, prefix, dims = 1:30, seed = 1234L) -> list`
  - `build_inclusive_eos_provenance(metadata, annotation_col, tier_call_col) -> data.frame`
  - `score_eos_state_programs(object, eos_gene_sets, assay = "Xenium") -> list(object, gene_summary, score_summary)`
  - `classify_eos_state_extremes(metadata, balance_col = "EosState_balance", quantile_probability = 0.10) -> factor`

- [ ] **Step 1: Extract Wang transfer from B1 into a documented function**

Normalize and PCA-transform the reference, intersect requested features with both assays, require at least 100 shared genes, run `FindTransferAnchors(..., reduction = "pcaproject", dims = 1:30)` and transfer main/subtype labels. Return predictions, anchors and the exact shared feature vector. Do not mutate the input query.

- [ ] **Step 2: Implement the exact inclusive Eos rule with provenance**

Assign precedence to original annotation Eos, then Tier 1, then Tier 2. Return `Eos_inclusive`, `Eos_origin`, `Eos_original_subtype` and the unchanged tier call. Assert that `Eos_inclusive` is identical to:

```r
metadata[[annotation_col]] == "Eosinophil" |
  metadata[[tier_call_col]] %in% c("REF_EOS_TIER1", "REF_EOS_TIER2")
```

- [ ] **Step 3: Extract Eos scoring without changing the B1 formulas**

Use all available 47 short-lived and 46 long-lived genes, derive the non-ribosomal short set using `!grepl("^Rp[ls]", gene)`, scale the state genes, compute short/long means, z-scores, balance/activity and detected-gene counts. Return gene availability and score summaries for immediate display.

- [ ] **Step 4: Add interpretation guards**

The extreme classifier must use q10/q90 plus positive component-score rules exactly as B1, but attach method text declaring the groups quantile-defined. Tests must confirm signature genes are tagged `DEFINITION_GENE` so later association tables cannot present them as independent validation.

- [ ] **Step 5: Run focused tests and commit**

```powershell
git add R/source.R tests/test_region3_complete479.R
git commit -m 'feat: add auditable reference and Eos helpers'
```

---

### Task 5: Implement corrected spatial and exploratory LR helpers

**Files:**
- Modify: `R/source.R`
- Test: `tests/test_region3_complete479.R`

**Interfaces:**
- Produces:
  - `build_disjoint_spatial_pools(coordinates, metadata, eos_col = "Eos_inclusive", cell_id_col = "cell", celltype_col) -> list(query, reference)`
  - `assert_disjoint_spatial_pools(query, reference, cell_id_col = "cell") -> invisible(TRUE)`
  - `compute_spatial_neighbours(query, reference, k = 15L, cell_id_col = "cell", celltype_col) -> list(closest, edges, by_type)`
  - `make_spatial_permutation_strata(coordinates, block_width, density, density_bins = 5L) -> factor`
  - `spatial_permutation_test(values, strata, statistic, n_permutations = 1000L, seed = 1234L) -> data.frame`
  - `score_spatial_lr_coexpression(edges, expression, interactions) -> data.frame`

- [ ] **Step 1: Implement disjoint-pool construction and hard stops**

Reference cells must be the exact complement of `Eos_inclusive`, not cells selected from an older annotation column. Stop on any intersecting ID, duplicated ID, non-finite coordinate or empty pool.

- [ ] **Step 2: Implement nearest-neighbour outputs with positive-distance validation**

Use `FNN::get.knnx()`. Return the single closest reference, 15-neighbour edges and nearest distance to every reference cell type. Stop if any query cell equals a matched reference cell or if any returned distance is `<= 0`.

- [ ] **Step 3: Implement deterministic spatially constrained nulls**

Define strata from fixed x/y blocks crossed with local-density quintiles. Permute labels or continuous state values only within strata containing at least two observations. Return observed statistic, null median, 2.5%/97.5% interval, empirical two-sided p-value `(1 + sum(abs(null) >= abs(observed))) / (1 + n_valid)` and valid permutation count.

- [ ] **Step 4: Implement LR co-expression scoring**

Filter interactions to simple ligand/receptor symbols present in the expression matrix. For each validated edge and direction, return ligand value, receptor value, both-detected indicator and product score. Name outputs `spatial_lr_coexpression`; do not use `communication_probability`, `CellChat_inference` or equivalent inferential labels.

- [ ] **Step 5: Test failures and deterministic results**

Assert overlap and zero-distance fixtures fail, identical seeds give identical null distributions, empirical p-values lie in `[0,1]`, and absent LR genes are reported rather than silently treated as zero.

- [ ] **Step 6: Run tests and commit**

```powershell
git add R/source.R tests/test_region3_complete479.R
git commit -m 'fix: validate Region 3 spatial analyses'
```

---

### Task 6: Build the new reader-facing B2 notebook

**Files:**
- Create: `scripts/build_region3_complete479_notebook.py`
- Create: `notebooks/B2_Region3_complete479_stepwise_reviewed.ipynb`
- Modify: `tests/test_notebook_contracts.py`

**Interfaces:**
- Consumes: functions from Tasks 3–5, existing config files, Region 3 QC bundle and Wang reference.
- Produces: one output-free R-kernel notebook with deterministic top-to-bottom execution and immediate methods/results/interpretation cells.

- [ ] **Step 1: Build the notebook skeleton**

Use Python JSON/notebook construction with alternating concise markdown and focused R code cells. Required top-level sections are exactly:

```text
0. Scope, decisions and reproducibility
1. Import Region 3 and align QC metadata
2. Review QC and create primary_include_revised cohort
3. Normalize all 479 genes and run PCA
4. Build PC1–30 graph, UMAP and clustering stability
5. Marker-based annotation
6. Wang-reference annotation
7. Consolidate annotations and uncertainty
8. Define inclusive Eosinophils with provenance
9. Score and review Eosinophil states
10. Exploratory Eosinophil expression associations
11. Corrected spatial-neighbour analysis
12. Exploratory spatial ligand–receptor co-expression
13. Save outputs, limitations and HPC handoff
```

- [ ] **Step 2: Retain every valid B1 output class**

Include corrected equivalents of all unique B1 outputs: input inventory, Seurat summaries, mask cross-tabs, spatial QC maps, flag burden, PCA elbow/loadings/heatmap/QC correlations, multi-resolution UMAP, marker coverage and score tables, marker feature plots, Wang-reference counts and transfer maps, consolidated annotation tables/UMAP/spatial maps, Eos tier performance and stacked proportions, Eos state histograms/scatters/detection diagnostics/dip test/mixture BIC/heatmap, exploratory association volcano and external-list concordance plots, continuous and extreme Eos spatial maps, nearest-type/distance/neighbourhood panels, suspicious-gene proximity checks and LR co-expression panels.

- [ ] **Step 3: Correct or caveat invalid B1 analyses in adjacent markdown**

Place markdown immediately before and after corrected cells. Mandatory statements include:

```markdown
The inclusive Eos definition is intentionally sensitive and includes Tier 1/2 rescue cells. `Eos_origin` must be used to distinguish annotation-defined from rescued cells.

The long-like and short-like groups are derived from the same signatures displayed below. Their association table is exploratory and is not independent validation or mouse-level differential expression.

Spatial query and reference pools are disjoint by cell ID. The notebook stops if any self-match or zero-distance match is detected.

Ligand–receptor values below are panel-limited spatial co-expression scores, not CellChat communication probabilities.
```

- [ ] **Step 4: Eliminate interactive-state defects**

Every object must be created before first use; every checkpoint filename must be defined once; no inline package installation is allowed; no empty cells remain; no cell reads `_spatial.rds` when the saved contract is `_spatial_passQC.rds`; and all plot options are set immediately before the relevant plot.

- [ ] **Step 5: Generate and validate the notebook**

Run the builder twice and confirm an unchanged hash on the second run. Validate JSON, R kernel metadata, output-free code cells, R syntax for every code cell, required headings and forbidden legacy strings.

- [ ] **Step 6: Commit**

```powershell
git add scripts/build_region3_complete479_notebook.py notebooks/B2_Region3_complete479_stepwise_reviewed.ipynb tests/test_notebook_contracts.py
git commit -m 'feat: add reviewed Region 3 complete-panel notebook'
```

---

### Task 7: Run bounded local verification and prepare HPC execution

**Files:**
- Create: `scripts/execute_region3_complete479_subset.ps1`
- Modify: `SCWAT_QUERY_AND_RESEARCH_PLAN.md`
- Modify: `README.md`

**Interfaces:**
- Consumes: B2 notebook, Region 3 subset input, active R source and local R installation.
- Produces: D:-local executed-subset notebook/test artifacts, complete test evidence and exact HPC commands.

- [ ] **Step 1: Add a D:-only bounded executor**

Set `R_LIBS_USER`, `TMPDIR`, `TMP`, `TEMP`, `PYTHONDONTWRITEBYTECODE` and Jupyter cache/config paths below the D: project. Inject subset paths and a reduced permutation count for smoke testing while leaving the notebook default at 1,000 full-HPC permutations.

- [ ] **Step 2: Execute a small Region 3 subset**

Run setup, import/alignment, revised mask, 479-gene validation, PCA1–30, clustering summaries, Eos provenance and a small spatial fixture. If local dependencies or memory prevent a stage, record the exact package/stage and retain the full HPC code without substituting a different method.

- [ ] **Step 3: Run the complete regression suite**

Run all R tests, Python notebook contracts, notebook JSON/R-cell parsing, deterministic builder check, `git diff --check` and a scan ensuring no output path references C:.

- [ ] **Step 4: Record verification and HPC commands**

Append actual pass/fail results to `SCWAT_QUERY_AND_RESEARCH_PLAN.md`. Add an HPC command block that sets the project root, activates the installed R kernel/environment, executes B2 top-to-bottom and writes the executed notebook beneath the declared Region 3 output root.

- [ ] **Step 5: Commit the verification handoff**

```powershell
git add scripts/execute_region3_complete479_subset.ps1 SCWAT_QUERY_AND_RESEARCH_PLAN.md README.md
git commit -m 'docs: add Region 3 HPC execution handoff'
```

---

### Task 8: Final review checkpoint

**Files:**
- Review: every file changed by Tasks 1–7

**Interfaces:**
- Consumes: implementation commits and verification logs.
- Produces: a release-ready local branch and an explicit list of checks deferred to full HPC execution.

- [ ] **Step 1: Reconcile the B1 audit against B2**

For every `RETAIN`, `CORRECT` or `RETAIN_WITH_CAVEAT` audit row, link the B2 section that implements it. Fail review if any unique valid B1 statistic or plot lacks a B2 counterpart.

- [ ] **Step 2: Review scientific claims**

Confirm the notebook never calls cell-level results biological replication, never treats signature-derived groups as independent validation, never calls LR products CellChat inference, and never hides Eos rescue provenance.

- [ ] **Step 3: Review Git and filesystem scope**

Confirm B1 is byte-identical to its baseline hash, no C: project/temp artifact was created, only intended D: repository files are staged, and HPC-returned data remain unmodified.

- [ ] **Step 4: Run final tests and commit any documentation-only corrections**

Expected: all locally runnable checks pass; heavy full-data notebook execution remains clearly marked `HPC_REQUIRED`.
