# scWAT QC Source Refactor and Fixed-Threshold Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Produce a documented, notebook-focused Xenium QC source, archive inactive functions, adopt fixed open-interval cell-complexity QC, and explicitly mark the invalid Region 3 zero-distance Eos neighbour analysis.

**Architecture:** `R/source.R` remains the only file sourced by active notebooks and contains every notebook-reachable function, transitive dependency, and active test/support entry point. Functions with no active notebook, test, or support-script consumer move to `R/source_bk.R`, which is archival and is not sourced by the QC notebooks. The four region notebooks generate identical region-specific QC contracts; the slide notebook only aggregates those contracts and never redefines masks.

**Tech Stack:** R 4.6.1, Matrix, ggplot2, Jupyter/IRkernel notebooks, Python `nbformat` for structure-preserving notebook edits, Git.

**Spec:** User-approved design in the Codex scWAT task on 2026-08-26.

## Global Constraints

- All implementation files, tests, temporary files, and generated artifacts must remain under `D:\Xiaonan\CODEX_projects\Yanan_Xenium`.
- Use `D:\Programs\R-4.6.1\bin\Rscript.exe` and `R_LIBS_USER=D:\Programs\R_library`.
- Preserve raw Xenium objects; masks are metadata fields and do not delete cells.
- Fixed primary complexity rule uses strict inequalities: `5 < nFeature_Xenium < 200` and `10 < nCount_Xenium < 1000`.
- `primary_include = qc_core_pass`; `strict_include = primary_include AND NOT qc_review_flag`.
- Region notebooks remain separate for Regions 1–4; Region 4 remains sensitivity/mapping-only in later downstream analysis.
- The slide summary consumes region-level masks and must not recompute or redefine them.
- The Eos nearest-cell result with overlapping query/reference identities is invalid and must be labelled at the exact notebook location; no biological conclusion may rely on it.

---

### Task 1: Contract tests for fixed QC masks

**Files:**
- Modify: `tests/test_source.R`

**Interfaces:**
- Consumes: `calculate_xenium_cell_qc(counts, cells, region_id, fixed_thresholds)` and `build_cell_downstream_masks(cell_metadata, spatial_hotspots, provenance)`.
- Produces: executable assertions for strict boundary behavior and mask nesting.

- [x] Add a literal sparse-matrix fixture containing cells at and immediately inside/outside all four fixed boundaries.
- [x] Assert exact boundary values fail and only open-interval values pass `qc_core_pass`.
- [x] Assert `primary_include` equals `qc_core_pass` and `strict_include` is nested within it.
- [x] Run the focused contract test and confirm the new assertions fail against the data-derived implementation for the intended reason.

### Task 2: Implement and document the fixed-threshold QC contract

**Files:**
- Modify: `R/source.R`
- Modify: `config/fixed_cell_qc_thresholds.tsv` or create it if absent

**Interfaces:**
- Consumes: sparse counts, aligned cell metadata, region identifier, and versioned fixed thresholds.
- Produces: cell metadata with `qc_core_pass`, review flags, threshold provenance, and downstream masks.

- [x] Add validated fixed-threshold loading/application helpers matching the Colon contract.
- [x] Change `calculate_xenium_cell_qc` to use fixed open bounds for `qc_core_pass` while retaining robust high-tail calculations only as diagnostic segmentation evidence.
- [x] Change `build_cell_downstream_masks` so primary and strict masks follow the approved definitions.
- [x] Add purpose/input/output documentation and explicit mask-rule provenance strings.
- [x] Run the new R contract assertions and confirm they pass.

### Task 3: Split active and archival R functions safely

**Files:**
- Modify: `R/source.R`
- Create: `R/source_bk.R`
- Modify: `tests/test_source.R`

**Interfaces:**
- Consumes: static notebook/test/support-script call inventory and transitive R-function dependencies.
- Produces: standalone active source and clearly labelled archival backup.

- [x] Add a source-contract test that sources `R/source.R` in a clean environment and exercises active public entry points.
- [x] Move only repo-inactive functions to `R/source_bk.R`; retain notebook-reachable functions, their dependencies, and active test/support functions.
- [x] Add section headers and function documentation describing purpose, required inputs, returned values, and written artifacts.
- [x] Verify active notebooks and tests do not source `source_bk.R`.
- [x] Source both files independently in controlled environments to check syntax and archive labelling.

### Task 4: Update five QC notebooks and mark the Eos error in place

**Files:**
- Modify: `notebooks/01_QC_Region1.ipynb`
- Modify: `notebooks/01_QC_Region2.ipynb`
- Modify: `notebooks/01_QC_Region3.ipynb`
- Modify: `notebooks/01_QC_Region4.ipynb`
- Modify: `notebooks/02_slide_QC_summary.ipynb`
- Modify: `notebooks/B1_Region3_primary_479.ipynb`
- Modify: `tests/test_notebook_contracts.py`

**Interfaces:**
- Consumes: documented fixed-threshold functions and existing region artifact bundles.
- Produces: separately runnable region notebooks, a non-redefining slide summary, and an explicit invalid-analysis warning beside the faulty neighbour code.

- [x] Add failing notebook-contract tests for the fixed threshold declaration, mask semantics, slide non-redefinition, and Eos overlap warning.
- [x] Edit notebooks structurally with `nbformat`, preserving cell order and existing outputs except where a changed explanatory cell requires replacement.
- [x] In `B1_Region3_primary_479.ipynb`, place a prominent warning immediately before the code that builds `reference_cells` from `Final_CellType_subtype` and calls `FNN::get.knnx()`.
- [x] Explain that rescued Eos were labelled through `Final_CellType_subtype_refined`, remained in the older-label non-Eos pool, and could self-match at distance zero; mark existing neighbour results invalid.
- [x] Run notebook structure and contract tests with temporary paths forced to D:.

### Task 5: Full verification, progress record, and Git integration

**Files:**
- Modify: `SCWAT_QUERY_AND_RESEARCH_PLAN.md`
- Modify: `README.md` if file descriptions or execution instructions changed

**Interfaces:**
- Consumes: all implementation and test changes.
- Produces: auditable test record and commit on `codex/notebook-qc-pipeline`.

- [x] Run R syntax/source checks, focused QC tests, extended tests, Python notebook tests, notebook JSON validation, and `git diff --check`.
- [x] Record all pass/fail results and any pre-existing environment discrepancies below.
- [x] Inspect the final diff for accidental notebook-output churn and any file outside D:.
- [x] Commit the isolated branch, merge it locally into `codex/notebook-qc-pipeline`, rerun verification, and push the updated branch to GitHub.

## Progress Log

### 2026-08-26 — Read-only inventory and approved design

- Confirmed clean starting branch `codex/notebook-qc-pipeline` at commit `86319f3`.
- Counted 144 functions in `R/source.R`; static analysis found 94 notebook-reachable functions before accounting for active test/support entry points.
- Confirmed Colon uses strict fixed bounds and defines `primary_include` only from the fixed core rule.
- Confirmed the zero-distance issue originates in `B1_Region3_primary_479.ipynb`, not inside an existing `source.R` nearest-neighbour function.
- User approved the recommended split and requested the warning at the exact faulty notebook location.

### 2026-08-26 — Isolated workspace and baseline

- Created D:-local worktree `.worktrees/source-qc-refactor` on branch `codex/source-qc-refactor`.
- Verified R 4.6.1 at `D:\Programs\R-4.6.1`.
- Baseline `test_source.R` failed an existing `CoordFixed` class assertion.
- Baseline `test_extended_qc.R` failed an existing Arrow `SKIP_ALLOWED` expectation.
- `python` was not on PATH; the Codex bundled Python executable was located for notebook validation, with all temp variables to be redirected to D:.

### 2026-08-26 — Fixed-threshold TDD and source split

- Added a focused boundary fixture and observed the expected red failure because the fixed-threshold contract/config did not exist.
- Added `config/fixed_cell_qc_thresholds.tsv`, validated threshold readers/application, and changed `qc_core_pass` to strict open intervals.
- Confirmed the focused green result: exact values 5/200 features and 10/1000 counts fail; only cells strictly inside both intervals pass.
- Changed masks to the approved Colon-consistent definitions: `primary_include = qc_core_pass`; `strict_include = primary_include & !qc_review_flag`.
- Archived 48 repo-inactive functions in `R/source_bk.R`; retained 98 active functions in `R/source.R`.
- Added function-level purpose/input/output comments throughout active `source.R` and a prominent Eos nearest-cell warning beside the active Eos refinement section.

### 2026-08-26 — Notebook updates and Eos error annotation

- Updated `01_QC_Region1.ipynb` through `01_QC_Region4.ipynb` to load the versioned fixed thresholds and pass them explicitly to `calculate_xenium_cell_qc()`.
- Updated region markdown to show strict inequalities, boundary behavior, and the new primary/strict mask definitions.
- Updated `02_slide_QC_summary.ipynb` to read but not redefine region masks and to assert that imported `primary_include` equals the fixed rule for every cell.
- Inserted a warning immediately before `reference_cells <- coords` in `B1_Region3_primary_479.ipynb`, plus code comments at the reference-pool and `FNN::get.knnx()` cells.
- Recorded the exact cause: refined-label rescued Eos remained in an unrefined-label non-Eos reference pool and could self-match at distance zero. Existing nearest-cell/neighbourhood results are not presently valid.
- Reworked the notebook updater to be transactional and idempotent after diagnosing newline-sensitive matching; running it twice produces one validation block and one Eos warning.

### 2026-08-26 — Local verification checkpoint

- `tests/test_fixed_cell_qc.R`: PASS.
- `tests/test_source_contract.R`: PASS.
- `tests/test_source.R`: PASS.
- `tests/test_extended_qc.R`: PASS.
- `tests/test_notebook_contracts.py`: 5/5 PASS using Python `unittest` (bundled Python lacks `pytest` and `nbformat`).
- Existing baseline discrepancies were reconciled accurately: ggplot2 4.x represents `coord_fixed()` as `CoordCartesian` with ratio 1, and installed Arrow correctly reports `PASS` rather than `SKIP_ALLOWED`.
- Notebook diff review found no stored-output churn: four region notebooks changed by approximately 45 lines each, the summary by 24 lines, and the large Region 3 Eos notebook by 14 inserted warning/comment lines.
- Full HPC execution remains required to regenerate section and slide artifacts under the new fixed primary cohort.

### 2026-08-26 — Executed local subset verification

- Updated `scripts/execute_local_subset.ps1` to use `D:\Programs\R-4.6.1\bin\Rscript.exe`, `R_LIBS_USER=D:\Programs\R_library`, and a run-specific D:-local temporary directory.
- Added fixed-threshold path injection for every region notebook and the summary notebook.
- Replaced stale expected-cell-count checks with direct rule equivalence checks on every saved mask row.
- The first Region 1 smoke test correctly failed at the obsolete C: R executable; after correction it exposed an existing LOCAL_SUBSET notebook bug where `transcript_gene_quality` was printed without being initialized.
- Initialized `transcript_gene_quality` as an empty data frame before the FULL_HPC/LOCAL_SUBSET branch in all four region notebooks; FULL_HPC still replaces it with Arrow-derived results.
- Region 1 rerun: all 117 R code cells executed and all section artifact checks passed.
- Full bounded run: Regions 1–4 and `02_slide_QC_summary.ipynb` executed successfully, including summary reload and mask-equivalence checks.
- Local 500-cell-per-region results: primary counts were Region 1 = 489, Region 2 = 492, Region 3 = 494, Region 4 = 492; strict counts were 460, 477, 470, and 482 respectively.
- Generated test outputs are under `D:\Xiaonan\CODEX_projects\Yanan_Xenium\adipose_analysis\scwat_qc_outputs\local_fixed_qc_all_sections` and are not part of the Git commit.

### 2026-08-26 — Final pre-commit gate

- Parsed both `R/source.R` and `R/source_bk.R` successfully with R 4.6.1.
- Re-ran all four R suites: fixed cell-QC, source contract, reusable function, and extended QC tests all passed.
- Re-ran six Python notebook-contract tests; all passed, including the point-of-use Eos warning assertion.
- Validated all four region notebooks and the slide summary with their appropriate notebook schemas.
- Parsed every R code cell in the four region notebooks, slide summary, and `B1_Region3_primary_479.ipynb`; no syntax errors were found.
- `git diff --check` reported no whitespace or conflict-marker errors. Line-ending notices reflect the repository's existing Windows checkout behavior.
- All test/cache/temp paths were redirected to D:. The bundled Python executable was read from C: but was configured not to write bytecode; no project, temporary, or test files were placed on C:.

### 2026-08-26 — Git integration

- Committed the isolated implementation as `958fa38` on `codex/source-qc-refactor`.
- Merged it into `codex/notebook-qc-pipeline` as merge commit `686feac`.
- Re-ran the complete verification suite from the merged target checkout; every gate passed.
- Pushed `codex/notebook-qc-pipeline` to the configured GitHub repository.

### 2026-08-27 — Colon-parity initial-QC notebook refactor

- User requested that the four scWAT `01_QC` notebooks and `02_slide_QC_summary.ipynb` use exactly the same concise reader flow as the Colon pipeline, with only scWAT paths, inputs, metadata, and four-region logic changed.
- Replaced the four previously executed 133-cell region notebooks with four output-free, region-locked 12-cell notebooks and replaced the 80-cell summary with an output-free 10-cell summary.
- Added a four-region scWAT bundle API with roxygen-style purpose, parameter, return, path-safety, sparse-alignment, and reload contracts.
- Kept initial technical readiness separate from later Region 3 anchor admission and Region 4 mapping-only decisions.
- A call-graph audit retained 51 actively reachable functions in `R/source.R` and moved 57 newly obsolete extended/evidence-only initial-QC functions to `R/source_bk.R`; previously archived downstream reference functions remain preserved there.
- Added roxygen-style documentation to every active named function and to the `%||%` helper.
- TDD red gates confirmed the previous 133/80-cell notebooks and missing scWAT bundle API; focused R and notebook contracts subsequently passed.
- Local D:-only Region 1 smoke execution passed. A complete bounded run of Regions 1-4 (500 cells each) plus the slide summary also passed and wrote test outputs under `adipose_analysis/scwat_qc_outputs/local_colon_parity_all`.
- Full-data execution is still required on HPC before using the initial-QC readiness results for downstream admission decisions.
- Committed the implementation as `565edcd` and merged it locally into `codex/notebook-qc-pipeline` as `b43ec07`; the complete post-merge suite passed.
- Two non-force GitHub push attempts failed because the current session could not connect to `github.com:443`. The local target branch remains clean and ahead of its remote; no force operation was attempted.

### 2026-09-01 — Region 3 complete-panel notebook review design

- Audited the executed `B1_Region3_primary_479.ipynb`, current `R/source.R`, September 1 HPC-returned QC tables and Region 3 downstream artifact inventory without changing them.
- Confirmed the user requires all 479 genes, `primary_include_revised`, PCs 1–30, clustering-stability assessment, marker plus Wang-reference annotation, and an inclusive Eos definition consisting of annotation Eos plus Tier 1/2 calls.
- Identified correctness risks to address in a new notebook: out-of-order execution, row-order mask recomputation, misleading PC12 names for PC1–30 analysis, forced Eos relabelling without provenance, circular cell-level Eos comparisons, overlapping spatial query/reference pools, and overinterpretation of ligand–receptor expression products.
- Wrote the proposed design to `docs/superpowers/specs/2026-09-01-region3-complete479-stepwise-review-design.md`. No implementation files or executed B1 outputs were changed at this checkpoint.

### 2026-09-01 — Approved design converted to implementation plan

- User approved the complete-479 design and required every valid B1 statistic, plot and analysis to be retained, with invalid methods corrected and commented at the point of use.
- Added a test-first, eight-task implementation plan at `docs/superpowers/plans/2026-09-01-region3-complete479-stepwise-review.md`.
- The plan preserves B1 as immutable provenance, adds a formal B1-to-B2 method audit, extracts reusable R functions, builds a deterministic B2 notebook, performs bounded D:-only subset tests, and leaves full-data execution to HPC.

### 2026-09-01 — September HPC baseline audit checkpoint

- Captured B1 SHA-256 `739BD6648CF48057E14699694556D33AD24E09110D0FDF31A6DD6109214D1037` and `R/source.R` SHA-256 `1A2AE63697CF87424A0570F6DA3DEF29928BD3B16FAE6822AFEC5A845E5F5E21` before new implementation.
- Added `docs/validation/2026-09-01-b1-region3-method-audit.md`, mapping every unique B1 analysis to `RETAIN`, `RETAIN_WITH_CAVEAT`, or `CORRECT` and specifying its B2 treatment.
- JSON parsing passed for all four HPC-returned QC notebooks, the slide summary and B1; none contains a stored Jupyter error object. `git diff --check` also passed.
- Local R 4.6.1 did not finish even a source-parse startup within 60 seconds when user/site startup files were disabled and all temporary paths were redirected to D:. The process was terminated and this is recorded as a local-environment validation gap; no alternative R installation was used.

### 2026-09-09 — Approved four-region tissue-branch implementation

- Reviewed the executed `B2_Region1_primary_479.ipynb` and classified the current full-panel workflow as useful exploratory analysis requiring reproducibility, confidence, spatial-null and tissue-boundary corrections before four-region reuse.
- Created the isolated D:-local worktree `.worktrees/split-region-domain-notebooks` on branch `codex/split-region-domain-notebooks`; the user checkout was not modified.
- The inherited baseline initially failed because five active helpers lacked an adjacent roxygen `@return` contract. Added accurate parameter/return documentation and restored all pre-existing R tests to PASS before feature implementation.
- Added tested reusable helpers for the revised primary mask, cell-ID-safe tissue partitions, branch selection, disjoint spatial pools, mandatory positive-distance checks, adjusted Rand index, cluster-stability summaries, PC-QC correlations, lymph-node candidate domains, and boundary-sensitivity Jaccard summaries.
- Added a deterministic builder and 12 output-free notebooks: all QC-passed, adipose-only, and lymph-node-only branches for each Region 1-4. Child notebooks consume a frozen domain manifest and never redefine the lymph-node boundary.
- Retained all 479 genes, `LogNormalize`, PC1-30 analysis, resolution plots, canonical-marker and Wang transfer, explicit `Uncertain` labels, Eosinophil annotation plus Tier1/2 evidence, continuous Eosinophil state, exploratory state-associated markers, corrected disjoint spatial-neighbour analysis, per-cell-type distances, and panel-limited ligand-receptor spatial co-expression.
- Corrected interpretation: 2.5-month Wang is primary and all-age is sensitivity; Eosinophil rescue does not overwrite principal cell type; mixture/tail results are descriptive; cell-level marker tests are not mouse-level DE; spatial summaries omit naive cell-level inferential p-values; Region 4 is visibly sensitivity-only.
- Replaced the convex-hull LN boundary with local confident lymphoid/DC enrichment plus bounded spatial expansion and five-setting sensitivity diagnostics. `LN_NOT_DETECTED` and minimum-cell gates prevent forced LN analysis.
- Structural tests for all 12 notebooks pass; every generated R code cell parses under R 4.6.1. A bounded Region 1 subset smoke test passed import, 500/500 mask matching, revised primary selection, 475 retained cells, all 479 genes, 30 PCs and three-seed Louvain clustering (median ARI 1). Maximum absolute PC-QC correlation was 0.900, confirming the intended technical-dominance review gate.
- Full Wang transfer, branch-specific annotation, complete spatial plots, optional dip/mixture diagnostics, exploratory LR scoring, and all full-data outputs remain HPC execution checkpoints.
- The legacy `tests/test_notebook_contracts.py` is already inconsistent with commit `d09892c`: it expects 12/10-cell initial-QC notebooks although the tracked files have 20/21 cells, and it expects `B1_Region3_primary_479.ipynb`, which that commit renamed to B2. This inherited suite remains a documented baseline discrepancy and was not rewritten as part of the tissue-branch implementation.
- Hardened two clean-kernel edge cases found during validation: optional diagnostic packages now yield typed skipped outputs, and path containment is validated before the first output/temp directory is created.
- Added reload-validated Seurat checkpoints with unique stage filenames for branch input, PCA, annotation, and final objects; each manifest records dimensions, feature/cell identity validation, byte size and MD5.
- Harmony is explicitly deferred: a single-section notebook has no defensible batch factor. It will be used only in the later admitted-region consensus, with section as technical batch and mouse retained as biological replicate.

### 2026-09-09 — Eosinophil composition-plot regression fix

- Reproduced the reported `could not find function "theme_cell"` failure in a fresh D:-local R 4.6.1 test using a real minimal Seurat object and the split-notebook package context.
- Root cause: `plot_eos_call_by_subtype()` retained one dependency on the legacy notebook-local `theme_cell()` helper after the reusable plotting theme had been standardized as `cell_style_theme()` in `R/source.R`.
- Replaced only that hidden theme dependency with `cell_style_theme(base_size = base_size)`; the plot's data preparation, ordering, labels, colours and return contract are unchanged.
- Added a regression test requiring a valid ggplot, the complete 3-subtype by 7-call plotting grid, and Eosinophil-first ordering without defining `theme_cell()`.
- The focused test now passes. The test also showed that this legacy helper uses the attached `%>%` operator; current split notebooks satisfy that declared setup dependency by attaching dplyr. Broader pipe refactoring was deliberately kept outside this targeted bug fix.

### 2026-09-09 — Eos spatial, CellChat and Wang-integration design

- User approved an architectural extension covering the mclust runtime error, consistent cell/macaron plots, downstream use of `Final_CellType_subtype`, biological cell-type/marker ordering, four explicit Eosinophil spatial-neighbour steps, spatial CellChat inference for top KNN-associated cell types, and a joint Wang–Xenium embedding.
- Added `docs/superpowers/specs/2026-09-09-eos-spatial-wang-integration-design.md` with explicit inputs, thresholds, statistical interpretation, outputs, skip/failure gates and local-versus-HPC verification boundaries.
- The design keeps CellChat group-based: continuous `EosState_balance` drives KNN association and prespecified state-enriched tails, while CellChat evaluates those groups. Section-level permutation results are not labelled mouse-level inference.
- Existing Wang label transfer remains primary annotation evidence; the new joint embedding is a supporting cross-modality concordance diagnostic and cannot overwrite Xenium labels or allow Region 4 to define the primary reference.

### 2026-09-09 — Test-first implementation plan

- Converted the approved design into seven independently testable tasks in `docs/superpowers/plans/2026-09-09-eos-spatial-wang-integration.md`.
- The plan sequences the namespace-safe mclust wrapper, cell/macaron ordering and plot contracts, four-stage Eosinophil KNN analysis, gated spatial CellChat adapter, Wang–Xenium integration diagnostics, deterministic regeneration of all 12 B2 notebooks, and bounded local/HPC verification.
- Each behavior change begins with a failing test and ends with a focused commit. Full-data CellChat and Wang integration remain HPC checkpoints rather than unverified local claims.

### 2026-09-09 — Extended Eosinophil/Wang notebook implementation checkpoint

- Implemented and tested a namespace-safe `mclust` diagnostic wrapper. It supplies the `mclustBIC` binding expected by `Mclust()` and returns typed `PASS`, skipped or runtime-failure evidence instead of terminating a notebook.
- Added shared scWAT biological label ordering, stable macaron palettes, cell-style ggplot formatting and marker-block ordering. `Final_CellType_subtype_with_uncertain` is retained for annotation review, while `Final_CellType_subtype` is explicitly used for downstream grouping and plots.
- Rebuilt spatial analysis as four visible stages: disjoint Eos/reference pools; k=1 and k=15 neighbour composition; distance to every adequately represented reference cell type; and continuous-state KNN association with deterministic top-three short-like and long-like neighbour types.
- Added a spatial CellChat adapter using normalized Xenium expression, aligned centroids and a physical scale derived from cell area. The lower and upper 30% of continuous Eos state define adequately sized CellChat groups; raw p<0.05 and BH FDR<0.10 are reported as exploratory within-section results only.
- Added an optional subtype-balanced 2.5-month Wang–Xenium Seurat CCA integration using shared panel genes and PCs 1–30. Label transfer remains primary; the joint embedding, cluster composition and cross-dataset Eos-neighbour diagnostics are supporting concordance evidence only.
- Regenerated all 12 Region 1–4 all-QC/adipose/LN notebooks deterministically. Structural contracts passed for every notebook and every R code cell parsed successfully.
- The D:-local Region 1 smoke test passed with 475 revised-primary cells, all 479 genes, 30 PCs and median three-seed clustering ARI=1. The observed maximum absolute PC–QC correlation was 0.900 and therefore remains a required visual review item, not an automatically regressed covariate.
- Local `mclust` and CellChat packages are unavailable, so package-gated skip behavior was verified locally. Full CellChat and Wang–Xenium integration, plots and output tables require HPC execution before scientific interpretation.
- Added a regression-tested non-fatal empty-pool gate so a small adipose/LN branch with no Eosinophils can finish and report `SKIPPED_EMPTY_EOS_OR_REFERENCE_POOL`; ID overlap, duplicate cross-pool coordinates and zero/negative neighbour distances remain hard failures.
- Final feature verification passed deterministic two-run notebook generation, 4 Python notebook contracts, R parsing of every generated code cell, seven focused/source/QC R suites, the D:-local subset smoke test and `git diff --check`.
- Completion-check review then added explicit optional-package version tables, saved PNG/PDF artifacts for the Wang, Steps 10.1–10.4, mclust and CellChat diagnostics, persisted mclust status/BIC evidence, a subtype and cross-dataset-distance Wang integration panel, optional Wang/CellChat checkpoints, and matching README HPC/download instructions. The subset smoke test now also exercises biological ordering, disjoint pools, exact k=1/k=15 edge counts, positive distances and typed optional-module skips.

### 2026-09-09 — Directional plots, CellChat API compatibility and safe exports

- Reworked all focal UMAP/spatial displays in the 12 generated B2 notebooks so contextual cells are drawn first with smaller pale points and the target cells are drawn last with larger, darker points. Figure subtitles were removed.
- Restored the canonical-marker DotPlot gradient to `#D9D9D9`–`#5A2F5E` and retained biological cell-type and marker-block ordering.
- Rebuilt Steps 10.1–10.4 around the requested Eosinophil-state tails. Intermediate cells are grey spatial context only; tail-comparison panels exclude them. Signed composition, relative proximity and continuous-state association consistently place long-like evidence to the right and short-like evidence to the left, ordered from strongest long-like to strongest short-like tendency.
- Increased the prespecified KNN recipient selection from three to five cell types per Eosinophil-state direction. The reported significant CellChat subset now contains only outgoing short-like/long-like Eosinophil interactions with the matching state-specific top-five neighbour set; the complete communication table remains available for audit.
- Diagnosed the HPC CellChat 2.2.0.9001 construction error as an API mismatch. Added a tested adapter that supplies current `spatial.factors` (`ratio`, `tol`) or legacy `scale.factors` only after inspecting `createCellChat()` formals. Xenium centroids are treated as micrometre coordinates (`ratio = 1`), with tolerance derived from median equivalent cell radius.
- Reproduced the final-save `length of 'dimnames' [2] not equal to array extent` error with an `AsIs` matrix/list column. Added TSV-safe rectangularization to both compressed and uncompressed writers; the exact regression fixture now round-trips successfully.
- Regenerated all 12 notebooks and passed Python structural contracts, R parsing of every code cell, seven focused/source/QC R suites and the D:-local Region 1 subset smoke test (475 cells, 479 genes, 30 PCs, median clustering ARI 1). CellChat itself is unavailable locally, so both constructor signatures were validated with injected API fixtures; full CellChat execution remains an HPC checkpoint.

### 2026-09-10 — All-age Wang reference, Eos-state propagation and largest-cluster LN rule

- Changed the complete all-age Wang object to the primary reference for label reconciliation, Eosinophil calibration and subtype-balanced Wang–Xenium CCA integration. The 2.5-month transfer remains an explicit age-matched sensitivity comparison.
- Rebuilt the Eosinophil-state histogram and signature scatter with separate layers: intermediate cells are drawn first in grey, and short-like/long-like cells are drawn last in dark blue/red. Both plot objects are explicitly printed.
- Traced the three `EosState_extreme not found` failures to metadata being dropped at table-construction boundaries. `calculate_eos_distance_by_cell_type()` now propagates the tail label into its cell-level result, `rank_eos_state_knn_associations()` propagates it into `per_eos`, and Step 10.1 uses the enriched `spatial_pools$all` table rather than the earlier coordinate table.
- Replaced lymph-node radius expansion with the requested largest-cluster rule. After local lymphoid-enrichment selection, DBSCAN runs on preliminary core coordinates with `eps=80` and `minPts=10`; cluster 0 is noise, only the largest non-zero cluster is eligible, deterministic ties use the smallest cluster ID, and the retained cluster must still contain at least 100 cells.
- Added `lymph_node_dbscan_cluster_sizes.tsv` and cell-level preliminary-core/DBSCAN-cluster fields so the selection can be reviewed directly.

### 2026-09-10 — Nested `x` final-save regression fix

- Reproduced the HPC error `Nested TSV column has incompatible row count: x` with the same structural condition: a two-row outer result containing a three-row nested `AsIs` data-frame column named `x`.
- Root cause: the shared TSV rectangularizer assumed every nested matrix/data-frame column had one nested row per outer row. Optional-package result objects can legally violate that assumption, so the final output checkpoint stopped before writing the remaining audit files.
- Preserved the existing behavior for aligned nested columns. A non-aligned nested object is now retained once as a deterministic serialized payload, with explicit `__nested_rows__` and `__nested_cols__` audit fields; no row-wise relationship is invented and no result column is silently discarded.
- Moved rectangularization inside the path-aware error handlers for both plain and gzipped TSV writers. Any future export failure will identify the exact destination file.
- Added an exact regression fixture to `tests/test_eos_extended_helpers.R`; the pre-fix test reproduced the reported error and the post-fix round trip retained two outer rows plus the complete nested payload and its 3-by-2 dimensions.
