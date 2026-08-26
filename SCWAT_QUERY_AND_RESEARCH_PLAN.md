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
