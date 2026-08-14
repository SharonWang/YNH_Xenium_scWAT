# scWAT Xenium Notebook QC Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build and locally validate one reusable R-kernel section notebook for Xenium Phases 0-2, run it independently over four scWAT sections, and summarize the resulting QC in a Cell-inspired slide-level notebook.

**Architecture:** Calculation, validation, I/O, aggregation, and plotting live in `R/source.R` with explicit inputs and returns. Thin R notebooks call those functions. Generated data, executed notebooks, figures, logs, libraries, and temporary files stay outside Git under the project root.

**Tech Stack:** R 4.3+, IRkernel/Jupyter, Matrix, jsonlite, ggplot2, Python nbformat/nbclient, Bash/Slurm.

## Global Constraints

- Scope is scWAT only; exactly four Xenium sections are expected.
- Canonical repository: `D:/Xiaonan/CODEX_projects/Yanan_Xenium/adipose_analysis/YNH_Xenium_scWAT`.
- Local project/data/temp root: `D:/Xiaonan/CODEX_projects/Yanan_Xenium`; no project or temporary writes to C:.
- HPC root: `/dssg/home/acct-svetoslav_chakarov/svetoslav_chakarov/Lab_members/Yanan_Hu`.
- Raw Xenium directories are read-only; local execution uses only the four 500-cell subsets.
- No automatic cell deletion; preserve raw sparse counts and every QC reason.
- Synthetic metadata blocks biology; poor-quality-cycle errors produce `HOLD` but diagnostic QC still writes when strict mode is false.
- Every code change updates repository tests/docs and the external progress markdown when project status changes.

## File map

- `.gitignore`: exclude generated artifacts.
- `R/source.R`: reusable Xenium functions.
- `config/eos_gene_sets.tsv`: verified 100-gene table.
- `tests/test_source.R`: unit tests.
- `notebooks/01_section_phase0_2_QC.ipynb`: one-section workflow.
- `notebooks/02_slide_QC_summary.ipynb`: four-section summary.
- `scripts/render_notebooks.py`: notebook generation/schema checks.
- `scripts/execute_local_subset.ps1`: D:-safe local execution.
- `shell/run_notebook_qc_hpc.sh`: HPC four-section loop.
- `slurm/scwat_notebook_qc.sbatch`: scheduler wrapper.
- `README.md`: input/output and run contract.

### Task 1: Repository contracts and panel configuration

**Files:** Create `.gitignore`, `config/eos_gene_sets.tsv`, `tests/test_source.R`.

**Interfaces:** Produces a 100-row `gene/gene_set` table and test harness that sources `R/source.R`.

- [ ] Write a failing test asserting 100 unique genes, group sizes `common=7`, `short_lived=47`, `long_lived=46`, no missing values, and existence of `R/source.R`.
- [ ] Run with D:-only temp storage; expect failure because `R/source.R` is absent:

```powershell
$env:TMPDIR='D:\Xiaonan\CODEX_projects\Yanan_Xenium\adipose_analysis\tmp'
$env:TEMP=$env:TMPDIR; $env:TMP=$env:TMPDIR
& 'C:\Program Files\R\R-4.3.3\bin\Rscript.exe' 'D:\Xiaonan\CODEX_projects\Yanan_Xenium\adipose_analysis\YNH_Xenium_scWAT\tests\test_source.R'
```

- [ ] Add the verified panel table and ignore `.ipynb_checkpoints/`, outputs, logs, temp, R libraries, executed notebooks, and R history/session files.
- [ ] Re-run panel assertions and confirm only the intended missing-source assertion remains red.
- [ ] Commit: `git commit -m "test: define Xenium QC repository contracts"`.

### Task 2: General path, metadata, discovery, and integrity functions

**Files:** Create `R/source.R`; modify `tests/test_source.R`.

**Interfaces:** Produce `assert_path_within()`, `validate_runtime_paths()`, `discover_xenium_sections()`, `discover_one_section()`, `create_synthetic_manifest()`, `validate_sample_manifest()`, `inventory_section_files()`, `read_mtx_dimensions()`, `validate_section_integrity()`, and `write_tsv()`.

- [ ] Add failing tests: descendants accepted; C:/siblings rejected; Regions 1-4 sorted; missing/duplicate sections rejected; synthetic mapping deterministic/balanced/non-interpretable; missing files reported; dimension mismatches rejected.
- [ ] Run tests and confirm the first missing-function failure.
- [ ] Implement explicit argument validation and informative errors. Require one exact region match. Every writer validates its destination before directory creation.
- [ ] Run tests; expect all configuration/path/metadata/discovery/integrity assertions to pass.
- [ ] Commit: `git commit -m "feat: add reusable Xenium configuration and integrity functions"`.

### Task 3: Alarm, panel, sparse import, and cell QC functions

**Files:** Modify `R/source.R`, `tests/test_source.R`.

**Interfaces:** Produce `extract_analysis_alarms()`, `read_custom_panel_genes()`, `reconcile_panel()`, `read_xenium_features()`, `read_xenium_barcodes()`, `import_xenium_mex()`, `robust_interval()`, and `calculate_xenium_cell_qc()`. QC returns `list(cell_metadata, thresholds, summary)` without subsetting.

- [ ] Add a tiny sparse-MEX fixture under the D: test temp root with aligned features/barcodes/cells, one alarm, and panel JSON. Test gene-only sparse dimensions, alarm level, missing/extra genes, deterministic thresholds, separate review flags, and unchanged cell count.
- [ ] Run and confirm failures at unimplemented functions.
- [ ] Implement with Matrix/jsonlite, median/MAD bounds with q01/q99 fallback, and control upper limit `max(0.05, q99.5)`.
- [ ] Re-run; expect `All reusable scWAT Xenium function tests passed.`
- [ ] Commit: `git commit -m "feat: add Xenium sparse import and section QC functions"`.

### Task 4: Artifacts and Cell-inspired plots

**Files:** Modify `R/source.R`, `tests/test_source.R`.

**Interfaces:** Produce `cell_style_theme()`, `section_palette()`, `plot_section_qc()`, `write_section_artifacts()`, `validate_section_artifacts()`, and `calculate_readiness_gates()`.

- [ ] Add failing tests for exact Region 1-4 palette names, equal spatial coordinates, complete artifact set, sparse-RDS reload, traversal rejection, HOLD for ERROR alarms, and PENDING for synthetic metadata.
- [ ] Run and confirm failure at `section_palette()`.
- [ ] Implement ggplot2-only helpers; save PDF and 300-dpi PNG; gzip cell TSV; preserve uncompressed raw-count RDS and verify after save.
- [ ] Re-run and confirm all outputs remain under D: temp and are cleaned afterward.
- [ ] Commit: `git commit -m "feat: add Xenium QC artifacts and Cell-inspired plots"`.

### Task 5: Parameterized section notebook

**Files:** Create `scripts/render_notebooks.py`, `notebooks/01_section_phase0_2_QC.ipynb`.

**Interfaces:** Consumes notebook parameters plus `R/source.R`; produces one complete section output bundle.

- [ ] Implement schema checks for R kernelspec, tagged parameter cell, parameters `PROJECT_ROOT/PIPELINE_REPO/INPUT_ROOT/REGION_ID/RUN_LABEL/METADATA_PATH/EXPECTED_SECTION_COUNT/SEED/STRICT_MODE`, ordered sections, bounded previews, and no C: path.
- [ ] Generate an incomplete notebook and observe custom validation fail on a missing required section.
- [ ] Add focused cells for setup, discovery, metadata, inventory, integrity, panel/alarms, sparse import, QC, figures, artifact writing, and reload checks.
- [ ] Run `python scripts/render_notebooks.py --validate`; expect nbformat and custom checks to pass.
- [ ] Commit: `git commit -m "feat: add parameterized section Phase 0-2 notebook"`.

### Task 6: Four-section aggregation and summary notebook

**Files:** Modify `R/source.R`, `tests/test_source.R`, `scripts/render_notebooks.py`; create `notebooks/02_slide_QC_summary.ipynb`.

**Interfaces:** Produce `validate_four_section_outputs()`, `read_slide_qc_outputs()`, `summarise_slide_qc()`, and `plot_slide_qc()`.

- [ ] Add four minimal output fixtures; test exact Region 1-4 coverage, rejection of missing/duplicates, combined row count, per-section proportions, worst-gate propagation, and four palette groups.
- [ ] Run and confirm failure at `validate_four_section_outputs()`.
- [ ] Implement aggregation plus cell yield, count/feature, review proportion, control, nucleus/segmentation, threshold, and readiness plots.
- [ ] Generate summary sections `Goal`, `Setup`, `Inputs`, `Completeness Checks`, `QC Results`, `Cell-style Figures`, `Readiness`, and `Outputs`; write combined TSV/RDS/figures/session info.
- [ ] Run unit/schema tests and commit: `git commit -m "feat: add four-section Xenium QC summary notebook"`.

### Task 7: Execute locally on all four subsets

**Files:** Create `scripts/execute_local_subset.ps1`; create/update `README.md`.

**Interfaces:** Produces executed notebooks and QC outputs under `adipose_analysis/scwat_qc_outputs/local_notebook_test`.

- [ ] Add an executor that creates D: temp before launching R/Jupyter, sets `TMPDIR/TMP/TEMP`, checks IRkernel/Jupyter, injects one Region at a time, and then executes the summary.
- [ ] Run reusable R tests.
- [ ] Execute Regions 1-4 separately; after each assert 500 cell rows, readiness file, reloadable RDS, and figures.
- [ ] Execute summary; assert 2,000 rows, four unique regions, reloadable summary, and preserved HOLD/PENDING status.
- [ ] Commit only code/docs: `git commit -m "test: validate notebook QC workflow on four scWAT subsets"`.

### Task 8: HPC loop and Slurm execution

**Files:** Create `shell/run_notebook_qc_hpc.sh`, `slurm/scwat_notebook_qc.sbatch`; modify `README.md`.

**Interfaces:** Consumes full `${PROJECT_ROOT}/adipose_data`; produces executed notebooks, four section bundles, slide summary, and logs outside Git.

- [ ] Add preflight for Rscript, Jupyter, IRkernel, Matrix/jsonlite/ggplot2, disk, paths, and one directory for each Region 1-4; set project-local temp/library/log/output paths.
- [ ] Loop through four parameter injections and stop if any section artifact validation fails; then execute the summary.
- [ ] Add Slurm wrapper using supplied `--chdir`, project-local logs, 8 CPUs, 96 GB, and 24 hours.
- [ ] On HPC run `bash -n shell/run_notebook_qc_hpc.sh` and `bash -n slurm/scwat_notebook_qc.sbatch`; expect zero exits.
- [ ] Document exact inputs/outputs/manual/Slurm/real-metadata commands and commit: `git commit -m "feat: add HPC execution for four-section notebook QC"`.

### Task 9: Final verification and progress record

**Files:** Modify `README.md`; modify external `D:/Xiaonan/CODEX_projects/Yanan_Xenium/PROJECT_RESEARCH_PLAN_AND_PROGRESS.md`.

**Interfaces:** Produces auditable validation evidence and a clean worktree.

- [ ] Run R tests, notebook schema validation, four-section execution, summary execution, reload checks, `git diff --check`, and `git status --short`.
- [ ] Reconcile with prior results: 2,000 cells, panel 100/100, ERROR/HOLD in Regions 1/2/4, no such error in Region 3, and zero deleted cells. Diagnose any discrepancy before completion.
- [ ] Record commits, commands, output root, observed per-section QC values, unresolved gates, and HPC status in README/progress markdown.
- [ ] Commit repository docs: `git commit -m "docs: record scWAT notebook QC validation"`.
- [ ] Handoff notebook/source paths, validation status, HPC copy/paste chunks, branch, and unpushed commits; do not mark full-data analysis complete before HPC review.
