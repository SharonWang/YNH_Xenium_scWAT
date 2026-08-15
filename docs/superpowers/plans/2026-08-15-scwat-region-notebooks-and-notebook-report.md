# scWAT Region Notebooks and Notebook-Only QC Report Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Provide four independently runnable scWAT section QC notebooks and make `02_slide_QC_summary.ipynb` the sole slide-level QC report.

**Architecture:** Keep `01_section_phase0_2_QC.ipynb` as the canonical generated template and deterministically render four parameter-locked copies. Keep calculations in `R/source.R`; use the notebooks only for orchestration and bounded display. Execute named section notebooks individually or sequentially on HPC, then execute the single summary notebook against their validated artifact bundles.

**Tech Stack:** R 4.3/IRkernel, Jupyter notebook JSON, Python 3 standard library, PowerShell local runner, Bash/Slurm HPC runner, existing R QC functions and artifact validators.

## Global Constraints

- Scope is scWAT only.
- All project inputs, outputs, logs, caches, and temporary files must stay below `D:/Xiaonan/CODEX_projects/Yanan_Xenium` locally.
- No implementation, build, test, or temporary artifact may be written to C:.
- HPC root is `/dssg/home/acct-svetoslav_chakarov/svetoslav_chakarov/Lab_members/Yanan_Hu/YNH_Xenium`.
- Full transcript data must be aggregated lazily and never collected in full into R memory.
- Full-data execution is HPC-only; local validation uses the deterministic subset or structural checks.
- Exact poor-cycle identity requires 10x diagnostics; affected genes remain `CANDIDATE_NOT_CONFIRMED` until external evidence is reviewed.
- Cells are not automatically deleted by QC.
- Sections 62308/62309 belong to Mouse 1 and 62310/62311 to Mouse 2; section is technical and mouse is biological.

---

### Task 1: Region-notebook generation contract

**Files:**
- Create: `tests/test_notebook_contracts.py`
- Modify: `scripts/render_notebooks.py`
- Create: `notebooks/01_section_phase0_2_QC_Region1.ipynb`
- Create: `notebooks/01_section_phase0_2_QC_Region2.ipynb`
- Create: `notebooks/01_section_phase0_2_QC_Region3.ipynb`
- Create: `notebooks/01_section_phase0_2_QC_Region4.ipynb`

**Interfaces:**
- Consumes: `section_cells()` and `inject_parameters(source_path, output_path, assignments)` from `scripts/render_notebooks.py`.
- Produces: `build_region_notebooks(repo: Path) -> list[Path]` and CLI flag `--build-regions`; four valid R-kernel notebooks whose tagged parameter cell contains exactly one matching `REGION_ID` assignment.

- [ ] **Step 1: Write the failing notebook-contract test**

```python
def test_four_region_notebooks_are_parameter_locked():
    expected = {f"Region_{i}" for i in range(1, 5)}
    observed = set()
    for index in range(1, 5):
        path = NOTEBOOKS / f"01_section_phase0_2_QC_Region{index}.ipynb"
        assert path.exists()
        notebook = json.loads(path.read_text(encoding="utf-8"))
        parameter_cell = next(cell for cell in notebook["cells"] if "parameters" in cell.get("metadata", {}).get("tags", []))
        match = re.search(r'^REGION_ID\s*<-\s*"(Region_[1-4])"$', "".join(parameter_cell["source"]), re.MULTILINE)
        assert match
        observed.add(match.group(1))
        assert notebook["metadata"]["kernelspec"]["name"] == "ir"
    assert observed == expected
```

- [ ] **Step 2: Run the contract test and confirm RED**

Run from the D: repository with D:-local temp variables:

```powershell
$env:TEMP='D:\Xiaonan\CODEX_projects\Yanan_Xenium\adipose_analysis\tmp'
$env:TMP=$env:TEMP
$env:PYTHONDONTWRITEBYTECODE='1'
& 'C:\Users\Xiaonan_Wang\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' -m unittest tests.test_notebook_contracts -v
```

Expected: FAIL because the four named notebooks do not exist.

- [ ] **Step 3: Implement deterministic region generation**

Add a generator that first writes the canonical section notebook, then injects only the matching region assignment into each committed copy:

```python
REGION_NOTEBOOKS = {
    "Region_1": "01_section_phase0_2_QC_Region1.ipynb",
    "Region_2": "01_section_phase0_2_QC_Region2.ipynb",
    "Region_3": "01_section_phase0_2_QC_Region3.ipynb",
    "Region_4": "01_section_phase0_2_QC_Region4.ipynb",
}

def build_region_notebooks(repo):
    source = repo / "notebooks" / "01_section_phase0_2_QC.ipynb"
    built = []
    for region_id, filename in REGION_NOTEBOOKS.items():
        output = repo / "notebooks" / filename
        inject_parameters(source, output, {"REGION_ID": region_id})
        validate_notebook(output, "section", expected_region_id=region_id)
        built.append(output)
    return built
```

Extend `validate_notebook()` with optional `expected_region_id` and add `--build-regions`.

- [ ] **Step 4: Generate notebooks and confirm GREEN**

Run `render_notebooks.py --build-section --build-regions`, then rerun the unit test. Expected: four named notebooks pass JSON, R-kernel, fixed-region, and canonical-structure checks.

- [ ] **Step 5: Commit the region-notebook contract**

```bash
git add tests/test_notebook_contracts.py scripts/render_notebooks.py notebooks/01_section_phase0_2_QC_Region*.ipynb
git commit -m "feat: add independently runnable region QC notebooks"
```

---

### Task 2: Notebook-only slide QC report

**Files:**
- Modify: `tests/test_notebook_contracts.py`
- Modify: `scripts/render_notebooks.py`
- Modify: `notebooks/02_slide_QC_summary.ipynb`

**Interfaces:**
- Consumes: `read_slide_qc_outputs()`, `summarise_slide_qc()`, `read_extended_slide_qc_outputs()`, and `summarise_extended_slide_qc()` from `R/source.R`.
- Produces: a reader-facing summary notebook with sections `## TL;DR and QC decision`, `## Inputs and validation`, `## Core QC distributions`, `## Alarm evidence and candidate genes`, `## Subset versus full-data burden`, `## Spatial QC`, `## Within-mouse concordance`, and `## Final QC decision and next actions`.

- [ ] **Step 1: Extend the contract test and confirm RED**

Add assertions that the summary notebook contains every required report heading, displays `CANDIDATE_NOT_CONFIRMED` and `REQUIRES_10X_DIAGNOSTICS`, and contains no reference to `build_report.R` or `report.html`.

Run the test. Expected: FAIL because the current summary notebook uses the earlier question/readiness headings and lacks the complete decision block.

- [ ] **Step 2: Restructure `summary_cells()`**

Use short markdown introductions and focused R cells. Build a bounded decision table from generated outputs:

```r
qc_decision <- merge(
  slide_summary$section_summary,
  slide_summary$readiness[, c("region_id", "status")],
  by = "region_id", all.x = TRUE, sort = FALSE
)
qc_decision$required_next_action <- ifelse(
  qc_decision$region_id %in% c("Region_1", "Region_2"),
  "Obtain 10x poor-cycle diagnostics",
  ifelse(qc_decision$region_id == "Region_3",
         "Review morphology in FDR-positive hotspot bins",
         "Obtain 10x diagnostics and review segmentation/cell area")
)
qc_decision
```

Display only bounded tables, existing Cell-style plots, and explicit interpretation limits. Keep all file writing in the final output cell through existing artifact writers.

- [ ] **Step 3: Regenerate the summary notebook and confirm GREEN**

Run `render_notebooks.py --build-summary`, validate it as type `summary`, then rerun `tests.test_notebook_contracts`. Expected: required headings and evidence labels pass.

- [ ] **Step 4: Run R unit tests**

Run `tests/test_source.R` and `tests/test_extended_qc.R` with R temporary paths on D:. Expected: PASS with no change to core/extended artifact calculations.

- [ ] **Step 5: Commit the notebook-only report**

```bash
git add tests/test_notebook_contracts.py scripts/render_notebooks.py notebooks/02_slide_QC_summary.ipynb
git commit -m "feat: make slide QC summary the review report"
```

---

### Task 3: Independent local and HPC execution paths

**Files:**
- Modify: `tests/test_notebook_contracts.py`
- Modify: `scripts/execute_local_subset.ps1`
- Modify: `shell/run_notebook_qc_hpc.sh`
- Modify: `slurm/scwat_notebook_qc.sbatch`
- Modify: `README.md`

**Interfaces:**
- Consumes: the four committed region notebooks and `02_slide_QC_summary.ipynb`.
- Produces: optional environment variable `REGION_ID=Region_1|Region_2|Region_3|Region_4|ALL`; `ALL` runs four sections then summary, while a single region runs and validates only its own notebook and bundle.

- [ ] **Step 1: Add runner-contract tests and confirm RED**

Assert that both runners reference the four committed filenames, the HPC runner accepts a single valid region or `ALL`, a single-region run does not run the slide summary, and the README provides one-region plus summary-only commands. Run the test and expect failure against the current template-only runners.

- [ ] **Step 2: Update the local subset runner**

Add a validated `RegionId` parameter defaulting to `ALL`. Resolve each region to its named notebook, inject runtime paths without changing the fixed identity, and run the summary only when all four bundles are requested.

- [ ] **Step 3: Update the HPC runner and Slurm wrapper**

Use:

```bash
REGION_ID="${REGION_ID:-ALL}"
case "${REGION_ID}" in ALL) regions=(Region_1 Region_2 Region_3 Region_4);; Region_[1-4]) regions=("${REGION_ID}");; *) exit 2;; esac
```

Map each region to `notebooks/01_section_phase0_2_QC_RegionN.ipynb`. Execute and validate one region when requested; execute the summary only for `ALL`. Keep all environment and preflight paths below `PROJECT_ROOT`.

- [ ] **Step 4: Update README copy/paste chunks**

Document inputs and outputs for each region notebook, show one-region HPC execution, four-region sequential execution, and summary-only rerun. State that local checks are subset/structural and that full-data evidence must be downloaded after HPC execution.

- [ ] **Step 5: Confirm GREEN and shell syntax**

Run the Python contract suite. If local Bash is unavailable, record `bash -n shell/run_notebook_qc_hpc.sh` and `bash -n slurm/scwat_notebook_qc.sbatch` as exact HPC validation commands rather than claiming they passed.

- [ ] **Step 6: Commit runner and documentation changes**

```bash
git add tests/test_notebook_contracts.py scripts/execute_local_subset.ps1 shell/run_notebook_qc_hpc.sh slurm/scwat_notebook_qc.sbatch README.md
git commit -m "feat: run scWAT QC one region at a time"
```

---

### Task 4: Remove duplicate report package and update project record

**Files:**
- Modify: `tests/test_notebook_contracts.py`
- Delete: `reports/2026-08-15_full_hpc_qc_acceptance/build_report.R`
- Delete: `reports/2026-08-15_full_hpc_qc_acceptance/report.html`
- Delete: `reports/2026-08-15_full_hpc_qc_acceptance/artifact.json`
- Delete: `reports/2026-08-15_full_hpc_qc_acceptance/candidate_cycle_affected_genes_affected_only.tsv`
- Modify: `D:/Xiaonan/CODEX_projects/Yanan_Xenium/PROJECT_RESEARCH_PLAN_AND_PROGRESS.md`

**Interfaces:**
- Consumes: the approved design and completed notebook pipeline.
- Produces: one report implementation (`02_slide_QC_summary.ipynb`) and a current progress entry with local/HPC validation boundaries.

- [ ] **Step 1: Add a failing retirement test**

Assert that the 2026-08-15 report-package directory is absent or contains no `build_report.R`, `report.html`, `artifact.json`, or duplicated candidate TSV. Run and expect failure because those tracked files currently exist.

- [ ] **Step 2: Delete only the duplicate generated report package**

Remove the four tracked report-package files with an explicit patch. Do not delete returned HPC outputs or the earlier historical report.

- [ ] **Step 3: Update the progress Markdown**

Record the four notebook names, the notebook-only summary report, tests actually run, full-HPC tests deferred, and current gates: Region 1 HOLD, Region 2 HOLD, Region 3 conditional morphology review, and Region 4 HOLD plus segmentation/cell-area review.

- [ ] **Step 4: Confirm GREEN and commit**

Run the contract suite and `git diff --check`. Commit the retirement and progress record with:

```bash
git add tests/test_notebook_contracts.py reports/2026-08-15_full_hpc_qc_acceptance D:/Xiaonan/CODEX_projects/Yanan_Xenium/PROJECT_RESEARCH_PLAN_AND_PROGRESS.md
git commit -m "docs: consolidate QC reporting in the summary notebook"
```

---

### Task 5: End-to-end bounded validation and handoff

**Files:**
- Modify only if a failing test reveals a defect in the files above.

**Interfaces:**
- Consumes: four named region notebooks, one summary notebook, local subset inputs, and existing validators.
- Produces: local test evidence and exact HPC commands for full-data validation.

- [ ] **Step 1: Run structural and R unit tests**

Run the Python notebook-contract suite, `tests/test_source.R`, and `tests/test_extended_qc.R` with D:-local temporary directories. Expected: all pass.

- [ ] **Step 2: Run the four-section local subset checkpoint**

Run `scripts/execute_local_subset.ps1 -RegionId ALL -RunLabel local_region_notebook_report_test`. Expected: 500 cells per section, 2,000 combined cells, unchanged core-pass/readiness values, explicit `NOT_RUN_LOCAL_SUBSET` transcript status, zero deleted cells, and valid section/summary artifact contracts.

- [ ] **Step 3: Validate notebook JSON and repository path safety**

Validate all six notebook files with `render_notebooks.py`, scan tracked implementation files for disallowed C: project/output/temp paths, and run `git diff --check`. Runtime executable paths in the local PowerShell launcher are the only permitted C: references and must not be used for project data, output, cache, or temporary storage.

- [ ] **Step 4: Record unexecuted full-HPC validation exactly**

If full data, Jupyter/IRkernel, Arrow, or Bash is unavailable locally, do not claim full execution. Provide:

```bash
REGION_ID=Region_1 bash "${PIPELINE_REPO}/shell/run_notebook_qc_hpc.sh"
REGION_ID=Region_2 bash "${PIPELINE_REPO}/shell/run_notebook_qc_hpc.sh"
REGION_ID=Region_3 bash "${PIPELINE_REPO}/shell/run_notebook_qc_hpc.sh"
REGION_ID=Region_4 bash "${PIPELINE_REPO}/shell/run_notebook_qc_hpc.sh"
REGION_ID=ALL bash "${PIPELINE_REPO}/shell/run_notebook_qc_hpc.sh"
```

The final `ALL` command may reuse valid region bundles or rerun them according to the implemented runner contract, then builds the summary.

- [ ] **Step 5: Verify Git state and push**

Confirm only intended changes exist, show the commit list, then push `codex/notebook-qc-pipeline` to its configured remote. Report the local validation level and request the returned HPC `RUN_ROOT` for full-data review.
