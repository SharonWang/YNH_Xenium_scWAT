# scWAT Evidence-Only QC Release Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Produce immutable, evidence-only downstream masks, section/gene decisions, release gates, and per-region downstream input bundles from the four Phase 0–2 Xenium QC notebooks.

**Architecture:** Section notebooks create cell-level masks and section-local provenance without changing raw sparse objects. The slide summary joins all four sections, freezes the cross-section gene tiers and Eos set, writes the six required decision tables, and emits one downstream RDS bundle per region; downstream biological-analysis gates remain explicitly pending until PCA, clustering, label transfer, and Eos stability analyses are run.

**Tech Stack:** R 4.x, Matrix, base R TSV/gzip/RDS I/O, ggplot2, Python 3 notebook generation/contract validation, Jupyter IRkernel, PowerShell and Slurm runners.

## Global Constraints

- Write every implementation file, build artifact, test output, cache, and temporary file below `D:/Xiaonan/CODEX_projects/Yanan_Xenium` locally.
- Use `/dssg/home/acct-svetoslav_chakarov/svetoslav_chakarov/Lab_members/Yanan_Hu/YNH_Xenium` as the HPC project root.
- Preserve raw sparse counts and all 479 panel genes; masks are annotations and never destructive filters.
- Region 1 is `PRIMARY_CONDITIONAL`, Region 2 is `PRIMARY_CONDITIONAL`, Region 3 is `PRIMARY`, and Region 4 is `SENSITIVITY_ONLY`.
- Never label a gene confirmed affected or confirmed unaffected; exact cycle-to-codeword mapping remains unavailable.
- Local execution is limited to the existing bounded subset; full data are run on HPC.
- Use test-first red/green cycles for every reusable function and notebook contract change.

---

### Task 1: Cell masks and immutable section decisions

**Files:**
- Modify: `R/source.R`
- Modify: `tests/test_source.R`

**Interfaces:**
- Consumes: QC cell metadata with `qc_core_pass`, `segmentation_multiplet_flag`, `high_control_flag`, `qc_review_flag`, `region_id`, and `cell_id`; spatial hotspot table with `grid_id` and `hotspot_status`.
- Produces: `section_downstream_status(region_id) -> character(1)` and `build_cell_downstream_masks(cell_metadata, spatial_hotspots, provenance) -> data.frame`.

- [ ] **Step 1: Write failing mask tests**

```r
mask_fixture <- data.frame(
  region_id = rep("Region_3", 4), cell_id = paste0("c", 1:4), grid_id = c("g1", "g1", "g2", "g3"),
  qc_core_pass = c(TRUE, FALSE, TRUE, TRUE), segmentation_multiplet_flag = c(FALSE, FALSE, TRUE, FALSE),
  high_control_flag = c(FALSE, FALSE, FALSE, TRUE), qc_review_flag = c(FALSE, TRUE, TRUE, TRUE)
)
hotspot_fixture <- data.frame(grid_id = "g1", hotspot_status = "MORPHOLOGY_REVIEW_REQUIRED")
masks <- build_cell_downstream_masks(mask_fixture, hotspot_fixture, provenance = "unit_fixture")
stopifnot(identical(masks$primary_include, c(TRUE, FALSE, FALSE, FALSE)))
stopifnot(identical(masks$strict_include, c(TRUE, FALSE, FALSE, FALSE)))
stopifnot(identical(masks$hotspot_sensitivity_include, c(FALSE, FALSE, FALSE, FALSE)))
stopifnot(section_downstream_status("Region_4") == "SENSITIVITY_ONLY")
```

- [ ] **Step 2: Run the source test and verify RED**

Run: `Rscript tests/test_source.R`

Expected: FAIL because `build_cell_downstream_masks` and `section_downstream_status` do not exist.

- [ ] **Step 3: Implement the minimal functions**

```r
section_downstream_status <- function(region_id) {
  status <- c(Region_1="PRIMARY_CONDITIONAL", Region_2="PRIMARY_CONDITIONAL", Region_3="PRIMARY", Region_4="SENSITIVITY_ONLY")
  if (length(region_id) != 1L || !region_id %in% names(status)) stop("Unknown scWAT region.", call. = FALSE)
  unname(status[[region_id]])
}

build_cell_downstream_masks <- function(cell_metadata, spatial_hotspots, provenance) {
  required <- c("region_id","cell_id","qc_core_pass","segmentation_multiplet_flag","high_control_flag","qc_review_flag")
  if (length(setdiff(required, names(cell_metadata)))) stop("Cell metadata lacks downstream-mask fields.", call. = FALSE)
  out <- cell_metadata
  out$primary_include <- out$qc_core_pass & !out$segmentation_multiplet_flag & !out$high_control_flag
  out$strict_include <- !out$qc_review_flag
  hotspot_ids <- if (nrow(spatial_hotspots)) spatial_hotspots$grid_id[spatial_hotspots$hotspot_status == "MORPHOLOGY_REVIEW_REQUIRED"] else character()
  out$hotspot_review_cell <- out$region_id == "Region_3" & "grid_id" %in% names(out) & out$grid_id %in% hotspot_ids
  out$hotspot_sensitivity_include <- out$primary_include & !out$hotspot_review_cell
  out$section_status <- vapply(out$region_id, section_downstream_status, character(1))
  out$provenance <- provenance
  out
}
```

- [ ] **Step 4: Run `Rscript tests/test_source.R` and verify GREEN**

- [ ] **Step 5: Commit**

```powershell
git add R/source.R tests/test_source.R
git commit -m "feat: add evidence-only cell masks"
```

### Task 2: Cross-section gene tiers and Eos decisions

**Files:**
- Modify: `R/source.R`
- Modify: `tests/test_extended_qc.R`

**Interfaces:**
- Consumes: candidate evidence returned by `rank_candidate_cycle_genes`, complete panel features, and `config/eos_gene_sets.tsv`.
- Produces: `build_gene_downstream_decision(candidates, panel_features, provenance)`, `build_eos_gene_decision(gene_decision, eos_sets, provenance)`.

- [ ] **Step 1: Write failing tier tests**

```r
gene_decision <- build_gene_downstream_decision(candidate_rank, data.frame(gene=c("A","B","C","D")), "unit_fixture")
tiers <- setNames(gene_decision$gene_tier, gene_decision$gene)
stopifnot(tiers[["A"]] == "PROVISIONAL_PRIMARY_FEATURES")
stopifnot(tiers[["B"]] == "TECHNICAL_RISK_SENSITIVITY_ONLY")
stopifnot(all(gene_decision$raw_panel_tier == "RAW_COMPLETE_PANEL"))
stopifnot(!any(grepl("CONFIRMED_(AFFECTED|UNAFFECTED)", unlist(gene_decision))))
```

- [ ] **Step 2: Run `Rscript tests/test_extended_qc.R` and verify RED**

- [ ] **Step 3: Implement recurrence counting and Eos intersection**

```r
build_gene_downstream_decision <- function(candidates, panel_features, provenance) {
  genes <- unique(as.character(panel_features$gene))
  recurrence <- vapply(genes, function(g) length(unique(candidates$region_id[candidates$gene == g & candidates$section_candidate_flag & candidates$region_id %in% c("Region_1","Region_2","Region_4")])), integer(1))
  no_signal <- vapply(genes, function(g) all(candidates$counts_per_10000[candidates$gene == g] == 0, na.rm = TRUE), logical(1))
  tier <- ifelse(recurrence >= 2L, "TECHNICAL_RISK_SENSITIVITY_ONLY", "PROVISIONAL_PRIMARY_FEATURES")
  tier[no_signal] <- "CONSERVATIVE_NO_SIGNAL_DETECTED"
  data.frame(gene=genes, raw_panel_tier="RAW_COMPLETE_PANEL", gene_tier=tier,
             alarm_positive_section_count=recurrence, provenance=provenance, stringsAsFactors=FALSE)
}
```

`build_eos_gene_decision` must retain rows whose `gene_tier` is `PROVISIONAL_PRIMARY_FEATURES` or `CONSERVATIVE_NO_SIGNAL_DETECTED`, preserve `common`, `short_lived`, and `long_lived`, and report both expected and retained counts. Full-HPC validation asserts 479 total, 67 conservative, 245 provisional-primary (including the conservative subset in the primary eligibility field), 234 technical-risk, and 53 retained Eos genes partitioned 4/27/22; local-subset validation records observed counts without pretending they equal full-data counts.

- [ ] **Step 4: Run `Rscript tests/test_extended_qc.R` and verify GREEN**

- [ ] **Step 5: Commit**

```powershell
git add R/source.R tests/test_extended_qc.R
git commit -m "feat: freeze evidence-only gene tiers"
```

### Task 3: Required tables, downstream bundles, and pending release gates

**Files:**
- Modify: `R/source.R`
- Modify: `tests/test_extended_qc.R`

**Interfaces:**
- Consumes: four section bundles, masks, gene/Eos decisions, section summary, and hotspot evidence.
- Produces: the six required TSV files, `downstream_inputs/Region_1.downstream_input.rds` through Region 4, and `downstream_input_manifest.tsv`.

- [ ] **Step 1: Write failing artifact-contract tests**

```r
required <- c("cell_downstream_masks.tsv.gz","section_downstream_decision.tsv","gene_downstream_decision.tsv",
              "eos_gene_decision_summary.tsv","hotspot_sensitivity_decision.tsv","evidence_only_qc_release.tsv")
stopifnot(all(required %in% evidence_only_required_artifacts()))
paths <- write_evidence_only_qc_artifacts(extended_slide_root, evidence_fixture)
stopifnot(all(file.exists(file.path(extended_slide_root, "slide_summary", required))))
stopifnot(all(file.exists(file.path(extended_slide_root, "downstream_inputs", paste0("Region_",1:4,".downstream_input.rds")))))
stopifnot(validate_evidence_only_qc_artifacts(extended_slide_root, stop_on_error=TRUE))
```

- [ ] **Step 2: Run the extended test and verify RED**

- [ ] **Step 3: Implement writers and validators**

Each table must include `run_label`, `execution_mode`, `generated_utc`, `source_artifact`, `decision_rule`, and applicable thresholds/counts. `evidence_only_qc_release.tsv` contains one row per gate with `gate_status="PENDING_DOWNSTREAM_ANALYSIS"` for PCA dominance, cluster mask stability, 245-vs-67 cell-type agreement, Eos stability, Region 3 hotspot sensitivity, Region 4 mapping confidence, and 234-gene dependence; it must not report an overall primary PASS before these gates are evaluated.

Each RDS bundle must contain:

```r
list(
  schema_version="evidence_only_qc_v1", region_id=region_id, section_status=section_status,
  counts=section_object$counts, features=section_object$features,
  cell_metadata=merge(section_object$qc$cell_metadata, masks, by=c("region_id","cell_id"), sort=FALSE),
  gene_sets=list(provisional_primary_245=primary_genes, conservative_67=conservative_genes,
                 technical_risk_234=risk_genes, raw_complete_panel_479=all_genes, eos_provisional_53=eos_genes),
  downstream_contract=if (region_id == "Region_4") "MAP_TO_REGION_1_3_REFERENCE_WITH_UNCERTAIN" else "REFERENCE_ELIGIBILITY_FROM_CELL_MASKS",
  provenance=provenance
)
```

- [ ] **Step 4: Run both R test suites and verify GREEN**

- [ ] **Step 5: Commit**

```powershell
git add R/source.R tests/test_extended_qc.R
git commit -m "feat: write evidence-only QC release artifacts"
```

### Task 4: Update the four section notebooks

**Files:**
- Modify: `scripts/render_notebooks.py`
- Modify: `tests/test_notebook_contracts.py`
- Regenerate: `notebooks/01_section_phase0_2_QC.ipynb`
- Regenerate: `notebooks/01_section_phase0_2_QC_Region1.ipynb`
- Regenerate: `notebooks/01_section_phase0_2_QC_Region2.ipynb`
- Regenerate: `notebooks/01_section_phase0_2_QC_Region3.ipynb`
- Regenerate: `notebooks/01_section_phase0_2_QC_Region4.ipynb`

**Interfaces:**
- Consumes: Task 1 functions and the existing section QC outputs.
- Produces: section-local `cell_downstream_masks.tsv.gz`, `section_downstream_decision.tsv`, and mask-count plots/tables.

- [ ] **Step 1: Add failing notebook contract assertions**

```python
for token in ["build_cell_downstream_masks", "cell_downstream_masks.tsv.gz", "section_downstream_decision.tsv",
              "PRIMARY_CONDITIONAL", "SENSITIVITY_ONLY", "raw objects are not modified"]:
    self.assertIn(token, notebook_text)
```

- [ ] **Step 2: Run `python tests/test_notebook_contracts.py` and verify RED**

- [ ] **Step 3: Add notebook cells and regenerate notebooks**

The generated section code calls `build_cell_downstream_masks(spatial_cells, spatial_hotspots, provenance)`, writes the two section-local tables under `${RUN_ROOT}/sections/<Region_ID>/`, displays counts for all three masks, and states that Region 3 hotspot exclusion is sensitivity-only and Region 4 cannot enter reference discovery.

Run:

```powershell
$env:PYTHONDONTWRITEBYTECODE='1'
python scripts/render_notebooks.py --all
```

- [ ] **Step 4: Run notebook contract tests and verify GREEN**

- [ ] **Step 5: Commit**

```powershell
git add scripts/render_notebooks.py tests/test_notebook_contracts.py notebooks
git commit -m "feat: add evidence-only masks to section notebooks"
```

### Task 5: Update slide summary and reader-facing evidence

**Files:**
- Modify: `scripts/render_notebooks.py`
- Modify: `tests/test_notebook_contracts.py`
- Regenerate: `notebooks/02_slide_QC_summary.ipynb`
- Modify: `README.md`

**Interfaces:**
- Consumes: all four completed section bundles.
- Produces: combined masks, fixed decisions, six required outputs, four downstream bundles, Cell-style supporting plots, and reload checks.

- [ ] **Step 1: Add failing summary contract assertions**

```python
for token in ["PROVISIONAL_PRIMARY_FEATURES", "CONSERVATIVE_NO_SIGNAL_DETECTED",
              "TECHNICAL_RISK_SENSITIVITY_ONLY", "RAW_COMPLETE_PANEL", "PENDING_DOWNSTREAM_ANALYSIS",
              "downstream_input_manifest.tsv", "Uncertain"]:
    self.assertIn(token, summary_text)
```

- [ ] **Step 2: Run notebook contract tests and verify RED**

- [ ] **Step 3: Replace obsolete HOLD/10x-required narrative and add release cells**

The summary notebook must show: fixed section decisions; per-section mask counts; gene-tier and Eos counts; Region 3 hotspot sensitivity counts; the pending downstream release-gate table; and output/reload validation. It must clearly distinguish current QC preparation from later PCA, cluster validation, Eos robustness, and Region 4 mapping, and must not claim those pending analyses passed.

- [ ] **Step 4: Regenerate notebooks, run contract tests, and verify GREEN**

- [ ] **Step 5: Commit**

```powershell
git add scripts/render_notebooks.py tests/test_notebook_contracts.py notebooks/02_slide_QC_summary.ipynb README.md
git commit -m "feat: report evidence-only downstream release"
```

### Task 6: Local subset execution, discrepancy report, and HPC commands

**Files:**
- Modify: `scripts/execute_local_subset.ps1`
- Modify: `shell/run_notebook_qc_hpc.sh`
- Modify: `slurm/scwat_notebook_qc.sbatch`
- Modify: `PROJECT_RESEARCH_PLAN_AND_PROGRESS.md` in the project root
- Create during validation only below D: `adipose_analysis/scwat_qc_outputs/evidence_only_subset_test/`

**Interfaces:**
- Consumes: completed Tasks 1–5.
- Produces: validated local subset artifacts and exact one-region-at-a-time HPC commands.

- [ ] **Step 1: Add failing runner assertions for the six outputs and four bundles**

The local runner reload check must assert all six slide files, four bundle files, four unique section decisions, no `TRUE` primary-release result while gates are pending, and zero paths beginning with `C:`.

- [ ] **Step 2: Run the local runner check and verify RED**

- [ ] **Step 3: Implement the runner checks and HPC validation command**

```bash
for region in Region_1 Region_2 Region_3 Region_4; do
  REGION_ID=${region} RUN_LABEL=evidence_only_qc_v1 bash "${PIPELINE_REPO}/shell/run_notebook_qc_hpc.sh"
done
REGION_ID=SUMMARY RUN_LABEL=evidence_only_qc_v1 bash "${PIPELINE_REPO}/shell/run_notebook_qc_hpc.sh"
Rscript -e "source(file.path(Sys.getenv('PIPELINE_REPO'),'R','source.R')); stopifnot(validate_evidence_only_qc_artifacts(file.path(Sys.getenv('PROJECT_ROOT'),'adipose_analysis','scwat_qc_outputs','evidence_only_qc_v1'), stop_on_error=TRUE))"
```

- [ ] **Step 4: Run all validation below D:**

```powershell
$env:TEMP='D:\Xiaonan\CODEX_projects\Yanan_Xenium\adipose_analysis\tmp'
$env:TMP=$env:TEMP
$env:R_USER='D:\Xiaonan\CODEX_projects\Yanan_Xenium\adipose_analysis\r_user'
$env:PYTHONDONTWRITEBYTECODE='1'
Rscript tests/test_source.R
Rscript tests/test_extended_qc.R
python tests/test_notebook_contracts.py
powershell -ExecutionPolicy Bypass -File scripts/execute_local_subset.ps1 -RegionId ALL -RunLabel evidence_only_subset_test
```

Expected: both R suites pass; notebook contracts pass; all five notebooks execute without error; local mode records full-data count checks as expected/not-applicable rather than fabricating 479/245/234/67/53 counts.

- [ ] **Step 5: Update progress with exact test results and discrepancies, then commit**

```powershell
git add scripts/execute_local_subset.ps1 shell/run_notebook_qc_hpc.sh slurm/scwat_notebook_qc.sbatch README.md
git commit -m "test: validate evidence-only QC subset workflow"
```

### Task 7: Final verification and handoff

**Files:**
- Verify all modified files.

**Interfaces:**
- Consumes: all completed implementation tasks.
- Produces: evidence-backed completion report and HPC commands; no push unless explicitly requested.

- [ ] **Step 1: Run `git diff --check`, both R tests, Python contract tests, and the bounded local execution once more**
- [ ] **Step 2: Inspect required TSV schemas and each RDS bundle’s dimensions, masks, gene vectors, and Region 4 mapping contract**
- [ ] **Step 3: Confirm `git status --short` contains no untracked test/temp outputs and scan tracked text for prohibited `C:` project paths**
- [ ] **Step 4: Report tests, expected local-vs-HPC discrepancies, output locations, and exact HPC order before any biological analysis**

