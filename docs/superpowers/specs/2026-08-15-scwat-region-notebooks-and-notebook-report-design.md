# scWAT region notebooks and notebook-only QC report design

**Date:** 2026-08-15  
**Status:** Approved design; implementation pending written-spec review  
**Scope:** scWAT only  
**Branch:** `codex/notebook-qc-pipeline`

## 1. Objective

Make the Phase 0-2 Xenium QC workflow match the way it is reviewed and run on HPC:

1. keep one committed Jupyter notebook for each of the four scWAT sections;
2. run and inspect each section independently;
3. use `notebooks/02_slide_QC_summary.ipynb` as the only slide-level QC summary and acceptance report;
4. keep reusable calculations in `R/source.R` so the region notebooks remain consistent;
5. remove the separate generated R/HTML acceptance-report package; and
6. preserve explicit local-subset versus full-HPC execution and validation states.

The scientific conclusions remain unchanged: Regions 1, 2, and 4 have direct poor-quality-cycle alarms; exact cycle/channel/codeword identity requires 10x diagnostics; affected-gene lists remain candidate evidence; and morphology labels require visual review.

## 2. Approved architecture

The repository retains `notebooks/01_section_phase0_2_QC.ipynb` as the canonical reusable template. A deterministic generator creates four committed, parameter-locked copies:

- `notebooks/01_section_phase0_2_QC_Region1.ipynb`
- `notebooks/01_section_phase0_2_QC_Region2.ipynb`
- `notebooks/01_section_phase0_2_QC_Region3.ipynb`
- `notebooks/01_section_phase0_2_QC_Region4.ipynb`

Each copy uses the same functions and cell structure but fixes its own `SECTION_ID`. The region notebooks may be executed one by one in Jupyter on HPC. Generated notebooks are committed so the user can open, run, and inspect the exact file for each region without manually changing parameters.

`notebooks/02_slide_QC_summary.ipynb` becomes the sole reader-facing slide report. It reads the four completed section artifact bundles, validates their contract, computes the slide-level diagnostics, displays concise tables and figures inline, and writes machine-readable summary artifacts below the configured run root. It does not depend on `reports/2026-08-15_full_hpc_qc_acceptance/build_report.R` or `report.html`.

## 3. Region notebook contract

Every region notebook contains the combined Phase 0-2 workflow and displays, in order:

1. project-root, environment, section, and run-label checks;
2. metadata contract and section identity;
3. provenance, file integrity, panel reconciliation, and QC gate;
4. Xenium import and cell-level QC without automatic cell deletion;
5. direct alarm evidence and explicit unresolved-cycle wording;
6. bounded per-gene transcript-quality summaries;
7. spatial edge, density, clustering, and hotspot diagnostics;
8. section-level QC tables and cell-style plots;
9. artifact validation and reload checks; and
10. an execution-status block identifying `LOCAL_SUBSET`, `FULL_HPC`, or an explicit skipped HPC-only diagnostic.

Every analytical chunk states its inputs, outputs, and output paths. Full-data transcript aggregation remains lazy and bounded; the complete transcript Parquet table must never be collected into R memory.

## 4. Slide-summary notebook contract

The summary notebook is intentionally simple to review and follows this reader-facing structure:

1. **TL;DR and decision table** — section status, direct alarms, review burden, and required next action.
2. **Inputs and validation** — run root, four expected section bundles, mouse mapping, row/identifier reconciliation, and artifact-contract checks.
3. **Core QC distributions** — cell yield, core-pass rate, review rate, transcript/gene distributions, controls, and cell area by section.
4. **Alarm evidence and candidate genes** — direct available alarms; candidate affected-gene tables and plots clearly labelled `CANDIDATE_NOT_CONFIRMED`; exact cycle identity marked as requiring 10x diagnostics.
5. **Subset versus full-data burden** — full-data ranking, subset reference ranking, rank changes, and descriptive agreement.
6. **Spatial QC** — edge and density effects, global clustering, hotspot counts and bounding boxes, and morphology-review status without automatic fold/tear labels.
7. **Within-mouse concordance** — 62308/62309 for Mouse 1 and 62310/62311 for Mouse 2, including the Region 3/4 cell-area difference.
8. **Final QC decision and next actions** — explicit section-level gates, unresolved review items, and criteria for release to downstream analysis.

The notebook writes its combined TSV/PNG/RDS artifacts to `${RUN_ROOT}/slide_summary/`. Its executed cells are the report; no duplicate HTML report is built or maintained.

## 5. Metadata and scientific interpretation

The fixed provisional metadata used until the user provides the final metadata is:

| Section | Xenium region | Mouse | Condition | Age |
|---|---|---|---|---|
| Region 1 | 62308 | Mouse 1 | WT, untreated | 8 weeks |
| Region 2 | 62309 | Mouse 1 | WT, untreated | 8 weeks |
| Region 3 | 62310 | Mouse 2 | WT, untreated | 8 weeks |
| Region 4 | 62311 | Mouse 2 | WT, untreated | 8 weeks |

Left/right labels are not used. Sections are technical tissue sections, not independent biological replicates. Mouse is the biological unit for later inference.

Current QC interpretation carried into the summary notebook:

- Regions 1, 2, and 4 remain on hold for 10x diagnostic clarification of the poor-quality-cycle alarms.
- Region 3 is a conditional pass pending morphology review of the 34 candidate hotspot bins.
- Region 4 additionally requires review of the large cell-area difference relative to Region 3.
- Slide-wide comparative gene-level conclusions remain on hold until the affected-gene decision is frozen.

## 6. Runner and documentation behavior

The HPC runner and README list the four named notebooks in Region 1 to Region 4 order, followed by `02_slide_QC_summary.ipynb`. The runner supports independent execution of any one region and sequential execution of all four. It checks the HPC project root:

`/dssg/home/acct-svetoslav_chakarov/svetoslav_chakarov/Lab_members/Yanan_Hu/YNH_Xenium`

It installs nothing, checks required R packages and inputs, and keeps notebooks, outputs, logs, caches, and temporary files below that root. Local validation uses only the deterministic subset under the D: project.

## 7. Removal scope

Remove the separate acceptance-report implementation and generated artifacts:

- `reports/2026-08-15_full_hpc_qc_acceptance/build_report.R`
- `reports/2026-08-15_full_hpc_qc_acceptance/report.html`
- report-only `artifact.json`
- the report-local duplicate candidate TSV

No full-HPC source outputs are deleted. The returned HPC run remains the evidence source and is read by the summary notebook when configured as its run root.

## 8. Testing and validation

Use test-driven development for generator, runner, and notebook-contract behavior.

Automated checks must verify:

- exactly four named region notebooks exist and each fixes the correct `SECTION_ID`;
- all four region notebooks share the same source template fingerprint apart from the approved parameter cell and title;
- the summary notebook contains the required report sections and has no dependency on `build_report.R` or `report.html`;
- notebook JSON and the R-kernel metadata are valid;
- no notebook, runner, or configuration contains a `C:` project/output/temp path;
- all generated outputs resolve under the supplied project/run root;
- local subset execution passes for all four region notebooks and the summary when dependencies permit;
- when local execution cannot exercise a full-data path, the exact HPC command and expected artifact contract are reported;
- current full-HPC counts and statuses reload consistently without recomputing heavy data locally; and
- the repository progress record is updated with the restructuring, validation level, and outstanding QC gates.

Notebook execution is considered fully validated only after the four region notebooks and summary notebook run top-to-bottom on HPC. Local structural tests and subset runs do not claim that full-data diagnostics passed.

## 9. Completion criteria

The change is complete when:

1. the four region notebooks are committed and independently runnable;
2. `02_slide_QC_summary.ipynb` is the sole QC summary/report notebook;
3. the separate R/HTML acceptance-report package is removed;
4. README and runners use the named region notebooks and summary order;
5. reusable logic remains in `R/source.R` with no divergent per-region calculations;
6. all available local structural/unit/subset checks pass, with full-data gaps explicitly delegated to HPC;
7. the project progress Markdown records the change and outstanding QC work; and
8. changes are committed to `codex/notebook-qc-pipeline` and pushed only after verification.
