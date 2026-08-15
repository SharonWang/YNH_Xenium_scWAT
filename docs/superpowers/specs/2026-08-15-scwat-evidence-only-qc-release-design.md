# scWAT evidence-only QC release design

**Date:** 2026-08-15  
**Status:** Approved scientific decisions; implementation pending written-spec review  
**Scope:** scWAT only  
**Branch:** `codex/notebook-qc-pipeline`

## 1. Objective

Convert the completed Phase 0-2 and extended Xenium QC evidence into reproducible downstream-use masks, gene tiers, section-use decisions, sensitivity inputs, and release gates without relying on confirmed cycle-to-codeword mapping from 10x Support.

The implementation must preserve raw sparse counts and raw metadata. It writes new analysis masks and downstream input objects; it never deletes cells, genes, transcripts, or source files.

This phase prepares inputs for downstream normalization, PCA, integration, clustering, reference construction, Region 4 mapping, and Eosinophil-state analysis. It does not claim that those downstream analyses have already passed their release checks.

## 2. Fixed section decisions

The section statuses are fixed by the approved evidence-only decision:

| Region | Section | Mouse | Status | Permitted use |
|---|---|---|---|---|
| Region 1 | 62308 | Mouse 1 | `PRIMARY_CONDITIONAL` | Reference construction with provisional primary genes and approved cell masks |
| Region 2 | 62309 | Mouse 1 | `PRIMARY_CONDITIONAL` | Reference construction with provisional primary genes and approved cell masks |
| Region 3 | 62310 | Mouse 2 | `PRIMARY` | Reference construction, complete-panel exploratory analysis, hotspot sensitivity analysis, and Eos robustness analysis |
| Region 4 | 62311 | Mouse 2 | `SENSITIVITY_ONLY` | Mapping to the finalized Region 1-3 reference and sensitivity-only summaries |

Region 4 must not contribute to PCA loading estimation, integration anchors, cluster discovery, marker definition, reference centroids, or primary gene-level results. Downstream mapping must allow `Uncertain` rather than force a low-confidence Region 4 cell into a reference label.

No complete section is removed from storage or from descriptive QC reporting.

## 3. Cell masks

Reusable functions create three logical columns for every cell:

```r
primary_include <- qc_core_pass &
  !segmentation_multiplet_flag &
  !high_control_flag

strict_include <- !qc_review_flag

hotspot_sensitivity_include <- primary_include &
  !(region_id == "Region_3" & morphology_review_hotspot)
```

The masks are nested where expected: `strict_include` must be a subset of `primary_include`. `hotspot_sensitivity_include` differs from `primary_include` only for Region 3 cells assigned to an FDR-positive morphology-review hotspot.

On the audited full dataset, the primary mask is expected to exclude 6,215 unique cells and retain 279,963 cells. The strict mask is expected to exclude 13,257 and retain 272,921. The Region 3 hotspot sensitivity is expected to remove 1,162 cells in 34 bins from the Region 3 primary mask. These are reconciliation expectations for the audited full run, not hard-coded filters. The functions derive counts from the current artifacts and fail if identifiers or hotspot joins do not reconcile.

The per-section output is:

- `cell_downstream_masks.tsv.gz`

It contains cell ID, region, mouse, section, all original QC flags, the three masks, hotspot/bin status, mask reasons, execution mode, run label, and source artifact paths.

## 4. Gene evidence and use tiers

Gene tiers are derived from the section-specific depletion and/or Q20-loss candidate flag in the three alarm-positive regions. No gene is described as confirmed affected or confirmed unaffected.

For each gene, calculate `alarm_positive_sections_flagged` across Regions 1, 2, and 4 and write non-mutually-exclusive use fields:

- `raw_panel_status = RAW_COMPLETE_PANEL` for all 479 genes;
- `conservative_evidence_status = CONSERVATIVE_NO_SIGNAL_DETECTED` when flagged in zero alarm-positive sections;
- `primary_feature_status = PROVISIONAL_PRIMARY_FEATURES` when flagged in zero or one alarm-positive section;
- `technical_risk_status = TECHNICAL_RISK_SENSITIVITY_ONLY` when flagged in at least two alarm-positive sections.

The audited full run must reconcile:

- 479 `RAW_COMPLETE_PANEL` genes;
- 67 `CONSERVATIVE_NO_SIGNAL_DETECTED` genes;
- 245 `PROVISIONAL_PRIMARY_FEATURES` genes;
- 234 `TECHNICAL_RISK_SENSITIVITY_ONLY` genes.

The 67 conservative genes are a subset of the 245 provisional primary features. The 234 technical-risk genes must not define primary PCA loadings, integration anchors, clusters, or reference labels, but remain available in the raw object and sensitivity outputs.

The slide-level output is:

- `gene_downstream_decision.tsv`

It contains the raw section-level metrics and flags, recurrence count, all terminology fields, Eos membership, thresholds, execution mode, run label, comparison caveats, and provenance paths.

## 5. Eosinophil gene coverage and robustness contract

The provisional primary Eos signature is the intersection of the 100-gene Eos list with `PROVISIONAL_PRIMARY_FEATURES`. The audited full run must reproduce:

- 4 common genes;
- 27 short-lived genes;
- 22 long-lived genes;
- 53 unique retained Eos genes.

The complete 100-gene set remains available for exploratory Region 3 analysis. The downstream Eos workflow must compare within Region 3:

1. 53-gene provisional scores versus complete 100-gene scores;
2. hotspot cells retained versus excluded;
3. primary versus strict cell masks; and
4. cell assignments and score direction across all combinations.

An Eos conclusion is releasable only when direction and cell assignments are stable. This QC phase writes the gene membership and downstream comparison contract. Before Eos scores exist, the corresponding release checks are `PENDING_DOWNSTREAM_ANALYSIS`, never `PASS`.

The output is:

- `eos_gene_decision_summary.tsv`

## 6. Cluster and reference validation contract

The downstream reference is discovered using Regions 1-3, `primary_include` cells, and the 245 provisional primary features. Region 3 and canonical markers must support every retained cluster.

Downstream validators must flag, merge, or reject a cluster that is:

- predominantly driven by one alarm-positive section;
- unsupported by canonical markers;
- unstable under the strict cell mask;
- inconsistent with the 67-gene conservative sensitivity analysis; or
- absent from Region 3 without a documented biological or morphological explanation.

The QC notebooks create the inputs and a versioned validation contract. They do not execute PCA, integration, clustering, or marker discovery in Phase 0-2. Their release table marks these checks `PENDING_DOWNSTREAM_ANALYSIS` until a later downstream notebook supplies the required results.

## 7. Region 4 mapping contract

Region 4 is mapped only after the Region 1-3 reference is finalized. Mapping uses the same 245 provisional features and reports, for every Region 4 cell:

- assigned reference label;
- mapping confidence;
- reference distance;
- second-best label and confidence margin;
- cell mask used; and
- `Uncertain` when the prespecified confidence/distance criteria are not met.

Threshold values must be trained or calibrated only from held-out/reference Region 1-3 cells and recorded in the output. They must not be chosen to maximize Region 4 assignment rate. Region 4 cell-type proportions and gene-level expression remain sensitivity-only regardless of mapping confidence.

This QC phase writes a mapping-input bundle and the release-check placeholder; it does not fabricate mapping outputs before the downstream reference exists.

## 8. Hotspot sensitivity decision

Region 3 hotspot cells remain in the primary analysis. The hotspot sensitivity mask excludes the cells falling in FDR-positive `MORPHOLOGY_REVIEW_REQUIRED` bins.

The output:

- `hotspot_sensitivity_decision.tsv`

records region, grid/bin identifiers, cell counts, flagged counts, review rate, FDR/effect statistics, bounding boxes, primary handling (`KEEP_IN_PRIMARY`), sensitivity handling (`EXCLUDE_IN_HOTSPOT_SENSITIVITY`), thresholds, execution mode, run label, and provenance.

The audited full run must reproduce 34 Region 3 bins and 1,162 cells. Local subset mode may produce no FDR-positive bins and must label the full-data reconciliation `NOT_ESTIMABLE_LOCAL_SUBSET`.

## 9. Release gates

The output:

- `evidence_only_qc_release.tsv`

contains one row per gate with `gate_id`, scope, status, measured value, threshold, evidence source, interpretation, run label, execution mode, and next required artifact.

Primary release automatically stops if any completed downstream check demonstrates:

- PCA or clustering dominated by section identity;
- unstable major cell-type assignments between primary and strict masks;
- contradictory major cell types between the 245-gene and 67-gene analyses;
- reversal or instability of the Eos-state conclusion;
- a changed main conclusion after Region 3 hotspot exclusion;
- poor or strongly cell-type-dependent Region 4 mapping confidence; or
- a primary conclusion driven mainly by the 234 technical-risk genes.

Gate statuses are `PASS`, `STOP`, or `PENDING_DOWNSTREAM_ANALYSIS`. Overall release is `STOP` if any gate is `STOP`, otherwise `PENDING_DOWNSTREAM_ANALYSIS` while any required downstream gate is pending, otherwise `PASS`.

QC-only gates verify section decisions, mask reconciliation, gene-tier reconciliation, Eos coverage, preservation of raw counts, and absence of Region 4 from the primary reference definition.

## 10. Downstream input bundles

Each region notebook writes its section-local cell masks and preserves the existing raw-count RDS. Because the gene recurrence tiers require evidence from all four section bundles, the summary notebook performs the final slide-level join and writes:

```text
${RUN_ROOT}/downstream_inputs/Region_1.downstream_input.rds
${RUN_ROOT}/downstream_inputs/Region_2.downstream_input.rds
${RUN_ROOT}/downstream_inputs/Region_3.downstream_input.rds
${RUN_ROOT}/downstream_inputs/Region_4.downstream_input.rds
${RUN_ROOT}/downstream_inputs/downstream_input_manifest.tsv
```

Each RDS contains or references:

- the unchanged raw sparse count matrix and feature metadata;
- cell metadata with the three masks;
- fixed section status and permitted-use fields;
- the shared 479-gene downstream decision table;
- the 245-gene primary feature vector;
- the 67-gene conservative sensitivity vector;
- the complete 479-gene vector;
- the 53-gene provisional Eos membership;
- source paths, file hashes where already available, run label, execution mode, and session information.

The summary notebook owns final bundle construction because gene tiers are slide-level evidence. Each bundle records the individual region notebook and section RDS from which it was assembled. No normalized matrix, PCA, cluster label, or Region 4 mapped label is invented at this stage.

## 11. Notebook presentation

All four named region notebooks gain a bounded **Evidence-only downstream cell masks** section showing:

- inputs and exact output path;
- section decision;
- primary/strict/hotspot-sensitivity counts;
- mask-overlap checks;
- a Cell-style mask-count plot; and
- artifact reload validation.

`notebooks/02_slide_QC_summary.ipynb` gains a concise **Evidence-only downstream release** section showing:

- section-use decision table and plot;
- combined cell-mask counts and overlap checks;
- gene-tier counts and plot;
- Eos coverage by gene set;
- Region 3 hotspot sensitivity summary;
- downstream input manifest;
- release-gate table distinguishing completed QC evidence from pending downstream validation; and
- the explicit Region 1-3 reference / Region 4 mapping boundary.

Every chunk states inputs, outputs, mode, thresholds, and interpretation limits.

## 12. Required outputs

The exact required filenames are:

- section-level `cell_downstream_masks.tsv.gz` in each section output directory;
- slide-level `section_downstream_decision.tsv`;
- slide-level `gene_downstream_decision.tsv`;
- slide-level `eos_gene_decision_summary.tsv`;
- slide-level `hotspot_sensitivity_decision.tsv`;
- slide-level `evidence_only_qc_release.tsv`.

Additional required handoff files are the four downstream RDS bundles and `downstream_input_manifest.tsv`.

## 13. Testing and execution modes

Use test-driven development for every reusable function and artifact contract.

Pure-R fixtures verify:

- exact mask logic, nesting, unique counts, and reasons;
- section decision mapping and Region 4 restrictions;
- zero/one versus at-least-two alarm-section gene recurrence;
- non-mutually-exclusive 67/245/234/479 terminology;
- Eos intersection counts;
- hotspot joins and Region 3-only sensitivity behavior;
- release status aggregation with `PASS`, `STOP`, and pending gates;
- preservation of raw sparse dimensions/counts;
- rejection of duplicate/missing cell or gene identifiers; and
- downstream bundle contents and Region 4 reference exclusions.

`LOCAL_SUBSET` executes all mask/bundle code but labels full-data gene/hotspot reconciliation checks as `NOT_ESTIMABLE_LOCAL_SUBSET` when full transcript evidence or hotspots are absent. Fixture tests still verify the full logic.

`FULL_HPC` requires four validated full section bundles and must reproduce the audited 67/245/234/479 gene counts, 53 Eos genes with 4/27/22 membership, 34 Region 3 hotspot bins, 1,162 hotspot cells, and mask counts unless a documented new run changes the evidence. Discrepancies stop downstream-input release and are reported rather than overwritten.

## 14. Completion criteria

Implementation is complete when:

1. failing-first tests cover every new reusable behavior;
2. all four region notebooks write and reload `cell_downstream_masks.tsv.gz`;
3. the summary writes and reloads all five remaining required TSV outputs;
4. four final downstream input bundles preserve raw counts and carry the approved masks/tier vectors;
5. Region 4 is absent from the primary-reference definition and marked sensitivity-only;
6. local subset execution completes for all four region notebooks and the summary with explicit full-data limitations;
7. exact full-HPC commands and expected reconciliation checks are documented;
8. the project progress record is updated; and
9. verified changes are committed to `codex/notebook-qc-pipeline` and pushed when GitHub connectivity is available.
