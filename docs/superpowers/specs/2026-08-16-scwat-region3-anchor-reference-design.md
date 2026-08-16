# scWAT Region 3 anchor and gate-first downstream reference design

**Date:** 2026-08-16
**Status:** Approved design; pending written-spec review
**Scope:** scWAT only
**Branch:** `codex/notebook-qc-pipeline`

## 1. Objective

Build the downstream scWAT Xenium reference through sequential, evidence-gated admission. Region 3 is the initial anchor because it is the only section without a direct poor-cycle alarm. Regions 1 and 2 are mapped to and validated against the frozen Region 3 anchor independently. Each passing region may join an expanded consensus; each failing region becomes sensitivity-only. Region 4 is always mapping-only and sensitivity-only.

The design preserves all raw sparse counts, raw metadata, existing Phase 0-2 masks, and all 479 panel genes. It adds labels, model objects, validation results, and release decisions without modifying raw objects.

The analysis remains descriptive and exploratory. The two mice are the biological units, the sections are technical units, and no downstream cell count is treated as an additional biological replicate.

## 2. Fixed Phase 0-2 inputs

The downstream workflow starts only from a Phase 0-2 run whose completed evidence-only QC gates pass. The expected full-data inputs are:

- four validated `downstream_inputs/Region_<n>.downstream_input.rds` bundles;
- `cell_downstream_masks.tsv.gz`;
- `section_downstream_decision.tsv`;
- `gene_downstream_decision.tsv`;
- `eos_gene_decision_summary.tsv`;
- `hotspot_sensitivity_decision.tsv`; and
- `evidence_only_qc_release.tsv`.

The frozen analytical contracts are:

- Region 1: `PRIMARY_CONDITIONAL`;
- Region 2: `PRIMARY_CONDITIONAL`;
- Region 3: `PRIMARY` and initial anchor;
- Region 4: `SENSITIVITY_ONLY` and mapping-only;
- primary clustering features: 245 `PROVISIONAL_PRIMARY_FEATURES`;
- strict gene sensitivity: 67 `CONSERVATIVE_NO_SIGNAL_DETECTED` genes;
- technical-risk sensitivity-only genes: 234 `TECHNICAL_RISK_SENSITIVITY_ONLY` genes;
- raw panel: all 479 genes; and
- provisional Eosinophil signature: 53 genes comprising 4 common, 27 short-lived, and 22 long-lived genes.

No gene is described as confirmed affected or confirmed unaffected.

## 3. Reference-state machine

The pipeline uses the following ordered state machine:

1. Build and validate a Region 3 anchor.
2. Freeze the passing anchor before examining Region 1 or Region 2 admission.
3. Map and validate Region 1 against the frozen anchor.
4. Map and validate Region 2 against the same frozen anchor.
5. Admit each section independently when it passes.
6. If at least one section is admitted, attempt an expanded consensus using Region 3 plus all admitted sections.
7. Validate the expanded consensus.
8. Use the expanded consensus if it passes; otherwise revert to the frozen Region 3 anchor.
9. Freeze the selected final reference.
10. Map Region 4 to the frozen final reference.
11. Run the Eosinophil robustness analysis with Region 3 as the primary evidence, the frozen eligible reference as validation evidence, and Region 4 as a separate sensitivity-only result.

The final reference state is exactly one of:

- `EXPANDED_CONSENSUS_R1_R2_R3`;
- `EXPANDED_CONSENSUS_R1_R3`;
- `EXPANDED_CONSENSUS_R2_R3`;
- `REGION3_ANCHOR_NO_ADMITTED_REGIONS`;
- `REGION3_ANCHOR_FALLBACK`; or
- `STOP_REGION3_ANCHOR_FAILED`.

Failure of Region 1, Region 2, or an expanded consensus does not by itself stop Region 3-based downstream analysis. Only failure of the Region 3 anchor's essential integrity, annotation, or stability gates causes `STOP_REGION3_ANCHOR_FAILED`.

## 4. Notebook 03: Region 3 anchor

Create `notebooks/03_region3_anchor_reference.ipynb`.

### Inputs

- Region 3 downstream input bundle;
- slide-level mask, gene, Eosinophil, and hotspot decision tables;
- canonical scWAT marker configuration; and
- versioned analysis thresholds and random seed.

### Primary construction

Use Region 3 `primary_include` cells and the 245 provisional primary features. Fit normalization, scaling, PCA, neighbor graph, and provisional clustering only from these cells and genes. The 234 technical-risk genes must not contribute to primary PCA loadings, neighbor construction, clustering, or reference-label definition.

All 479 genes remain available for exploratory Region 3 marker inspection. Canonical markers outside the 245-gene set may support biological interpretation but may not rescue a cluster whose primary definition is unstable.

### Internal validation branches

Run the same reusable workflow for:

1. primary cells with 245 genes;
2. strict cells with 245 genes;
3. hotspot-sensitivity cells with 245 genes; and
4. primary cells with 67 conservative genes.

Match clusters across branches using shared-cell overlap and marker agreement. For every provisional anchor cluster, report size, canonical markers, branch-to-branch correspondence, label stability, hotspot dependence, and whether marker support depends primarily on technical-risk genes.

Major labels are labels representing at least 1% of Region 3 primary cells or at least 500 Region 3 primary cells. Small populations remain reportable, but their uncertainty must not be hidden by aggregate major-label metrics.

### Anchor decisions

- `PASS_ANCHOR`: essential data checks pass; every retained major cluster has canonical support; major labels are stable across the primary, strict, hotspot, and conservative branches; and no primary cluster depends mainly on technical-risk genes.
- `REVIEW_ANCHOR`: computation completes but one or more non-major populations require manual marker or morphology review. This state does not permit automatic Region 1-2 admission or final release until resolved.
- `STOP_ANCHOR`: integrity fails, a major label is unsupported or contradictory across sensitivity branches, hotspot exclusion changes the principal reference interpretation, or clustering is driven mainly by technical-risk genes.

Only `PASS_ANCHOR` permits automatic admission testing. The passing anchor object and its labels are immutable inputs to the Region 1 and Region 2 notebooks.

### Outputs

- `region3_anchor_reference.rds`;
- `region3_anchor_cell_labels.tsv.gz`;
- `region3_anchor_cluster_markers.tsv.gz`;
- `region3_anchor_cluster_stability.tsv`;
- `region3_anchor_qc_gate.tsv`; and
- a compact provenance/session manifest.

## 5. Notebooks 04 and 05: independent Region 1 and Region 2 admission

Create:

- `notebooks/04_region1_anchor_validation.ipynb`; and
- `notebooks/05_region2_anchor_validation.ipynb`.

The notebooks use identical functions, thresholds, output schemas, and seeds. Each candidate section is evaluated independently against the same frozen Region 3 anchor. Neither candidate may change anchor PCA loadings, clusters, centroids, or labels during admission testing.

### Mapping calibration

Calibrate label-confidence and reference-distance thresholds using stratified held-out Region 3 cells. Thresholds are derived globally and per major anchor label where sample size permits. They are frozen before candidate mapping and are never selected to maximize Region 1 or Region 2 assignment rates.

Every mapped candidate cell receives:

- best reference label;
- mapping confidence;
- reference distance;
- second-best label;
- confidence margin;
- calibration threshold used;
- mask status; and
- `Uncertain` when confidence or distance is insufficient.

### Admission diagnostics

For each candidate section, assess:

- mapping coverage and `Uncertain` fraction relative to held-out Region 3 calibration;
- canonical-marker agreement within each mapped cell type;
- within-label candidate-versus-Region 3 pseudobulk profile concordance;
- label stability for overlapping cells under primary and strict masks;
- compatibility of major labels under the 245-gene and 67-gene branches;
- candidate-enriched branches in an admission-only joint diagnostic embedding;
- candidate-specific diagnostic clusters;
- dependence of mapped labels or candidate-specific structure on technical-risk genes; and
- all results by cell type so poor performance in a single lineage is not concealed by a slide-wide average.

The admission-only joint embedding may diagnose section-driven structure but is never used as the anchor or final reference. A candidate-enriched population is not rejected solely because its abundance differs from Region 3. It becomes technical evidence against admission when section enrichment is accompanied by weak canonical support, poor mapping, unstable labels, technical-risk-gene dependence, or disappearance under strict/conservative sensitivity analysis.

Data-adaptive mapping thresholds come from Region 3 held-out distributions. Non-mapping stability thresholds are prespecified in the implementation configuration and printed in every output; they must not be changed after inspecting candidate outcomes without creating a new versioned run.

### Admission decisions

- `PASS_ADMIT_TO_CONSENSUS`: mapping and marker evidence are compatible with the Region 3 anchor, major labels are stable, and no unsupported section-driven cluster is introduced.
- `FAIL_SENSITIVITY_ONLY`: the section fails one or more biological/technical admission gates but its data remain available for mapped or section-local sensitivity analysis.
- `STOP_DATA_INTEGRITY`: identifiers, matrices, feature contracts, or required inputs fail structural validation. This stops processing of that candidate section but does not invalidate the frozen Region 3 anchor.

Region 1 and Region 2 decisions do not rescue or veto one another. A passing region is admitted even if the other region fails.

### Outputs per candidate

- `region<n>_anchor_mapping.tsv.gz`;
- `region<n>_celltype_concordance.tsv`;
- `region<n>_section_cluster_diagnostics.tsv`;
- `region<n>_admission_gate.tsv`; and
- a candidate-specific plots and provenance manifest.

## 6. Notebook 06: expanded consensus and fallback

Create `notebooks/06_region1_3_consensus_reference.ipynb`.

The eligible input is Region 3 plus every Region 1-2 section with `PASS_ADMIT_TO_CONSENSUS`:

| Region 1 | Region 2 | Eligible reference input |
|---|---|---|
| Pass | Pass | Regions 1, 2, and 3 |
| Pass | Fail | Regions 1 and 3 |
| Fail | Pass | Regions 2 and 3 |
| Fail | Fail | Region 3 only; no expanded fit is required |

When at least one candidate is admitted, rebuild normalization, PCA, any integration model, neighbor graph, clustering, and marker decisions using eligible sections only. Use the 245 provisional genes and primary masks. Fit the reference on a documented section-balanced training set, stratified by preliminary anchor label, so a section with more cells cannot dominate model estimation. Project all eligible primary cells back into the fitted reference for final labelling.

Any integration is restricted to eligible reference sections. Its unintegrated and integrated diagnostics must both be retained. Corrected values may support reference geometry and labels but must not replace raw counts for gene-level reporting.

### Consensus validation

Report:

- section composition of every cluster;
- Region 3 representation and marker support;
- canonical-marker agreement across admitted sections;
- primary-versus-strict label stability;
- 245-versus-67 gene compatibility;
- Region 3 hotspot sensitivity;
- unintegrated and integrated section-association diagnostics;
- section predictability from the retained latent space;
- technical-risk-gene dependence; and
- clusters absent from Region 3 with their biological and morphological evidence.

An expanded consensus passes only when it improves or preserves reference stability without creating unsupported section-driven clusters or contradictory major cell types. Clusters absent from Region 3 require canonical-marker support and stability across masks; otherwise they are merged, rejected, or held for review.

### Fallback behavior

If the expanded consensus passes, freeze it as the final reference. If it fails any essential consensus gate, record the failure and freeze the original Region 3 anchor as `REGION3_ANCHOR_FALLBACK`. The failed expanded object may be retained for audit but must not supply released labels.

Fallback is a valid primary-reference outcome, not a global analysis failure.

### Outputs

- `eligible_reference_membership.tsv`;
- `expanded_consensus_reference.rds` when attempted;
- `expanded_consensus_cell_labels.tsv.gz` when attempted;
- `consensus_cluster_markers.tsv.gz`;
- `consensus_cluster_decision.tsv`;
- `consensus_reference_gate.tsv`;
- `final_frozen_reference.rds`;
- `final_frozen_reference_labels.tsv.gz`; and
- `final_reference_decision.tsv`.

## 7. Notebook 07: Region 4 mapping-only sensitivity

Create `notebooks/07_region4_consensus_mapping.ipynb`.

Region 4 is processed only after `final_frozen_reference.rds` exists. Region 4 must not alter normalization parameters, integration anchors, feature selection, PCA loadings, neighbor graphs, clusters, marker selection, cell-type definitions, centroids, thresholds, or reference labels.

Calibrate mapping thresholds only from held-out cells belonging to the frozen eligible reference. Apply those frozen thresholds to Region 4 and assign `Uncertain` rather than forcing insufficient-confidence cells.

All Region 4 cell-type proportions, gene-level summaries, spatial summaries, and Eosinophil results are labelled `SENSITIVITY_ONLY`. Poor overall or cell-type-specific mapping produces `REGION4_MAPPING_NOT_RELEASABLE`. That state blocks Region 4 interpretation only and cannot invalidate an acceptable Region 3 anchor or expanded consensus.

Outputs:

- `region4_consensus_mapping.tsv.gz`;
- `region4_mapping_by_celltype.tsv`;
- `region4_mapping_qc_gate.tsv`; and
- a separate Region 4 sensitivity plot manifest.

## 8. Notebook 08: Eosinophil robustness

Create `notebooks/08_region3_eosinophil_robustness.ipynb`.

The evidence hierarchy is fixed:

1. Region 3 provides the primary Eosinophil evidence.
2. The frozen eligible reference provides supporting validation.
3. Region 4 provides a separate sensitivity-only result after mapping.

### Region 3 primary comparisons

Compare:

- the 53-gene provisional score with the complete 100-gene score;
- common, short-lived, and long-lived components;
- primary and strict masks;
- Region 3 hotspots retained and excluded; and
- Region 3 anchor labels with the final frozen-reference labels.

### Eligible-reference validation

For any admitted Regions 1-2, test score direction, mapped Eosinophil assignments, canonical markers, and agreement with the Region 3-defined state. Identical abundance is not required because tissue composition may differ. Candidate-section evidence supports but does not override the Region 3 primary result.

### Region 4 sensitivity

Calculate a Region 4 Eosinophil result only among adequately mapped cells and report it in separate tables and plot panels. Region 4 cannot establish, reverse, or rescue the primary Eosinophil conclusion. If Region 4 mapping is not releasable, its Eosinophil result remains descriptive and uninterpretable.

### Eosinophil release

Use `EOS_PRIMARY_RELEASED` only when score direction and cell assignments remain stable across Region 3 gene sets, cell masks, hotspot handling, and anchor/frozen-reference labels, and when the conclusion does not depend primarily on the excluded/higher-risk Eosinophil genes. Otherwise use `EOS_PRIMARY_NOT_STABLE` without stopping unrelated cell-type reference release.

## 9. Notebook 09: final downstream release summary

Create `notebooks/09_downstream_release_summary.ipynb` as the reader-facing downstream QC and release report.

It consolidates all gate tables without recomputing upstream models and reports:

- the Region 3 anchor decision;
- independent Region 1 and Region 2 admission decisions;
- the eligible expanded-consensus composition;
- expanded-consensus pass or Region 3 fallback;
- final frozen-reference identity and provenance;
- retained, merged, rejected, and review clusters;
- cell counts by section, mask, label, and reference role;
- Region 4 mapping quality and sensitivity release state;
- Eosinophil primary and sensitivity decisions; and
- unresolved limitations.

The release hierarchy is:

- `PRIMARY_REFERENCE_RELEASED_EXPANDED`;
- `PRIMARY_REFERENCE_RELEASED_REGION3_FALLBACK`;
- `PRIMARY_REFERENCE_STOP_REGION3_FAILED`;
- `EOS_PRIMARY_RELEASED` or `EOS_PRIMARY_NOT_STABLE`; and
- `REGION4_SENSITIVITY_RELEASED` or `REGION4_SENSITIVITY_NOT_RELEASABLE`.

Region 4 failure cannot stop primary-reference release. Eosinophil instability blocks the Eosinophil-state conclusion but does not stop unrelated cell-type reference release. Expanded-consensus failure triggers Region 3 fallback. Region 3 anchor failure is the only reference-level global stop.

## 10. Reusable implementation boundaries

Add reusable functions to `R/source.R`, grouped by purpose and testable without notebook execution:

- input and Phase 0-2 contract validation;
- normalization and fixed-feature preparation;
- Region 3 anchor fitting;
- sensitivity-branch fitting and cluster matching;
- held-out mapping calibration;
- candidate mapping and admission metrics;
- section-driven structure diagnostics;
- eligible-section selection;
- section-balanced consensus fitting;
- consensus validation and Region 3 fallback;
- frozen-reference mapping;
- Eosinophil scoring and stability comparisons;
- release-state aggregation; and
- provenance-aware artifact writers/reload validators.

Notebook cells orchestrate these functions and present results. Scientific calculations must not be duplicated independently across notebooks.

All stochastic steps use a versioned seed. Every output records input paths, file hashes where available, feature set, cell mask, thresholds, package versions, execution mode, run label, reference state, and upstream gate identifiers.

## 11. Testing and execution

Use failing-first tests before implementation changes. Local tests use the existing bounded subset and write only beneath the D: project. Full matrices and heavy models run on HPC.

Tests must cover:

- Region 3-only anchor construction and failure states;
- immutability of the frozen anchor during Region 1-2 mapping;
- held-out threshold calibration and `Uncertain` assignment;
- independent Region 1 and Region 2 admission;
- all four admission combinations in the eligibility table;
- section-balanced consensus fitting;
- expanded-consensus pass and Region 3 fallback;
- proof that Region 4 cannot enter reference fitting or label definition;
- Region 4 mapping failure isolated from primary release;
- separate primary, eligible-reference, and Region 4 Eosinophil evidence;
- Eosinophil instability isolated from general reference release;
- raw-count preservation and identifier alignment;
- deterministic results under the fixed seed; and
- artifact reload and provenance reconciliation.

Each notebook contains explicit input, output, checkpoint, and interpretation-limit cells. HPC execution proceeds one notebook at a time and stops at each gate for review. A later notebook must refuse to run when its required upstream gate or frozen artifact is absent.

## 12. Completion criteria

Implementation is complete when:

1. all new behavior has failing-first unit or contract tests;
2. Notebook 03 produces and reload-validates a frozen Region 3 anchor;
3. Notebooks 04 and 05 independently map and gate Regions 1 and 2;
4. Notebook 06 implements all admission combinations and verified Region 3 fallback;
5. Notebook 07 proves Region 4 is mapping-only and isolates its release decision;
6. Notebook 08 separates Region 3 primary, eligible-reference validation, and Region 4 sensitivity Eosinophil evidence;
7. Notebook 09 reports the final reference and all branch-specific release states without recomputation;
8. local bounded tests pass with limitations stated explicitly;
9. copy-and-paste HPC commands and expected outputs are documented notebook by notebook;
10. full-HPC outputs reload and reconcile before any scientific release;
11. the project progress document is updated at each checkpoint; and
12. verified changes are committed and pushed to `codex/notebook-qc-pipeline` when connectivity is available.
