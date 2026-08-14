# scWAT Xenium extended QC diagnostics design

**Date:** 2026-08-14  
**Status:** Approved design; implementation pending written-spec review  
**Scope:** scWAT only  
**Pipeline shape:** one combined Phase 0-2 notebook run once per section, followed by one slide-level QC summary notebook

## 1. Objective

Extend the existing two-notebook QC workflow to answer four technical questions on the complete HPC dataset:

1. Which alarm evidence is directly available, and which panel genes are plausible candidates for effects from poor imaging cycles?
2. Does the full-data section ranking of QC review burden reproduce the validated subset ranking?
3. Are QC review flags spatially enriched at tissue edges, in dense cell aggregates, or in statistically defined hotspots that merit morphology review?
4. Are the two technical sections from each mouse concordant after full-data QC?

The extension remains diagnostic. It must not convert candidate genes into confirmed cycle-affected genes, automatically label spatial hotspots as folds or tears, treat cells or sections as biological replicates, or clear the existing Xenium `HOLD` gates.

## 2. Approved approach

Add reusable functions to `R/source.R`, additional code chunks to both existing R-kernel notebooks, new tests, and matching local/HPC runner checks. The section notebook computes section-local summaries and spatial diagnostics. The slide notebook combines all four section bundles, ranks candidate genes, compares subset and full-data results, and evaluates within-mouse technical concordance.

The local subset remains a pipeline smoke test. Transcript-Parquet and morphology-heavy work is HPC-first. A subset run may emit explicit `NOT_RUN_LOCAL_SUBSET` records when large inputs were intentionally excluded. A full-data run must fail clearly when a required input or package is absent.

## 3. Inputs

### 3.1 Existing required inputs per section

- `experiment.xenium`
- `analysis_summary.html`
- `metrics_summary.csv`
- `gene_panel.json`
- `cells.csv.gz` or the corresponding cell table
- cell-feature sparse matrix
- verified `config/scwat_sample_manifest.tsv`

### 3.2 Extended full-data inputs

- `transcripts.parquet`, read by projected columns and summarized without collecting the complete table;
- `morphology.ome.tif`, used only for optional downsampled review panels when a supported reader is available;
- cell centroids and segmentation metadata already imported during Phase 2;
- subset reference table containing the validated review rates and fixed subset order `Region_1 > Region_3 > Region_2 > Region_4`.

### 3.3 Runtime requirements

- Existing R requirements: `Matrix`, `jsonlite`, and `ggplot2`.
- Full transcript diagnostics: R package `arrow`.
- Scalable nearest-neighbour calculations: R package `RANN`.
- Optional morphology rendering: a preinstalled OME-TIFF-capable reader. Its absence does not prevent coordinate-based hotspot diagnostics, but it produces an explicit manual-review manifest instead of an image overlay.

The pipeline installs nothing. Local and HPC preflight checks report missing requirements before analysis begins. All temporary and output paths remain below the supplied project root.

## 4. Section-notebook extensions

### 4.1 Direct alarm evidence

Preserve the alarm ID, level, title, full message, raised state, section identifier, and source filename in `cycle_alarm_evidence.tsv`.

Add an explicit interpretation field:

- `DIRECT_EVIDENCE`: the Xenium alarm is present in the supplied output;
- `NO_ALARM_REPORTED`: no corresponding alarm is present;
- `CYCLE_IDENTITY_UNRESOLVED_REQUIRES_10X`: standard outputs do not identify the exact poor cycle;
- `GENE_EFFECT_UNCONFIRMED`: any downstream gene list is a candidate diagnostic list, not a confirmed mapping.

### 4.2 Per-gene transcript-quality summary

Project only the required Parquet columns after validating the actual schema. Accept documented aliases through one centralized schema resolver; fail on ambiguous or missing gene/QV fields in full mode.

For every panel gene and section, write `gene_transcript_quality.tsv` containing:

- gene name and panel gene-set membership;
- raw gene-expression counts from the sparse cell matrix;
- counts per 10,000 gene-expression transcripts;
- number and fraction of cells with at least one transcript;
- transcript rows observed in `transcripts.parquet`;
- mean transcript QV;
- fraction of transcripts with QV at least 20;
- the number of represented codeword indices when available;
- diagnostic execution status and any missing-input reason.

No per-transcript table is written to the QC output. Arrow aggregation must remain lazy/streamed until the bounded gene-level summary is collected.

### 4.3 Spatial QC diagnostics

Use section centroids and the binary `qc_review_flag` without loading a full pairwise distance matrix.

- Estimate local cell density from the configurable `k`th-nearest-neighbour distance using `RANN`.
- Define dense aggregates operationally as the highest configured local-density quantile, defaulting to the top 10%. This is a density class, not a cell-type annotation.
- Rasterize centroids onto a configurable spatial grid. Occupied bins adjacent to empty bins form the computational edge proxy. Report the grid size and edge definition in every output.
- Compare review rates at edge versus interior and dense versus non-dense locations using counts, proportions, absolute differences, and continuity-corrected risk ratios.
- Test global spatial clustering of the binary review flag with a k-nearest-neighbour Moran-style statistic and within-section label permutations.
- Identify review hotspots on the spatial grid using empirical permutation probabilities, Benjamini-Hochberg adjustment, minimum cell-count rules, and effect-size reporting.

Write:

- `spatial_qc_global.tsv`;
- `spatial_qc_edge_density.tsv`;
- `spatial_qc_hotspots.tsv`;
- `spatial_qc_cell_annotations.tsv.gz` with edge, density, and hotspot annotations;
- `spatial_manual_review_manifest.tsv` containing hotspot bounding boxes and the matching morphology/Xenium Explorer source path;
- Cell-inspired spatial QC figures showing all cells, review flags, edge class, density, and hotspot boxes.

Hotspots are labelled `MORPHOLOGY_REVIEW_REQUIRED`. The pipeline must never infer “fold” or “tear” from coordinates alone.

## 5. Slide-summary extensions

### 5.1 Candidate affected-gene ranking

Region 3 is the only alarm-negative section in the current slide. Candidate ranking therefore uses transparent evidence tiers rather than a formal differential-expression test.

For each gene and alarm-positive section, compute abundance and transcript-quality contrasts against Region 3. Preserve the distinction between the paired comparison within Mouse 2 (`Region_4` versus `Region_3`) and cross-mouse comparisons (`Region_1`/`Region_2` versus `Region_3`).

Configurable advisory thresholds are applied to absolute effect sizes and empirical within-panel ranks. The default candidate logic requires depletion and/or quality degradation; thresholds and intermediate values are written to the output so candidates can be reproduced and sensitivity-tested.

- **Tier A:** concordant abundance depletion and poorer QV/Q20 evidence in the paired Region 4 versus Region 3 comparison.
- **Tier B:** concordant abundance/quality evidence recurring in at least two alarm-positive sections relative to Region 3.
- **Tier C:** depletion-only or quality-only exploratory evidence.
- **Unranked:** insufficient transcript-quality data or no advisory evidence.

Write `candidate_cycle_affected_genes.tsv` for all panel genes, not only flagged genes. Include raw metrics, contrasts, empirical ranks, evidence tier, comparison type, caveats, and the statement `CANDIDATE_NOT_CONFIRMED`. Write a focused heatmap and ranked-effect plot. The report must state that exact cycle identity and definitive cycle-to-gene mapping require 10x diagnostic information.

### 5.2 Subset-versus-full review-burden ranking

Store the validated subset rates in a versioned configuration table rather than reading a local-only output path. On the full HPC run:

- calculate full-data review proportions and ranks;
- report rank changes and absolute rate differences;
- calculate descriptive Spearman and Kendall rank agreement across four sections;
- plot subset and full rates together with section labels.

These four-section correlations are descriptive. Do not interpret their p-values as biological evidence. Write `subset_full_qc_ranking.tsv` and `subset_full_qc_rank_agreement.tsv`.

### 5.3 Spatial review summary

Aggregate the section spatial tables without pooling coordinates across sections. Report:

- edge/interior and dense/non-dense effects by section;
- global clustering statistics and permutation results;
- the number, size, effect magnitude, and bounding box of review hotspots;
- whether an optional morphology overlay was created or manual review is required.

Write `combined_spatial_qc.tsv`, `combined_spatial_hotspots.tsv`, a faceted spatial hotspot figure, and a concise review-priority table. Statistical clustering supports prioritization; only visual morphology review may classify a hotspot as an edge artifact, fold, tear, or dense aggregate.

### 5.4 Within-mouse technical concordance

Compare the verified pairs:

- Mouse 1: section 62308 (`Region_1`) versus 62309 (`Region_2`);
- Mouse 2: section 62310 (`Region_3`) versus 62311 (`Region_4`).

For each pair, report:

- cell yield and QC review-rate differences;
- ratios of median transcript counts, detected genes, cell area, and control fractions;
- one-dimensional quantile/Wasserstein-style distribution distances for cell-level QC variables, calculated from bounded quantile grids;
- Spearman correlation of per-gene log counts per 10,000;
- Spearman correlation of per-gene detection fractions;
- scatterplots with the identity line and labelled outlying genes.

Advisory defaults are configurable and are not universal Xenium acceptance standards: gene-profile Spearman correlation at least 0.90, absolute review-rate difference at most 5 percentage points, median count ratio between 0.67 and 1.50, and median detected-gene ratio between 0.80 and 1.25. The output reports each criterion separately and assigns `CONCORDANT`, `REVIEW`, or `NOT_ESTIMABLE`. It does not change Phase 0-2 readiness gates.

Write `within_mouse_section_concordance.tsv`, `within_mouse_gene_concordance.tsv`, and pair-specific concordance figures. With two mice, no mouse-level inferential model is attempted.

## 6. Notebook presentation

The section notebook gains clearly separated chunks for extended-input preflight, alarm evidence, gene transcript quality, spatial diagnostics, figures, and artifact reload checks. The summary notebook gains chunks matching the four scientific questions, followed by an updated readiness interpretation.

Every chunk states:

- inputs;
- outputs and exact paths;
- execution mode (`FULL_HPC`, `LOCAL_SUBSET`, or `NOT_RUN_LOCAL_SUBSET`);
- checks and failure conditions;
- scientific interpretation limits.

The notebooks remain copy/paste-friendly R notebooks. Reusable calculations stay in `R/source.R`; notebook cells orchestrate functions and display bounded results.

## 7. Execution modes and failure behavior

### `LOCAL_SUBSET`

- Run all pure-R summary and spatial functions on the deterministic subset.
- If transcript Parquet or morphology data is intentionally absent, write explicit skip-status tables and continue.
- Never imply that skipped local diagnostics passed.

### `FULL_HPC`

- Require all four complete Xenium section directories, verified metadata, `transcripts.parquet`, `arrow`, and `RANN`.
- Require projected/lazy transcript aggregation; reject code paths that collect the full transcript table.
- Fail the affected section with a clear message when required schemas or packages are absent.
- Preserve partial diagnostic evidence under the run directory but do not mark the extended diagnostics complete.

Optional morphology rendering may be unavailable. In that case the full run remains valid for coordinate-based spatial statistics and writes manual-review bounding boxes; it must not claim morphology inspection was completed.

## 8. Testing and verification

Use test-driven development for every new reusable behavior.

Pure-R fixtures will verify:

- direct alarm classification and unresolved-cycle wording;
- candidate tiers from hand-calculated gene contrasts;
- subset/full rank changes and descriptive rank agreement;
- edge and density enrichment counts/effect sizes;
- spatial clustering/hotspot detection on known clustered and randomized fixtures;
- within-mouse pair construction, metric contrasts, advisory flags, and `NOT_ESTIMABLE` behavior;
- path guards and deterministic permutation seeds;
- zero cell deletion and preservation of existing Phase 0-2 gates.

Integration checks will verify:

- current local subset execution produces valid skip records for unavailable large inputs;
- every section writes and reloads the complete extended artifact contract;
- the summary rejects missing or duplicate section diagnostics;
- output row counts and section/mouse identifiers reconcile;
- notebooks contain no local `C:` path and all temporary/output paths remain project-root constrained;
- the HPC launcher checks `arrow`, `RANN`, input files, memory-safe mode, and output roots before execution.

Because no local Parquet reader is installed, Arrow integration is an explicit HPC checkpoint. The repository will provide exact preflight and execution chunks; it will not install packages automatically.

## 9. Outputs and reporting

All new section outputs remain below:

`${RUN_ROOT}/sections/<Region_ID>/`

All combined outputs remain below:

`${RUN_ROOT}/slide_summary/`

The portable QC report and project progress record will be updated after local smoke tests and again after the full HPC results are available. The local report must distinguish tested calculations, skipped HPC-only diagnostics, candidate evidence, and unresolved 10x cycle identity.

## 10. Completion criteria

Implementation is complete only when:

1. all new pure-R tests have been observed failing before implementation and then pass;
2. both notebook templates validate and the local subset run completes with explicit skip states where necessary;
3. the extended output contract reloads without ambiguity for all four sections and the slide summary;
4. existing QC counts, panel reconciliation, metadata mapping, and readiness results remain unchanged unless a documented new full-data run provides new evidence;
5. the HPC runner contains exact full-data preflight and execution instructions;
6. the progress record and report identify the poor-cycle candidates as unconfirmed and exact cycle identity as requiring 10x diagnostics.
