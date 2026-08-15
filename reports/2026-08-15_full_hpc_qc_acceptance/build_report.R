options(stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: Rscript build_report.R <full_hpc_run_root> <repository_root>", call. = FALSE)
run_root <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
repo_root <- normalizePath(args[[2]], winslash = "/", mustWork = TRUE)
report_dir <- file.path(repo_root, "reports", "2026-08-15_full_hpc_qc_acceptance")
summary_root <- file.path(run_root, "slide_summary")
template_path <- file.path(repo_root, "reports", "2026-08-14_scwat_qc_summary", "artifact.json")
if (!requireNamespace("jsonlite", quietly = TRUE)) stop("jsonlite is required.", call. = FALSE)
read_tsv <- function(name) utils::read.delim(file.path(summary_root, name), check.names = FALSE)
as_rows <- function(data) lapply(seq_len(nrow(data)), function(index) as.list(data[index, , drop = FALSE]))
rate <- function(numerator, denominator) ifelse(denominator > 0, numerator / denominator, NA_real_)

artifact <- jsonlite::fromJSON(template_path, simplifyVector = FALSE)
qc <- read_tsv("combined_qc_summary.tsv")
readiness <- read_tsv("combined_readiness.tsv")
alarm <- read_tsv("combined_cycle_alarm_evidence.tsv")
candidates <- read_tsv("candidate_cycle_affected_genes.tsv")
ranking <- read_tsv("subset_full_qc_ranking.tsv")
rank_agreement <- read_tsv("subset_full_qc_rank_agreement.tsv")
spatial <- read_tsv("combined_spatial_qc.tsv")
hotspots <- read_tsv("combined_spatial_hotspots.tsv")
concordance <- read_tsv("within_mouse_section_concordance.tsv")
gene_qc <- read_tsv("combined_gene_transcript_quality.tsv")
cells <- utils::read.delim(gzfile(file.path(summary_root, "combined_cell_qc_metadata.tsv.gz")), check.names = FALSE)
gene_sets <- utils::read.delim(file.path(repo_root, "config", "eos_gene_sets.tsv"), check.names = FALSE)

if (!"section_candidate_flag" %in% names(candidates)) {
  candidates$section_candidate_flag <- candidates$abundance_depletion | candidates$quality_degradation
}
if (!"section_evidence_status" %in% names(candidates)) {
  candidates$section_evidence_status <- ifelse(
    candidates$abundance_depletion & candidates$quality_degradation, "DEPLETION_AND_Q20_LOSS",
    ifelse(candidates$abundance_depletion, "DEPLETION_ONLY", ifelse(candidates$quality_degradation, "Q20_LOSS_ONLY", "NO_SECTION_LEVEL_SIGNAL"))
  )
}

qc$core_pass_rate <- rate(qc$core_qc_pass, qc$input_cells)
qc$review_rate <- rate(qc$review_flagged, qc$input_cells)
qc$section_id <- c("62308", "62309", "62310", "62311")[match(qc$region_id, paste0("Region_", 1:4))]
qc$mouse_id <- c("Mouse_1", "Mouse_1", "Mouse_2", "Mouse_2")[match(qc$region_id, paste0("Region_", 1:4))]
readiness$xenium_alarm <- ifelse(readiness$status == "HOLD", "poor_quality_cycles_detected (ERROR)", "None reported")
readiness$biological_interpretation_allowed <- ifelse(readiness$status == "PASS", "Section only, with morphology caveat", "No")

affected <- candidates[candidates$abundance_depletion | candidates$quality_degradation, , drop = FALSE]
affected <- merge(affected, gene_sets, by = "gene", all.x = TRUE, sort = FALSE)
candidate_burden <- do.call(rbind, lapply(split(affected, affected$region_id), function(d) data.frame(
  region_id = d$region_id[[1]], affected_genes = length(unique(d$gene)),
  depletion_genes = length(unique(d$gene[d$abundance_depletion])),
  q20_loss_genes = length(unique(d$gene[d$quality_degradation])),
  both_signals = length(unique(d$gene[d$abundance_depletion & d$quality_degradation])),
  eos_common = length(unique(d$gene[d$gene_set == "common" & !is.na(d$gene_set)])),
  eos_short_lived = length(unique(d$gene[d$gene_set == "short_lived" & !is.na(d$gene_set)])),
  eos_long_lived = length(unique(d$gene[d$gene_set == "long_lived" & !is.na(d$gene_set)])),
  comparison_type = unique(d$comparison_type)[[1]], stringsAsFactors = FALSE
)))
rownames(candidate_burden) <- NULL
candidate_top <- affected[!is.na(affected$gene_set), , drop = FALSE]
candidate_top <- candidate_top[order(candidate_top$region_id, !(candidate_top$abundance_depletion & candidate_top$quality_degradation), candidate_top$log2_count_ratio, candidate_top$q20_difference), ]
candidate_top <- do.call(rbind, lapply(split(candidate_top, candidate_top$region_id), head, 8L))
candidate_top <- candidate_top[, c("region_id", "gene", "gene_set", "comparison_type", "log2_count_ratio", "q20_difference", "evidence_tier", "candidate_status")]

ranking$subset_review_fraction <- ranking$review_fraction
ranking$review_rate_change_pp <- 100 * ranking$review_fraction_difference
ranking <- ranking[, c("region_id", "subset_review_fraction", "subset_rank", "full_review_fraction", "full_rank", "rank_change", "review_rate_change_pp", "comparison_status")]
global_spatial <- spatial[spatial$diagnostic_type == "global_clustering", c("region_id", "statistic", "empirical_p", "permutations", "k")]
hotspot_counts <- aggregate(list(hotspot_bins = hotspots$hotspot_status == "MORPHOLOGY_REVIEW_REQUIRED"), list(region_id = hotspots$region_id), sum)
global_spatial <- merge(global_spatial, hotspot_counts, by = "region_id", all.x = TRUE, sort = FALSE)
global_spatial$interpretation <- ifelse(global_spatial$hotspot_bins > 0, "Global clustering plus candidate morphology-review bins", "Small but significant global clustering; no grid hotspot passed FDR")
spatial_effect <- spatial[spatial$diagnostic_type == "edge_density_enrichment" & ((spatial$class_type == "edge_proxy" & spatial$class == "edge") | (spatial$class_type == "local_density" & spatial$class == "dense")), c("region_id", "class_type", "cells", "flagged", "review_rate", "risk_ratio", "absolute_rate_difference")]
spatial_effect$contrast <- ifelse(spatial_effect$class_type == "edge_proxy", "Edge vs interior", "Dense vs non-dense")
spatial_effect <- spatial_effect[, c("region_id", "contrast", "cells", "flagged", "review_rate", "risk_ratio", "absolute_rate_difference")]

cell_distribution <- do.call(rbind, lapply(split(cells, cells$region_id), function(d) data.frame(
  region_id = d$region_id[[1]], median_counts = stats::median(d$nCount_Xenium),
  p05_counts = unname(stats::quantile(d$nCount_Xenium, 0.05)), p95_counts = unname(stats::quantile(d$nCount_Xenium, 0.95)),
  median_features = stats::median(d$nFeature_Xenium), p05_features = unname(stats::quantile(d$nFeature_Xenium, 0.05)),
  p95_features = unname(stats::quantile(d$nFeature_Xenium, 0.95)), median_cell_area = stats::median(d$cell_area),
  p95_cell_area = unname(stats::quantile(d$cell_area, 0.95)), stringsAsFactors = FALSE
)))
rownames(cell_distribution) <- NULL
gene_quality <- do.call(rbind, lapply(split(gene_qc, gene_qc$region_id), function(d) data.frame(
  region_id = d$region_id[[1]], panel_genes = nrow(d), zero_count_genes = sum(d$raw_counts == 0),
  median_fraction_q20 = stats::median(d$fraction_q20, na.rm = TRUE),
  genes_fraction_q20_below_0_8 = sum(d$fraction_q20 < 0.8, na.rm = TRUE),
  median_mean_qv = stats::median(d$mean_qv, na.rm = TRUE), stringsAsFactors = FALSE
)))
rownames(gene_quality) <- NULL
concordance$sections <- paste(concordance$section_a, concordance$section_b, sep = " / ")
concordance$qualified_status <- ifelse(concordance$median_area_ratio > 1.5, "PRESPECIFIED_CONCORDANT; AREA_REVIEW", concordance$concordance_status)
total_cells <- sum(qc$input_cells); total_core <- sum(qc$core_qc_pass); total_review <- sum(qc$review_flagged)
headline <- data.frame(total_cells = total_cells, sections = 4L, core_pass_cells = total_core, core_pass_rate = total_core / total_cells, review_flagged_cells = total_review, review_rate = total_review / total_cells, panel_genes_matched = 100L, alarm_positive_sections = 3L, cells_deleted = 0L)

source_root <- "adipose_analysis/hpc_return/2026-08-15_full_qc/scwat_qc_outputs/full_notebook_qc_v1/slide_summary"
for (index in seq_along(artifact$manifest$sources)) {
  item <- artifact$manifest$sources[[index]]
  item$path <- switch(item$id,
    "methods-source" = "R/source.R", "metadata-source" = "config/scwat_sample_manifest.tsv",
    "panel-source" = "adipose_analysis/hpc_return/2026-08-15_full_qc/scwat_qc_outputs/full_notebook_qc_v1/sections/Region_1/panel_reconciliation.tsv",
    "qc-summary-source" = paste0(source_root, "/combined_qc_summary.tsv"), "readiness-source" = paste0(source_root, "/combined_readiness.tsv"),
    "alarm-source" = paste0(source_root, "/combined_analysis_alerts.tsv"), "extended-alarm-source" = paste0(source_root, "/combined_cycle_alarm_evidence.tsv"),
    "candidate-source" = paste0(source_root, "/candidate_cycle_affected_genes.tsv"), "ranking-source" = paste0(source_root, "/subset_full_qc_ranking.tsv"),
    "spatial-source" = paste0(source_root, "/combined_spatial_qc.tsv"), "concordance-source" = paste0(source_root, "/within_mouse_section_concordance.tsv"), item$path)
  if (item$id == "qc-summary-source") {
    item$label <- "Full-data combined section QC summary"
    item$query$description <- "Section-specific summaries generated from all 286,178 segmented cells."
    item$query$filters <- list("scWAT", "Regions 1-4", "complete full-HPC sections")
  }
  if (item$id == "candidate-source") {
    item$label <- "Full-data provisional candidate affected-gene diagnostics"
    item$query$description <- "Cross-section count-depletion and transcript-QV diagnostics; candidates are not cycle-confirmed."
    item$query$filters <- list("alarm-positive Regions 1, 2, and 4", "prespecified depletion/Q20 thresholds")
  }
  if (item$id == "ranking-source") {
    item$label <- "Subset versus full-data QC ranking"
    item$query$description <- "Descriptive comparison of fixed subset and complete-section review burden."
  }
  if (item$id == "spatial-source") {
    item$label <- "Full-data spatial QC diagnostics"
    item$query$description <- "Complete-section kNN permutation, edge/density proxy, and grid-hotspot diagnostics."
    item$query$filters <- list("FULL_HPC", "k=15", "999 permutations", "100-micron bins")
  }
  if (item$id == "concordance-source") {
    item$label <- "Within-mouse full-data technical concordance"
    item$query$description <- "Advisory complete-section comparison for the two verified technical pairs."
  }
  if (item$id == "panel-source") {
    item$query$description <- "The 100 expected Eos genes reconciled with the installed custom Xenium panel."
    item$query$filters <- list("7 common", "47 short-lived", "46 long-lived")
  }
  artifact$manifest$sources[[index]] <- item
}
artifact$surface <- "report"
artifact$manifest$title <- "scWAT Xenium Full-Data QC Acceptance Report"
artifact$manifest$description <- "Technical acceptance assessment of four full scWAT Xenium sections returned from HPC on 2026-08-15."
artifact$manifest$generatedAt <- "2026-08-15T18:00:00+08:00"
artifact$snapshot$generatedAt <- artifact$manifest$generatedAt
artifact$snapshot$status <- "ready"
artifact$sources <- artifact$manifest$sources
artifact$manifest$cards[[1]]$description <- "All segmented cells assessed across the four complete sections."
artifact$manifest$cards[[2]]$description <- "Cells within section-specific robust transcript-count and detected-feature intervals."
artifact$manifest$cards[[3]]$description <- "Cells carrying at least one overlapping QC review flag."
artifact$manifest$cards[[4]]$description <- "Expected Eos target genes reconciled with every installed section panel."
artifact$manifest$cards[[5]]$description <- "The QC pipeline flags but does not automatically delete cells."
artifact$manifest$charts[[1]]$title <- "Full-data QC review-flag rate by section"
artifact$manifest$charts[[1]]$subtitle <- "Complete sections; denominator is all segmented cells in each section"
artifact$manifest$charts[[2]]$title <- "Global spatial clustering statistic by section"
artifact$manifest$charts[[2]]$subtitle <- "Complete sections; k=15 and 999 within-section label permutations"

artifact$manifest$tables[[1]]$subtitle <- "Complete sections; exact counts and rates"
artifact$manifest$tables[[2]]$subtitle <- "Worst-case integrity, panel, metadata, and Xenium-alarm gate"
artifact$manifest$tables[[4]]$title <- "Selected Eos candidate affected-gene comparisons"
artifact$manifest$tables[[4]]$subtitle <- "Up to eight strongest Eos-target rows per alarm-positive section; complete filtered list is saved separately"
artifact$manifest$tables[[4]]$columns <- list(
  list(field="region_id",label="Section",type="text"), list(field="gene",label="Gene",type="text"), list(field="gene_set",label="Eos set",type="text"),
  list(field="comparison_type",label="Comparison",type="text"), list(field="log2_count_ratio",label="log2 count ratio",format="decimal"),
  list(field="q20_difference",label="Q20 fraction difference",format="decimal"), list(field="evidence_tier",label="Cross-section tier",type="text"),
  list(field="candidate_status",label="Confirmation",type="text"))
artifact$manifest$tables[[5]]$title <- "Subset versus full-data review burden"
artifact$manifest$tables[[5]]$subtitle <- sprintf("Descriptive rank agreement across four sections: Spearman rho %.2f; Kendall tau %.2f", rank_agreement$spearman_rho[[1]], rank_agreement$kendall_tau[[1]])
artifact$manifest$tables[[5]]$columns <- list(
  list(field="region_id",label="Section",type="text"), list(field="subset_review_fraction",label="Subset review rate",format="percent"),
  list(field="subset_rank",label="Subset rank",format="number"), list(field="full_review_fraction",label="Full review rate",format="percent"),
  list(field="full_rank",label="Full rank",format="number"), list(field="rank_change",label="Rank change",format="number"),
  list(field="review_rate_change_pp",label="Rate change, pp",format="decimal"))
artifact$manifest$tables[[6]]$title <- "Full-data spatial clustering and hotspot screen"
artifact$manifest$tables[[6]]$subtitle <- "Global permutation test plus FDR-screened 100-micron candidate bins"
artifact$manifest$tables[[6]]$columns <- list(
  list(field="region_id",label="Section",type="text"), list(field="statistic",label="kNN statistic",format="decimal"),
  list(field="empirical_p",label="Empirical p",format="decimal"), list(field="hotspot_bins",label="Candidate bins",format="number"),
  list(field="interpretation",label="Interpretation",type="text"))
artifact$manifest$tables[[7]]$subtitle <- "Prespecified criteria plus a separate non-gated cell-area caution"
artifact$manifest$tables[[7]]$columns <- list(
  list(field="mouse_id",label="Mouse",type="text"), list(field="sections",label="Sections",type="text"),
  list(field="review_rate_difference",label="Review-rate difference",format="percent"), list(field="median_count_ratio",label="Median count ratio",format="decimal"),
  list(field="median_feature_ratio",label="Median feature ratio",format="decimal"), list(field="median_area_ratio",label="Median area ratio",format="decimal"),
  list(field="gene_count_spearman",label="Gene count rho",format="decimal"), list(field="gene_detection_spearman",label="Gene detection rho",format="decimal"),
  list(field="qualified_status",label="Qualified assessment",type="text"))

table_like <- function(id,title,subtitle,dataset,source_id,columns) list(id=id,title=title,subtitle=subtitle,dataset=dataset,sourceId=source_id,defaultSort=list(field="region_id",direction="asc"),density="spacious",layout="full",columns=columns)
artifact$manifest$tables[[8]] <- table_like("candidate-burden-table","Candidate affected-gene burden by section","Threshold-crossing genes; Regions 1/2 are cross-mouse diagnostics and Region 4 is within-mouse","candidate_burden","candidate-source",list(
  list(field="region_id",label="Section",type="text"),list(field="affected_genes",label="Affected",format="number"),list(field="depletion_genes",label="Depleted",format="number"),
  list(field="q20_loss_genes",label="Q20 loss",format="number"),list(field="both_signals",label="Both",format="number"),list(field="eos_common",label="Eos common",format="number"),
  list(field="eos_short_lived",label="Eos short-lived",format="number"),list(field="eos_long_lived",label="Eos long-lived",format="number"),list(field="comparison_type",label="Comparison",type="text")))
artifact$manifest$tables[[9]] <- table_like("gene-quality-table","Per-gene transcript-quality distribution","All 479 panel genes per section; Q20 values are descriptive and not cycle identities","gene_quality","candidate-source",list(
  list(field="region_id",label="Section",type="text"),list(field="panel_genes",label="Panel genes",format="number"),list(field="zero_count_genes",label="Zero-count genes",format="number"),
  list(field="median_fraction_q20",label="Median fraction Q20",format="percent"),list(field="genes_fraction_q20_below_0_8",label="Genes below 0.8",format="number"),list(field="median_mean_qv",label="Median mean QV",format="decimal")))
artifact$manifest$tables[[10]] <- table_like("cell-distribution-table","Cell-level QC distributions","Medians and 5th/95th percentiles from complete sections","cell_distribution","qc-summary-source",list(
  list(field="region_id",label="Section",type="text"),list(field="median_counts",label="Median counts",format="number"),list(field="p05_counts",label="Count P05",format="number"),
  list(field="p95_counts",label="Count P95",format="number"),list(field="median_features",label="Median features",format="number"),list(field="p05_features",label="Feature P05",format="number"),
  list(field="p95_features",label="Feature P95",format="number"),list(field="median_cell_area",label="Median area",format="decimal"),list(field="p95_cell_area",label="Area P95",format="decimal")))
artifact$manifest$tables[[11]] <- table_like("spatial-effect-table","Edge and dense-aggregate proxy effect sizes","Coordinate proxies do not establish morphology","spatial_effect","spatial-source",list(
  list(field="region_id",label="Section",type="text"),list(field="contrast",label="Contrast",type="text"),list(field="cells",label="Proxy cells",format="number"),
  list(field="flagged",label="Flagged",format="number"),list(field="review_rate",label="Proxy rate",format="percent"),list(field="risk_ratio",label="Risk ratio",format="decimal"),
  list(field="absolute_rate_difference",label="Absolute difference",format="percent")))

artifact$manifest$blocks <- list(
  list(id="title",type="markdown",body="# scWAT Xenium Full-Data QC Acceptance Report",layout="full"),
  list(id="technical-summary",type="markdown",sourceId="qc-summary-source",layout="full",body=sprintf("## Cell-level QC passes, but the slide is not yet acceptable for unrestricted gene-level biology\n\nAcross %s cells, %.2f%% pass the section-specific core count/feature gate and %.2f%% carry at least one broader review flag; no cells were deleted. This supports continued technical processing and morphology review. It does **not** override the separate Xenium alarm gate described below.",format(total_cells,big.mark=","),100*total_core/total_cells,100*total_review/total_cells)),
  list(id="acceptance-summary",type="markdown",sourceId="readiness-source",layout="full",body="### Acceptance decision\n\n**Overall: HOLD for comparative biological analysis.** Region 3 is conditionally acceptable for section-level exploratory work with the morphology caveat documented below. Regions 1, 2, and 4 remain HOLD until 10x diagnostics identify the affected cycles/codewords and the resulting gene scope is documented. Raw data and all cells should be retained."),
  list(id="headline-metrics",type="metric-strip",cardIds=lapply(artifact$manifest$cards,`[[`,"id"),layout="full"),
  list(id="cell-heading",type="markdown",sourceId="qc-summary-source",layout="full",body="## All four sections have strong cell-level QC\n\nCore-pass rates are 99.06% to 99.62%, while review-flag rates are 3.22% to 5.62%. Region 3 has the highest full-data review burden and Region 4 the lowest. These rates quantify segmentation/count review burden; they do not override the independent imaging-cycle alarm gate."),
  list(id="review-chart-block",type="chart",chartId="review-rate-chart",layout="full"),list(id="section-table-block",type="table",tableId="section-qc-table",layout="full"),list(id="cell-distribution-block",type="table",tableId="cell-distribution-table",layout="full"),
  list(id="alarm-heading",type="markdown",sourceId="extended-alarm-source",layout="full",body="## Direct 10x evidence identifies three alarm-positive sections, not the exact cycles\n\nRegions 1, 2, and 4 each report `poor_quality_cycles_detected`; Region 3 reports no poor-cycle alarm. The supplied Xenium summaries do not identify cycle number, fluorophore, codeword, or a confirmed gene list. Exact cycle identity therefore requires the corresponding 10x diagnostic package/support output."),list(id="alarm-evidence-table-block",type="table",tableId="alarm-evidence-table",layout="full"),
  list(id="candidate-heading",type="markdown",sourceId="candidate-source",layout="full",body="## Region 4 shows broad within-mouse depletion and Q20 loss; Regions 1 and 2 remain cross-mouse diagnostics\n\nAgainst its paired section Region 3, Region 4 has 353/479 genes crossing at least one candidate threshold, including 200 with both depletion and Q20 loss. This broad pattern is consistent with a substantial section-level technical effect, but it is not a cycle assignment. Regions 1 and 2 have 198 and 174 candidates respectively; because their reference is from the other mouse, those counts can mix technical and inter-mouse biological variation. Every row remains `CANDIDATE_NOT_CONFIRMED`."),list(id="candidate-burden-block",type="table",tableId="candidate-burden-table",layout="full"),list(id="candidate-table-block",type="table",tableId="candidate-table",layout="full"),list(id="gene-quality-block",type="table",tableId="gene-quality-table",layout="full"),
  list(id="ranking-heading",type="markdown",sourceId="ranking-source",layout="full",body=sprintf("## Full data broadly reproduces, but does not exactly preserve, the subset ranking\n\nThe full-data order is Region 3, Region 1, Region 2, Region 4 versus Region 1, Region 3, Region 2, Region 4 in the subset. Spearman rho is %.2f and Kendall tau is %.2f across only four sections. The subset correctly identified the two lower-burden sections, but Regions 1 and 3 swapped; full-data values should be used for decisions.",rank_agreement$spearman_rho[[1]],rank_agreement$kendall_tau[[1]])),list(id="ranking-table-block",type="table",tableId="ranking-status-table",layout="full"),
  list(id="spatial-heading",type="markdown",sourceId="spatial-source",layout="full",body="## Review flags are weakly globally clustered, with localized candidate bins only in Region 3\n\nAll four global kNN tests have empirical p=0.001, but the statistics are small (0.0197 to 0.0496) and the very large cell counts make statistical significance easy to obtain. Edge proxies show modest enrichment in Regions 1, 2, and 4 (risk ratios 1.49, 1.35, and 1.81) but essentially none in Region 3 (1.03). Dense-cell proxies are not enriched in any section. Only Region 3 contains FDR-positive grid bins: 34 bins with median review rate 29.5% versus 5.62% globally. Coordinates alone cannot label edges, folds, tears, or aggregates; these bins require image review."),list(id="spatial-chart-block",type="chart",chartId="spatial-clustering-chart",layout="full"),list(id="spatial-table-block",type="table",tableId="spatial-table",layout="full"),list(id="spatial-effect-block",type="table",tableId="spatial-effect-table",layout="full"),
  list(id="concordance-heading",type="markdown",sourceId="concordance-source",layout="full",body="## Both mouse pairs meet the prespecified concordance criteria, but Mouse 2 has a large cell-area shift\n\nMouse 1 (62308/62309) and Mouse 2 (62310/62311) meet all five prespecified advisory criteria for review rate, median counts, median features, gene-count correlation, and gene-detection correlation. Mouse 2 nevertheless has a 1.82-fold median cell-area ratio and area-distribution distance 0.665; cell area was reported but not part of the original gate. Treat Mouse 2 as concordant for the prespecified metrics with a separate segmentation/morphology review requirement."),list(id="concordance-table-block",type="table",tableId="concordance-table",layout="full"),list(id="readiness-table-block",type="table",tableId="readiness-table",layout="full"),
  list(id="scope-definitions",type="markdown",layout="full",body="## Scope, data, and metric definitions\n\n**Population:** all 286,178 segmented cells in four WT, untreated, normal 8-week scWAT sections. **Biological replicate:** mouse (n=2). **Technical unit:** section; 62308/62309 are Mouse 1 and 62310/62311 are Mouse 2. **Core pass:** count and detected-feature values within robust section-specific bounds. **Review flagged:** one or more overlapping core, nucleus, cell-area, high-control, or segmentation-multiplet flags. **Candidate gene:** a section/reference comparison crossing log2 depletion <= -0.5 and/or Q20-fraction difference <= -0.05. Candidate status is diagnostic, not confirmed."),
  list(id="methodology",type="markdown",sourceId="methods-source",layout="full",body="## Methodology and robustness\n\nEach section was processed independently. Full transcript tables were projected to gene/QV/codeword fields and summarized per gene; cell matrices remained sparse. Spatial clustering used k=15 and 999 deterministic within-section label permutations; hotspots used 100-micron bins, at least 20 cells, and BH FDR <=0.05. Subset/full agreement is descriptive because there are only four sections. No cell was automatically removed, and no cell is treated as a biological replicate."),
  list(id="limitations",type="markdown",layout="full",body="## Limitations and uncertainty\n\nThe poor-cycle alarm is direct evidence, but cycle identity and affected-gene mapping are absent. Region 3 is used as the diagnostic reference; only Region 4 is a within-mouse comparison, so Region 1/2 gene candidates are confounded by mouse. Spatial proxies are coordinate-based and cannot diagnose morphology without DAPI/morphology images. Two mice are insufficient for population-level biological inference. The verified 100-gene Eos list contains 7 common, 47 short-lived, and 46 long-lived genes; this resolves the earlier arithmetic inconsistency with the stated six common genes."),
  list(id="next-steps",type="markdown",layout="full",body="## Recommended next steps\n\n1. Request/export the 10x cycle/codeword diagnostics for Regions 1, 2, and 4; keep their biological gate at HOLD until reviewed.\n2. Inspect Region 3's 34 candidate bins and the edge-enriched regions against morphology/DAPI images in Xenium Explorer.\n3. Review the Mouse 2 cell-area shift, especially segmentation boundaries and tissue composition in Regions 3/4.\n4. If 10x identifies affected genes, exclude or sensitivity-test those genes before normalization, annotation, Eos-state scoring, or spatial differential analysis.\n5. Use mouse as the replicate for downstream inference; with n=2 and no treatment contrast, keep analyses descriptive/exploratory unless more animals are added."),
  list(id="further-questions",type="markdown",layout="full",body="## Further questions\n\n- Which exact cycles, channels, and codewords were impaired according to 10x diagnostics?\n- Do Region 3 hotspots correspond to tissue edge, fold, tear, segmentation crowding, or true anatomical structure?\n- Does the Region 3/4 cell-area shift persist after morphology-guided segmentation review?\n- After confirmed affected genes are handled, do within-mouse expression profiles remain concordant?")
)
artifact$snapshot$datasets <- list(headline_metrics=as_rows(headline),qc_by_section=as_rows(qc),readiness_by_section=as_rows(readiness),alarm_evidence=as_rows(alarm),candidate_top=as_rows(candidate_top),candidate_burden=as_rows(candidate_burden),ranking_status=as_rows(ranking),spatial_global=as_rows(global_spatial),spatial_effect=as_rows(spatial_effect),within_mouse_concordance=as_rows(concordance),gene_quality=as_rows(gene_quality),cell_distribution=as_rows(cell_distribution))
utils::write.table(
  affected,
  file.path(report_dir, "candidate_cycle_affected_genes_affected_only.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
artifact_path <- file.path(report_dir, "artifact.json")
jsonlite::write_json(artifact, artifact_path, auto_unbox=TRUE, pretty=TRUE, null="null", digits=NA)
cat(artifact_path, "\n")
