options(stringsAsFactors = FALSE)

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (!length(script_arg)) stop("Run this test with Rscript.", call. = FALSE)
test_path <- normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = TRUE)
repo_root <- dirname(dirname(test_path))
source(file.path(repo_root, "R", "source.R"))

expect_error <- function(expr, pattern = NULL) {
  error <- tryCatch({ force(expr); NULL }, error = identity)
  stopifnot(inherits(error, "error"))
  if (!is.null(pattern)) stopifnot(grepl(pattern, conditionMessage(error), ignore.case = TRUE))
  invisible(error)
}

# A broken config parser would silently change scientific thresholds.
config_path <- file.path(repo_root, "config", "extended_qc_defaults.tsv")
config <- read_extended_qc_config(config_path)
stopifnot(config$qv_threshold == 20, config$spatial_k == 15L, config$permutations == 999L)
stopifnot(config$seed == 20260814L, is.integer(config$seed), is.numeric(config$dense_quantile))

bad_config <- utils::read.delim(config_path, check.names = FALSE)
bad_config <- rbind(bad_config, bad_config[1, , drop = FALSE])
bad_path <- file.path(tempdir(), "duplicate_extended_qc_config.tsv")
utils::write.table(bad_config, bad_path, sep = "\t", quote = FALSE, row.names = FALSE)
expect_error(read_extended_qc_config(bad_path), "duplicate")

# AUTO must select full mode only when the large transcript input exists.
mode_root <- file.path(tempdir(), "extended_qc_modes")
unlink(mode_root, recursive = TRUE, force = TRUE)
subset_region_dir <- file.path(mode_root, "subset")
full_region_dir <- file.path(mode_root, "full")
dir.create(subset_region_dir, recursive = TRUE)
dir.create(full_region_dir, recursive = TRUE)
invisible(file.create(file.path(full_region_dir, "transcripts.parquet")))
stopifnot(resolve_extended_qc_mode("AUTO", subset_region_dir) == "LOCAL_SUBSET")
stopifnot(resolve_extended_qc_mode("AUTO", full_region_dir) == "FULL_HPC")
stopifnot(resolve_extended_qc_mode("local_subset", full_region_dir) == "LOCAL_SUBSET")
expect_error(resolve_extended_qc_mode("BAD", subset_region_dir), "AUTO")

# Local mode may skip large inputs but must report the skip explicitly.
preflight <- extended_qc_preflight("LOCAL_SUBSET", subset_region_dir, config)
stopifnot(preflight$status[preflight$check == "transcripts_parquet"] == "SKIP_ALLOWED")
stopifnot(preflight$status[preflight$check == "arrow"] == "SKIP_ALLOWED")
full_preflight <- extended_qc_preflight("FULL_HPC", full_region_dir, config)
if (!requireNamespace("arrow", quietly = TRUE)) stopifnot(full_preflight$status[full_preflight$check == "arrow"] == "FAIL")
if (!requireNamespace("RANN", quietly = TRUE)) stopifnot(full_preflight$status[full_preflight$check == "RANN"] == "FAIL")
stopifnot("dplyr" %in% full_preflight$check)
if (!requireNamespace("dplyr", quietly = TRUE)) stopifnot(full_preflight$status[full_preflight$check == "dplyr"] == "FAIL")

# Alarm evidence must never imply that the exact cycle or gene mapping is known.
alarm_fixture <- data.frame(
  raw_value = TRUE, formatted_value = "True", raised = TRUE,
  title = "Poor quality imaging cycles detected", message = "Diagnostic message",
  level = "ERROR", id = "poor_quality_cycles_detected", stringsAsFactors = FALSE
)
evidence <- build_cycle_alarm_evidence(alarm_fixture, "Region_1")
stopifnot(evidence$evidence_status == "DIRECT_EVIDENCE")
stopifnot(evidence$cycle_identity_status == "CYCLE_IDENTITY_UNRESOLVED_REQUIRES_10X")
stopifnot(evidence$gene_effect_status == "GENE_EFFECT_UNCONFIRMED")
no_evidence <- build_cycle_alarm_evidence(alarm_fixture[0, ], "Region_3")
stopifnot(no_evidence$evidence_status == "NO_ALARM_REPORTED", nrow(no_evidence) == 1L)

# Transcript summaries use literal hand-calculated QV and codeword results.
transcript_fixture <- data.frame(
  feature_name = c("A", "A", "B"), qv = c(30, 10, 25), codeword_index = c(1L, 1L, 2L),
  stringsAsFactors = FALSE
)
schema <- resolve_transcript_schema(names(transcript_fixture))
stopifnot(schema$gene == "feature_name", schema$qv == "qv", schema$codeword == "codeword_index")
expect_error(resolve_transcript_schema(c("gene", "feature_name", "qv")), "ambiguous")
expect_error(resolve_transcript_schema(c("feature_name", "x")), "QV")
transcript_qc <- summarise_transcript_quality_table(transcript_fixture, "Region_1", 20)
stopifnot(transcript_qc$transcript_rows[transcript_qc$gene == "A"] == 2L)
stopifnot(transcript_qc$mean_qv[transcript_qc$gene == "A"] == 20)
stopifnot(transcript_qc$fraction_q20[transcript_qc$gene == "A"] == 0.5)
stopifnot(transcript_qc$represented_codewords[transcript_qc$gene == "A"] == 1L)

# Arrow-compatible queries aggregate before collection and avoid unsupported n_distinct().
require_package("dplyr")
projected_fixture <- data.frame(
  gene = c("A", "A", "B", "B"), qv = c(30, 10, 25, NA),
  codeword = c(1L, 1L, 2L, 3L), stringsAsFactors = FALSE
)
transcript_queries <- build_transcript_quality_queries(projected_fixture, qv_threshold = 20, has_codeword = TRUE)
bounded_fixture <- dplyr::collect(transcript_queries$summary)
codeword_fixture <- dplyr::collect(transcript_queries$codewords)
stopifnot(nrow(bounded_fixture) == 2L, nrow(codeword_fixture) == 2L)
stopifnot(bounded_fixture$transcript_rows[bounded_fixture$gene == "A"] == 2L)
stopifnot(bounded_fixture$fraction_q20[bounded_fixture$gene == "A"] == 0.5)
stopifnot(!anyNA(bounded_fixture$mean_qv), !anyNA(bounded_fixture$fraction_q20))

# Matrix summaries preserve all panel genes, including zero-count genes.
require_package("Matrix")
count_fixture <- Matrix::Matrix(matrix(c(2, 0, 3, 0, 0, 0), nrow = 3L, byrow = TRUE), sparse = TRUE)
rownames(count_fixture) <- c("A", "B", "C")
colnames(count_fixture) <- c("cell1", "cell2")
gene_sets <- c(A = "common", B = "short_lived", C = "long_lived")
matrix_qc <- summarise_gene_matrix_qc(count_fixture, "Region_1", gene_sets)
stopifnot(matrix_qc$raw_counts[matrix_qc$gene == "A"] == 2)
stopifnot(matrix_qc$detection_fraction[matrix_qc$gene == "B"] == 0.5)
stopifnot(matrix_qc$raw_counts[matrix_qc$gene == "C"] == 0)
combined_qc <- combine_gene_quality(matrix_qc, transcript_qc)
stopifnot(nrow(combined_qc) == 3L, combined_qc$transcript_status[combined_qc$gene == "C"] == "NO_TRANSCRIPTS")

# Spatial edge logic must identify the 20 perimeter cells of a complete 6x6 grid.
spatial_fixture <- expand.grid(x_centroid = seq(5, 55, 10), y_centroid = seq(5, 55, 10))
spatial_fixture$cell_id <- sprintf("cell_%02d", seq_len(nrow(spatial_fixture)))
spatial_fixture$region_id <- "Region_1"
spatial_fixture$qc_review_flag <- spatial_fixture$x_centroid < 20 & spatial_fixture$y_centroid < 20
spatial_fixture <- spatial_fixture[, c("region_id", "cell_id", "x_centroid", "y_centroid", "qc_review_flag")]
spatial_grid <- assign_spatial_grid(spatial_fixture, grid_size_um = 10)
stopifnot(sum(spatial_grid$edge_proxy) == 20L)
edge_summary <- summarise_spatial_enrichment(spatial_grid)
stopifnot(edge_summary$flagged[edge_summary$class == "edge"] == 3L)
stopifnot(edge_summary$flagged[edge_summary$class == "interior"] == 1L)
stopifnot(all(c("risk_ratio", "absolute_rate_difference") %in% names(edge_summary)))
stopifnot(length(unique(edge_summary$risk_ratio)) == 1L, is.finite(unique(edge_summary$risk_ratio)))
stopifnot(length(unique(edge_summary$absolute_rate_difference)) == 1L, is.finite(unique(edge_summary$absolute_rate_difference)))
stopifnot(abs(unique(edge_summary$risk_ratio) - ((3.5 / 21) / (1.5 / 17))) < 1e-12)
stopifnot(abs(unique(edge_summary$absolute_rate_difference) - (3 / 20 - 1 / 16)) < 1e-12)
all_edge_fixture <- spatial_grid
all_edge_fixture$edge_proxy <- TRUE
all_edge_summary <- summarise_spatial_enrichment(all_edge_fixture)
stopifnot(all(is.na(all_edge_summary$risk_ratio)), all(is.na(all_edge_summary$absolute_rate_difference)))

# A known four-cell corner cluster should produce a morphology-review hotspot.
hotspot_grid <- assign_spatial_grid(spatial_fixture, grid_size_um = 20)
hotspots_1 <- find_spatial_qc_hotspots(hotspot_grid, permutations = 199L, min_bin_cells = 4L, fdr = 0.10, seed = 20260814L)
hotspots_2 <- find_spatial_qc_hotspots(hotspot_grid, permutations = 199L, min_bin_cells = 4L, fdr = 0.10, seed = 20260814L)
stopifnot(identical(hotspots_1, hotspots_2))
stopifnot(any(hotspots_1$hotspot_status == "MORPHOLOGY_REVIEW_REQUIRED"))
stopifnot(all(c("fdr_threshold", "min_bin_cells_threshold") %in% names(hotspots_1)))
stopifnot(all(hotspots_1$fdr_threshold == 0.10), all(hotspots_1$min_bin_cells_threshold == 4L))

# kNN calculations are deterministic locally and return explicit non-estimable states.
density <- calculate_knn_density(spatial_fixture, k = 4L, mode = "LOCAL_SUBSET")
stopifnot(length(density) == nrow(spatial_fixture), all(is.finite(density)), all(density > 0))
clustered <- test_spatial_flag_clustering(spatial_fixture, k = 4L, permutations = 199L, seed = 20260814L, mode = "LOCAL_SUBSET")
stopifnot(clustered$status == "ESTIMATED", clustered$permutations == 199L)
no_flags <- spatial_fixture; no_flags$qc_review_flag <- FALSE
not_estimable <- test_spatial_flag_clustering(no_flags, k = 4L, permutations = 19L, seed = 20260814L, mode = "LOCAL_SUBSET")
stopifnot(not_estimable$status == "NOT_ESTIMABLE")

# Candidate tiers distinguish paired, recurring, partial, and insufficient evidence.
candidate_fixture <- expand.grid(region_id = paste0("Region_", 1:4), gene = c("A", "B", "C", "D"), stringsAsFactors = FALSE)
candidate_fixture$counts_per_10000 <- 100
candidate_fixture$detection_fraction <- 0.5
candidate_fixture$fraction_q20 <- 0.90
candidate_fixture$counts_per_10000[candidate_fixture$gene == "A" & candidate_fixture$region_id == "Region_4"] <- 40
candidate_fixture$fraction_q20[candidate_fixture$gene == "A" & candidate_fixture$region_id == "Region_4"] <- 0.80
candidate_fixture$counts_per_10000[candidate_fixture$gene == "B" & candidate_fixture$region_id %in% c("Region_1", "Region_2")] <- 40
candidate_fixture$fraction_q20[candidate_fixture$gene == "B" & candidate_fixture$region_id %in% c("Region_1", "Region_2")] <- 0.80
candidate_fixture$counts_per_10000[candidate_fixture$gene == "C" & candidate_fixture$region_id == "Region_1"] <- 40
candidate_fixture$fraction_q20[candidate_fixture$gene == "D"] <- NA_real_
candidate_rank <- rank_candidate_cycle_genes(candidate_fixture, config)
tier <- setNames(candidate_rank$evidence_tier[!duplicated(candidate_rank$gene)], candidate_rank$gene[!duplicated(candidate_rank$gene)])
stopifnot(tier[["A"]] == "Tier_A", tier[["B"]] == "Tier_B", tier[["C"]] == "Tier_C", tier[["D"]] == "Unranked")
stopifnot(all(candidate_rank$candidate_status == "CANDIDATE_NOT_CONFIRMED"))
stopifnot(all(candidate_rank$exact_cycle_status == "REQUIRES_10X_DIAGNOSTICS"))
stopifnot(candidate_rank$comparison_type[candidate_rank$gene == "A" & candidate_rank$region_id == "Region_4"] == "WITHIN_MOUSE_TECHNICAL_PAIR")
stopifnot(candidate_rank$section_candidate_flag[candidate_rank$gene == "A" & candidate_rank$region_id == "Region_4"])
stopifnot(!candidate_rank$section_candidate_flag[candidate_rank$gene == "A" & candidate_rank$region_id == "Region_1"])
stopifnot(candidate_rank$section_evidence_status[candidate_rank$gene == "A" & candidate_rank$region_id == "Region_4"] == "DEPLETION_AND_Q20_LOSS")
stopifnot(candidate_rank$section_evidence_status[candidate_rank$gene == "C" & candidate_rank$region_id == "Region_1"] == "DEPLETION_ONLY")

# Evidence-only gene decisions are overlapping status fields: the conservative
# zero-alarm set is a subset of provisional primary features.
gene_decision_fixture <- build_gene_downstream_decision(
  candidate_rank, panel_genes = c("A", "B", "C", "D"),
  run_label = "unit_run", execution_mode = "LOCAL_SUBSET", provenance = "unit_fixture"
)
decision_by_gene <- gene_decision_fixture[match(c("A", "B", "C", "D"), gene_decision_fixture$gene), ]
stopifnot(identical(decision_by_gene$alarm_positive_section_count, c(1L, 2L, 1L, 0L)))
stopifnot(identical(decision_by_gene$primary_feature_status,
                    c("PROVISIONAL_PRIMARY_FEATURES", "EXCLUDED_FROM_PRIMARY_FEATURES",
                      "PROVISIONAL_PRIMARY_FEATURES", "PROVISIONAL_PRIMARY_FEATURES")))
stopifnot(decision_by_gene$conservative_evidence_status[[4]] == "CONSERVATIVE_NO_SIGNAL_DETECTED")
stopifnot(decision_by_gene$technical_risk_status[[2]] == "TECHNICAL_RISK_SENSITIVITY_ONLY")
stopifnot(all(decision_by_gene$raw_panel_status == "RAW_COMPLETE_PANEL"))
stopifnot(!any(grepl("CONFIRMED_(AFFECTED|UNAFFECTED)", unlist(decision_by_gene))))

eos_fixture <- data.frame(gene = c("A", "B", "D", "X"), gene_set = c("common", "common", "short_lived", "long_lived"))
eos_decision_fixture <- build_eos_gene_decision(
  gene_decision_fixture, eos_fixture, "unit_run", "LOCAL_SUBSET", "unit_fixture"
)
stopifnot(identical(eos_decision_fixture$retained_provisional, c(TRUE, FALSE, TRUE, FALSE)))
stopifnot(all(eos_decision_fixture$complete_signature_status == "RAW_COMPLETE_EOS_100"))

# Full-data ranks are compared descriptively with the fixed subset reference.
subset_reference_fixture <- utils::read.delim(file.path(repo_root, "config", "subset_qc_reference.tsv"), check.names = FALSE)
full_fixture <- data.frame(region_id = paste0("Region_", 1:4), input_cells = 1000L, review_flagged = c(80L, 20L, 60L, 40L))
rank_comparison <- compare_subset_full_qc(full_fixture, subset_reference_fixture)
stopifnot(identical(rank_comparison$ranking$full_rank, c(1L, 4L, 2L, 3L)))
stopifnot(is.numeric(rank_comparison$agreement$spearman_rho), is.numeric(rank_comparison$agreement$kendall_tau))
stopifnot(rank_comparison$agreement$interpretation == "DESCRIPTIVE_FOUR_SECTIONS")

# Verified metadata must produce exactly one technical pair per mouse.
manifest_fixture <- utils::read.delim(file.path(repo_root, "config", "scwat_sample_manifest.tsv"), check.names = FALSE)
pairs <- build_section_pairs(manifest_fixture)
stopifnot(identical(pairs$section_a, c("62308", "62310")))
stopifnot(identical(pairs$section_b, c("62309", "62311")))
stopifnot(identical(pairs$region_a, c("Region_1", "Region_3")))
stopifnot(identical(pairs$region_b, c("Region_2", "Region_4")))
incomplete_manifest <- manifest_fixture[manifest_fixture$region_id != "Region_2", , drop = FALSE]
incomplete_pairs <- build_section_pairs(incomplete_manifest)
stopifnot(incomplete_pairs$pair_status[incomplete_pairs$mouse_id == "Mouse_1"] == "NOT_ESTIMABLE")

# Mouse 1 is concordant; Mouse 2 deliberately breaches review/count/gene criteria.
concordance_cells <- do.call(rbind, list(
  data.frame(region_id="Region_1", nCount_Xenium=c(90,100,110), nFeature_Xenium=c(45,50,55), cell_area=c(190,200,210), control_fraction=c(.01,.01,.02)),
  data.frame(region_id="Region_2", nCount_Xenium=c(95,100,105), nFeature_Xenium=c(48,50,52), cell_area=c(195,200,205), control_fraction=c(.01,.01,.02)),
  data.frame(region_id="Region_3", nCount_Xenium=c(90,100,110), nFeature_Xenium=c(45,50,55), cell_area=c(190,200,210), control_fraction=c(.01,.01,.02)),
  data.frame(region_id="Region_4", nCount_Xenium=c(180,200,220), nFeature_Xenium=c(72,80,88), cell_area=c(380,400,420), control_fraction=c(.03,.04,.05))
))
concordance_summary <- data.frame(
  region_id=paste0("Region_",1:4), input_cells=1000L, review_flagged=c(20L,30L,10L,100L), stringsAsFactors=FALSE
)
concordance_genes <- do.call(rbind, lapply(paste0("Region_",1:4), function(region) {
  values <- if (region == "Region_4") 20:1 else 1:20
  data.frame(region_id=region, gene=sprintf("Gene%02d",1:20), counts_per_10000=values,
             detection_fraction=values/25, stringsAsFactors=FALSE)
}))
concordance <- calculate_within_mouse_concordance(
  manifest_fixture, concordance_summary, concordance_cells, concordance_genes, config
)
stopifnot(concordance$summary$concordance_status[concordance$summary$mouse_id == "Mouse_1"] == "CONCORDANT")
stopifnot(concordance$summary$concordance_status[concordance$summary$mouse_id == "Mouse_2"] == "REVIEW")
stopifnot(abs(concordance$summary$review_rate_difference[concordance$summary$mouse_id == "Mouse_1"] - 0.01) < 1e-12)
stopifnot(abs(concordance$summary$median_count_ratio[concordance$summary$mouse_id == "Mouse_2"] - 2) < 1e-12)
stopifnot(concordance$summary$gene_count_spearman[concordance$summary$mouse_id == "Mouse_1"] == 1)
stopifnot(concordance$summary$gene_count_spearman[concordance$summary$mouse_id == "Mouse_2"] == -1)
stopifnot(nrow(concordance$genes) == 40L)
incomplete_concordance <- calculate_within_mouse_concordance(
  incomplete_manifest, concordance_summary, concordance_cells, concordance_genes, config
)
stopifnot(incomplete_concordance$summary$concordance_status[incomplete_concordance$summary$mouse_id == "Mouse_1"] == "NOT_ESTIMABLE")

# Extended section artifacts must preserve matrix metrics and record local transcript skips.
extended_root <- file.path(tempdir(), "extended_section_artifacts")
unlink(extended_root, recursive = TRUE, force = TRUE)
dir.create(extended_root, recursive = TRUE)
local_gene_quality <- matrix_qc
local_gene_quality$transcript_rows <- NA_integer_
local_gene_quality$mean_qv <- NA_real_
local_gene_quality$fraction_q20 <- NA_real_
local_gene_quality$represented_codewords <- NA_integer_
local_gene_quality$transcript_status <- "NOT_RUN_LOCAL_SUBSET"
spatial_annotations <- spatial_grid
spatial_annotations$local_density <- calculate_knn_density(spatial_annotations, k = 4L, mode = "LOCAL_SUBSET")
spatial_annotations$dense_aggregate <- spatial_annotations$local_density >= stats::quantile(spatial_annotations$local_density, 0.90)
spatial_global <- test_spatial_flag_clustering(spatial_annotations, k = 4L, permutations = 19L, seed = 20260814L, mode = "LOCAL_SUBSET")
spatial_enrichment <- summarise_spatial_enrichment(spatial_annotations)
manual_review <- hotspots_1[hotspots_1$hotspot_status == "MORPHOLOGY_REVIEW_REQUIRED", , drop = FALSE]
extended_plots <- plot_extended_spatial_qc(spatial_annotations, spatial_enrichment, hotspots_1, "Region_1")
stopifnot(all(c("review_map", "edge_density", "hotspots") %in% names(extended_plots)))
stopifnot(all(vapply(extended_plots, inherits, logical(1), what = "ggplot")))

extended_paths <- write_extended_section_artifacts(
  project_root = extended_root, output_dir = file.path(extended_root, "Region_1"), region_id = "Region_1",
  mode = "LOCAL_SUBSET", preflight = preflight, cycle_alarm_evidence = evidence,
  gene_quality = local_gene_quality, spatial_global = spatial_global,
  spatial_edge_density = spatial_enrichment, spatial_hotspots = hotspots_1,
  spatial_cells = spatial_annotations, manual_review_manifest = manual_review,
  plots = extended_plots
)
required_extended <- extended_section_required_artifacts("Region_1", "LOCAL_SUBSET")
stopifnot(all(file.exists(file.path(extended_root, "Region_1", required_extended))))
stopifnot(all(file.exists(extended_paths)))
stopifnot(validate_extended_section_artifacts(file.path(extended_root, "Region_1"), "Region_1", "LOCAL_SUBSET"))
reloaded_extended <- read_extended_section_artifacts(file.path(extended_root, "Region_1"), "Region_1", "LOCAL_SUBSET")
stopifnot(all(reloaded_extended$gene_quality$transcript_status == "NOT_RUN_LOCAL_SUBSET"))
stopifnot(all(reloaded_extended$gene_quality$raw_counts == local_gene_quality$raw_counts))
expect_error(
  validate_extended_section_artifacts(file.path(extended_root, "Region_1"), "Region_1", "FULL_HPC", stop_on_error = TRUE),
  "NOT_RUN_LOCAL_SUBSET"
)

# Exactly four same-mode extended section bundles are required for slide aggregation.
extended_slide_root <- file.path(tempdir(), "extended_slide_run")
unlink(extended_slide_root, recursive = TRUE, force = TRUE)
for (region_index in seq_len(4L)) {
  region <- paste0("Region_", region_index)
  region_cells <- concordance_cells[concordance_cells$region_id == region, , drop = FALSE]
  region_cells$cell_id <- paste0(region, "_cell_", seq_len(nrow(region_cells)))
  region_cells$x_centroid <- seq_len(nrow(region_cells)) * 10
  region_cells$y_centroid <- seq_len(nrow(region_cells)) * 5 + region_index
  region_cells$qc_review_flag <- seq_len(nrow(region_cells)) == 1L
  region_cells$qc_core_pass <- !region_cells$qc_review_flag
  region_cells$segmentation_multiplet_flag <- FALSE
  region_cells$high_control_flag <- FALSE
  region_cells <- assign_spatial_grid(region_cells, 20)
  region_cells$local_density <- calculate_knn_density(region_cells, 2L, "LOCAL_SUBSET")
  region_cells$dense_aggregate <- region_cells$local_density >= stats::quantile(region_cells$local_density, 0.90)
  region_global <- test_spatial_flag_clustering(region_cells, 2L, 19L, 20260814L, "LOCAL_SUBSET")
  region_global$region_id <- region
  region_enrichment <- summarise_spatial_enrichment(region_cells); region_enrichment$region_id <- region
  region_hotspots <- find_spatial_qc_hotspots(region_cells, 19L, 20L, 0.05, 20260814L)
  region_hotspots$region_id <- rep(region, nrow(region_hotspots))
  region_gene_quality <- concordance_genes[concordance_genes$region_id == region, , drop = FALSE]
  region_gene_quality$raw_counts <- region_gene_quality$counts_per_10000
  region_gene_quality$detected_cells <- round(region_gene_quality$detection_fraction * nrow(region_cells))
  region_gene_quality$matrix_cells <- nrow(region_cells)
  region_gene_quality$transcript_rows <- NA_integer_
  region_gene_quality$mean_qv <- NA_real_
  region_gene_quality$fraction_q20 <- NA_real_
  region_gene_quality$represented_codewords <- NA_integer_
  region_gene_quality$transcript_status <- "NOT_RUN_LOCAL_SUBSET"
  region_alarm <- build_cycle_alarm_evidence(if (region == "Region_3") alarm_fixture[0, ] else alarm_fixture, region)
  region_plots <- plot_extended_spatial_qc(region_cells, region_enrichment, region_hotspots, region)
  write_extended_section_artifacts(
    project_root = dirname(extended_slide_root), output_dir = file.path(extended_slide_root, "sections", region),
    region_id = region, mode = "LOCAL_SUBSET", preflight = preflight,
    cycle_alarm_evidence = region_alarm, gene_quality = region_gene_quality,
    spatial_global = region_global, spatial_edge_density = region_enrichment,
    spatial_hotspots = region_hotspots, spatial_cells = region_cells,
    manual_review_manifest = region_hotspots, plots = region_plots
  )
  require_package("Matrix")
  section_counts <- Matrix::Matrix(
    matrix(seq_len(20L * nrow(region_cells)), nrow = 20L,
           dimnames = list(paste0("Gene", seq_len(20L)), region_cells$cell_id)),
    sparse = TRUE
  )
  saveRDS(
    list(counts = section_counts, cells = region_cells,
         features = data.frame(gene = rownames(section_counts)), region_id = region,
         raw_counts_preserved = TRUE),
    file.path(extended_slide_root, "sections", region, paste0(region, ".phase0_2_qc.rds")),
    compress = FALSE
  )
}
extended_coverage <- validate_four_extended_section_outputs(extended_slide_root)
stopifnot(identical(extended_coverage$region_id, paste0("Region_", 1:4)))
stopifnot(length(unique(extended_coverage$mode)) == 1L, extended_coverage$mode[[1]] == "LOCAL_SUBSET")
extended_slide_data <- read_extended_slide_qc_outputs(extended_slide_root)
stopifnot(nrow(extended_slide_data$gene_quality) == 80L)
stopifnot(nrow(extended_slide_data$spatial_cells) == 12L)
stopifnot(identical(unique(extended_slide_data$gene_quality$region_id), paste0("Region_", 1:4)))

extended_slide_summary <- summarise_extended_slide_qc(
  extended_slide_data = extended_slide_data,
  section_summary = concordance_summary,
  manifest = manifest_fixture,
  config = config,
  subset_reference = subset_reference_fixture
)
stopifnot(all(extended_slide_summary$ranking$comparison_status == "NOT_RUN_LOCAL_SUBSET"))
stopifnot(nrow(extended_slide_summary$concordance$summary) == 2L)
stopifnot(all(extended_slide_summary$candidates$candidate_status == "CANDIDATE_NOT_CONFIRMED"))
extended_slide_plots <- plot_extended_slide_qc(extended_slide_data, extended_slide_summary)
stopifnot(all(c("alarm_evidence", "candidate_genes", "ranking", "spatial", "concordance") %in% names(extended_slide_plots)))
extended_slide_paths <- write_extended_slide_qc_artifacts(
  project_root = dirname(extended_slide_root), run_root = extended_slide_root,
  extended_slide_data = extended_slide_data, extended_slide_summary = extended_slide_summary,
  plots = extended_slide_plots
)
stopifnot(all(file.exists(extended_slide_paths)), validate_extended_slide_qc_artifacts(extended_slide_root))
affected_only <- utils::read.delim(
  file.path(extended_slide_root, "slide_summary", "candidate_cycle_affected_genes_affected_only.tsv"),
  check.names = FALSE
)
stopifnot(nrow(affected_only) > 0L, all(affected_only$section_candidate_flag))

# Evidence-only slide artifacts and downstream bundles have an explicit,
# reloadable contract while biological-analysis gates remain pending.
eos_release_fixture <- data.frame(
  gene = paste0("Gene", 1:4),
  gene_set = c("common", "common", "short_lived", "long_lived"),
  stringsAsFactors = FALSE
)
evidence_release <- summarise_evidence_only_qc(
  extended_slide_data = extended_slide_data,
  candidates = extended_slide_summary$candidates,
  eos_gene_sets = eos_release_fixture,
  run_label = "unit_evidence_release", execution_mode = "LOCAL_SUBSET",
  provenance = "unit_fixture"
)
stopifnot(all(c("section_input_cells", "section_primary_include_cells", "qc_threshold_source") %in% names(evidence_release$cell_masks)))
stopifnot(all(c("expected_full_raw_count", "expected_full_provisional_count", "expected_full_technical_risk_count") %in% names(evidence_release$genes)))
stopifnot(all(c("expected_full_retained_total", "expected_full_gene_set_count") %in% names(evidence_release$eos)))
stopifnot(all(c("fdr_threshold", "min_bin_cells_threshold", "permutations", "seed") %in% names(evidence_release$hotspots)))
required_release <- c(
  "cell_downstream_masks.tsv.gz", "section_downstream_decision.tsv",
  "gene_downstream_decision.tsv", "eos_gene_decision_summary.tsv",
  "hotspot_sensitivity_decision.tsv", "evidence_only_qc_release.tsv"
)
stopifnot(all(required_release %in% evidence_only_required_artifacts()))
release_paths <- write_evidence_only_qc_artifacts(
  project_root = dirname(extended_slide_root), run_root = extended_slide_root,
  evidence_summary = evidence_release
)
stopifnot(all(file.exists(release_paths)))
stopifnot(all(file.exists(file.path(
  extended_slide_root, "downstream_inputs", paste0("Region_", 1:4, ".downstream_input.rds")
))))
stopifnot(validate_evidence_only_qc_artifacts(extended_slide_root, stop_on_error = TRUE))
stopifnot(all(evidence_release$release$gate_status %in% c("PASS", "STOP", "PENDING_DOWNSTREAM_ANALYSIS")))
stopifnot(evidence_release$release$gate_status[evidence_release$release$gate_id == "overall_primary_release"] == "PENDING_DOWNSTREAM_ANALYSIS")
region4_bundle <- readRDS(file.path(extended_slide_root, "downstream_inputs", "Region_4.downstream_input.rds"))
stopifnot(region4_bundle$section_status == "SENSITIVITY_ONLY")
stopifnot(region4_bundle$downstream_contract == "MAP_TO_REGION_1_3_REFERENCE_WITH_UNCERTAIN")
stopifnot(ncol(region4_bundle$counts) == nrow(region4_bundle$cell_metadata))

status_path <- file.path(extended_slide_root, "sections", "Region_4", "extended_qc_status.tsv")
status_fixture <- utils::read.delim(status_path, check.names = FALSE)
status_fixture$mode <- "FULL_HPC"
utils::write.table(status_fixture, status_path, sep = "\t", quote = FALSE, row.names = FALSE)
expect_error(validate_four_extended_section_outputs(extended_slide_root), "mixed")
status_fixture$mode <- "LOCAL_SUBSET"
utils::write.table(status_fixture, status_path, sep = "\t", quote = FALSE, row.names = FALSE)
dir.create(file.path(extended_slide_root, "sections", "Region_5"), recursive = TRUE)
expect_error(validate_four_extended_section_outputs(extended_slide_root), "exactly four")

cat("All extended scWAT Xenium QC tests passed.\n")
