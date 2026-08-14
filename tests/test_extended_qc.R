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

# A known four-cell corner cluster should produce a morphology-review hotspot.
hotspot_grid <- assign_spatial_grid(spatial_fixture, grid_size_um = 20)
hotspots_1 <- find_spatial_qc_hotspots(hotspot_grid, permutations = 199L, min_bin_cells = 4L, fdr = 0.10, seed = 20260814L)
hotspots_2 <- find_spatial_qc_hotspots(hotspot_grid, permutations = 199L, min_bin_cells = 4L, fdr = 0.10, seed = 20260814L)
stopifnot(identical(hotspots_1, hotspots_2))
stopifnot(any(hotspots_1$hotspot_status == "MORPHOLOGY_REVIEW_REQUIRED"))

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

cat("All extended scWAT Xenium QC tests passed.\n")
