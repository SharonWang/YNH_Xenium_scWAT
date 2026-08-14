#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-/dssg/home/acct-svetoslav_chakarov/svetoslav_chakarov/Lab_members/Yanan_Hu}"
PIPELINE_REPO="${PIPELINE_REPO:-${PROJECT_ROOT}/adipose_analysis/YNH_Xenium_scWAT}"
INPUT_ROOT="${INPUT_ROOT:-${PROJECT_ROOT}/adipose_data}"
RUN_LABEL="${RUN_LABEL:-full_notebook_qc_$(date +%Y%m%d_%H%M%S)}"
METADATA_PATH="${METADATA_PATH:-${PIPELINE_REPO}/config/scwat_sample_manifest.tsv}"
RUN_ROOT="${PROJECT_ROOT}/adipose_analysis/scwat_qc_outputs/${RUN_LABEL}"
TEMP_ROOT="${PROJECT_ROOT}/adipose_analysis/tmp"
R_LIBRARY_ROOT="${PROJECT_ROOT}/adipose_analysis/R_libs"
EXECUTED_ROOT="${RUN_ROOT}/executed_notebooks"

for path in "${PIPELINE_REPO}" "${INPUT_ROOT}" "${RUN_ROOT}" "${TEMP_ROOT}" "${R_LIBRARY_ROOT}"; do
  case "${path}" in "${PROJECT_ROOT}"|"${PROJECT_ROOT}"/*) ;; *) echo "Unsafe path outside PROJECT_ROOT: ${path}" >&2; exit 2;; esac
done
if [[ -n "${METADATA_PATH}" ]]; then
  case "${METADATA_PATH}" in "${PROJECT_ROOT}"/*) ;; *) echo "METADATA_PATH must be below PROJECT_ROOT" >&2; exit 2;; esac
  [[ -f "${METADATA_PATH}" ]] || { echo "Missing metadata manifest: ${METADATA_PATH}" >&2; exit 4; }
fi

mkdir -p "${TEMP_ROOT}" "${R_LIBRARY_ROOT}" "${EXECUTED_ROOT}" "${PROJECT_ROOT}/adipose_analysis/scwat_qc_logs"
export TMPDIR="${TEMP_ROOT}"
export TMP="${TEMP_ROOT}"
export TEMP="${TEMP_ROOT}"
export R_LIBS_USER="${R_LIBRARY_ROOT}"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export PYTHONDONTWRITEBYTECODE=1

command -v Rscript >/dev/null 2>&1 || { echo "Rscript is not available; load the site R module." >&2; exit 3; }
command -v python3 >/dev/null 2>&1 || { echo "python3 is required for notebook parameter injection." >&2; exit 3; }
command -v jupyter >/dev/null 2>&1 || { echo "jupyter is required for notebook execution." >&2; exit 3; }
[[ -d "${INPUT_ROOT}" ]] || { echo "Missing Xenium input root: ${INPUT_ROOT}" >&2; exit 4; }
[[ -f "${PIPELINE_REPO}/R/source.R" ]] || { echo "Missing repository source.R" >&2; exit 4; }
echo "Project filesystem capacity:"
df -h "${PROJECT_ROOT}"

Rscript -e 'p <- c("IRkernel","Matrix","jsonlite","ggplot2","arrow","dplyr","RANN"); ok <- vapply(p, requireNamespace, logical(1), quietly=TRUE); if(!all(ok)) stop("Missing R packages: ", paste(p[!ok], collapse=", ")); cat("Full-HPC R package preflight passed:\n", paste(sprintf("%s=%s", p, vapply(p, function(x) as.character(packageVersion(x)), character(1))), collapse="\n"), "\n")'
jupyter kernelspec list | grep -E '(^|[[:space:]])ir([[:space:]]|$)' >/dev/null || { echo "The Jupyter R kernelspec 'ir' is not registered." >&2; exit 5; }

for region in Region_1 Region_2 Region_3 Region_4; do
  mapfile -t region_dirs < <(find "${INPUT_ROOT}" -mindepth 1 -maxdepth 1 -type d -name "*__${region}__*" -print)
  [[ "${#region_dirs[@]}" -eq 1 ]] || { echo "Expected one input directory for ${region}; found ${#region_dirs[@]}." >&2; exit 6; }
  region_dir="${region_dirs[0]}"
  [[ -f "${region_dir}/transcripts.parquet" ]] || { echo "FULL_HPC requires ${region_dir}/transcripts.parquet" >&2; exit 7; }
  echo "${region} projected input sizes:"
  du -sh "${region_dir}/transcripts.parquet" "${region_dir}/cell_feature_matrix" "${region_dir}/cells.csv.gz"
done

RENDERER="${PIPELINE_REPO}/scripts/render_notebooks.py"
SECTION_TEMPLATE="${PIPELINE_REPO}/notebooks/01_section_phase0_2_QC.ipynb"
SUMMARY_TEMPLATE="${PIPELINE_REPO}/notebooks/02_slide_QC_summary.ipynb"
python3 "${RENDERER}" --validate "${SECTION_TEMPLATE}" --type section
python3 "${RENDERER}" --validate "${SUMMARY_TEMPLATE}" --type summary

for region in Region_1 Region_2 Region_3 Region_4; do
  executed="${EXECUTED_ROOT}/${region}.executed.ipynb"
  python3 "${RENDERER}" --inject "${SECTION_TEMPLATE}" "${executed}" \
    --set "PROJECT_ROOT=${PROJECT_ROOT}" --set "PIPELINE_REPO=${PIPELINE_REPO}" \
    --set "INPUT_ROOT=${INPUT_ROOT}" --set "REGION_ID=${region}" --set "RUN_LABEL=${RUN_LABEL}" \
    --set "METADATA_PATH=${METADATA_PATH}" --set "EXPECTED_SECTION_COUNT=4L" \
    --set "SEED=20260814L" --set "STRICT_MODE=FALSE" \
    --set "EXTENDED_QC_MODE=FULL_HPC" \
    --set "EXTENDED_QC_CONFIG_PATH=${PIPELINE_REPO}/config/extended_qc_defaults.tsv"
  jupyter nbconvert --execute --to notebook --inplace \
    --ExecutePreprocessor.kernel_name=ir --ExecutePreprocessor.timeout=-1 "${executed}"
  Rscript -e "source('${PIPELINE_REPO}/R/source.R'); out <- '${RUN_ROOT}/sections/${region}'; stopifnot(validate_section_artifacts(out, '${region}')); stopifnot(validate_extended_section_artifacts(out, '${region}', 'FULL_HPC', stop_on_error=TRUE)); e <- read_extended_section_artifacts(out, '${region}', 'FULL_HPC'); stopifnot(!any(e[['gene_quality']][['transcript_status']] == 'NOT_RUN_LOCAL_SUBSET'), nrow(e[['gene_quality']]) > 0L, nrow(e[['spatial_cells']]) > 0L, e[['status']][['cells_deleted']] == 0L)"
done

summary_executed="${EXECUTED_ROOT}/slide_summary.executed.ipynb"
python3 "${RENDERER}" --inject "${SUMMARY_TEMPLATE}" "${summary_executed}" \
  --set "PROJECT_ROOT=${PROJECT_ROOT}" --set "PIPELINE_REPO=${PIPELINE_REPO}" \
  --set "RUN_LABEL=${RUN_LABEL}" --set "EXPECTED_SECTION_COUNT=4L" \
  --set "METADATA_PATH=${METADATA_PATH}" \
  --set "EXTENDED_QC_CONFIG_PATH=${PIPELINE_REPO}/config/extended_qc_defaults.tsv" \
  --set "SUBSET_REFERENCE_PATH=${PIPELINE_REPO}/config/subset_qc_reference.tsv"
jupyter nbconvert --execute --to notebook --inplace \
  --ExecutePreprocessor.kernel_name=ir --ExecutePreprocessor.timeout=-1 "${summary_executed}"
Rscript -e "source('${PIPELINE_REPO}/R/source.R'); x <- readRDS('${RUN_ROOT}/slide_summary/slide_qc_summary.rds'); stopifnot(nrow(x[['data']][['cell_metadata']]) > 0L, length(unique(x[['data']][['cell_metadata']][['region_id']])) == 4L); stopifnot(validate_extended_slide_qc_artifacts('${RUN_ROOT}', stop_on_error=TRUE)); root <- '${RUN_ROOT}/slide_summary'; alarm <- read.delim(file.path(root, 'combined_cycle_alarm_evidence.tsv'), check.names=FALSE); gene <- read.delim(file.path(root, 'combined_gene_transcript_quality.tsv'), check.names=FALSE); candidate <- read.delim(file.path(root, 'candidate_cycle_affected_genes.tsv'), check.names=FALSE); ranking <- read.delim(file.path(root, 'subset_full_qc_ranking.tsv'), check.names=FALSE); concordance <- read.delim(file.path(root, 'within_mouse_section_concordance.tsv'), check.names=FALSE); direct <- alarm[['region_id']][alarm[['evidence_status']] == 'DIRECT_EVIDENCE']; stopifnot(identical(direct, c('Region_1','Region_2','Region_4')), nrow(gene) > 0L, !any(gene[['transcript_status']] == 'NOT_RUN_LOCAL_SUBSET'), nrow(candidate) > 0L, all(candidate[['candidate_status']] == 'CANDIDATE_NOT_CONFIRMED'), all(candidate[['exact_cycle_status']] == 'REQUIRES_10X_DIAGNOSTICS'), all(ranking[['comparison_status']] == 'FULL_DATA_COMPARISON'), nrow(concordance) == 2L)"
echo "Completed scWAT notebook QC: ${RUN_ROOT}"
