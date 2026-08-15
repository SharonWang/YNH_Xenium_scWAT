#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-/dssg/home/acct-svetoslav_chakarov/svetoslav_chakarov/Lab_members/Yanan_Hu/YNH_Xenium}"
PIPELINE_REPO="${PIPELINE_REPO:-${PROJECT_ROOT}/adipose_analysis/YNH_Xenium_scWAT}"
INPUT_ROOT="${INPUT_ROOT:-${PROJECT_ROOT}/adipose_data}"
RUN_LABEL="${RUN_LABEL:-full_notebook_qc_$(date +%Y%m%d_%H%M%S)}"
REGION_ID="${REGION_ID:-ALL}"
METADATA_PATH="${METADATA_PATH:-${PIPELINE_REPO}/config/scwat_sample_manifest.tsv}"
RUN_ROOT="${PROJECT_ROOT}/adipose_analysis/scwat_qc_outputs/${RUN_LABEL}"
TEMP_ROOT="${PROJECT_ROOT}/adipose_analysis/tmp"
R_LIBRARY_ROOT="${PROJECT_ROOT}/adipose_analysis/R_libs"
EXECUTED_ROOT="${RUN_ROOT}/executed_notebooks"

RUN_SUMMARY=0
case "${REGION_ID}" in
  ALL)
    regions=(Region_1 Region_2 Region_3 Region_4)
    RUN_SUMMARY=1
    ;;
  SUMMARY)
    regions=()
    RUN_SUMMARY=1
    ;;
  Region_1|Region_2|Region_3|Region_4)
    regions=("${REGION_ID}")
    ;;
  *)
    echo "REGION_ID must be ALL, SUMMARY, or Region_1 through Region_4; received ${REGION_ID}." >&2
    exit 2
    ;;
esac

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

if [[ "${REGION_ID}" == "SUMMARY" ]]; then
  Rscript -e 'p <- c("IRkernel","Matrix","jsonlite","ggplot2"); ok <- vapply(p, requireNamespace, logical(1), quietly=TRUE); if(!all(ok)) stop("Missing R packages: ", paste(p[!ok], collapse=", ")); cat("Summary R package preflight passed:\n", paste(sprintf("%s=%s", p, vapply(p, function(x) as.character(packageVersion(x)), character(1))), collapse="\n"), "\n")'
else
  Rscript -e 'p <- c("IRkernel","Matrix","jsonlite","ggplot2","arrow","dplyr","RANN"); ok <- vapply(p, requireNamespace, logical(1), quietly=TRUE); if(!all(ok)) stop("Missing R packages: ", paste(p[!ok], collapse=", ")); cat("Full-HPC R package preflight passed:\n", paste(sprintf("%s=%s", p, vapply(p, function(x) as.character(packageVersion(x)), character(1))), collapse="\n"), "\n")'
fi
jupyter kernelspec list | grep -E '(^|[[:space:]])ir([[:space:]]|$)' >/dev/null || { echo "The Jupyter R kernelspec 'ir' is not registered." >&2; exit 5; }

for region in "${regions[@]}"; do
  mapfile -t region_dirs < <(find "${INPUT_ROOT}" -mindepth 1 -maxdepth 1 -type d -name "*__${region}__*" -print)
  [[ "${#region_dirs[@]}" -eq 1 ]] || { echo "Expected one input directory for ${region}; found ${#region_dirs[@]}." >&2; exit 6; }
  region_dir="${region_dirs[0]}"
  [[ -f "${region_dir}/transcripts.parquet" ]] || { echo "FULL_HPC requires ${region_dir}/transcripts.parquet" >&2; exit 7; }
  echo "${region} projected input sizes:"
  du -sh "${region_dir}/transcripts.parquet" "${region_dir}/cell_feature_matrix" "${region_dir}/cells.csv.gz"
done

RENDERER="${PIPELINE_REPO}/scripts/render_notebooks.py"
SUMMARY_TEMPLATE="${PIPELINE_REPO}/notebooks/02_slide_QC_summary.ipynb"
for region in "${regions[@]}"; do
  python3 "${RENDERER}" --validate "$(python3 "${RENDERER}" --region-notebook "${region}")" --type section
done
if [[ "${RUN_SUMMARY}" -eq 1 ]]; then
  python3 "${RENDERER}" --validate "${SUMMARY_TEMPLATE}" --type summary
fi

for region in "${regions[@]}"; do
  section_notebook="$(python3 "${RENDERER}" --region-notebook "${region}")"
  executed="${EXECUTED_ROOT}/${region}.executed.ipynb"
  python3 "${RENDERER}" --inject "${section_notebook}" "${executed}" \
    --set "PROJECT_ROOT=${PROJECT_ROOT}" --set "PIPELINE_REPO=${PIPELINE_REPO}" \
    --set "INPUT_ROOT=${INPUT_ROOT}" --set "REGION_ID=${region}" --set "RUN_LABEL=${RUN_LABEL}" \
    --set "METADATA_PATH=${METADATA_PATH}" --set "EXPECTED_SECTION_COUNT=4L" \
    --set "SEED=20260814L" --set "STRICT_MODE=FALSE" \
    --set "EXTENDED_QC_MODE=FULL_HPC" \
    --set "EXTENDED_QC_CONFIG_PATH=${PIPELINE_REPO}/config/extended_qc_defaults.tsv"
  jupyter nbconvert --execute --to notebook --inplace \
    --ExecutePreprocessor.kernel_name=ir --ExecutePreprocessor.timeout=-1 "${executed}"
  Rscript -e "source('${PIPELINE_REPO}/R/source.R'); out <- '${RUN_ROOT}/sections/${region}'; stopifnot(validate_section_artifacts(out, '${region}')); stopifnot(validate_extended_section_artifacts(out, '${region}', 'FULL_HPC', stop_on_error=TRUE)); e <- read_extended_section_artifacts(out, '${region}', 'FULL_HPC'); m <- read.delim(gzfile(file.path(out,'cell_downstream_masks.tsv.gz')), check.names=FALSE); d <- read.delim(file.path(out,'section_downstream_decision.tsv'), check.names=FALSE); stopifnot(!any(e[['gene_quality']][['transcript_status']] == 'NOT_RUN_LOCAL_SUBSET'), nrow(e[['gene_quality']]) > 0L, nrow(e[['spatial_cells']]) > 0L, e[['status']][['cells_deleted']] == 0L, nrow(m) == nrow(e[['spatial_cells']]), d[['section_status']] == section_downstream_status('${region}'))"
done

if [[ "${RUN_SUMMARY}" -eq 1 ]]; then
  summary_executed="${EXECUTED_ROOT}/slide_summary.executed.ipynb"
  python3 "${RENDERER}" --inject "${SUMMARY_TEMPLATE}" "${summary_executed}" \
    --set "PROJECT_ROOT=${PROJECT_ROOT}" --set "PIPELINE_REPO=${PIPELINE_REPO}" \
    --set "RUN_LABEL=${RUN_LABEL}" --set "EXPECTED_SECTION_COUNT=4L" \
    --set "METADATA_PATH=${METADATA_PATH}" \
    --set "EXTENDED_QC_CONFIG_PATH=${PIPELINE_REPO}/config/extended_qc_defaults.tsv" \
    --set "SUBSET_REFERENCE_PATH=${PIPELINE_REPO}/config/subset_qc_reference.tsv"
  jupyter nbconvert --execute --to notebook --inplace \
    --ExecutePreprocessor.kernel_name=ir --ExecutePreprocessor.timeout=-1 "${summary_executed}"
  Rscript -e "source('${PIPELINE_REPO}/R/source.R'); x <- readRDS('${RUN_ROOT}/slide_summary/slide_qc_summary.rds'); stopifnot(nrow(x[['data']][['cell_metadata']]) > 0L, length(unique(x[['data']][['cell_metadata']][['region_id']])) == 4L); stopifnot(validate_extended_slide_qc_artifacts('${RUN_ROOT}', stop_on_error=TRUE), validate_evidence_only_qc_artifacts('${RUN_ROOT}', stop_on_error=TRUE)); root <- '${RUN_ROOT}/slide_summary'; genes <- read.delim(file.path(root,'gene_downstream_decision.tsv'), check.names=FALSE); eos <- read.delim(file.path(root,'eos_gene_decision_summary.tsv'), check.names=FALSE); release <- read.delim(file.path(root,'evidence_only_qc_release.tsv'), check.names=FALSE); sections <- read.delim(file.path(root,'section_downstream_decision.tsv'), check.names=FALSE); stopifnot(nrow(genes)==479L, sum(genes[['conservative_evidence_status']]=='CONSERVATIVE_NO_SIGNAL_DETECTED')==67L, sum(genes[['primary_feature_status']]=='PROVISIONAL_PRIMARY_FEATURES')==245L, sum(genes[['technical_risk_status']]=='TECHNICAL_RISK_SENSITIVITY_ONLY')==234L, sum(eos[['retained_provisional']])==53L, identical(sections[['section_status']],c('PRIMARY_CONDITIONAL','PRIMARY_CONDITIONAL','PRIMARY','SENSITIVITY_ONLY')), release[['gate_status']][release[['gate_id']]=='overall_primary_release']=='PENDING_DOWNSTREAM_ANALYSIS')"
fi
echo "Completed scWAT notebook QC target ${REGION_ID}: ${RUN_ROOT}"
