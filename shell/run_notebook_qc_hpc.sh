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

command -v Rscript >/dev/null 2>&1 || { echo "Rscript is not available; load the site R module." >&2; exit 3; }
command -v python3 >/dev/null 2>&1 || { echo "python3 is required for notebook parameter injection." >&2; exit 3; }
command -v jupyter >/dev/null 2>&1 || { echo "jupyter is required for notebook execution." >&2; exit 3; }
[[ -d "${INPUT_ROOT}" ]] || { echo "Missing Xenium input root: ${INPUT_ROOT}" >&2; exit 4; }
[[ -f "${PIPELINE_REPO}/R/source.R" ]] || { echo "Missing repository source.R" >&2; exit 4; }
df -h "${PROJECT_ROOT}"

Rscript -e 'p <- c("IRkernel","Matrix","jsonlite","ggplot2"); ok <- vapply(p, requireNamespace, logical(1), quietly=TRUE); if(!all(ok)) stop("Missing R packages: ", paste(p[!ok], collapse=", "))'
jupyter kernelspec list | grep -E '(^|[[:space:]])ir([[:space:]]|$)' >/dev/null || { echo "The Jupyter R kernelspec 'ir' is not registered." >&2; exit 5; }

for region in Region_1 Region_2 Region_3 Region_4; do
  count=$(find "${INPUT_ROOT}" -mindepth 1 -maxdepth 1 -type d -name "*__${region}__*" | wc -l)
  [[ "${count}" -eq 1 ]] || { echo "Expected one input directory for ${region}; found ${count}." >&2; exit 6; }
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
    --set "SEED=20260814L" --set "STRICT_MODE=FALSE"
  jupyter nbconvert --execute --to notebook --inplace \
    --ExecutePreprocessor.kernel_name=ir --ExecutePreprocessor.timeout=-1 "${executed}"
  Rscript -e "source('${PIPELINE_REPO}/R/source.R'); stopifnot(validate_section_artifacts('${RUN_ROOT}/sections/${region}', '${region}'))"
done

summary_executed="${EXECUTED_ROOT}/slide_summary.executed.ipynb"
python3 "${RENDERER}" --inject "${SUMMARY_TEMPLATE}" "${summary_executed}" \
  --set "PROJECT_ROOT=${PROJECT_ROOT}" --set "PIPELINE_REPO=${PIPELINE_REPO}" \
  --set "RUN_LABEL=${RUN_LABEL}" --set "EXPECTED_SECTION_COUNT=4L"
jupyter nbconvert --execute --to notebook --inplace \
  --ExecutePreprocessor.kernel_name=ir --ExecutePreprocessor.timeout=-1 "${summary_executed}"
Rscript -e "x <- readRDS('${RUN_ROOT}/slide_summary/slide_qc_summary.rds'); stopifnot(nrow(x[['data']][['cell_metadata']]) > 0L, length(unique(x[['data']][['cell_metadata']][['region_id']])) == 4L)"
echo "Completed scWAT notebook QC: ${RUN_ROOT}"
