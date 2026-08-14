#!/usr/bin/env bash
#SBATCH --job-name=scwat_xenium_qc
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=24:00:00
#SBATCH --chdir=/dssg/home/acct-svetoslav_chakarov/svetoslav_chakarov/Lab_members/Yanan_Hu
#SBATCH --output=/dssg/home/acct-svetoslav_chakarov/svetoslav_chakarov/Lab_members/Yanan_Hu/adipose_analysis/scwat_qc_logs/%x_%j.out
#SBATCH --error=/dssg/home/acct-svetoslav_chakarov/svetoslav_chakarov/Lab_members/Yanan_Hu/adipose_analysis/scwat_qc_logs/%x_%j.err

set -euo pipefail
export PROJECT_ROOT=/dssg/home/acct-svetoslav_chakarov/svetoslav_chakarov/Lab_members/Yanan_Hu
bash "${PROJECT_ROOT}/adipose_analysis/YNH_Xenium_scWAT/shell/run_notebook_qc_hpc.sh"
