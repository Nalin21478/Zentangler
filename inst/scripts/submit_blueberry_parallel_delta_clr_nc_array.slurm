#!/bin/bash
#SBATCH --job-name=Blueberry_Delta_CLR_NC
#SBATCH --partition=scu-cpu
#SBATCH --output=Blueberry_Delta_CLR_NC_%A_%a.out
#SBATCH --error=Blueberry_Delta_CLR_NC_%A_%a.err
#SBATCH --time=150:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1

source /athena/mallicklab/scratch/naa4050/miniconda3/etc/profile.d/conda.sh
conda activate r-env

PROJECT_DIR=${PROJECT_DIR:-/home/fs01/naa4050/Zentangler_Blueberry}
RESULTS_RDATA=${RESULTS_RDATA:-${PROJECT_DIR}/Parallel_delta_CLR.RData}
GRID_FILE=${GRID_FILE:-${PROJECT_DIR}/blueberry_parallel_delta_clr_nc_grid.csv}
OUT_DIR=${OUT_DIR:-${PROJECT_DIR}/blueberry_parallel_delta_clr_nc_results}
MAE_DIR=${MAE_DIR:-${PROJECT_DIR}/data}
TASK_OFFSET=${TASK_OFFSET:-0}
N_PERM=${N_PERM:-100}
Q_THRESHOLD=${Q_THRESHOLD:-0.25}
SEED=${SEED:-1}
GRID_SCRIPT=${GRID_SCRIPT:-${PROJECT_DIR}/inst/scripts/make_blueberry_parallel_delta_clr_nc_grid.R}
TASK_SCRIPT=${TASK_SCRIPT:-${PROJECT_DIR}/inst/scripts/run_blueberry_parallel_delta_clr_nc_hpc_task.R}

which R
R --version
which Rscript

echo "SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"
echo "TASK_OFFSET=${TASK_OFFSET}"
echo "Effective grid row=$((SLURM_ARRAY_TASK_ID + TASK_OFFSET))"
echo "RESULTS_RDATA=${RESULTS_RDATA}"
echo "GRID_FILE=${GRID_FILE}"
echo "OUT_DIR=${OUT_DIR}"
echo "MAE_DIR=${MAE_DIR}"
echo "N_PERM=${N_PERM}"
echo "Q_THRESHOLD=${Q_THRESHOLD}"
echo "SEED=${SEED}"
echo "GRID_SCRIPT=${GRID_SCRIPT}"
echo "TASK_SCRIPT=${TASK_SCRIPT}"

mkdir -p "${OUT_DIR}"

if [ ! -f "${GRID_FILE}" ]; then
  Rscript "${GRID_SCRIPT}" "${RESULTS_RDATA}" "${GRID_FILE}"
fi

Rscript "${TASK_SCRIPT}" \
  --grid "${GRID_FILE}" \
  --task-id "${SLURM_ARRAY_TASK_ID}" \
  --task-offset "${TASK_OFFSET}" \
  --out-dir "${OUT_DIR}" \
  --mae-dir "${MAE_DIR}" \
  --n-perm "${N_PERM}" \
  --q-threshold "${Q_THRESHOLD}" \
  --seed "${SEED}"
