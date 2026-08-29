#!/bin/bash
#SBATCH --job-name=Blueberry_NC
#SBATCH --partition=scu-cpu
#SBATCH --output=Blueberry_NC_%A_%a.out
#SBATCH --error=Blueberry_NC_%A_%a.err
#SBATCH --time=150:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1

source /athena/mallicklab/scratch/naa4050/miniconda3/etc/profile.d/conda.sh
conda activate r-env

PROJECT_DIR=${PROJECT_DIR:-/home/fs01/naa4050/Zentangler_Blueberry}
GRID_FILE=${GRID_FILE:-${PROJECT_DIR}/blueberry_parallel_negative_control_grid.csv}
OUT_DIR=${OUT_DIR:-${PROJECT_DIR}/blueberry_parallel_negative_control_results}
MAE_DIR=${MAE_DIR:-${PROJECT_DIR}/data}
TASK_OFFSET=${TASK_OFFSET:-0}
N_PERM=${N_PERM:-500}
Q_THRESHOLD=${Q_THRESHOLD:-0.25}
SEED=${SEED:-1}
TASK_SCRIPT=${TASK_SCRIPT:-${PROJECT_DIR}/inst/scripts/run_blueberry_parallel_negative_control_hpc_task.R}

which R
R --version
which Rscript

echo "SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"
echo "TASK_OFFSET=${TASK_OFFSET}"
echo "Effective grid row=$((SLURM_ARRAY_TASK_ID + TASK_OFFSET))"
echo "GRID_FILE=${GRID_FILE}"
echo "OUT_DIR=${OUT_DIR}"
echo "MAE_DIR=${MAE_DIR}"
echo "N_PERM=${N_PERM}"
echo "Q_THRESHOLD=${Q_THRESHOLD}"
echo "SEED=${SEED}"
echo "TASK_SCRIPT=${TASK_SCRIPT}"

mkdir -p "${OUT_DIR}"

Rscript "${TASK_SCRIPT}" \
  --grid "${GRID_FILE}" \
  --task-id "${SLURM_ARRAY_TASK_ID}" \
  --task-offset "${TASK_OFFSET}" \
  --out-dir "${OUT_DIR}" \
  --mae-dir "${MAE_DIR}" \
  --n-perm "${N_PERM}" \
  --q-threshold "${Q_THRESHOLD}" \
  --seed "${SEED}"
