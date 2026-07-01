#!/bin/bash
# =============================================================================
# run_pipeline.sh
# SLURM job script for the MACSima MICS pipeline
#
# Usage (submit from the project root directory):
#   sbatch scripts/run_pipeline.sh              # full pipeline
#   sbatch scripts/run_pipeline.sh --from 3     # restart from step 3
#   sbatch scripts/run_pipeline.sh --steps 5,6,7
#
# Adjust #SBATCH directives to match your cluster's partition and resources.
# The conda environment name and activation path may also need to be adapted.
# =============================================================================

#SBATCH --job-name=MACSima_pipeline
#SBATCH --output=logs/slurm_%j.log
#SBATCH --error=logs/slurm_%j.log
#SBATCH --time=8-00:00:00        # walltime: adjust to your cluster's partition limit
#SBATCH --mem=192G               # RAM: FindNeighbors on ~5M cells requires >=128G
#SBATCH --cpus-per-task=8        # Harmony and UMAP benefit from multi-threading
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-type=BEGIN,END,FAIL
# #SBATCH --mail-user=your.email@institution.edu   # uncomment and set your email

# --------------------------------------------------------------------------
# Conda environment (R + Seurat installed in seurat_env)
# Adapt the path below to your cluster's conda installation.
# --------------------------------------------------------------------------

# source /path/to/conda/etc/profile.d/conda.sh
source /lustre/software/anaconda/2022.10_all/etc/profile.d/conda.sh
conda activate seurat_env

# --------------------------------------------------------------------------
# Variables — derived from submission directory, no hardcoded paths
# --------------------------------------------------------------------------

# SLURM_SUBMIT_DIR is set automatically to the directory from which sbatch
# was called. Always submit from the project root: cd MACSima_pipeline && sbatch ...
PROJECT_DIR="${SLURM_SUBMIT_DIR}"
SCRIPT="$PROJECT_DIR/scripts/run_pipeline.R"
LOG_DIR="$PROJECT_DIR/logs"

mkdir -p "$LOG_DIR"

# Optional: pass MACSIMA_DIR to R so config.R picks it up via Sys.getenv()
export MACSIMA_DIR="$PROJECT_DIR"

# Additional arguments forwarded to run_pipeline.R (e.g. --from 3)
EXTRA_ARGS="$@"

# --------------------------------------------------------------------------
# Environment diagnostics
# --------------------------------------------------------------------------

echo "============================================="
echo "MACSima pipeline — SLURM job $SLURM_JOB_ID"
echo "  Date:         $(date)"
echo "  Node:         $SLURMD_NODENAME"
echo "  CPUs:         $SLURM_CPUS_PER_TASK"
echo "  RAM (requested): ${SLURM_MEM_PER_NODE}MB"
echo "  Project dir:  $PROJECT_DIR"
echo "============================================="

R --version | head -1

# --------------------------------------------------------------------------
# Multi-threaded library settings
# --------------------------------------------------------------------------

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export BLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK

# --------------------------------------------------------------------------
# Launch pipeline
# --------------------------------------------------------------------------

echo "Running: Rscript $SCRIPT $EXTRA_ARGS"
echo ""

Rscript "$SCRIPT" $EXTRA_ARGS
EXIT_CODE=$?

echo ""
echo "============================================="
if [ $EXIT_CODE -eq 0 ]; then
  echo "Pipeline completed successfully."
else
  echo "Pipeline FAILED (exit code: $EXIT_CODE)."
  echo "Check: $LOG_DIR/slurm_${SLURM_JOB_ID}.log"
fi
echo "Date: $(date)"
echo "============================================="

exit $EXIT_CODE
