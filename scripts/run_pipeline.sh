#!/bin/bash
# =============================================================================
# run_pipeline.sh
# SLURM job script per la pipeline MACSima
#
# Uso:
#   sbatch scripts/run_pipeline.sh           # pipeline completa
#   sbatch scripts/run_pipeline.sh --from 3  # riparte dallo step 3
#   sbatch scripts/run_pipeline.sh --steps 5,6,7
#
# Adatta le righe #SBATCH in base alla coda/partizione del tuo cluster.
# =============================================================================

#SBATCH --job-name=MACSima_pipeline
#SBATCH --output=/lustre/home/gfiscon/projects/MACSima_pipeline/logs/slurm_%j.log
#SBATCH --error=/lustre/home/gfiscon/projects/MACSima_pipeline/logs/slurm_%j.log
#SBATCH --time=12:00:00          # walltime: aggiusta in base al cluster
#SBATCH --mem=192G               # RAM: FindNeighbors su ~4M cellule richiede >=128G
#SBATCH --cpus-per-task=8        # Harmony e UMAP parallelizzano su thread
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=federicadannunzio@gmail.com

# --------------------------------------------------------------------------
# Ambiente conda (R 4.4 + Seurat installati in seurat_env)
# --------------------------------------------------------------------------

source /lustre/software/anaconda/2022.10_all/etc/profile.d/conda.sh
conda activate seurat_env

# --------------------------------------------------------------------------
# Variabili
# --------------------------------------------------------------------------

PROJECT_DIR="/lustre/home/gfiscon/projects/MACSima_pipeline"
SCRIPT="$PROJECT_DIR/scripts/run_pipeline.R"
LOG_DIR="$PROJECT_DIR/logs"

mkdir -p "$LOG_DIR"

# Passa eventuali argomenti aggiuntivi (es. --from 3) a run_pipeline.R
EXTRA_ARGS="$@"

# --------------------------------------------------------------------------
# Diagnostica ambiente
# --------------------------------------------------------------------------

echo "============================================="
echo "MACSima pipeline — SLURM job $SLURM_JOB_ID"
echo "  Data:        $(date)"
echo "  Nodo:        $SLURMD_NODENAME"
echo "  CPUs:        $SLURM_CPUS_PER_TASK"
echo "  RAM (richiesta): ${SLURM_MEM_PER_NODE}MB"
echo "  Working dir: $PROJECT_DIR"
echo "============================================="

R --version | head -1

# --------------------------------------------------------------------------
# Imposta numero di thread per librerie multi-threaded
# (Harmony, data.table, openblas)
# --------------------------------------------------------------------------

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export BLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK

# --------------------------------------------------------------------------
# Lancia la pipeline
# --------------------------------------------------------------------------

echo "Avvio: Rscript $SCRIPT $EXTRA_ARGS"
echo ""

Rscript "$SCRIPT" $EXTRA_ARGS
EXIT_CODE=$?

echo ""
echo "============================================="
if [ $EXIT_CODE -eq 0 ]; then
  echo "Pipeline completata con successo."
else
  echo "Pipeline FALLITA (exit code: $EXIT_CODE)."
  echo "Controlla: $LOG_DIR/slurm_${SLURM_JOB_ID}.err"
fi
echo "Fine: $(date)"
echo "============================================="

exit $EXIT_CODE
