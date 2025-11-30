#!/bin/bash
#SBATCH --job-name=BioHabs_Pipeline
#SBATCH --output=logs/hpc_%j.log
#SBATCH --error=logs/hpc_%j.err
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G

# ==============================================================================
# BioHabs HPC Deployment Script
# ==============================================================================
# Usage: sbatch BioHabs/deploy_hpc.sh
#
# Prerequisites:
# - R (>= 4.0.0)
# - GDAL (>= 3.0)
# - renv (initialized)
# ==============================================================================

# 1. Load Modules
# module load R/4.3.0
# module load gdal/3.6.2
# module load java/11  # Required for H2O

# 2. Set Environment Variables
export R_LIBS_USER=~/R/library
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# 3. Navigate to Project Root
# Assumes script is submitted from project root or BioHabs subdir
if [ -d "BioHabs" ]; then
    PROJECT_ROOT="."
elif [ -f "run_pipeline.R" ]; then
    PROJECT_ROOT=".."
    cd ..
else
    echo "Error: Could not find project root. Submit from project root."
    exit 1
fi

echo "Running BioHabs Pipeline in $PROJECT_ROOT"
echo "Mode: REPRO"

# 4. Restore Environment (renv)
Rscript -e "renv::restore()"

# 5. Execute Pipeline
Rscript BioHabs/run_pipeline.R --mode=REPRO

echo "Pipeline execution finished."
