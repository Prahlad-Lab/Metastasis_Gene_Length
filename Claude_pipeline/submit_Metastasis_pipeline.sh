#!/bin/bash
#SBATCH --job-name=metastasis_pipeline
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=logs/metastasis_pipeline_%j.out
#SBATCH --error=logs/metastasis_pipeline_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --partition=general-compute
#SBATCH --cluster=ub-hpc
#SBATCH --qos=general-compute
#SBATCH --mail-user=Johnny.CruzCorchado@RoswellPark.org
#SBATCH --mail-type=END,FAIL

# ==============================================================================
# 1. Setup Environment
# ==============================================================================
# Create the logs directory if it doesn't exist so Slurm doesn't crash
mkdir -p logs

# Activate your Conda/Mamba environment that contains Nextflow and your Python packages
# (Make sure to initialize your shell first if running in a batch script)
source $(conda info --base)/etc/profile.d/conda.sh
conda activate nextflow.25

echo "Nextflow Manager script completed."
WORK_DIR="/vscratch/grp-vprahlad/Metastasis_Gene_Length/work"
mkdir -p "$WORK_DIR"


# ==============================================================================
# 2. Execute Pipeline
# ==============================================================================
echo "Starting Nextflow Metastasis Pipeline..."

# Tell Nextflow to run main.nf, use Slurm to submit sub-tasks, and resume if interrupted
nextflow run main.nf -profile mixed \
	-w "$WORK_DIR" \
	-with-report logs/report_STAR_${SLURM_JOB_ID}.html \
	-with-trace logs/trace_STAR_${SLURM_JOB_ID}.txt \
	-with-timeline logs/timeline_STAR_${SLURM_JOB_ID}.html 

echo "Nextflow pipeline completed."
