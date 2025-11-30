#!/bin/bash
#SBATCH --job-name=GENESPACE
#SBATCH --partition=pibu_el8
#SBATCH --time=1-00:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --output=/data/users/%u/annotation_of_eukaryotic_genome/annotation/log/GENESPACE_%j.out
#SBATCH --error=/data/users/%u/annotation_of_eukaryotic_genome/annotation/log/GENESPACE_%j.err

set -Eeuo pipefail

#----- CONFIG -----
WD="/data/users/${USER}/annotation_of_eukaryotic_genome/genespace_workingDirectory"
RSCRIPT="/data/users/${USER}/annotation_of_eukaryotic_genome/scripts/08c_genespace_run.R"
COURSEDIR="/data/courses/assembly-annotation-course/CDS_annotation"
IMG="${COURSEDIR}/containers/genespace_latest.sif"

[ -d "$WD" ] || { echo "[ERROR] Working dir missing: $WD"; exit 2; }
[ -s "$RSCRIPT" ] || { echo "[ERROR] R script missing: $RSCRIPT"; exit 3; }
[ -s "$IMG" ] || { echo "[ERROR] Container missing: $IMG"; exit 4; }

export SLURM_CPUS_PER_TASK="${SLURM_CPUS_PER_TASK:-16}"
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export LANG=C.UTF-8
export LC_ALL=C.UTF-8

#----- RUN R SCRIPT GENESPACE -----
echo "[INFO] Running GENESPACE in container..."
apptainer exec --bind /data "$IMG" Rscript "$RSCRIPT" "$WD" diamond orthofinder "/usr/local/bin"

echo "[OK] Finished."
