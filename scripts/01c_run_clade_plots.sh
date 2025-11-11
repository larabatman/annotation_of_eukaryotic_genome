#!/usr/bin/env bash
#SBATCH --job-name=tesorter_plot
#SBATCH --partition=pibu_el8
#SBATCH --time=00:15:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

set -euo pipefail

#----- CONFIG ------
RESULTS_DIR="/data/users/lland/annotation_of_eukaryotic_genome/EDTA_annotation/istisu1.bp.p_ctg.fa.LTR_analysis"
RSCRIPT_PATH="/data/users/lland/annotation_of_eukaryotic_genome/scripts"

GFF="${RESULTS_DIR}/genomic.fna.mod.LTR.intact.raw.gff3"
CLS="${RESULTS_DIR}/genomic.fna.mod.LTR.raw.fa.rexdb-plant.cls.tsv"

#Sanity checks
[ -s "$GFF" ] || { echo "[ERR] Missing GFF: $GFF"; exit 2; }
[ -s "$CLS" ] || { echo "[ERR] Missing TEsorter TSV: $CLS"; exit 3; }
[ -s "$RSCRIPT_PATH" ] || { echo "[ERR] Missing R script: $RSCRIPT_PATH"; exit 4; }

# R env
module purge
module load R/4.2.1-foss-2021a || true
command -v Rscript >/dev/null || { echo "[ERR] Rscript not on PATH"; exit 5; }

# Run
mkdir -p "${RESULTS_DIR}/plots"
cd "$RESULTS_DIR"
echo "[INFO] Running plot script in ${RESULTS_DIR} ..."
Rscript "$RSCRIPT_PATH"
echo "[OK] Plots written to: ${RESULTS_DIR}/plots"
