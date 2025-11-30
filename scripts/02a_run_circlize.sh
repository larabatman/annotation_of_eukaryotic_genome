#!/usr/bin/env bash
#SBATCH --job-name=circlize_te_all
#SBATCH --partition=pibu_el8
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --output=%x_%j.out

set -euo pipefail

#----- CONFIG -----
WORKDIR="/data/users/lland/annotation_of_eukaryotic_genome"
FASTA="${WORKDIR}/assembly/hifiasm/istisu1.bp.p_ctg.fa"
FAI="${FASTA}.fai"
GFF="${WORKDIR}/EDTA_annotation/istisu1.bp.p_ctg.fa.mod.EDTA.TEanno.gff3"
OUTDIR="${WORKDIR}/EDTA_annotation/circlize_all"
RSCRIPT_PATH="${WORKDIR}/scripts/02b_circlize_te_density.R"

mkdir -p "$OUTDIR"

#----- OPTIONAL SAMtools faidx FOR FAI -----
module purge || true
module load SAMtools || true
if ! command -v samtools >/dev/null 2>&1; then
  echo "[ERROR] samtools not available." >&2; exit 2
fi
[ -s "$FAI" ] || samtools faidx "$FASTA"

#----- LAUNCH RSCRIPT -----
module load R/4.2.1-foss-2021a || true
command -v Rscript >/dev/null || { echo "[ERROR] Rscript not found."; exit 3; }

Rscript "$RSCRIPT_PATH"
echo "[OK] circlize done -> $OUTDIR"
