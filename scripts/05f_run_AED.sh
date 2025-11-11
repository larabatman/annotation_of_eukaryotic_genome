#!/usr/bin/env bash
#SBATCH --job-name=AED
#SBATCH --partition=pibu_el8
#SBATCH --time=00:10:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --output=/data/users/%u/annotation_of_eukaryotic_genome/annotation/log/AED_%j.out
#SBATCH --error=/data/users/%u/annotation_of_eukaryotic_genome/annotation/log/AED_%j.err
set -euo pipefail

FINAL_DIR="/data/users/${USER}/annotation_of_eukaryotic_genome/annotation/final"
GFF_IN="${FINAL_DIR}/istisu1.bp.p_ctg.all.noseq.gff.renamed.functional.gff"
AED_TXT="${FINAL_DIR}/istisu1.bp.p_ctg.functional.AED.txt"
AED_VALS="${FINAL_DIR}/istisu1.bp.p_ctg.functional.AED.values.tsv"

COURSE="/data/courses/assembly-annotation-course/CDS_annotation"
IMG="${COURSE}/containers/MAKER_3.01.03.sif"
MAKERBIN_HOST="${COURSE}/softwares/Maker_v3.01.03/src/bin"

# AED CDF (bin=0.025)
apptainer exec --bind /data --bind "$COURSE" "$IMG" \
  perl "$MAKERBIN_HOST/AED_cdf_generator.pl" -b 0.025 "$GFF_IN" > "$AED_TXT"

# Flat list of AED per mRNA (handles _AED= or AED=)
awk -F'\t' '$3=="mRNA"{ if (match($9,/(_)?AED=([0-9.]+)/,m)) print m[2] }' "$GFF_IN" > "$AED_VALS"

# Quick % ≤0.5
awk '{n++; if($1<=0.5) k++} END{printf("AED<=0.5: %.1f%% (%d/%d)\n",(n?100*k/n:0),k,n)}' "$AED_VALS"

echo "[OK] Wrote: $AED_TXT and $AED_VALS"
