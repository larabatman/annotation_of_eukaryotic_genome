#!/bin/bash
#SBATCH --job-name=AGAT_stats
#SBATCH --partition=pibu_el8
#SBATCH --time=00:30:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=/data/users/%u/annotation_of_eukaryotic_genome/annotation/log/AGAT_%j.out
#SBATCH --error=/data/users/%u/annotation_of_eukaryotic_genome/annotation/log/AGAT_%j.err

set -Eeuo pipefail

FINAL_DIR="/data/users/${USER}/annotation_of_eukaryotic_genome/annotation/final"
GFF_IN="${FINAL_DIR}/filtered.genes.renamed.gff3"

OUT_DIR="/data/users/${USER}/annotation_of_eukaryotic_genome/annotation/agat_stats"
OUT_TXT="${OUT_DIR}/annotation_statistics.txt"

APPTAINER_IMG="/data/courses/assembly-annotation-course/CDS_annotation/containers/agat_1.5.1--pl5321hdfd78af_0.sif"

mkdir -p "$OUT_DIR"

[ -s "$GFF_IN" ] || { echo "[ERROR] Missing GFF: $GFF_IN"; exit 2; }
[ -s "$APPTAINER_IMG" ] || { echo "[ERROR] Missing AGAT container: $APPTAINER_IMG"; exit 3; }

echo "[INFO] Running AGAT statistics -> $OUT_TXT"
/usr/bin/apptainer exec --bind /data "$APPTAINER_IMG" agat_sp_statistics.pl -i "$GFF_IN" -o "$OUT_TXT"

[ -s "$OUT_TXT" ] || { echo "[ERROR] Empty AGAT output: $OUT_TXT"; exit 4; }

echo "[OK] Done. See $OUT_TXT"
echo "[INFO] Quick peek:"
grep -E 'number_of_genes|number_of_transcripts|number_of_CDS' -n "$OUT_TXT" || true
