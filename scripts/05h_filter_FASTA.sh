#!/bin/bash
#SBATCH --job-name=FilterFASTA
#SBATCH --partition=pibu_el8
#SBATCH --time=00:10:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --output=/data/users/%u/annotation_of_eukaryotic_genome/annotation/log/FilterFASTA_%j.out
#SBATCH --error=/data/users/%u/annotation_of_eukaryotic_genome/annotation/log/FilterFASTA_%j.err

set -Eeuo pipefail

#----- CONTIG -----
FINAL_DIR="/data/users/${USER}/annotation_of_eukaryotic_genome/annotation/final"
GFF="${FINAL_DIR}/filtered.genes.renamed.gff3"
TX_ALL="${FINAL_DIR}/istisu1.bp.p_ctg.all.maker.all.maker.transcripts.fasta.renamed.fasta"
PR_ALL="${FINAL_DIR}/istisu1.bp.p_ctg.all.maker.all.maker.proteins.fasta.renamed.fasta"

LIST="${FINAL_DIR}/list.mRNA.ids.txt"
TX_OUT="${FINAL_DIR}/istisu1.bp.p_ctg.transcripts.renamed.filtered.fasta"
PR_OUT="${FINAL_DIR}/istisu1.bp.p_ctg.proteins.renamed.filtered.fasta"

module load SAMtools/1.13-GCC-10.3.0

[[ -s "$GFF" && -s "$TX_ALL" && -s "$PR_ALL" ]]

#----- HANDLING mRNA IDs ------
awk -F'\t' '$3=="mRNA" && match($9,/(^|;)ID=([^;]+)/,m){print m[2]}' "$GFF" \
  | sort -u > "$LIST"
echo "[INFO] mRNA IDs: $(wc -l < "$LIST")"

#----- INDEXING FASTAs-----
samtools faidx "$TX_ALL"
samtools faidx "$PR_ALL"

#----- SUBSTETTING BY ID LIST ------
samtools faidx "$TX_ALL" -r "$LIST" > "$TX_OUT"
samtools faidx "$PR_ALL" -r "$LIST" > "$PR_OUT"

echo "[INFO] kept transcripts: $(grep -c '^>' "$TX_OUT")"
echo "[INFO] kept proteins  : $(grep -c '^>' "$PR_OUT")"

#report IDs that had no FASTA match (expected to be 0)
comm -23 <(sort "$LIST") <(grep '^>' "$TX_OUT" | sed 's/^>//' | cut -d' ' -f1 | sort) \
  | awk 'NR==1{print "[WARN] IDs missing from FASTA:"} {print "  "$0}' || true

echo "[OK] Done."
