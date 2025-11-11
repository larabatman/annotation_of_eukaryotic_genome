#!/bin/bash
#SBATCH --job-name=UniProt_QC
#SBATCH --partition=pibu_el8
#SBATCH --time=00:15:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=/data/users/%u/annotation_of_eukaryotic_genome/annotation/log/UniProt_QC_%j.out
#SBATCH --error=/data/users/%u/annotation_of_eukaryotic_genome/annotation/log/UniProt_QC_%j.err

set -Eeuo pipefail

FINAL_DIR="/data/users/${USER}/annotation_of_eukaryotic_genome/annotation/final"
PROT_IN="${FINAL_DIR}/istisu1.bp.p_ctg.proteins.renamed.filtered.fasta"
BLAST_BEST="${FINAL_DIR}/blastp_vs_uniprot.besthits"

[ -s "$PROT_IN" ] || { echo "[ERROR] Missing proteins FASTA: $PROT_IN"; exit 2; }

echo "[INFO] Counting total protein entries ..."
N_ALL=$(grep -c '^>' "$PROT_IN" || true)
echo "[INFO] Total proteins: $N_ALL"

if [ -s "$BLAST_BEST" ]; then
  N_HIT=$(cut -f1 "$BLAST_BEST" | sort -u | wc -l | awk '{print $1}')
  echo "[INFO] Proteins with UniProt best hit: $N_HIT ($(( 100 * N_HIT / (N_ALL>0?N_ALL:1) ))%)"
else
  echo "[WARN] No UniProt besthits file; cannot compute % with hits."
fi

echo "[INFO] Crude length bias check (aa lengths) ..."
# Build an ID->length table for proteins
awk '
  BEGIN{RS=">"; ORS=""; FS="\n"}
  NR>1{
    header=$1; split(header,a," "); id=a[1];
    seq="";
    for(i=2;i<=NF;i++){ gsub(/[ \t\r]/,"",$i); seq=seq $i }
    print id "\t" length(seq) "\n"
  }
' "$PROT_IN" > "${FINAL_DIR}/prot.lengths.tsv"

if [ -s "$BLAST_BEST" ]; then
  cut -f1 "$BLAST_BEST" | sort -u > "${FINAL_DIR}/ids.hit.txt"
  awk 'NR==FNR{hit[$1]=1; next} {print $2 > ( ($1 in hit) ? "'"${FINAL_DIR}"'/len.hit.txt" : "'"${FINAL_DIR}"'/len.nohit.txt" )}' "${FINAL_DIR}/ids.hit.txt" "${FINAL_DIR}/prot.lengths.tsv"
  echo -n "[INFO] median length (hit):   "
  sort -n "${FINAL_DIR}/len.hit.txt"   | awk ' {a[NR]=$1} END{ if(NR%2){print a[(NR+1)/2]} else {print (a[NR/2]+a[NR/2+1])/2} }'
  echo -n "[INFO] median length (nohit): "
  sort -n "${FINAL_DIR}/len.nohit.txt" | awk ' {a[NR]=$1} END{ if(NR==0){print "NA"} else if(NR%2){print a[(NR+1)/2]} else {print (a[NR/2]+a[NR/2+1])/2} }'
else
  echo "[WARN] Skipping length bias check (no besthits)."
fi

echo "[OK] QC done."
