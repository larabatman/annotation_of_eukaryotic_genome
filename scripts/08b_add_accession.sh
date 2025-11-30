#!/bin/bash
# Usage:
#   bash 10a_add_accession.sh ACC_LABEL /path/to/filtered.gff3 /path/to/longest_proteins.fasta

set -Eeuo pipefail

#----- VARIABLE INPUT -----
ACC="$1" #Accession label
GFF="$2" #GFF3
PFA="$3" #

BASE="/data/users/${USER}/annotation_of_eukaryotic_genome"
GS_ROOT="${BASE}/genespace_workingDirectory"
PEP_DIR="${GS_ROOT}/peptide"
BED_DIR="${GS_ROOT}/bed"

[ -n "${ACC:-}" ] || { echo "Usage: $0 ACC_LABEL filtered.gff3 longest_proteins.fasta"; exit 1; }
[ -s "$GFF" ] || { echo "[ERROR] Missing GFF: $GFF"; exit 2; }
[ -s "$PFA" ] || { echo "[ERROR] Missing proteins FASTA: $PFA"; exit 3; }

mkdir -p "${PEP_DIR}" "${BED_DIR}"

OUT_BED="${BED_DIR}/${ACC}.bed"
OUT_FA="${PEP_DIR}/${ACC}.fa"

awk -F'\t' '$3=="gene"{split($9,a,";");id=""; for(i=1;i<=length(a);i++){ if(a[i] ~ /^ID=/){split(a[i],b,"="); id=b[2]} } if(id!=""){s=$4-1; if(s<0)s=0; print $1, s, $5, id}}' OFS='\t' "$GFF" \
| sort -k1,1 -k2,2n > "$OUT_BED"

awk 'BEGIN{RS=">"; ORS=""} NR>1{ n=split($1,t,/[ \t\r\n]/); id=t[1]; sub(/-R.*/,"",id); seq=""; for(i=2;i<=NF;i++){ gsub(/[ \t\r]/,"",$i); seq=seq $i } if(length(seq)>0){ print ">" id "\n" seq "\n" }}' "$PFA" > "$OUT_FA"

[ -s "$OUT_BED" ] || { echo "[ERROR] Empty BED for $ACC"; exit 4; }
[ -s "$OUT_FA" ]  || { echo "[ERROR] Empty peptide FASTA for $ACC"; exit 5; }

echo "[OK] Added $ACC:"
echo "  - $OUT_BED"
echo "  - $OUT_FA"
