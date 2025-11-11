#!/usr/bin/env bash
#SBATCH --job-name=extract_CG
#SBATCH --partition=pibu_el8
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=%x_%j.out

set -euo pipefail

#----- CONFIG -----
WORKDIR=/data/users/lland/annotation_of_eukaryotic_genome
GENOME=$WORKDIR/assembly/hifiasm/istisu1.bp.p_ctg.fa
EDTA_DIR=$WORKDIR/EDTA_annotation
TELIB=$EDTA_DIR/$(basename "$GENOME").mod.EDTA.TElib.fa
OUT=$WORKDIR/TEsorter_refine
mkdir -p "$OUT"

#Sanity check for the EDTA TE library
[ -s "$TELIB" ] || { echo "ERROR: missing EDTA library: $TELIB"; exit 1; }

#----- OUTPUT PATHS -----
COP="$OUT/Copia_sequences.fa"
GYP="$OUT/Gypsy_sequences.fa"

#----- SUBSETTING WITH SEQKIT ------
# Prefer seqkit if available; otherwise use an awk fallback
if command -v seqkit >/dev/null 2>&1; then
  seqkit grep -r -p "Copia" "$TELIB" > "$COP"
  seqkit grep -r -p "Gypsy" "$TELIB" > "$GYP"
else
  # Fallback: stream FASTA records and keep those whose header contains the token
  awk -v tok="Copia" 'BEGIN{RS=">"; ORS=""} NR>1{hdr=$0; n=substr(hdr,1,index(hdr,"\n")-1); if (hdr ~ tok) print ">" hdr}' "$TELIB" > "$COP"
  awk -v tok="Gypsy" 'BEGIN{RS=">"; ORS=""} NR>1{hdr=$0; n=substr(hdr,1,index(hdr,"\n")-1); if (hdr ~ tok) print ">" hdr}' "$TELIB" > "$GYP"
fi
# Report count headers
echo "[OK] Wrote:"
echo "  $COP  ($(grep -c '^>' "$COP") seqs)"
echo "  $GYP  ($(grep -c '^>' "$GYP") seqs)"
