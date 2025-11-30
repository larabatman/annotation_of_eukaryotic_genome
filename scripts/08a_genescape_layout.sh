#!/bin/bash
#SBATCH --job-name=GS_layout
#SBATCH --partition=pibu_el8
#SBATCH --time=00:05:00
#SBATCH --mem=1G
#SBATCH --cpus-per-task=1
#SBATCH --output=/data/users/%u/annotation_of_eukaryotic_genome/annotation/log/GS_layout_%j.out
#SBATCH --error=/data/users/%u/annotation_of_eukaryotic_genome/annotation/log/GS_layout_%j.err

set -Eeuo pipefail

#----- CONFIG -----
ACC="ISTISU1"  # short label for your accession
BASE="/data/users/${USER}/annotation_of_eukaryotic_genome"
FINAL_DIR="${BASE}/annotation/final"

#----- INPUTS -----
ACC_BED_SRC="${FINAL_DIR}/filtered.genes.renamed.gff3"              #we will convert to BED below if needed
ACC_PEPTIDE_SRC="${FINAL_DIR}/istisu1.bp.p_ctg.proteins.renamed.filtered.longest.fasta"  #longest peptide per gene, headers like ISTxxxxx-RA
COURSEDIR="/data/courses/assembly-annotation-course/CDS_annotation/data"
TAIR_FA_SRC="${COURSEDIR}/TAIR10.fa"
TAIR_BED_SRC="${COURSEDIR}/TAIR10.bed"

#----- OUTPUTS -----
GS_ROOT="${BASE}/genespace_workingDirectory"
PEP_DIR="${GS_ROOT}/peptide"
BED_DIR="${GS_ROOT}/bed"

mkdir -p "${PEP_DIR}" "${BED_DIR}"

#----- BUILD BED -----
#Filter rows with feature gene, parse ID=... from column 9 (attributes)
#Then convert GFF 1-based start ($4) to BED 0-based, clamp at 0 and output chr start0 end geneID sorted by chr/ start
ACC_BED_TMP="${BED_DIR}/${ACC}.bed"
awk -F'\t' '
  $3=="gene" {
    n=split($9,a,";"); id=""
    for(i=1;i<=n;i++){ if(a[i] ~ /^ID=/){ split(a[i],b,"="); id=b[2] } }
    if(id!=""){
      s=$4-1; if(s<0) s=0
      print $1, s, $5, id
    }
  }
' OFS='\t' "${ACC_BED_SRC}" | sort -k1,1 -k2,2n > "${ACC_BED_TMP}"
[ -s "${ACC_BED_TMP}" ] || { echo "[ERROR] BED is empty: ${ACC_BED_TMP}"; exit 2; }

#----- BUILD FASTA -----
#Read each FASTA record, take the first token header as the ID and strop isoform suffic to reduce gene ID 
#Concatenate sequence lines, remove whitespace and emit geneID + sequence
#This should result in peptide headers matching BED 
ACC_PEP_TMP="${PEP_DIR}/${ACC}.fa"
awk '
  BEGIN{RS=">"; ORS=""}
  NR>1{
    n=split($1,tok,/[ \t\r\n]/); id=tok[1]; sub(/-R.*/,"",id)
    seq=""
    for(i=2;i<=NF;i++){ gsub(/[ \t\r]/,"",$i); seq=seq $i }
    if(length(seq)>0){ print ">" id "\n" seq "\n" }
  }
' "${ACC_PEPTIDE_SRC}" > "${ACC_PEP_TMP}"
[ -s "${ACC_PEP_TMP}" ] || { echo "[ERROR] Peptide FASTA is empty: ${ACC_PEP_TMP}"; exit 3; }

#Sanity checks 
cut -f4 "${ACC_BED_TMP}" | sort > "${GS_ROOT}/_ids.bed.txt"
grep '^>' "${ACC_PEP_TMP}" | sed 's/^>//' | sort > "${GS_ROOT}/_ids.fa.txt"
M1=$(comm -23 "${GS_ROOT}/_ids.bed.txt" "${GS_ROOT}/_ids.fa.txt" | wc -l | awk '{print $1}')
M2=$(comm -13 "${GS_ROOT}/_ids.bed.txt" "${GS_ROOT}/_ids.fa.txt" | wc -l | awk '{print $1}')
echo "[INFO] ${ACC} BED entries: $(wc -l < "${ACC_BED_TMP}")"
echo "[INFO] ${ACC} peptide entries: $(grep -c '^>' "${ACC_PEP_TMP}")"
echo "[INFO] IDs in BED not in peptide: ${M1}"
echo "[INFO] IDs in peptide not in BED: ${M2}"
