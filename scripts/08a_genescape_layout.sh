#!/bin/bash
#SBATCH --job-name=GS_layout
#SBATCH --partition=pibu_el8
#SBATCH --time=00:05:00
#SBATCH --mem=1G
#SBATCH --cpus-per-task=1
#SBATCH --output=/data/users/%u/annotation_of_eukaryotic_genome/annotation/log/GS_layout_%j.out
#SBATCH --error=/data/users/%u/annotation_of_eukaryotic_genome/annotation/log/GS_layout_%j.err

set -Eeuo pipefail

# ---- Edit these names ----
ACC="ISTISU1"  # short label for your accession
BASE="/data/users/${USER}/annotation_of_eukaryotic_genome"
FINAL_DIR="${BASE}/annotation/final"

# These two were produced earlier:
ACC_BED_SRC="${FINAL_DIR}/filtered.genes.renamed.gff3"              # we will convert to BED below if needed
ACC_PEPTIDE_SRC="${FINAL_DIR}/istisu1.bp.p_ctg.proteins.renamed.filtered.longest.fasta"  # headers like ISTxxxxx-RA

# If you already ran the previous prep and have gene-only names:
#   ${BASE}/genespace_input/${ACC}.bed
#   ${BASE}/genespace_input/${ACC}.fa
# You can point ACC_BED_SRC and ACC_PEPTIDE_SRC directly to them and skip the conversion step.

COURSEDIR="/data/courses/assembly-annotation-course/CDS_annotation/data"
TAIR_FA_SRC="${COURSEDIR}/TAIR10.fa"
TAIR_BED_SRC="${COURSEDIR}/TAIR10.bed"

# Destination (GENESPACE workingDirectory)
GS_ROOT="${BASE}/genespace_workingDirectory"
PEP_DIR="${GS_ROOT}/peptide"
BED_DIR="${GS_ROOT}/bed"

mkdir -p "${PEP_DIR}" "${BED_DIR}"

# ---- Build BED from your filtered GFF (gene features only; 0-based start) ----
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

# ---- Build peptide FASTA with headers == gene IDs (strip -R.*) ----
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

# ---- Minimal consistency check ----
cut -f4 "${ACC_BED_TMP}" | sort > "${GS_ROOT}/_ids.bed.txt"
grep '^>' "${ACC_PEP_TMP}" | sed 's/^>//' | sort > "${GS_ROOT}/_ids.fa.txt"
M1=$(comm -23 "${GS_ROOT}/_ids.bed.txt" "${GS_ROOT}/_ids.fa.txt" | wc -l | awk '{print $1}')
M2=$(comm -13 "${GS_ROOT}/_ids.bed.txt" "${GS_ROOT}/_ids.fa.txt" | wc -l | awk '{print $1}')
echo "[INFO] ${ACC} BED entries: $(wc -l < "${ACC_BED_TMP}")"
echo "[INFO] ${ACC} peptide entries: $(grep -c '^>' "${ACC_PEP_TMP}")"
echo "[INFO] IDs in BED not in peptide: ${M1}"
echo "[INFO] IDs in peptide not in BED: ${M2}"
