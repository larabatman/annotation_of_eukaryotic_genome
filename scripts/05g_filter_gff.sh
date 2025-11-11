#!/bin/bash
#SBATCH --job-name=Maker_FilterGFF
#SBATCH --partition=pibu_el8
#SBATCH --time=00:20:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --output=/data/users/%u/annotation_of_eukaryotic_genome/annotation/log/FilterGFF_%j.out
#SBATCH --error=/data/users/%u/annotation_of_eukaryotic_genome/annotation/log/FilterGFF_%j.err

set -Eeuo pipefail

#----- CONFIG -----
FINAL_DIR="/data/users/${USER}/annotation_of_eukaryotic_genome/annotation/final"
GFF_IN="${FINAL_DIR}/istisu1.bp.p_ctg.all.noseq.gff.renamed.functional.gff" #Input to filter

COURSEDIR="/data/courses/assembly-annotation-course/CDS_annotation"
IMG="${COURSEDIR}/containers/MAKER_3.01.03.sif"
MAKERBIN_HOST="${COURSEDIR}/softwares/Maker_v3.01.03/src/bin"

#----- OUTPUTS -----
GFF_QF="${FINAL_DIR}/istisu1.bp.p_ctg.all.noseq.gff.renamed.functional_iprscan_quality_filtered.gff"
GFF_GENES="${FINAL_DIR}/filtered.genes.renamed.gff3"

#Sanity checks
[[ -s "$GFF_IN" ]] || { echo "[ERROR] Missing input GFF: $GFF_IN"; exit 2; }

#----- QUALITY FILTER -----
echo "[INFO] Running quality_filter.pl (-s: AED<1 and/or Pfam present) ..."

# Runs the MAKER quality filter inside the container; writes to GFF_QF
/usr/bin/apptainer exec --bind /data --bind "${COURSEDIR}" "${IMG}" perl "${MAKERBIN_HOST}/quality_filter.pl" -s "${GFF_IN}" > "${GFF_QF}"
# -s keeps transcripts with AED < 1 and/or Pfam tag if present in the GFF 
[[ -s "$GFF_QF" ]] || { echo "[ERROR] quality_filter.pl produced empty output: $GFF_QF"; exit 3; }

echo "[INFO] Keeping only gene-structure features (gene,mRNA,exon,CDS,UTRs) ..."
#Keeping the FEATURE subset 
FEATURES=$'\tgene\t|\tmRNA\t|\texon\t|\tCDS\t|\tfive_prime_UTR\t|\tthree_prime_UTR\t'
grep -P "$FEATURES" "${GFF_QF}" > "${GFF_GENES}"

[[ -s "$GFF_GENES" ]] || { echo "[ERROR] Feature-filtered GFF is empty: $GFF_GENES"; exit 4; }

echo "[INFO] Feature types present in ${GFF_GENES}:"
cut -f3 "${GFF_GENES}" | sort | uniq -c | sort -nr

echo "[INFO] Quick counts:"
echo -n "  genes : ";  grep -P $'\tgene\t'  "${GFF_GENES}" | wc -l
echo -n "  mRNA  : ";  grep -P $'\tmRNA\t'  "${GFF_GENES}" | wc -l
echo -n "  CDS   : ";  grep -P $'\tCDS\t'   "${GFF_GENES}" | wc -l

echo "[OK] Done."
echo "Outputs:"
echo "  - ${GFF_QF}"
echo "  - ${GFF_GENES}"
