#!/usr/bin/env bash
#SBATCH --job-name=iprscan
#SBATCH --partition=pibu_el8
#SBATCH --time=12:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=2
#SBATCH --output=/data/users/%u/annotation_of_eukaryotic_genome/annotation/log/InterProScan_%j.out
#SBATCH --error=/data/users/%u/annotation_of_eukaryotic_genome/annotation/log/InterProScan_%j.err

set -euo pipefail
#----- CONFIG ------
FINAL_DIR="/data/users/${USER}/annotation_of_eukaryotic_genome/annotation/final" #with renamed GFF and FASTA 
PROT="${FINAL_DIR}/istisu1.bp.p_ctg.all.maker.proteins.fasta.renamed.fasta"
GFF_IN="${FINAL_DIR}/istisu1.bp.p_ctg.all.noseq.gff.renamed.gff"
TSV_OUT="${FINAL_DIR}/output.iprscan.tsv"
GFF_OUT="${FINAL_DIR}/istisu1.bp.p_ctg.all.noseq.gff.renamed.functional.gff"

COURSE="/data/courses/assembly-annotation-course/CDS_annotation"
IPR_IMG="${COURSE}/containers/interproscan_latest.sif"
IPR_DATA="${COURSE}/data/interproscan-5.70-102.0/data"
MAKER_IMG="${COURSE}/containers/MAKER_3.01.03.sif"

TMP="${FINAL_DIR}/ipr_tmp_${SLURM_JOB_ID}"
mkdir -p "$TMP"

#----- RUN InterProScan -----
#Pfam only, with GO and InterProf xrefs
#Run Pfam HMMs for domain calls, avoiding cached DB artifacts for reproducible run 
# --seqtype p: scanning proteins
apptainer exec \
  --bind "${IPR_DATA}:/opt/interproscan/data" \
  --bind "${FINAL_DIR}:/work" \
  --bind "${TMP}:/temp" \
  "${IPR_IMG}" /opt/interproscan/interproscan.sh \
  -appl Pfam --disable-precalc --goterms --iprlookup \
  -f tsv --seqtype p -cpu "${SLURM_CPUS_PER_TASK:-2}" -T /temp \
  -i "/work/$(basename "$PROT")" \
  -o "/work/$(basename "$TSV_OUT")"

#----- INJECT DOMAINS INTO GFF -----
#merge domain and GO annotations back into feature attributes in GFF 
apptainer exec --bind /data "${MAKER_IMG}" \
  ipr_update_gff "${GFF_IN}" "${TSV_OUT}" > "${GFF_OUT}"
