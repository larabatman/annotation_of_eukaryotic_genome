#!/bin/bash
#SBATCH --job-name=TAIR10_hits
#SBATCH --partition=pibu_el8
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=10
#SBATCH --output=/data/users/%u/annotation_of_eukaryotic_genome/annotation/log/TAIR10_%j.out
#SBATCH --error=/data/users/%u/annotation_of_eukaryotic_genome/annotation/log/TAIR10_%j.err

set -Eeuo pipefail

#----- CONFIG ------
FINAL_DIR="/data/users/${USER}/annotation_of_eukaryotic_genome/annotation/final"
DBDIR="/data/users/${USER}/annotation_of_eukaryotic_genome/annotation/blastdbs"

#------ INPUTS ------
PROT_IN="${FINAL_DIR}/istisu1.bp.p_ctg.proteins.renamed.filtered.fasta"
TAIR_FASTA="${DBDIR}/TAIR10_pep_20110103_representative_gene_model.fa"
TAIR_DB="${DBDIR}/TAIR10_rep"

#----- OUTPUTS ------
BLAST_OUT="${FINAL_DIR}/blastp_vs_TAIR10.outfmt6"
BLAST_BEST="${FINAL_DIR}/blastp_vs_TAIR10.besthits"

#----- MODULES ------
module load BLAST+/2.15.0-gompi-2021a

#CPU variable
CPUS="${SLURM_CPUS_PER_TASK}"

#Sanity check
[ -s "$PROT_IN" ]   || { echo "[ERROR] Missing proteins FASTA: $PROT_IN"; exit 2; }
[ -s "$TAIR_FASTA" ] || { echo "[ERROR] Missing TAIR FASTA (expected $TAIR_FASTA). Build/copy it first."; exit 3; }
#Building our own TAIR10 DB as the one for the course was inaccessible. 
if [ ! -s "${TAIR_DB}.pin" ] || [ ! -s "${TAIR_DB}.psq" ] || [ ! -s "${TAIR_DB}.phr" ]; then
  echo "[INFO] Local TAIR DB missing, building with makeblastdb ..."
  makeblastdb -in "$TAIR_FASTA" -dbtype prot -parse_seqids -out "$TAIR_DB"
fi

#------ RUNNING BLASTP ------
echo "[INFO] BLASTP vs local TAIR10 DB ..."
#protein-protein search keeping up to 10 hits per query
blastp -query "$PROT_IN" -db "$TAIR_DB" -num_threads "$CPUS" -outfmt 6 -evalue 1e-5 -max_target_seqs 10 -out "$BLAST_OUT"

#Fail check
[ -s "$BLAST_OUT" ] || { echo "[ERROR] No TAIR10 hits or BLAST failed: $BLAST_OUT"; exit 4; }

#----- BEST HITS ------
#Order by qseqid then bitscore descending, keeping the first best hit per query
echo "[INFO] Best hit per query (highest bitscore) ..."
LC_ALL=C sort -k1,1 -k12,12gr "$BLAST_OUT" | awk '!seen[$1]++' > "$BLAST_BEST"
[ -s "$BLAST_BEST" ] || { echo "[ERROR] Empty besthits file."; exit 5; }

#------ LOOKING FOR FLC ------
echo "[INFO] Searching for FLC (AT5G10140) among best hits ..."
grep -w 'AT5G10140' "$BLAST_BEST" || echo "FLC (AT5G10140) not found in best hits."

echo "[OK] Done."
echo "Outputs:"
echo "  - $BLAST_OUT"
echo "  - $BLAST_BEST"
