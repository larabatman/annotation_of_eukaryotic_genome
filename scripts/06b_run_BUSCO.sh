#!/bin/bash
#SBATCH --job-name=BUSCO
#SBATCH --partition=pibu_el8
#SBATCH --time=06:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --output=/data/users/%u/annotation_of_eukaryotic_genome/annotation/log/BUSCO_%j.out
#SBATCH --error=/data/users/%u/annotation_of_eukaryotic_genome/annotation/log/BUSCO_%j.err

set -Eeuo pipefail

#----- CONFIG -----
FINAL_DIR="/data/users/${USER}/annotation_of_eukaryotic_genome/annotation/final"
PROT_LONGEST="${FINAL_DIR}/istisu1.bp.p_ctg.proteins.renamed.filtered.longest.fasta"
TRAN_LONGEST="${FINAL_DIR}/istisu1.bp.p_ctg.transcripts.renamed.filtered.longest.fasta"

LINEAGE="brassicales_odb10"

module load BUSCO/5.4.2-foss-2021a

#Sanity checks
[ -s "$PROT_LONGEST" ] || { echo "[ERROR] Missing longest proteins: $PROT_LONGEST"; exit 2; }
[ -s "$TRAN_LONGEST" ] || { echo "[ERROR] Missing longest transcripts: $TRAN_LONGEST"; exit 3; }

#----- RUN BUSCO -----
#-i: protein inputs, one isoform per gene
#-l: lineage dataset, here brassicales_odn10
#-o: output folders
#-m proteins: protein mode
echo "[INFO] BUSCO proteins mode on $PROT_LONGEST"
busco -i "$PROT_LONGEST" -l "$LINEAGE" -o busco_proteins -m proteins -c "${SLURM_CPUS_PER_TASK}"
#-m transcriptome for transcript mode
echo "[INFO] BUSCO transcriptome mode on $TRAN_LONGEST"
busco -i "$TRAN_LONGEST" -l "$LINEAGE" -o busco_transcripts -m transcriptome -c "${SLURM_CPUS_PER_TASK}"

echo "[OK] BUSCO runs launched."
echo "  - busco_proteins/short_summary*.txt"
echo "  - busco_transcripts/short_summary*.txt"
