#!/bin/bash
#SBATCH --job-name=UniProt_annot
#SBATCH --partition=pibu_el8
#SBATCH --time=04:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=10
#SBATCH --output=/data/users/%u/annotation_of_eukaryotic_genome/annotation/log/UniProt_%j.out
#SBATCH --error=/data/users/%u/annotation_of_eukaryotic_genome/annotation/log/UniProt_%j.err

set -Eeuo pipefail

# -------------------- Directories --------------------
FINAL_DIR="/data/users/${USER}/annotation_of_eukaryotic_genome/annotation/final"
COURSEDIR="/data/courses/assembly-annotation-course/CDS_annotation"
MAKERBIN="${COURSEDIR}/softwares/Maker_v3.01.03/src/bin"

# -------------------- Inputs --------------------
PROT_IN="${FINAL_DIR}/istisu1.bp.p_ctg.proteins.renamed.filtered.fasta"
GFF_IN="${FINAL_DIR}/filtered.genes.renamed.gff3"

# Curated UniProt (Viridiplantae reviewed) FASTA; BLAST DB already built
UNIPROT_FA="${COURSEDIR}/data/uniprot/uniprot_viridiplantae_reviewed.fa"
DB="${UNIPROT_FA}"

# -------------------- Outputs --------------------
BLAST_OUT="${FINAL_DIR}/blastp_vs_uniprot.outfmt6"
BLAST_BEST="${FINAL_DIR}/blastp_vs_uniprot.besthits"
FASTA_UNI="${FINAL_DIR}/istisu1.bp.p_ctg.proteins.renamed.filtered.fasta.UniProt"
GFF_UNI="${FINAL_DIR}/filtered.genes.renamed.gff3.UniProt.gff3"

# -------------------- Tools --------------------
module load BLAST+/2.15.0-gompi-2021a

# -------------------- Checks --------------------
[ -s "$PROT_IN" ]     || { echo "[ERROR] Missing proteins FASTA: $PROT_IN"; exit 2; }
[ -s "$GFF_IN" ]      || { echo "[ERROR] Missing filtered GFF: $GFF_IN"; exit 3; }
[ -s "$UNIPROT_FA" ]  || { echo "[ERROR] Missing UniProt FASTA: $UNIPROT_FA"; exit 4; }
# Ensure CPU count even if SLURM var is absent
CPUS="${SLURM_CPUS_PER_TASK:-10}"

# -------------------- BLASTP vs UniProt --------------------
echo "[INFO] Running BLASTP vs UniProt reviewed (Viridiplantae) ..."
# -query: your proteins
# -db: UniProt reviewed Viridiplantae
# -num_threads: CPU threads
# -outfmt 6: tabular (12 cols; 11=evalue, 12=bitscore)
# -evalue 1e-5: significance cutoff
# -max_target_seqs 10: keep up to 10 hits per query
blastp -query "$PROT_IN" -db "$DB" -num_threads "$CPUS" -outfmt 6 -evalue 1e-5 -max_target_seqs 10 -out "$BLAST_OUT"

[ -s "$BLAST_OUT" ] || { echo "[ERROR] BLAST produced no hits or failed: $BLAST_OUT"; exit 5; }

# -------------------- Best hit per query --------------------
echo "[INFO] Keeping best hit per query (highest bitscore) ..."
# Sort by qseqid asc, bitscore desc; keep first per qseqid
LC_ALL=C sort -k1,1 -k12,12gr "$BLAST_OUT" | awk '!seen[$1]++' > "$BLAST_BEST"
[ -s "$BLAST_BEST" ] || { echo "[ERROR] No best hits written (empty $BLAST_BEST)"; exit 6; }

# -------------------- Inject UniProt into FASTA --------------------
echo "[INFO] Writing UniProt names into FASTA headers ..."
# maker_functional_fasta <uniprot.fasta> <blast.besthits> <maker_proteins.fasta>
"$MAKERBIN/maker_functional_fasta" "$UNIPROT_FA" "$BLAST_BEST" "$PROT_IN" > "$FASTA_UNI"
[ -s "$FASTA_UNI" ] || { echo "[ERROR] maker_functional_fasta produced empty output: $FASTA_UNI"; exit 7; }

# -------------------- Inject UniProt into GFF --------------------
echo "[INFO] Writing UniProt names into GFF attributes ..."
# maker_functional_gff <uniprot.fasta> <blast.outfmt6> <maker.gff3>
"$MAKERBIN/maker_functional_gff" "$UNIPROT_FA" "$BLAST_OUT" "$GFF_IN" > "$GFF_UNI"
[ -s "$GFF_UNI" ] || { echo "[ERROR] maker_functional_gff produced empty output: $GFF_UNI"; exit 8; }

# -------------------- Summary --------------------
echo "[OK] Done."
echo "Outputs:"
echo "  - $BLAST_OUT"
echo "  - $BLAST_BEST"
echo "  - $FASTA_UNI"
echo "  - $GFF_UNI"
