#!/bin/bash
#SBATCH --job-name=trinity2genome
#SBATCH --partition=pibu_el8
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --output=/data/users/%u/annotation_of_eukaryotic_genome/assembly/Trinity/log/trinity2genome_%j.out
#SBATCH --error=/data/users/%u/annotation_of_eukaryotic_genome/assembly/Trinity/log/trinity2genome_%j.err

set -Eeuo pipefail

# --- paths ---
GENOME="/data/users/${USER}/annotation_of_eukaryotic_genome/assembly/hifiasm/istisu1.bp.p_ctg.fa"
TRINITY="/data/users/${USER}/annotation_of_eukaryotic_genome/assembly/Trinity/trinity.Trinity.fasta"   # adjust if named differently
OUTDIR="/data/users/${USER}/annotation_of_eukaryotic_genome/assembly/Trinity/aln_to_genome"
IDX="$OUTDIR/genome.mmi"
BAM="$OUTDIR/trinity_vs_genome.spliced.sorted.bam"

mkdir -p "$OUTDIR" "/data/users/${USER}/annotation_of_eukaryotic_genome/assembly/Trinity/log"

# --- modules ---
module load minimap2/2.20-GCCcore-10.3.0
module load SAMtools/1.13-GCC-10.3.0 

# --- sanity ---
[ -s "$GENOME" ]  || { echo "[ERROR] Missing genome: $GENOME"; exit 2; }
[ -s "$TRINITY" ] || { echo "[ERROR] Missing Trinity transcripts: $TRINITY"; exit 3; }

# --- index genome for minimap2 (once) ---
if [ ! -s "$IDX" ]; then
  echo "[INFO] Building minimap2 index -> $IDX"
  minimap2 -d "$IDX" "$GENOME"
else
  echo "[INFO] Reusing index: $IDX"
fi

# --- splice-aware alignment (cDNA/transcripts) -> sorted BAM ---
T="${SLURM_CPUS_PER_TASK:-8}"
echo "[INFO] Aligning Trinity to genome with $T threads"
minimap2 -t "$T" -ax splice -uf -k14 --secondary=no "$IDX" "$TRINITY" \
  | samtools view -@ "$T" -b - \
  | samtools sort -@ "$T" -o "$BAM"

samtools index "$BAM"

echo "[OK] BAM ready:"
echo "  $BAM"
echo "  ${BAM}.bai"

# Next in Geneious: import the genome FASTA, group into Sequence List, import this BAM,
# then drag the BAM track above the gene track in the side panel.
