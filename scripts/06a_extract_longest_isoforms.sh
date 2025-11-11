#!/bin/bash
#SBATCH --job-name=longest_iso
#SBATCH --partition=pibu_el8
#SBATCH --time=00:20:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=/data/users/%u/annotation_of_eukaryotic_genome/annotation/log/longest_%j.out
#SBATCH --error=/data/users/%u/annotation_of_eukaryotic_genome/annotation/log/longest_%j.err

set -Eeuo pipefail

#----- CONFIG -----
cd "/data/users/${USER}/annotation_of_eukaryotic_genome/annotation/final"
#Input filtered FASTA file from previous step
IN_PROT="istisu1.bp.p_ctg.proteins.renamed.filtered.fasta"
IN_TRAN="istisu1.bp.p_ctg.transcripts.renamed.filtered.fasta"

#Output longest FASTA sequence: avoid isomers, one sequence per transcript
OUT_PROT="istisu1.bp.p_ctg.proteins.renamed.filtered.longest.fasta"
OUT_TRAN="istisu1.bp.p_ctg.transcripts.renamed.filtered.longest.fasta"

#Load module
module load SAMtools/1.13-GCC-10.3.0

#Sanity check
[ -s "$IN_PROT" ] || { echo "Missing $IN_PROT"; exit 2; }
[ -s "$IN_TRAN" ] || { echo "Missing $IN_TRAN"; exit 3; }

#Helper function using samtools to index FASTA and extract IDs of longest isoforms
pick_longest() {
  local in="$1" out="$2"
  samtools faidx "$in"
  samtools faidx "$in" -r <(
    awk '{
      # .fai: $1=name, $2=length
      if (FILENAME ~ /\.fai$/) {
        g=$1; sub(/-R.*/, "", g)            # base gene id (drop -R*)
        if ($2 > max[g]) { max[g]=$2; best[g]=$1 }
      }
    } END { for (g in best) print best[g] }' "${in}.fai"
  ) > "$out"
}

pick_longest "$IN_PROT" "$OUT_PROT"
pick_longest "$IN_TRAN" "$OUT_TRAN"

echo "proteins:   $(grep -c '^>' "$IN_PROT")  → longest: $(grep -c '^>' "$OUT_PROT")"
echo "transcripts: $(grep -c '^>' "$IN_TRAN") → longest: $(grep -c '^>' "$OUT_TRAN")"
