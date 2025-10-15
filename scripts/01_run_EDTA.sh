#!/usr/bin/env bash
#SBATCH --time=2-00:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=20
#SBATCH --job-name=EDTA
#SBATCH --partition=pibu_el8
#SBATCH --output=EDTA_%j.out

# Paths
WORKDIR=/data/users/lland/annotation_of_eukaryotic_genome
GENOME=$WORKDIR/assembly/hifiasm/istisu1.bp.p_ctg.fa
OUTDIR=$WORKDIR/EDTA_annotation
IMG=/data/courses/assembly-annotation-course/CDS_annotation/containers/EDTA2.2.sif

# (Optional) CDS for filtering TE candidates overlapping genes:
CDS=/data/courses/assembly-annotation-course/CDS_annotation/data/TAIR10_cds_20110103_representative_gene_model_updated.fasta

mkdir -p "$OUTDIR"
cd "$OUTDIR"


apptainer exec \
  --bind "$WORKDIR","/data/courses/assembly-annotation-course" \
  "$IMG" \
  EDTA.pl \
    --genome "$GENOME" \                 # Input genome FASTA
    --species others \                   # Parameter preset (generic)
    --step all \                         # Run both library building and annotation
    --sensitive 1 \                      # More sensitive TE discovery (slower)
    --cds "$CDS" \                       # CDS to help filter gene-like sequences
    --anno 1 \                           # Perform genome annotation/masking w/ final TE lib
    --threads "$SLURM_CPUS_PER_TASK"     # Match SLURM threads automatically