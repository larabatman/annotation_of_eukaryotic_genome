#!/bin/bash
#SBATCH --job-name=MakerMerge
#SBATCH --partition=pibu_el8
#SBATCH --time=00:20:00
#SBATCH --time-min=00:10:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --output=/data/users/%u/annotation_of_eukaryotic_genome/annotation/log/MakerMerge_%j.out
#SBATCH --error=/data/users/%u/annotation_of_eukaryotic_genome/annotation/log/MakerMerge_%j.err

set -euo pipefail

#----- CONFIG -----
MAKER_IMG="/data/courses/assembly-annotation-course/CDS_annotation/containers/MAKER_3.01.03.sif"

#Master index telling where each contig's outputs is
DATASTORE_INDEX="/data/users/${USER}/annotation_of_eukaryotic_genome/annotation/maker_run1.maker.output/maker_run1_master_datastore_index.log"
OUTDIR="/data/users/${USER}/annotation_of_eukaryotic_genome/annotation"
#Target filenames
GFF_WITH_FASTA="${OUTDIR}/istisu1.bp.p_ctg.all.gff"
GFF_NO_FASTA="${OUTDIR}/istisu1.bp.p_ctg.all.noseq.gff"
FA_PREFIX="${OUTDIR}/istisu1.bp.p_ctg.all.maker"   # -> .transcripts.fasta / .proteins.fasta

#Sanity check
[ -s "$DATASTORE_INDEX" ] || { echo "ERROR: datastore index not found: $DATASTORE_INDEX"; exit 2; }
mkdir -p "$OUTDIR"

#Merge GFF with embedded FASTA
/usr/bin/apptainer exec --bind /data "$MAKER_IMG" \
  gff3_merge -s -d "$DATASTORE_INDEX" > "$GFF_WITH_FASTA"

#Merge GFF without FASTA
/usr/bin/apptainer exec --bind /data "$MAKER_IMG" \
  gff3_merge -n -s -d "$DATASTORE_INDEX" > "$GFF_NO_FASTA"

#Merge transcripts & proteins FASTA
/usr/bin/apptainer exec --bind /data "$MAKER_IMG" \
  fasta_merge -d "$DATASTORE_INDEX" -o "$FA_PREFIX"

#Report
ls -lh "$GFF_WITH_FASTA" "$GFF_NO_FASTA" "${FA_PREFIX}.transcripts.fasta" "${FA_PREFIX}.proteins.fasta"
