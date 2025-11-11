#!/usr/bin/env bash
#SBATCH --job-name=MakerRenameIDs
#SBATCH --partition=pibu_el8
#SBATCH --time=00:10:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=/data/users/%u/annotation_of_eukaryotic_genome/annotation/log/MakerRename_%j.out
#SBATCH --error=/data/users/%u/annotation_of_eukaryotic_genome/annotation/log/MakerRename_%j.err

set -euo pipefail

#----- CONFIG ----
IMG="/data/courses/assembly-annotation-course/CDS_annotation/containers/MAKER_3.01.03.sif"

#inputs from merge step
NOSEQ_GFF="../annotation/istisu1.bp.p_ctg.all.noseq.gff"
ALL_GFF="../annotation/istisu1.bp.p_ctg.all.gff"
PROT="../annotation/istisu1.bp.p_ctg.all.maker.proteins.fasta"
TRAN="../annotation/istisu1.bp.p_ctg.all.maker.transcripts.fasta"

PREFIX="IST"                          # 3–4 letter accession tag
FINAL="../annotation/final"
mkdir -p "$FINAL"
cd "$FINAL"

#stage copies to modify in place
cp -f "$NOSEQ_GFF" "$(basename "$NOSEQ_GFF").renamed.gff"
cp -f "$ALL_GFF"   "$(basename "$ALL_GFF").renamed.gff"
cp -f "$PROT"      "$(basename "$PROT").renamed.fasta"
cp -f "$TRAN"      "$(basename "$TRAN").renamed.fasta"

GFF_NOSEQ_R="$(basename "$NOSEQ_GFF").renamed.gff"
GFF_ALL_R="$(basename "$ALL_GFF").renamed.gff"
PROT_R="$(basename "$PROT").renamed.fasta"
TRAN_R="$(basename "$TRAN").renamed.fasta"

#----- CREATE ID MAP from the NOSEQ GFF -----
#--justify 7 yields IDs with 7 digits 
apptainer exec --bind /data "$IMG" maker_map_ids --prefix "$PREFIX" --justify 7 "$GFF_NOSEQ_R" > id.map
apptainer exec --bind /data "$IMG" map_gff_ids   id.map "$GFF_NOSEQ_R"
apptainer exec --bind /data "$IMG" map_gff_ids   id.map "$GFF_ALL_R"
apptainer exec --bind /data "$IMG" map_fasta_ids id.map "$PROT_R"
apptainer exec --bind /data "$IMG" map_fasta_ids id.map "$TRAN_R"

echo "[OK] Renamed: $GFF_NOSEQ_R, $GFF_ALL_R, $PROT_R, $TRAN_R"
