#!/usr/bin/env bash
#SBATCH --job-name=maker_ctl_setup
#SBATCH --partition=pibu_el8
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=%x_%j.out

set -euo pipefail

#----- CONFIG -----
WORKDIR=/data/users/lland/annotation_of_eukaryotic_genome
ANNODIR=$WORKDIR/annotation
GENOME=$WORKDIR/assembly/hifiasm/istisu1.bp.p_ctg.fa
IMG=/data/courses/assembly-annotation-course/CDS_annotation/containers/MAKER_3.01.03.sif

#EDTA outputs
EDTA_DIR=$WORKDIR/EDTA_annotation
PREFIX=$(basename "$GENOME")
TELIB=$EDTA_DIR/${PREFIX}.mod.EDTA.TElib.fa

#Evidence
COURSE=/data/courses/assembly-annotation-course/CDS_annotation/data
#TAIR peptides
TAIR_PEP=$COURSE/TAIR10_pep_20110103_representative_gene_model
#Reviewed UniProt Viridiplantae
UNIPROT_PEP=$COURSE/uniprot_viridiplantae_reviewed.fa
#PTERP20
PTREP=$COURSE/PTREP20
#Trinity FASTA
TRINITY=$WORKDIR/assembly/Trinity/trinity.Trinity.fasta

mkdir -p "$ANNODIR"
cd "$ANNODIR"

#----- GENERATE CONTROL FILES -----
apptainer exec --bind "$WORKDIR","/data/courses/assembly-annotation-course" "$IMG" maker -CTL
[ -s maker_opts.ctl ] || { echo "maker_opts.ctl missing"; exit 1; }

#----- SCRATCH -----
SCR="${SCRATCH:-$ANNODIR/tmp}"
mkdir -p "$SCR"

#Escape '&' for sed
esc() { printf '%s' "$1" | sed 's/[&]/\\&/g'; }

#----- PATCH maker_opts.ctl -----
#disable DFam, use EDTA TElib, evidence paths, cpus=1
sed -i \
  -e "s|^genome=.*|genome=$(esc "$GENOME")|" \
  -e "s|^est=.*|est=$(esc "$TRINITY")|" \
  -e "s|^protein=.*|protein=$(esc "$TAIR_PEP"),$(esc "$UNIPROT_PEP")|" \
  -e "s|^model_org=.*|model_org=|" \
  -e "s|^rmlib=.*|rmlib=$(esc "$TELIB")|" \
  -e "s|^repeat_protein=.*|repeat_protein=$(esc "$PTREP")|" \
  -e "s|^augustus_species=.*|augustus_species=arabidopsis|" \
  -e "s|^est2genome=.*|est2genome=1|" \
  -e "s|^protein2genome=.*|protein2genome=1|" \
  -e "s|^cpus=.*|cpus=1|" \
  -e "s|^alt_splice=.*|alt_splice=1|" \
  -e "s|^TMP=.*|TMP=$(esc "$SCR")|" \
  maker_opts.ctl

echo "==== maker_opts.ctl key fields ===="
grep -E '^(genome|est=|protein=|model_org=|rmlib=|repeat_protein=|augustus_species=|est2genome=|protein2genome=|cpus=|alt_splice=|TMP=)' maker_opts.ctl
echo "[OK] Control files ready in: $ANNODIR"
