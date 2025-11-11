#!/usr/bin/env bash
#SBATCH --job-name=parseRM
#SBATCH --partition=pibu_el8
#SBATCH --time=00:15:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=%x_%j.out

set -euo pipefail

#----- CONFIG ------
WORKDIR=/data/users/lland/annotation_of_eukaryotic_genome
GENOME=$WORKDIR/assembly/hifiasm/istisu1.bp.p_ctg.fa
EDTA=$WORKDIR/EDTA_annotation
OUT=$WORKDIR/TE_age/parseRM
mkdir -p "$OUT"

PREFIX=$(basename "$GENOME")
MODOUT="$EDTA/${PREFIX}.mod.EDTA.anno/${PREFIX}.mod.out" #RepeatMasker .out protuced by EDTA
#Sanity check for the .mod.out file
[ -s "$MODOUT" ] || { echo "ERROR: not found: $MODOUT"; exit 1; }

#Script given by the course
PARSERM=/data/courses/assembly-annotation-course/CDS_annotation/scripts/05-parseRM.pl
[ -s "$PARSERM" ] || { echo "ERROR: parseRM.pl missing at $PARSERM"; exit 2; }

module purge
#Load BioPerl to run the .pl script
module load BioPerl/1.7.8-GCCcore-10.3.0 || true

cd "$OUT"

#Symlink the .out into OUT and run parseRM on the local symlink
ln -sf "$MODOUT" "$(basename "$MODOUT")"

#Run with binning (-l 50,1): bin size 1% across 0-50% divergence 
#Call parseRM.pl on the symlinked .out ans save stdout and stderr separately 
perl "$PARSERM" -i "$(basename "$MODOUT")" -l 50,1 -v \
  > parseRM_stdout.log 2> parseRM_stderr.log

echo "[INFO] Files in $OUT after parseRM:"
ls -lh "$OUT" | sed -n '1,200p'
