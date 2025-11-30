#!/usr/bin/env bash
#SBATCH --job-name=EDTA
#SBATCH --partition=pibu_el8
#SBATCH --time=3-00:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem=128G
#SBATCH --output=%x_%j.out

# Abort on error -e, undefined variable -u and pipe failures with pipefile
set -euo pipefail

#----- INPUTS -----
WORKDIR=/data/users/lland/annotation_of_eukaryotic_genome # project root
GENOME=$WORKDIR/assembly/hifiasm/istisu1.bp.p_ctg.fa # assembly FASTA to annotate
IMG=/data/courses/assembly-annotation-course/CDS_annotation/containers/EDTA2.2.sif # Apptainer image with EDTA
CDS_LOCAL=$WORKDIR/ref/TAIR10_cds.fa # coding sequences to improve classification

#----- WORKING DIR -----
RUN_DIR="$WORKDIR/EDTA_run_annotation"
mkdir -p "$RUN_DIR"
cd "$RUN_DIR"

#----- CHECKS -----
# sanity checks: genome and containers must exist
[ -s "$GENOME" ]   || { echo "ERROR: genome not found: $GENOME" >&2; exit 1; }
[ -s "$IMG" ]      || { echo "ERROR: EDTA image not found: $IMG" >&2; exit 1; }

# Making an array of arguments to avoid backslashing 
args=(
  --genome "$GENOME"  # input assembly
  --species others  # generic model 
  --step all  # run full EDTA pipeline with structure discoveryl classification, masking, annotation
  --sensitive 1 # include RepeatModeler for better TE discovery
  --anno 1  # produce final TE annotation with GFF and masked sequences files 
  --threads "${SLURM_CPUS_PER_TASK}"  # use SLURM CPUs
  --cds "$CDS_LOCAL"
)

#------ RUN EDTA -----
echo "[INFO] Running: EDTA.pl ${args[*]}" # exact command log
# Run in container, bind project and course data
apptainer exec --bind "$WORKDIR","/data/courses/assembly-annotation-course" \
  "$IMG" EDTA.pl "${args[@]}" # call EDTA with the array of arguments 

#----- WRAP-UP ------
echo "[OK] EDTA completed."
echo "[INFO] Key outputs in: $RUN_DIR"
# Peak at the script 
ls -lh "$RUN_DIR" | sed -n '1,200p'