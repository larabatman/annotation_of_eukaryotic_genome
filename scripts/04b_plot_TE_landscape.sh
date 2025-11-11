#!/usr/bin/env bash
#SBATCH --job-name=plot_div
#SBATCH --partition=pibu_el8
#SBATCH --time=00:25:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --output=%x_%j.out

set -euo pipefail

#----- COBNFIG -----
WORKDIR=/data/users/lland/annotation_of_eukaryotic_genome
PARSE_DIR=$WORKDIR/TE_age/parseRM          # parseRM outputs live here
OUTPLOT=$WORKDIR/TE_age/landscape          # plots go here
PLOT_R=$WORKDIR/scripts/04b_plot_div.R

mkdir -p "$OUTPLOT"

#Sanity check for input files
[ -d "$PARSE_DIR" ] || { echo "ERROR: missing $PARSE_DIR"; exit 1; }
[ -s "$PLOT_R" ]    || { echo "ERROR: missing plot script $PLOT_R"; exit 2; }

module purge
module load R/4.2.1-foss-2021a

#Course script
Rscript "$PLOT_R" "$PARSE_DIR" "$OUTPLOT"

echo "[OK] Landscape written to: $OUTPLOT"
ls -lh "$OUTPLOT" | sed -n '1,100p'
