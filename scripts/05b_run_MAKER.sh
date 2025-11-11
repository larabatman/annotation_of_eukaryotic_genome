#!/usr/bin/env bash
#SBATCH --job-name=maker_mpi
#SBATCH --partition=pibu_el8
#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=50
#SBATCH --mem=120G
#SBATCH --output=/data/users/%u/annotation_of_eukaryotic_genome/annotation/log/maker_%j.out
#SBATCH --error=/data/users/%u/annotation_of_eukaryotic_genome/annotation/log/maker_%j.err

set -euo pipefail

#----- CONFIG -----
WORKDIR=/data/users/${USER}/annotation_of_eukaryotic_genome
ANNODIR=$WORKDIR/annotation
COURSEDIR=/data/courses/assembly-annotation-course/CDS_annotation
IMG=$COURSEDIR/containers/MAKER_3.01.03.sif

module load OpenMPI/4.1.1-GCC-10.3.0

cd "$ANNODIR"

#tmp in the container
SCR="${SCRATCH:-$ANNODIR/tmp}"
mkdir -p "$SCR"

#Sanity check for the control files
for f in maker_opts.ctl maker_bopts.ctl maker_exe.ctl; do
  [ -s "$f" ] || { echo "Missing $f"; exit 1; }
done

NTASKS=${SLURM_NTASKS:-${SLURM_NTASKS_PER_NODE:-50}}
echo "[INFO] $(date)  tasks=$NTASKS  wd=$ANNODIR"

#binding necessary files
BIND=(--bind /data --bind "$SCR":/TMP --bind "$COURSEDIR")
# expose host AUGUSTUS config only if user set it
if [ -n "${AUGUSTUS_CONFIG_PATH:-}" ]; then
  BIND+=(--bind "$AUGUSTUS_CONFIG_PATH" --env AUGUSTUS_CONFIG_PATH="$AUGUSTUS_CONFIG_PATH")
fi

mpiexec -n "$NTASKS" \
  apptainer exec "${BIND[@]}" "$IMG" \
  maker -mpi -TMP /TMP -base maker_run1 \
        maker_opts.ctl maker_bopts.ctl maker_exe.ctl

echo "[INFO] $(date)  done"
