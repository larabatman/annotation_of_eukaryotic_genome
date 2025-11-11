#!/usr/bin/env bash
#SBATCH --job-name=tesorter_results
#SBATCH --partition=pibu_el8
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --output=%x_%j.out

set -euo pipefail

#----- CONFIG -----
GENOME="${GENOME:-/data/users/lland/annotation_of_eukaryotic_genome/assembly/hifiasm/istisu1.bp.p_ctg.fa}" #path to FASTA assembly 
OUTDIR="${OUTDIR:-/data/users/lland/annotation_of_eukaryotic_genome/EDTA_annotation}" #EDTA output root 
THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-8}}" #default threads 
IMG="${IMG:-/data/courses/assembly-annotation-course/CDS_annotation/containers/TEsorter_1.3.0.sif}" #apptainer image for TEsorter 1.3.0
BIND="${BIND:-/data}"   # bind root visible to container
RESULTS_DIR="${RESULTS_DIR:-}"  # if empty, will default below

#----- LOCATING EDTA LTR CANDIDATES AND SANITY CHECKS ------
PREFIX="$(basename "$GENOME")" #istisu1.bp.p_ctg.fa for example
RAW="${OUTDIR}/${PREFIX}.mod.EDTA.raw" #EDTA's raw subfolder
LTRFA="${RAW}/${PREFIX}.mod.LTR.raw.fa" #EDTA's raw LTR candidates, the TEsorter input
# Existence checks for genome, outdir, image and LTR FASTA:
[ -s "$GENOME" ] || { echo "ERROR: GENOME not found: $GENOME" >&2; exit 1; }
[ -d "$OUTDIR" ] || { echo "ERROR: OUTDIR not found: $OUTDIR" >&2; exit 1; }
[ -s "$IMG" ]    || { echo "ERROR: TEsorter image not found: $IMG" >&2; exit 1; }
[ -s "$LTRFA" ]  || { echo "ERROR: Missing LTR FASTA: $LTRFA (EDTA LTR step not done?)" >&2; exit 2; }

#----- RESULTS DIR -----
RESULTS_DIR="${RESULTS_DIR:-$OUTDIR/${PREFIX}.LTR_analysis}"
mkdir -p "$RESULTS_DIR"
echo "[INFO] Writing results in: $RESULTS_DIR"

#----- RUN TESORTER -----
cd "$RESULTS_DIR"
echo "[INFO] Running TEsorter on: $LTRFA  (threads=$THREADS)"
#run TEstorter on the EDTA LTR candidates against REXdb-plamt HMMs, assigning superfamily and clades
apptainer exec --bind "$BIND" "$IMG" \
  TEsorter "$LTRFA" -db rexdb-plant -p "$THREADS"

#Output sanity check
CLS="${RESULTS_DIR}/$(basename "$LTRFA").rexdb-plant.cls.tsv"
[ -s "$CLS" ] || { echo "ERROR: expected TSV missing: $CLS"; exit 3; }
GFF_INTACT_SRC="$RAW/${PREFIX}.mod.LTR.intact.raw.gff3"
[ -s "$GFF_INTACT_SRC" ] || { echo "ERROR: intact GFF missing: $GFF_INTACT_SRC"; exit 4; }

#Symlinks to match the TSV and GFF3 files to the R plotting script 
ln -sf "$(basename "$CLS")" "${RESULTS_DIR}/genomic.fna.mod.LTR.raw.fa.rexdb-plant.cls.tsv"
echo "[INFO] Linked: genomic.fna.mod.LTR.raw.fa.rexdb-plant.cls.tsv -> $(basename "$CLS")"

ln -sf "$GFF_INTACT_SRC" "${RESULTS_DIR}/genomic.fna.mod.LTR.intact.raw.gff3"
echo "[INFO] Linked GFF -> genomic.fna.mod.LTR.intact.raw.gff3"


#----- SUMMARY STATISTICS -----
CLS_TSV="$CLS" #shorthand
GFF_INTACT="$RAW/${PREFIX}.mod.LTR.intact.raw.gff3" #EDTA's intact LTR retrotransposons GFF
N_FASTA=$(grep -c '^>' "$LTRFA" || echo 0) #number of input sequences to TEsorter

#if there is an intact GFF file, count features that represent intect elements by skipping wrapper features to get N_INTACT
if [ -s "$GFF_INTACT" ]; then
  N_INTACT=$(awk -F'\t' '!/^#/ && $3!="long_terminal_repeat" && $3!="repeat_region" && $3!="target_site_duplication"{c++} END{print c+0}' "$GFF_INTACT")
else
  N_INTACT="NA"
fi

#Prepare output filenames
SF_OUT="${RESULTS_DIR}/tesorter_counts_per_superfamily.tsv"
CL_OUT="${RESULTS_DIR}/tesorter_counts_per_clade.tsv"
SUM_OUT="${RESULTS_DIR}/tesorter_summary.txt"

#Detect header columns; for rows where the TE name has _INT for intact label, strip and keep unique intact elements per superfamily and clade. 
# Write two count table and a summary with number of input sequences, number of unique _INT-classified elements and classification rate as well as number of intact LTR-RTs from EDTA's GFF as cross-check 
awk -F'\t' -v sf_out="$SF_OUT" -v cl_out="$CL_OUT" -v n_fa="$N_FASTA" -v n_intact="$N_INTACT" '
  BEGIN{OFS="\t"}
  NR==1{
    for(i=1;i<=NF;i++){
      h=tolower($i)
      if(h=="te"||h=="name") te=i
      if(h=="superfamily") sf=i
      if(h=="clade") cl=i
    }
    next
  }
  {
    teval=$te
    if (teval ~ /_INT/) {
      name=teval; sub(/#.*$/,"",name); sub(/_INT$/,"",name)
      if (!(name in seen)) {
        seen[name]=1; total++
        s=($sf==""?"NA":$sf); c=($cl==""?"Unclassified":$cl)
        sf_cnt[s]++; cl_cnt[c]++
      }
    }
  }
  END{
    print "Superfamily","count" > sf_out
    for(k in sf_cnt) print k,sf_cnt[k] >> sf_out
    print "Clade","count" > cl_out
    for(k in cl_cnt) print k,cl_cnt[k] >> cl_out
    print "TEsorter input (FASTA seqs):", n_fa > "'$SUM_OUT'"
    print "Unique INT-classified TEs  :", (total+0) >> "'$SUM_OUT'"
    print "Classification rate        :", (n_fa>0? sprintf("%.2f%%",100*total/n_fa):"NA") >> "'$SUM_OUT'"
    print "Intact LTR-RTs in GFF      :", n_intact >> "'$SUM_OUT'"
  }' "$CLS_TSV"

echo "[INFO] Wrote:"
echo "  $SUM_OUT"
echo "  $SF_OUT"
echo "  $CL_OUT"
