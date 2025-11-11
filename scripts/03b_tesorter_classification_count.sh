#!/usr/bin/env bash
#SBATCH --job-name=tesorter_CG
#SBATCH --partition=pibu_el8
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --output=%x_%j.out

set -euo pipefail

#----- CONFIG -----
WORKDIR=/data/users/lland/annotation_of_eukaryotic_genome
GENOME=$WORKDIR/assembly/hifiasm/istisu1.bp.p_ctg.fa
EDTA_DIR=$WORKDIR/EDTA_annotation
TEANNO=$EDTA_DIR/$(basename "$GENOME").mod.EDTA.TEanno.gff3
REFINE=$WORKDIR/TEsorter_refine
COP="$REFINE/Copia_sequences.fa"
GYP="$REFINE/Gypsy_sequences.fa"

#----- TEsorter APPTAINER -----
IMG=/data/courses/assembly-annotation-course/CDS_annotation/containers/TEsorter_1.3.0.sif
THREADS=${SLURM_CPUS_PER_TASK:-8}

# Sanity checks for inputs
mkdir -p "$REFINE"
[ -s "$COP" ] || { echo "ERROR: missing $COP (run step 3a)"; exit 1; }
[ -s "$GYP" ] || { echo "ERROR: missing $GYP (run step 3a)"; exit 1; }

cd "$REFINE"

#----- RUN TESORTER ON COPIA AND GYPSY ------
# REXdb-plant DB
apptainer exec --bind "$WORKDIR" "$IMG" TEsorter "$COP" -db rexdb-plant -p "$THREADS"
apptainer exec --bind "$WORKDIR" "$IMG" TEsorter "$GYP" -db rexdb-plant -p "$THREADS"

#----- CHECK NAMING -----
COP_B=$(basename "$COP"); GYP_B=$(basename "$GYP")
COP_CLS="$REFINE/${COP_B}.rexdb-plant.cls.tsv";  [ -s "$COP_CLS" ] || COP_CLS="$REFINE/${COP_B%.fa}.cls.tsv"
GYP_CLS="$REFINE/${GYP_B}.rexdb-plant.cls.tsv";  [ -s "$GYP_CLS" ] || GYP_CLS="$REFINE/${GYP_B%.fa}.cls.tsv"
[ -s "$COP_CLS" ] || { echo "ERROR: Copia cls tsv not found in $REFINE"; exit 2; }
[ -s "$GYP_CLS" ] || { echo "ERROR: Gypsy cls tsv not found in $REFINE"; exit 3; }

#----- NORMALIZE TESORTER TSV VARIANTS ----
# Find columns by name; or parse from a single Classification field; deduplicate (TE, SF, Clade)
ALL_CLS="$REFINE/tesorter_all_class.tsv"
awk -F'\t' 'BEGIN{OFS="\t"}
FNR==1{
  delete idx
  for(i=1;i<=NF;i++){
    h=tolower($i)
    if(h=="te" || h=="name" || h=="query") idx["te"]=i
    if(h=="superfamily") idx["sf"]=i
    if(h=="clade") idx["cl"]=i
    if(h=="classification") idx["cls"]=i
  }
  if (idx["te"] || idx["sf"] || idx["cl"] || idx["cls"]) next
}
{
  te = (idx["te"] ? $idx["te"] : $1)
  sf = (idx["sf"] ? $idx["sf"] : "")
  cl = (idx["cl"] ? $idx["cl"] : "")
  if (sf=="" && idx["cls"]) {
    cls=$idx["cls"]
    gsub(/[[:space:]]+/,"",cls)
    gsub(/\|/,"/",cls)    # allow LTR|Copia|Angela
    n=split(cls, a, "/")
    if (n>=2) {
      tl=tolower(a[1])
      if (tl=="ltr" || tl=="dna" || tl=="rc") sf=a[2]; else sf=a[1]
    } else if (n==1) { sf=a[1] }
    if (n>=3) cl=a[3]
  }
  # Clean TE id: drop trailing #..., remove _INT suffix
  sub(/#.*/,"",te); sub(/_INT$/,"",te)
  if (sf=="") sf="Unknown"
  if (cl=="") cl="Unclassified"
  print te, sf, cl
}' "$COP_CLS" "$GYP_CLS" \
  | LC_ALL=C sort -u > "$ALL_CLS"

echo "[OK] Wrote: $ALL_CLS  (TE,Superfamily,Clade)"

#----- COUNT CONSENSUS SEQUENCES (library entries) -----
LIB_COUNTS="$REFINE/counts_library_by_clade.tsv"
awk -F'\t' 'BEGIN{OFS="\t"} {c[$3]++} END{print "clade","n_library"; for(k in c) print k,c[k]}' "$ALL_CLS" \
  | LC_ALL=C sort -k2,2nr > "$LIB_COUNTS"
echo "[OK] Wrote: $LIB_COUNTS"

# Also library counts by (Superfamily, Clade)
PAIR_COUNTS_LIB="$REFINE/counts_library_by_superfamily_clade.tsv"
awk -F'\t' 'BEGIN{OFS="\t"} {key=$2"/"$3; c[key]++}
  END{print "superfamily","clade","n_library";
      for(k in c){split(k,a,"/"); print a[1],a[2],c[k]}}' \
  "$ALL_CLS" | LC_ALL=C sort -k1,1 -k3,3nr > "$PAIR_COUNTS_LIB"
echo "[OK] Wrote: $PAIR_COUNTS_LIB"

#----- MAP CLADES BACK TO EDTA TEanno.gff3 (copy counts in genome) -----
if [ -s "$TEANNO" ]; then
  # Extract clean TE names from GFF, drop meta features and normalize common suffixes
  LC_ALL=C awk -F'\t' 'BEGIN{OFS="\t"}
    /^#/ {next}
    $3=="long_terminal_repeat" || $3=="repeat_region" || $3=="target_site_duplication" {next}
    {
      name=""
      if (match($9, /(^|;)Name=([^;]+)/, m)) name=m[2]
      if (name!="") {
        # Normalize common EDTA/TEsorter naming quirks so join hits:
        gsub(/(%2[Dd])/, "-", name)                    # decode "-"" if percent-encoded
        sub(/#.*/,"",name)                             # strip trailing #...
        gsub(/(_|-)?LTR[53]?$/,"",name)                # remove -LTR, -LTR5/3, _LTR
        gsub(/(_|-)?I(NT)?$/,"",name)                  # remove _I, -I, _INT, -INT
        print name
      }
    }' "$TEANNO" | LC_ALL=C sort -u > "$REFINE/__names_in_gff.txt"

  # Join TE names in GFF with classification table
  LC_ALL=C sort -k1,1 "$ALL_CLS" > "$REFINE/__cls_sorted.tsv"
  join -t $'\t' -1 1 -2 1 "$REFINE/__names_in_gff.txt" "$REFINE/__cls_sorted.tsv" \
    > "$REFINE/__name_sf_cl.tsv" || true

  # Stats: how many names mapped?
  GFF_N=$(wc -l < "$REFINE/__names_in_gff.txt" || echo 0)
  MAP_N=$(wc -l < "$REFINE/__name_sf_cl.tsv"  || echo 0)
  echo "[INFO] GFF names: $GFF_N ; mapped with class: $MAP_N"

  # Counts per clade (elements, not library)
  ANN_COUNTS="$REFINE/counts_TEanno_by_clade.tsv"
  awk -F'\t' 'BEGIN{OFS="\t"} {cl=$3; c[cl]++}
       END{print "clade","n_TEanno"; for(k in c) print k,c[k]}' \
    "$REFINE/__name_sf_cl.tsv" | LC_ALL=C sort -k2,2nr > "$ANN_COUNTS"
  echo "[OK] Wrote: $ANN_COUNTS"


# Counts per (Superfamily, Clade) — header first, then sorted rows
  ANN_PAIR_COUNTS="$REFINE/counts_TEanno_by_superfamily_clade.tsv"
  {
    echo -e "superfamily\tclade\tn_TEanno"
    awk -F'\t' 'BEGIN{OFS="\t"} {key=$2"/"$3; c[key]++}
        END{for(k in c){split(k,a,"/"); print a[1],a[2],c[k]}}' \
      "$REFINE/__name_sf_cl.tsv" \
    | LC_ALL=C sort -k1,1 -k3,3nr
  } > "$ANN_PAIR_COUNTS"
  echo "[OK] Wrote: $ANN_PAIR_COUNTS"

  # Totals per superfamily (overall Copia vs Gypsy elements)
  ANN_SF_COUNTS="$REFINE/counts_TEanno_by_superfamily.tsv"
  awk -F'\t' 'NR>1{sf[$1]+=$3} END{print "superfamily","n_TEanno"; for(s in sf) print s,sf[s]}' \
    "$ANN_PAIR_COUNTS" | LC_ALL=C sort -k2,2nr > "$ANN_SF_COUNTS"
  echo "[OK] Wrote: $ANN_SF_COUNTS"

  # Quick readable top lists
  echo
  echo "Top clades by elements (EDTA):"
  head -n 15 "$ANN_COUNTS" || true

  echo
  echo "[Copia] top clades:"
  awk -F'\t' 'NR>1 && $1=="Copia"{print $2,$3}' "$ANN_PAIR_COUNTS" | head -n 10 || true

  echo
  echo "[Gypsy] top clades:"
  awk -F'\t' 'NR>1 && $1=="Gypsy"{print $2,$3}' "$ANN_PAIR_COUNTS" | head -n 10 || true

  # Cleanup temps
  rm -f "$REFINE/__names_in_gff.txt" "$REFINE/__cls_sorted.tsv"
  # Keep __name_sf_cl.tsv for debugging; comment next line if you want to inspect it
  # rm -f "$REFINE/__name_sf_cl.tsv"

else
  echo "[WARN] TEanno GFF not found: $TEANNO — skipping per-clade counts from annotation."
fi

echo
echo "Top clades (library entries):"
head -n 10 "$LIB_COUNTS" || true
