#!/bin/bash
#SBATCH --job-name=UniProt_QC
#SBATCH --partition=pibu_el8
#SBATCH --time=00:15:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=/data/users/%u/annotation_of_eukaryotic_genome/annotation/log/UniProt_QC_%j.out
#SBATCH --error=/data/users/%u/annotation_of_eukaryotic_genome/annotation/log/UniProt_QC_%j.err

set -Eeuo pipefail

#------ CONFIG -------
FINAL_DIR="/data/users/${USER}/annotation_of_eukaryotic_genome/annotation/final"

#Protein FASTA from MAKER (renamed + filtered)
PROT_IN="${FINAL_DIR}/istisu1.bp.p_ctg.proteins.renamed.filtered.fasta"

#Best hit table (tabular BLASTP) vs UniProt (one best row per query)
BLAST_BEST="${FINAL_DIR}/blastp_vs_uniprot.besthits"

#Sanity check: protein FASTA must exist and be non-empty
[ -s "$PROT_IN" ] || { echo "[ERROR] Missing proteins FASTA: $PROT_IN"; exit 2; }

echo "[INFO] Counting total protein entries ..."
# Count number of FASTA headers ("proteins") = lines starting with '>'
N_ALL=$(grep -c '^>' "$PROT_IN" || true)
echo "[INFO] Total proteins: $N_ALL"

#If we have a UniProt besthits table, count how many unique proteins have a hit
if [ -s "$BLAST_BEST" ]; then
  #column 1 = qseqid, cut + uniq + count
  N_HIT=$(cut -f1 "$BLAST_BEST" | sort -u | wc -l | awk '{print $1}')
  #crude percentage of proteins with at least one UniProt best hit
  echo "[INFO] Proteins with UniProt best hit: $N_HIT ($(( 100 * N_HIT / (N_ALL>0?N_ALL:1) ))%)"
else
  echo "[WARN] No UniProt besthits file; cannot compute % with hits."
fi

echo "[INFO] Crude length bias check (aa lengths) ..."
#Build an ID -> length table for proteins from the FASTA
awk '
  BEGIN{RS=">"; ORS=""; FS="\n"}
  # each record is one FASTA entry (RS=">")
  NR>1{
    header=$1; split(header,a," "); id=a[1];  # get ID from header (first token)
    seq="";
    # concatenate all sequence lines, strip whitespace
    for(i=2;i<=NF;i++){ gsub(/[ \t\r]/,"",$i); seq=seq $i }
    # print: ID <tab> length
    print id "\t" length(seq) "\n"
  }
' "$PROT_IN" > "${FINAL_DIR}/prot.lengths.tsv"

#If we have UniProt besthits, split lengths into hit vs no-hit files and compare medians
if [ -s "$BLAST_BEST" ]; then
  #list of protein IDs that have a UniProt best hit
  cut -f1 "$BLAST_BEST" | sort -u > "${FINAL_DIR}/ids.hit.txt"

  #Read ids.hit.txt first (NR==FNR), mark hit[ID]=1
  #Then read prot.lengths.tsv and route length to len.hit.txt or len.nohit.txt
  awk 'NR==FNR{
           hit[$1]=1; next
       }
       {
           # $1 = ID, $2 = length
           print $2 > ( ($1 in hit) ? "'"${FINAL_DIR}"'/len.hit.txt"
                                     : "'"${FINAL_DIR}"'/len.nohit.txt" )
       }' "${FINAL_DIR}/ids.hit.txt" "${FINAL_DIR}/prot.lengths.tsv"

  #Median length for hit proteins
  echo -n "[INFO] median length (hit):   "
  sort -n "${FINAL_DIR}/len.hit.txt" | awk '{
      a[NR]=$1
    }
    END{
      if(NR%2){           # odd number of entries
        print a[(NR+1)/2]
      } else {            # even -> mean of the two middle values
        print (a[NR/2]+a[NR/2+1])/2
      }
    }'

  #Median length for no-hit proteins
  echo -n "[INFO] median length (nohit): "
  sort -n "${FINAL_DIR}/len.nohit.txt" | awk '{
      a[NR]=$1
    }
    END{
      if(NR==0){
        print "NA"        # no no-hit proteins at all
      } else if(NR%2){
        print a[(NR+1)/2]
      } else {
        print (a[NR/2]+a[NR/2+1])/2
      }
    }'
else
  echo "[WARN] Skipping length bias check (no besthits)."
fi

echo "[OK] QC done."
