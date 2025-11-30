#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

#----- CONFIG ------
base      <- "/data/users/lland/annotation_of_eukaryotic_genome/annotation/final"
prot_fa   <- file.path(base, "istisu1.bp.p_ctg.proteins.renamed.filtered.fasta")

uni_all   <- file.path(base, "blastp_vs_uniprot.outfmt6")
uni_best  <- file.path(base, "blastp_vs_uniprot.besthits")
tair_all  <- file.path(base, "blastp_vs_TAIR10.outfmt6")
tair_best <- file.path(base, "blastp_vs_TAIR10.besthits")

#UniProt-annotated files (for sanity peek)
fasta_uni <- file.path(base, "istisu1.bp.p_ctg.proteins.renamed.filtered.fasta.UniProt")
gff_uni   <- file.path(base, "filtered.genes.renamed.gff3.UniProt.gff3")

out_dir   <- file.path(base, "func_annot_report")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

#TAIR IDs / symbols to sanity check (FLC, RCA, etc.)
targets <- c("AT5G10140", "FLC", "RCA")

#------ HELPER ------

#BLAST outfmt6 column names
cols12 <- c("qseqid","sseqid","pident","length","mismatch","gapopen",
            "qstart","qend","sstart","send","evalue","bitscore")

#read a BLAST table
read_blast <- function(p) {
  read_tsv(p,
           col_names = cols12,
           col_types = cols(.default = col_guess()),
           progress  = FALSE)
}

#parse FASTA and return named vector: id -> length(aa)
fasta_lengths <- function(p) {
  x     <- read_lines(p)
  is_hd <- startsWith(x, ">")
  ids   <- sub("^>(\\S+).*","\\1", x[is_hd])

  lens <- integer(length(ids))
  cur  <- 0L
  k    <- 0L

  for (i in seq_along(x)) {
    if (is_hd[i]) {
      if (i > 1) lens[k] <- cur  # close previous seq
      k   <- k + 1L
      cur <- 0L
    } else {
      cur <- cur + nchar(gsub("\\s", "", x[i]))
    }
  }
  lens[k] <- cur
  stats::setNames(lens, ids)
}

pct <- function(x) sprintf("%.2f%%", 100 * x)

#grep a few lines matching a pattern in a text file (sanity peeks)
peek <- function(path, pattern, n = 5) {
  if (!file.exists(path)) return(character())
  y <- grep(pattern,
            read_lines(path, n_max = 5000),
            value = TRUE,
            ignore.case = TRUE)
  head(y, n)
}

#load a BLAST pair (all + besthits) as list
load_set <- function(allp, bestp, label) {
  if (!file.exists(allp) || !file.exists(bestp)) return(NULL)
  all  <- read_blast(allp)  %>% mutate(set = label)
  best <- read_blast(bestp) %>% mutate(set = label)
  list(all = all, best = best)
}

#----- LOAD CORE DATA ------

if (!file.exists(prot_fa)) {
  stop("Protein FASTA not found: ", prot_fa)
}

qlen  <- fasta_lengths(prot_fa)          # id -> length
total <- length(qlen)                    # total # proteins

UNI <- load_set(uni_all,  uni_best,  "UniProt")
TAI <- load_set(tair_all, tair_best, "TAIR10")

#----- COVERGAE per DBs -----

cov_tbl <- tibble(set = character(), covered = integer(),
                  total = integer(), pct = numeric())

if (!is.null(UNI)) {
  cov_tbl <- bind_rows(
    cov_tbl,
    tibble(
      set     = "UniProt",
      covered = n_distinct(UNI$best$qseqid),  # proteins with best hit
      total   = total,
      pct     = n_distinct(UNI$best$qseqid) / total
    )
  )
}

if (!is.null(TAI)) {
  cov_tbl <- bind_rows(
    cov_tbl,
    tibble(
      set     = "TAIR10",
      covered = n_distinct(TAI$best$qseqid),
      total   = total,
      pct     = n_distinct(TAI$best$qseqid) / total
    )
  )
}

#----- SET TALLIES ------

ids_all  <- names(qlen)   # all proteins
ids_uni  <- if (!is.null(UNI)) unique(UNI$best$qseqid) else character()
ids_tai  <- if (!is.null(TAI)) unique(TAI$best$qseqid) else character()

ids_inter     <- intersect(ids_uni, ids_tai)              # hit in both
ids_tair_only <- setdiff(ids_tai, ids_uni)                # only TAIR
ids_unip_only <- setdiff(ids_uni, ids_tai)                # only UniProt
ids_nohit     <- setdiff(ids_all, union(ids_uni, ids_tai))# no hit

set_tbl <- tibble(
  category = c("UniProt ∩ TAIR",
               "TAIR-only",
               "UniProt-only",
               "No hit (neither)"),
  n   = c(length(ids_inter),
          length(ids_tair_only),
          length(ids_unip_only),
          length(ids_nohit)),
  pct = c(length(ids_inter),
          length(ids_tair_only),
          length(ids_unip_only),
          length(ids_nohit)) / total
)

#----- LENGTH STATISTICS: HIT vs NO-HIT ------

#build a data frame of lengths, labelled hit/no_hit per set
mk_len_df <- function(best, label) {
  ids_hit <- unique(best$qseqid)
  tibble(
    id  = names(qlen),
    len = as.integer(qlen),
    hit = if_else(names(qlen) %in% ids_hit, "hit", "no_hit"),
    set = label
  )
}

len_both <- bind_rows(
  if (!is.null(UNI)) mk_len_df(UNI$best, "UniProt") else NULL,
  if (!is.null(TAI)) mk_len_df(TAI$best, "TAIR10") else NULL
)

len_stats  <- tibble()
wilcox_tbl <- tibble()

if (nrow(len_both) > 0) {
  #per set × hit: n, mean, median, sd, range
  len_stats <- len_both %>%
    group_by(set, hit) %>%
    summarise(
      n          = n(),
      mean_len   = mean(len),
      median_len = median(len),
      sd_len     = sd(len),
      min_len    = min(len),
      max_len    = max(len),
      .groups    = "drop"
    )
}

#----- SANITY PEEKS ------

fasta_heads <- peek(fasta_uni, "^>")                         # UniProt-annotated FASTA headers
gff_anno    <- peek(gff_uni, "Note=|Dbxref=|product=")       # GFF attributes with function

#------ TARGET GENE CHECKS ------

targets_re  <- paste0("(", paste(targets, collapse = "|"), ")")
target_hits <- list()

#find TAIR best hits whose subject matches any of the target IDs/symbols
if (file.exists(tair_best) && file.info(tair_best)$size > 0) {
  tb <- read_blast(tair_best)
  tb$match <- grepl(targets_re, tb$sseqid, ignore.case = TRUE)
  th <- tb %>%
    filter(match) %>%
    select(qseqid, sseqid, pident, length, evalue, bitscore) %>%
    distinct()
  target_hits$tair_best <- th
}

##ook for targets in annotated FASTA headers / GFF attributes
target_hits$fasta_uni <- peek(fasta_uni, targets_re, n = 20)
target_hits$gff_uni   <- peek(gff_uni,   targets_re, n = 20)

#----- SUMMARY ------

sum_txt <- file.path(out_dir, "summary.txt")
con     <- file(sum_txt, "wt")

#total proteins
writeLines(sprintf("Total proteins (from FASTA): %d", total), con)

#coverage
if (nrow(cov_tbl) > 0) {
  writeLines("\nCoverage (best hit per protein):", con)
  for (i in seq_len(nrow(cov_tbl))) {
    r <- cov_tbl[i, ]
    line <- sprintf(
      "  %s: %d/%d (%s)",
      r$set, r$covered, r$total, pct(r$pct)
    )
    writeLines(line, con)
  }
} else {
  writeLines("\nNo BLAST besthits found (UniProt/TAIR10).", con)
}

#set tallies (intersection/only/no-hit)
writeLines("\nSet tallies (exact):", con)
for (i in seq_len(nrow(set_tbl))) {
  r <- set_tbl[i, ]
  line <- sprintf("  %-22s %6d (%s)", r$category, r$n, pct(r$pct))
  writeLines(line, con)
}

#length stats
writeLines("\nLength statistics by hit status (per set):", con)
if (nrow(len_stats) > 0) {
  for (i in seq_len(nrow(len_stats))) {
    r <- len_stats[i, ]
    line <- sprintf(
      "  %s | %s: n=%d, mean=%.1f aa, median=%.1f aa, sd=%.1f, range=[%d, %d] aa",
      r$set, r$hit, r$n,
      r$mean_len, r$median_len, r$sd_len,
      r$min_len, r$max_len
    )
    writeLines(line, con)
  }
} else {
  writeLines("  (no length stats available)", con)
}

# sanity peek of annotated FASTA / GFF
writeLines("\nSanity peek (UniProt-annotated FASTA headers):", con)
if (length(fasta_heads)) {
  writeLines(paste0("  ", fasta_heads), con)
} else {
  writeLines("  (missing or empty)", con)
}

writeLines("\nSanity peek (UniProt-annotated GFF attributes):", con)
if (length(gff_anno)) {
  writeLines(paste0("  ", gff_anno), con)
} else {
  writeLines("  (missing or no Note/Dbxref/product)", con)
}

# target genes
writeLines(
  paste0("\nTargets searched (TAIR IDs/symbols): ",
         paste(targets, collapse = ", ")),
  con
)

writeLines("\nTarget matches (TAIR besthits):", con)
if (!is.null(target_hits$tair_best) && nrow(target_hits$tair_best) > 0) {
  utils::write.table(
    target_hits$tair_best,
    file      = con,
    append    = TRUE,
    quote     = FALSE,
    row.names = FALSE,
    sep       = "\t"
  )
} else {
  writeLines("  none found in TAIR besthits", con)
}

writeLines("\nTarget matches in annotated FASTA headers:", con)
if (length(target_hits$fasta_uni)) {
  writeLines(paste0("  ", target_hits$fasta_uni), con)
} else {
  writeLines("  none", con)
}

writeLines("\nTarget matches in annotated GFF attributes:", con)
if (length(target_hits$gff_uni)) {
  writeLines(paste0("  ", target_hits$gff_uni), con)
} else {
  writeLines("  none", con)
}

close(con)

message("[DONE] Summary → ", sum_txt)
