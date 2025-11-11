#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(ggplot2) })

# -------- paths (hardcoded) --------
in_txt  <- "/data/users/lland/annotation_of_eukaryotic_genome/annotation/agat_stats/annotation_statistics.txt"
out_dir <- "/data/users/lland/annotation_of_eukaryotic_genome/annotation/agat_stats/plots"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# -------- parse AGAT report (two blocks) --------
x <- readLines(in_txt, warn = FALSE)

i_mrna <- grep("^[- ]+mrna[- ]+$", x, ignore.case = TRUE)
i_coll <- grep("^mrna have isoforms!", x, ignore.case = TRUE)
if (length(i_mrna) != 1L || length(i_coll) != 1L) {
  stop("Could not find expected AGAT section headers in: ", in_txt)
}

block1 <- x[(i_mrna + 1):(i_coll - 1)]           # all isoforms
block2 <- x[(i_coll + 1):length(x)]              # collapsed isoforms

parse_pairs <- function(lines, set_name) {
  # keep lines that end with a number
  keep <- grepl("\\d\\s*$", lines)
  y <- lines[keep]
  if (!length(y)) return(data.frame())
  metric <- sub("\\s+[0-9.]+\\s*$", "", y)
  value  <- as.numeric(sub(".*?([0-9.]+)\\s*$", "\\1", y))
  data.frame(set = set_name,
             metric = tolower(trimws(metric)),
             value = value,
             stringsAsFactors = FALSE)
}

d1 <- parse_pairs(block1, "all_isoforms")
d2 <- parse_pairs(block2, "collapsed_isoforms")
df <- rbind(d1, d2)

getv <- function(d, set, key) {
  v <- d$value[d$set == set & d$metric == tolower(key)]
  if (length(v)) v[1] else NA_real_
}

# -------- derived summaries --------
# Counts
want_counts <- c("number of gene","number of mrna","number of exon","number of cds")
counts <- subset(df, metric %in% tolower(want_counts))
nice_name <- c("number of gene"="Genes","number of mrna"="mRNAs","number of exon"="Exons","number of cds"="CDS")
counts$metric <- factor(nice_name[counts$metric], levels = c("Genes","mRNAs","Exons","CDS"))

# UTR coverage (both / one / none)
mk_utr <- function(set_name){
  n_mrna   <- getv(df, set_name, "Number of mrna")
  n_both   <- getv(df, set_name, "Number of mrnas with utr both sides")
  n_atleast<- getv(df, set_name, "Number of mrnas with at least one utr")
  data.frame(
    set = set_name,
    category = factor(c("both UTRs","one UTR","no UTR"), levels = c("no UTR","one UTR","both UTRs")),
    count = c(n_both, max(n_atleast - n_both, 0), max(n_mrna - n_atleast, 0)),
    total = n_mrna
  )
}
utr <- rbind(mk_utr("all_isoforms"), mk_utr("collapsed_isoforms"))
utr$perc <- 100 * utr$count / utr$total

# Single-exon percentages
mk_single <- function(set_name){
  data.frame(
    set = set_name,
    entity = c("Genes","mRNAs"),
    count_single = c(getv(df, set_name, "Number of single exon gene"),
                     getv(df, set_name, "Number of single exon mrna")),
    total = c(getv(df, set_name, "Number of gene"),
              getv(df, set_name, "Number of mrna"))
  )
}
single <- rbind(mk_single("all_isoforms"), mk_single("collapsed_isoforms"))
single$pct <- 100 * single$count_single / single$total

# Structure means
struct <- rbind(
  data.frame(set="all_isoforms",
             metric=c("Exons per mRNA","Introns (CDS) per mRNA"),
             value=c(getv(df,"all_isoforms","mean exons per mrna"),
                     getv(df,"all_isoforms","mean introns in cdss per mrna"))),
  data.frame(set="collapsed_isoforms",
             metric=c("Exons per mRNA","Introns (CDS) per mRNA"),
             value=c(getv(df,"collapsed_isoforms","mean exons per mrna"),
                     getv(df,"collapsed_isoforms","mean introns in cdss per mrna")))
)
struct$metric <- factor(struct$metric, levels=c("Exons per mRNA","Introns (CDS) per mRNA"))

# -------- plotting helpers --------
save_plot <- function(p, file, w=8, h=5){
  ggsave(filename=file.path(out_dir, file), plot=p, width=w, height=h, units="in", dpi=150)
  message("Saved: ", file.path(out_dir, file))
}

# -------- P1: counts --------
p_counts <- ggplot(counts, aes(x=metric, y=value, fill=set)) +
  geom_col(position=position_dodge(width=0.75), width=0.7) +
  scale_y_continuous(labels = scales::label_number(big.mark=",")) +
  labs(title="AGAT counts", x=NULL, y="Count") +
  scale_fill_manual(values=c(all_isoforms="#999999", collapsed_isoforms="#377eb8"),
                    name=NULL, labels=c("All isoforms","Collapsed")) +
  theme_minimal(base_size=12)
save_plot(p_counts, "agat_counts.png")

# -------- P2: UTR coverage --------
p_utr <- ggplot(utr, aes(x=set, y=perc, fill=category)) +
  geom_col(position="stack", width=0.7) +
  scale_fill_manual(values=c("no UTR"="#e41a1c","one UTR"="#ff7f00","both UTRs"="#4daf4a"), name=NULL) +
  scale_y_continuous(limits=c(0,100), breaks=seq(0,100,20)) +
  labs(title="UTR coverage per mRNA", x=NULL, y="% of mRNAs") +
  theme_minimal(base_size=12) +
  coord_flip()
save_plot(p_utr, "agat_UTR_coverage.png")

# -------- P3: single-exon percentages --------
p_single <- ggplot(single, aes(x=entity, y=pct, fill=set)) +
  geom_col(position=position_dodge(width=0.7), width=0.65) +
  labs(title="Single-exon fraction", x=NULL, y="% single-exon") +
  scale_fill_manual(values=c(all_isoforms="#999999", collapsed_isoforms="#377eb8"),
                    name=NULL, labels=c("All isoforms","Collapsed")) +
  theme_minimal(base_size=12)
save_plot(p_single, "agat_single_exon_pct.png")

# -------- P4: structure means --------
p_struct <- ggplot(struct, aes(x=metric, y=value, fill=set)) +
  geom_col(position=position_dodge(width=0.7), width=0.65) +
  labs(title="Structure means", x=NULL, y="Mean per mRNA") +
  scale_fill_manual(values=c(all_isoforms="#999999", collapsed_isoforms="#377eb8"),
                    name=NULL, labels=c("All isoforms","Collapsed")) +
  theme_minimal(base_size=12)
save_plot(p_struct, "agat_structure_means.png")

message("[OK] Plots in: ", out_dir)
