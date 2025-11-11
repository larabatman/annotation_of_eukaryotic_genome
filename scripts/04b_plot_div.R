#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(reshape2)
  library(tidyverse)
  library(data.table)
})

## --- hardcoded paths (same as course divergence) ---
data   <- "/data/users/lland/annotation_of_eukaryotic_genome/TE_age/parseRM/istisu1.bp.p_ctg.fa.mod.out.landscape.Div.Rname.tab"
outdir <- "/data/users/lland/annotation_of_eukaryotic_genome/TE_age/landscape"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
stopifnot(file.exists(data))

## --- read from the TRUE header line; avoid the “0,0,0” duplicate-name trap ---
hdr <- readLines(data, warn = FALSE, encoding = "UTF-8")
idx <- which(grepl("(^|\\t)(\ufeff)?Rname(\\t|$)", hdr, useBytes = TRUE))[1]
if (is.na(idx)) stop("Header with 'Rname' not found in file: ", data)

rep <- fread(
  data,
  header = TRUE,
  sep = "\t",
  skip = idx - 1,          # start exactly at header
  check.names = FALSE,
  strip.white = TRUE
)
# if header still weird, uniquify to prevent setnames error
if (anyDuplicated(names(rep))) setnames(rep, make.unique(names(rep), sep = "_"))

stopifnot(ncol(rep) >= 4)
setnames(rep, old = names(rep)[1:3], new = c("Rname","Rclass","Rfam"))

## --- rename bracketed bins -> "1","2",... preserving order (like divergence script) ---
bin_cols <- grep("^\\[", names(rep), value = TRUE)
if (length(bin_cols) == 0) {
  # fallback: assume columns 4..n are bins already
  bin_cols <- names(rep)[4:ncol(rep)]
}
left <- suppressWarnings(as.integer(sub("^\\[(\\d+);.*", "\\1", bin_cols)))
ord  <- order(left, na.last = NA)
if (length(ord) != length(bin_cols)) ord <- seq_along(bin_cols)
setnames(rep, old = bin_cols[ord], new = as.character(seq_along(bin_cols)))

## --- build families EXACTLY like divergence plot, drop unknown + Helitron (course note) ---
rep <- rep %>% filter(tolower(Rfam) != "unknown")
rep$fam <- paste(rep$Rclass, rep$Rfam, sep = "/")
rep <- rep %>% filter(fam != "DNA/Helitron")

## --- melt long, drop artificial youngest bin [0;1[, keep positives only ---
rep_m <- reshape2::melt(rep, id.vars = c("Rname","Rclass","Rfam","fam"),
                        variable.name = "bin", value.name = "value")
rep_m$bin <- as.integer(as.character(rep_m$bin))
rep_m <- rep_m %>% filter(bin != 1, !is.na(value), value > 0)

## --- AGE from bin CENTER; T = K/(2r), r = 8.22e-9 subs/site/year ---
r <- 8.22e-9
rep_m <- rep_m %>%
  mutate(age_Myr = (((bin - 0.5)/100) / (2*r)) / 1e6)

## --- enforce SAME categories + order as divergence (so colors match with Brewer::Paired) ---
levels_course <- c("LTR/Copia","LTR/Gypsy",
                   "DNA/DTA","DNA/DTC","DNA/DTH","DNA/DTM","DNA/DTT",
                   "MITE/DTA","MITE/DTC","MITE/DTH","MITE/DTM")
rep_m$fam <- factor(rep_m$fam, levels = levels_course)
rep_m <- rep_m %>% filter(!is.na(fam))
stopifnot(nrow(rep_m) > 0)

## --- PLOT: identical styling; only x becomes age_Myr ---
p <- ggplot(rep_m, aes(fill = fam, x = age_Myr, weight = value/1e6)) +
  geom_bar() +
  cowplot::theme_cowplot() +
  scale_fill_brewer(palette = "Paired", na.translate = FALSE) +
  xlab("Insertion age (Myr)") +
  ylab("Sequence (Mbp)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1, size = 9))

ggsave(file.path(outdir, "TE_landscape_age.pdf"), p, width = 10, height = 5, useDingbats = FALSE)

## --- minimal CSVs (same categories/order) ---
summary_mbp <- rep_m %>%
  group_by(fam) %>%
  summarize(total_Mbp = sum(value, na.rm = TRUE)/1e6, .groups = "drop") %>%
  arrange(match(fam, levels(rep_m$fam)))
write.csv(summary_mbp, file.path(outdir, "summary_mbp_by_family.csv"), row.names = FALSE)

agg <- rep_m %>% group_by(fam, bin) %>% summarize(Mbp = sum(value, na.rm = TRUE)/1e6, .groups = "drop")
peaks <- agg %>%
  group_by(fam) %>%
  slice_max(Mbp, n = 1, with_ties = FALSE) %>%
  mutate(peak_age_Myr = (((bin - 0.5)/100) / (2*r)) / 1e6) %>%
  arrange(match(fam, levels(rep_m$fam)))
write.csv(peaks, file.path(outdir, "peaks_by_family.csv"), row.names = FALSE)

message("[OK] TE_landscape_age.pdf written to: ", outdir)
