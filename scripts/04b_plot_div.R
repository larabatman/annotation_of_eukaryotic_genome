#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(reshape2)
  library(dplyr)
  library(ggplot2)
  library(cowplot)
})

#----- CONFIG -----
data   <- "istisu1.bp.p_ctg.fa.mod.out.landscape.Div.Rname.tab"  # parseRM output
outdir <- "TE_age/landscape"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

#----- READ TABLE (one row per TE family; many divergence bins) -----
rep <- fread(data, sep = "\t", header = TRUE, check.names = FALSE)

#first 3 columns = TE name, class, family
setnames(rep, old = names(rep)[1:3], new = c("Rname","Rclass","Rfam"))

#identify divergence-bin columns (the ones that start with '[') and rename them as 1,2,3,... so they are easy to treat as numeric bins
bin_cols <- grep("^\\[", names(rep), value = TRUE)
if (length(bin_cols) == 0) bin_cols <- names(rep)[4:ncol(rep)]
setnames(rep, old = bin_cols, new = as.character(seq_along(bin_cols)))

#----- TE FAMILY KEY AND FILTER -----
#drop "unknown" families; build fam = "Rclass/Rfam" (LTR/Copia)
rep <- rep %>% filter(tolower(Rfam) != "unknown")
rep <- rep %>% mutate(fam = paste(Rclass, Rfam, sep = "/"))

#remove Helitrons for this plot
rep <- rep %>% filter(fam != "DNA/Helitron")

#keep only the TE families we care about and fix their plotting order
levels_course <- c("LTR/Copia","LTR/Gypsy",
                   "DNA/DTA","DNA/DTC","DNA/DTH","DNA/DTM","DNA/DTT",
                   "MITE/DTA","MITE/DTC","MITE/DTH","MITE/DTM")
rep$fam <- factor(rep$fam, levels = levels_course)
rep      <- rep %>% filter(!is.na(fam))

#----- wide -> long: one row per family × bin -----
rep_m <- melt(rep,
              id.vars      = c("Rname","Rclass","Rfam","fam"),
              variable.name = "bin",
              value.name    = "value")

#convert bin label back to integer (1,2,3,..)
rep_m$bin <- as.integer(as.character(rep_m$bin))

#drop artificial youngest bin (bin 1) and rows with 0 bp
rep_m <- rep_m %>%
  filter(bin != 1, !is.na(value), value > 0)

#----- CONVERT DIVERGENCE TO AGE ESTIMATION (Myr) -----
#bin k represents divergence K ≈ (k - 0.5)% = (k-0.5)/100
#T = K / (2r) with r = 8.22e-9 subs/site/year
r <- 8.22e-9
rep_m <- rep_m %>%
  mutate(age_Myr = (((bin - 0.5)/100) / (2 * r)) / 1e6)

#----- PLOT STACKED TE VS AGE -----
p <- ggplot(rep_m, aes(x = age_Myr,
                       fill = fam,
                       weight = value / 1e6)) +   # convert bp -> Mbp
  geom_bar() +
  cowplot::theme_cowplot() +
  scale_fill_brewer(palette = "Paired", na.translate = FALSE) +
  xlab("Insertion age (Myr)") +
  ylab("Sequence (Mbp)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1, size = 9))

ggsave(file.path(outdir, "TE_landscape_age.pdf"),
       p, width = 10, height = 5, useDingbats = FALSE)

#----- summaries: total Mbp per family -----
summary_mbp <- rep_m %>%
  group_by(fam) %>%
  summarise(total_Mbp = sum(value, na.rm = TRUE) / 1e6, .groups = "drop") %>%
  arrange(match(fam, levels(rep_m$fam)))

write.csv(summary_mbp,
          file.path(outdir, "summary_mbp_by_family.csv"),
          row.names = FALSE)

#----- per-family peak age (bin with max Mbp) -----
agg <- rep_m %>%
  group_by(fam, bin) %>%
  summarise(Mbp = sum(value, na.rm = TRUE) / 1e6, .groups = "drop")

peaks <- agg %>%
  group_by(fam) %>%
  slice_max(Mbp, n = 1, with_ties = FALSE) %>%
  mutate(peak_age_Myr = (((bin - 0.5)/100) / (2 * r)) / 1e6) %>%
  arrange(match(fam, levels(rep_m$fam)))

write.csv(peaks,
          file.path(outdir, "peaks_by_family.csv"),
          row.names = FALSE)
