#!/usr/bin/env Rscript

library(ggplot2)

#----- CONFIG ------
FINAL_DIR <- "/data/users/lland/annotation_of_eukaryotic_genome/annotation/final"

AED_VALS_FILE <- file.path(FINAL_DIR, "istisu1.bp.p_ctg.functional.AED.values.tsv")

#Read AED scores (one value per line)
aedscores <- read.table(AED_VALS_FILE, header = FALSE)
colnames(aedscores) <- "AED"

#----- ECDF PLOT -----
p_ecdf <- ggplot(aedscores, aes(x = AED)) +
  stat_ecdf(geom = "step") +
  labs(
    x = "AED",
    y = "Cumulative fraction of genes"
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = c(0, 0.25, 0.5, 0.75, 1.0)
  ) +
  geom_vline(xintercept = 0.25, linetype = "dotted") +
  geom_vline(xintercept = 0.50, linetype = "dotted") +
  theme_bw()

ggsave(
  filename = file.path(FINAL_DIR, "istisu1.bp.p_ctg.functional.AED.ecdf.ggplot.pdf"),
  plot = p_ecdf,
  width = 5,
  height = 4
)

#----- HISTOGRAM AED ------
p_hist <- ggplot(aedscores, aes(x = AED)) +
  geom_histogram(binwidth = 0.05, boundary = 0, closed = "left") +
  labs(
    x = "AED",
    y = "Number of genes"
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = c(0, 0.25, 0.5, 0.75, 1.0)
  ) +
  geom_vline(xintercept = 0.25, linetype = "dotted") +
  geom_vline(xintercept = 0.50, linetype = "dotted") +
  theme_bw()

ggsave(
  filename = file.path(FINAL_DIR, "istisu1.bp.p_ctg.functional.AED.hist.ggplot.pdf"),
  plot = p_hist,
  width = 5,
  height = 4
)
