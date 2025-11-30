#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

#----- CONFIG -----
FINAL_DIR <- "/data/users/lland/annotation_of_eukaryotic_genome/annotation/final"

len_hit_file   <- file.path(FINAL_DIR, "len.hit.txt")
len_nohit_file <- file.path(FINAL_DIR, "len.nohit.txt")

len_hit   <- scan(len_hit_file)
len_nohit <- scan(len_nohit_file)

cat("[INFO] N hit   =", length(len_hit), "\n")
cat("[INFO] N nohit =", length(len_nohit), "\n")

df <- bind_rows(
  data.frame(length = len_hit,   status = "UniProt hit"),
  data.frame(length = len_nohit, status = "No UniProt hit")
)

#Global stats
medians <- df %>%
  group_by(status) %>%
  summarise(
    n          = n(),
    mean_len   = mean(length),
    median_len = median(length)
  )

print(medians)

#Cap visuals at 99th percentile to avoid extreme tail
upper99 <- as.numeric(quantile(df$length, 0.99))

#----- DENSITY PLOT ------
p_density <- ggplot(df, aes(x = length, colour = status)) +
  geom_density() +
  coord_cartesian(xlim = c(0, upper99)) +
  labs(
    x = "Protein length (aa)",
    y = "Density"
  ) +
  theme_bw() +
  theme(
    legend.title = element_blank()
  )

ggsave(
  filename = file.path(FINAL_DIR, "uniprot_length_density.ggplot.pdf"),
  plot     = p_density,
  width    = 5,
  height   = 4
)

#----- BOXPLOT ------
p_box <- ggplot(df, aes(x = status, y = length)) +
  geom_boxplot(outlier.size = 0.4) +
  coord_cartesian(ylim = c(0, upper99)) +
  labs(
    x = "",
    y = "Protein length (aa)"
  ) +
  theme_bw() +
  theme(
    axis.text.x  = element_text(angle = 20, hjust = 1),
    legend.position = "none"
  )

ggsave(
  filename = file.path(FINAL_DIR, "uniprot_length_boxplot.ggplot.pdf"),
  plot     = p_box,
  width    = 5,
  height   = 4
)

#------ ECDF -----
p_ecdf <- ggplot(df, aes(x = length, colour = status)) +
  stat_ecdf(geom = "step") +
  coord_cartesian(xlim = c(0, upper99)) +
  labs(
    x = "Protein length (aa)",
    y = "Cumulative fraction of proteins"
  ) +
  theme_bw() +
  theme(
    legend.title = element_blank()
  )

ggsave(
  filename = file.path(FINAL_DIR, "uniprot_length_ecdf.ggplot.pdf"),
  plot     = p_ecdf,
  width    = 5,
  height   = 4
)
