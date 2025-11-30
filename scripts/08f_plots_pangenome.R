#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(tidyverse); library(scales) })

#----- CONFIG ------
#gene_counts_per_genome.csv came from the previous script
csv    <- "/data/users/lland/annotation_of_eukaryotic_genome/genespace_workingDirectory/gene_counts_per_genome.csv"
outdir <- "/data/users/lland/annotation_of_eukaryotic_genome/genespace_workingDirectory"

#focal genome for the cheese plot
focal  <- "ISTISU1"

#minimum % threshold for drawing labels on stacked % bars
LABEL_MIN_PCT <- 2

stopifnot(file.exists(csv))
dir.create(file.path(outdir, "donuts"), showWarnings = FALSE)

#read per-genome core/accessory/specific counts + percentages
tab <- readr::read_csv(csv, show_col_types = FALSE)
tab$genome <- factor(tab$genome, levels = tab$genome)

#colours for categories
pal <- c(Core="#2C7FB8", Accessory="#7570B3", `Species-specific`="#D95F02")

#----- long format for counts -----
counts_long <- tab |>
  select(genome, gene_core, gene_accessory, gene_specific) |>
  pivot_longer(-genome, names_to = "category", values_to = "count") |>
  mutate(
    #rename columns to human-readable categories
    category = recode(
      category,
      gene_core      = "Core",
      gene_accessory = "Accessory",
      gene_specific  = "Species-specific"
    ),
    category = factor(category, levels = c("Core","Accessory","Species-specific"))
  )

#----- long format for percentages -----
pct_long <- tab |>
  select(genome, percent_core, percent_accessory, percent_specific) |>
  pivot_longer(-genome, names_to = "category", values_to = "pct") |>
  mutate(
    category = recode(
      category,
      percent_core     = "Core",
      percent_accessory= "Accessory",
      percent_specific = "Species-specific"
    ),
    category = factor(category, levels = c("Core","Accessory","Species-specific"))
  )

#----- stacked bar: raw gene counts per genome ------
p_counts <- ggplot(counts_long, aes(genome, count, fill = category)) +
  geom_col(width = 0.7, color = "white", linewidth = 0.2) +
  scale_fill_manual(values = pal, drop = FALSE) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, .05))) +
  labs(
    title = "Gene set counts per genome",
    x     = NULL,
    y     = "Genes (count)",
    fill  = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position   = "top",
    axis.text.x       = element_text(angle = 45, hjust = 1),
    plot.title        = element_text(face = "bold")
  )

ggsave(file.path(outdir, "bar_counts_per_genome.pdf"),
       p_counts, width = 9.5, height = 5.5)

#------ stacked 100% bar: percentages per genome, with labels above threshold ------
labels_df <- pct_long |>
  group_by(genome) |>
  arrange(genome, category) |>
  mutate(
    ypos = cumsum(pct) - pct / 2,
    lab  = if_else(pct >= LABEL_MIN_PCT, paste0(round(pct,1), "%"), "")
  )

p_pct <- ggplot(pct_long, aes(genome, pct / 100, fill = category)) +
  geom_col(width = 0.7, color = "white", linewidth = 0.2, position = "fill") +
  geom_text(
    aes(label = ifelse(pct >= LABEL_MIN_PCT, paste0(round(pct,1), "%"), "")),
    position = position_fill(vjust = 0.52),
    color    = "white",
    fontface = "bold",
    size     = 3,
    na.rm    = TRUE
  ) +
  scale_fill_manual(values = pal, drop = FALSE) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = c(0, 0)) +
  theme_minimal(base_size = 12)

ggsave(file.path(outdir, "bar_percent_per_genome.pdf"),
       p_pct, width = 9.5, height = 5.5)

#----- cheese for focal genome -----
frow <- tab[tab$genome == focal, , drop = FALSE]
stopifnot(nrow(frow) == 1)

# values for focal genome
vals <- c(
  Core              = as.numeric(frow$gene_core),
  Accessory         = as.numeric(frow$gene_accessory),
  `Species-specific`= as.numeric(frow$gene_specific)
)

don_df <- data.frame(
  category = factor(names(vals), levels = c("Core","Accessory","Species-specific")),
  value    = as.numeric(vals)
)

RING_WIDTH <- 0.60

p_cheese <- ggplot(don_df, aes(x = 1, y = value, fill = category)) +
  geom_col(width = RING_WIDTH, color = "white", size = 0.2) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = pal, drop = FALSE) +
  xlim(0.5, 1.6) +
  theme_void(base_size = 12) +
  theme(legend.position = "top")   # donut used as legend-like summary

ggsave(file.path(outdir, paste0("cheese_", focal, "_legend_only.pdf")),
       p_cheese, width = 6, height = 6)

#----- PANGENOME FREQUENCY PLOTS (OG vs genes vs n_present) -----
freq_csv <- file.path(outdir, "pangenome_frequency_summary.csv")
if (!file.exists(freq_csv)) stop("Missing: ", freq_csv)

freq <- readr::read_csv(freq_csv, show_col_types = FALSE)

#long format: Orthogroups and Genes stacked
freq_long <- freq |>
  mutate(n_present = as.integer(n_present)) |>
  pivot_longer(cols = c(Orthogroups, Genes),
               names_to = "type",
               values_to = "count") |>
  mutate(
    type        = factor(type, levels = c("Orthogroups","Genes")),
    n_present_f = factor(n_present, levels = rev(sort(unique(n_present))))
  )

#compute % share within each series (Orthogroups vs Genes)
freq_pct <- freq_long |>
  group_by(type) |>
  mutate(pct = 100 * count / sum(count)) |>
  ungroup()

pal2 <- c("Genes" = "#0072B2", "Orthogroups" = "#D55E00")

#1) Faceted barplot: counts vs n_present, separate panel for OGs and genes
p_counts_fac <- ggplot(freq_long, aes(x = n_present_f, y = count, fill = type)) +
  geom_col(width = 0.7, color = "white", size = 0.2, show.legend = FALSE) +
  facet_wrap(~ type, scales = "free_y", nrow = 1) +
  scale_y_continuous(labels = comma) +
  scale_fill_manual(values = pal2) +
  labs(
    x = "Genomes with orthogroup (n_present; 9 = core)",
    y = "Count"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    strip.text       = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(outdir, "pangenome_freq_counts_faceted.pdf"),
       p_counts_fac, width = 9.5, height = 5.0)

#2) Line plot: % of OGs vs % of genes as a function of n_present
p_pct_lines <- ggplot(freq_pct, aes(x = n_present, y = pct, color = type)) +
  geom_line(size = 0.9) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = 1:max(freq$n_present)) +
  scale_y_continuous(labels = percent_format(scale = 1)) +
  scale_color_manual(values = pal2, name = NULL) +
  labs(
    x     = "Genomes with orthogroup (n_present)",
    y     = "Share of series (%)",
    title = "Pangenome frequency: Orthogroups vs Genes"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position   = "top",
    panel.grid.minor  = element_blank(),
    plot.title        = element_text(face = "bold")
  )

ggsave(file.path(outdir, "pangenome_freq_percent_lines.pdf"),
       p_pct_lines, width = 9.5, height = 5.0)

message("[OK] wrote: ",
        file.path(outdir, "bar_counts_per_genome.pdf"), ", ",
        file.path(outdir, "bar_percent_per_genome.pdf"), ", ",
        file.path(outdir, paste0("cheese_", focal, "_legend_only.pdf")), ", ",
        file.path(outdir, "pangenome_freq_counts_faceted.pdf"), ", ",
        file.path(outdir, "pangenome_freq_percent_lines.pdf"))
