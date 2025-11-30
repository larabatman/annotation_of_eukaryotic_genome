#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(tidyverse))

#C------ CONFIG ------
args <- commandArgs(trailingOnly = TRUE)

#working directory where pangenome_matrix.rds lives
wd <- if (length(args) >= 1) args[1] else
  "/data/users/lland/annotation_of_eukaryotic_genome/genespace_workingDirectory"

#focal genome
focal_genome <- if (length(args) >= 2) args[2] else "ISTISU1"

stopifnot(dir.exists(wd))
pg_rds <- file.path(wd, "pangenome_matrix.rds")
stopifnot(file.exists(pg_rds))

#------ LOAD PANGENOME ------
#GENESPACE pangenome matrix: each row = orthogroup (pgID),
#columns for genomes are list-columns of gene IDs
pangenome <- readRDS(pg_rds)

#genome columns are the list-cols (one list per genome)
genome_cols <- names(pangenome)[sapply(pangenome, is.list)]
if (!(focal_genome %in% genome_cols)) {
  stop("focal_genome '", focal_genome, "' not found. Available: ",
       paste(genome_cols, collapse = ", "))
}

pg <- as_tibble(pangenome)

#----- cleaning: drop out-of-synteny genes --------
#helper: turn a list entry into a clean character vector, remove NAs & those ending in '*'
clean_gene_list <- function(v) {
  if (is.null(v)) return(character(0))
  v <- as.character(v)
  v <- trimws(v[!is.na(v)])
  unique(v[!grepl("\\*$", v)])     # drop out-of-synteny genes (GENESPACE marks them with '*')
}

#apply cleaning to all genome list-columns
pg <- pg %>%
  mutate(across(all_of(genome_cols), ~ lapply(.x, clean_gene_list)))

#drop orthogroups that are empty in ALL genomes after cleaning
pg <- pg %>%
  rowwise() %>%
  filter(sum(sapply(c_across(all_of(genome_cols)), length)) > 0) %>%
  ungroup()

n_genomes <- length(genome_cols)

#------ ORTHOGROUP: PRESENCE/ ABSENCE ------
#presence_tbl: for each orthogroup (pgID) and each genome: TRUE if ≥1 gene in that OG
presence_tbl <- pg %>%
  transmute(pgID, across(all_of(genome_cols), ~ lengths(.x) > 0))

#----- CORE, ACCESSORY, SPECIFIC CATEGORIES -----
pg_flags <- presence_tbl %>%
  mutate(
    n_present = select(., all_of(genome_cols)) %>% rowSums(),
    category  = case_when(
      n_present == n_genomes ~ "core",             # present in all genomes
      n_present == 1         ~ "species_specific", # present in exactly one genome
      TRUE                   ~ "accessory"         # present in 2..(n_genomes-1)
    )
  )

#----- GENE COUNTS PER ORTHOGROUP and GENOME -----
# count unique gene IDs in each OG × genome cell
count_genes <- function(g) {
  if (is.null(g)) return(0L)
  length(unique(as.character(g)))
}

gene_counts_matrix <- pg %>%
  select(pgID, all_of(genome_cols)) %>%
  mutate(across(all_of(genome_cols), ~ sapply(.x, count_genes)))

#----- PER-GENOME TOTALS BY CATEGORY -----
#attach OG category to gene counts
gene_counts_w_cat <- pg_flags %>%
  select(pgID, category) %>%
  left_join(gene_counts_matrix, by = "pgID")

#sum gene counts per genome × category
gene_by_cat <- gene_counts_w_cat %>%
  pivot_longer(cols = all_of(genome_cols),
               names_to = "genome",
               values_to = "gene_count") %>%
  group_by(genome, category) %>%
  summarise(gene_count = sum(gene_count), .groups = "drop")

#total number genes per genome (sum over all categories)
gene_totals <- gene_by_cat %>%
  group_by(genome) %>%
  summarise(gene_total = sum(gene_count), .groups = "drop")

#wide table: gene_core, gene_accessory, gene_specific and percentages
gene_counts_per_genome <- gene_by_cat %>%
  pivot_wider(names_from = category,
              values_from = gene_count,
              values_fill = 0) %>%
  rename(
    gene_core      = core,
    gene_accessory = accessory,
    gene_specific  = species_specific
  ) %>%
  left_join(gene_totals, by = "genome") %>%
  mutate(
    percent_core     = round(100 * gene_core     / pmax(gene_total, 1), 2),
    percent_specific = round(100 * gene_specific / pmax(gene_total, 1), 2)
  )

#------ FREQUENCY: ORTHOGROUPS VS GENES BY NUMBER OF GENOMES  ------
#og_freq: how many OGs are present in exactly k genomes (k = 1..n)
og_freq <- pg_flags %>%
  count(n_present, name = "count") %>%
  mutate(type = "Orthogroups")

#gene_freq: for each gene, attach its OG's n_present and count genes per n_present
all_genes_with_presence <- pg %>%
  select(pgID, all_of(genome_cols)) %>%
  pivot_longer(cols = all_of(genome_cols),
               names_to = "genome",
               values_to = "genes_list") %>%
  #drop empty list entries
  filter(!sapply(genes_list, is.null) & sapply(genes_list, length) > 0) %>%
  unnest_longer(genes_list) %>%
  rename(gene = genes_list) %>%
  mutate(gene = as.character(gene)) %>%
  distinct(pgID, genome, gene) %>%       # prevent double-counting the same gene
  left_join(select(pg_flags, pgID, n_present), by = "pgID")

gene_freq <- all_genes_with_presence %>%
  distinct(pgID, gene, n_present) %>%    # one row per unique gene
  count(n_present, name = "count") %>%
  mutate(type = "Genes")

#combined frequency table: OG counts + gene counts as a function of n_present
freq_data <- bind_rows(og_freq, gene_freq)

#----- PLOT: PANGENOME FREQUENCY -----
p <- ggplot(freq_data, aes(x = n_present, y = count, fill = type)) +
  geom_col(position = "dodge", alpha = 0.8, width = 0.7) +
  scale_x_continuous(breaks = 1:n_genomes, labels = 1:n_genomes) +
  scale_y_continuous(labels = scales::comma) +
  scale_fill_manual(values = c("Genes" = "#0072B2", "Orthogroups" = "#D55E00")) +
  labs(
    x = "Number of genomes",
    y = "Count",
    fill = NULL,
    title = "Pangenome composition: distribution across genomes",
    subtitle = paste("Total genomes:", n_genomes)
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position    = "top",
    panel.grid.minor   = element_blank(),
    plot.title         = element_text(face = "bold")
  )

ggsave(file.path(wd, "pangenome_frequency_plot.pdf"),
       p, width = 10, height = 6)

#----- SUMMARY -----
#summary of counts per n_present for both OGs and genes
freq_summary <- freq_data %>%
  pivot_wider(names_from = type,
              values_from = count,
              values_fill = 0) %>%
  mutate(
    label = case_when(
      n_present == n_genomes ~ "Core (all genomes)",
      n_present == 1         ~ "Species-specific (1 genome)",
      TRUE                   ~ paste0("Shared (", n_present, " genomes)")
    )
  ) %>%
  arrange(desc(n_present)) %>%
  select(n_present, label, Orthogroups, Genes)

write_csv(gene_counts_per_genome,
          file.path(wd, "gene_counts_per_genome.csv"))

write_csv(freq_summary,
          file.path(wd, "pangenome_frequency_summary.csv"))

write_csv(presence_tbl,
          file.path(wd, "pangenome_presence_matrix.csv"))

## ----------- CONSOLE PEEK -----------
cat("Genomes:", paste(genome_cols, collapse = ", "), "\n")
cat("Focal:", focal_genome, "\n")
print(head(freq_summary, 10))
