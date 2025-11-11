######################################
# BUSCO comparison figures (reads short_summary*.txt)
# Keeps the original plotting block style; only data-loading is generalized
######################################

library(ggplot2)
library(grid)

# --- Hardcoded inputs ---
hifiasm_sum <- "/data/users/lland/annotation_of_eukaryotic_genome/assembly/hifiasm/BUSCO/short_summary.specific.brassicales_odb10.hifiasm_busco.txt"
trinity_sum <- "/data/users/lland/annotation_of_eukaryotic_genome/assembly/Trinity/BUSCO/short_summary.specific.brassicales_odb10.trinity_busco.txt"
prot_sum    <- "/data/users/lland/annotation_of_eukaryotic_genome/annotation/BUSCO/busco_proteins/short_summary.specific.brassicales_odb10.busco_proteins.txt"
tran_sum    <- "/data/users/lland/annotation_of_eukaryotic_genome/annotation/BUSCO/busco_transcripts/short_summary.specific.brassicales_odb10.busco_transcripts.txt"

# --- Outputs ---
out_dir <- "/data/users/lland/annotation_of_eukaryotic_genome/annotation/BUSCO"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_genome <- file.path(out_dir, "busco_figure_genome.png")
out_tx     <- file.path(out_dir, "busco_figure_transcripts.png")
out_all    <- file.path(out_dir, "busco_figure_all.png")

# --- Plot config (same spirit as template) ---
my_colors <- c("#56B4E9", "#3492C7", "#F0E442", "#F04442")  # S, D, F, M
my_bar_height <- 0.75
my_title <- "BUSCO Assessment Results"
my_family <- "sans"
my_size_ratio <- 1
my_width <- 20; my_height <- 15; my_unit <- "cm"

# --- Parser for BUSCO short summary (percent or counts) ---
parse_busco <- function(path, label){
  x <- readLines(path, warn = FALSE)
  line <- x[grepl("^\\s*C:", x)]
  if (!length(line)) stop("No 'C:' line in: ", path, call. = FALSE)

  # Percent form
  rx_pct <- "C:\\s*([0-9.]+)%\\s*\\[\\s*S:\\s*([0-9.]+)%,\\s*D:\\s*([0-9.]+)%\\s*\\]\\s*,\\s*F:\\s*([0-9.]+)%\\s*,\\s*M:\\s*([0-9.]+)%\\s*,\\s*n:\\s*(\\d+)"
  m <- regexec(rx_pct, line); mm <- regmatches(line, m)[[1]]
  if (length(mm)) {
    n <- as.integer(mm[7])
    S <- round(n * as.numeric(mm[3]) / 100)
    D <- round(n * as.numeric(mm[4]) / 100)
    F <- round(n * as.numeric(mm[5]) / 100)
    M <- n - (S + D + F)  # keep totals consistent
    return(data.frame(species=label, S=S, D=D, F=F, M=M, n=n, check.names=FALSE))
  }

  # Counts form
  rx_cnt <- "C:\\s*(\\d+)\\s*\\[\\s*S:\\s*(\\d+)\\s*,\\s*D:\\s*(\\d+)\\s*\\]\\s*,\\s*F:\\s*(\\d+)\\s*,\\s*M:\\s*(\\d+)\\s*,\\s*n:\\s*(\\d+)"
  m <- regexec(rx_cnt, line); mm <- regmatches(line, m)[[1]]
  if (length(mm)) {
    S <- as.integer(mm[3]); D <- as.integer(mm[4]); F <- as.integer(mm[5]); M <- as.integer(mm[6]); n <- as.integer(mm[7])
    return(data.frame(species=label, S=S, D=D, F=F, M=M, n=n, check.names=FALSE))
  }

  stop("Could not parse BUSCO summary in: ", path, call. = FALSE)
}

# --- Template-like plotting helper (expects a 'busco' data.frame: species,S,D,F,M,n) ---
plot_busco <- function(busco, outfile, title_suffix=""){
  my_output <- outfile
  dir.create(dirname(my_output), recursive = TRUE, showWarnings = FALSE)

  # Expand into long format (as in the template)
  my_species <- rep(busco$species, each=4)
  my_values <- c(rbind(busco$S, busco$D, busco$F, busco$M))
  my_percentage <- with(busco, c(rbind(S/n*100, D/n*100, F/n*100, M/n*100)))
  my_species <- factor(my_species)
  my_species <- factor(my_species, levels(my_species)[c(length(levels(my_species)):1)]) # reverse order

  labsize <- if (length(levels(my_species)) > 10) 0.66 else 1
  category <- factor(rep(c("S","D","F","M"), c(1)))
  category <- factor(category, levels(category)[c(4,1,2,3)])  # order M,S,D,F (bottom->top)
  df <- data.frame(my_species, my_percentage, my_values, category)

  figure <- ggplot() +
    geom_bar(aes(y = my_percentage, x = my_species, fill = category),
             position = position_stack(reverse = TRUE),
             data = df, stat="identity", width=my_bar_height) +
    coord_flip() +
    theme_gray(base_size = 8) +
    scale_y_continuous(labels = c("0","20","40","60","80","100"), breaks = c(0,20,40,60,80,100)) +
    scale_fill_manual(values = my_colors,labels =c(" Complete (C) and single-copy (S)  ",
                                                   " Complete (C) and duplicated (D)",
                                                   " Fragmented (F)  ",
                                                   " Missing (M)")) +
    ggtitle(paste0(my_title, if(nchar(title_suffix)) paste0(" — ", title_suffix) else "")) +
    xlab("") + ylab("\n%BUSCOs") +
    theme(plot.title = element_text(family=my_family, hjust=0.5, colour = "black",
                                    size = rel(2.2)*my_size_ratio, face = "bold")) +
    theme(legend.position="top",legend.title = element_blank()) +
    theme(legend.text = element_text(family=my_family, size = rel(1.2)*my_size_ratio)) +
    theme(panel.background = element_rect(color="#FFFFFF", fill="white")) +
    theme(panel.grid.minor = element_blank()) +
    theme(panel.grid.major = element_blank()) +
    theme(axis.text.y = element_text(family=my_family, colour = "black", size = rel(1.66)*my_size_ratio)) +
    theme(axis.text.x = element_text(family=my_family, colour = "black", size = rel(1.66)*my_size_ratio)) +
    theme(axis.line = element_line(size=1*my_size_ratio, colour = "black")) +
    theme(axis.ticks.length = unit(.85, "cm")) +
    theme(axis.ticks.y = element_line(colour="white", size = 0)) +
    theme(axis.ticks.x = element_line(colour="#222222")) +
    theme(axis.ticks.length = unit(0.4, "cm")) +
    theme(axis.title.x = element_text(family=my_family, size=rel(1.2)*my_size_ratio)) +
    guides(fill = guide_legend(override.aes = list(colour = NULL))) +
    guides(fill=guide_legend(nrow=2,byrow=TRUE))

  # Inline text labels per bar
  for(i in rev(seq_along(levels(my_species)))){
    detailed_values <- my_values[my_species==my_species[my_species==levels(my_species)[i]]]
    total_buscos <- sum(detailed_values)
    figure <- figure +
      annotate("text",
               label=paste("C:", detailed_values[1] + detailed_values[2],
                           " [S:", detailed_values[1], ", D:", detailed_values[2],
                           "], F:", detailed_values[3], ", M:", detailed_values[4],
                           ", n:", total_buscos, sep=""),
               y=3, x = i, size = labsize*4*my_size_ratio, colour = "black", hjust=0, family=my_family)
  }

  ggsave(figure, file=my_output, width = my_width, height = my_height, unit = my_unit)
  message("Saved: ", my_output)
}

# --- Build data and plot ---

# Genome comparison
busco_genome <- rbind(
  parse_busco(hifiasm_sum, "assembly (hifiasm)"),
  parse_busco(prot_sum,    "annotation (proteins)")
)
plot_busco(busco_genome, out_genome, "Genome: assembly vs annotation (proteins)")

# Transcriptome comparison
busco_tx <- rbind(
  parse_busco(trinity_sum, "transcriptome (Trinity)"),
  parse_busco(tran_sum,    "annotation (transcripts)")
)
plot_busco(busco_tx, out_tx, "Transcriptome: assembly vs annotation (transcripts)")

# All four together (nice overview)
busco_all <- rbind(
  parse_busco(hifiasm_sum, "assembly (hifiasm)"),
  parse_busco(prot_sum,    "annotation (proteins)"),
  parse_busco(trinity_sum, "transcriptome (Trinity)"),
  parse_busco(tran_sum,    "annotation (transcripts)")
)
plot_busco(busco_all, out_all, "All datasets")
