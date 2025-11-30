#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(circlize)
})

#----- CONFIG -----
fai      <- "/data/users/lland/annotation_of_eukaryotic_genome/assembly/hifiasm/istisu1.bp.p_ctg.fa.fai"
te_gff   <- "/data/users/lland/annotation_of_eukaryotic_genome/EDTA_annotation/istisu1.bp.p_ctg.fa.mod.EDTA.TEanno.gff3"
gene_gff <- "/data/users/lland/annotation_of_eukaryotic_genome/annotation/final/filtered.genes.renamed.gff3"
outd     <- "/data/users/lland/annotation_of_eukaryotic_genome/EDTA_annotation/circlize_all"

#----- PARAMETERS ------
n_scaff      <- 5L        # top N scaffolds (same as other script)
win_size     <- 100000L    # window size for densities
smooth_lines <- TRUE       # filled ribbons + outline

#TE families to show in subset plots
te_subset_cg    <- c("Gypsy")
te_subset_cgtir <- c("Copia", "Gypsy", "TIR")

dir.create(outd, showWarnings = FALSE, recursive = TRUE)

#----- HELPER FUNCITONS ------
read_fai <- function(p)
  fread(p, header = FALSE)[, .(chr = V1, len = as.numeric(V2))]

read_gff <- function(p){
  dt <- fread(cmd = sprintf("awk '!/^#/' %s", shQuote(p)),
              sep = "\t", header = FALSE, quote = "")
  setnames(dt, c("seqid","source","type","start","end",
                 "score","strand","phase","attr"))
  dt[]
}

get_classification <- function(attr){
  m <- regexpr("(^|;)[Cc]lassification=([^;]+)", attr, perl = TRUE)
  x <- ifelse(m > 0,
              sub(".*[Cc]lassification=([^;]+).*", "\\1", regmatches(attr, m)),
              NA_character_)
  gsub("\\s+","", x)
}

#Same mapping as in the TE-only circos script
map_super <- function(cls){
  parts <- tstrsplit(ifelse(is.na(cls), "", cls), "/", fixed = TRUE)
  a <- toupper(ifelse(is.na(parts[[1]]), "", parts[[1]]))
  b <- toupper(ifelse(is.na(parts[[2]]), "", parts[[2]]))
  ifelse(a=="LTR" & b=="COPIA","Copia",
  ifelse(a=="LTR" & b=="GYPSY","Gypsy",
  ifelse(a=="LTR","LTR-Other",
  ifelse(a=="DNA" & b %in% c("DTA","DTC","DTM","DTH","DTE","DTP","DTT","DTX"),"TIR",
  ifelse(a=="DNA" & b=="HELITRON","Helitron",
  ifelse(a=="DNA","DNA-Other",
  ifelse(a=="LINE","LINE",
  ifelse(a=="SINE","SINE","Unknown"))))))))
}

sanitize_tracks <- function(lst){
  if (length(lst)==0) return(lst)
  lst <- lst[!is.na(names(lst))]
  keep <- vapply(lst, function(d) is.data.frame(d) && nrow(d) > 0, logical(1))
  lst[keep]
}

draw_density <- function(seqlens_df, track_list, cols,
                         title_txt, outfile_prefix,
                         win_size, smooth_lines = TRUE) {

  track_list <- sanitize_tracks(track_list)
  if (!length(track_list)) {
    message("[skip] ", title_txt, " : no intervals")
    return(invisible())
  }

  max_radial <- 0.85
  #choose height so all tracks fit within max_radial
  track_height <- max_radial / length(track_list)
  #cap it so single-track plots don't become huge
  if (track_height > 0.12) track_height <- 0.12

  cols <- cols[names(track_list)]
  cols[is.na(cols)] <- "#999999"

  ascii_title <- sprintf(
    "%s (window=%s bp)",
    gsub("[^ -~]", "-", title_txt),
    format(win_size, big.mark = ",")
  )

  n_chr   <- nrow(seqlens_df)
  gap_deg <- min(2, 360 / max(1, n_chr) / 10)

  dens_list <- NULL
  ymax <- 1
  if (smooth_lines) {
    dens_list <- lapply(track_list, function(df) {
      genomicDensity(df, window.size = win_size)
    })
    ymax <- max(unlist(lapply(dens_list, function(df) df[[4]])), na.rm = TRUE)
    if (!is.finite(ymax) || ymax <= 0) ymax <- 1
  }

  for (fmt in c("pdf","png")) {
    fn <- file.path(outd, sprintf("%s_%dw.%s",
                                  outfile_prefix, win_size %/% 1000, fmt))
    if (fmt == "pdf") {
      pdf(fn, 10, 10)
    } else {
      png(fn, 1200, 1200, res = 150)
    }

    #LEFT PANEL: circos
    par(fig = c(0, 0.85, 0, 1), mar = c(1, 1, 4, 1), xpd = NA)
    circos.clear()
    circos.par(
      gap.degree   = gap_deg,
      start.degree = 90,
      track.height = track_height,
      track.margin = c(0.001, 0.001),
      cell.padding = c(0, 0, 0, 0)
    )
    circos.genomicInitialize(seqlens_df, plotType = c("axis", "labels"))

    if (smooth_lines) {
      for (nm in names(dens_list)) {
        d <- dens_list[[nm]]
        col_current <- cols[[nm]]
        #more opaque + thicker line
        fill_col    <- grDevices::adjustcolor(col_current, alpha.f = 0.8)
        line_lwd    <- 2

        circos.genomicTrack(
          d,
          ylim = c(0, ymax),
          track.height = track_height,
          panel.fun = function(region, value, ...) {
            v2 <- cbind(value, bottom = 0)
            circos.genomicRect(
              region, v2,
              ytop.column    = 1,
              ybottom.column = 2,
              col    = fill_col,
              border = NA,
              ...
            )
            circos.genomicLines(region, value, col = col_current, lwd = line_lwd, ...)
          }
        )
      }
    } else {
      for (nm in names(track_list)) {
        circos.genomicDensity(
          track_list[[nm]],
          window.size  = win_size,
          col          = cols[[nm]],
          track.height = track_height
        )
      }
    }

    title(main = ascii_title)

    #RIGHT PANEL: legend
    par(fig = c(0.85, 1, 0, 1), mar = c(1, 1, 4, 2), new = TRUE, xpd = NA)
    plot.new()
    labs <- names(track_list)
    nc   <- if (length(labs) > 12) 2 else 1
    legend(
      "left",
      legend = labs,
      fill   = cols[labs],
      border = NA,
      bty    = "n",
      cex    = 0.9,
      ncol   = nc,
      xpd    = NA
    )

    dev.off()
    message("Saved: ", fn)
  }
}

#----- LOAD GENOME -----
fai_dt  <- read_fai(fai)
setorder(fai_dt, -len)
fai_top <- fai_dt[1:min(n_scaff, .N)]
seqlens <- data.frame(chr = fai_top$chr, start = 0, end = fai_top$len)

#----- LOAD TEs -----
te <- read_gff(te_gff)
te <- te[seqid %in% fai_top$chr]
te[, classification := get_classification(attr)]
te[, superfamily    := map_super(classification)]

te_dt <- te[!is.na(superfamily) & superfamily != ""]
sf_te_all <- te_dt[, .(chr = seqid, start, end, superfamily)]
sf_list_te_all <- split(sf_te_all[, .(chr, start, end)], sf_te_all$superfamily)
sf_list_te_all <- sanitize_tracks(sf_list_te_all)
families_all <- sort(names(sf_list_te_all))
#collapse all TE superfamilies into one track ---
sf_te_merged <- sf_te_all[, .(chr, start, end)]
sf_list_te_merged <- list(TEs = sf_te_merged)
cols_te_merged <- c(TEs = "#e41a1c")

#----- LOAD GENES -----
g <- read_gff(gene_gff)
g <- g[type == "gene" & seqid %in% fai_top$chr]
genes_df <- g[, .(chr = seqid, start, end)]
sf_list_genes <- list(Genes = genes_df)

#colours
pal_te <- c(
  Gypsy      = "#d95f02",
  Copia      = "#1b9e77",
  `DNA-Other`= "#7570b3",
  TIR        = "#e7298a",
  Helitron   = "#66a61e",
  LINE       = "#e6ab02",
  SINE       = "#a6761d",
  `LTR-Other`= "#666666",
  Unknown    = "#999999"
)

te_cols <- pal_te[families_all]
miss <- is.na(te_cols)
if (any(miss)) {
  extras <- families_all[miss]
  te_cols[miss] <- grDevices::rainbow(length(extras))
  names(te_cols) <- families_all
}

cols_genes <- c(Genes = "#023858")

#----- PLOTS -----

#1) Genes only
draw_density(
  seqlens,
  sf_list_genes,
  cols_genes,
  sprintf("Gene density — top %d scaffolds", nrow(seqlens)),
  sprintf("GENE_density_top%dscaf", nrow(seqlens)),
  win_size,
  smooth_lines = smooth_lines
)

#2) Genes vs all TEs merged into one track
sf_list_genes_te <- c(sf_list_genes, sf_list_te_merged)
cols_genes_te    <- c(cols_genes, cols_te_merged)

draw_density(
  seqlens,
  sf_list_genes_te,
  cols_genes_te,
  sprintf("Genes vs all TEs — top %d scaffolds", nrow(seqlens)),
  sprintf("GENE_vs_allTE_top%dscaf", nrow(seqlens)),
  win_size,
  smooth_lines = smooth_lines
)

#3) Genes + ALL TE superfamilies
sf_list_comb_all <- c(sf_list_genes, sf_list_te_all)
cols_comb_all    <- c(cols_genes, te_cols)
draw_density(
  seqlens,
  sf_list_comb_all,
  cols_comb_all,
  sprintf("Genes + all TE superfamilies — top %d scaffolds", nrow(seqlens)),
  sprintf("GENE_TE_ALL_density_top%dscaf", nrow(seqlens)),
  win_size,
  smooth_lines = smooth_lines
)

#4) Genes + Copia + Gypsy
cg <- intersect(te_subset_cg, families_all)
if (length(cg) > 0) {
  sf_list_cg  <- sf_list_te_all[cg]
  cols_cg     <- te_cols[cg]
  sf_list_cg_comb <- c(sf_list_genes, sf_list_cg)
  cols_cg_comb    <- c(cols_genes, cols_cg)

  draw_density(
    seqlens,
    sf_list_cg_comb,
    cols_cg_comb,
    sprintf("Genes + Copia + Gypsy — top %d scaffolds", nrow(seqlens)),
    sprintf("GENE_TE_Copia_Gypsy_top%dscaf", nrow(seqlens)),
    win_size,
    smooth_lines = smooth_lines
  )
} else {
  message("[copia-gypsy] none present in TE annotation")
}

#5) Genes + Copia + Gypsy + TIR
cgt <- intersect(te_subset_cgtir, families_all)
if (length(cgt) > 0) {
  sf_list_cgt   <- sf_list_te_all[cgt]
  cols_cgt      <- te_cols[cgt]
  sf_list_cgt_c <- c(sf_list_genes, sf_list_cgt)
  cols_cgt_c    <- c(cols_genes, cols_cgt)

  draw_density(
    seqlens,
    sf_list_cgt_c,
    cols_cgt_c,
    sprintf("Genes + Copia + Gypsy + TIR — top %d scaffolds", nrow(seqlens)),
    sprintf("GENE_TE_Copia_Gypsy_TIR_top%dscaf", nrow(seqlens)),
    win_size,
    smooth_lines = smooth_lines
  )
} else {
  message("[copia-gypsy-TIR] none present in TE annotation")
}

message("[OK] Wrote plots to: ", outd)
