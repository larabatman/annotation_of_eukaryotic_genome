#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(circlize)
})

#---- CONFIG----
gff  <- "/data/users/lland/annotation_of_eukaryotic_genome/EDTA_annotation/istisu1.bp.p_ctg.fa.mod.EDTA.TEanno.gff3"
fai  <- "/data/users/lland/annotation_of_eukaryotic_genome/assembly/hifiasm/istisu1.bp.p_ctg.fa.fai"
outd <- "/data/users/lland/annotation_of_eukaryotic_genome/EDTA_annotation/circlize_all"

#---- PARAMETERS FOR SEGMENTS ----
n_scaff      <- 10L       # always top 10 scaffolds
win_size     <- 100000L   # density window
smooth_lines <- TRUE      # use smoothed density, not bars

dir.create(outd, showWarnings = FALSE, recursive = TRUE)

#----- READ FILES -----
read_fai <- function(p) {
  fread(p, header = FALSE)[, .(chr = V1, len = as.numeric(V2))]
}

read_gff <- function(p){
  dt <- fread(cmd = sprintf("awk '!/^#/' %s", shQuote(p)),
              sep = "\t", header = FALSE, quote = "")
  setnames(dt, c("seqid","source","type","start","end",
                 "score","strand","phase","attr"))
  dt
}

#----- CLASSIFY TEs-----
get_classification <- function(attr){
  m <- regexpr("(^|;)[Cc]lassification=([^;]+)", attr, perl = TRUE)
  x <- ifelse(m > 0,
              sub(".*[Cc]lassification=([^;]+).*", "\\1", regmatches(attr, m)),
              NA_character_)
  gsub("\\s+","", x)
}

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

draw_density <- function(seqlens_df, sf_list, cols,
                         title_txt, outfile_prefix,
                         win_size, smooth_lines = TRUE) {

  #keep only non-empty data frames
  sf_list <- sf_list[vapply(sf_list, function(d) is.data.frame(d) && nrow(d) > 0, logical(1))]
  if (!length(sf_list)) {
    message("[skip] ", title_txt, " : no intervals")
    return(invisible())
  }

  max_radial   <- 0.85
  track_height <- min(0.08, max_radial / length(sf_list))

  cols <- cols[names(sf_list)]
  cols[is.na(cols)] <- "#999999"

  ascii_title <- sprintf(
    "%s (window=%s bp)",
    gsub("[^ -~]", "-", title_txt),
    format(win_size, big.mark = ",")
  )

  #dynamic gap just to be safe, but with 10 scaffolds it will stay at 2°
  n_chr   <- nrow(seqlens_df)
  gap_deg <- min(2, 360 / max(1, n_chr) / 10)

  #precompute densities if using smooth lines
  dens_list <- NULL
  ymax <- 1
  if (smooth_lines) {
    dens_list <- lapply(sf_list, function(df) {
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

    #---- LEFT PANEL: circos ----
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
      #one track per TE superfamily, filled area + outline
      for (nm in names(dens_list)) {
        d <- dens_list[[nm]]
        col_current <- cols[[nm]]
        fill_col    <- grDevices::adjustcolor(col_current, alpha.f = 0.4)

        circos.genomicTrack(
          d,
          ylim = c(0, ymax),
          track.height = track_height,
          panel.fun = function(region, value, ...) {
            # value: first column = density
            v2 <- cbind(value, bottom = 0)
            # filled area
            circos.genomicRect(
              region, v2,
              ytop.column = 1,
              ybottom.column = 2,
              col = fill_col,
              border = NA,
              ...
            )
            #line on top
            circos.genomicLines(region, value, col = col_current, ...)
          }
        )
      }
    } else {
      #fallback: original bar-style density if the smoothed line does not work
      for (nm in names(sf_list)) {
        circos.genomicDensity(
          sf_list[[nm]],
          window.size  = win_size,
          col          = cols[[nm]],
          track.height = track_height
        )
      }
    }

    title(main = ascii_title)

    #---- RIGHT PANEL: legend ----
    par(fig = c(0.85, 1, 0, 1), mar = c(1, 1, 4, 2), new = TRUE, xpd = NA)
    plot.new()
    labs <- names(sf_list)
    nc   <- if (length(labs) > 10) 2 else 1
    legend(
      "left",
      legend = labs,
      fill   = cols[labs],
      border = NA,
      bty    = "n",
      cex    = 0.85,
      ncol   = nc,
      xpd    = NA
    )

    dev.off()
  }
}

#---- LOAD ANNOTATION ----
fai_dt <- read_fai(fai)
setorder(fai_dt, -len)
fai_top <- fai_dt[1:min(n_scaff, .N)]
seqlens <- data.frame(chr = fai_top$chr, start = 0, end = fai_top$len)

gff_dt <- read_gff(gff)
gff_dt <- gff_dt[seqid %in% fai_top$chr]

gff_dt[, classification := get_classification(attr)]
gff_dt[, superfamily    := map_super(classification)]
gff_dt[, bp             := pmax(0L, end - start + 1L)]

sf_dt <- gff_dt[, .(chr = seqid, start, end, superfamily)]
sf_dt <- sf_dt[!is.na(superfamily) & superfamily != ""]

sf_list_all <- split(sf_dt[, .(chr, start, end)], sf_dt$superfamily)
sf_list_all <- sanitize_tracks(sf_list_all)
families_all <- sort(names(sf_list_all))

pal_fixed <- c(
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

cols <- pal_fixed[families_all]
miss <- is.na(cols)
if (any(miss)) {
  extras <- families_all[miss]
  cols[miss] <- grDevices::rainbow(length(extras))
  names(cols) <- families_all
}

message("[all-superfamilies] families: ", paste(families_all, collapse = ", "))

#1) All superfamilies, 10 scaffolds
draw_density(
  seqlens,
  sf_list_all,
  cols,
  sprintf("TE density - all superfamilies - top %d scaffolds", nrow(seqlens)),
  sprintf("TE_density_ALL_top%dscaf", nrow(seqlens)),
  win_size,
  smooth_lines = smooth_lines
)

#2) Copia vs Gypsy only
cg <- intersect(c("Copia","Gypsy"), families_all)
if (length(cg) >= 1) {
  lst_cg  <- sf_list_all[cg]
  cols_cg <- cols[cg]
  message("[copia-gypsy] families: ", paste(cg, collapse=", "))
  draw_density(
    seqlens,
    lst_cg,
    cols_cg,
    sprintf("TE density - Copia vs Gypsy - top %d scaffolds", nrow(seqlens)),
    sprintf("TE_density_Copia_Gypsy_top%dscaf", nrow(seqlens)),
    win_size,
    smooth_lines = smooth_lines
  )
} else {
  message("[copia-gypsy] families: none")
}

#3) Copia + Gypsy + TIR
cgt <- intersect(c("Copia","Gypsy","TIR"), families_all)
if (length(cgt) >= 1) {
  lst_cgt  <- sf_list_all[cgt]
  cols_cgt <- cols[cgt]
  message("[copia-gypsy-TIR] families: ", paste(cgt, collapse=", "))
  draw_density(
    seqlens,
    lst_cgt,
    cols_cgt,
    sprintf("TE density - Copia, Gypsy, TIR - top %d scaffolds", nrow(seqlens)),
    sprintf("TE_density_Copia_Gypsy_TIR_top%dscaf", nrow(seqlens)),
    win_size,
    smooth_lines = smooth_lines
  )
} else {
  message("[copia-gypsy-TIR] families: none")
}

message("[OK] Wrote plots to: ", outd)
