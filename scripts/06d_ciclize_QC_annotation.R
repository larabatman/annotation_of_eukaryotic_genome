#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(data.table); library(circlize) })

# ---- hard-coded paths ----
fai   <- "/data/users/lland/annotation_of_eukaryotic_genome/assembly/hifiasm/istisu1.bp.p_ctg.fa.fai"
te_gff<- "/data/users/lland/annotation_of_eukaryotic_genome/EDTA_annotation/istisu1.bp.p_ctg.fa.mod.EDTA.TEanno.gff3"
gene_gff <- "/data/users/lland/annotation_of_eukaryotic_genome/annotation/final/filtered.genes.renamed.gff3"
outd  <- "/data/users/lland/annotation_of_eukaryotic_genome/EDTA_annotation/circlize_all"

# ---- knobs ----
n_scaff  <- 10L       # top N scaffolds
win_size <- 100000L   # window size for densities

dir.create(outd, showWarnings = FALSE, recursive = TRUE)

# -------- helpers ----------
read_fai <- function(p) fread(p, header=FALSE)[, .(chr=V1, len=as.numeric(V2))]
read_gff <- function(p){
  dt <- fread(cmd=sprintf("awk '!/^#/' %s", shQuote(p)),
              sep="\t", header=FALSE, quote="")
  setnames(dt, c("seqid","source","type","start","end","score","strand","phase","attr"))
  dt[]
}
get_field <- function(attr, key) {
  m <- regexpr(sprintf("(^|;)%s=([^;]+)", key), attr, perl=TRUE, ignore.case=TRUE)
  x <- ifelse(m>0, sub(sprintf(".*%s=([^;]+).*", key), "\\1", regmatches(attr,m)), NA_character_)
  trimws(x)
}
# classification looks like "LTR/Gypsy/..." or "LTR/Copia/..."
map_super <- function(classification, type, family){
  parts <- tstrsplit(ifelse(is.na(classification),"",classification), "/", fixed=TRUE)
  a <- toupper(ifelse(is.na(parts[[1]]),"",parts[[1]]))
  b <- toupper(ifelse(is.na(parts[[2]]),"",parts[[2]]))
  fam <- toupper(ifelse(is.na(family),"",family))
  # prefer explicit Copia/Gypsy, else infer from type/family
  out <- ifelse(a=="LTR" & b=="GYPSY","Gypsy",
         ifelse(a=="LTR" & b=="COPIA","Copia",
         ifelse(grepl("GYPSY", fam), "Gypsy",
         ifelse(grepl("COPIA", fam), "Copia",
         ifelse(a=="LTR", "LTR-Other", NA_character_)))))
  # as a last resort: type == LTR_retrotransposon but no subfamily
  out[is.na(out) & toupper(type)=="LTR_RETROTRANSPOSON"] <- "LTR-Other"
  out
}
sanitize_tracks <- function(lst){
  if (length(lst)==0) return(lst)
  lst <- lst[!is.na(names(lst))]
  keep <- vapply(lst, function(d) is.data.frame(d) && nrow(d)>0, logical(1))
  lst[keep]
}
draw_density <- function(seqlens_df, sf_list, cols, title_txt, outfile_prefix, win_size){
  sf_list <- sanitize_tracks(sf_list)
  if (!length(sf_list)) { message("[skip] ", title_txt, " : no intervals"); return(invisible()) }

  max_radial   <- 0.85
  track_height <- min(0.08, max_radial / length(sf_list))
  cols <- cols[names(sf_list)]; cols[is.na(cols)] <- "#999999"
  ascii_title <- sprintf("%s (window=%s bp)", gsub("[^ -~]", "-", title_txt), format(win_size, big.mark=","))

  for (fmt in c("pdf","png")){
    fn <- file.path(outd, sprintf("%s_%dw.%s", outfile_prefix, win_size %/% 1000, fmt))
    if (fmt == "pdf") pdf(fn, 10, 10) else png(fn, 1400, 1400, res = 150)

    # left: circos
    par(fig = c(0, 0.83, 0, 1), mar = c(1, 1, 4, 1), xpd = NA)
    circos.clear()
    circos.par(gap.degree = 2, start.degree = 90,
               track.height = track_height,
               track.margin = c(0.001, 0.001),
               cell.padding = c(0, 0, 0, 0))
    circos.genomicInitialize(seqlens_df, plotType = c("axis", "labels"))

    for(nm in names(sf_list)){
      circos.genomicDensity(sf_list[[nm]],
                            window.size  = win_size,
                            col          = cols[[nm]],
                            track.height = track_height)
    }
    title(main = ascii_title)

    # right: legend
    par(fig = c(0.83, 1, 0, 1), mar = c(1, 1, 4, 2), new = TRUE, xpd = NA)
    plot.new()
    labs <- names(sf_list)
    nc   <- if (length(labs) > 12) 2 else 1
    legend("left", legend = labs, fill = cols[labs], border = NA,
           bty = "n", cex = 0.9, ncol = nc, xpd = NA)

    dev.off()
    message("Saved: ", fn)
  }
}

# -------- load genome/scaffolds --------
fai_dt <- read_fai(fai); setorder(fai_dt, -len)
fai_top <- fai_dt[1:min(n_scaff, .N)]
seqlens <- data.frame(chr=fai_top$chr, start=0, end=fai_top$len)

# -------- load TEs (EDTA) --------
te <- read_gff(te_gff)
te <- te[seqid %in% fai_top$chr]
te[, classification := get_field(attr, "classification")]
te[, family := get_field(attr, "family")]
te[, superfamily := map_super(classification, type, family)]

# keep only Copia & Gypsy
te_cg <- te[superfamily %in% c("Copia","Gypsy")]
sf_te <- te_cg[, .(chr=seqid, start, end, superfamily)]
sf_list_te <- split(sf_te[, .(chr,start,end)], sf_te$superfamily)
sf_list_te <- sanitize_tracks(sf_list_te)

# -------- load genes --------
g <- read_gff(gene_gff)
g <- g[type == "gene" & seqid %in% fai_top$chr]
genes_df <- g[, .(chr=seqid, start, end)]
sf_list_genes <- list(Genes = genes_df)

# -------- colors --------
cols_te <- c(Copia="#1b9e77", Gypsy="#d95f02")
cols_genes <- c(Genes="#023858")
cols_combined <- c(cols_genes, cols_te)

# -------- plots --------
# 1) Genes only
draw_density(seqlens, sf_list_genes, cols_genes,
             sprintf("Gene density — top %d scaffolds", nrow(seqlens)),
             sprintf("GENE_density_top%dscaf", nrow(seqlens)),
             win_size)

# 2) LTR-TEs (Copia & Gypsy) only
draw_density(seqlens, sf_list_te, cols_te,
             sprintf("TE density (LTR Copia/Gypsy) — top %d scaffolds", nrow(seqlens)),
             sprintf("TE_density_Copia_Gypsy_top%dscaf", nrow(seqlens)),
             win_size)

# 3) Combined: Genes + Copia + Gypsy
sf_list_combined <- c(sf_list_genes, sf_list_te)
draw_density(seqlens, sf_list_combined, cols_combined,
             sprintf("Genes + LTR-TEs (Copia/Gypsy) — top %d scaffolds", nrow(seqlens)),
             sprintf("GENE_TE_density_combined_top%dscaf", nrow(seqlens)),
             win_size)

message("[OK] Wrote plots to: ", outd)
