#!/usr/bin/env Rscript

#------ CONFIG ------
suppressPackageStartupMessages(library(GENESPACE))  # container provides this old-API build

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Provide the GENESPACE workingDirectory as arg 1.")
wd <- normalizePath(args[1], mustWork = TRUE)

diamond.path     <- if (length(args) >= 2) args[2] else "diamond"
orthofinder.path <- if (length(args) >= 3) args[3] else "orthofinder"
path2mcscanx <- if (length(args) >= 4) args[4] else "/usr/local/bin"

message("[INFO] wd: ", wd)
message("[INFO] path2mcscanx (dir): ", path2mcscanx)

#Genomes that are present:
pepf <- list.files(file.path(wd, "peptide"),
                   pattern="\\.(fa(sta)?|faa)$", full.names=TRUE, ignore.case=TRUE)
bedf <- list.files(file.path(wd, "bed"),
                   pattern="\\.bed$", full.names=TRUE, ignore.case=TRUE)
stopifnot(length(pepf) > 0, length(bedf) > 0)
#We need at least two genomes with matching names in both folders 
stem <- function(x) sub("\\.[^.]+$", "", basename(x))
pepNames <- stem(pepf); bedNames <- stem(bedf)
common <- intersect(pepNames, bedNames)
if (length(common) < 2) stop("Need ≥2 genomes with BOTH peptide + bed.")

pick <- function(tbl, namesWanted){
  setNames(vapply(namesWanted, \(nm) normalizePath(tbl[stem(tbl)==nm], TRUE), ""),
           namesWanted)
}
pepMap <- pick(pepf, common)
bedMap <- pick(bedf, common)

sp <- data.frame(
  genome   = common,
  path_pep = unname(pepMap[common]),
  path_bed = unname(bedMap[common]),
  stringsAsFactors = FALSE
)
write.csv(sp, file.path(wd, "species_table.csv"), row.names=FALSE)

#Pick TAIR10 as the reference genome 
refGenome <- if ("TAIR10" %in% sp$genome) "TAIR10" else sp$genome[1]
message("[INFO] Genomes: ", paste(sp$genome, collapse=", "))
message("[INFO] Ref genome: ", refGenome)

#Initialize and run the pipeline 
gpar <- init_genespace(
  wd = wd,
  path2mcscanx = "/usr/local/bin",
  verbose = TRUE
)

# Threads from Slurm
thr <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "8"))
gpar$diamond$threads     <- thr
gpar$orthofinder$threads <- thr

#To prevent overwriting on completed steps
out <- run_genespace(gpar, overwrite = FALSE)

#Build and save the pangenome matrix (presence/absence) in reference coordinates
pg <- query_pangenes(
  out,
  bed = NULL,
  refGenome = refGenome,
  transform = TRUE,
  showArrayMem = TRUE,
  showNSOrtho = TRUE,
  maxMem2Show = Inf
)
saveRDS(out, file.path(wd, "genespace_results.rds"))
saveRDS(pg,  file.path(wd, "pangenome_matrix.rds"))
message("[OK] Saved genespace_results.rds and pangenome_matrix.rds")
