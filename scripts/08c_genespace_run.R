#!/usr/bin/env Rscript
# Usage:
#   Rscript genespace_run.R /path/to/genespace_workingDirectory diamond orthofinder /usr/local/bin

suppressPackageStartupMessages(library(GENESPACE))  # container provides this old-API build

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Provide the GENESPACE workingDirectory as arg 1.")
wd <- normalizePath(args[1], mustWork = TRUE)

# Optional: diamond & orthofinder args (kept for completeness; old API finds them on PATH)
diamond.path     <- if (length(args) >= 2) args[2] else "diamond"
orthofinder.path <- if (length(args) >= 3) args[3] else "orthofinder"

# IMPORTANT: pass the **directory** that contains MCScanX and MCScanX_h
path2mcscanx <- if (length(args) >= 4) args[4] else "/usr/local/bin"

message("[INFO] wd: ", wd)
message("[INFO] path2mcscanx (dir): ", path2mcscanx)

# Discover inputs
pepf <- list.files(file.path(wd, "peptide"),
                   pattern="\\.(fa(sta)?|faa)$", full.names=TRUE, ignore.case=TRUE)
bedf <- list.files(file.path(wd, "bed"),
                   pattern="\\.bed$", full.names=TRUE, ignore.case=TRUE)
stopifnot(length(pepf) > 0, length(bedf) > 0)

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

refGenome <- if ("TAIR10" %in% sp$genome) "TAIR10" else sp$genome[1]
message("[INFO] Genomes: ", paste(sp$genome, collapse=", "))
message("[INFO] Ref genome: ", refGenome)

# Init (old API)
gpar <- init_genespace(
  wd = wd,
  path2mcscanx = path2mcscanx   # <-- DIRECTORY with MCScanX & MCScanX_h
)

# Threads from Slurm
thr <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "8"))
gpar$diamond$threads     <- thr
gpar$orthofinder$threads <- thr

# Old API: single call runs pipeline; do not overwrite completed steps
out <- run_genespace(gpar, overwrite = FALSE)

# Pangenome matrix
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
