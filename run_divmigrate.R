#!/usr/bin/env Rscript

#Import packages without printing package messages
suppressPackageStartupMessages({
  library(diveRsity)
  library(qgraph)
  library(reshape2)
})

#Reads command line arguments from Pop_Script_2.py into R
args <- commandArgs(trailingOnly = TRUE)

#Check that all required files to run analyses are included in the run code
if (!(length(args) %in% c(5, 6))) {
  stop("Usage: Rscript run_divmigrate.R <multi_snp_genepop> <outdir> <stat> <boots> <threads> [node_names_csv]")
}

input_file <- args[1] #names multi_snp_genepop file as input
outdir <- args[2] #outdir is the included path
stat <- args[3] #String argument
boots <- as.integer(args[4]) #Assign bootstrap value (integer)
threads <- as.integer(args[5]) #number of threads defined in command
node_names_arg <- if (length(args) == 6) args[6] else NA #Make sure number of arguments pass is 6, otherwise report NA

#Create output directory for this script
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

#Check point for threads to be a number greater than 1
if (is.na(threads) || threads < 1) {
  stop("threads must be an integer >= 1")
}


cat("Running divMigrate...\n")
cat(paste0("Statistic: ", stat, "\n"))
cat(paste0("Bootstrap replicates: ", boots, "\n"))
cat(paste0("Threads requested: ", threads, "\n"))

# divMigrate accepts genepop .gen/.txt, so create a temporary .gen copy
temp_gen_file <- tempfile(fileext = ".gen")
file.copy(input_file, temp_gen_file, overwrite = TRUE)

# The diveRsity API only exposes parallel on/off through para
use_parallel <- threads > 1
cat(paste0("Parallel enabled: ", use_parallel, "\n"))

results <- divMigrate(
  infile = temp_gen_file,
  stat = stat,
  plot_network = FALSE,
  boots = boots,
  para = use_parallel
)

if (!is.null(results$nmRelMig)) {
  write.csv(
    results$nmRelMig,
    file.path(outdir, "divmigrate_nmRelMig.csv"),
    row.names = TRUE
  )
}

if (!is.null(results$nmRelMigSig)) {
  write.csv(
    results$nmRelMigSig,
    file.path(outdir, "divmigrate_nmRelMigSig.csv"),
    row.names = TRUE
  )
}

if (!is.null(results$gRelMig)) {
  gmat <- results$gRelMig
} else if (!is.null(results$nmRelMig)) {
  gmat <- results$nmRelMig
} else {
  stop("Could not find a relative migration matrix in the divMigrate result.")
}

diag(gmat) <- NA

if (!is.na(node_names_arg)) {
  labs <- trimws(unlist(strsplit(node_names_arg, ",")))
  if (length(labs) == nrow(gmat)) {
    rownames(gmat) <- colnames(gmat) <- labs
  }
}

write.csv(
  round(gmat, 3),
  file = file.path(outdir, "divmigrate_relative_migration_matrix.csv"),
  na = ""
)

gmat_long <- na.omit(
  melt(gmat, varnames = c("From", "To"), value.name = "RelMigration")
)

write.csv(
  gmat_long,
  file = file.path(outdir, "divmigrate_relative_migration_long.csv"),
  row.names = FALSE
)

pdf(
  file.path(outdir, "divmigrate_network.pdf"),
  width = 8,
  height = 6,
  useDingbats = FALSE
)
par(mar = c(4, 4, 4, 4), xpd = NA)
qgraph(
  gmat,
  labels = if (!is.null(colnames(gmat))) colnames(gmat) else NULL,
  edge.labels = TRUE,
  curve = 0.8
)
dev.off()

png(
  filename = file.path(outdir, "divmigrate_network.png"),
  width = 2400,
  height = 1800,
  res = 300
)
par(mar = c(4, 4, 4, 4), xpd = NA)
qgraph(
  gmat,
  labels = if (!is.null(colnames(gmat))) colnames(gmat) else NULL,
  edge.labels = TRUE,
  curve = 0.8
)
dev.off()

capture.output(
  results,
  file = file.path(outdir, "divmigrate_summary.txt")
)

saveRDS(results, file = file.path(outdir, "divmigrate_result.rds"))

cat("divMigrate complete.\n")
