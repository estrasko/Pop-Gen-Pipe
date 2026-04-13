#!/usr/bin/env Rscript

#Import packages without printing package messages
suppressPackageStartupMessages({
  library(diveRsity)
  library(qgraph)
  library(reshape2)
})

#Reads command line arguments from Pop_Script_2.py into R
args <- commandArgs(trailingOnly = TRUE)

if (!(length(args) %in% c(5, 6, 7))) {
  stop("Usage: Rscript run_divmigrate.R <multi_snp_genepop> <outdir> <stat> <boots> <threads> [node_names_csv] [pop_colors_csv]")
}

input_file <- args[1]
outdir <- args[2]
stat <- args[3]
boots <- as.integer(args[4])
threads <- as.integer(args[5])
node_names_arg <- if (length(args) >= 6 && args[6] != "__NONE__") args[6] else NA
pop_colors_file <- if (length(args) == 7) args[7] else NA

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

# If a population color file is provided, validate it and use the first color
# as a single default node color for the whole divMigrate network.
node_color <- "lightblue"

if (!is.na(pop_colors_file)) {
  cat("Reading population colors...\n")
  color_df <- read.csv(pop_colors_file, stringsAsFactors = FALSE)

  if (!all(c("Population", "Color") %in% names(color_df))) {
    stop("Population colors CSV must contain columns named 'Population' and 'Color'.")
  }

  if (nrow(color_df) < 1) {
    stop("Population colors CSV must contain at least one row.")
  }

  node_color <- color_df$Color[1]

  write.csv(
    color_df,
    file = file.path(outdir, "divmigrate_population_colors.csv"),
    row.names = FALSE
  )
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
  curve = 0.8,
  color = "white",
  edge.color = node_color,        
  label.color = "black"
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
  curve = 0.8,
  color = "white",
  edge.color = node_color,
  label.color = "black"
)
dev.off()

capture.output(
  results,
  file = file.path(outdir, "divmigrate_summary.txt")
)

saveRDS(results, file = file.path(outdir, "divmigrate_result.rds"))

cat("divMigrate complete.\n")
