#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(diveRsity)
  library(qgraph)
  library(reshape2)
})

args <- commandArgs(trailingOnly = TRUE)

if (!(length(args) %in% c(4, 5))) {
  stop("Usage: Rscript run_divmigrate.R <multi_snp_genepop> <outdir> <stat> <boots> [node_names_csv]")
}

input_file <- args[1]
outdir <- args[2]
stat <- args[3]
boots <- as.integer(args[4])
node_names_arg <- if (length(args) == 5) args[5] else NA

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

cat("Running divMigrate...\n")
results <- divMigrate(
  infile = input_file,
  stat = stat,
  plot_network = FALSE,
  boots = boots
)

if (!is.null(results$nmRelMig)) {
  write.csv(results$nmRelMig, file.path(outdir, "divmigrate_nmRelMig.csv"), row.names = TRUE)
}

if (!is.null(results$nmRelMigSig)) {
  write.csv(results$nmRelMigSig, file.path(outdir, "divmigrate_nmRelMigSig.csv"), row.names = TRUE)
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
par(mar = c(6, 10, 6, 6), xpd = NA)
qgraph(
  gmat,
  nodeNames = if (!is.null(colnames(gmat))) colnames(gmat) else NULL,
  legend = TRUE,
  edge.labels = TRUE,
  curve = 2.5,
  mar = c(2, 10, 6, 6)
)
dev.off()

png(
  filename = file.path(outdir, "divmigrate_network.png"),
  width = 2400,
  height = 1800,
  res = 300
)
par(mar = c(6, 10, 6, 6), xpd = NA)
qgraph(
  gmat,
  nodeNames = if (!is.null(colnames(gmat))) colnames(gmat) else NULL,
  legend = TRUE,
  edge.labels = TRUE,
  curve = 2.5,
  mar = c(2, 10, 6, 6)
)
dev.off()

capture.output(
  results,
  file = file.path(outdir, "divmigrate_summary.txt")
)

saveRDS(results, file = file.path(outdir, "divmigrate_result.rds"))

cat("divMigrate complete.\n")
