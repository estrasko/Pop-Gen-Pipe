#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(diveRsity)
  library(qgraph)
})

args <- commandArgs(trailingOnly = TRUE)

input <- args[1]
outdir <- args[2]
stat <- args[3]
boots <- as.integer(args[4])

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

res <- divMigrate(input, stat=stat, boots=boots, plot_network=FALSE)

mat <- res$nmRelMig
diag(mat) <- NA

png(file.path(outdir, "divmigrate.png"), width=2000, height=1500, res=300)
qgraph(mat)
dev.off()
