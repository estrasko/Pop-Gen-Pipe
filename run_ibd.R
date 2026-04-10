#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(vegan)
  library(ecodist)
})

args <- commandArgs(trailingOnly = TRUE)

fst <- as.dist(read.csv(args[1], header=FALSE))
geo <- as.dist(read.csv(args[2], header=FALSE))
outdir <- args[3]

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

mantel_res <- mantel.randtest(geo, fst)
mrm_res <- MRM(fst ~ geo)

png(file.path(outdir, "ibd.png"))
plot(geo, fst)
dev.off()
