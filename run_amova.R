#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(adegenet)
  library(poppr)
})

args <- commandArgs(trailingOnly = TRUE)

gen <- read.genepop(args[1])
popmap <- read.csv(args[2])
outdir <- args[3]

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

if ("Sample" %in% names(popmap)) {
  if (!identical(indNames(gen), popmap$Sample)) {
    stop("Sample mismatch")
  }
}

genclone <- as.genclone(gen)
strata(genclone) <- data.frame(Population=popmap$Population)

res <- poppr.amova(genclone, ~Population)
rand <- randtest(res)

png(file.path(outdir, "amova.png"))
plot(rand)
dev.off()
