#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(adegenet)
})

args <- commandArgs(trailingOnly = TRUE)

genind_obj <- read.genepop(args[1])
popmap <- read.csv(args[2])
outdir <- args[3]

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ---- VALIDATION ----
gen_names <- indNames(genind_obj)

if ("Sample" %in% names(popmap)) {
  if (!identical(gen_names, popmap$Sample)) {
    stop("Sample names do not match between genepop and popmap.")
  }
}

pop(genind_obj) <- popmap$Population

# ---- CROSS VALIDATION ----
xval <- xvalDapc(
  tab(genind_obj, NA.method="mean"),
  pop(genind_obj),
  n.pca.max=100,
  training.set=0.9,
  result="groupMean",
  n.rep=30,
  xval.plot=FALSE
)

best_n_pca <- xval$best.n.pca

dapc_res <- dapc(genind_obj, n.pca=best_n_pca, n.da=2)

png(file.path(outdir, "dapc.png"), width=2000, height=1500, res=300)
scatter(dapc_res)
dev.off()

saveRDS(dapc_res, file=file.path(outdir, "dapc.rds"))
