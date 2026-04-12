#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(adegenet)
  library(ggplot2)
})

# Capture the command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Missing arguments. Expected 3, got ", length(args),
     ".\nUsage: Rscript run_dapc.R <multi_snp_genepop> <popmap_file> <outdir>")
}

input_file <- args[1]
popmap_file <- args[2]
outdir <- args[3]

# Create directory to store output files
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Create a temporary .gen copy of the input file for use with read.genepop()
# adegenet::read.genepop() expects a .gen extension
temp_gen_file <- tempfile(fileext = ".gen")
file.copy(input_file, temp_gen_file, overwrite = TRUE)

cat("Reading multi-SNP genepop file...\n")
genind_obj <- read.genepop(temp_gen_file)

cat("Reading popmap...\n")
popmap_df <- read.csv(popmap_file, stringsAsFactors = FALSE)

if (!("Population" %in% names(popmap_df))) {
  stop("Popmap file must contain a column named 'Population'.")
}

# Extract individual sample names from the genind object 
genepop_samples <- indNames(genind_obj)

# Check that sample names in the popmap file match the names from the genind object
if (!is.null(genepop_samples) && "Sample" %in% names(popmap_df)) {
  popmap_samples <- as.character(popmap_df$Sample)

  if (!setequal(genepop_samples, popmap_samples)) {
    missing_in_popmap <- setdiff(genepop_samples, popmap_samples)
    missing_in_genepop <- setdiff(popmap_samples, genepop_samples)

    stop(
      paste0(
        "Sample mismatch between genepop and popmap.\n",
        "Missing in popmap: ", paste(missing_in_popmap, collapse = ", "), "\n",
        "Missing in genepop: ", paste(missing_in_genepop, collapse = ", ")
      )
    )
  }

  popmap_df <- popmap_df[match(genepop_samples, popmap_samples), ]
} else {
  if (nInd(genind_obj) != nrow(popmap_df)) {
    stop(
      paste0(
        "Mismatch between number of individuals in genepop (", nInd(genind_obj),
        ") and number of rows in popmap (", nrow(popmap_df), ")."
      )
    )
  }

  warning(
    "Popmap does not contain a 'Sample' column. Validation is limited to row count only."
  )
}

# Label samples in genind object by population using the popmap file
pop(genind_obj) <- as.factor(popmap_df$Population)

# Determine maximum number of PCA components based on your dataset
max_pca_allowed <- min(
  100,
  nInd(genind_obj) - 1,
  ncol(tab(genind_obj, NA.method = "mean"))
)

if (max_pca_allowed < 2) {
  stop("Not enough data to run DAPC cross-validation.")
}

cat("Running DAPC cross-validation...\n")

# Set random seed to ensure reproducible results across runs
set.seed(123)

xval_obj <- xvalDapc(
  x = tab(genind_obj, NA.method = "mean"),
  grp = pop(genind_obj),
  n.pca.max = max_pca_allowed,
  training.set = 0.9,
  result = "groupMean",
  center = TRUE,
  scale = FALSE,
  n.rep = 30,
  xval.plot = FALSE
)

best_n_pca <- NULL

if (!is.null(xval_obj$`Number of PCs Achieving Highest Mean Success`)) {
  best_n_pca <- xval_obj$`Number of PCs Achieving Highest Mean Success`
} else if (!is.null(xval_obj$best.n.pca)) {
  best_n_pca <- xval_obj$best.n.pca
}

if (is.null(best_n_pca)) {
  stop("Could not determine best n.pca from xvalDapc output.")
}

best_n_pca <- as.integer(best_n_pca)
n_groups <- nPop(genind_obj)
best_n_da <- min(2, n_groups - 1)

if (best_n_da < 1) {
  stop("DAPC requires at least two populations.")
}

cat(paste0("Selected n.pca = ", best_n_pca, "\n"))
cat(paste0("Selected n.da = ", best_n_da, "\n"))

dapc_result <- dapc(genind_obj, n.pca = best_n_pca, n.da = best_n_da)

capture.output(
  list(
    selected_n_pca = best_n_pca,
    selected_n_da = best_n_da,
    dapc_summary = dapc_result
  ),
  file = file.path(outdir, "dapc_summary.txt")
)

saveRDS(xval_obj, file = file.path(outdir, "dapc_xval_result.rds"))
saveRDS(dapc_result, file = file.path(outdir, "dapc_result.rds"))

png(
  filename = file.path(outdir, "dapc_scatter.png"),
  width = 2200,
  height = 1800,
  res = 300
)
scatter(
  dapc_result,
  scree.da = TRUE,
  posi.da = "bottomleft",
  cell = 0,
  cstar = 0,
  clab = 0
)
dev.off()

png(
  filename = file.path(outdir, "dapc_assignplot.png"),
  width = 2400,
  height = 1600,
  res = 300
)
assignplot(dapc_result)
dev.off()

cat("DAPC complete.\n")
