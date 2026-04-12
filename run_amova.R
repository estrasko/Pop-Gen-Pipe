#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(adegenet)
  library(poppr)
})

# Capture the command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Missing arguments. Expected 3, got ", length(args),
     ".\nUsage: Rscript run_amova.R <haps_genepop> <popmap_file> <outdir>")
}

input_file <- args[1]
popmap_file <- args[2]
outdir <- args[3]

# Create directory to store output files
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Create a temporary .gen copy of the input file for use with read.genepop()
# adegenet::read.genepop() requires a .gen extension
temp_gen_file <- tempfile(fileext = ".gen")
file.copy(input_file, temp_gen_file, overwrite = TRUE)

cat("Reading haplotype genepop file...\n")
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

# Create a genclone object for use with poppr.amova()
genclone_obj <- as.genclone(genind_obj)

# Label samples in genclone object by population using the popmap file
strata(genclone_obj) <- data.frame(Population = popmap_df$Population)

cat("Running AMOVA...\n")
amova_result <- poppr.amova(
  genclone_obj,
  ~Population,
  cutoff = 0.5,
  method = "ade4"
)

cat("Running randomization test...\n")
amova_randtest <- randtest(amova_result, nrepet = 999)

# Save results and randtest plot to files
capture.output(
  amova_result,
  file = file.path(outdir, "amova_result.txt")
)

capture.output(
  amova_randtest,
  file = file.path(outdir, "amova_randtest.txt")
)

saveRDS(amova_result, file = file.path(outdir, "amova_result.rds"))
saveRDS(amova_randtest, file = file.path(outdir, "amova_randtest.rds"))

png(
  filename = file.path(outdir, "amova_randtest_plot.png"),
  width = 2200,
  height = 1800,
  res = 300
)
plot(amova_randtest)
dev.off()

cat("AMOVA complete.\n")
