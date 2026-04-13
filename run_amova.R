#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(adegenet)
  library(poppr)
})

# Capture the command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (!(length(args) %in% c(3, 4))) {
  stop("Missing arguments. Expected 3 or 4, got ", length(args),
     ".\nUsage: Rscript run_amova.R <haps_genepop> <popmap_file> <outdir> [pop_colors_csv]")
}

input_file <- args[1]
popmap_file <- args[2]
outdir <- args[3]
pop_colors_file <- if (length(args) == 4) args[4] else NA

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

# If a population color file is provided, validate it and save a copy with the AMOVA output.
# The AMOVA randtest plot is not population-specific, so the colors are documented rather than plotted.
if (!is.na(pop_colors_file)) {
  cat("Reading population colors...\n")
  color_df <- read.csv(pop_colors_file, stringsAsFactors = FALSE)

  if (!all(c("Population", "Color") %in% names(color_df))) {
    stop("Population colors CSV must contain columns named 'Population' and 'Color'.")
  }

  if (anyDuplicated(color_df$Population)) {
    stop("Population colors CSV contains duplicate Population entries.")
  }

  pop_levels <- sort(unique(as.character(popmap_df$Population)))
  color_lookup <- setNames(color_df$Color, color_df$Population)

  missing_colors <- setdiff(pop_levels, names(color_lookup))
  if (length(missing_colors) > 0) {
    stop(
      paste0(
        "Missing colors for populations: ",
        paste(missing_colors, collapse = ", ")
      )
    )
  }

  write.csv(
    color_df,
    file = file.path(outdir, "amova_population_colors.csv"),
    row.names = FALSE
  )

  writeLines(
    c(
      "Population colors were validated successfully.",
      "The AMOVA randtest plot is not population-specific, so colors were documented but not applied to the plot."
    ),
    con = file.path(outdir, "amova_color_info.txt")
  )
}

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
