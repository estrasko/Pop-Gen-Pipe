#!/usr/bin/env Rscript

# Load packages
suppressPackageStartupMessages({
  library(vegan)
  library(ecodist)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

# Check number of arguments
if (!(length(args) %in% c(3, 4))) {
  stop("Usage: Rscript run_ibd.R <fst_csv> <distance_csv> <outdir> [summary_stats_csv]")
}

# Assign inputs
fst_file <- args[1]
distance_file <- args[2]
outdir <- args[3]
summary_stats_file <- if (length(args) == 4) args[4] else NA

# Create output directory if it doesn't exist
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Read Fst matrix and geographic distance matrix
cat("Reading Fst matrix...\n")
fst_df <- read.csv(fst_file, header = FALSE, check.names = FALSE)
fst_mat <- as.matrix(fst_df)

cat("Reading geographic distance matrix...\n")
distance_df <- read.csv(distance_file, header = FALSE, check.names = FALSE)
geo_mat <- as.matrix(distance_df)

#Checking to make sure there's enough populations to run IBD
n_pop <- nrow(fst_mat)

if (n_pop <3) {
   message("Not enough populations for IBD analysis (need >=3. Skipping Mantel and MRM analyses)")

  writeLines(
    paste("Fst:", fst_mat[1,2], "Distance:", geo_mat[1,2]),
    con = file.path(outdir, "ibd_skipped.txt")
  )
  
  quit(save = "no", status = 0)
}


# Basic validation
if (nrow(fst_mat) != ncol(fst_mat)) {
  stop("fst.csv must be a square matrix.")
}

if (nrow(geo_mat) != ncol(geo_mat)) {
  stop("geo.csv must be a square matrix.")
}

if (!all(dim(fst_mat) == dim(geo_mat))) {
  stop("fst.csv and geo.csv must have the same dimensions.")
}

fst_dist <- as.dist(fst_mat)
geo_dist <- as.dist(geo_mat)

cat("Running Mantel test...\n")
mantel_result <- vegan::mantel(
  geo_dist,
  fst_dist,
  method = "pearson",
  permutations = 1000
)

cat("Running MRM...\n")
mrm_result <- MRM(fst_dist ~ geo_dist, method = "linear", mrank = FALSE)

capture.output(
  mantel_result,
  file = file.path(outdir, "mantel_result.txt")
)

capture.output(
  mrm_result,
  file = file.path(outdir, "mrm_result.txt")
)

saveRDS(mantel_result, file = file.path(outdir, "mantel_result.rds"))
saveRDS(mrm_result, file = file.path(outdir, "mrm_result.rds"))

# Simple scatter plot of geographic distance vs Fst
pdf(file.path(outdir, "fst_vs_distance.pdf"), width = 7, height = 4.5)
plot(
  x = as.vector(geo_dist),
  y = as.vector(fst_dist),
  xlab = "Geographic distance",
  ylab = "Fst"
)
abline(lm(as.vector(fst_dist) ~ as.vector(geo_dist)))
dev.off()

png(
  filename = file.path(outdir, "fst_vs_distance.png"),
  width = 2200,
  height = 1800,
  res = 300
)
plot(
  x = as.vector(geo_dist),
  y = as.vector(fst_dist),
  xlab = "Geographic distance",
  ylab = "Fst"
)
abline(lm(as.vector(fst_dist) ~ as.vector(geo_dist)))
dev.off()

if (!is.na(summary_stats_file)) {
  cat("Reading summary statistics file...\n")
  stats_df <- read.csv(summary_stats_file, stringsAsFactors = FALSE)

  required_cols <- c("Distance")
  if (!all(required_cols %in% names(stats_df))) {
    stop("Summary stats CSV must contain a 'Distance' column if provided.")
  }

  metric_map <- list(
    Ho = "observed_heterozygosity_vs_distance.png",
    He = "expected_heterozygosity_vs_distance.png",
    Pi = "pi_vs_distance.png",
    AR = "allelic_richness_vs_distance.png"
  )

  regression_summaries <- list()

  for (metric in names(metric_map)) {
    if (metric %in% names(stats_df)) {
      model <- lm(stats_df[[metric]] ~ stats_df$Distance)
      regression_summaries[[metric]] <- summary(model)

      p <- ggplot(stats_df, aes(x = Distance, y = .data[[metric]])) +
        geom_point(shape = 1) +
        geom_smooth(method = "lm", se = FALSE) +
        expand_limits(y = 0) +
        theme_bw()

      ggsave(
        filename = file.path(outdir, metric_map[[metric]]),
        plot = p,
        width = 7,
        height = 4.5,
        dpi = 300
      )
    }
  }

  capture.output(
    regression_summaries,
    file = file.path(outdir, "distance_regression_summaries.txt")
  )
}

cat("IBD complete.\n")
