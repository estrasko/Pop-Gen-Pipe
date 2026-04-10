#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(vegan)
  library(ecodist)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

if (!(length(args) %in% c(3, 4))) {
  stop("Usage: Rscript run_ibd.R <fst_csv> <distance_csv> <outdir> [summary_stats_csv]")
}

fst_file <- args[1]
distance_file <- args[2]
outdir <- args[3]
summary_stats_file <- if (length(args) == 4) args[4] else NA

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

cat("Reading Fst matrix...\n")
fst_df <- read.csv(fst_file, header = FALSE, check.names = FALSE)
fst_dist <- quasieuclid(as.dist(as.matrix(fst_df)))

cat("Reading geographic distance matrix...\n")
distance_df <- read.csv(distance_file, header = FALSE, check.names = FALSE)
geo_dist <- quasieuclid(as.dist(as.matrix(distance_df)))

cat("Running Mantel test...\n")
mantel_result <- mantel.randtest(geo_dist, fst_dist, nrepet = 1000)

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

png(
  filename = file.path(outdir, "mantel_randtest_plot.png"),
  width = 2200,
  height = 1800,
  res = 300
)
plot(mantel_result)
dev.off()

pdf(file.path(outdir, "fst_vs_distance.pdf"), width = 7, height = 4.5)
plot(geo_dist, fst_dist)
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
