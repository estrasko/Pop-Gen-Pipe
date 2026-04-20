#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(adegenet)
  library(ggplot2)
})

# Capture the command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (!(length(args) %in% c(3, 4))) {
  stop("Usage: Rscript run_dapc.R <multi_snp_genepop> <popmap_file> <outdir> [pop_colors_csv]")
}

input_file <- args[1]
popmap_file <- args[2]
outdir <- args[3]
pop_colors_file <- if (length(args) == 4) args[4] else NA

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

# Optional color support
custom_colors <- NULL
ind_colors <- NULL
grp_colors <- NULL

if (!is.na(pop_colors_file)) {
  cat("Reading population colors...\n")
  color_df <- read.csv(pop_colors_file, stringsAsFactors = FALSE)

  if (!all(c("Population", "Color") %in% names(color_df))) {
    stop("Population colors CSV must contain columns named 'Population' and 'Color'.")
  }

  if (anyDuplicated(color_df$Population)) {
    stop("Population colors CSV contains duplicate Population entries.")
  }

  pop_levels <- levels(pop(genind_obj))
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

  custom_colors <- color_lookup
  grp_colors <- unname(custom_colors[pop_levels])
  ind_colors <- unname(custom_colors[as.character(pop(genind_obj))])

  cat("Using custom population colors:\n")
  print(data.frame(Population = pop_levels, Color = grp_colors))
}

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

# Find the optimal number of PCs for DAPC analysis
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

# Specify the number of DA axes to keep
n_groups <- nPop(genind_obj)
best_n_da <- min(2, n_groups - 1)

if (best_n_da < 1) {
  stop("DAPC requires at least two populations.")
}

cat(paste0("Selected n.pca = ", best_n_pca, "\n"))
cat(paste0("Selected n.da = ", best_n_da, "\n"))

cat("Running DAPC...\n")
dapc_result <- dapc(genind_obj, n.pca = best_n_pca, n.da = best_n_da)

# Save results and plots to files
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

# Added code accounts for one discriminant function axis
coords <- as.matrix(dapc_result$ind.coord)
eig <- dapc_result$eig
percent_var <- eig / sum(eig) * 100

# Plotting for 1 axis
if (best_n_da == 1 || ncol(coords) == 1) {

  cat("Generating 1-axis density plot...\n")

  xlab_text <- paste0(
    "Discriminant Axis 1 (",
    round(percent_var[1], 1),
    "%)"
  )

  png(
    filename = file.path(outdir, "dapc_axis1_density.png"),
    width = 2200,
    height = 1800,
    res = 300
  )

  plot(
    density(coords[, 1]),
    main = xlab_text,
    xlab = "DAPC score",
    ylab = "Density",
    lwd = 2
  )

  dev.off()

# Plotting for 2 axes
} else {

  cat("Generating 2-axis scatter plot...\n")

  xlab_text <- paste0(
    "Discriminant Axis 1 (",
    round(percent_var[1], 1),
    "%)"
  )

  ylab_text <- paste0(
    "Discriminant Axis 2 (",
    round(percent_var[2], 1),
    "%)"
  )

  dapc_plot_df <- data.frame(
    LD1 = coords[, 1],
    LD2 = coords[, 2],
    Population = as.factor(pop(genind_obj))
  )

  p <- ggplot(dapc_plot_df, aes(x = LD1, y = LD2, color = Population)) +
    geom_point(size = 6, alpha = 0.8) +
    theme_bw() +
    labs(
      x = xlab_text,
      y = ylab_text,
      title = "DAPC Scatter Plot"
    )

  if (!is.null(custom_colors)) {
    p <- p + scale_color_manual(values = custom_colors)
  }

  ggsave(
    filename = file.path(outdir, "dapc_scatter.png"),
    plot = p,
    width = 8,
    height = 6,
    dpi = 300
  )
}

png(
  filename = file.path(outdir, "dapc_assignplot.png"),
  width = 2400,
  height = 1600,
  res = 300
)

assignplot(dapc_result)

dev.off()

cat("DAPC complete.\n")
