#!/usr/bin/env Rscript

# generate dataset with certain seed
set.seed(1)
data <- dyntoy::generate_dataset(
  id = "specific_example/scorpius",
  num_cells = 99,
  num_features = 101,
  model = "linear",
  normalise = FALSE
)

# add method specific args (if needed)
data$parameters <- list()

data$seed <- 1L

# write example dataset to file
file <- commandArgs(trailingOnly = TRUE)
if (length(file) > 0) {
  file <- file[[1]]
  dynutils::write_h5(data, file)
}

data
