#!/usr/local/bin/Rscript

requireNamespace("dyncli", quietly = TRUE)
task <- dyncli::main()

library(SCORPIUS, warn.conflicts = FALSE)

output <- SCORPIUS:::ti_scorpius_run_fun(
  expression = task$expression,
  priors = task$priors,
  parameters = task$parameters,
  seed = task$seed,
  verbose = task$verbose
)

dyncli::write_output(output, task$output)
