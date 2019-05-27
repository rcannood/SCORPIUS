context("Testing extract_modules.R")

dataset <- generate_dataset(num_genes = 101, num_samples = 91, num_groups = 4)
expression <- dataset$expression
time <- seq(-1, 1, length.out = nrow(expression))

check_modules <- function(modules, expression) {
  expect_is( modules, "data.frame" )
  expect_equal( colnames(modules), c("feature", "orig_index", "module", "within_module_ordering") )
  expect_true( all(colnames(expression) %in% modules$feature) )
  expect_equal( nrow(modules), ncol(expression) )
}

test_that("With generated data", {
  modules <- extract_modules(expression)

  check_modules(modules, expression)

  modules <- extract_modules(expression, time)

  check_modules(modules, expression)
})

dataset <- generate_dataset(num_genes = 10, num_samples = 100, num_groups = 4)
expression <- dataset$expression
time <- seq(-1, 1, length.out = nrow(expression))

test_that("With generated data with edge case", {
  modules <- extract_modules(expression)

  check_modules(modules, expression)

  modules <- extract_modules(expression, time)

  check_modules(modules, expression)
})

test_that("it fails gracefully", {
  expect_error(extract_modules(list(), time), "must be a numeric matrix")
  expect_error(extract_modules(expression, list()), "must be a numeric vector")

  expect_warning(extract_modules(rbind(runif(1001), runif(1001))), "has more than 1000 features")
})
