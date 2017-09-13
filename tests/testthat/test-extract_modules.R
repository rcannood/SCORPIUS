context("Extract modules")

test_that("With generated data", {
  dataset <- generate_dataset(type = "poly", num_genes = 101, num_samples = 91, num_groups = 4)
  expression <- dataset$expression

  time <- seq(-1, 1, length.out = nrow(expression))

  modules <- extract_modules(expression)

  expect_is( modules, "data.frame" )
  expect_equal( colnames(modules), c("feature", "orig_index", "module", "within_module_ordering") )
  expect_true( all(colnames(expression) %in% modules$feature) )
  expect_equal( nrow(modules), ncol(expression) )

  modules <- extract_modules(expression, time)

  expect_is( modules, "data.frame" )
  expect_equal( colnames(modules), c("feature", "orig_index", "module", "within_module_ordering") )
  expect_true( all(colnames(expression) %in% modules$feature) )
  expect_equal( nrow(modules), ncol(expression) )
})

test_that("With generated data with edge case", {
  dataset <- generate_dataset(type = "poly", num_genes = 10, num_samples = 100, num_groups = 4)
  expression <- dataset$expression

  time <- seq(-1, 1, length.out = nrow(expression))

  modules <- extract_modules(expression)

  expect_is( modules, "data.frame" )
  expect_equal( colnames(modules), c("feature", "orig_index", "module", "within_module_ordering") )
  expect_true( all(colnames(expression) %in% modules$feature) )
  expect_equal( nrow(modules), ncol(expression) )

  modules <- extract_modules(expression, time)

  expect_is( modules, "data.frame" )
  expect_equal( colnames(modules), c("feature", "orig_index", "module", "within_module_ordering") )
  expect_true( all(colnames(expression) %in% modules$feature) )
  expect_equal( nrow(modules), ncol(expression) )
})
