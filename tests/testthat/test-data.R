context("Dataset")

test_that("Testing ginhoux dataset", {
  expect_error(data(ginhoux), NA)

  expect_is(ginhoux, "list")

  expect_named(ginhoux)

  expect_equal(sort(names(ginhoux)), c("expression", "sample_info"))

  expect_is(ginhoux$expression, c("matrix", "data.frame"))
  expect_is(ginhoux$sample_info, "data.frame")

  expression <- ginhoux$expression
  sample_info <- ginhoux$sample_info

  expect_equal(nrow(expression), nrow(sample_info))

  expect_equal(rownames(sample_info), rownames(expression))

  expect_named(sample_info)
  expect_equal(colnames(sample_info), "group_name")
})

test_that("Testing generate_dataset", {
  expect_error(dataset <- generate_dataset(), NA)

  expect_is(dataset, "list")

  expect_named(dataset)

  expect_equal(sort(names(dataset)), c("expression", "sample_info"))

  expect_is(dataset$expression, c("matrix", "data.frame"))
  expect_is(dataset$sample_info, "data.frame")

  expression <- dataset$expression
  sample_info <- dataset$sample_info

  expect_equal(nrow(expression), nrow(sample_info))

  expect_equal(rownames(sample_info), rownames(expression))

  expect_named(sample_info)
  expect_equal(colnames(sample_info), "group_name")
})
