context("Gene importances")

test_that("With generated data", {
  dataset <- generate_dataset(type = "poly", num_genes = 400, num_samples = 101, num_groups = 4)
  expression <- dataset$expression %>% scale_quantile(0)

  time <- seq(-1, 1, length.out = nrow(expression))

  amount_noise <- runif(ncol(expression))
  for (i in seq_len(ncol(expression))) {
    pct <- amount_noise[[i]]
    expression[,i] <- pct * expression[,i] + (1-pct) * runif(nrow(expression))
  }

  gimp1 <- gene_importances(expression, time, ntree = 5) %>% slice(match(colnames(expression), gene))
  gimp2 <- gene_importances(expression, time, ntree = 300) %>% slice(match(colnames(expression), gene))
  gimp3 <- gene_importances(expression, time, ntree = 10000) %>% slice(match(colnames(expression), gene))

  cor_gimp1 <- cor(gimp1$importance, amount_noise)
  cor_gimp2 <- cor(gimp2$importance, amount_noise)
  cor_gimp3 <- cor(gimp3$importance, amount_noise)
  expect_lt(cor_gimp1, cor_gimp2 + .05) # allow for some error
  expect_lt(cor_gimp2, cor_gimp3 + .05) # allow for some error
  expect_gt(cor_gimp3, .35)

  gimp <- gene_importances(expression, time, num_permutations = 5) %>% slice(match(colnames(expression), gene))

  cor_gimp <- cor(gimp$importance, amount_noise)
  expect_gt(cor_gimp, .35)

  cor_gimp <- cor(-gimp$pvalue, amount_noise)
  expect_gt(cor_gimp, .35)
})
