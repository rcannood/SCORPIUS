context("Testing evaluation.R")

test_that("evaluate_trajectory", {
  time <- runif(1000)
  groups <- cut(time, breaks = 4, labels = F)
  expect_equal(evaluate_trajectory(time, groups), 1)
  expect_equal(evaluate_trajectory(time, rank(time)), 1)

  ix <- which(groups %in% c(2,3))
  six <- sample(ix, length(ix)/2)
  groups[six] <- sample(groups[six])
  badscore <- evaluate_trajectory(time, groups)
  expect_lt(badscore, 1)

  ix <- seq_along(time)
  six <- sample(ix, length(ix)*3/4)
  groups[six] <- sample(groups[six])
  worsescore <- evaluate_trajectory(time, groups)
  expect_lt(worsescore, badscore)

  expect_error(evaluate_trajectory(list(), groups), "must be a numeric vector")
  expect_error(evaluate_trajectory(time, list()), "must be a numeric vector or a factor")
  expect_error(evaluate_trajectory(time, c(groups, groups)), "must have equal lengths")
  expect_error(evaluate_trajectory(c(time, time), groups), "must have equal lengths")
})

test_that("evaluate_dim_red", {
  space <- matrix(runif(1000), ncol = 2)
  groups <- kmeans(space, 5)$cluster
  goodscore <- evaluate_dim_red(space, groups)
  expect_gt(goodscore, .9)
  expect_gt(evaluate_dim_red(space, as.integer(groups)), .9)

  ix <- seq_along(groups)
  six <- sample(ix, length(ix)*.2)
  groups[six] <- sample(groups[six])
  badscore <- evaluate_dim_red(space, groups)
  expect_lt(badscore, goodscore)

  ix <- seq_along(groups)
  six <- sample(ix, length(ix)*.75)
  groups[six] <- sample(groups[six])
  worsescore <- evaluate_dim_red(space, groups)
  expect_lt(worsescore, badscore)

  expect_error(evaluate_dim_red(list(), groups), "must be a numeric matrix or, a data frame")
  expect_error(evaluate_dim_red(space, list()), "must be a numeric vector or a factor")
  expect_error(evaluate_dim_red(space, c(groups, groups)), "nrow.*length.*must be the same")
  expect_error(evaluate_dim_red(rbind(space, space), groups), "nrow.*length.*must be the same")
})
