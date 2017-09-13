context("Evaluation functions")

test_that("evaluate_trajectory", {
  time <- runif(1000)
  groups <- cut(time, breaks = 4, labels = F)
  expect_equal(evaluate_trajectory(time, groups), 1)

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
})

test_that("evaluate_dim_red", {
  space <- matrix(runif(1000), ncol = 2)
  groups <- kmeans(space, 5)$cluster
  goodscore <- evaluate_dim_red(space, groups)
  expect_gt(goodscore, .9)

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
})
