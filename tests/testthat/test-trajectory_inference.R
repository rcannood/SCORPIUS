context("Dimensionality reduction")

test_that("With generated data", {
  dataset <- generate_dataset(num_genes = 500, num_samples = 300, num_groups = 4)
  expression <- dataset$expression

  # testing for 1 dimension
  space <- reduce_dimensionality(expression, "spearman", ndim = 1)
  traj <- infer_trajectory(space)

  expect_is(traj, "list")
  expect_is(traj, "SCORPIUS::trajectory")

  expect_equal(names(traj), c("path", "time"))
  expect_equal( ncol(traj$path), ncol(space) )
  expect_equal( colnames(traj$path), colnames(space) )
  expect_equal( length(traj$time), nrow(expression) )
  expect_equal( names(traj$time), rownames(expression) )

  # testing for 2 dimensions
  space <- reduce_dimensionality(expression, "spearman", ndim=2)
  traj <- infer_trajectory(space)

  expect_is(traj, "list")
  expect_is(traj, "SCORPIUS::trajectory")

  expect_equal(names(traj), c("path", "time"))
  expect_equal( ncol(traj$path), ncol(space) )
  expect_equal( colnames(traj$path), colnames(space) )
  expect_equal( length(traj$time), nrow(expression) )
  expect_equal( names(traj$time), rownames(expression) )

  # reverse_trajectory
  rev_traj <- reverse_trajectory(traj)

  expect_is(rev_traj, "list")
  expect_is(rev_traj, "SCORPIUS::trajectory")

  expect_equal(names(rev_traj), c("path", "time"))
  expect_equal( ncol(rev_traj$path), ncol(space) )
  expect_equal( colnames(rev_traj$path), colnames(space) )
  expect_equal( length(rev_traj$time), nrow(expression) )
  expect_equal( names(rev_traj$time), rownames(expression) )
})


test_that("With generated data and edge case", {
  dataset <- generate_dataset(num_genes = 500, num_samples = 10, num_groups = 4)
  expression <- dataset$expression

  # testing for 1 dimension
  space <- reduce_dimensionality(expression, "spearman", ndim=1)
  traj <- infer_trajectory(space)

  expect_is(traj, "list")
  expect_is(traj, "SCORPIUS::trajectory")

  expect_equal(names(traj), c("path", "time"))
  expect_equal( ncol(traj$path), ncol(space) )
  expect_equal( colnames(traj$path), colnames(space) )
  expect_equal( length(traj$time), nrow(expression) )
  expect_equal( names(traj$time), rownames(expression) )

  # testing for 2 dimensions
  space <- reduce_dimensionality(expression, "spearman", ndim = 2)
  traj <- infer_trajectory(space)

  expect_is(traj, "list")
  expect_is(traj, "SCORPIUS::trajectory")

  expect_equal(names(traj), c("path", "time"))
  expect_equal( ncol(traj$path), ncol(space) )
  expect_equal( colnames(traj$path), colnames(space) )
  expect_equal( length(traj$time), nrow(expression) )
  expect_equal( names(traj$time), rownames(expression) )

  # reverse_trajectory
  rev_traj <- reverse_trajectory(traj)

  expect_is(rev_traj, "list")
  expect_is(rev_traj, "SCORPIUS::trajectory")

  expect_equal(names(rev_traj), c("path", "time"))
  expect_equal( ncol(rev_traj$path), ncol(space) )
  expect_equal( colnames(rev_traj$path), colnames(space) )
  expect_equal( length(rev_traj$time), nrow(expression) )
  expect_equal( names(rev_traj$time), rownames(expression) )
})
