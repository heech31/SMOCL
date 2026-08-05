test_that("smocl validates its inputs", {
  expect_error(
    smocl(
      x = matrix("a", nrow = 3, ncol = 1),
      k = 2,
      lambda = 0.5
    ),
    "`x` must be a numeric matrix or data frame"
  )

  expect_error(
    smocl(
      x = matrix(1, nrow = 3, ncol = 1),
      k = 2,
      lambda = 0.5
    ),
    "constant columns"
  )

  x <- matrix(rnorm(20), nrow = 10, ncol = 2)

  expect_error(
    smocl(x = x, k = 1, lambda = 0.5),
    "`k` must be an integer"
  )

  expect_error(
    smocl(x = x, k = 2, lambda = 0),
    "`lambda` must contain positive finite values"
  )
})

test_that("smocl returns a valid fitted object", {
  baseball_path <- system.file(
    "extdata",
    "baseball.rds",
    package = "SMOCL"
  )

  expect_true(nzchar(baseball_path))

  baseball <- readRDS(baseball_path)

  set.seed(20260726)

  fit <- smocl(
    x = baseball$x,
    k = 4,
    lambda = c(0.1, 0.5, 1),
    max_iter = 100,
    n_basis = 5,
    degree = 2,
    tol = 1e-6
  )

  expect_s3_class(fit, "smocl")
  expect_length(fit$cluster, nrow(baseball$x))
  expect_equal(dim(fit$fx), dim(baseball$x))
  expect_equal(fit$k, 4L)
  expect_true(fit$lam_wcss %in% c(0.1, 0.5, 1))
  expect_equal(sum(table(fit$cluster)), nrow(baseball$x))
})

test_that("print.smocl prints a concise fit summary", {
  fit <- structure(
    list(
      cluster = c(1L, 1L, 2L),
      fx = matrix(0, nrow = 3, ncol = 2),
      k = 2L,
      lam_wcss = 0.5
    ),
    class = "smocl"
  )

  output <- capture.output(returned <- print(fit))

  expect_identical(returned, fit)
  expect_true(any(grepl("Sparse Monotone Clustering fit", output)))
  expect_true(any(grepl("Number of observations: 3", output)))
  expect_true(any(grepl("Number of clusters: 2", output)))
  expect_true(any(grepl("Selected lambda: 0.5", output)))
})
