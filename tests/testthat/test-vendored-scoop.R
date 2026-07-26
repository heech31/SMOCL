test_that("vendored cooperative lasso reproduces the reference fit", {

  set.seed(20260726)

  n <- 80L

  Z <- matrix(
    rnorm(n * 6L),
    nrow = n,
    ncol = 6L
  )

  Z <- scale(
    Z,
    center = TRUE,
    scale = TRUE
  )

  groups <- rep(
    1:2,
    each = 3L
  )

  beta_true <- c(
    0.8, 0.5, 0.2,
    -0.6, -0.4, -0.1
  )

  y <- 0.3 +
    drop(Z %*% beta_true) +
    rnorm(n, sd = 0.5)

  coop_lasso_internal <- getFromNamespace(
    "coop.lasso",
    "SMOCL"
  )

  fit <- coop_lasso_internal(
    x = Z,
    y = y,
    group = groups,
    family = "gaussian",
    intercept = TRUE,
    normalize = FALSE,
    wk = sqrt(tabulate(groups)),
    lambda = 5,
    eps = 1e-8
  )

  estimated <- as.numeric(
    methods::slot(
      fit,
      "coefficients"
    )
  )

  expected <- c(
     0.264763491,
     0.666765007,
     0.334062659,
     0.163646441,
    -0.326364985,
    -0.369945207,
    -0.006391538
  )

  expect_s4_class(
    fit,
    "scoopfit"
  )

  expect_length(
    estimated,
    7L
  )

  expect_equal(
    estimated,
    expected,
    tolerance = 1e-7
  )
})
