#' Fit Sparse Monotone Clustering
#'
#' Fits the Sparse MOnotone CLustering (SMOCL) model to identify
#' ordinally interpretable clusters and variables associated with the
#' estimated cluster ordering.
#'
#' @param x A numeric matrix or data frame. Rows are observations and
#'   columns are variables.
#' @param k A positive integer giving the number of clusters.
#' @param lambda A positive numeric vector of candidate penalty parameters.
#' @param max_iter Maximum number of alternating updates for each candidate
#'   value of `lambda`.
#' @param n_basis Number of spline basis functions for each variable.
#' @param degree Degree of the I-spline basis.
#' @param tol Convergence tolerance.
#'
#' @return An object of class `"smocl"` containing the ordered cluster
#'   assignments, estimated monotone functions, cooperative-lasso
#'   coefficients, selected penalty parameter, and tuning paths.
#'
#' @examples
#' baseball_path <- system.file(
#'   "extdata",
#'   "baseball.rds",
#'   package = "SMOCL"
#' )
#'
#' baseball <- readRDS(baseball_path)
#'
#' set.seed(20260726)
#'
#' fit <- smocl(
#'   x = baseball$x,
#'   k = 4,
#'   lambda = c(0.1, 0.5, 1),
#'   max_iter = 100,
#'   n_basis = 5,
#'   degree = 2
#' )
#'
#' print(fit)
#'
#' @export
smocl <- function(
    x,
    k,
    lambda = seq(0.1, 2, length.out = 50),
    max_iter = 500,
    n_basis = 3,
    degree = 2,
    tol = 1e-6) {

  if (is.data.frame(x)) {
    x <- data.matrix(x)
  }

  if (!is.matrix(x) || !is.numeric(x)) {
    stop("`x` must be a numeric matrix or data frame.", call. = FALSE)
  }

  if (nrow(x) < 3L || ncol(x) < 1L) {
    stop("`x` must contain at least three rows and one column.", call. = FALSE)
  }

  if (anyNA(x) || any(!is.finite(x))) {
    stop("`x` must not contain missing or non-finite values.", call. = FALSE)
  }

  column_sd <- apply(x, 2L, stats::sd)

  if (any(!is.finite(column_sd)) || any(column_sd == 0)) {
    stop("`x` must not contain constant columns.", call. = FALSE)
  }

  if (
    length(k) != 1L ||
      !is.numeric(k) ||
      !is.finite(k) ||
      k != as.integer(k) ||
      k < 2L ||
      k >= nrow(x)
  ) {
    stop(
      "`k` must be an integer between 2 and `nrow(x) - 1`.",
      call. = FALSE
    )
  }

  if (
    !is.numeric(lambda) ||
      length(lambda) < 1L ||
      any(!is.finite(lambda)) ||
      any(lambda <= 0)
  ) {
    stop("`lambda` must contain positive finite values.", call. = FALSE)
  }

  if (
    length(max_iter) != 1L ||
      max_iter != as.integer(max_iter) ||
      max_iter < 1L
  ) {
    stop("`max_iter` must be a positive integer.", call. = FALSE)
  }

  if (
    length(n_basis) != 1L ||
      n_basis != as.integer(n_basis) ||
      n_basis < 1L
  ) {
    stop("`n_basis` must be a positive integer.", call. = FALSE)
  }

  if (
    length(degree) != 1L ||
      degree != as.integer(degree) ||
      degree < 1L
  ) {
    stop("`degree` must be a positive integer.", call. = FALSE)
  }

  if (
    length(tol) != 1L ||
      !is.numeric(tol) ||
      !is.finite(tol) ||
      tol <= 0
  ) {
    stop("`tol` must be a positive finite number.", call. = FALSE)
  }

  fit <- MOCL_wcss_scale(
    Xf = x,
    K = as.integer(k),
    max.iter = as.integer(max_iter),
    df = as.integer(n_basis),
    degree = as.integer(degree),
    lambda_seq = lambda,
    tol = tol
  )

  fit$call <- match.call()
  fit$k <- as.integer(k)
  fit$variable_names <- colnames(x)

  class(fit) <- "smocl"
  fit
}
