#' Print a SMOCL Fit
#'
#' Prints a concise summary of a fitted SMOCL object.
#'
#' @param x An object of class `"smocl"`.
#' @param ... Additional arguments. Currently unused.
#'
#' @return The input object, invisibly.
#'
#' @method print smocl
#' @export
print.smocl <- function(x, ...) {
  cluster_sizes <- table(x$cluster)

  cat("Sparse Monotone Clustering fit\n")
  cat("--------------------------------\n")
  cat("Number of observations:", length(x$cluster), "\n")
  cat("Number of variables:", ncol(x$fx), "\n")
  cat("Number of clusters:", x$k, "\n")
  cat("Selected lambda:", x$lam_wcss, "\n")

  cat("Cluster sizes:\n")
  print(cluster_sizes)

  invisible(x)
}
