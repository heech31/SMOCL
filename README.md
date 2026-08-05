# SMOCL: Sparse Monotone Clustering with Ordinal Interpretation

`SMOCL` provides an R implementation of Sparse MOnotone CLustering
with ordinal interpretation. The method estimates ordered clusters,
monotone variable transformations, and a sparse set of variables
associated with the cluster ordering.

## Installation

The current package-development branch can be installed from GitHub
using `remotes`.

```r
# install.packages("remotes")

remotes::install_github(
  "Mos-start2092/SMOCL",
  ref = "package-dev-scoop"
)
```

The cooperative-lasso implementation required by SMOCL is included in
the package. The archived `scoop` package does not need to be installed
separately.

## Quick start


the package. The archived `scoop` package does not need to be installed
Load the package and the Baseball data included with SMOCL.

```r
library(SMOCL)

baseball_path <- system.file(
  "extdata",
  "baseball.rds",
  package = "SMOCL"
)

baseball <- readRDS(baseball_path)
mlb_x <- baseball$x
```

Fit SMOCL with four clusters and three candidate penalty parameters.

```r
set.seed(20260726)

fit <- smocl(
  x = mlb_x,
  k = 4,
  lambda = c(0.1, 0.5, 1),
  max_iter = 100,
  n_basis = 5,
  degree = 2,
  tol = 1e-6
)

fit
```

The fitted object contains the ordered cluster assignments, estimated
monotone functions, cooperative-lasso coefficients, selected penalty
parameter, and tuning paths.

```r
fit$cluster
fit$lam_wcss
fit$coef
fit$fx
```

Variables with nonzero estimated monotone functions can be identified
as follows.

```r
selected <- colSums(abs(fit$fx)) > 0

colnames(mlb_x)[selected]
```

## Cluster visualization

The following optional example visualizes the fitted clusters using
t-SNE. It requires the `Rtsne` package.

```r
# install.packages("Rtsne")

set.seed(252)

tsne <- Rtsne::Rtsne(
  mlb_x,
  dims = 2,
  perplexity = 30,
  max_iter = 500
)

plot(
  tsne$Y[, 1],
  tsne$Y[, 2],
  type = "n",
  xlab = "Y1",
  ylab = "Y2"
)

text(
  tsne$Y[, 1],
  tsne$Y[, 2],
  labels = fit$cluster,
  col = fit$cluster,
  cex = 1.2
)
```

![](README_files/figure-gfm/tsne-1.png)

## Estimated monotone functions

The estimated monotone functions for the selected variables can be
displayed using base R graphics.

```r
group <- fit$cluster
fx <- fit$fx
selected_index <- which(colSums(abs(fx)) > 0)
variable_names <- colnames(mlb_x)

if (length(selected_index) > 0L) {
  old_par <- par(no.readonly = TRUE)
  par(mfrow = c(1, length(selected_index)))

  for (j in selected_index) {
    mean_x <- tapply(mlb_x[, j], group, mean)
    mean_fx <- tapply(fx[, j], group, mean)

    plot(
      mlb_x[, j],
      fx[, j],
      col = group,
      pch = group,
      cex = 0.8,
      xlab = variable_names[j],
      ylab = "Estimated monotone function",
      main = variable_names[j]
    )

    points(
      mean_x,
      mean_fx,
      pch = 20,
      cex = 2
    )
  }

  par(old_par)
}
```

![](README_files/figure-gfm/functions-1.png)

## Legacy research scripts

The original research functions are retained in the `Functions/`
directory for reproducibility. New analyses should ordinarily use the
exported `smocl()` interface.

## Reference

Chung, H. C., An, S., Cha, H., and Ahn, J.  
“Monotone Clusters with Ordinal Interpretation.”  
*Statistical Analysis and Data Mining* (to appear).
