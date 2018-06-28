## ---- include=FALSE, cache=FALSE-----------------------------------------
library("knitr")
knitr::opts_knit$set(self.contained = FALSE)
knitr::opts_chunk$set(tidy = TRUE, collapse=TRUE, comment = "#>",
                      tidy.opts=list(blank=FALSE, width.cutoff=55))

## ------------------------------------------------------------------------
X <- as.matrix(NSW[, 2:10])
## weight matrix
Ahalf <- diag(c(0.15, 0.6, 2.5, 2.5, 2.5, 0.5, 0.5, 0.1, 0.1))
d <- NSW$treated
y <- NSW$re78

## Compute matrix of distances between treated and control units
D0 <- distMat(X, d, Ahalf, method="manhattan", FullMatrix=FALSE)

## ------------------------------------------------------------------------
Avar <- chol(solve(cov(X)))
DMvar <- distMat(X, d, Avar, method="euclidean",
                 FullMatrix=TRUE)
sigma2 <- nnvar(DMvar, d, y, J=30)

