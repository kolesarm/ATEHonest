## ---- include=FALSE, cache=FALSE----------------------------------------------
library("knitr")
knitr::opts_knit$set(self.contained = FALSE)
knitr::opts_chunk$set(tidy = TRUE, collapse=TRUE, comment = "#>",
                      tidy.opts=list(blank=FALSE, width.cutoff=55))

## -----------------------------------------------------------------------------
library("ATEHonest")
X <- as.matrix(NSWexper[, 2:10])
d <- NSWexper$treated
y <- NSWexper$re78
DMvar <- distMat(X, chol(solve(cov(X))), method="euclidean")
sigma2 <- nnvar(DMvar, d, y, J=3)

## -----------------------------------------------------------------------------
Ahalf <- diag(c(0.15, 0.6, 2.5, 2.5, 2.5, 0.5, 0.5, 0.1, 0.1))
D0 <- distMat(X, Ahalf, method="manhattan", d)

## -----------------------------------------------------------------------------
op <- ATTOptPath(y, d, D0, maxsteps=120)
ATTOptEstimate(op, mean(sigma2), C=1, sigma2final=sigma2, opt.criterion="RMSE")
ATTOptEstimate(op, mean(sigma2), C=1, sigma2final=sigma2, opt.criterion="FLCI")

## -----------------------------------------------------------------------------
ATTEffBounds(op, mean(sigma2), C=1)

## -----------------------------------------------------------------------------
op <- ATTOptPath(path=op, maxsteps=290)
ATTEffBounds(op, mean(sigma2), C=1)

## -----------------------------------------------------------------------------
ATTMatchEstimate(ATTMatchPath(y, d, D0, M=1, tol=1e-12), mean(sigma2),
C=1, sigma2final=sigma2)

## -----------------------------------------------------------------------------
mp <- ATTMatchPath(y, d, D0, M=1:10, tol=1e-12)
ATTMatchEstimate(mp, mean(sigma2), C=1, sigma2final=sigma2, opt.criterion="FLCI")
ATTMatchEstimate(mp, mean(sigma2), C=1, sigma2final=sigma2, opt.criterion="RMSE")

