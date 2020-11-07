## ---- include=FALSE, cache=FALSE----------------------------------------------
library("knitr")
knitr::opts_knit$set(self.contained = FALSE)
knitr::opts_chunk$set(tidy = TRUE, collapse=TRUE, comment = "#>",
                      tidy.opts=list(blank=FALSE, width.cutoff=55))
oldoptions <- options(digits=5)

## -----------------------------------------------------------------------------
library("ATEHonest")
X <- as.matrix(NSWexper[, 2:10])
d <- NSWexper$treated
y <- NSWexper$re78

## -----------------------------------------------------------------------------
Ahalf <- diag(c(0.15, 0.6, 2.5, 2.5, 2.5, 0.5, 0.5, 0.1, 0.1))
D0 <- distMat(X, Ahalf, method="manhattan", d)

## -----------------------------------------------------------------------------
DM <- distMat(X, chol(solve(cov(X))), method="euclidean")

## -----------------------------------------------------------------------------
c1 <- ATTOptEstimate(y=y, d=d, D0=D0, C=1, DM=DM, opt.criterion="RMSE")
print(c1)

## -----------------------------------------------------------------------------
ATTOptEstimate(path=c1$path, C=1, DM=DM, opt.criterion="FLCI")

## -----------------------------------------------------------------------------
ATTEffBounds(c1$path, DM=DM, C=1)

## -----------------------------------------------------------------------------
op <- ATTOptPath(path=c1$path, maxsteps=290)
ATTEffBounds(op, DM=DM, C=1)

## -----------------------------------------------------------------------------
ATTMatchEstimate(ATTMatchPath(y, d, D0, M=1, tol=1e-12), C=1, DM=DM)

## -----------------------------------------------------------------------------
mp <- ATTMatchPath(y, d, D0, M=1:10, tol=1e-12)
ATTMatchEstimate(mp, C=1, DM=DM, opt.criterion="FLCI")
ATTMatchEstimate(mp, C=1, DM=DM, opt.criterion="RMSE")

## ----cleanup, include=FALSE---------------------------------------------------
options(oldoptions)

