context("Check solution path for ATT in NSW")

test_that("Solution path for ATT in NSW", {

## main specification (based on 1980 census regressions)
X <- as.matrix(NSW[, 2:10])
Ahalf <- diag(c(0.15, 0.6, 2.5, 2.5, 2.5, 0.5, 0.5, 0.1, 0.1))
d <- NSW$treated
y <- NSW$re78

D0 <- distMat(X, d, Ahalf, method="manhattan", FullMatrix=FALSE)

## Distance matrix for variance estimation
Avar <- chol(solve(cov(X)))
DMvar <- distMat(X, d, Avar, method="euclidean",
                 FullMatrix=TRUE)

sigma2 <- nnvar(DMvar, d, y, J=30)


## Check matching against Tim's code
expect_equal(ATTmatch(M=2, y, d, D0, sigma2)$att, 1.69911476)
expect_equal(ATTmatch(M=10, y, d, D0, sigma2)$att, 1.64632730)
expect_equal(ATTmatch(M=1, y, d, D0, sigma2)$att, 1.391618501)
expect_equal(ATTmatch(M=4, y, d, DMvar[d==1, d==0], sigma2)$att, 1.42077963)



## r <- ATTh(D0, maxiter=500)
## res <- r$res[r$res[, "delta"]>0, ]
## rr <- as.data.frame(ATTpath(res, d, y, C=1, mean(sigma2)))
## rr$delta <- rr$delta/sqrt(mean(sigma2))
## ## Plot the results
## d1 <- reshape2::melt(rr[, c(1, 2, 3, 4, 7, 8)],
##                           id.vars="delta")
## qplot(x=delta, y=value, data=d1, geom="line")+facet_grid( variable ~ ., scales="free_y")

## ## Optimal results
## rr[c(which.min(rr$rmse), which.min(rr$hl)), ]




})
