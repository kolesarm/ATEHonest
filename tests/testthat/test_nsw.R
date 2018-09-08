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
attm <- function(M, D0)
    unname(unlist(ATTMatchPath(y, d, sigma2, D0, C=1, M)$ep[, c("att", "maxbias")]))

expect_equal(attm(M= 1, D0), c(1.391618501, 1.48333929))
expect_equal(attm(M= 2, D0), c(1.69911476, 1.61376836))
expect_equal(attm(M=10, D0), c(1.64632730, 1.98280592))
expect_equal(attm(M=4, DMvar[d==1, d==0]), c(1.42077963, 0.96131554))
})
