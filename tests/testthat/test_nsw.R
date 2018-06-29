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
})
