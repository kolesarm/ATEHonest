context("Check solution path for ATT in NSW")

test_that("Solution path for ATT in NSW", {

## main specification (based on 1980 census regressions)
X <- as.matrix(NSW[, 2:10])
Ahalf <- diag(c(0.15, 0.6, 2.5, 2.5, 2.5, 0.5, 0.5, 0.1, 0.1))
d <- NSW$treated
y <- NSW$re78

D0 <- distMat(X, Ahalf, method="manhattan", d)

## Distance matrix for variance estimation
DMvar <- distMat(X, chol(solve(cov(X))), method="euclidean", d)

## Check matching against Tim's code
attm <- function(M, D0)
    unname(unlist(ATTMatchPath(y, d, D0, M)$ep[, c("att", "maxbias")]))

expect_equal(attm(M= 1, D0), c(1.391618501, 1.48333929))
expect_equal(attm(M= 2, D0), c(1.69911476, 1.61376836))
expect_equal(attm(M=10, D0), c(1.64632730, 1.98280592))
expect_equal(attm(M=4, DMvar), c(1.42077963, 0.96131554))
})


context("Matching estimates ")

test_that("Check we match Abadie and Imbens (2011, JBES) estimates", {

    X <- as.matrix(ATEHonest::NSWexper[, 2:10])
    d <- ATEHonest::NSWexper$treated
    y <- ATEHonest::NSWexper$re78
    D2 <- ATEHonest::distMat(X, diag(1/sqrt(diag(var(X)))),
                             method="euclidean", d)

    mp2 <- ATEHonest::ATTMatchPath(y, d, D2, M=c(1:4, 16, 64, Inf), tol=1e-12)

    expect_equal(round(mp2$ep$att[c(1, 4:7)], 4),
                 c(1.2232, 1.9946, 1.7533, 2.2049, 1.7943))

    X <- as.matrix(ATEHonest::NSW[, 2:10])
    d <- ATEHonest::NSW$treated
    y <- ATEHonest::NSW$re78
    D2 <- ATEHonest::distMat(X, diag(1/sqrt(diag(var(X)))),
                             method="euclidean", d)

    mp2 <- ATEHonest::ATTMatchPath(y, d, D2, M=c(1:4, 16, 64), tol=1e-12)

    expect_equal(round(mp2$ep$att[c(1, 4:6)], 4),
                 c(2.0735,  1.6187,  0.4692, -0.1114))

})
