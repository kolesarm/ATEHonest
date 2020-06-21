context("Check computation of distance matrix")

test_that("Various ways of computing the distance matrix", {
    Ahalf <- diag(c(0.15, 0.6, 2.5, 2.5, 2.5, 0.5, 0.5, 0.1, 0.1))
    d <- NSW$treated
    X <- as.matrix(NSW[, 2:10])
    expect_equal(sum(distMat(X, chol(solve(cov(X))),
                             method="euclidean"))-28130262.545790588, 0)
    expect_equal(sum(distMat(X, Ahalf, method="manhattan", d))-
                 11844361.9818461447, 0)

    d <- NSWexper$treated
    ## Add a few more lines to X, so it's above 500 cutoff
    X <- as.matrix(NSWexper[, 2:10])
    X <- rbind(X, X[1:50, ])
    d <- c(d, d[1:50])
    expect_equal(distMat(X, var(X), method="maximum"),
                 unname(as.matrix(stats::dist(X %*% var(X),
                                              method="maximum"))))
    expect_equal(distMat(X, Ahalf, method="euclidean"),
                 unname(as.matrix(stats::dist(X %*% Ahalf,
                                              method="euclidean"))))
    expect_equal(distMat(X, Ahalf, method="manhattan"),
                 unname(as.matrix(stats::dist(X %*% Ahalf,
                                              method="manhattan"))))
})
