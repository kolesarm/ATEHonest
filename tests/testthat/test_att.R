context("Check solution path for ATT")

test_that("Solution path for ATT in small examples", {

    distMat <- function(x0, x1) {
        Dm <- unname(as.matrix(dist(c(x0, x1))))
        Dm[(length(x0)+1):(length(x0)+length(x1)), 1:length(x0)]
    }

    x0 <- c(3, 4, 7, 7, 9)
    x1 <- c(1, 4, 5, 7)
    expect_silent(t0 <- ATTh(distMat(x0, x1), check=TRUE))
    expect_lt(nrow(t0), 10)

    x0 <- c(1, 7, 7, 8, 9)
    x1 <- c(1, 2, 4, 8)
    expect_silent(t1 <- ATTh(distMat(x0, x1), check=TRUE))
    expect_lt(nrow(t1), 7)

    x0 <- c(2, 3, 7, 7, 7, 7, 8, 10, 10, 10)
    x1 <- c(1, 1, 1, 1, 1, 2, 5, 8, 8, 9, 10)
    expect_silent(t2 <- ATTh(distMat(x0, x1), check=TRUE))
    expect_lt(nrow(t2), 9)

    x0 <- c(1,  1,  3,  3,  4,  6, 17, 17, 21, 22, 23, 23, 26, 30, 36, 37, 47,
            53, 58, 61)
    x1 <- c(9, 15, 26, 27, 28, 31, 32, 40, 51, 52, 56)
    expect_silent(t3 <- ATTh(distMat(x0, x1), check=TRUE))



    set.seed(42)
    x0 <- rnorm(100)
    x1 <- rnorm(100)
    expect_silent(ATTh(distMat(x0, x1), check=TRUE))



    ## CVXR::installed_solvers()
})
