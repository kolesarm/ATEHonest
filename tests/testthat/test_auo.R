context("Check solution path for average untreated outcome")

test_that("Solution path for AUO in small examples", {

    distMat <- function(x0, x1) {
        Dm <- unname(as.matrix(dist(c(x0, x1))))
        Dm[(length(x0)+1):(length(x0)+length(x1)), 1:length(x0)]
    }

    x0 <- c(3, 4, 7, 7, 9)
    x1 <- c(1, 4, 5, 7)
    D0 <- distMat(x0, x1)
    r <- AUOh(D0, check=TRUE)
    expect_silent(AUOh(D0, check=TRUE))

    x0 <- c(1, 7, 7, 8, 9)
    x1 <- c(1, 2, 4, 8)
    expect_silent(AUOh(distMat(x0, x1), check=TRUE))

    x0 <- c(2, 3, 7, 7, 7, 7, 8, 10, 10, 10)
    x1 <- c(1, 1, 1, 1, 1, 2, 5, 8, 8, 9, 10)
    expect_silent(AUOh(distMat(x0, x1), check=TRUE))


    x0 <- c(1,  1,  3,  3,  4,  6, 17, 17, 21, 22, 23, 23, 26, 30, 36, 37, 47,
            53, 58, 61)
    x1 <- c(9, 15, 26, 27, 28, 31, 32, 40, 51, 52, 56)
    ## r <- AUOh(distMat(x0, x1), check=TRUE)


    ## TODO: check we're not taking too many steps

    ## CVXR::installed_solvers()
})
