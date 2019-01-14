context("Check solution path for ATT")

test_that("Solution path for ATT in small examples", {

    dd <- function(x0, x1) c(rep(FALSE, length(x0)),
                             rep(TRUE, length(x1)))

    x0 <- list(c(3, 4, 7, 7, 9),
               c(1, 7, 7, 8, 9),
               c(2, 3, 7, 7, 7, 7, 8, 10, 10, 10),
               c(1,  1,  3,  3,  4,  6, 17, 17, 21, 22, 23, 23, 26, 30, 36, 37,
                 47, 53, 58, 61),
               c(-2.19632667892222866, -2.15997868323220610,
                 -1.31097131989365789, -1.16421297762310583,
                 -0.61199606702408638, -0.37223289750184541,
                 -0.35269074554765034, -0.22176143627504877,
                 -0.00749056392427201, 0.12011828315494359, 0.18267388345099736,
                 0.18593160956753849, 0.22915318333472903, 0.29678161157067251,
                 0.42253423096032311, 0.67552340069946715, 1.14904473037269139,
                 1.68623728943529416, 1.85267571248565455, 3.17444522610536684))

    x1 <- list(c(1, 4, 5, 7),
               c(1, 2, 4, 8),
               c(1, 1, 1, 1, 1, 2, 5, 8, 8, 9, 10),
               c(9, 15, 26, 27, 28, 31, 32, 40, 51, 52, 56),
               c(-0.417160506610267, -0.400400493634091, -0.323482341050380,
                 0.191550650274630, 0.250676837337075,  0.388565760568126,
                 0.478185868404733,  0.525553657035341, 1.044897968992634,
                 1.181192580509525))

    tt <- list()
    for (j in seq_along(x0)) {
        d <- dd(x0[[j]], x1[[j]])
        D0 <- distMat(c(x0[[j]], x1[[j]]), d=d)
        if (j != 4)
            expect_silent(tt[[j]] <- ATTh(D0, check=TRUE)$res)
        ## For j==4 the solution paths on xps and travis do not match officepc
        ## (where there is no message)
        ## expect_message(tt[[j]] <- ATTh(D0, check=TRUE)$res)

        y <- d
        ## Check maximum bias matches that from LP
        res <- tt[[j]][tt[[j]][, "delta"]!=0, ]
        r0 <- ATTOptEstimate(ATTOptPath(res, y, d), sigma2=1, C=0.5)
        expect_equal(0.5*ATTbias(r0$w, D0), r0$e$maxbias)
        r1 <- ATTOptEstimate(ATTOptPath(res, y, d), sigma2=1, C=1)
        expect_equal(ATTbias(r1$w, D0), r1$e$maxbias)
        ## Check optimum matches CVX
        rmse <- function(delta, C)
            ATTOptEstimate(ATTOptPath(
                res=matrix(ATTbrute(delta2=delta^2, D0), nrow=1), y, d),
                sigma2=1, C=C)$e$rmse

        cvx0 <- unlist(optimize(function(de)
            rmse(de, 0.5), c(min(res[, "delta"]), max(res[, "delta"]))))
        cvx1 <- unlist(optimize(function(de)
            rmse(de, 1), c(min(res[, "delta"]), max(res[, "delta"]))))
        hom0 <- unlist(r0$e[, c("delta", "rmse")])
        hom1 <- unlist(r1$e[, c("delta", "rmse")])

        expect_lt(abs(cvx1[2]-hom1[2]), 1e-3)
        expect_lt(abs(cvx0[2]-hom0[2]), 1e-3)
        ## The deltas sometimes differ
        ## expect_lt(abs(cvx1[1]-hom1[1]), 1e-4)
        ## expect_lt(abs(cvx0[1]-hom0[1]), 1e-4)
    }

    expect_lt(nrow(tt[[1]]), 8)
    expect_lt(nrow(tt[[2]]), 9)
    expect_lt(nrow(tt[[3]]), 9)
    expect_lt(nrow(tt[[4]]), 25)

    ## TODO: CVX solver fails at more iterations
    set.seed(42)
    x0 <- sort(rnorm(100))
    x1 <- sort(rnorm(100))
    expect_silent(t5 <- ATTh(distMat(c(x0, x1), dd(x0, x1)), check=TRUE,
                                     maxiter=100)$res)
    ## CVXR::installed_solvers()
})
