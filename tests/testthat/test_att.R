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
        expect_silent(tt[[j]] <- ATTOptPath(d, d, D0, check=TRUE))

        ## Check maximum bias matches that from LP
        r0 <- ATTOptEstimate(op=tt[[j]], sigma2=1, C=0.5)
        expect_equal(0.5*ATTbias(r0$k[d==0], D0), unname(r0$e["maxbias"]))
        r1 <- ATTOptEstimate(op=tt[[j]], sigma2=1, C=1)
        expect_equal(1*ATTbias(r1$k[d==0], D0), unname(r1$e["maxbias"]))
        ## Check optimum matches CVX
        rmse <- function(delta, C, mvar) {
            op <- tt[[j]]
            op$res <- matrix(ATTbrute(delta2=delta^2, D0), nrow=1)
            ## Use previous mvar, since CVX path is not exact, determining
            ## binding constraints is tricky.
            unlist(ATTOptEstimate(op=ATTOptPath(path=op, check=FALSE),
                                  sigma2=1, C=C, mvar=mvar)$e)
        }
        expect_lt(max(abs(rmse(r0$res[1], C=0.5,
                               r0$e["usd"]^2-r0$e["rsd"]^2)-r0$e)), 1e-4)
        expect_lt(max(abs(rmse(r1$res[1], C=1.0,
                               r1$e["usd"]^2-r1$e["rsd"]^2)/r1$e-1)), 1e-3)

        ## Check delta and omega
        check_equiv <- function(res) {
            ## mu0*n1=sum(m0)
            n <- nrow(D0)+ncol(D0)
            expect_equal(rowSums(res[, 2:(ncol(D0)+1), drop=FALSE]),
                         res[, n+2]*nrow(D0))
            ## delta
            expect_equal(2*sqrt(rowSums(res[, 2:(ncol(D0)+1), drop=FALSE]^2)+
                                res[, n+2]^2*nrow(D0)),
                         res[, 1])
        }
        check_equiv(tt[[j]]$res)
        check_equiv(matrix(r0$res, nrow=1))
        check_equiv(matrix(r1$res, nrow=1))
        ## delta and omega in the $e
        check_d0 <- function(r0, C) {
            expect_equal(unname(r0$res[1])*C/1, unname(r0$e["delta"]))
            expect_equal(C*unname(2*sum(r0$res[2:(length(d)+1)])/sum(d)),
                         unname(r0$e["omega"]))
        }
        check_d0(r0, C=0.5)
        check_d0(r1, C=1)
    }
    expect_lt(nrow(tt[[1]]$res), 8)
    expect_lt(nrow(tt[[2]]$res), 9)
    expect_lt(nrow(tt[[3]]$res), 9)
    expect_lt(nrow(tt[[4]]$res), 25)

    ## Make path longer if necessary
    d <- dd(x0[[5]], x1[[5]])
    D0 <- distMat(c(x0[[5]], x1[[5]]), d=d)
    op <- ATTOptPath(d, d, D0, check=FALSE, maxsteps=10)
    expect_message(r0 <- ATTOptEstimate(op=op, sigma2=1, C=0.5))
    ## Solution at end of path
    expect_message(r0 <- ATTOptEstimate(d=d, y=d, D0=D0, sigma2=1, C=0.1))
    expect_equal(nrow(r0$op$ep), 45)

})

context("Efficiency calculations")

test_that("Alternative way of computing modulus efficiency", {
    dt <- NSWexper[c(1:20, 431:445), ]
    X <- as.matrix(dt[, 2:10])
    d <- dt$treated
    D0 <- distMat(X, diag(c(0.15, 0.6, 2.5, 2.5, 2.5, 0.5, 0.5, 0.1, 0.1)),
                  method="manhattan", d)
    sigma2 <- 40

    op <- ATTOptPath(y=d, d=d, D0=D0, maxsteps=10)
    expect_warning(ATTEffBounds(op, sigma2=sigma2, C=1))
    op <- ATTOptPath(path=op, maxsteps=50)
    eb <- ATTEffBounds(op, sigma2=sigma2, C=1)
    expect_equal(eb$onesided, 0.993176377)

    ATTEffBounds2 <- function(op, sigma2, C=1, beta=0.8, alpha=0.05) {
        res <- op$res
        d <- op$d
        n <- ncol(res)-3
        n1 <- sum(d)
        n0 <- n-n1
        ## normalize delta by standard deviation:
        del0 <- res[, "delta"]
        mu0 <- res[, "mu"]
        m0 <- res[, 2:(n0+1), drop=FALSE]
        r0 <- res[, (n0+2):(n+1), drop=FALSE]

        ## Modulus when sigma2=1 and C=1
        mod11 <- function(del) {
            idx <- which.max(del0>=del)
            if (idx==1 || del==del0[idx])
                return(list(omega=2*(mu0[idx]+mean(r0[idx, ])),
                            domega=0.5*del/(n1 * mu0[idx])))

            fn <- function(w)
                2*sqrt((n1*((1-w)*mu0[idx-1]+w*mu0[idx])^2 +
                        sum(((1-w)*m0[idx-1, ]+w*m0[idx, ])^2)))-del

            w <- stats::uniroot(fn, interval=c(0, 1))$root
            list(omega=2*((1-w)*mu0[idx-1]+w*mu0[idx] +
                          mean((1-w)*r0[idx-1, ]+w*r0[idx, ])),
                 domega=0.5*del/(n1 * ((1-w)*mu0[idx-1]+w*mu0[idx])))
        }

        ## One-sided
        zal <- stats::qnorm(1-alpha)
        ## Rescaling to modulus(C, sigma)
        sig <- sqrt(mean(sigma2)) / C

        integrand <- function(z)
            vapply(z, function(z)
                mod11(2*(zal-z) * sig)$omega * stats::dnorm(z), numeric(1))
        ## Maximum integrable point
        lbar <- zal-max(del0)/(2*sig)

        lo <- -zal                          # lower endpoint
        while(integrand(lo)>1e-8) lo <- max(lo-2, lbar)
        num <- stats::integrate(integrand, lo, zal, abs.tol=1e-6)$value
        den <- 2*ATTOptEstimate(op=op,
                                sigma2init=mean(sigma2), C=C,
                                sigma2=mean(sigma2), alpha,
                                opt.criterion="FLCI")$e["hl"]
        C*num/den
    }
    expect_lt(abs(eb$twosided-ATTEffBounds2(op, sigma2, C=1)), 1e-5)
    expect_lt(abs(ATTEffBounds(op, sigma2=sigma2, C=4)$twosided-
                 ATTEffBounds2(op, sigma2, C=4)), 1e-5)
})
