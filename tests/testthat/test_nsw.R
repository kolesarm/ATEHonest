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

context("Efficiency calculations")

test_that("Alternative way of computing modulus efficiency", {

    X <- as.matrix(ATEHonest::NSWexper[, 2:10])
    d <- ATEHonest::NSWexper$treated
    D0 <- distMat(X, diag(c(0.15, 0.6, 2.5, 2.5, 2.5, 0.5, 0.5, 0.1, 0.1)),
                  method="manhattan", d)
    sigma2 <- 40

    res <- ATEHonest::ATTh(D0, maxiter=300)$res

    eb <- ATTEffBounds(res, d, sigma2, C=1)
    expect_equal(eb$onesided, 0.991823078)
    ## Alternative modulus calculation
    expect_warning(ATTEffBounds(res[1:100, ], d, sigma2, C=1))

    ATTEffBounds2 <- function(res, d, sigma2, C=1, beta=0.8, alpha=0.05) {
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
                 domega=0.5*del/(n1 * ((1-w)*mu0[idx-1]+w*mu0[idx]) ))
        }

        ## One-sided
        zal <- stats::qnorm(1-alpha)
        ## Rescaling to modulus(C, sigma)
        sig <- sqrt(mean(sigma2)) / C

        integrand <- function(z)
            sapply(z, function(z)
                mod11(2*(zal-z) * sig)$omega * stats::dnorm(z))
        ## Maximum integrable point
        lbar <- zal-max(del0)/(2*sig)

        lo <- -zal                          # lower endpoint
        while(integrand(lo)>1e-8) lo <- max(lo-2, lbar)
        num <- stats::integrate(integrand, lo, zal, abs.tol=1e-6)$value
        den <- 2*ATTOptEstimate(ATTOptPath(res, d, d), mean(sigma2), C=C,
                                sigma2final=mean(sigma2), alpha,
                                opt.criterion="FLCI")$e$hl
        C*num/den
    }
    expect_equal(eb$twosided, ATTEffBounds2(res, d, sigma2, C=1))
    expect_equal(ATTEffBounds(res, d, sigma2, C=4)$twosided,
                 ATTEffBounds2(res, d, sigma2, C=4))
})
