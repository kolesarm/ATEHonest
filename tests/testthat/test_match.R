context("Check matching solution path in NSW")

test_that("Matching path for ATT in NSW", {

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


context("Matching estimates")

test_that("Check we match Abadie and Imbens (2011, JBES) estimates", {
    X <- as.matrix(NSWexper[, 2:10])
    d <- NSWexper$treated
    y <- NSWexper$re78
    D2 <- distMat(X, diag(1/sqrt(diag(var(X)))),
                  method="euclidean", d)

    mp2 <- ATTMatchPath(y, d, D2, M=c(1:4, 16, 64, Inf), tol=1e-12)

    expect_equal(round(mp2$ep$att[c(1, 4:7)], 4),
                 c(1.2232, 1.9946, 1.7533, 2.2049, 1.7943))

    X <- as.matrix(NSW[, 2:10])
    d <- NSW$treated
    y <- NSW$re78
    D2 <- distMat(X, diag(1/sqrt(diag(var(X)))), method="euclidean", d)

    mp2 <- ATTMatchPath(y, d, D2, M=c(1:4, 16, 64), tol=1e-12)

    expect_equal(round(mp2$ep$att[c(1, 4:6)], 4),
                 c(2.0735,  1.6187,  0.4692, -0.1114))

})

context("Matching estimator for ATT in NSWexper")
test_that("Standard errors for matching estimator are correct", {
    Ahalf <- diag(c(0.15, 0.6, 2.5, 2.5, 2.5, 0.5, 0.5, 0.1, 0.1))
    D0 <- distMat(NSWexper[, 2:10], Ahalf, method="manhattan", NSWexper$treated)
    mp <- ATTMatchPath(NSWexper$re78, NSWexper$treated, D0, M=2, tol=1e-12)
    DM <- distMat(NSWexper[, 2:10], Ahalf, method="manhattan")
    r <- ATTMatchEstimate(mp, C=1, DM=DM)
    expect_equal(unname(r$e[c("sd", "hl")]), c(0.692625208, 2.437109605))
    o1 <- utils::capture.output(print(r, digits=6))
    e1 <- c("", "",
            "|     | Estimate| Max. bias|       SE|CI                   |  M|",
            "|:----|--------:|---------:|--------:|:--------------------|--:|",
            "|CATT |  1.92327|   1.29784| 0.709765|(-0.542033, 4.38857) |  2|",
            "|PATT |  1.92327|          | 0.738429|(-0.821869, 4.66841) |   |")
    expect_equal(o1, e1)
})

context("Variance for the ATE")
test_that("Direct calculation of unconditional variance", {

    nnUncondVar <- function(D0, M, tol, d, y, sigma2) {
        if (length(sigma2)==1)
            sigma2 <- rep(sigma2, length(d))
        Lam <- D0*0
        Y0 <- vector(length=nrow(D0))
        for (i in seq_len(nrow(D0))) {
            ## Find NN of i, within tolerance
            idx <- D0[i, ] <= sort(D0[i, ])[M]+tol
            Y0[i] <- mean(y[!d][idx])
            Lam[i, idx] <- 1/sum(idx)^2
        }
        tau <- mean(y[d]-Y0)
        v2 <- tcrossprod(Lam %*% diag(sigma2[!d]), Lam)

        ## Unconditional and Conditional
        c((sum(v2-diag(diag(v2)))+sum((y[d]-Y0-tau)^2))/nrow(D0)^2,
        (sum(v2)+sum(sigma2[d]))/nrow(D0)^2)
    }

    x0 <- c(-2.19632667892222866, -2.15997868323220610,
                 -1.31097131989365789, -1.16421297762310583,
                 -0.61199606702408638, -0.37223289750184541,
                 -0.35269074554765034, -0.22176143627504877,
                 -0.00749056392427201, 0.12011828315494359, 0.18267388345099736,
                 0.18593160956753849, 0.22915318333472903, 0.29678161157067251,
                 0.42253423096032311, 0.67552340069946715, 1.14904473037269139,
                 1.68623728943529416, 1.85267571248565455, 3.17444522610536684)
    x1 <- c(-0.417160506610267, -0.400400493634091, -0.323482341050380,
                 0.191550650274630, 0.250676837337075,  0.388565760568126,
                 0.478185868404733,  0.525553657035341, 1.044897968992634,
                 1.181192580509525)
    d <- c(rep(FALSE, length(x0)), rep(TRUE, length(x1)))
    D0 <- distMat(c(x0, x1), d=d)
    sigma2 <- 40
    mp1 <- ATTMatchPath(d, d, D0, M=1, tol=1e-12)
    r1 <- ATTMatchEstimate(mp1, sigma2=sigma2, C=1)
    expect_equal(nnUncondVar(D0, M=1, tol=1e-12, d, y=d, sigma2),
                 unname(r1$e[c("usd", "rsd")])^2)

    mp3 <- ATTMatchPath(d, d, D0, M=3, tol=1e-12)
    r3 <- ATTMatchEstimate(mp1, C=1, sigma2=4)
    expect_equal(nnUncondVar(D0, M=1, tol=1e-12, d, y=d, 4),
                 unname(r3$e[c("usd", "rsd")])^2)


})
