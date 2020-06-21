## Matching weights for untreated
ATTMatchW <- function(D0, M, tol) {
    n0 <- ncol(D0)
    w0 <- rep(0, n0)
    if (M==Inf)
        return(-(w0+1)/n0)

    for (i in seq_len(nrow(D0))) {
        ## Find NN of i, within tolerance
        idx <- D0[i, ] <= sort(D0[i, ])[M]+tol
        w0[idx] <- w0[idx] + 1/sum(idx)
    }
    -w0/nrow(D0)
}


## update estimation path with new C or new variance
UpdatePath <- function(ep, resw, Cratio, sigma2, alpha=0.05, beta=0.8) {
    if (length(sigma2)==1) sigma2 <- rep(sigma2, ncol(resw))
    ep$maxbias <- Cratio*ep$maxbias
    ep$sd <- sqrt(drop(resw^2 %*% sigma2))
    ep$hl <- cv(ep$maxbias/ep$sd, alpha) * ep$sd
    ep$lower <- ep$att - ep$maxbias - stats::qnorm(1-alpha)*ep$sd
    ep$upper <- ep$att + ep$maxbias + stats::qnorm(1-alpha)*ep$sd
    ## worst-case quantile of excess length
    ep$maxel <- 2*ep$maxbias +
        ep$sd * (stats::qnorm(1-alpha)+stats::qnorm(beta))
    ep$rmse <- sqrt(ep$sd^2+ep$maxbias^2)
    ep
}


#' Compute the matching estimator for the CATT
#'
#' Computes the matching estimating and matching weights for a range of matches
#' \code{M}. The output of this function is used as an input for
#' \code{\link{ATTMatchEstimate}} for inference on the CATT.
#' @param M a vector determining the number of matches. If \code{Inf}, then use
#'     the simple difference in means estimator.
#' @template D0
#' @template data
#' @param tol numerical tolerance for determining nearest neighbors in
#'     constructing matches
#' @examples
#' Ahalf <- diag(c(0.15, 0.6, 2.5, 2.5, 2.5, 0.5, 0.5, 0.1, 0.1))
#' D0 <- distMat(NSWexper[, 2:10], Ahalf, method="manhattan", NSWexper$treated)
#' ATTMatchPath(NSWexper$re78, NSWexper$treated, D0, M=c(1, 2), tol=1e-12)
#' @export
ATTMatchPath <- function(y, d, D0, M=1:25, tol=1e-12) {
    n1 <- nrow(D0)
    maxbias <- vector(length=length(M))
    resw <- matrix(nrow=length(M), ncol=(ncol(D0)+n1))
    resw[, d==1] <- 1/n1

    for (j in seq_along(M)) {
        resw[j, d==0] <- ATTMatchW(D0, M[j], tol)
        maxbias[j] <- ATTbias(resw[j, d==0], D0)
    }
    list(ep=data.frame(att=drop(resw %*% y), maxbias=maxbias, M=M),
         resw=resw, d=d)
}


#' Inference on the CATT using the matching estimator
#'
#' Computes matching estimator and confidence intervals (CIs) for the CATT. If
#' \code{ATTMatchPath} used a single \code{M}, the estimator and CIs are based
#' on a matching estimator with this number of matches. Otherwise, optimize the
#' number of matches from the set in \code{M} according to \code{opt.criterion}.
#' @param mp Output of \code{ATTMatchPath}
#' @param alpha Level of confidence interval, \code{1-alpha}.
#' @param beta The quantile \code{beta} of excess length for determining
#'     performance of one-sided CIs.
#' @param sigma2 Estimate of the conditional variance of the outcome, used to
#'     optimize the tuning parameter.
#' @param C Lipschitz smoothness constant
#' @param sigma2final vector of variance estimates with length{n} for
#'     determining standard error of the optimal estimators. In contrast,
#'     \code{sigma2} is used only for determining the optimal tuning parameter.
#' @param opt.criterion One of \code{"RMSE"} (root mean squared error),
#'     \code{"OCI"} (one-sided confidence intervals), \code{"FLCI"}
#'     (fixed-length two-sided confidence intervals)
#' @examples
#' Ahalf <- diag(c(0.15, 0.6, 2.5, 2.5, 2.5, 0.5, 0.5, 0.1, 0.1))
#' D0 <- distMat(NSWexper[, 2:10], Ahalf, method="manhattan", NSWexper$treated)
#' mp <- ATTMatchPath(NSWexper$re78, NSWexper$treated, D0, M=c(1, 2), tol=1e-12)
#' ## Distance matrix for variance estimation
#' DM <- distMat(NSWexper[, 2:10], Ahalf, method="manhattan")
#' sigma2 <- nnvar(DM, NSWexper$treated, NSWexper$re78, J=3)
#' ## Estimator and CI based on a single match
#' ATTMatchEstimate(mp, mean(sigma2), C=1, sigma2final=sigma2)
#' @export
ATTMatchEstimate <- function(mp, sigma2, C=1, sigma2final=sigma2, alpha=0.05,
                             beta=0.8, opt.criterion="RMSE") {
    ## Update estimate path with new value of C
    mp$ep <- UpdatePath(mp$ep, mp$resw, C, sigma2, alpha, beta)

    ## Index of criterion to optimize
    idx <- if (opt.criterion=="RMSE") {
               which.max(names(mp$ep)=="rmse")
           } else if (opt.criterion=="OCI") {
               which.max(names(mp$ep)=="maxel")
           } else if (opt.criterion=="FLCI") {
               which.max(names(mp$ep)=="hl")
           }
    i <- which.min(mp$ep[[idx]])
    if (i==nrow(mp$ep) & nrow(mp$ep)>1)
        warning("Optimum found at end of path")

    ## Robust se, C=1 to keep bias the same
    er <- UpdatePath(mp$ep[i, ], mp$resw[i, ], Cratio=1, sigma2final, alpha,
                     beta)

    structure(list(e=cbind(mp$ep[i, ],
                           data.frame(rsd=er$sd, rlower=er$lower,
                                      rupper=er$upper, rhl=er$hl, rrmse=er$rmse,
                                      rmaxel=er$maxel, C=C)),
                   w=mp$resw[i, mp$d==0]), class="ATTEstimate")
}


## Compute worst-case bias
ATTbias <- function(w, D0) {
    ## If w receives zero weight, we can drop it
    idx <- w != 0
    w0 <- w[idx]
    D0 <- D0[, idx]

    n1 <- nrow(D0)
    n0 <- ncol(D0)

    ## Use LP
    f.obj <- c(w0, rep(1, n1)/n1)
    f.rhs <- as.vector(D0)
    f.dir <- rep("<=", length(f.rhs))
    ## Constraint number, column number, value
    f.con <- cbind(rep(seq_along(f.rhs), each=2),
                   as.vector(t(expand.grid((n0+1):(n0+n1), 1:n0))),
                   rep(c(1, -1), length(f.rhs)))

    r <- lpSolve::lp("max", f.obj, , f.dir, f.rhs, dense.const=f.con) #nolint
    r$objval
}


## Weights from solution path
ATTOptW <- function(res, d) {
    n0 <- length(d)-sum(d)
    m <- res[, 2:(n0+1), drop=FALSE]
    resw <- matrix(1/sum(d), nrow=nrow(res), ncol=length(d))
    resw[, d==0] <- -m/rowSums(m)
    resw
}


#' Compute the class of optimal linear estimators for the CATT
#'
#' Computes the class of optimal linear estimators that minimize variance
#' subject to a bound on bias, and the optimal weights. The output of this
#' function is used by \code{\link{ATTOptEstimate}} for optimal estimation and
#' inference on the CATT.
#' @param res The \code{res} element of the output of \code{ATTh} on which to
#'     base the estimate
#' @template data
#' @export
ATTOptPath <- function(res, y, d) {
    n <- length(d)
    n0 <- n-sum(d)
    ## Vector or matrix?
    if (!is.matrix(res))
        stop("res needs to be a matrix")
    rmean <- rowMeans(res[, (n0+2):(n+1), drop=FALSE])
    m <- res[, 2:(n0+1), drop=FALSE]
    maxbias <- rmean - apply(m, 1, function(x) sum(x^2)) / rowSums(m)
    resw <- ATTOptW(res, d)

    list(ep=data.frame(att=drop(resw %*% y), maxbias=maxbias,
                       delta=unname(res[, 1]),
                       omega=unname(2*(res[, n+2]+rmean))),
         resw=resw, y=y, d=d, res=res)
}

#' Optimal estimation and inference for the CATT
#'
#' Computes the estimator and confidence intervals (CIs) for the CATT. The
#' tuning parameter is chosen to optimize \code{opt.criterion} criterion.
#' @param op Output of \code{ATTOptPath}.
#' @inheritParams ATTMatchEstimate
#' @examples
#' Ahalf <- diag(c(0.15, 0.6, 2.5, 2.5, 2.5, 0.5, 0.5, 0.1, 0.1))
#' D0 <- distMat(NSWexper[, 2:10], Ahalf, method="manhattan", NSWexper$treated)
#' ## Distance matrix for variance estimation
#' DM <- distMat(NSWexper[, 2:10], Ahalf, method="manhattan")
#' sigma2 <- nnvar(DM, NSWexper$treated, NSWexper$re78, J=3)
#' ## Compute homotopy/solution path, and the class of optimal estimators based
#' ## on the solution path
#' attpath <- ATTh(D0, maxiter=200)$res
#' op <- ATTOptPath(attpath, NSWexper$re78, NSWexper$treated)
#' ATTOptEstimate(op, mean(sigma2), C=1, sigma2final=sigma2,
#'                opt.criterion="RMSE")
#' ATTOptEstimate(op, mean(sigma2), C=1, sigma2final=sigma2,
#'                opt.criterion="FLCI")
#' @export
ATTOptEstimate <- function(op, sigma2, C=1, sigma2final=sigma2, alpha=0.05,
                           beta=0.8, opt.criterion="RMSE") {
    ## Drop delta=0, back out weights
    keep <- op$ep$delta > 0
    res <- op$res[keep, , drop=FALSE] # nolint
    ep <- UpdatePath(op$ep[keep, ], op$resw[keep, , drop=FALSE], # nolint
                     C, sigma2, alpha, beta)
    ## update delta
    if (length(sigma2)>1)
        warning(paste("The solution path is assuming homoskedastic variance,",
                      "but supplied sigma2 implies it's not homoskedastic"))
    ep$delta <- res[, 1] <- ep$delta/sqrt(mean(sigma2))

    up <- function(res) {
        r <- ATTOptPath(res, op$y, op$d)
        UpdatePath(r$ep, r$resw, C, sigma2, alpha, beta)
    }

    if (nrow(res)==1) {
        resopt <- res
    } else {

    ## Index of criterion to optimize
    idx <- if (opt.criterion=="RMSE") {
               which.max(names(ep)=="rmse")
           } else if (opt.criterion=="OCI") {
               which.max(names(ep)=="maxel")
           } else if (opt.criterion=="FLCI") {
               which.max(names(ep)=="hl")
           }
    i <- which.min(ep[[idx]])
    ## Assume minimum is either in [i, i+1] or [i-1, i], which is true if
    ## criterion is convex. Sometimes solution can be stuck, so not enough to
    ## take i+1 (but since which.min takes first minimum, we can always take
    ## i-1)
    ip <- i+1
    while(ip<=nrow(res) & ep$delta[i]==ep$delta[ip])
        ip <- ip+1

    if (ip<=nrow(res)) {
        f1 <- function(w)
            up((1-w)*res[i, , drop=FALSE]+w*res[i+1, , drop=FALSE])[[idx]] # nolint
        opt1 <- stats::optimize(f1, interval=c(0, 1))
    } else {
        opt1 <- list(minimum=0, objective=min(ep[[idx]]))
    }

    if (i>1) {
        f0 <- function(w)
            up((1-w)*res[i-1, , drop=FALSE]+w*res[i, , drop=FALSE])[[idx]] # nolint
        opt0 <- stats::optimize(f0, interval=c(0, 1))
    } else {
        opt0 <- list(minimum=1, objective=min(ep[[idx]]))
    }

    if (opt1$objective < opt0$objective) {
        resopt <- (1-opt1$minimum)*res[i, , drop=FALSE] +  # nolint
            opt1$minimum*res[min(i+1, nrow(res)), , drop=FALSE]  # nolint
    } else {
        resopt <- (1-opt0$minimum)*res[max(i-1, 1), , drop=FALSE] +  # nolint
            opt0$minimum*res[i, , drop=FALSE]  # nolint
    }
    }
    r1 <- up(resopt)
    resw1 <- ATTOptW(resopt, op$d)
    ## C=1 to keep bias the same
    r2 <- UpdatePath(r1, resw1, Cratio=1, sigma2final, alpha, beta)

    if (r1$delta==max(ep$delta) & nrow(res)>1)
        warning("Optimum found at end of path")

    structure(list(e=cbind(r1, data.frame(rsd=r2$sd, rlower=r2$lower,
                                          rupper= r2$ upper, rhl=r2$hl,
                                          rrmse=r2$rmse, rmaxel=r2$maxel, C=C)),
                   res=resopt, w=resw1[, op$d==0]),
              class="ATTEstimate")
}

#' @export
print.ATTEstimate <- function(x, digits = getOption("digits"), ...) {

    fmt <- function(x) format(x, digits=digits, width=digits+1)


    r <- x$e[, c("att", "maxbias", "rsd", "rhl")]
    r <- cbind(r, l=r$att-r$rhl, u=r$att+r$rhl)
    r <- fmt(r)
    r$ci <- paste0("(", r$l, ", ",  r$u, ")")
    r$rhl <- r$l <- r$u <- NULL

    names(r) <- c("Estimate", "Max. bias", "SE", "CI")

    if ("M" %in% names(x$e))
        r$M <- c("M"=x$e$M)
    else
        r$delta <- c("delta"=x$e$delta)

    print(knitr::kable(r))
    cat("\n")
    invisible(x)
}


#' Efficiency bounds for confidence intervals
#'
#' Computes the asymptotic efficiency of two-sided fixed-length confidence
#' intervals at smooth functions, as well as the efficiency of one-sided
#' confidence intervals that optimize a given \code{beta} quantile of excess
#' length, using the formula described in Appendix A of Armstrong and Kolesár
#' (2018)
#' @inheritParams ATTMatchEstimate
#' @param res The \code{res} element of the output of \code{ATTh}.
#' @param d Vector of treatment indicators with length \code{n}
#' @references \cite{Armstrong, T. B., and M. Kolesár (2018): Finite-Sample
#'     Optimal Estimation and Inference on Average Treatment Effects Under
#'     Unconfoundedness, Unpublished manuscript}
#' @export
ATTEffBounds <- function(res, d, sigma2, C=1, beta=0.8, alpha=0.05) {
    if (length(sigma2)>1)
        warning(paste("The solution path is assuming homoskedastic variance,",
                      "but supplied sigma2 implies it's not homoskedastic"))

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
    d0 <- (zal + stats::qnorm(beta)) * sig
    if (2*d0>max(del0)) {
        warning("Path too short to compute one-sided efficiency")
        eff1 <- NaN
    } else
        eff1 <- mod11(2*d0)$omega/sum(unlist(mod11(d0))*c(1, d0))


    integrand <- function(z)
        vapply(z, function(z)
            mod11(2*(zal-z) * sig)$omega * stats::dnorm(z), numeric(1))
    ## Maximum integrable point
    lbar <- zal-max(del0)/(2*sig)

    if (integrand(lbar)>1e-8) {
        warning("Path too short to compute two-sided efficiency")
        num <- den <- NaN
    } else {
        lo <- -zal                          # lower endpoint
        while(integrand(lo)>1e-8) lo <- max(lo-2, lbar)
        num <- stats::integrate(integrand, lo, zal, abs.tol=1e-6)$value



        len <- function(del) {
            if (del==0)
                return(Inf)
            m <- mod11(del*sig)
            2*sig*cv((m$omega/(m$domega*sig)-del)/2, alpha)*m$domega
        }
        den <- stats::optimize(len, interval=c(1e-6, max(del0/sig)))$objective
    }
    list(onesided=eff1, twosided=num/den)
}
