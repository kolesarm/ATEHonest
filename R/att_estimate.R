## update estimation path with new C or new variance
UpdatePath <- function(ep, K, Cratio, sigma2, alpha=0.05, beta=0.8, ucse=NULL) {
    if (length(sigma2)==1)
        sigma2 <- rep(sigma2, ncol(K))

    ep$maxbias <- Cratio*ep$maxbias
    if (is.null(ucse)) {
        ep$sd <- sqrt(drop(K^2 %*% sigma2))
        ep$hl <- cv(ep$maxbias/ep$sd, alpha) * ep$sd
    } else {
        ## PATT half-length uses usual CVs
        ep$sd <- ucse
        ep$hl <- ep$maxbias + stats::qnorm(1-alpha/2) * ep$sd
    }

    ep$lower <- ep$att - ep$maxbias - stats::qnorm(1-alpha)*ep$sd
    ep$upper <- ep$att + ep$maxbias + stats::qnorm(1-alpha)*ep$sd
    ## worst-case quantile of excess length
    ep$maxel <- 2*ep$maxbias +
        ep$sd * (stats::qnorm(1-alpha)+stats::qnorm(beta))
    ep$rmse <- sqrt(ep$sd^2+ep$maxbias^2)
    if(!is.null(ep$omega)) {
        ep$omega <- Cratio*ep$omega
        ep$delta <- Cratio*ep$delta/sqrt(mean(sigma2))
    }
    ep
}


## Weights from homotopy path
ATTOptK <- function(res, d) {
    n0 <- length(d)-sum(d)
    m <- res[, 2:(n0+1), drop=FALSE]
    K <- matrix(1/sum(d), nrow=nrow(res), ncol=length(d))
    K[, d==0] <- -m/rowSums(m)
    K
}


#' Class of optimal linear estimators for the ATT
#'
#' Use a LASSO-like algorithm to compute the solution path
#' \eqn{\{\hat{L}_{\delta}:\delta>0\}}{{hatL_delta}_{delta>0}} tracing out the
#' class of optimal linear estimators that minimize variance subject to a bound
#' on bias. The output of this function is used by \code{\link{ATTOptEstimate}}
#' for optimal estimation and inference on the CATT and PATT
#' @param maxsteps maximum number of steps in the solution path. If the full
#'     solution path is shorter than \code{maxsteps}, compute the whole path.
#' @param tol numerical tolerance for rounding error when finding the nearest
#'     neighbors. All observations with effective distance within \code{tol} of
#'     the closest are considered to be active.
#' @param check check at each step that the solution matches that obtained by
#'     direct optimization using \code{\link[CVXR]{CVXR-package}} (generic
#'     convex optimizer package).
#' @param path Optionally, supply previous output of \code{ATTOptPath}. If not
#'     provided, the path is started at the beginning (at \eqn{\delta=0}). If
#'     provided, it starts at the step where the previous call to
#'     \code{ATTOptPath} ended.
#' @template D0
#' @template data
#' @return A list with the following elements: \describe{
#'
#' \item{y}{Output vector, as supplied by \code{y}}
#'
#' \item{d}{Vector of treatment indicators, as supplied by \code{d}}
#'
#' \item{D0}{Matrix of distances, as supplied by \code{D0}}
#'
#' \item{res}{A matrix with rows corresponding to steps in the solution path, so
#'   that the maximum number of rows is \code{maxsteps}, and columns
#'   corresponding to the state variables \eqn{\delta}, \eqn{m}, \eqn{r},
#'   \eqn{\mu}, and \code{drop}.}
#'
#' \item{K}{Matrix of weights \eqn{k} associated with the optimal estimator at
#' each step}
#'
#' \item{ep}{A data frame with columns \code{delta}, \code{omega},
#' \code{maxbias}, and \code{att}, corresponding to \eqn{\delta},
#' \eqn{\omega(\delta)}, the scaled worst-case bias, and the ATT estimate.}
#'
#' }
#'
#' The remaining elements are state variables at the last step of the solution
#' path (see Appendix A in Armstrong and Kolesár (2018) for details and
#' notation): \describe{
#'
#' \item{m0}{A vector of length \code{n0} of corresponding to \eqn{m}.}
#'
#' \item{r0}{A vector of length \code{n1} of corresponding to \eqn{r}.}
#'
#' \item{mu}{A scalar corresponding to the Lagrange multiplier \eqn{\mu}.}
#'
#' \item{D}{A matrix of effective distances with dimension \code{n1} by
#' \code{n0}.}
#'
#' \item{Lam}{A sparse matrix of Lagrange multipliers \eqn{\Lambda} with
#'             dimension \code{n1} by \code{n0}.}
#'
#' \item{N0}{A sparse matrix of effective nearest neighbors with dimension
#'             \code{n1} by \code{n0}.}
#'
#' \item{drop}{An indicator if an observation has been dropped from an active
#'   set, or added.}
#'
#' }
#' @references \cite{Armstrong, T. B., and M. Kolesár (2018): Finite-Sample
#'     Optimal Estimation and Inference on Average Treatment Effects Under
#'     Unconfoundedness, \url{https://arxiv.org/abs/1712.04594}}
#' @examples
#' x0 <- c(0, 1, 2, 3)
#' x1 <- c(1, 4, 5)
#' d <- c(rep(FALSE, length(x0)), rep(TRUE, length(x1)))
#' D0 <- distMat(c(x0, x1), d=d)
#' ## Compute first three steps
#' p1 <- ATTOptPath(d, d, D0, maxsteps=3)
#' ## Compute the remaining steps, checking them against CVX solution
#' ATTOptPath(path=p1, maxsteps=50, check=TRUE)
#' @export
ATTOptPath <- function(y, d, D0, maxsteps=50, tol, path=NULL, check=FALSE) {
    if (is.null(path))
        path <- list(y=y, d=d, D0=D0)
    if (missing(tol))
        tol <- .Machine$double.eps*ncol(path$D0)*nrow(path$D0)
    ## Construct or update homotopy
    path <- ATTh(path, maxsteps, check,  tol)

    n <- length(path$d)
    n0 <- n-sum(path$d)
    ## Recompute $ep
    rmean <- rowMeans(path$res[, (n0+2):(n+1), drop=FALSE])
    m <- path$res[, 2:(n0+1), drop=FALSE]
    maxbias <- rmean - apply(m, 1, function(x) sum(x^2)) / rowSums(m)
    path$K <- ATTOptK(path$res, path$d)

    path$ep <- data.frame(att=drop(path$K %*% path$y), maxbias=maxbias,
                     delta=unname(path$res[, 1]),
                     omega=2*unname((path$res[, n+2]+rmean)),
                     lindw=apply(path$K^2, 1, max)/rowSums(path$K^2))
    path
}

#' Optimal estimation and inference for the CATT and the PATT
#'
#' Computes the estimator and confidence intervals (CIs) for the CATT and the
#' PATT. The tuning parameter \eqn{\delta}{delta} is chosen to optimize
#' \code{opt.criterion} criterion.
#' @param op Output of \code{ATTOptPath}.
#' @param sigma2init estimate of the conditional variance of the outcome, used
#'     to choose the optimal smoothing parameter \eqn{\delta}{delta}. If not
#'     supplied, use homoskedastic variance estimate based on a nearest neighbor
#'     variance estimator.
#' @param extrasteps If the optimal smoothing parameter \eqn{\delta}{delta} is
#'     attained at the end of the solution path, compute additional
#'     \code{extrasteps} steps in the solution path and estimate it again. If
#'     \code{extrasteps==0}, then do not compute extra steps.
#' @inheritParams ATTMatchEstimate
#' @return Returns an object of class \code{"ATTEstimate"}. An object of class
#'     \code{"ATTEstimate"} is a list containing the following components:
#'     \describe{
#'
#'     \item{e}{Data frame with columns TODO}
#'
#'     \item{k}{weights TODO}
#'
#'     \item{res}{TODO}
#'
#'     \item{op}{TODO}
#'
#'     }
#' @references \cite{Armstrong, T. B., and M. Kolesár (2018): Finite-Sample
#'     Optimal Estimation and Inference on Average Treatment Effects Under
#'     Unconfoundedness, \url{https://arxiv.org/abs/1712.04594}}
#' @examples
#' ## Use NSW experimental subsample with 25 treated and untreated units
#' dt <- NSWexper[c(1:25, 421:445), ]
#' Ahalf <- diag(c(0.15, 0.6, 2.5, 2.5, 2.5, 0.5, 0.5, 0.1, 0.1))
#' D0 <- distMat(dt[, 2:10], Ahalf, method="manhattan", dt$treated)
#' ## Distance matrix for variance estimation
#' DM <- distMat(dt[, 2:10], Ahalf, method="manhattan")
#' ## Compute the solution path, first 50 steps will be sufficient
#' op <- ATTOptPath(dt$re78, dt$treated, D0, maxsteps=50)
#' ATTOptEstimate(op, C=1, DM=DM, opt.criterion="RMSE")
#' ATTOptEstimate(op, C=1, DM=DM, opt.criterion="FLCI")
#' @export
ATTOptEstimate <- function(op, C=1, opt.criterion="RMSE", sigma2init, sigma2,
                           mvar, DM, alpha=0.05, beta=0.8, J=3, extrasteps=50) {
    if (missing(sigma2))
        sigma2 <- nnvar(DM, op$d, op$y, J=J)
    if (missing(sigma2init))
        sigma2init <- mean(sigma2)

    ## Drop delta=0
    keep <- op$ep$delta > 0
    res <- op$res[keep, , drop=FALSE]
    ep <- UpdatePath(op$ep[keep, ], op$K[keep, , drop=FALSE],
                     C, sigma2init, alpha, beta)
    ## update delta
    if (length(sigma2init)>1)
        warning(paste("The solution path is assuming homoskedastic variance,",
                      "but supplied sigma2init implies it's not homoskedastic"))
    up <- function(res) {
        op$res <- res
        r <- ATTOptPath(path=op, maxsteps=0, check=FALSE)
        UpdatePath(r$ep, r$K, C, sigma2init, alpha, beta)
    }

    if (nrow(res)==1) {
        resopt <- res
    } else {
        ## Index of criterion to optimize
        idx <- which.max(names(ep)==switch(opt.criterion, RMSE = "rmse",
                                           OCI = "maxel", FLCI = "hl"))
        i <- which.min(ep[[idx]])
        ## Assume minimum is either in [i, i+1] or [i-1, i], which is true if
        ## criterion is convex. Sometimes solution can be stuck, so not enough
        ## to take i+1 (but since which.min takes first minimum, we can always
        ## take i-1)
        ip <- i+1
        while(ip<=nrow(res) & ep$delta[i]==ep$delta[ip])
            ip <- ip+1

        if (ip<=nrow(res)) {
            f1 <- function(w)
                up((1-w)*res[i, , drop=FALSE]+w*res[i+1, , drop=FALSE])[[idx]]
            opt1 <- stats::optimize(f1, interval=c(0, 1))
        } else {
            opt1 <- list(minimum=0, objective=min(ep[[idx]]))
        }

        if (i>1) {
            f0 <- function(w)
                up((1-w)*res[i-1, , drop=FALSE]+w*res[i, , drop=FALSE])[[idx]]
            opt0 <- stats::optimize(f0, interval=c(0, 1))
        } else {
            opt0 <- list(minimum=1, objective=min(ep[[idx]]))
        }
        if (opt1$objective < opt0$objective) {
            resopt <- (1-opt1$minimum)*res[i, , drop=FALSE] +
                opt1$minimum*res[min(i+1, nrow(res)), , drop=FALSE]
        } else {
            resopt <- (1-opt0$minimum)*res[max(i-1, 1), , drop=FALSE] +
                opt0$minimum*res[i, , drop=FALSE]
        }
    }
    ## Fix delta
    resopt[, 1] <- 2*sqrt(sum(op$d)*resopt[, length(op$d)+2]^2 +
                          sum(resopt[, 2:(sum(1-op$d)+1)]^2))
    oh <- up(resopt) # homoskedastic

    K <- ATTOptK(resopt, op$d)

    if (oh$delta==max(ep$delta) & nrow(res)>1) {
        if (extrasteps>0) {
            path_len <- nrow(op$ep)+extrasteps
            message("Increasing length of solution path to ", path_len)
            op <- ATTOptPath(path=op, maxsteps=path_len)
            return(ATTOptEstimate(op, C, opt.criterion, sigma2init, sigma2,
                           mvar, DM, alpha, beta, J, extrasteps))
        } else
            warning("Optimum found at end of path")
    }


    ## Robust variance, C=1 to keep bias the same
    or <- UpdatePath(oh, K, Cratio=1, sigma2, alpha, beta)
    if (missing(mvar))
        mvar <- nnMarginalVarO(op$D0, drop(resopt), op$d, op$y, sigma2,
                               tol=1e-12)
    ## Marginal variance could be negative in small samples
    ou <- UpdatePath(or, K, Cratio=1, sigma2, alpha, beta,
                     ucse=sqrt(or$sd^2+mvar))
    structure(list(e=c(unlist(oh), rsd=or$sd, rlower=or$lower, rupper=or$upper,
                       rhl=or$hl, rrmse=or$rmse, rmaxel=or$maxel, usd=ou$sd,
                       ulower=ou$lower, uupper=ou$upper, uhl=ou$hl, C=C),
              res=drop(resopt), k=drop(K), op=op),
              class="ATTEstimate")
}

#' @export
print.ATTEstimate <- function(x, digits = getOption("digits"), ...) {

    fmt <- function(x) format(x, digits=digits, width=digits+1)

    r <- data.frame(rbind(x$e[c("att", "maxbias", "rsd", "rhl")],
               c(x$e["att"], NA, unlist(x$e[c("usd", "uhl")]))))
    r <- cbind(r, l=r$att-r$rhl, u=r$att+r$rhl)
    r[, c("l", "u")] <- fmt(r[, c("l", "u")])

    r$ci <- paste0("(", r$l, ", ",  r$u, ")")
    r$rhl <- r$l <- r$u <- NULL

    colnames(r) <- c("Estimate", "Max. bias", "SE", "CI")
    rownames(r) <- c("CATT", "PATT")

    if ("M" %in% names(x$e))
        r$M <- c("M"=x$e["M"], NA)
    else
        r$delta <- c("delta"=x$e["delta"], NA)

    backup_options <- options(knitr.kable.NA = "", digits=digits)
    print(knitr::kable(r))
    options(backup_options)

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
#' @param op The output of \code{ATTOptPath}.
#' @param sigma2 estimate of the conditional variance of the outcome (assuming
#'     homoskedasticity). If not supplied, use homoskedastic variance estimate
#'     based on a nearest neighbor variance estimator.
#' @param DM distance matrix with dimension \code{n} by \code{n} to determine
#'     nearest neighbors when when estimating \code{sigma2}.
#' @param J number of nearest neighbors to use when estimating \code{sigma2}.
#' @return A list with two elements, \code{onesided} and \code{twosided}, for
#'     one- and two-sided efficiency.
#' @references \cite{Armstrong, T. B., and M. Kolesár (2018): Finite-Sample
#'     Optimal Estimation and Inference on Average Treatment Effects Under
#'     Unconfoundedness, \url{https://arxiv.org/abs/1712.04594}}
#' @examples
#' ## Use NSW experimental subsample with 25 treated and untreated units
#' dt <- NSWexper[c(1:25, 421:445), ]
#' Ahalf <- diag(c(0.15, 0.6, 2.5, 2.5, 2.5, 0.5, 0.5, 0.1, 0.1))
#' D0 <- distMat(dt[, 2:10], Ahalf, method="manhattan", dt$treated)
#' ## Distance matrix for variance estimation
#' DM <- distMat(dt[, 2:10], Ahalf, method="manhattan")
#' ## Compute the solution path, first 50 steps will be sufficient
#' op <- ATTOptPath(dt$re78, dt$treated, D0, maxsteps=50)
#' ATTEffBounds(op, C=1, DM=DM)
#' @export
ATTEffBounds <- function(op, C=1, beta=0.8, alpha=0.05, sigma2, J=3, DM) {
    if (missing(sigma2))
        sigma2 <- mean(nnvar(DM, op$d, op$y, J=J))
    if (length(sigma2)>1)
        warning(paste("The solution path is assuming homoskedastic variance,",
                      "but supplied sigma2 implies it's not homoskedastic"))
    n <- length(op$d)
    n1 <- sum(op$d)
    n0 <- n-n1
    del0 <- op$res[, "delta"]
    mu0 <- op$res[, "mu"]
    m0 <- op$res[, 2:(n0+1), drop=FALSE]
    r0 <- op$res[, (n0+2):(n+1), drop=FALSE]

    ## Modulus and its derivative when sigma2=1 and C=1
    mod11 <- function(del) {
        idx <- which.max(del0>=del)
        if (abs(del-del0[idx]) < .Machine$double.eps^(2/3))
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
    ## Minimum integrable point, make it > -100 to prevent numerical issues
    lbar <- zal-min(max(del0)/(2*sig)-1e-8, 100)
    if (integrand(lbar)>1e-6) {
        warning("Path too short to compute two-sided efficiency")
        num <- den <- NaN
    } else {
        lo <- -zal                          # lower endpoint
        while(integrand(lo)>1e-6) lo <- max(lo-2, lbar)
        num <- stats::integrate(integrand, lo, zal, abs.tol=1e-6)$value

        len <- function(del) {
            if (del==0 | del*2*sig> 1e5) # again, to prevent numerical issues
                return(Inf)
            m <- mod11(del)
            2*sig*cv((m$omega/(m$domega*sig)-del/sig)/2, alpha)*m$domega
        }
        idx <- min(which.min(vapply(del0, len, numeric(1))), length(del0)-1)
        den <- stats::optimize(len, interval=del0[c(idx-1, (idx+1))])$objective
    }
    list(onesided=eff1, twosided=num/den)
}
