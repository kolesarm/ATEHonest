## Matching weights for the untreated
ATTMatchW <- function(D0, M, tol) {
    n0 <- ncol(D0)
    w0 <- rep(0, n0)
    if (M==Inf)
        return(rep(-1/n0, n0))

    for (i in seq_len(nrow(D0))) {
        ## Find NN of i, within tolerance
        idx <- D0[i, ] <= sort(D0[i, ])[M]+tol
        w0[idx] <- w0[idx] + 1/sum(idx)
    }
    -w0/nrow(D0)
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
    f.con <- cbind(rep(seq_len(n1*n0), each=2),
                   as.vector(t(expand.grid((n0+1):(n0+n1), 1:n0))),
                   rep(c(1, -1), length(f.rhs)))

    r <- lpSolve::lp("max", f.obj, , f.dir, f.rhs, dense.const=f.con) #nolint
    r$objval
}


#' Compute the matching estimator for the ATT
#'
#' Computes the matching estimator and the matching weights for a range of
#' matches \code{M}. The output of this function is used as an input for
#' \code{\link{ATTMatchEstimate}} for inference on the CATT.
#' @param M a vector of integers determining the number of matches. If
#'     \code{Inf}, then use the simple difference in means estimator.
#' @template D0
#' @template data
#' @param tol numerical tolerance for determining nearest neighbors in
#'     constructing matches
#' @return List with the following components
#' \describe{
#'
#' \item{ep}{A data frame with columns \code{M}, \code{maxbias}, and \code{att},
#' corresponding to the number of matches, the scaled worst-case bias, and the
#' CATT estimate.}
#'
#' \item{K}{A matrix where each row \code{j} corresponds to the linear weights
#' \eqn{k} used to form the matching estimator with \code{M[j]} matches.}
#'
#' \item{d}{Vector of treatment indicators, as supplied by \code{d}}
#' }
#' @examples
#' Ahalf <- diag(c(0.15, 0.6, 2.5, 2.5, 2.5, 0.5, 0.5, 0.1, 0.1))
#' D0 <- distMat(NSWexper[, 2:10], Ahalf, method="manhattan", NSWexper$treated)
#' mp <- ATTMatchPath(NSWexper$re78, NSWexper$treated, D0, M=1:2, tol=1e-12)
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
         K=resw, d=d)
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
#' @param sigma2final vector of variance estimates with length \code{n} for
#'     determining the standard error of the optimal estimator. In contrast,
#'     \code{sigma2} is used only for determining the optimal tuning parameter.
#' @param opt.criterion One of \code{"RMSE"} (root mean squared error),
#'     \code{"OCI"} (one-sided confidence intervals), \code{"FLCI"}
#'     (fixed-length two-sided confidence intervals)
#' @return Returns an object of class \code{"ATTEstimate"}. An object of class
#'     \code{"ATTEstimate"} is a list containing the following components:
#' \describe{
#' \item{e}{Data frame with columns TODO}
#' \item{w}{weights}
#' }
#' @examples
#' Ahalf <- diag(c(0.15, 0.6, 2.5, 2.5, 2.5, 0.5, 0.5, 0.1, 0.1))
#' D0 <- distMat(NSWexper[, 2:10], Ahalf, method="manhattan", NSWexper$treated)
#' mp <- ATTMatchPath(NSWexper$re78, NSWexper$treated, D0, M=c(1, 2), tol=1e-12)
#' ## Distance matrix for variance estimation
#' DM <- distMat(NSWexper[, 2:10], Ahalf, method="manhattan")
#' sigma2 <- nnvar(DM, NSWexper$treated, NSWexper$re78, J=3)
#' ## Estimator and CI based on a single match is better than with 2 matches for
#' ## RMSE
#' ATTMatchEstimate(mp, mean(sigma2), C=1, sigma2final=sigma2)
#' @export
ATTMatchEstimate <- function(mp, sigma2, C=1, sigma2final=sigma2, alpha=0.05,
                             beta=0.8, opt.criterion="RMSE") {
    ## Update estimate path with new value of C
    mp$ep <- UpdatePath(mp$ep, mp$K, C, sigma2, alpha, beta)

    ## Index of criterion to optimize
    colidx <- switch(opt.criterion, RMSE = "rmse", OCI = "maxel",
                     FLCI = "hl")
    idx <- which.max(names(mp$ep)==colidx)
    i <- which.min(mp$ep[[idx]])
    if (i==nrow(mp$ep) & nrow(mp$ep)>1)
        warning("Optimum found at end of path")

    ## Robust se, C=1 to keep bias the same
    er <- UpdatePath(mp$ep[i, ], mp$K[i, ], Cratio=1, sigma2final, alpha,
                     beta)

    structure(list(e=cbind(mp$ep[i, ],
                           data.frame(rsd=er$sd, rlower=er$lower,
                                      rupper=er$upper, rhl=er$hl, rrmse=er$rmse,
                                      rmaxel=er$maxel, C=C)),
                   w=mp$K[i, mp$d==0]), class="ATTEstimate")
}
