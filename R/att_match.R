## Matching weights for the untreated
ATTMatchK0 <- function(D0, M, tol) {
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

    r <- lpSolve::lp("max", f.obj, , f.dir, f.rhs, dense.const=f.con)
    r$objval
}


#' Matching estimator for the ATT
#'
#' Computes the matching estimator and the matching weights for a range of
#' matches \code{M}. The output of this function is used as an input for
#' \code{\link{ATTMatchEstimate}} for inference on the conditional average
#' treatment effect for the treated (CATT) and population average treatment
#' effect for the treated (PATT).
#' @param M a vector of integers determining the number of matches. If
#'     \code{Inf}, then use the simple difference in means estimator.
#' @template D0
#' @template data
#' @param tol numerical tolerance for determining nearest neighbors in
#'     constructing matches
#' @return List with the following components: \describe{
#'
#' \item{ep}{A data frame with columns \code{M}, \code{maxbias}, \code{att}, and
#' \code{lindw} corresponding to the number of matches, the scaled worst-case
#' bias, the ATT estimate, and the largest Lindeberg weight.}
#'
#' \item{K}{A matrix where each row \code{j} corresponds to the linear weights
#' \eqn{k} used to form the matching estimator with \code{M[j]} matches.}
#'
#' \item{d}{Vector of treatment indicators, as supplied by \code{d}}
#'
#' \item{y}{Vector of outcomes, as supplied by \code{y}}
#'
#' \item{tol}{The tolerance parameter \code{tol}, as supplied by \code{tol}}
#'
#' \item{D0}{The distance matrix, as supplied by \code{D0}}
#' }
#' @examples
#' ## Construct distance matrix
#' Ahalf <- diag(c(0.15, 0.6, 2.5, 2.5, 2.5, 0.5, 0.5, 0.1, 0.1))
#' D0 <- distMat(NSWexper[, 2:10], Ahalf, method="manhattan", NSWexper$treated)
#' mp <- ATTMatchPath(NSWexper$re78, NSWexper$treated, D0, M=1:2, tol=1e-12)
#' @export
ATTMatchPath <- function(y, d, D0, M=1:25, tol=1e-12) {
    n1 <- nrow(D0)
    maxbias <- vector(length=length(M))
    K <- matrix(nrow=length(M), ncol=(ncol(D0)+n1))
    K[, d==1] <- 1/n1

    for (j in seq_along(M)) {
        K[j, d==0] <- ATTMatchK0(D0, M[j], tol)
        maxbias[j] <- ATTbias(K[j, d==0], D0)
    }
    list(ep=data.frame(att=drop(K %*% y), maxbias=maxbias, M=M,
                       lindw=apply(K^2, 1, max)/rowSums(K^2)),
         K=K, d=d, y=y, tol=tol, D0=D0)
}


#' Inference on the CATT and the PATT using the matching estimator
#'
#' Computes the matching estimator and confidence intervals (CIs) for the
#' conditional average treatment effect for the treated (CATT) and population
#' average treatment effect for the treated (PATT). If \code{ATTMatchPath} used
#' a single \code{M}, the estimator and CIs are based on a matching estimator
#' with this number of matches. Otherwise, optimize the number of matches
#' according to \code{opt.criterion}.
#' @param mp output of \code{ATTMatchPath}
#' @param C Lipschitz smoothness constant
#' @param opt.criterion criterion to optimize. One of \code{"RMSE"} (root mean
#'     squared error), \code{"OCI"} (one-sided confidence intervals),
#'     \code{"FLCI"} (fixed-length two-sided confidence intervals)
#' @param alpha determines confidence level, \code{1-alpha}.
#' @param beta quantile \code{beta} of excess length for determining performance
#'     of one-sided CIs.
#' @param sigma2init estimate of the conditional variance of the outcome, used
#'     to optimize the number of matches. If not supplied, use homoskedastic
#'     variance estimate based on a nearest neighbor variance estimator.
#' @param sigma2 vector of variance estimates with length \code{n} for
#'     determining the conditional standard error of the optimal estimator. In
#'     contrast, \code{sigma2init} is used only for determining the optimal
#'     tuning parameter. If not supplied, use the nearest neighbor variance
#'     estimator.
#' @param mvar Marginal variance estimate (variance of the CATT) used to
#'     construct CIs for the PATT. If not supplied use the matching estimator of
#'     Abadie and Imbens (2006).
#' @param DM distance matrix with dimension \code{n} by \code{n} to determine
#'     nearest neighbors when when estimating \code{sigma2init} and
#'     \code{sigma2}.
#' @param J number of nearest neighbors to use when estimating \code{sigma2init}
#'     and \code{sigma2}.
#' @references \cite{Abadie, A. and G. W. Imbens (2006):
#'     "Large sample properties of matching estimators for
#'      average treatment effects,"
#'     Econometrica, 74, 235–267.}
#'
#' \cite{Armstrong, T. B., and M. Kolesár (2020): Finite-Sample
#'     Optimal Estimation and Inference on Average Treatment Effects Under
#'     Unconfoundedness, \url{https://arxiv.org/abs/1712.04594}}
#' @examples
#' Ahalf <- diag(c(0.15, 0.6, 2.5, 2.5, 2.5, 0.5, 0.5, 0.1, 0.1))
#' D0 <- distMat(NSWexper[, 2:10], Ahalf, method="manhattan", NSWexper$treated)
#' mp <- ATTMatchPath(NSWexper$re78, NSWexper$treated, D0, M=c(1, 2), tol=1e-12)
#' ## Distance matrix for variance estimation
#' DM <- distMat(NSWexper[, 2:10], Ahalf, method="manhattan")
#' ## Estimator based on a single match is better than with 2 matches for RMSE
#' ATTMatchEstimate(mp, C=1, DM=DM)
#' @return Returns an object of class \code{"ATTEstimate"}. An object of class
#'     \code{"ATTEstimate"} is a list containing the following components:
#'     \describe{
#'
#'     \item{e}{Vector with elements \code{"att"} (value of ATT estimate),
#'     \code{"maxbias"} (worst-case bias), \code{"M"} (number of matches),
#'     \code{"lindw"} (maximal Lindeberg weight \eqn{Lind(k)}), \code{"sd"},
#'     \code{"hl"}, \code{"lower"}, \code{"upper"}, \code{"maxel"},
#'     \code{"rmse"} (standard deviation, half-length of two-sided CI, lower and
#'     upper endpoint of one-sided CIs, worst-case excess length of one-sided CI
#'     at quantile \code{beta} and RMSE of estimator, assuming conditional
#'     variance equals \code{sigma2inint}) \code{"rsd"}, \code{"rlower"},
#'     \code{"rupper"}, \code{"rhl"}, \code{"rrmse"}, \code{"rmaxel"} (same
#'     quantities, but calculated using \code{sigma2}), \code{"usd"},
#'     \code{"ulower"}, \code{"uupper"}, \code{"uhl"} (standard deviation,
#'     endpoints for one-sided CIs and half-length of two-sided CI for PATE),
#'     and \code{"C"} (Value of Lipschitz constant \code{C})}
#'
#'     \item{k}{Vector of matching weights \eqn{k(x_i, d_i)}}
#'     }
#' @export
ATTMatchEstimate <- function(mp, C=1, opt.criterion="RMSE", sigma2init,
                             sigma2, mvar, DM, alpha=0.05, beta=0.8, J=3) {
    if (missing(sigma2))
        sigma2 <- nnvar(DM, mp$d, mp$y, J=J)
    if (missing(sigma2init))
        sigma2init <- mean(sigma2)

    ## Update estimate path with new value of C
    mp$ep <- UpdatePath(mp$ep, mp$K, C, sigma2init, alpha, beta)

    ## Index of criterion to optimize
    idx <- which.max(names(mp$ep) ==
                     switch(opt.criterion, RMSE = "rmse", OCI = "maxel",
                            FLCI = "hl"))
    i <- which.min(mp$ep[[idx]])
    if (i==nrow(mp$ep) & nrow(mp$ep)>1)
        warning("Optimum found at end of path")

    ## Robust se, C=1 to keep bias the same
    er <- UpdatePath(mp$ep[i, ], mp$K[i, , drop=FALSE], Cratio=1, sigma2,
                     alpha, beta)
    ## SE for PATE
    if (missing(mvar))
        mvar <- nnMarginalVar(mp$D0, er$M, mp$d, mp$y, sigma2, mp$tol)

    eu <- UpdatePath(mp$ep[i, ], mp$K[i, , drop=FALSE], Cratio=1, sigma2,
                     alpha, beta, ucse=sqrt(er$sd^2+mvar))
    structure(list(e=c(unlist(mp$ep[i, ]), rsd=er$sd, rlower=er$lower,
                       rupper=er$upper, rhl=er$hl, rrmse=er$rmse,
                       rmaxel=er$maxel, usd=eu$sd, ulower=eu$lower,
                       uupper=eu$upper, uhl=eu$hl, C=C),
                   k=mp$K[i, ]), class="ATTEstimate")
}
