#' Critical values for CIs around a biased Gaussian estimator.
#'
#' Computes the critical value \eqn{cv_{1-alpha}(B)} that is needed to make the
#' confidence interval \eqn{X\pm cv_{1-alpha}(B)} have coverage \eqn{1-alpha} if
#' \eqn{X} is normally distributed with variance one and maximum bias at most
#' \eqn{B}.
#'
#' @param B Maximum bias, a non-negative vector.
#' @param alpha Determines CI level, \code{1-alpha}. Needs to be between 0 and
#'     1. Can be a vector of values.
#' @return  Critical value
#' @examples
#' # 90% critical value:
#' cv(B = 1, alpha = 0.1)
#' # 95% critical values
#' cv(B = c(0, 1, 3), alpha = 0.05)
#' # 96 and 90% critical values
#' cv(B = 1, alpha = c(0.05, 0.1))
#' @export
cv <- function(B, alpha=0.05)
    sqrt(stats::qchisq(1-alpha, df = 1, ncp = B^2))


#' Matrix of distances between observations
#'
#' Compute and return the distance matrix using the distance measure in
#' \code{method}.
#' @param X Design matrix for control variables with dimension \code{[n p]} (if
#'     \eqn{p=1} and \code{X} is a vector, then it is converted to a matrix)
#' @param d Vector of treatment indicators, needed only when
#'     \code{FullMatrix==FALSE}
#' @param Ahalf \eqn{A^{1/2}} weight matrix so that the distances are computed
#'     between \eqn{Ax_{0,j}} and \eqn{Ax_{1,i}}.
#' @param FullMatrix Indicator for whether to return the full
#'     \code{n0+n1}-by-\code{n0+n1} distance matrix, or only the
#'     \code{n1}-by-\code{n0} distance matrix between treated and untreated
#'     observations
#' @return \code{[n1 n0]} or \code{[n n]} matrix of distances.
#' @inheritParams stats::dist
#' @export
distMat <- function(X, d, Ahalf=diag(NCOL(X)),
                    method="euclidean", FullMatrix=FALSE, p=2) {
    if (class(X)!= "matrix")
        X <- as.matrix(X)

    if (nrow(X)>500) {
        Dm <- matrix(nrow=nrow(X), ncol=nrow(X))
        ## retrieve method in case full name not given
        method <- attr(stats::dist(c(0, 1), method=method), "method")
        p <- if (method=="euclidean") 2L
             else if (method=="manhattan") 1L
             else if (method=="minkowski") p
             else if (method=="maximum") Inf
             else stop("Method '", method, "' not yet implemented for large X")
        for (i in 1:nrow(X)) {
            Dm[i, ] <- if (p != Inf)
                           colSums((abs(Ahalf %*% (X[i, ] - t(X))))^p)^(1/p)
                       else
                           apply(abs(Ahalf %*% (X[i, ] - t(X))), 2, max)
        }
    } else {
        Dm <- unname(as.matrix(stats::dist(X, method=method)))
    }

    if (!FullMatrix)
        Dm <- Dm[d==1, d==0]
    Dm
}


#' Nearest-neighbor variance estimator
#'
#' Calculates an \code{n}-vector sigma2 of estimates of the variance of \code{y}
#' using a nearest-neighbor esitmator among observations with the same treatment
#' status \code{d}.
#' @param DM distance matrix with dimension \code{[n n]}.
#' @param d vector of treatment status indicators with length \code{n}
#' @param y outcome vector with length \code{n}
#' @param J number of nearest neighbors to average over
#' @param tol numerical tolerance for determining neighbors
#' @export
nnvar <- function(DM, d, y, J=3, tol=0) {
    ehat2 <- vector(length=length(d))

    for (i in 1:nrow(DM)) {
        ## distance from i to other people with same treatment, including
        ## oneself
        di <- sort(DM[i, d==d[i]])[J+1]
        idx <- DM[i, ]<=di+tol & d==d[i]     # in case of ties
        ehat2[i] <- (y[i]- mean(y[idx]))^2 * (sum(idx)+1)/sum(idx)
    }
    ehat2
}
