#' Critical values for inference based on a biased Gaussian estimator.
#'
#' Critical value \eqn{cv_{1-\alpha}(B)} such that the confidence interval
#' \eqn{X\pm cv_{1-\alpha}(B)}{X +- cv_{1-alpha}(B)} will have coverage
#' \eqn{1-\alpha}, where \eqn{X} is normally distributed with variance equal to
#' \eqn{1} and bias bounded by \eqn{B} in absolute value.
#' @param B Maximum bias, a non-negative vector.
#' @param alpha Scalar between 0 and 1 determining the confidence level,
#'     \code{1-alpha}
#' @return Critical value \eqn{cv_{1-\alpha}(B)}{cv_{1-alpha}(B)}
#' @examples
#' # 90% critical value:
#' cv(B = 1, alpha = 0.1)
#' # 95% critical values
#' cv(B = c(0, 1, 3), alpha = 0.05)
#' @export
cv <- function(B, alpha=0.05) {
    idx <- B<10 & !is.na(B)
    r <- B+stats::qnorm(1-alpha)
    r[idx] <- sqrt(stats::qchisq(1-alpha, df = 1, ncp = B[idx]^2))
    r
}

#' Matrix of distances between observations
#'
#' Compute a matrix of distances between \code{n} observations using the
#' distance measure in \code{method}.
#' @param X Design matrix of covariates with dimension \code{n} by \code{p}, or
#'     else a vector of length \code{n} if there is a single covariate.
#' @param d Vector of treatment indicators with length \code{n}. If supplied,
#'     return the \code{n1} by \code{n0} submatrix corresponding to distances
#'     between treated and untreated observations. Otherwise return the full
#'     \code{n} by \code{n} matrix
#' @param Ahalf Weight matrix with dimension \code{p} by \code{p} so that the
#'     distances are computed between \code{Ahalf \%*\%  X[i, ]}.
#' @return Matrix of distances with dimension \code{n} by \code{n} or else
#'     \code{n1} by \code{n0}
#' @inheritParams stats::dist
#' @examples
#' ## 4 units, unit 1 and 3 are treated.
#' distMat(X=c(1, 2, 3, 4), d=c(TRUE, FALSE, TRUE, FALSE))
#' @export
distMat <- function(X, Ahalf=diag(NCOL(X)),
                    method="euclidean", d=NULL, p=2) {
    X <- as.matrix(X) %*% t(Ahalf)

    if (nrow(X)>500) {
        Dm <- matrix(nrow=nrow(X), ncol=nrow(X))
        ## retrieve method in case full name not given
        method <- attr(stats::dist(c(0, 1), method=method), "method")
        p <- switch(method, euclidean=2L, manhattan=1L, minkowski=p,
                    maximum=Inf,
                    stop("Method '", method, "' not implemented for large X"))
        for (i in seq_len(nrow(X))) {
            Dm[i, ] <- if (p != Inf)
                           colSums(abs(X[i, ] - t(X))^p)^(1/p)
                       else
                           apply(abs(X[i, ] - t(X)), 2, max)
        }
    } else {
        Dm <- unname(as.matrix(stats::dist(X, method=method)))
    }

    if (!is.null(d))
        Dm <- Dm[d==1, d==0]
    Dm
}


## Nearest-neighbor variance estimator
## @param DM distance matrix with dimension \code{n} by \code{n}.
## @template data
## @param J number of nearest neighbors to average over
## @param tol numerical tolerance for determining nearest neighbors in
##     constructing matches
## @export
nnvar <- function(DM, d, y, J=3, tol=0) {
    ehat2 <- vector(length=length(d))

    for (i in seq_len(nrow(DM))) {
        ## distance from i to other people with same treatment, including
        ## oneself
        di <- sort(DM[i, d==d[i]])[J+1]
        idx <- DM[i, ]<=di+tol & d==d[i]     # in case of ties
        ehat2[i] <- (y[i]- mean(y[idx]))^2 * (sum(idx)+1)/sum(idx)
    }
    ehat2
}

## Matching Marginal variance estimate
nnMarginalVar <- function(D0, M, d, y, sigma2, tol) {
    if (length(sigma2)==1)
        sigma2 <- rep(sigma2, length(d))

    w0 <- rep(0, ncol(D0))
    Y0 <- vector(length=nrow(D0))

    for (i in seq_len(nrow(D0))) {
        ## Find NN of i, within tolerance
        idx <- D0[i, ] <= sort(D0[i, ])[min(M, ncol(D0))]+tol
        Y0[i] <- mean(y[!d][idx])
        w0[idx] <- w0[idx] + 1/sum(idx)^2
    }
    tau <- mean(y[d]-Y0)

    (sum((y[d]-Y0-tau)^2)-sum(sigma2[d])-sum(w0*sigma2[!d]))/nrow(D0)^2
}

## Marginal variance for optimal estimator
nnMarginalVarO <- function(D0, res, d, y, sigma2, tol=1e-12) {
    if (length(sigma2)==1)
        sigma2 <- rep(sigma2, length(d))

    n0 <- ncol(D0)
    n1 <- nrow(D0)
    ## Which constraints are binding
    m0 <- res[2:(n0+1)]
    r0 <- res[(n0+2):(n0+n1+1)]
    N0 <- Matrix::Matrix(outer(r0, m0, "-") >= D0-tol, sparse=TRUE)
    N <- cbind(rbind(Matrix::Matrix(0L, nrow=n0, ncol=n0), N0),
               Matrix::Matrix(0L, nrow=(n0+n1), ncol=n1))
    cl <- igraph::components(igraph::graph.adjacency(N))
    W <- solveLam(cl, N0, m0/res[(n0+n1+2)])

    Y0 <-  Matrix::drop(W %*% y[d==0])
    tau <- mean(y[d]-Y0)
    w0 <- Matrix::colSums(W^2)

    (sum((y[d]-Y0-tau)^2)-sum(sigma2[d])-sum(w0*sigma2[!d])) / n1^2

}
