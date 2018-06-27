#' Matrix of distances between observations
#'
#' Compute and return the distance matrix using the distance measure in
#' \code{method}.
#' @param x0,x1 Design matrices for treated and untreated units with
#'     dimensions \code{[n0 p]} and \code{[n1 p]} (if \eqn{p=1} and these are
#'     vectors, they are converted to a matrix)
#' @param Ahalf \eqn{A^{1/2}} weight matrix so that the distances are computed
#'     between \eqn{Ax_{0,j}} and \eqn{Ax_{1,i}}.
#' @param FullMatrix Indicator for whether to return the full
#'     \code{n0+n1}-by-\code{n0+n1} distance matrix, or only the
#'     \code{n1}-by-\code{n0} distance matrix between treated and untreated
#'     observations
#' @inheritParams stats::dist
#' @export
distMat <- function(x0, x1, Ahalf=diag(ncol(as.matrix(x1))),
                    method="euclidean", FullMatrix=FALSE) {
    X <- rbind(as.matrix(x0), as.matrix(x1))
    if (nrow(X)>500) {
        Dm <- matrix(nrow=nrow(X), ncol=nrow(X))
        ## retrieve method in case full name not given
        method <- attr(stats::dist(c(0, 1), method=method), "method")
        p <- if (method=="euclidean") 2L else if (method=="manhattan") 1L
        for (i in 1:nrow(X)) {
            Dm[i, ] <- colSums((Ahalf %*% (abs(X[i, ] - t(X))))^p)^(1/p)
        }
    } else {
        Dm <- unname(as.matrix(stats::dist(X, method=method)))
    }

    if (!FullMatrix)
        Dm <- Dm[(nrow(as.matrix(x0))+1):nrow(X), 1:nrow(as.matrix(x0))]
    Dm
}


#' Matching estimator for the ATT
#' @param M number of matches
#' @template D0
#' @param y0,y1 outcomes for treated and untreated units
ATTmatch <- function(M, y0, y1, D0, tol=1e-13) {
    w1 <- rep(1, length(y1))/length(y1)
    w0 <- rep(0, length(y0))
    for (i in 1:length(y1)) {
        ## Find NN of i
        idx <- D0[i, ] <= sort(D0[i, ])[M]+tol
        w0[idx] <- w0[idx] + 1/sum(idx)
    }
    w0 <- -w0/length(y1)
    list(w0=w0, w1=w1, att=sum(w1*y1)+sum(w0*y0))
}


#' Nearest-neighbor variance estimator
#' @template D0
#' @param d treatment status
#' @param y outcome
#' @param J number of nearest neighbors
#' @param tol distance tolerance
nnvar <- function(DM, d, y, J=3, tol=1e-13) {
    ehat2 <- vector(length=length(d))

    for (i in 1:nrow(DM)) {
        idx <- d==d[i]
## TODO:fix Ji adjustment
## idx[i] <- FALSE
        ## distance from i to other people with same treatment
        di <- sort(DM[i, idx])[J]
        idx <- DM[i, ]<=di+tol & idx     # in case of ties

        ehat2[i] <- (y[i]- mean(y[idx]))^2 #* Ji/(Ji+1)
    }

    ehat2
}


ATTEstimate <- function(res, n0, y) {
    n <- length(y)
    mean(y[((n0+1):n)]) -
        apply(res, 1, function(x) sum(x[2:(n0+1)]*y[1:n0])/sum(x[2:(n0+1)]))
}
