#' Matrix of distances between observations
#'
#' Compute and return the distance matrix using the distance measure in
#' \code{method}.
#' @param X Design matrix for control variables with dimension \code{[n p]} (if
#'     \eqn{p=1} and \code{X} is a vector, then it is converted to a matrix)
#' @param d Vector of treatment indicators
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
                    method="euclidean", FullMatrix=FALSE) {
    if (class(X)!= "matrix")
        X <- as.matrix(X)

    if (nrow(X)>500) {
        Dm <- matrix(nrow=nrow(X), ncol=nrow(X))
        ## retrieve method in case full name not given
        method <- attr(stats::dist(c(0, 1), method=method), "method")
        p <- if (method=="euclidean") 2L else if (method=="manhattan") 1L
        for (i in 1:nrow(X)) {
            Dm[i, ] <- colSums((abs(Ahalf %*% (X[i, ] - t(X))))^p)^(1/p)
        }
    } else {
        Dm <- unname(as.matrix(stats::dist(X, method=method)))
    }

    if (!FullMatrix)
        Dm <- Dm[d==1, d==0]
    Dm
}


#' Nearest-neighbor variance estimator
#' @param DM distance matrix
#' @param d treatment status
#' @param y outcome
#' @param J number of nearest neighbors
#' @param tol numerical tolerance for finding neighbors
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


#' Matching estimator for the ATT
#' @param M number of matches
#' @param y outcome vector
#' @param d vector of treatment indicators
#' @template D0
#' @param sigma2 vector of variances
ATTmatch <- function(M, y, d, D0, sigma2) {
    n1 <- nrow(D0)
    w <- rep(NA, length(y))
    w0 <- rep(0, ncol(D0))
    for (i in 1:sum(d)) {
        ## Find NN of i
        idx <- D0[i, ] <= sort(D0[i, ])[M]
        w0[idx] <- w0[idx] + 1/sum(idx)
    }
    w[d==0] <- -w0/n1
    w[d==1] <- 1/n1

    list(w=w, att=sum(w*y), sd=sqrt(sum(w^2*sigma2)))
}




ATTpath <- function(res, d, y, C=1, sigma2, alpha=0.05) {
    n <- length(y)
    n0 <- n-sum(d)
    if (length(sigma2)==1) sigma2 <- sigma2*rep(1, n)

    ms <- res[, 2:(n0+1)]
    rs <- res[, (n0+2): (n+1)]

    att <- mean(y[d==1]) - drop(ms %*% y[d==0]) / rowSums(ms)
    maxbias <- rowMeans(rs) - apply(ms, 1, function(x) sum(x^2)) / rowSums(ms)
    sd <- sqrt(mean(sigma2[d==1])/sum(d) +
               apply(ms, 1, function(x) sum(x^2*sigma2[d==0])) / rowSums(ms)^2)
    hl <- CVb(maxbias/sd, alpha)$cv * sd
    lower <- att - maxbias - stats::qnorm(1-alpha)*sd
    upper <- att + maxbias + stats::qnorm(1-alpha)*sd
    omega <- 2*(res[, "mu"]+rowMeans(rs))


    list(att=att, maxbias=C*maxbias, sd=sd, rmse=sqrt(sd^2+C*maxbias^2),
         lower=lower, upper=upper, hl=hl, delta=res[, "delta"],
         omega=omega)
}
