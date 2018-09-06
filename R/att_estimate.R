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


#' Matching estimator for the ATT
#' @param M number of matches
#' @param y outcome vector
#' @param d vector of treatment indicators
#' @template D0
#' @param sigma2 vector of variances
#' @export
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

#' build optimal estimator given a solution path
#' @param res output of \code{ATTh}
#' @param y outcome vector
#' @param d vector of treatment indicators
#' @param C smoothness constant
#' @param sigma2 vector of variances
#' @param alpha CI coverage
#' @param beta quantile of excess length
#' @export
ATTEstimatePath <- function(res, d, y, C=1, sigma2,
                            alpha=0.05, beta=0.8) {
    n <- length(y)
    n0 <- n-sum(d)
    if (length(sigma2)==1) sigma2 <- sigma2*rep(1, n)

    if (length(dim(res)) > 1L) {
        m <- res[, 2:(n0+1)]
        r <- res[, (n0+2):(n+1)]
        mu <- res[, "mu"]
        delta <- res[, "delta"]
        att <- mean(y[d==1]) - drop(m %*% y[d==0]) / rowSums(m)
        maxbias <- rowMeans(r) - apply(m, 1, function(x) sum(x^2)) / rowSums(m)
        sd <- sqrt(mean(sigma2[d==1])/sum(d) + apply(m, 1, function(x)
            sum(x^2*sigma2[d==0])) / rowSums(m)^2)
        omega <- 2*(mu+rowMeans(r))
    } else {
        ## It's a vector
        m <- res[2:(n0+1)]
        r <- res[(n0+2):(n+1)]
        mu <- res["mu"]
        delta <- res["delta"]
        att <- mean(y[d==1]) - sum(m * y[d==0]) / sum(m)
        maxbias <- mean(r) -  sum(m^2) / sum(m)
        sd <- sqrt(mean(sigma2[d==1])/sum(d) + sum(m^2*sigma2[d==0]) / sum(m)^2)
        omega <- 2*(mu+mean(r))
    }
    maxbias <- C*maxbias
    hl <- CVb(maxbias/sd, alpha)$cv * sd
    lower <- att - maxbias - stats::qnorm(1-alpha)*sd
    upper <- att + maxbias + stats::qnorm(1-alpha)*sd
    ## worst-case quantile of excess length
    maxel <- 2*maxbias + sd * (stats::qnorm(1-alpha)+stats::qnorm(beta))

    data.frame(att=att, maxbias=maxbias, sd=sd, lower=lower, upper=upper, hl=hl,
               rmse=sqrt(sd^2+maxbias^2), maxel=maxel, omega=unname(omega),
               delta=unname(delta))
}

#' build optimal estimator given a solution path
#' @param res output of \code{ATTh}
#' @param ep output of \code{ATTEstimatePath}
#' @param y outcome vector
#' @param d vector of treatment indicators
#' @param C smoothness constant
#' @param sigma2 vector of variances
#' @param sigma2final vector of variances for final variance.
#' @param alpha CI coverage
#' @param beta quantile of excess length
#' @param opt.criterion One of \code{"RMSE"}, \code{"OCI"},  \code{"FLCI"}
#' @param UpdateC Update C that's assumed in \code{ep}?
#' @export
ATTOptEstimate <- function(res, ep, d, y, C=1, sigma2,
                           sigma2final=sigma2, alpha=0.05, beta=0.8,
                           opt.criterion="RMSE", UpdateC=TRUE) {
    ## Update estimate path with new value of C
    ep$maxbias <- C*ep$maxbias
    ep[, c("rmse", "lower", "upper", "hl", "maxel")] <-
        cbind(sqrt(ep$sd^2+ep$maxbias^2),
              ep$att - ep$maxbias - stats::qnorm(1-alpha)*ep$sd,
              ep$att + ep$maxbias + stats::qnorm(1-alpha)*ep$sd,
              CVb(ep$maxbias/ep$sd, alpha)$cv * ep$sd,
              2*ep$maxbias + ep$sd * (stats::qnorm(1-alpha)+stats::qnorm(beta)))

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
            ATTEstimatePath((1-w)*res[i, ]+w*res[i+1, ], d, y, C,
                            sigma2, alpha, beta)[[idx]]
        opt1 <- stats::optimize(f1, interval=c(0, 1))
    } else {
        opt1 <- list(minimum=0, objective=min(ep[[idx]]))
    }

    if (i>1) {
        f0 <- function(w)
            ATTEstimatePath((1-w)*res[i-1, ]+w*res[i, ], d, y, C,
                            sigma2, alpha, beta)[[idx]]
        opt0 <- stats::optimize(f0, interval=c(0, 1))
    } else {
        opt0 <- list(minimum=1, objective=min(ep[[idx]]))
    }

    if (opt1$objective < opt0$objective) {
        resopt <- (1-opt1$minimum)*res[i, ] +
            opt1$minimum*res[min(i+1, nrow(res)), ]
    } else {
        resopt <- (1-opt0$minimum)*res[max(i-1, 1), ]+opt0$minimum*res[i, ]
    }

    r1 <- ATTEstimatePath(resopt, d, y, C, sigma2, alpha, beta)
    r2 <- ATTEstimatePath(resopt, d, y, C, sigma2final, alpha, beta)
    if (r1$delta==max(ep$delta))
        warning("Optimum found at end of path")
    cbind(r1, data.frame(rsd=r2$sd, rlower=r2$lower, rupper=r2$upper,
         rhl=r2$hl, rrmse=r2$rmse, rmaxel=r2$maxel, C=C))
}
