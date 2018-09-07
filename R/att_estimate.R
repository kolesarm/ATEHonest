#' Matching estimator for the ATT
#' @param M number of matches
#' @param y outcome vector
#' @param d vector of treatment indicators
#' @param C smoothness constant
#' @param alpha CI coverage
#' @template D0
#' @param sigma2 vector of variances
#' @export
ATTMatchPath <- function(M, y, d, D0, sigma2, C, alpha=0.05, beta=0.8) {
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

    att <- sum(w*y)
    maxbias <- C*ATTbias(w[d==0], D0)
    sd <- sqrt(sum(w^2*sigma2))
    hl <- cv(maxbias/sd, alpha) * sd
    lower <- att - maxbias - stats::qnorm(1-alpha)*sd
    upper <- att + maxbias + stats::qnorm(1-alpha)*sd
    ## worst-case quantile of excess length
    maxel <- 2*maxbias + sd * (stats::qnorm(1-alpha)+stats::qnorm(beta))

    data.frame(att=att, maxbias=maxbias, sd=sd, lower=lower, upper=upper, hl=hl,
               rmse=sqrt(sd^2+maxbias^2), maxel=maxel, M=M)

}

#' Compute worst-case bias
#' @keywords internal
ATTbias <- function(w, D0) {
    ## If w receives zero weight, we can drop it
    idx <- w != 0
    w0 <- w[idx]
    D0 <- D0[, idx]

    n1 <- nrow(D0)
    n0 <- ncol(D0)

    ## CVX method
    ## mb <- CVXR::Variable(n0)
    ## rb <- CVXR::Variable(n1)
    ## ob <- CVXR::Maximize(mean(rb) + sum(w0*mb))
    ## ## outer(r, m, "-")<=D0 using kronecker
    ## con <- list(kronecker(t(rep(1, n0)), rb) -
    ##             t(kronecker(t(rep(1, n1)), mb)) <= D0)
    ## s <- solve(CVXR::Problem(ob, con))
    ## return(s$value)

    ## Use LP
    f.obj <- c(w0, rep(1, n1)/n1)
    f.rhs <- as.vector(D0)
    f.dir <- rep("<=", length(f.rhs))
    ## Constraint number, column number, value
    f.con.d <- cbind(rep(1:length(f.rhs), each=2),
                      as.vector(t(expand.grid((n0+1):(n0+n1), 1:n0))),
                      rep(c(1, -1), length(f.rhs)))

    r <- lpSolve::lp("max", f.obj, , f.dir, f.rhs, dense.const=f.con.d)
    r$objval
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
ATTOptPath <- function(res, d, y, C=1, sigma2,
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
    hl <- cv(maxbias/sd, alpha)$cv * sd
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
              cv(ep$maxbias/ep$sd, alpha)$cv * ep$sd,
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
