#' Matching weights for untreated
#' @keywords internal
ATTMatchW <- function(D0, M, tol) {
    n0 <- ncol(D0)
    w0 <- rep(0, n0)
    if (M==Inf)
        return(-(w0+1)/n0)

    for (i in 1:nrow(D0)) {
        ## Find NN of i, within tolerance
        idx <- D0[i, ] <= sort(D0[i, ])[M]+tol
        w0[idx] <- w0[idx] + 1/sum(idx)
    }
    -w0/nrow(D0)
}


#' update estimation path with new C or new variance
#' @keywords internal
UpdatePath <- function(ep, resw, Cratio, sigma2, alpha, beta) {
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


#' Matching estimator for the ATT
#' @param M number of matches, vector. If \code{Inf}, then use simple difference
#'     in means
#' @template D0
#' @template data
#' @param tol numerical tolerance for determining neighbors
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


#' build optimal estimator given a solution path
#' @param mp Output of \code{ATTMatchPath}
#' @template data
#' @template D0
#' @param sigma2final vector of variance estimates with length{n} for
#'     determining standard error of the optimal estimators. In contrast,
#'     \code{sigma2} is used only for determining the optimal tuning parameter.
#' @param opt.criterion One of \code{"RMSE"}, \code{"OCI"}, \code{"FLCI"}
#' @param Mrange Range of Ms
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


#' Compute worst-case bias
#' @keywords internal
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
    f.con <- cbind(rep(1:length(f.rhs), each=2),
                   as.vector(t(expand.grid((n0+1):(n0+n1), 1:n0))),
                   rep(c(1, -1), length(f.rhs)))

    r <- lpSolve::lp("max", f.obj, , f.dir, f.rhs, dense.const=f.con)
    r$objval
}


#' Weights from solution path
#' @keywords internal
ATTOptW <- function(res, d) {
    n0 <- length(d)-sum(d)
    m <- res[, 2:(n0+1), drop=FALSE]
    resw <- matrix(1/sum(d), nrow=nrow(res), ncol=length(d))
    resw[, d==0] <- -m/rowSums(m)
    resw
}


#' build optimal estimator given a solution path
#' @param res output of \code{ATTh}
#' @template data
#' @export
ATTOptPath <- function(res, y, d, sigma2, C=1, alpha=0.05, beta=0.8) {
    n <- length(y)
    n0 <- n-sum(d)
    ## Vector or matrix?
    if (!is.matrix(res))
        stop("res needs to be a matrix")
    rmean <- rowMeans(res[, (n0+2):(n+1), drop=FALSE])
    m <- res[, 2:(n0+1), drop=FALSE]
    maxbias <- rmean - apply(m, 1, function(x) sum(x^2)) / rowSums(m)
    resw <- ATTOptW(res, d)

    UpdatePath(data.frame(att=drop(resw %*% y), maxbias=maxbias,
                          delta=unname(res[, 1]),
                          omega=unname(2*(res[, n+2]+rmean))),
               resw, C, sigma2, alpha, beta)
}

#' build optimal estimator given a solution path
#' @param res The \code{res} element of the output of \code{ATTh} on which to
#'     base the estimate
#' @param ep Output of \code{ATTOptPath} at \code{C=1}. This parameter is
#'     optional, if supplied, it will speed up the calculation.
#' @template data
#' @param sigma2final vector of variance estimates with length{n} for
#'     determining standard error of the optimal estimators. In contrast,
#'     \code{sigma2} is used only for determining the optimal tuning parameter.
#' @param opt.criterion One of \code{"RMSE"}, \code{"OCI"}, \code{"FLCI"}
#' @export
ATTOptEstimate <- function(res, ep=NULL, y, d, sigma2, C=1,
                           sigma2final=sigma2, alpha=0.05, beta=0.8,
                           opt.criterion="RMSE") {
    ## Drop delta=0, back out weights
    keep <- res[, "delta"]>0
    res <- res[keep, , drop=FALSE]

    ep <- if (is.null(ep)) {
              ATTOptPath(res, y, d, sigma2, C=C, alpha, beta)
          } else {
              UpdatePath(ep[keep, ], ATTOptW(res, d), C, sigma2, alpha, beta)
          }

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
            ATTOptPath((1-w)*res[i, , drop=FALSE]+w*res[i+1, , drop=FALSE],
                       y, d, sigma2, C, alpha, beta)[[idx]]
        opt1 <- stats::optimize(f1, interval=c(0, 1))
    } else {
        opt1 <- list(minimum=0, objective=min(ep[[idx]]))
    }

    if (i>1) {
        f0 <- function(w)
            ATTOptPath((1-w)*res[i-1, , drop=FALSE]+w*res[i, , drop=FALSE],
                       y, d, sigma2, C, alpha, beta)[[idx]]
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

    r1 <- ATTOptPath(resopt, y, d, sigma2, C, alpha, beta)
    resw1 <- ATTOptW(resopt, d)
    ## C=1 to kep bias the same
    r2 <- UpdatePath(r1, resw1, C=1, sigma2final, alpha, beta)

    if (r1$delta==max(ep$delta))
        warning("Optimum found at end of path")

    structure(list(e=cbind(r1, data.frame(rsd=r2$sd, rlower=r2$lower,
                                          rupper= r2$ upper, rhl=r2$hl,
                                          rrmse=r2$rmse, rmaxel=r2$maxel, C=C)),
                   res=resopt, w=resw1[, d==0]),
              class="ATTEstimate")
}
