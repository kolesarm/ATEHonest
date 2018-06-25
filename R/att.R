#' Check solution validity
#' @param m values for untreated observations
#' @param r values for treated observations
#' @inheritParams ATTh
ATTcheck <- function(m, r, mu, D0) {
    n <- nrow(D0) + ncol(D0)

    ## Check constraints hold, up to tolerance
    tol <- 10*n*.Machine$double.eps
    if (max(outer(r, m, "-")-D0) > tol) {
        warning("Lipschitz constraints do not hold")
        return(1)
    }
    ## Brute force solution should yield lower modulus
    modulus <- 2*(mean(r)+mu)
    delta2 <- 4*(sum(m^2)+nrow(D0)*mu^2)
    cvxs <- ATTbrute(delta2, D0)
    cvxmod <- cvxs[1]
    m1 <- cvxs[2:(ncol(D0)+1)]
    mu1 <- cvxs[n+2]
    cvxdelta2 <- 4*(sum(m1^2)+nrow(D0)*mu1^2)

    if (cvxmod-modulus > sqrt(max(cvxdelta2-delta2, 0))) {
        warning("CVX found higher modulus, ", cvxmod, ". Homotopy got ",
                modulus, ".\n Difference ", cvxmod-modulus)
        return(1)
    }

    ## TODO: solutions should also be close together

    0
}

#' Find modulus using CVX
#' @param delta2 square of delta for modulus
#' @param D0 matrix of distances with dimension \code{[n1 n0]} between untreated
#'     and treated units
ATTbrute <- function(delta2, D0) {
    n1 <- nrow(D0)
    n0 <- ncol(D0)

    mb <- CVXR::Variable(n0)
    rb <- CVXR::Variable(n1)
    mu <- CVXR::Variable(1)

    ob <- 2*CVXR::Maximize(mean(rb) + mu)
    ## outer(r, m, "-")<=D0 using kronecker
    con <- list(kronecker(t(rep(1, n0)), rb) -
                t(kronecker(t(rep(1, n1)), mb)) <= D0,
                4*(sum(mb^2)+n1*mu^2) <= delta2)
    s <- solve(CVXR::Problem(ob, con))
    c(s$value, s$getValue(mb), s$getValue(rb), s$getValue(mu))
}


#' Calculate next step in ATT homotopy
#' @param s list of state variables
#' @keywords internal
ATTstep <- function(s) {
    n0 <- ncol(s$D)
    n1 <- nrow(s$D)
    n <- n0+n1

    N <- cbind(rbind(Matrix::Matrix(0L, nrow=n0, ncol=n0), s$N0),
               Matrix::Matrix(0L, nrow=n, ncol=n1))
    cl <- igraph::components(igraph::graph.adjacency(N))
    le <- 1:cl$no
    ## avg number of matches per cluster
    am <- as.vector((table(factor(cl$membership[(n0+1):n], levels=le))) /
                    (table(factor(cl$membership[1:n0], levels=le))))

    ## Direction for m
    delta <- am[cl$membership[1:n0]]

    ## New distance matrix is D[1, ]+delta
    d1 <- (s$D-s$r0)/outer(apply(Matrix::t(s$N0) * delta, 2, max), delta, "-")
    d1[is.nan(d1) | d1<0] <- Inf
    if (!identical(s$minind, 0L))
        d1[s$minind] <- Inf

    ## Direction for Lambda
    ## DL <- Matrix::Matrix(0, nrow=n1, ncol=n0)
    ## DL[s$N0!=0] <- drop(MASS::ginv(as.matrix(A[, as.vector(s$N0)])) %*%
    ##                              c(delta, rep(1, n1)))

    DL <- Matrix::Matrix(0, nrow=n1, ncol=n0)
    ## R_k x M_k block of Lambda
    for (k in 1:cl$no) {
        if (cl$csize[k]>1) {
            ## Extract block of Lambda
            ind <- which(cl$membership==k)
            M0 <- ind[ind<=n0]
            M1 <- ind[ind>n0]-n0
            if (length(M0)==1) {
                DL[M1, M0] <- 1L
            } else if (length(M1)==1) {
                DL[M1, M0] <- 1/length(M0)
            } else {
                l0 <- length(M0)
                l1 <- length(M1)
                lN0 <- s$N0[M1, M0]
                A <- rbind(kronecker(Matrix::Diagonal(l0),
                                     Matrix::Matrix(1, nrow=1, ncol=l1)),
                           kronecker(Matrix::Matrix(1, nrow=1, ncol=l0),
                                     Matrix::Diagonal(l1)))
                A <- as.matrix(A[, as.vector(lN0)])
                lN0 <- as(lN0, "dgCMatrix")
                ## In same some of the sparse elements are zero
                lN0@x[lN0@x!=0] <- drop(MASS::ginv(A) %*%
                                        c(delta[M0], rep(1, l1)))
                DL[M1, M0] <- lN0
            }
        }
    }
    ## TODO: Add inequality constraint that DL>=0 if Lam==0

    if (min(DL)>=0) {
        d2 <- Inf
    } else {
        d2 <- -s$Lam/DL
        d2[DL>=0] <- Inf
    }
    s$DL <- DL
    print(c(min(d1), min(d2)))


    d <- min(d2, d1)

    ## TODO: this is probably not right
    if (cl$no==1) d <- Inf              # all connected

    ## Update
    s$m0 <- s$m0+d*delta
    s$Lam <- s$Lam+d*DL
    s$D <- t(t(s$D)+d*delta)
    s$r0 <- apply(s$D, 1, min)
    s$mu <- sum(s$m0) / n1

    if (min(d2) < min(d1)) {
        s$minind <- Matrix::which(d2==min(d2), arr.ind=TRUE)
        ## In case multiple drop, arrayInd(which(as.matrix(signif(d2, 10) ==
                                             ## min(signif(d2, 12)))),

        s$N0[s$minind] <- FALSE
        s$Lam[s$minind] <- 0L
        s$dmin <- 2L
    } else {
        ## s$N0[which.min(d1)] <- TRUE
        ## ## ## In case multiple join
        s$N0[which(d1==min(d1))] <- TRUE

        s$dmin <- 1L
    }
    s
}


#' Average treatment effect for the treated homotopy
#' @param D0 matrix of distances with dimension \code{[n1 n0]} between untreated
#'     and treated units
#' @param maxiter maximum number of steps in the homotopy
#' @param check check solution using CVX at each step
#' @return Returns solution path
ATTh <- function(D0, maxiter=50, check=FALSE) {
    n0 <- ncol(D0)
    n1 <- nrow(D0)
    n <- n0+n1

    ## Initialize state variables
    r0 <- apply(D0, 1, min)
    s <- list(m0=rep(0, n0),
              Lam=Matrix::Matrix(0, nrow=n1, ncol=n0),
              D=D0,
              r0=r0,
              N0=Matrix::Matrix(D0==r0),
              minind=0L,
              mu=0)
    res <- matrix(c(0, s$m0, r0, s$mu, NA), nrow=1)
    colnames(res) <- c("delta", 1:n, "mu", "drop")

    ## A <- rbind(kronecker(Matrix::Diagonal(n0),
    ##                      Matrix::Matrix(1, nrow=1, ncol=n1)),
    ##            kronecker(Matrix::Matrix(1, nrow=1, ncol=n0),
    ##                      Matrix::Diagonal(n1)))

    while (sum(s$m0^2)<Inf && nrow(res) <= maxiter) {
        ## work out directions
        s <- ATTstep(s)
        res <- rbind(res, c(2*sqrt(sum(s$m0^2)), s$m0, s$r0, s$mu, s$dmin))
        if (check && max(s$m0)<Inf && ATTcheck(s$m0, s$r0, s$mu, D0))
            stop("Solution doesn't agree with CVX")
    }
    res[-nrow(res), ]
}
