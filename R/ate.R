#' Check solution validity
#' @param m values for untreated observations
#' @param r values for treated observations
#' @inheritParams AUOh
AUOcheck <- function(m, r, D0) {
    n <- nrow(D0) + ncol(D0)

    ## Check constraints hold, up to tolerance
    tol <- 10*n*.Machine$double.eps
    if (max(outer(r, m, "-")-D0) > tol) {
        warning("Lipschitz constraints do not hold")
        return(1)
    }
    ## Brute force solution should yield lower modulus
    modulus <- 2*sum(c(m, r))/n
    cvxmod <- AUObrute(2*sqrt(sum(m^2)), D0)[1]
    if (modulus < cvxmod) {
        warning("CVX found higher modulus, ", cvxmod, ". Homotopy got ",
                modulus, ".\n Difference ", cvxmod-modulus)
        if (cvxmod - modulus < 1e-10 ) {
            return(0)
        } else {
            return(1)
        }
    }

    ## TODO: solutions should also be close together

    0
}

#' Find modulus using CVX
#' @param delta delta for modulus
#' @param D0 matrix of distances with dimension \code{[n1 n0]} between untreated
#'     and treated units
AUObrute <- function(delta, D0) {
    n1 <- nrow(D0)
    n0 <- ncol(D0)
    n <- n0+n1

    mb <- CVXR::Variable(n0)
    rb <- CVXR::Variable(n1)

    ob <- 2*CVXR::Maximize(sum(mb)+sum(rb))/n
    ## outer(r, m, "-")<=D0 using kronecker
    con <- list(kronecker(t(rep(1, n0)), rb) -
                t(kronecker(t(rep(1, n1)), mb)) <= D0,
                CVXR::p_norm(mb, p=2) <= delta/2)
    s <- solve(CVXR::Problem(ob, con))
    c(s$value, s$getValue(mb), s$getValue(rb))
}


#' Calculate next step in AUO homotopy
#' @param s list of state variables
#' @keywords internal
AUOstep <- function(s, A) {
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
    delta <- 1+am[cl$membership[1:n0]]
    ## New distance matrix is D[1, ]+delta
    d1 <- (s$D-s$d0)/outer(apply(Matrix::t(s$N0) * delta, 2, max), delta, "-")
    d1[is.nan(d1) | d1<0] <- Inf
    if (!identical(s$minind, 0L))
        d1[s$minind] <- Inf

    ## Direction for Lambda
    DL <- Matrix::Matrix(0, nrow=n1, ncol=n0)
    DL[s$N0!=0] <- drop(MASS::ginv(as.matrix(A[, as.vector(s$N0)])) %*%
                                 c(delta-1, rep(1, n1)))
    ## TODO: Add inequality constraint that DL>=0 if Lam==0

    if (min(DL)>=0) {
        d2 <- Inf
    } else {
        d2 <- -s$Lam/DL
        d2[DL>=0] <- Inf
    }

    d <- min(d2, d1)

    ## TODO: this is probably not right
    if (cl$no==1) d <- Inf              # all connected

    ## Update
    s$m0 <- s$m0+d*delta
    s$Lam <- s$Lam+d*DL
    s$D <- t(t(s$D)+d*delta)
    s$d0 <- apply(s$D, 1, min)

    if (min(d2)<min(d1)) {
        s$minind <- Matrix::which(d2==min(d2), arr.ind=TRUE)
        ## In case multiple drop, arrayInd(which(as.matrix(signif(d2, 10) ==
                                             ## min(signif(d2, 12)))),

        s$N0[s$minind] <- FALSE
        s$Lam[s$minind] <- 0L
        s$dmin <- 2L
    } else {
        s$N0[which.min(d1)] <- TRUE
        ## In case multiple join, it should be the case that min(d1)=0 in the
        ## next step, rather than s$N0 <- Matrix::Matrix(s$D==s$d0), which has
        ## numerical issues

        s$dmin <- 1L
    }
    s
}


#' Average untreated effect homotopy
#' @param D0 matrix of distances with dimension \code{[n1 n0]} between untreated
#'     and treated units
#' @param maxiter maximum number of steps in the homotopy
#' @param check check solution using CVX at each step
#' @return Returns solution path
AUOh <- function(D0, maxiter=50, check=FALSE) {
    n0 <- ncol(D0)
    n1 <- nrow(D0)
    n <- n0+n1

    ## Initialize state variables
    d0 <- apply(D0, 1, min)
    s <- list(m0=rep(0, n0),
              Lam=Matrix::Matrix(0, nrow=n1, ncol=n0),
              D=D0,
              d0=d0,
              N0=Matrix::Matrix(D0==d0),
              minind=0L)
    res <- matrix(c(0, s$m0, d0, NA), nrow=1)
    colnames(res) <- c("delta", 1:n, "drop")

    A <- rbind(kronecker(Matrix::Diagonal(n0),
                         Matrix::Matrix(1, nrow=1, ncol=n1)),
               kronecker(Matrix::Matrix(1, nrow=1, ncol=n0),
                         Matrix::Diagonal(n1)))

    while (sum(s$m0^2)<Inf && nrow(res) <= maxiter) {
        ## work out directions
        s <- AUOstep(s, A)
        res <- rbind(res, c(2*sqrt(sum(s$m0^2)), s$m0, s$d0, s$dmin))
        if (check && max(s$m0)<Inf && AUOcheck(s$m0, s$d0, D0))
            stop("Solution doesn't agree with CVX")
    }
    res[-nrow(res), ]
}
