## Check solution validity using CVX
## @param m values for untreated observations
## @param r values for treated observations
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
    cvx_modulus <- cvxs[length(cvxs)]

    if (cvx_modulus-modulus > sqrt(max(cvxs[1]^2-delta2, 0))) {
        warning("CVX found higher modulus, ", cvx_modulus, ". Homotopy got ",
                modulus, ".\n Difference ", cvx_modulus-modulus, " at delta=",
                sqrt(delta2))
        return(1)
    }

    ## Solutions should also be close together
    diff <- max(abs(c(cvxs[(ncol(D0)+2):(n+1)]-r, cvxs[2:(ncol(D0)+1)]-m)))
    if (diff > 1e-1)
        message("CVX solution differs from homotopy, by ", round(diff, 3),
                " at delta=", sqrt(delta2), ":\n", paste(cvxs, collapse="; "))

    0
}

## Find modulus using CVX
## @param delta2 square of delta for modulus
## @param D0 matrix of distances with dimension \code{[n1 n0]} between untreated
##     and treated units
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
    s <- CVXR::solve(CVXR::Problem(ob, con))
    c(delta=2*sqrt(sum(s$getValue(mb)^2)+nrow(D0)*s$getValue(mu)^2),
      s$getValue(mb), s$getValue(rb), mu=s$getValue(mu), omega=s$value)
}


## Calculate next step in ATT homotopy
## @param s list of state variables
## @return Updated list \code{s} of state variables
ATTstep <- function(s, tol=.Machine$double.eps*n0*n1) {
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

    ## Direction for m and r
    delta <- am[cl$membership[1:n0]]
    deltar <- am[cl$membership[(n0+1):(n0+n1)]]

    den <- outer(deltar, delta, "-")
    d1 <- (s$D-s$r0)/den
    d1[is.nan(d1) | den<=0] <- Inf

    DL <- Matrix::Matrix(0, nrow=n1, ncol=n0)
    ## R_k x M_k block of Lambda
    for (k in 1:cl$no) {
        if (cl$csize[k]>1) {
            ## Extract block of Lambda
            ind <- which(cl$membership==k)
            M0 <- ind[ind<=n0]
            M1 <- ind[ind>n0]-n0
            l0 <- length(M0)
            l1 <- length(M1)

            if (l0==1) {
                DL[M1, M0] <- 1L
            } else if (l1==1) {
                DL[M1, M0] <- 1/length(M0)
            } else {
                lN0 <- s$N0[M1, M0]
                A <- rbind(kronecker(Matrix::Diagonal(l0),
                                     Matrix::Matrix(1, nrow=1, ncol=l1)),
                           kronecker(Matrix::Matrix(1, nrow=1, ncol=l0),
                                     Matrix::Diagonal(l1)))
                A <- as.matrix(A[, as.vector(lN0)])
                lN0 <- methods::as(lN0, "dgCMatrix")
                ## If some of the sparse elements are zero
                lN0@x[lN0@x!=0] <- drop(MASS::ginv(A) %*%
                                        c(delta[M0], rep(1, l1)))
                DL[M1, M0] <- lN0
            }
        }
    }
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
    if(d>0) {
        s$m0 <- s$m0+d*delta
        s$Lam <- s$Lam+d*DL
        s$D <- t(t(s$D)+d*delta)
        s$r0 <- s$r0+d*deltar
        s$mu <- sum(s$m0) / n1
    }

    if (min(d2) < min(d1)) {
        ## Sometimes there is numerical rounding error
        minind <- Matrix::which(d2<=min(d2)+tol, arr.ind=TRUE)
        s$N0[minind] <- FALSE
        s$Lam[minind] <- 0L
        s$drop <- TRUE
    } else {
        ## In case multiple join
        s$N0[which(d1<=min(d1)+tol)] <- TRUE
        s$drop <- FALSE
    }
    s
}


#' Homotopy for average treatment effect for the treated
#'
#' Calculates optimal weights \eqn{m} and \eqn{r} for the control and treated
#' observations as a function of \eqn{\delta}, or equivalently \eqn{\mu}, using
#' the algorithm described in the appendix to Armstrong and Kolesár (2018)
#' @template D0
#' @param maxiter maximum number of steps in the homotopy. If the homotopy has
#'     less steps than \code{maxiter}, returns the whole solution path.
#' @param check check at each step that the solution matches that in the
#'     \code{\link[CVXR]{CVXR-package}} (generic convex optimizer package).
#' @param tol numerical tolerance for rounding error when finding the nearest
#'     neighbor. All observations with effective distance within \code{tol} of
#'     the closest are considered to be active.
#' @param h Optionally, supply previous output of \code{ATTh}. If not provided,
#'     the homotopy is started at the beginning. If provided, it starts at the
#'     step where the previous call to \code{ATTh} ended.
#' @return A list with two elements:
#'
#' \describe{
#'
#'   \item{res}{A matrix with rows corresponding to steps in the homotopy, so
#'   that the maximum number of rows is \code{maxiter} (if homotopy started at
#'   the beginning), and columns corresponding to \eqn{\delta}, \eqn{m},
#'   \eqn{r}, \eqn{\mu}, and \code{drop}, an indicator if an observations has
#'   been dropped from an active set, or added }
#'
#'   \item{m0}{A vector of length \code{n0} of corresponding to \eqn{m} at the
#'   last step.}
#'
#'   \item{r0}{A vector of length \code{n1} of corresponding to \eqn{r} at the
#'   last step.}
#'
#'   \item{mu}{A scalar corresponding to \eqn{\mu} at the last step.}
#'
#'   \item{D}{A matrix of effective distances with dimension \code{[n1 n0]} at
#'   the last step.}
#'
#'   \item{Lam}{A sparse matrix of Lagrange multipliers with dimension
#'             \code{[n1 n0]} at the last step.}
#'
#'   \item{N0}{A sparse matrix of nearest neighbors with dimension
#'             \code{[n1 n0]} at the last step.}
#' }
#' @examples
#' x0 <- c(0, 1, 2, 3)
#' x1 <- c(1, 4, 5)
#' d <- c(rep(FALSE, length(x0)), rep(TRUE, length(x1)))
#' D0 <- distMat(c(x0, x1), d=d)
#' ## Compute first three steps
#' h <- ATTh(D0, maxiter=3)
#' ## Compute the remaining steps, checking them against CVX solution
#' h2 <- ATTh(D0, h=h, check=TRUE)
#' ## Check
#' x0 <- c(1,  1,  3,  3,  4,  6, 17, 17, 21, 22, 23, 23, 26, 30, 36, 37,
#'         47, 53, 58, 61)
#' x1 <- c(9, 15, 26, 27, 28, 31, 32, 40, 51, 52, 56)
#' d <- c(rep(FALSE, length(x0)), rep(TRUE, length(x1)))
#' D0 <- distMat(c(x0, x1), d=d)
#' h <- ATTh(D0, check=TRUE, maxiter=20)
#' h$res
#' @references \cite{Armstrong, T. B., and M. Kolesár (2018): Finite-Sample
#'     Optimal Estimation and Inference on Average Treatment Effects Under
#'     Unconfoundedness, Unpublished manuscript}
#' @export
ATTh <- function(D0, h, maxiter=50, check=FALSE,
                 tol=.Machine$double.eps*ncol(D0)*nrow(D0)) {
    n0 <- ncol(D0)
    n1 <- nrow(D0)
    n <- n0+n1

    ## Initialize state variables
    if (missing(h)) {
        r0 <- apply(D0, 1, min)
        h <- list(m0=rep(0, n0), Lam=Matrix::Matrix(0, nrow=n1, ncol=n0), D=D0,
                  r0=r0, N0=Matrix::Matrix(D0<=r0), # without tol
                  mu=0, drop=NA, res=matrix(nrow=0, ncol=n+3))
        colnames(h$res) <- c("delta", 1:n, "mu", "drop")
    }

    stopme <- FALSE
    while (sum(h$m0^2) < Inf && nrow(h$res) < maxiter && !isTRUE(stopme)) {
        if (check && ATTcheck(h$m0, h$r0, h$mu, D0))
            stop("Solution doesn't agree with CVX")
        h$res <- rbind(h$res, c(2*sqrt(n1*h$mu^2 + sum(h$m0^2)), h$m0,
                                h$r0, h$mu, h$drop))
        stopme <- tryCatch({
            h <- ATTstep(h, tol)
        }, error = function(e) { # nolint
            message(conditionMessage(e))
            cat("Stopping ATTh at step ", nrow(h$res), "\n")
            return(TRUE)
        })
    }
    h
}
