#' Compute matrix of distances between treated and untreated observations
#' @param x0,x1 Design matrices for treated and control variables with
#'     dimensions \code{[n0 p]} and \code{n1 p}
#' @param Ahalf \eqn{A^{1/2}} weight matrix
#' @inheritParams stats::dist
distMat <- function(x0, x1, Ahalf=diag(ncol(x1)), method="euclidean") {
    x0 <- x0 %*% Ahalf
    x1 <- x1 %*% Ahalf
    Dm <- unname(as.matrix(stats::dist(c(x0, x1), method=method)))
    Dm[(length(x0)+1):(length(x0)+length(x1)), 1:length(x0)]
}

#' Nearest-neighbor
#' @param DM distance matrix
#' @param d treatment status
#' @param y outcome
#' @param J number of nearest neighbors
nnvar <- function(DM, d, y, J=3) {
    ehat2 <- vector(length=length(d))

    for (i in 1:nrow(DM)) {
        idx <- d==d[i]
        idx[i] <- FALSE
        ## distance from i to other people with same treatment
        di <- sort(DM[i, idx])[J]
        Ji <- sum(DM[i, ]<=di & idx)     # in case of ties
        ehat2[i] <- (y[i]- mean(y[DM[i, ]<=di & idx]))^2 * Ji/(Ji+1)
    }

    ehat2
}


ATTEstimate <- function (res, n0, y) {
    n <- length(y)
    mean(y[((n0+1):n)]) -
        apply(res, 1, function(x) sum(x[2:(n0+1)]*y[1:n0])/sum(x[2:(n0+1)]))
}
