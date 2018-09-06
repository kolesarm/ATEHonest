#' Critical values for CIs based on a biased Gaussian estimator.
#'
#' Computes the critical value \eqn{cv_{1-alpha}(B)} that is needed to make the
#' confidence interval \eqn{X\pm cv} have coverage \eqn{1-alpha} if \eqn{X} is
#' normally distributed with variance one and maximum bias at most \eqn{B}.
#'
#' @param B Maximum bias, a non-negative number.
#' @param alpha Determines CI level, \code{1-alpha}. Needs to be between 0 and
#'     1. Can be a vector of values.
#' @return  Critical value
#' @examples
#' # 90% critical value:
#' CVb(B = 1, alpha = 0.1)
#' CVb(B = 1, alpha = c(0.05, 0.1))
#' @export
CVb <- function(B, alpha=0.05)
    sqrt(stats::qchisq(1-alpha, df = 1, ncp = B^2))
