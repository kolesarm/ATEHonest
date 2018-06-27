context("Check solution path for ATT in NSW")

test_that("Solution path for ATT in NSW", {

## main specification (based on 1980 census regressions)
X <- as.matrix(NSW[, 2:10])
Ahalf <- diag(c(0.15, 0.6, 2.5, 2.5, 2.5, 0.5, 0.5, 0.1, 0.1))
d <- NSW$treated
y <- NSW$re78
DM <- distMat(X, Ahalf, method="manhattan")


## Use same distance matrix for variance estimation, alternative is to use
## Mahalanobis, distmat(X, t(chol(solve(cov(X)))), p=1)
## sigma2 <- nnvar(DM, d, y, J=30)
Cvals <- 0.01


D0 <- DM[d==1, d==0]
r <- ATTh(D0, maxiter=70)
n1 <- nrow(D0)
n0 <- ncol(D0)

ATTEstimate(r$res, n0, y)

expect_equal(0, 0)

## compute nearest neighbor residuals for variance estimates
## uhat=nn_resid(Mvarest,x,d,outcome,dist_weight_mat,dist_p);
## standard deviation of errors in homoskedastic normal model for
## computing optimal kernels (will compute optimal kernels using this
## value and homoskedasticity, then compute confidence intervals that
## are valid under heteroskedasticity based on this kernel)
## sigmahat=sqrt(mean(uhat.^2));

## nonrobust_se = sigmahat*sqrt(sum((optweight.^2)));
## robust_se = sqrt(sum((optweight.^2).*(uhat.^2)));



})
