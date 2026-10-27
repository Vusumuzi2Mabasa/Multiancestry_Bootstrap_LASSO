suppressPackageStartupMessages(library(jsonlite)); suppressPackageStartupMessages(library(glmnet))
args <- commandArgs(TRUE); root <- args[1]; source(file.path(root,'src/lasso.R'))
set.seed(731593); errors <- numeric(); kkt <- numeric()
for(p in c(2,5,8)) for(rho in c(0,.5,.9)) for(lambda in c(.5,1.5,3)) {
  A <- outer(1:p,1:p,function(i,j)rho^abs(i-j)); z <- rnorm(p)*3
  fit <- solve_lasso(A,matrix(z,p,1),lambda,tol=1e-10)
  # X'X/m=A, X'y/m=z, with no centering or implicit standardization.
  X <- sqrt(p)*chol(A); y <- solve(t(X),p*z)
  g <- glmnet(X,y,alpha=1,lambda=lambda,standardize=FALSE,intercept=FALSE,thresh=1e-18,maxit=100000)
  errors <- c(errors,max(abs(as.numeric(coef(g))[-1]-fit$beta)))
  kkt <- c(kkt,fit$kkt)
}
print(max(errors));stopifnot(max(errors)<1e-5)
z <- c(-4,-.5,0,2,5); f <- solve_lasso(diag(5),matrix(z),1.5)
stopifnot(max(abs(f$beta-soft(z,1.5)))<1e-12)
A <- outer(1:8,1:8,function(i,j).7^abs(i-j)); Z <- matrix(rnorm(800),8,100)
f <- solve_lasso(A,Z,1.5); signs <- c(-1,1,-1,1,1,1,-1,1)
flip <- solve_lasso(A*outer(signs,signs),Z*signs,1.5)
stopifnot(max(abs(f$beta*signs-flip$beta))<1e-10)
per <- sample(8); perm <- solve_lasso(A[per,per],Z[per,],1.5)
stopifnot(max(abs(f$beta[per,]-perm$beta))<1e-5)
E <- t(chol(A))%*%matrix(rnorm(8*50000),8)
cov_error <- max(abs(tcrossprod(E)/ncol(E)-A));stopifnot(cov_error<.025)
# Direct individual fixed-design residual bootstrap has the specified score law.
X <- matrix(rnorm(1200*8),1200,8)%*%chol(A);X <- scale(X,TRUE,FALSE)
d <- colSums(X^2);R <- crossprod(X)/sqrt(outer(d,d))
Zres <- crossprod(X,matrix(rnorm(1200*2000),1200,2000))/sqrt(d)
res_cov_error <- max(abs(tcrossprod(Zres)/2000-R));stopifnot(res_cov_error<.10)
write_json(list(passed=TRUE,glmnet_max_difference=max(errors),max_kkt=max(kkt),
 covariance_max_difference=cov_error,individual_residual_covariance_difference=res_cov_error,
 R=as.character(getRversion()),glmnet=as.character(packageVersion('glmnet'))),
 file.path(root,'results/lasso_controls.json'),pretty=TRUE,auto_unbox=TRUE,digits=16)
