# Ordinary LASSO, independently fitted population coefficients; no group penalty.
soft <- function(x, lambda) sign(x)*pmax(abs(x)-lambda,0)
solve_lasso <- function(A,Z,lambda,tol=1e-7,maxit=10000L) {
  p <- nrow(A); Z <- as.matrix(Z)
  if(p==0) return(list(beta=Z,kkt=0,sweeps=0L))
  T <- matrix(0,p,ncol(Z)); E <- Z
  for(it in seq_len(maxit)) {
    change <- 0
    for(j in seq_len(p)) {
      new <- soft(E[j,]+A[j,j]*T[j,],lambda)/A[j,j]
      delta <- new-T[j,]; T[j,] <- new
      E <- E - outer(A[,j],delta); change <- max(change,max(abs(delta)))
    }
    if(change<tol) break
  }
  G <- A%*%T-Z
  violation <- ifelse(abs(T)>1e-8,abs(G+lambda*sign(T)),pmax(abs(G)-lambda,0))
  kkt <- max(violation)
  if(!is.finite(kkt)||kkt>1e-6||it==maxit) stop('LASSO KKT/iteration failure')
  list(beta=T,kkt=kkt,sweeps=it)
}
bootstrap_lasso <- function(context,B,lambdas,seed) {
  set.seed(seed); K <- length(context); p <- length(context[[1]]$z)
  if(p==0) return(list(frequencies=matrix(numeric(),0,length(lambdas)),kkt=0))
  draws <- lapply(context,function(c) {
    A <- do.call(rbind,lapply(c$A,unlist)); z <- unlist(c$z)
    list(A=A,Z=z+t(chol(A))%*%matrix(rnorm(p*B),p,B))
  })
  F <- matrix(0,p,length(lambdas)); worst <- 0
  for(l in seq_along(lambdas)) {
    selected <- matrix(FALSE,p,B)
    for(k in seq_len(K)) {
      fit <- solve_lasso(draws[[k]]$A,draws[[k]]$Z,lambdas[l])
      selected <- selected | abs(fit$beta)>1e-8; worst <- max(worst,fit$kkt)
    }
    F[,l] <- rowMeans(selected)
  }
  list(frequencies=F,kkt=worst)
}
if(sys.nframe()==0L) {
  suppressPackageStartupMessages(library(jsonlite))
  a <- commandArgs(TRUE); input <- fromJSON(a[1],simplifyVector=FALSE)
  result <- lapply(input$cases,function(c) {
    started <- proc.time()[3]
    f <- bootstrap_lasso(c$context,input$B,unlist(input$lambdas),c$seed)
    list(id=c$id,frequencies=unname(f$frequencies),kkt=f$kkt,seconds=unname(proc.time()[3]-started))
  })
  write_json(result,a[2],auto_unbox=TRUE,digits=16,pretty=FALSE)
}
