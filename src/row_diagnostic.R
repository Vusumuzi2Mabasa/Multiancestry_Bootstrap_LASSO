suppressPackageStartupMessages(library(jsonlite))
args <- commandArgs(TRUE);root <- args[1];source(file.path(root,'src/lasso.R'))
input <- fromJSON(file.path(root,'work/row_input.json'),simplifyVector=FALSE)
out <- list();set.seed(512711)
for(c in input) {
  matrices <- lapply(c$rows,function(x) list(g=do.call(rbind,lapply(x$g,unlist)),y=unlist(x$y)))
  sufficient <- function(g,y) {
    X <- scale(g,TRUE,FALSE);yc <- y-mean(y);d <- colSums(X^2)
    list(A=crossprod(X)/sqrt(outer(d,d)),z=as.numeric(crossprod(X,yc))/sqrt(d))
  }
  contexts <- lapply(matrices,function(x)sufficient(x$g,x$y))
  B <- 500L;p<-8L;row_selected<-matrix(FALSE,p,B);par_selected<-row_selected
  for(k in seq_along(matrices)) {
    data<-matrices[[k]];co<-contexts[[k]]
    Z<-co$z+t(chol(co$A))%*%matrix(rnorm(p*B),p,B)
    par_selected<-par_selected | abs(solve_lasso(co$A,Z,1.5)$beta)>1e-8
    for(b in seq_len(B)) {
      ix<-sample.int(nrow(data$g),replace=TRUE)
      s<-sufficient(data$g[ix,],data$y[ix])
      row_selected[,b]<-row_selected[,b] | abs(solve_lasso(s$A,matrix(s$z),1.5)$beta[,1])>1e-8
    }
  }
  for(j in seq_len(p))out[[length(out)+1L]]<-data.frame(region=c$region,scenario=c$scenario,index=j-1L,
    row_frequency=mean(row_selected[j,]),parametric_frequency=mean(par_selected[j,]),B=B,
    variance='known_noise_1',lambda=1.5)
}
write.table(do.call(rbind,out),file.path(root,'results/row_bootstrap_diagnostic.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
