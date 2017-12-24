


funmap.curve<-function(dat1,pheno){
  ym <- as.numeric(colMeans(pheno))
  value <- c()
  allpar <- list()
  for(ii in 1:10){
    par <- c(24.6554536,0.9656415,0.1592335)*runif(3,0.8,1.2)
    res <- optim(par,LC.mle,y=ym,times=dat1$sample_times,method="BFGS",control=list(trace=T,maxit=10000))
    value <- c(value,res$value)
    allpar[[ii]] <- res$par
  }
  p1 <- allpar[[which(value==min(value))]]
  parin1 <- c(0.5,sd(pheno[,8]))
  parin <- c(parin1,p1)
  r0<- optim( parin, fun.curve.mlefunc, y = pheno, time.std=dat1$sample_times,
              method ="BFGS",control=list(maxit=2000,trace=F))
  
  cat("Estimated parameters:", "log(L)=", r0$value, "PAR=", r0$par, "\n");
  return(c(r0$par,r0$value));
}



fun.curve.mlefunc<-function( par, y, time.std)
{
  len.cov <- 2;
  par.covar <- par[1:len.cov];
  sig <- SAD1.get_mat(par.covar, times=time.std,traits = 1);
  
  par1 <- par[-c(1:len.cov)]
  mu0 <- LC.get_mu(par=par1,times=time.std)
  
  fy0<- dmvnorm(y,mu0,sig)
  A <- -sum(log(fy0));
  
  return (A);
}

