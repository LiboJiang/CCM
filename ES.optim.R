

ES.H1.map <- function(dat,pheno,fsnp,H0,select.snp){
  
  ng <- dim(fsnp[,select.snp])[2]
  phenos <- as.matrix( pheno )
  allparin <- matrix(NA,ng,50)
  for(i in select.snp){
    snp.type <- names(table(as.character(unlist(c(fsnp[,i])))))
    g.par <- c()
    snp.index <- list()
    for(j in 1:length(snp.type)){
      index <- which(fsnp[,i]==snp.type[j])
      yy <- as.numeric(colMeans(phenos[index,]))
      #r1 <- optim(H0[6:11],s.mle,s.t=dat1$sample_times,s.y=yy,x0=yy[1],y0=yy[16],method="BFGS",
                  #control=list(maxit=2000,trace=T))
      snp.index[[j]] <- index
      g.par <- c(g.par,H0[6:11])
    }
    parin <- c(H0[1:5],g.par)
    
    res <- try(optim( parin, ES.mlefunc, y = phenos, time.std=dat1$sample_times, snp.index=snp.index,method ="BFGS",
                  control=list(maxit=5000,trace=F)),TRUE)
    if(class(res)=="try-error"){
      tmp <- NA
      cat("Estimated parameters:","SNP=",i,"LR=",NA, "\n");
      allparin[i,1:length(tmp)] <- tmp
    }else{
      LR <- 2*(H0[12]-res$value)
      tmp <- c(LR,res$value,res$par)
      cat("Estimated parameters:","SNP=",i,"LR=",LR, "log(L)=", res$value, "PAR=", res$par, "\n");
      allparin[i,1:length(tmp)] <- tmp
    }
    

  }
  return(allparin);
}










ES.mlefunc<-function( par, y, time.std,snp.index )
{
  len.cov <- 5;
  par.covar <- par[1:len.cov];
  sig <- SAD3.get_mat(par.covar, time.std,traits=2);
  
  m  <- length(y[1,]);
  n  <- length(y[,1]);
  len <- 0
  A <- 0
  for(k in 1:length(snp.index)){
    gen_par<- par[(len.cov+1+len):(len.cov+ 6+len)];
    yy0 <- y[snp.index[[k]],]
    myy0 <- as.numeric(colMeans(yy0))
    mu0 <- com.get_mu( gen_par, time.std,x0=myy0[1],y0=myy0[16]);
    fy0 <- dmvnorm(yy0,mu0,sig)
    A <- A -sum(log(fy0));
    len <- len + 6
  }
  return (A);
}