require("deSolve")

s.mle <- function(s.par,s.y,s.t,x0,y0){
  A <- sum((s.y - com.get_mu(s.par,s.t,x0,y0))^2 )
  A
}
s.mlecom1 <- function(s.par,s.y,s.t,x0,y0){
  A <- sum((s.y - com.get_mu1(s.par,s.t,x0,y0))^2 )
  A
}
com.get_mu <- function(par, times, x0,y0)
{
  par0 <- par;
  if (class(par0)!="list")
  {
    par0 <- list(
      r1 = par[1],
      a12 = par[2],
      k1 = par[3],
      r2 = par[4],
      a21 = par[5],
      k2 = par[6]);
  }
  
  state0 <- c(X=x0, Y=y0);
  y <- COMP.f( par0, state0, times );
  
  #if (class(y)=="try-error" )
  #  return(rep(NaN, length(times)*2));
  index.time <- 1:length(times)
  return ( c(y[index.time, 2:3] ) );
}


COMP.f <-function( parameters, state, times ){
  Lorenz<-function(t, state, parameters) 
  {
    with( as.list(c(state, parameters)),
          {
            dX <- r1*X*(1-((X+a12*Y)/k1))
            dY <- r2*Y*(1-((a21*X+Y)/k2))
            
            list(c(dX, dY))
          }
          
    ) # end with(as.list ...
  }
  
  out <- try(rk(y = state, times = times, func = Lorenz, parms = parameters,method="rk4") );
  out;
}

s.mle1 <- function(s.par,s.y,s.t,x0,y0,Nmax,Nmin1,Nmin2){
  A <- sum((s.y - com.get_mu1(s.par,s.t,x0,y0,Nmax,Nmin1,Nmin2))^2 )
  A
}

com.get_mu1 <- function(par, times,x0,y0)
{
  par0 <- par;
  if (class(par0)!="list")
  {
    par0 <- list(
      r1 = par[1],
      m1 = par[2],
      a12= par[3],
      k1 = par[4],
      n1 = par[5],
      r2 = par[6],
      m2 = par[7],
      a21= par[8],
      k2 = par[9],
      n2 = par[10]);
  }
  
  state0 <- c(X=x0, Y=y0);
  y <- COMP.f1( par0, state0, times );
  index.time <- 1:length(times)
  return ( c(y[index.time, 2:3] ) );
}


COMP.f1 <-function( parameters, state, times ){
  Lorenz<-function(t, state, parameters) 
  {
    with( as.list(c(state, parameters)),
          {
            dX <- r1*X*((1-((X+a12*Y)/k1)))^m1*(1-(8.52/X)^n1)
            dY <- r2*Y*((1-((a21*X+Y)/k2)))^m2*(1-(8.52/Y)^n2)
            
            list(c(dX, dY))
          }
          
    ) # end with(as.list ...
  }
  
  out <- try(rk(y = state, times = times, func = Lorenz, parms = parameters,method="rk4") );
  out;
}


LC.get_mu <- function(par, times, options=list())
{
  A <- par[1]/(1+par[2]*exp(-par[3]*times))
  return (A);
}


LC.mle <- function(par,y,times){
  
  sum((y-LC.get_mu(par,times))^2)
}




curve.get_est_param<-function(dat1,pheno1,pheno2)
{
  phenos <- cbind(pheno1,pheno2)
  ym <- as.numeric(colMeans(phenos))
  value <- c()
  allpar <- list()
  for(ii in 1:10){
    s.par <- c(0.1585815,4.5978148,125.2705139,0.1233422,2.8023220,91.2431834)*runif(6,0.85,1.15)
    res <- optim(s.par,s.mle,s.y=ym,s.t=dat1$sample_times,x0=ym[1],y0=ym[16],
                 method="BFGS",control=list(trace=F,maxit=10000))
    value <- c(value,res$value)
    allpar[[ii]] <- res$par
  }
  p1 <- allpar[[which(value==min(value))]]
  parin1 <- c(runif(1),sd(pheno1),runif(1),sd(pheno2),runif(1))
  parin <- c(parin1,p1)
  r0<- optim( parin, curve.mlefunc, y = phenos, time.std=dat1$sample_times,x0=ym[1],y0=ym[16],
              method ="BFGS",control=list(maxit=2000,trace=F))
  
  cat("Estimated parameters:", "log(L)=", r0$value, "PAR=", r0$par, "\n");
  return(c(r0$par,r0$value));
}

curve.mlefunc<-function( par, y, time.std,x0,y0)
{
  len.cov <- 5;
  par.covar <- par[1:len.cov];
  sig <- SAD3.get_mat(par.covar, times=1:length(time.std),traits = 2);
  
  m  <- length(y[1,]);
  n  <- length(y[,1]);
  
  par1 <- par[-c(1:len.cov)]
  mu0 <- com.get_mu(par1,time.std,x0,y0)
  
  fy0<- dmvnorm(y,mu0,sig)
  A <- -sum(log(fy0));
  
  return (A);
}
