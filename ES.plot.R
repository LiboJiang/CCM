ES.plot.ind <- function(dat){
  
  times <- dat$sample_times
  N <- dat$sample_N
  
  pdf("Figure-all-ind.pdf",height=8,width=8)
  par(fig=c(0,0.5,0.5,1))
  par(mar=c(2,3.4,4.2,0.2),las=1)
  plot(times,dat$ep.pi[1,],type="l",lwd=2,col="#757575",xlab=" ",
       ylab="Abundance",cex.lab=1.5,mgp = c(2.2, 1, 0),xaxt="n",yaxt="n",
       xlim=c(times[1],times[15]),ylim=c(min(dat$ep.pi)-3,max(dat$ep.pi)))
  axis(2,c(10,15,20,25,30),c(10,15,20,25,30),lwd=2,cex.axis=1.2)
  axis(1,times,times,lwd=2,cex.axis=1.2)
  for(i in 2:N){
    lines(times,dat$ep.pi[i,],lwd=2,col="#757575")
  }
  mtext(expression(italic("E. coli")),side=3, line=0.4,cex=2.0,adj=0.5)
  mtext("A", side=3, line=0.4,cex=2.2,adj=-0.2)
  par(fig=c(0.5,1,0.5,1),new=TRUE)
  plot(times,dat$sp.pi[1,],type="l",lwd=2,col="#757575",xlab=" ",
       ylab=" ",cex.lab=1.5,mgp = c(2.2, 1, 0),xaxt="n",yaxt="n",
       xlim=c(times[1],times[15]),ylim=c(min(dat$sp.pi)-3,max(dat$sp.pi)))
  axis(2,c(10,15,20,25,30),c(10,15,20,25,30),lwd=2,cex.axis=1.2)
  axis(1,times,times,lwd=2,cex.axis=1.2)
  for(i in 2:N){
    lines(times,dat$sp.pi[i,],lwd=2,col="#757575")
  }
  mtext(expression(italic("S. aureus")),side=3, line=0.4,cex=2.0,adj=0.5)
  
  par(fig=c(0,0.5,0,0.5),new=TRUE)
  par(mar=c(4,3.4,2.2,0.2),las=1)
  plot(times,dat$ep.p[1,],type="l",lwd=2,col="#757575",xlab="Time (h) ",
       ylab="Abundance",cex.lab=1.5,mgp = c(2.2, 1, 0),xaxt="n",yaxt="n",
       xlim=c(times[1],times[15]),ylim=c(min(dat$ep.p)-3,max(dat$ep.p)))
  axis(2,c(10,15,20,25,30),c(10,15,20,25,30),lwd=2,cex.axis=1.2)
  axis(1,times,times,lwd=2,cex.axis=1.2)
  for(i in 2:N){
    lines(times,dat$ep.p[i,],lwd=2,col="#757575")
  }
  mtext("B", side=3, line=0.4,cex=2.2,adj=-0.2)
  par(fig=c(0.5,1,0,0.5),new=TRUE)
  plot(times,dat$sp.p[1,],type="l",lwd=2,col="#757575",xlab="Time (h) ",
       ylab=" ",cex.lab=1.5,mgp = c(2.2, 1, 0),xaxt="n",yaxt="n",
       xlim=c(times[1],times[15]),ylim=c(min(dat$sp.p)-3,max(dat$sp.p)))
  axis(2,c(10,15,20,25,30),c(10,15,20,25,30),lwd=2,cex.axis=1.2)
  axis(1,times,times,lwd=2,cex.axis=1.2)
  for(i in 2:N){
    lines(times,dat$sp.p[i,],lwd=2,col="#757575")
  }
  dev.off()
}


ES.qq.plot <- function(dat){
  
  
  Asp <- qqnorm(as.numeric(unlist(c(dat$sp.p))),plot.it=F)
  Aep <- qqnorm(as.numeric(unlist(c(dat$ep.p))),plot.it=F)
  Asp.i <- qqnorm(as.numeric(unlist(c(dat$sp.pi))),plot.it=F)
  Aep.i <- qqnorm(as.numeric(unlist(c(dat$ep.pi))),plot.it=F)
  
  
  pdf("Figure-normal-test.pdf",height=8,width=8)
  par(fig=c(0,0.5,0.5,1))
  par(mar=c(2,3.4,4.2,0.2),las=1)
  plot(Aep.i$x,Aep.i$y,pch=16,lwd=1,col="#757575",xlab=" ",
       ylab="Sample Quantiles",cex.lab=1.5,mgp = c(2.2, 1, 0),xaxt="n",yaxt="n")
  axis(1,seq(-3,3,1),seq(-3,3,1),lwd=2,cex.axis=1.2)
  axis(2,seq(10,30,5),seq(10,30,5),lwd=2,cex.axis=1.2)
  qqline(as.numeric(unlist(c(dat$ep.pi))),col=2,lwd=2)
  
  mtext(expression(italic("E. coli")),side=3, line=0.4,cex=2.0,adj=0.5)
  mtext("A", side=3, line=0.4,cex=2.2,adj=-0.2)
  
  par(fig=c(0.5,1,0.5,1),new=TRUE)
  plot(Asp.i$x,Asp.i$y,pch=16,lwd=1,col="#757575",xlab=" ",
       ylab=" ",cex.lab=1.5,mgp = c(2.2, 1, 0),xaxt="n",yaxt="n")
  axis(1,seq(-3,3,1),seq(-3,3,1),lwd=2,cex.axis=1.2)
  axis(2,seq(10,25,5),seq(10,25,5),lwd=2,cex.axis=1.2)
  qqline(as.numeric(unlist(c(dat$sp.pi))),col=2,lwd=2)
  
  mtext(expression(italic("S. aureus")),side=3, line=0.4,cex=2.0,adj=0.5)
  par(fig=c(0,0.5,0,0.5),new=TRUE)
  par(mar=c(4,3.4,2.2,0.2),las=1)
  plot(Aep$x,Aep$y,pch=16,lwd=1,col="#757575",xlab=" Theoretial Quantiles",
       ylab="Sample Quantiles",cex.lab=1.5,mgp = c(2.2, 1, 0),xaxt="n",yaxt="n")
  axis(1,seq(-3,3,1),seq(-3,3,1),lwd=2,cex.axis=1.2)
  axis(2,seq(15,30,5),seq(15,30,5),lwd=2,cex.axis=1.2)
  qqline(as.numeric(unlist(c(dat$ep.p))),col=2,lwd=2)
  mtext("B", side=3, line=0.4,cex=2.2,adj=-0.2)
  par(fig=c(0.5,1,0,0.5),new=TRUE)
  plot(Asp$x,Asp$y,pch=16,lwd=1,col="#757575",xlab=" Theoretial Quantiles",
       ylab=" ",cex.lab=1.5,mgp = c(2.2, 1, 0),xaxt="n",yaxt="n")
  axis(1,seq(-3,3,1),seq(-3,3,1),lwd=2,cex.axis=1.2)
  axis(2,seq(12,24,4),seq(12,24,4),lwd=2,cex.axis=1.2)
  qqline(as.numeric(unlist(c(dat$sp.p))),col=2,lwd=2)
  dev.off()
  
  
}


ES.ind.adjust.plot <- function(alldat){
  
  mape <- c(alldat$estr)
  maps <- c(alldat$sstr)
  eresii <- c(alldat$eresii)
  eresi <- c(alldat$eresi)
  sresii <- c(alldat$sresii)
  sresi <- c(alldat$sresi)
  pdf("Figure-adjust-residual.pdf",height=8,width=8)
  par(fig=c(0,0.5,0.5,1))
  par(mar=c(2,3.4,4.2,0.2),las=1)
  plot(mape,eresii,pch=16,lwd=1,col="#757575",xlab=" ",
       ylab="Residuals",cex.lab=1.5,mgp = c(2.2, 1, 0),xaxt="n",yaxt="n",xlim=c(0.5,5.5))
  axis(1,seq(1,5,1),seq(1,5,1),lwd=2,cex.axis=1.2)
  axis(2,seq(round(min(eresii),0),round(max(eresii),0),2),
       seq(round(min(eresii),0),round(max(eresii),0),2),lwd=2,cex.axis=1.2)
  abline(h=0,lwd=2,lty=2)
  
  mtext(expression(italic("E. coli")),side=3, line=0.4,cex=2.0,adj=0.5)
  mtext("A", side=3, line=0.4,cex=2.2,adj=-0.2)
  par(fig=c(0.5,1,0.5,1),new=TRUE)
  plot(maps,sresii,pch=16,lwd=1,col="#757575",xlab=" ",
       ylab=" ",cex.lab=1.5,mgp = c(2.2, 1, 0),xaxt="n",yaxt="n",xlim=c(0.5,8.5))
  axis(1,seq(1,8,1),seq(1,8,1),lwd=2,cex.axis=1.2)
  axis(2,seq(round(min(sresii),0),round(max(sresii),0),2),
       seq(round(min(sresii),0),round(max(sresii),0),2),lwd=2,cex.axis=1.2)
  abline(h=0,lwd=2,lty=2)
  
  mtext(expression(italic("S. aureus")),side=3, line=0.4,cex=2.0,adj=0.5)
  
  par(fig=c(0,0.5,0,0.5),new=TRUE)
  par(mar=c(4,3.4,2.2,0.2),las=1)
  plot(mape,eresi,pch=16,lwd=1,col="#757575",xlab="Population",
       ylab="Residuals",cex.lab=1.5,mgp = c(2.2, 1, 0),xaxt="n",yaxt="n",xlim=c(0.5,5.5))
  axis(1,seq(1,5,1),seq(1,5,1),lwd=2,cex.axis=1.2)
  axis(2,seq(round(min(eresi),0),round(max(eresi),0),2),
       seq(round(min(eresi),0),round(max(eresi),0),2),lwd=2,cex.axis=1.2)
  abline(h=0,lwd=2,lty=2)
  
  mtext("B", side=3, line=0.4,cex=2.2,adj=-0.2)
  par(fig=c(0.5,1,0,0.5),new=TRUE)
  plot(maps,sresi,pch=16,lwd=1,col="#757575",xlab="Population",
       ylab=" ",cex.lab=1.5,mgp = c(2.2, 1, 0),xaxt="n",yaxt="n",xlim=c(0.5,8.5))
  axis(1,seq(1,8,1),seq(1,8,1),lwd=2,cex.axis=1.2)
  axis(2,seq(round(min(sresi),0),round(max(sresi),0),2),
       seq(round(min(sresi),0),round(max(sresi),0),2),lwd=2,cex.axis=1.2)
  abline(h=0,lwd=2,lty=2)
  dev.off()
  
}


ES.pvalue.qq.plot <- function(Pei,Pe,Psi,Ps){
  
  allpei <- c(Pei)
  allpe <- c(Pe)
  allpsi <- c(Psi)
  allps <- c(Ps)
  
  
  pdf("Figure-Pvalue-QQ.pdf",height=8,width=8)
  par(fig=c(0,0.5,0.5,1))
  par(mar=c(2,4,4.2,0.2),las=1)
  o=-log10(sort(allpei,decreasing=F))  
  e=-log10(ppoints(length(allpei))) 
  plot(e,o,pch=16,col="#757575",xlab=" ",
       ylab=expression(Observed ~ ~-log[10](italic(p))),cex.lab=1.5,mgp = c(2.5, 1, 0),xaxt="n",yaxt="n")
  axis(1,seq(0,5,1),seq(0,5,1),lwd=2,cex.axis=1.2)
  axis(2,seq(0,6,1),seq(0,6,1),lwd=2,cex.axis=1.2)
  abline(0,1,col= "red ",lwd=4)
  
  mtext(expression(italic("E. coli")),side=3, line=0.4,cex=2.0,adj=0.5)
  mtext("A", side=3, line=0.4,cex=2.2,adj=-0.2)
  
  par(fig=c(0.5,1,0.5,1),new=TRUE)
  o=-log10(sort(allpsi,decreasing=F))  
  e=-log10(ppoints(length(allpsi))) 
  plot(e,o,pch=16,col="#757575",xlab=" ",
       ylab=" ",cex.lab=1.5,mgp = c(2.2, 1, 0),xaxt="n",yaxt="n")
  axis(1,seq(0,5,1),seq(0,5,1),lwd=2,cex.axis=1.2)
  axis(2,seq(0,6,1),seq(0,6,1),lwd=2,cex.axis=1.2)
  abline(0,1,col= "red ",lwd=4)
  
  mtext(expression(italic("S. aureus")),side=3, line=0.4,cex=2.0,adj=0.5)
  
  par(fig=c(0,0.5,0,0.5),new=TRUE)
  o=-log10(sort(allpe,decreasing=F))  
  e=-log10(ppoints(length(allpe))) 
  par(mar=c(4,4,2.2,0.2),las=1)
  plot(e,o,pch=16,col="#757575",xlab=expression(Expected ~ ~-log[10](italic(p))),
       ylab=expression(Observed ~ ~-log[10](italic(p))),cex.lab=1.5,mgp = c(2.5, 1, 0),xaxt="n",yaxt="n")
  axis(1,seq(0,5,1),seq(0,5,1),lwd=2,cex.axis=1.2)
  axis(2,seq(0,5,1),seq(0,5,1),lwd=2,cex.axis=1.2)
  abline(0,1,col= "red ",lwd=4)
  mtext("B", side=3, line=0.4,cex=2.2,adj=-0.2)
  
  par(fig=c(0.5,1,0,0.5),new=TRUE)
  o=-log10(sort(allps,decreasing=F))  
  e=-log10(ppoints(length(allps))) 
  plot(e,o,pch=16,col="#757575",xlab=expression(Expected ~ ~-log[10](italic(p))),
       ylab=" ",cex.lab=1.5,mgp = c(2.5, 1, 0),xaxt="n",yaxt="n")
  axis(1,seq(0,5,1),seq(0,5,1),lwd=2,cex.axis=1.2)
  axis(2,seq(0,8,2),seq(0,8,2),lwd=2,cex.axis=1.2)
  abline(0,1,col= "red ",lwd=4)
  dev.off()
  
  
}


ES.single.man <- function (P,seg=9,nchr=12,dist=0.4,filename="Figure-E-man.pdf",
                           filename1="E-sig.txt",thre1=NULL) 
{
  logP <- P$logp
  if(is.null(thre1)){
    thre <- P$thre1
  }else{
    thre <- thre1
  }
  
  pos <- as.numeric(colnames(logP))
  n.pos <- pos/10^6
  Seg <- round(seq(0,round(max(n.pos),2),length=seg),1)
  chr <- c()
  col.i <- rep(c("#008B00","#EE9A00"),round(nchr/2))
  allcol <- c()
  for(i in 1:nchr){
    
    A <- length(which(n.pos<i*0.4))- length(which(n.pos<(i-1)*0.4))
    chr <- c(chr,rep(i,A))
    allcol <- c(allcol,rep(col.i[i],A))
  }
  sig.snp <- c()
  pdf(filename,height=12,width=6)
  for(i in 1:5){
    par(fig=c(0,1,0.6666,0.9999))
    par(mar=c(3.2,4,2.2,0.2),las=1)
    plot(0,0,pch=16,type="n",col="#757575",xlab="SNP Position (Mb)",xlim=c(0,max(n.pos)),ylim=c(0,max(c(logP[3*i-2,],thre*1.05))),
         ylab=expression(-log[10](italic(p))),cex.lab=1.5,mgp = c(2.3, 1, 0),xaxt="n",yaxt="n")
    axis(1,Seg,Seg,lwd=2,cex.axis=1.2)
    axis(2,seq(0,round(max(c(logP[3*i-2,],thre*1.05))),1),seq(0,round(max(c(logP[3*i-2,],thre*1.05))),1),lwd=2,cex.axis=1.2)
    abline(h=thre,col= "blue ",lwd=2,lty=2)
    points(n.pos,logP[3*i-2,],pch=16,col=allcol,cex=0.7)
    f1 <- paste("Time point ",3*i-1,sep="")
    mtext(f1,side=3, line=0.4,cex=2.0,adj=0.5)
    sig.snp <- c(sig.snp,pos[which(logP[3*i-2,]>thre1)])
    par(fig=c(0,1,0.3333,0.6666),new=TRUE)
    plot(0,0,pch=16,type="n",col="#757575",xlab="SNP Position (Mb)",xlim=c(0,max(n.pos)),ylim=c(0,max(c(logP[3*i-1,],thre*1.05))),
         ylab=expression(-log[10](italic(p))),cex.lab=1.5,mgp = c(2.3, 1, 0),xaxt="n",yaxt="n")
    axis(1,Seg,Seg,lwd=2,cex.axis=1.2)
    axis(2,seq(0,round(max(c(logP[3*i-1,],thre*1.05))),1),seq(0,round(max(c(logP[3*i-1,],thre*1.05))),1),lwd=2,cex.axis=1.2)
    abline(h=thre,col= "blue ",lwd=2,lty=2)
    points(n.pos,logP[3*i-1,],pch=16,col=allcol,cex=0.7)
    f2 <- paste("Time point ",3*i,sep="")
    mtext(f2,side=3, line=0.4,cex=2.0,adj=0.5)
    sig.snp <- c(sig.snp,pos[which(logP[3*i-1,]>thre1)])
    par(fig=c(0,1,0,0.3333),new=TRUE)
    plot(0,0,pch=16,type="n",col="#757575",xlab="SNP Position (Mb)",xlim=c(0,max(n.pos)),ylim=c(0,max(c(logP[3*i,],thre*1.05))),
         ylab=expression(-log[10](italic(p))),cex.lab=1.5,mgp = c(2.3, 1, 0),xaxt="n",yaxt="n")
    axis(1,Seg,Seg,lwd=2,cex.axis=1.2)
    axis(2,seq(0,round(max(c(logP[3*i,],thre*1.05))),1),seq(0,round(max(c(logP[3*i,],thre*1.05))),1),lwd=2,cex.axis=1.2)
    abline(h=thre,col= "blue ",lwd=2,lty=2)
    points(n.pos,logP[3*i,],pch=16,col=allcol,cex=0.7)
    f3 <- paste("Time point ",3*i+1,sep="")
    mtext(f3,side=3, line=0.4,cex=2.0,adj=0.5)
    sig.snp <- c(sig.snp,pos[which(logP[3*i,]>thre1)])
  }
  write.table(sig.snp,file=filename1)
  dev.off()
}

fp.plot_ind_curves<-function( dat, pheno1,pheno2,rows=NA, cols=NA, 
                              max_curves = NULL, selected=NULL,filename="1.pdf" )
{
  max_log <- 0;
  min_log <- 0;
  pheno <- cbind(pheno1,pheno2)
  if (is.null(max_curves))
    max_curves <- min(100,dat$sample_N);
  
  if ( is.null(selected) )
    selected <- c(1:max_curves)
  
  selected.greater <- which( selected > dat$sample_N); 
  if (length (selected.greater) )
    selected <- selected[- selected.greater ];
  
  if ( is.na (rows) )
  {
    cols <- ceiling(  (length(selected))^0.5 );
    rows <- ceiling( (length(selected))/cols);
  }
  
  minv <- min(as.matrix(pheno[selected,]), na.rm=TRUE)*0.9;
  maxv <- max(as.matrix(pheno[selected,]), na.rm=TRUE)*1.1;
  
  px<-dat$sample_times;
  xmax <- max(px, na.rm=TRUE);
  xmin <- min(px, na.rm=TRUE);
  
  py <- length ( pheno[,1] );
  rx <- ceiling( py/8 );
  xunit <- c(1, 2, 3, 5, 8, 10, 20, 50, 100, 200, 500, 1000, 2000);
  rx <- xunit[min(which(xunit>rx))];
  if (rx>2000) rx<-as.integer(py/8);
  
  ry<- ceiling((maxv - minv)/8)
  yunit <- c(1, 2, 3, 5, 8, 10, 20, 50, 100, 200, 500, 1000, 2000);
  ry <- yunit[min(which(yunit>ry))];
  if (ry>2000) ry<-as.integer((maxv - minv)/8);
  
  p.width <- (xmax-xmin)*cols/9*10;
  p.height <- (maxv-minv)*rows/9*10;
  pdf(filename,height=14,width=9)
  par(mar=c(1,1,1,1))
  plot(c(0,0), c(0,0), type="n",xaxt="n",yaxt="n",frame=T, xlab="",ylab="", xlim=c(0, p.width), ylim=c(0, p.height) );
  for (i in 0:(rows-1) )
  {
    for (j in 0:(cols-1) )
    {
      sub_fig<- c( 1/cols*j, 1/cols*(j+1), 1/rows*i,  1/rows*(i+1) )*0.98+0.02;
      sub_rc <- c( p.width * sub_fig[1], p.height* sub_fig[3], 
                   p.width * sub_fig[2], p.height* sub_fig[4] );
      rect(sub_rc[1], sub_rc[2], sub_rc[3], sub_rc[4], border="black", lwd=1);
      sub_fc <- c( sub_rc[1]-xmin, sub_rc[2]-minv );
      
      if ( (i*cols+j+1) <= length(selected) )
      {
        idx <- i*cols+j+1;
        px<-as.numeric( dat$sample_times);
        py<-as.numeric( pheno[selected[idx], ] );
        
        lines(px+sub_fc[1], py[1:15]+sub_fc[2], lwd=1,col="red");
        
        for (k in 1:length(px))
          points(px[k]+sub_fc[1], py[1:15][k]+sub_fc[2] ,type="o", pch=19, cex=0.5,col="red")
        
        lines(px+sub_fc[1], py[16:30]+sub_fc[2], lwd=1,col="blue");
        
        for (k in 1:length(px))
          points(px[k]+sub_fc[1], py[16:30][k]+sub_fc[2] ,type="o", pch=19, cex=0.5,col="blue")
        
        h<-( maxv - minv )
        if (i==0)
          for (k in as.integer(xmin/rx):as.integer(xmax/rx))
          {
            if ( k*rx+sub_fc[1]+xmin <= sub_rc[1] || k*rx+sub_fc[1]+xmin >= sub_rc[3]) next;
            segments(k*rx+sub_fc[1]+xmin, minv+sub_fc[2], k*rx+sub_fc[1]+xmin, minv+h/40+sub_fc[2]);
          }
        
        if (j==0)
          for (k in as.integer(minv/ry):as.integer(maxv/ry) )
          {
            if ( k*ry+sub_fc[2] <= sub_rc[2] || k*ry+sub_fc[2] >= sub_rc[4]) next;
            segments(xmin+sub_fc[1], k*ry+sub_fc[2], xmin+(px-1)/40+sub_fc[1], k*ry+sub_fc[2]);
          }
      }
      else
      {
        #plot(c(1:2), c(1:2),xlim=c(1, 10), ylim=c(minv, maxv), 
        #		type="n", main="", xlab="", ylab="", xaxt="n", yaxt="n",xaxs="i", yaxs="i");
      }
      
      if (i*cols+j+1 <= dat$sample_N)
        text( 1+(xmax-1)/5 + sub_fc[1], h*8.8/10+minv + sub_fc[2], paste(i*cols+j+1, "", sep=""), font=4);
      
      if (j==0)
      {
        for (k in as.integer(minv/ry):(as.integer(maxv/ry)) )
        {
          if ((k*ry)<maxv && (k*ry)>minv && (k*ry-minv)/(maxv-minv)<0.95 && (k*ry-minv)/(maxv-minv)>0.01)
            text(0.5+sub_fc[1], k*ry+ sub_fc[2], k*ry, adj=c(1, 0.5),cex=1);
        }
      }
      
      if (i==0)
      {
        for (k in as.integer(xmin/rx):as.integer(xmax/rx))
        {	
          if ((k*rx)/xmax<0.9 && (k*rx)/xmax>0.05)
            text(k*rx+ sub_fc[1]+xmin, -1+ sub_fc[2]+minv, (k*rx), adj=c(0.5,1),cex=1)
        }
      }
    }
  }
  dev.off()
}


fp.plot_ind_curves.com<-function( dat, pheno1,pheno2,rows=NA, cols=NA, 
                                  max_curves = NULL, selected=NULL,filename="1.pdf" )
{
  max_log <- 0;
  min_log <- 0;
  
  if (is.null(max_curves))
    max_curves <- min(100,dat$sample_N);
  
  if ( is.null(selected) )
    selected <- c(1:max_curves)
  
  selected.greater <- which( selected > dat$sample_N); 
  if (length (selected.greater) )
    selected <- selected[- selected.greater ];
  
  if ( is.na (rows) )
  {
    cols <- ceiling(  (length(selected))^0.5 );
    rows <- ceiling( (length(selected))/cols);
  }
  
  pheno <- cbind(pheno1,pheno2)
  minv <- min(as.matrix(pheno[selected,]), na.rm=TRUE)*0.9;
  maxv <- max(as.matrix(pheno[selected,]), na.rm=TRUE)*1.1;
  NN <- dat$sample_N
  fpar <- c()
  for(i in 1:NN){
    
    e <- as.numeric((pheno1[i,]))
    s <- as.numeric((pheno2[i,]))
    value <- c()
    allpar <- list()
    for(ii in 1:10){
      s.par <- c(0.1585815,4.5978148,125.2705139,0.1233422,2.8023220,91.2431834)*runif(6,0.85,1.15)
      res <- try(optim(s.par,s.mle,s.y=c(e,s),s.t=dat1$sample_times,x0=e[1],y0=s[1],
                       method="BFGS",control=list(trace=T,maxit=10000)),TRUE)
      if(class(res)=="try-error")
        next;
      value <- c(value,res$value)
      allpar[[ii]] <- res$par
    }
    parin <- allpar[[which(value==min(value))]]
    fpar <- rbind(fpar,parin)
  }
  
  pheno11 <- c()
  for(i in 1:NN){
    e <- as.numeric((pheno1[i,]))
    s <- as.numeric((pheno2[i,]))
    r1 <- com.get_mu(fpar[i,],times=dat1$sample_times,x0=e[1],y0=s[1]) 
    pheno11 <- rbind(pheno11,r1)
  }
  px<-dat$sample_times;
  xmax <- max(px, na.rm=TRUE);
  xmin <- min(px, na.rm=TRUE);
  
  py <- length ( pheno[,1] );
  rx <- ceiling( py/8 );
  xunit <- c(1, 2, 3, 5, 8, 10, 20, 50, 100, 200, 500, 1000, 2000);
  rx <- xunit[min(which(xunit>rx))];
  if (rx>2000) rx<-as.integer(py/8);
  
  ry<- ceiling((maxv - minv)/8)
  yunit <- c(1, 2, 3, 5, 8, 10, 20, 50, 100, 200, 500, 1000, 2000);
  ry <- yunit[min(which(yunit>ry))];
  if (ry>2000) ry<-as.integer((maxv - minv)/8);
  
  p.width <- (xmax-xmin)*cols/9*10;
  p.height <- (maxv-minv)*rows/9*10;
  pdf(filename,height=14,width=9)
  par(mar=c(1,1,1,1))
  plot(c(0,0), c(0,0), type="n",xaxt="n",yaxt="n",frame=T, xlab="",ylab="", xlim=c(0, p.width), ylim=c(0, p.height) );
  for (i in 0:(rows-1) )
  {
    for (j in 0:(cols-1) )
    {
      sub_fig<- c( 1/cols*j, 1/cols*(j+1), 1/rows*i,  1/rows*(i+1) )*0.98+0.02;
      sub_rc <- c( p.width * sub_fig[1], p.height* sub_fig[3], 
                   p.width * sub_fig[2], p.height* sub_fig[4] );
      rect(sub_rc[1], sub_rc[2], sub_rc[3], sub_rc[4], border="black", lwd=1);
      sub_fc <- c( sub_rc[1]-xmin, sub_rc[2]-minv );
      
      if ( (i*cols+j+1) <= length(selected) )
      {
        idx <- i*cols+j+1;
        px<-as.numeric( dat$sample_times);
        py<-as.numeric( pheno[selected[idx], ] );
        py1<-as.numeric( pheno11[selected[idx], ] );
        if (length(which(py==-1))>0)
        {
          px<-px[-(which(py==-1))];
          py<-py[-(which(py==-1))];
          py1<-py1[-(which(py1==-1))];
        }
        
        if (length(which(is.na(py)))>0)
        {
          px<-px[-(which(is.na(py)))];
          py<-py[-(which(is.na(py)))];
          py1<-py1[-(which(is.na(py1)))];
        }
        
        
        for (k in 1:length(px))
          points(px[k]+sub_fc[1], py[1:15][k]+sub_fc[2] ,type="o", pch=19, cex=0.5,col="red")
        
        lines(px+sub_fc[1], py1[1:15]+sub_fc[2], lwd=2,col="red");
        for (k in 1:length(px))
          points(px[k]+sub_fc[1], py[16:30][k]+sub_fc[2] ,type="o", pch=19, cex=0.5,col="blue")
        
        lines(px+sub_fc[1], py1[16:30]+sub_fc[2], lwd=2,col="blue");
        h<-( maxv - minv )
        if (i==0)
          for (k in as.integer(xmin/rx):as.integer(xmax/rx))
          {
            if ( k*rx+sub_fc[1]+xmin <= sub_rc[1] || k*rx+sub_fc[1]+xmin >= sub_rc[3]) next;
            segments(k*rx+sub_fc[1]+xmin, minv+sub_fc[2], k*rx+sub_fc[1]+xmin, minv+h/40+sub_fc[2]);
          }
        
        if (j==0)
          for (k in as.integer(minv/ry):as.integer(maxv/ry) )
          {
            if ( k*ry+sub_fc[2] <= sub_rc[2] || k*ry+sub_fc[2] >= sub_rc[4]) next;
            segments(xmin+sub_fc[1], k*ry+sub_fc[2], xmin+(px-1)/40+sub_fc[1], k*ry+sub_fc[2]);
          }
      }
      else
      {
        #plot(c(1:2), c(1:2),xlim=c(1, 10), ylim=c(minv, maxv), 
        #		type="n", main="", xlab="", ylab="", xaxt="n", yaxt="n",xaxs="i", yaxs="i");
      }
      
      if (i*cols+j+1 <= dat$sample_N)
        text( 1+(xmax-1)/5 + sub_fc[1], h*8.8/10+minv + sub_fc[2], paste(i*cols+j+1, "", sep=""), font=4);
      
      if (j==0)
      {
        for (k in as.integer(minv/ry):(as.integer(maxv/ry)) )
        {
          if ((k*ry)<maxv && (k*ry)>minv && (k*ry-minv)/(maxv-minv)<0.95 && (k*ry-minv)/(maxv-minv)>0.01)
            text(0.5+sub_fc[1], k*ry+ sub_fc[2], k*ry, adj=c(1, 0.5),cex=1);
        }
      }
      
      if (i==0)
      {
        for (k in as.integer(xmin/rx):as.integer(xmax/rx))
        {	
          if ((k*rx)/xmax<0.9 && (k*rx)/xmax>0.05)
            text(k*rx+ sub_fc[1]+xmin, -1+ sub_fc[2]+minv, (k*rx), adj=c(0.5,1),cex=1)
        }
      }
    }
  }
  dev.off()
}
ES.plot.mean <- function(dat,allpar){
  
  times <- dat$sample_times
  nt <- seq(min(times),max(times),1)
  N <- dat$sample_N
  x0 <- mean(dat$ep.p[,1])
  y0 <- mean(dat$sp.p[,1])
  ny <- com.get_mu(allpar,times=nt,x0=x0,y0=y0)
  pdf("Figure-all-mean.pdf",height=4,width=8)
  par(fig=c(0,0.5,0,1))
  par(mar=c(2,3.4,4.2,0.2),las=1)
  plot(times,dat$ep.p[1,],type="l",lwd=2,col="#757575",xlab=" ",
       ylab="Abundance",cex.lab=1.5,mgp = c(2.2, 1, 0),xaxt="n",yaxt="n",
       xlim=c(times[1],times[15]),ylim=c(min(dat$ep.p)-3,max(dat$ep.p)))
  axis(2,c(10,15,20,25,30),c(10,15,20,25,30),lwd=2,cex.axis=1.2)
  axis(1, seq(0,36,6), seq(0,36,6),lwd=2,cex.axis=1.2)
  for(i in 2:N){
    lines(times,dat$ep.p[i,],lwd=2,col="#757575")
  }
  lines(nt,ny[1:length(nt)],lwd=3,col="red")
  mtext(expression(italic("E. coli")),side=3, line=0.4,cex=2.0,adj=0.5)
  #mtext("A", side=3, line=0.4,cex=2.2,adj=-0.2)
  par(fig=c(0.5,1,0,1),new=TRUE)
  plot(times,dat$sp.p[1,],type="l",lwd=2,col="#757575",xlab=" ",
       ylab=" ",cex.lab=1.5,mgp = c(2.2, 1, 0),xaxt="n",yaxt="n",
       xlim=c(times[1],times[15]),ylim=c(min(dat$sp.p)-3,max(dat$sp.p)))
  axis(2,c(10,15,20,25,30),c(10,15,20,25,30),lwd=2,cex.axis=1.2)
  axis(1, seq(0,36,6), seq(0,36,6),lwd=2,cex.axis=1.2)
  for(i in 2:N){
    lines(times,dat$sp.p[i,],lwd=2,col="#757575")
  }
  lines(nt,ny[(length(nt)+1):(2*length(nt))],lwd=3,col="red")
  mtext(expression(italic("S. aureus")),side=3, line=0.4,cex=2.0,adj=0.5)
  dev.off()
}


manh1 <- function(snpinfo,snp,LR,file="manh-1.pdf"){
  
  df <- c()
  for(i in 1:dim(snp)[2]){
    
    ll <- length(table(snp[,i]))
    df <- c(df,ll)
  }
  
  z <- -log(pchisq(LR,df*6,lower.tail = F))
  
  x <- as.numeric(snpinfo[1,])
  y <- as.numeric(snpinfo[2,])
  
  colM <- rep("black",length(x))
  colM[which(z>-log(1e-6/length(x)))] <- "red"
  
  resdata <- data.frame(x=x,y=y,z=z,c=colM)
  colnames(resdata) <- c("x","y","-log(P-value)","cc")
  
  sig.snp <- snp[,which(z>-log(1e-6/length(x)))]
  sig.snp.info <- snpinfo[,which(z>-log(1e-6/length(x)))]
  
  write.csv(sig.snp.info,file="sig-snp.csv")
  p <- ggplot(resdata, aes(x = x, y = y, size = `-log(P-value)`,colour=as.character(cc)))
  p <- p + geom_point()
  p <- p + scale_size_continuous(range = c(0,5))
  p <- p + scale_x_continuous(breaks=c(0,1e6,2e6,3e6,4e6),labels=c(0:4),expand = c(0.02,0.02))
  p <- p + scale_y_continuous(breaks=c(0,5e5,1e6,1.5e6,2e6),labels=c(0,0.5,1,1.5,2),expand = c(0.02,0.02))
  p <- p + scale_color_manual(" ",breaks=(c("Metabolism","Organismal Systems")),values=c("black","red"))
  p <- p + xlab(expression(paste(italic("E. coli"),"(Mb)"))) 
  p <- p + ylab(expression(paste(italic("S. aureus"),"(Mb)")))
  p <- p + theme_zg_m()
  ggsave(file,height=6,width=10)
}






theme_zg_m <- function(..., bg='white'){
  require(grid)
  theme_classic(...,base_family="serif") +
    theme(rect=element_rect(fill=bg),
          plot.margin=unit(c(1,0.5,0.5,0.5), 'lines'),
          panel.background=element_rect(fill='transparent', color='black'),
          panel.border=element_rect(fill='transparent', color='transparent'),
          panel.grid=element_blank(),
          plot.title=element_text(size=25,vjust=1,hjust=0),
          axis.title = element_text(color='black',size=25),
          axis.title.x = element_text(vjust=0),
          axis.title.y = element_text(vjust=1.5),
          axis.ticks.length = unit(0.4,"lines"),
          axis.text.x=element_text(size=20),
          axis.text.y=element_text(size=20),
          axis.ticks = element_line(color='black'),
          axis.ticks.margin = unit(0.8,"lines"),
          #legend.direction = "horizontal", 
          #legend.justification = "top", 
          legend.position="top",
          legend.text=element_text(color='black',size=15,vjust=3),
          legend.title=element_text(color='black',size=25),
          strip.background=element_rect(fill='transparent', color='transparent'),
          strip.text=element_text(color='black',vjust=-3,hjust=0.06,size=16))
  
}

genotype.plot <- function(allpar,snp,index,times=dat$sample_times,epp=dat$ep.p,spp=dat$sp.p){
  
  snp.type <- names(table(as.character(unlist(c(snp[,index])))))
  ty <- cbind(epp,spp)
  glist <- list()
  respar <- allpar[index,8:31]
  k <- 1
  for(i in 1:length(snp.type)){
    index1 <- which(snp[,index]==snp.type[i])
    cm <- as.numeric(colMeans(ty[index1,]))
    fitc <- com.get_mu(respar[(i*6-5):(6*i)],times,x0=cm[1],y0=cm[16])
    trait1 <- fitc[1:15]
    trait2 <- fitc[16:30]
    phe <- c();phs <- c();
    for(ii in 1:length(index1)){
      tmp1 <- cbind(rep(ii,length(times)),times,epp[index1[ii],])
      tmp2 <- cbind(rep(ii,length(times)),times,spp[index1[ii],])
      phe <- rbind(phe,tmp1)
      phs <- rbind(phs,tmp2)
    }
    phef <- data.frame(time=times,ph=trait1)
    phsf <- data.frame(time=times,ph=trait2)
    phe <- as.data.frame(phe);phs <- as.data.frame(phs)
    colnames(phe) <- c("ind","time","ph")
    colnames(phs) <- c("ind","time","ph")
    g1 <- ggplot(phe)
    g1 <- g1 + geom_line(aes(x=time,y=ph,group=ind),colour="#B0B0B0",aphla=0.5)
    g1 <- g1 + geom_line(data=phef,aes(x=time,y=ph),colour="red",size=1)
    g1 <- g1 + scale_x_continuous(breaks=c(seq(6,36,6)),labels=seq(6,36,6))
    g1 <- g1 + scale_y_continuous(limits=c(min(dat$ep.p),max(dat$ep.p)))
    g1 <- g1 + annotate("text",x=6,y=max(dat$ep.p)*0.97,label=snp.type[i],size=8)
    g1 <- g1 + xlab("Time (h)")+ ylab(expression(paste(italic("E. coli"),"(Abundance)")))
    g1 <- g1 + theme_zg_m()
    g2 <- ggplot(phs)
    g2 <- g2 + geom_line(aes(x=time,y=ph,group=ind),colour="#B0B0B0",aphla=0.5)
    g2 <- g2 + geom_line(data=phsf,aes(x=time,y=ph),colour="red",size=1)
    g2 <- g2 + scale_x_continuous(breaks=c(seq(6,36,6)),labels=seq(6,36,6))
    g2 <- g2 + scale_y_continuous(limits=c(min(dat$sp.p),max(dat$sp.p)))
    g2 <- g2 + annotate("text",x=6,y=max(dat$sp.p)*0.97,label=snp.type[i],size=8)
    g2 <- g2 + xlab("Time (h)")+ ylab(expression(paste(italic("S. aureus"),"(Abundance)")))
    g2 <- g2 + theme_zg_m()
    
    glist[[k]] <- g1
    glist[[k+1]] <- g2
    k <- k + 2
  }
  if(length(glist)==8){
    pdf(paste("SNP",index,".pdf",sep=""),heigh=9,width=20)
    multiplot(plotlist = glist,cols=4)
    dev.off()
  }
}


multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

ES.plot.interaction <- function(allpar,dat,snp,index){
  
  snp.type <- names(table(as.character(unlist(c(snp[,index])))))
  
  respar <- allpar[index,8:31]
  parin <- c()
  for(i in 1:length(snp.type)){
    tmppar <- respar[(i*6-5):(6*i)]
    tmppar1 <- c(tmppar[3],tmppar[6]/tmppar[5],tmppar[6],tmppar[3]/tmppar[2])
    parin <- rbind(parin,tmppar1)
  }  
  maxx <- max(parin[,1:2])+60;
  maxy <- max(parin[,3:4])+30;
  file <- paste("Figure-interaction",index,".pdf",sep="")
  pdf(file,height=7.5,width=9)
  par(fig=c(0,0.5,0.5,1))
  par(mar=c(0,0,0,0),las=1, bty="n" )
  plot( 1:2, 1:2, xlim=c( -45, maxx), ylim=c(-20,maxy), 
        type="n", main="", xlab="", ylab="", xaxt="n", yaxt="n",xaxs="i", yaxs="i")
  segments(parin[1,1],0,0,parin[1,4],col="#00FF00",lwd=3)
  segments(parin[1,2],0,0,parin[1,3],col="#CD2626",lwd=3)
  arrows(0,0,0,maxy-15,col="#104E8B",lwd=3,angle=30,length=0.14)
  arrows(0,0,maxx-45,0,col="#104E8B",lwd=3,angle=30,length=0.14)
  text(-5,-5,0,col="black",cex=1.4)
  text(parin[1,1],-7,expression(italic(K[1])),col="black",cex=1.4)
  text(-25,parin[1,4],expression(italic(K[1]/alpha[12])),col="black",cex=1.4)
  text(parin[1,2],-7,expression(italic(K[2]/alpha[21])),col="black",cex=1.4)
  text(-15,parin[1,3],expression(italic(K[2])),col="black",cex=1.4)
  text(maxx-22,0,expression(italic("E. coli")),col="black",cex=1.4)
  text(0,maxy-9,expression(italic("S. aureus")),col="black",cex=1.4)
  text((maxx-45)/2,maxy-15,snp.type[1],col="black",cex=1.4)
  inter.type <- "Competition"
  text((maxx-45)/2,maxy-25,inter.type,col="red",cex=1.4)
  text((maxx-180),maxy-25,"A",col="black",cex=1.4)
  par(fig=c(0.5,1,0.5,1),new=TRUE)
  plot( 1:2, 1:2, xlim=c( -45, maxx), ylim=c(-20,maxy), 
        type="n", main="", xlab="", ylab="", xaxt="n", yaxt="n",xaxs="i", yaxs="i")
  segments(parin[2,1],0,0,parin[2,4],col="#00FF00",lwd=3)
  segments(parin[2,2],0,0,parin[2,3],col="#CD2626",lwd=3)
  arrows(0,0,0,maxy-15,col="#104E8B",lwd=3,angle=30,length=0.14)
  arrows(0,0,maxx-45,0,col="#104E8B",lwd=3,angle=30,length=0.14)
  text(-5,-5,0,col="black",cex=1.4)
  text(parin[2,1],-7,expression(italic(K[1])),col="black",cex=1.4)
  text(-25,parin[2,4],expression(italic(K[1]/alpha[12])),col="black",cex=1.4)
  text(parin[2,2],-7,expression(italic(K[2]/alpha[21])),col="black",cex=1.4)
  text(-15,parin[2,3],expression(italic(K[2])),col="black",cex=1.4)
  text(maxx-22,0,expression(italic("E. coli")),col="black",cex=1.4)
  text(0,maxy-9,expression(italic("S. aureus")),col="black",cex=1.4)
  text((maxx-45)/2,maxy-15,snp.type[2],col="black",cex=1.4)
  inter.type <- "Competition"
  text((maxx-45)/2,maxy-25,inter.type,col="red",cex=1.4)
  text((maxx-180),maxy-25,"B",col="black",cex=1.4)
  par(fig=c(0,0.5,0,0.5),new=TRUE)
  plot( 1:2, 1:2, xlim=c( -45, maxx), ylim=c(-20,maxy), 
        type="n", main="", xlab="", ylab="", xaxt="n", yaxt="n",xaxs="i", yaxs="i")
  segments(parin[3,1],0,0,parin[3,4],col="#00FF00",lwd=3)
  segments(parin[3,2],0,0,parin[3,3],col="#CD2626",lwd=3)
  arrows(0,0,0,maxy-15,col="#104E8B",lwd=3,angle=30,length=0.14)
  arrows(0,0,maxx-45,0,col="#104E8B",lwd=3,angle=30,length=0.14)
  text(-5,-5,0,col="black",cex=1.4)
  text(parin[3,1],-7,expression(italic(K[1])),col="black",cex=1.4)
  text(-25,parin[3,4],expression(italic(K[1]/alpha[12])),col="black",cex=1.4)
  text(parin[3,2],-7,expression(italic(K[2]/alpha[21])),col="black",cex=1.4)
  text(-15,parin[3,3],expression(italic(K[2])),col="black",cex=1.4)
  text(maxx-22,0,expression(italic("E. coli")),col="black",cex=1.4)
  text(0,maxy-9,expression(italic("S. aureus")),col="black",cex=1.4)
  text((maxx-45)/2,maxy-15,snp.type[3],col="black",cex=1.4)
  inter.type <- "Competition"
  text((maxx-45)/2,maxy-25,inter.type,col="red",cex=1.4)
  text((maxx-180),maxy-25,"C",col="black",cex=1.4)
  par(fig=c(0.5,1,0,0.5),new=TRUE)
  plot( 1:2, 1:2, xlim=c( -45, maxx), ylim=c(-20,maxy), 
        type="n", main="", xlab="", ylab="", xaxt="n", yaxt="n",xaxs="i", yaxs="i")
  segments(parin[4,1],0,0,parin[4,4],col="#00FF00",lwd=3)
  segments(parin[4,2],0,0,parin[4,3],col="#CD2626",lwd=3)
  arrows(0,0,0,maxy-15,col="#104E8B",lwd=3,angle=30,length=0.14)
  arrows(0,0,maxx-45,0,col="#104E8B",lwd=3,angle=30,length=0.14)
  text(-5,-5,0,col="black",cex=1.4)
  text(parin[4,1],-7,expression(italic(K[1])),col="black",cex=1.4)
  text(-25,parin[4,4],expression(italic(K[1]/alpha[12])),col="black",cex=1.4)
  text(parin[4,2],-7,expression(italic(K[2]/alpha[21])),col="black",cex=1.4)
  text(-15,parin[4,3],expression(italic(K[2])),col="black",cex=1.4)
  text(maxx-22,0,expression(italic("E. coli")),col="black",cex=1.4)
  text(0,maxy-9,expression(italic("S. aureus")),col="black",cex=1.4)
  text((maxx-45)/2,maxy-15,snp.type[4],col="black",cex=1.4)
  inter.type <- "Competition"
  text((maxx-45)/2,maxy-25,inter.type,col="red",cex=1.4)
  text((maxx-180),maxy-25,"D",col="black",cex=1.4)
  dev.off()
}
one.mu.plot <- function(allpar,dat,snp,index,epp=dat$ep.p,spp=dat$sp.p){
  
  require(ggplot2)
  snp.type <- names(table(as.character(unlist(c(snp[,index])))))
  respar <- allpar[index,8:31]
  ty <- cbind(epp,spp)
  times <- seq(dat$sample_times[1],dat$sample_times[15],0.1)
  trait1 <- c();trait2 <- c();
  for(i in 1:length(snp.type)){
    index1 <- which(snp[,index]==snp.type[i])
    cm <- as.numeric(colMeans(ty[index1,]))
    fitc <- com.get_mu(respar[(i*6-5):(6*i)],times,x0=cm[1],y0=cm[16])
    tmp1 <- cbind(rep(i,length(times)),times,fitc[1:length(times)])
    tmp2 <- cbind(rep(i,length(times)),times,fitc[(1+length(times)):(2*length(times))])
    trait1 <- rbind(trait1,tmp1)
    trait2 <- rbind(trait2,tmp2)
  }
  trait1 <- as.data.frame(trait1)
  trait2 <- as.data.frame(trait2)
  colnames(trait1) <- c("ind","times","effect")
  colnames(trait2) <- c("ind","times","effect")
  coll <- c("black","blue","red", "darkgreen")
  p <- ggplot(trait1)
  p <- p + geom_line(aes(x=times,y=effect,colour=as.factor(ind)),linetype=1,size=1)
  p <- p + scale_x_continuous(limits=c(2,36),breaks=seq(6,36,6),labels=seq(6,36,6))
  p <- p + scale_y_continuous(limits=c(14,26),breaks=seq(15,25,5),labels=seq(15,25,5))
  p <- p + scale_colour_manual(values = c("black","blue","red", "darkgreen"),breaks=c(1:4),labels=snp.type)
  p <- p + annotate("text", x = 4.5, y = 26*0.99, label=c("A"),size=8)
  p <- p + xlab("Time(h)") + ylab(expression(paste(italic("E. coli"),"(Abundance)")))+ theme_zg_one_mu(base_size = 20)
  #p
  
  p1 <- ggplot(trait2)
  p1 <- p1 + geom_line(aes(x=times,y=effect,colour=as.factor(ind)),linetype=1,size=1)
  p1 <- p1 + scale_x_continuous(limits=c(2,36),breaks=seq(6,36,6),labels=seq(6,36,6))
  p1 <- p1 + scale_y_continuous(limits=c(12,24),breaks=seq(13,23,4),labels=seq(13,23,4))
  p1 <- p1 + scale_colour_manual(values = c("black","blue","red", "darkgreen"),breaks=c(1:4),labels=snp.type)
  p1 <- p1 + annotate("text", x = 4.5, y = 24*0.99, label=c("B"),size=8)
  p1 <- p1 + xlab("Time(h)") + ylab(expression(paste(italic("S. aureus"),"(Abundance)")))+ theme_zg_one_mu(base_size = 20)
  #p1
  
  trait3 <- data.frame(ind=trait1[,1],x=trait1[,3],y=trait2[,3])
  p2 <- ggplot(trait3)
  p2 <- p2 + geom_point(aes(x=x,y=y,colour=as.factor(ind)), linetype=1,size=1)
  p2 <- p2 + scale_colour_manual(values = c("black","blue","red", "darkgreen"),breaks=c(1:4),labels=snp.type)
  p2 <- p2 + scale_x_continuous(limits=c(14,26),breaks=seq(15,25,5),labels=seq(15,25,5))
  p2 <- p2 + scale_y_continuous(limits=c(12,24),breaks=seq(13,23,4),labels=seq(13,23,4))
  p2 <- p2 + annotate("text", x =15, y = 24*0.99, label=c("C"),size=8)
  p2 <- p2 + xlab(expression(paste(italic("S. aureus"),"(Abundance)"))) + ylab(expression(paste(italic("E. coli"),"(Abundance)")))+ theme_zg_one_mu(base_size = 20)
  file <- paste("Figure-genotype",index,".pdf",sep="")
  pdf(file,height=11,width=5.5)
  multiplot(p,p1,p2,cols=1)
  dev.off()
  #ggsave("a.pdf",height=3.5,width=5)
  
}



theme_zg_one_mu <- function(..., bg='white'){
  require(grid)
  theme_classic(...,base_family="serif") +
    theme(rect=element_rect(fill=bg),
          plot.margin=unit(rep(0.5,4), 'lines'),
          panel.background=element_rect(fill='transparent', color='black'),
          panel.border=element_rect(fill='transparent', color='transparent'),
          panel.grid=element_blank(),
          plot.title=element_text(size=12,hjust=0.5),
          axis.title = element_text(color='black',hjust=0.5,vjust=0.3,size=20),
          axis.ticks.length = unit(-0.4,"lines"),
          axis.text.x=element_text(),
          axis.text.y=element_text(angle = 90,hjust=0.5),
          axis.title.x=element_text(vjust=-0.2),
          axis.title.y=element_text(vjust=1.4),
          axis.ticks = element_line(color='black'),
          axis.ticks.margin = unit(0.8,"lines"),
          legend.position='right',
          legend.title=element_blank(),
          #legend.title = element_text(size = rel(0.8), face = "bold", hjust = 0),
          legend.key=element_rect(fill='transparent', color='transparent'),
          strip.background=element_rect(fill='transparent', color='transparent'),
          strip.text=element_text(color='black',vjust=-1.8,hjust=0.02,size=20))
  
}

one.abi.plot <- function(allpar,dat,snp,index,epp=dat$ep.p,spp=dat$sp.p){
  
  require(ggplot2)
  snp.type <- names(table(as.character(unlist(c(snp[,index])))))
  respar <- allpar[index,8:31]
  ty <- cbind(epp,spp)
  times <- seq(dat$sample_times[1],dat$sample_times[15],0.1)
  fitres <- c()
  for(i in 1:length(snp.type)){
    index1 <- which(snp[,index]==snp.type[i])
    cm <- as.numeric(colMeans(ty[index1,]))
    fitc <- com.get_mu(respar[(i*6-5):(6*i)],times,x0=cm[1],y0=cm[16])
    fitres <- rbind(fitres,fitc)
  }
  effect <- one.abi.effect(mu=fitres)
  a11 <- data.frame(x=times,y=effect[1,1:341])
  a12 <- data.frame(x=times,y=effect[2,1:341])
  i1 <- data.frame(x=times,y=effect[3,1:341])
  a21 <- data.frame(x=times,y=effect[1,342:682])
  a22 <- data.frame(x=times,y=effect[2,342:682])
  i2 <- data.frame(x=times,y=effect[3,342:682])
  
  p <- ggplot(a11)
  p <- p + geom_line(aes(x=x,y=y),linetype=1,size=1,colour="red")
  p <- p + geom_line(data=a12,aes(x=x,y=y),linetype=1,size=1,colour="darkgreen")
  p <- p + geom_line(data=i1,aes(x=x,y=y),linetype=1,size=1,colour="blue")
  p <- p + scale_x_continuous(limits=c(2,40),breaks=seq(6,36,6),labels=seq(6,36,6))
  p <- p + annotate("text", x = 39, y = (a11[341,2]),label = paste("a","[","1","<","-","1","]"),parse = TRUE,size=6,colour="red") 
  p <- p + annotate('text', x = 39, y = (a12[341,2]),label = paste("a","[","1","<","-","2","]"),parse = TRUE,,size=6,colour="darkgreen") 
  p <- p + annotate('text', x = 39, y = (i1[341,2]),label = "I[1]",parse = TRUE,size=6,colour="blue") 
  p <- p + annotate("text", x = 4.5, y =max(effect[,1:341])*0.95, label=c("A"),size=8)
  p <- p + xlab("Time(h)") + ylab(expression(paste(italic("E. coli"),"(Abundance)")))+ theme_zg_one_abi(base_size = 20)
  #p
  
  p1 <- ggplot(a21)
  p1 <- p1 + geom_line(aes(x=x,y=y),linetype=1,size=1,colour="red")
  p1 <- p1 + geom_line(data=a22,aes(x=x,y=y),linetype=1,size=1,colour="darkgreen")
  p1 <- p1 + geom_line(data=i2,aes(x=x,y=y),linetype=1,size=1,colour="blue")
  p1 <- p1 + scale_x_continuous(limits=c(2,40),breaks=seq(6,36,6),labels=seq(6,36,6))
  p1 <- p1 + annotate("text", x = 39, y = (a21[341,2]),label = paste("a","[","2","<","-","1","]"),parse = TRUE,size=6,colour="red") 
  p1 <- p1 + annotate('text', x = 39, y = (a22[341,2]),label = paste("a","[","2","<","-","2","]"),parse = TRUE,,size=6,colour="darkgreen") 
  p1 <- p1 + annotate('text', x = 39, y = (i2[341,2]),label = "I[2]",parse = TRUE,size=6,colour="blue") 
  p1 <- p1 + annotate("text", x = 4.5, y =max(effect[,342:682])*0.95, label=c("B"),size=8)
  p1 <- p1 + xlab("Time(h)") + ylab(expression(paste(italic("S. aureus"),"(Abundance)")))+ theme_zg_one_abi(base_size = 20)
  #p1
  
  a <- data.frame(ind=rep(1,341),x=effect[1,1:341],y=effect[1,342:682])
  b <- data.frame(ind=rep(1,341),x=effect[2,1:341],y=effect[2,342:682])
  i <- data.frame(ind=rep(1,341),x=effect[3,1:341],y=effect[3,342:682])
  p2 <- ggplot(a)
  p2 <- p2 + geom_point(aes(x=x,y=y),colour="red", linetype=1,size=1)
  p2 <- p2 + geom_point(data=b,aes(x=x,y=y),colour="darkgreen", linetype=1,size=1)
  p2 <- p2 + geom_point(data=i,aes(x=x,y=y),colour="blue", linetype=1,size=1)
  p2 <- p2 + annotate("text", x =-0.93, y = max(effect[,342:682])*0.95, label=c("C"),size=8)
  p2 <- p2 + annotate('text', x = a$x[341]*0.8, y = a$y[341]*1.3,label = "a[1]",parse = TRUE,size=6,colour="red") 
  p2 <- p2 + annotate('text', x = b$x[341]*0.8, y = b$y[341]*1.3,label = "a[2]",parse = TRUE,size=6,colour="darkgreen") 
  p2 <- p2 + annotate('text', x = i$x[341]*0.8, y =  i$y[341]*1.3,label = "I",parse = TRUE,size=6,colour="blue")
  p2 <- p2 + xlab(expression(paste(italic("S. aureus"),"(Abundance)"))) + ylab(expression(paste(italic("E. coli"),"(Abundance)")))+ theme_zg_one_abi(base_size = 20)
  #p2
  file <- paste("Figure-genetic-effect",index,".pdf",sep="")
  pdf(file,height=11,width=5.5)
  multiplot(p,p1,p2,cols=1)
  dev.off()
  #ggsave("a.pdf",height=3.5,width=5)
  
}




one.abi.effect <- function(mu){
  
  a1 <- (mu[4,]+mu[3,]-mu[2,]-mu[1,])/4
  a2 <- (mu[4,]+mu[2,]-mu[3,]-mu[1,])/4
  i  <- (mu[4,]+mu[1,]-mu[3,]-mu[2,])/4
  effect <- rbind(a1,a2,i)
  return(effect)
}


theme_zg_one_abi <- function(..., bg='white'){
  require(grid)
  theme_classic(...,base_family="serif") +
    theme(rect=element_rect(fill=bg),
          plot.margin=unit(rep(0.5,4), 'lines'),
          panel.background=element_rect(fill='transparent', color='black'),
          panel.border=element_rect(fill='transparent', color='transparent'),
          panel.grid=element_blank(),
          plot.title=element_text(size=12,hjust=0.5),
          axis.title = element_text(color='black',hjust=0.5,vjust=0.3,size=20),
          axis.ticks.length = unit(-0.4,"lines"),
          axis.text.x=element_text(),
          axis.text.y=element_text(angle = 90,hjust=0.5),
          axis.title.x=element_text(vjust=-0.2),
          axis.title.y=element_text(vjust=1.4),
          axis.ticks = element_line(color='black'),
          axis.ticks.margin = unit(0.8,"lines"),
          legend.position='none',
          legend.title=element_blank(),
          #legend.title = element_text(size = rel(0.8), face = "bold", hjust = 0),
          legend.key=element_rect(fill='transparent', color='transparent'),
          strip.background=element_rect(fill='transparent', color='transparent'),
          strip.text=element_text(color='black',vjust=-1.8,hjust=0.02,size=20))
  
}


