ES.cov.test <- function(pheno,snp){
  y <- pheno
  snps <- snp
  nsnp <- dim(snps)[1]
  time <- dim(y)[2]
  p.matrix <- matrix(rep(0,nsnp*time),nrow=time)
  fdr.matrix <- matrix(rep(0,nsnp*time),nrow=time)
  for(i in 1:time){
    for(j in 1:nsnp){
      tmpsnp <- snps[j,]
      miss <- which(is.na(tmpsnp))
      if(length(miss)>0){
        newsnp <- tmpsnp[-miss]
        yy <- y[-miss,i]
      }else{
        newsnp <- tmpsnp
        yy <- y[,i]
      }
      symbol <- names(table(as.character(unlist(c(newsnp)))))
      index1 <- which(newsnp==symbol[1])
      y1 <- yy[index1]
      index0 <- which(newsnp==symbol[2])
      y0 <- yy[index0]
      #snpdata <- data.frame(X=c(y1,y0),A=factor(c(rep(1,length(index1)),
                                                  #rep(2,length(index0)))))
      #snp.aov <- aov(X ~ A, data=snpdata)
      var.t <- var.test(y0,y1)
      if(var.t$p.value>0.05)
        var.i <- TRUE
      else
        var.i <- FALSE
      p.value <- t.test(y0,y1,var.equal = var.i)$p.value#summary(snp.aov)[[1]][[1,"Pr(>F)"]] 
      p.matrix[i,j] <- p.value
    }
    fdr.matrix[i,] <- p.adjust(p.matrix[i,],method="fdr")      
  }
  
  colnames(p.matrix)<- rownames(snps)
  colnames(fdr.matrix)<- rownames(snps)
  
  logp <- -log10(p.matrix)
  thre1 <- -log10(0.05/nsnp)
  thre2 <- -log10(0.01/nsnp)
  return(list(p.value=p.matrix,fdr=fdr.matrix,logp=logp,thre1=thre1,thre2=thre2))
}