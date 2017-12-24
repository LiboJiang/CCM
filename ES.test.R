
source("ES.load.R")
source("ES.util.R")
source("ES.covar.R")
source("ES.plot.R")
source("ES.util.R")
source("ES.ODE.R")
source("ES.cov.R")
source("ES.optim.R")
source("ES.fun.curve.R")
source("ES.fun.optim.R")
library(mvtnorm)
############################### data load and process ###################################

dat <- SE.load(e.file="./data/coculture-Ecoli.csv",s.file="./data/coculture-Saureus.csv",
                    e.snp.file="./data/gen_ecoli.txt",s.snp.file="./data/gen_saureus.txt",
                    ei.file="data/monocultrue-Ecoli.csv",si.file="./data/monocultrue-Saureus.csv")

#ES.plot.ind(dat) #each individual growth curve
#ES.plink(SNP=dat$s.snp,fil="S") ####snps were tranlated to plink format
#ES.qq.plot(dat) ##data normal distribution test
dat1 <- ES.adjust.s(dat,eefile="Structure/E.coli/E/Eput_simple.5.meanQ",ssfile="Structure/S.aureus/S/Sput_simple.8.meanQ")
#phenotype adjust exclude population effect
#ES.plot.ind(dat1)
#ES.ind.adjust.plot(alldat=dat1) #residual distribution plot

###########################T-test analysis for single point ######################

E.i <- ES.cov.test(pheno=dat1$ep.pi,snp=dat1$e.snp)
E <- ES.cov.test(pheno=dat1$ep.p,snp=dat1$e.snp)
S.i <- ES.cov.test(pheno=dat1$sp.pi,snp=dat1$s.snp)
S <- ES.cov.test(pheno=dat1$sp.p,snp=dat1$s.snp)

save(E.i,file="E-i.RData"); save(E,file="E.RData");
save(S.i,file="S-i.RData"); save(S,file="S.RData");

#ES.pvalue.qq.plot(Pei=E.i$p.value,Pe=E$p.value,Psi=S.i$p.value,Ps=S$p.value) #Q-Q plot for p-values
#ES.single.man(P=E,seg=9,nchr=12,dist=0.4,filename="Figure-E-man.pdf",filename1="E-sig.txt",thre1=4.2)
#ES.single.man(P=E.i,seg=9,nchr=12,dist=0.4,filename="Figure-Ei-man.pdf",filename1="Ei-sig.txt",thre1=4.2) 
#ES.single.man(P=S,seg=9,nchr=12,dist=0.2,filename="Figure-S-man.pdf",filename1="S-sig.txt",thre1=3.8) 
#ES.single.man(P=S.i,seg=9,nchr=12,dist=0.2,filename="Figure-Si-man.pdf",filename1="Si-sig.txt",thre1=3.8) 


########################### system mapping ######################


fp.plot_ind_curves(dat1, pheno1=dat1$ep.p,pheno2=dat1$sp.p,rows=9, cols=5,filename="Figure-ES-raw.pdf" )

fp.plot_ind_curves.com( dat1, pheno1=dat1$ep.p,pheno2=dat1$sp.p,rows=9, cols=5,filename="Figure-ES-fit.pdf" )


H0 <- curve.get_est_param(dat1,pheno1=dat1$ep.p,pheno2=dat1$sp.p) # parameter initialization for system mapping

#ES.plot.mean(dat=dat1,allpar=H0[6:11])

load("E-i.RData")
load("E.RData")
load("S-i.RData")
load("S.RData")


dat <- snp.select(dat,E.i,E,S.i,S) #  combining snps
#load("dat.RData")


ret <- ES.H1.map(dat,pheno=cbind(dat$ep.p,dat$sp.p),fsnp=dat$fsnp,H0,select.snp=c(1:20)) # Genome-genome interaction scan

#####functional mapping#######################


H0ep <- funmap.curve(dat1,pheno=dat1$ep.p) #parameter initialization for function mapping 
H0epi <- funmap.curve(dat1,pheno=dat1$ep.pi)
H0sp <- funmap.curve(dat1,pheno=dat1$sp.p)
H0spi <- funmap.curve(dat1,pheno=dat1$sp.pi)


H1ep <- funmap(dat1,pheno=dat1$ep.p,H0=H0ep,SNP=dat1$e.snp) ## Genome scan
H1epi <- funmap(dat1,pheno=dat1$ep.pi,H0=H0epi,SNP=dat1$e.snp)
H1sp <- funmap(dat1,pheno=dat1$sp.p,H0=H0sp,SNP=dat1$s.snp)
H1spi <- funmap(dat1,pheno=dat1$sp.pi,H0=H0spi,SNP=dat1$s.snp)
