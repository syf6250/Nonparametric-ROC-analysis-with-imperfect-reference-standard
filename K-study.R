#K-study: single biomarker simualtion(normal settings)
library(mgcv)
library(spatstat.univar)
library(rootSolve)
source('FUN-simu.R')
source('FUN-np.R')
source('FUN.R')
library(doParallel)
source('FUN-para.R')
library(openxlsx2)


nseq <- 500#c(100,300,500)
senseq <- 0.9#c(0.95,0.9,0.75)
Kseq <- c(10,30,50,70,90)
nu.seq <- 3*10^seq(from=-1.5,to=2.5,length.out=10)
s <- 0.2#where to evaluate ROC
s.range <- c(0.1,0.3)#where to compute pAUC

RUN <- 1000

perf.allcase <- NULL;tune.allcase <- NULL
ptm <- proc.time()
for(j1 in 1:length(senseq)){
  pi0 <- computePi(senseq[j1],senseq[j1],0.4)$pi0
  pi1 <- computePi(senseq[j1],senseq[j1],0.4)$pi1
  for(j2 in 1:length(nseq)){
    n <- m <- nseq[j2]
    set.mat <- data.frame(se=rep(senseq[j1],length(Kseq)),n=rep(n,length(Kseq)),K=Kseq)
    #a unified function to conduct EM for each K
    cluster <- makeCluster(min(77,detectCores()-2))
    registerDoParallel(cluster)
    res.l <- list()
    res.l <- foreach(run = 1:RUN,
                     .packages = c("mgcv","spatstat.univar","rootSolve"))%dopar%{
                       pack.Kstudy(run,nu.seq,n,m,pi0,pi1,
                                   Kseq,knots="q",datatype="norm",
                                   s=s,s.range=s.range)
                     }
    stopCluster(cl = cluster)
    sumout <- sum.pack.Kstudy(res.l,Kseq)
    #save(sumout,file=paste0("n",n,"_RUN",RUN,"_se",senseq[j1],"_norm-Kstudy.Rdata"))
    perf.allcase <- rbind(perf.allcase,cbind(set.mat,sumout$perf))
    tune.allcase <- rbind(tune.allcase,cbind(set.mat,sumout$id))
  }
}
proc.time() - ptm

perf.allcase
#tune.allcase

write_xlsx(perf.allcase,"perf-allcase-normal-Kstudy.xlsx") 
#write_xlsx(tune.allcase,"tune-allcase-normal-Kstudy.xlsx") 











#K-study: single biomarker simualtion(gamma settings)
library(mgcv)
library(spatstat.univar)
library(rootSolve)
source('FUN-simu.R')
source('FUN-np.R')
source('FUN.R')
library(doParallel)
source('FUN-para.R')
library(openxlsx2)


nseq <- 500#c(100,300,500)
senseq <- 0.9#c(0.95,0.9,0.75)
Kseq <- c(10,30,50,70,90)
nu.seq <- 3*10^seq(from=-2,to=2.5,length.out=10)
s <- 0.2#where to evaluate ROC
s.range <- c(0.1,0.3)#where to compute pAUC

RUN <- 1000

perf.allcase <- NULL;tune.allcase <- NULL
ptm <- proc.time()
for(j1 in 1:length(senseq)){
  pi0 <- computePi(senseq[j1],senseq[j1],0.4)$pi0
  pi1 <- computePi(senseq[j1],senseq[j1],0.4)$pi1
  for(j2 in 1:length(nseq)){
    n <- m <- nseq[j2]
    set.mat <- data.frame(se=rep(senseq[j1],length(Kseq)),n=rep(n,length(Kseq)),K=Kseq)
    #a unified function to conduct EM for each K
    cluster <- makeCluster(min(77,detectCores()-2))
    registerDoParallel(cluster)
    res.l <- list()
    res.l <- foreach(run = 1:RUN,
                     .packages = c("mgcv","spatstat.univar","rootSolve"))%dopar%{
                       pack.Kstudy(run,nu.seq,n,m,pi0,pi1,
                                   Kseq,knots="q",datatype="gamma",
                                   s=s,s.range=s.range)
                     }
    stopCluster(cl = cluster)
    sumout <- sum.pack.Kstudy(res.l,Kseq)
    #save(sumout,file=paste0("n",n,"_RUN",RUN,"_se",senseq[j1],"_gamma-Kstudy.Rdata"))
    perf.allcase <- rbind(perf.allcase,cbind(set.mat,sumout$perf))
    tune.allcase <- rbind(tune.allcase,cbind(set.mat,sumout$id))
  }
}
proc.time() - ptm

perf.allcase
#tune.allcase

write_xlsx(perf.allcase,"perf-allcase-gamma-Kstudy.xlsx") 
#write_xlsx(tune.allcase,"tune-allcase-gamma-Kstudy.xlsx") 


