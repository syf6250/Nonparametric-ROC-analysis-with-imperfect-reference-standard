#single biomarker simualtion(normal settings)
library(mgcv)
library(spatstat.univar)
library(rootSolve)
source('FUN-simu.R')
source('FUN-np.R')
source('FUN.R')
library(doParallel)
source('FUN-para.R')
source('FUN-boot.R')
library(openxlsx2)

nseq <- 300#c(100,300,500)
senseq <- c(0.95,0.9,0.75)
nu.seq <- 3*10^seq(from=-1.5,to=2.5,length.out=10)
s <- 0.2#where to evaluate ROC
s.range <- c(0.1,0.3)#where to compute pAUC

RUN <- 1000

perf.allcase <- NULL;cover.allcase <- NULL
ptm <- proc.time()
for(j1 in 1:length(senseq)){
  pi0 <- computePi(senseq[j1],senseq[j1],0.4)$pi0
  pi1 <- computePi(senseq[j1],senseq[j1],0.4)$pi1
  for(j2 in 1:length(nseq)){
    n <- m <- nseq[j2]
    set.mat <- data.frame(se=rep(senseq[j1],3),n=rep(n,3),
                          method=c("proposed","NP","naive"))
    
    cluster <- makeCluster(min(77,detectCores()-2))
    registerDoParallel(cluster)
    res.l <- list()
    res.l <- foreach(run = 1:RUN,
                     .packages = c("mgcv","spatstat.univar","rootSolve"))%dopar%{
                       pack.uni.boot(run,nu.seq,n,m,pi0,pi1,
                                     K=50,knots="q",datatype="norm",
                                     s=s,s.range=s.range,B=500,btune=F)
                     }
    stopCluster(cl = cluster)
    save(res.l,file=paste0("n",n,"_RUN",RUN,"_se",senseq[j1],"_normBoot.Rdata"))
    sumout <- sum.pack.uni.boot(res.l)
    print(sumout$conv.prop.Boot)
    perf.allcase <- rbind(perf.allcase,cbind(set.mat,sumout$perf))
    cover.allcase <- rbind(cover.allcase,cbind(set.mat,sumout$cover))
  }
}
proc.time() - ptm

cover.allcase 
#write_xlsx(perf.allcase,"perf-allcase-normal.xlsx")
write_xlsx(cover.allcase,"cover-allcase-normal.xlsx")
