#two biomarkers simualtion-sep(normal settings)
library(mvtnorm)
library(mgcv)
library(spatstat.univar)
library(rootSolve)
source('FUN-simu.R')
source('FUN-np.R')
source('FUN.R')
library(doParallel)
source('FUN-para.R')
library(openxlsx2)

nseq <- c(100,300,500)
senseq <- c(0.95,0.9,0.75)
nu.seq <- 3*10^seq(from=-1.5,to=2.5,length.out=10)
s <- 0.2#where to evaluate ROC
s.range <- c(0.1,0.3)#where to compute pAUC

rho <- 0.2#corr=0.2

RUN <- 1000
perf.allcase <- NULL;tune.allcase <- NULL
perf.allcase1 <- NULL;perf.allcase2 <- NULL
ptm <- proc.time()
for(j1 in 1:length(senseq)){
  pi0 <- computePi(senseq[j1],senseq[j1],0.4)$pi0
  pi1 <- computePi(senseq[j1],senseq[j1],0.4)$pi1
  for(j2 in 1:length(nseq)){
    n <- m <- nseq[j2]
    set.mat <- data.frame(se=rep(senseq[j1],3),n=rep(n,3),
                          method=c("CV5","NP","naive"))
    
    cluster <- makeCluster(min(77,detectCores()-2))
    registerDoParallel(cluster)
    res.l <- list()
    res.l <- foreach(run = 1:RUN,
                     .packages = c("mgcv","spatstat.univar","rootSolve","mvtnorm"))%dopar%{
                       pack.bi(run,nu.seq,n,m,pi0,pi1,rho=rho,
                               K=50,knots="q",datatype="norm",
                               s=s,s.range=s.range)
                     }
    stopCluster(cl = cluster)
    sumout <- sum.pack.bi(res.l)
    print(sumout$conv)
    # save(sumout,file=paste0("n",n,"_RUN",RUN,"_se",senseq[j1],"_norm2.Rdata"))
    perf.allcase <- rbind(perf.allcase,cbind(set.mat,sumout$perf))
    # tune.allcase <- rbind(tune.allcase,cbind(set.mat[1:2,1:2],
    #                                          rbind(sumout$id.cv5.1,sumout$id.cv5.2)))
    
    perf.allcase1 <- rbind(perf.allcase1,cbind(set.mat,sumout$perf1))
    perf.allcase2 <- rbind(perf.allcase2,cbind(set.mat,sumout$perf2))
  }
}
proc.time() - ptm

perf.allcase

write_xlsx(perf.allcase,"perf-allcase-normal2.xlsx") 

write_xlsx(perf.allcase1,"perf-1-h-normal2.xlsx") 
write_xlsx(perf.allcase2,"perf-2-h-normal2.xlsx") 


