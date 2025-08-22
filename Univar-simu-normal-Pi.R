#single biomarker simualtion(normal settings: sensitivity analysis for misspecified pi0 pi1)
library(mgcv)
library(spatstat.univar)
library(rootSolve)
source('FUN-simu.R')
source('FUN-np.R')
source('FUN.R')
library(doParallel)
source('FUN-para.R')
source('FUN-sensi.R')
library(openxlsx2)


nseq <- c(100,300,500)
Pi0seq <- Pi1seq <- c(1,0.95,0.9,0.85,0.8)
nu.seq <- 3*10^seq(from=-1.5,to=2.5,length.out=10)
s <- 0.2#where to evaluate ROC
s.range <- c(0.1,0.3)#where to compute pAUC

RUN <- 1000

perf.allcase <- NULL;tune.allcase <- NULL
ptm <- proc.time()
for(j1 in 1:length(Pi0seq)){
  pi0 <- Pi0seq[j1]
  pi1 <- Pi1seq[j1]
  for(j2 in 1:length(nseq)){
    n <- m <- nseq[j2]
    set.mat <- data.frame(pi0=rep(pi0,4),pi1=rep(pi1,4),n=rep(n,4),
                          method=c("AIC","CV5","NP","naive"))
    
    cluster <- makeCluster(min(77,detectCores()-2))
    registerDoParallel(cluster)
    res.l <- list()
    res.l <- foreach(run = 1:RUN,
                     .packages = c("mgcv","spatstat.univar","rootSolve"))%dopar%{
                       pack.uni.asmp2(run,nu.seq,n,m,pi0,pi1,
                                      K=50,knots="q",datatype="norm",
                                      s=s,s.range=s.range,
                                      pi0.true=0.9,pi1.true=0.9)
                     }
    stopCluster(cl = cluster)
    sumout <- sum.pack.uni.assump1(res.l)
    print(sumout$conv)
    #save(sumout,file=paste0("n",n,"_RUN",RUN,"_pi",pi0,pi1,"_norm-sensiPi.Rdata"))
    perf.allcase <- rbind(perf.allcase,cbind(set.mat,sumout$perf))
    tune.allcase <- rbind(tune.allcase,cbind(set.mat[1:2,],rbind(sumout$id.aic,sumout$id.cv5)))
  }
}
proc.time() - ptm

perf.allcase
#tune.allcase
write_xlsx(perf.allcase,paste0("pi",pi0,pi1,"-perf-normal-sensi.xlsx")) 





