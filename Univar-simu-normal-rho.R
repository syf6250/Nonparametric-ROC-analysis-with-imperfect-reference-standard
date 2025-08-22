#single biomarker simualtion(normal settings: sensitivity analysis)
library(mgcv)
library(spatstat.univar)
library(rootSolve)
source('FUN-simu.R')
source('FUN-np.R')
source('FUN.R')
source('FUN-sensi.R')
library(doParallel)
source('FUN-para.R')
library(openxlsx2)

rhoseq <- c(0,0.1,0.2,0.3)
varX1=0.5;varX2=0.5;vareps=0.5
prev=0.4
coef=list(bX1=1,bGT=1,bX2=1,bGR=4.78,aR=-2.39)

sp <- integrate(function(x){1/(1+exp(coef$bX2*x+coef$aR))*dnorm(x,sd=sqrt(varX2))},
                lower = -Inf,upper = Inf)$value
se <- integrate(function(x){1/(1+exp(-coef$bX2*x-coef$aR-coef$bGR))*dnorm(x,sd=sqrt(varX2))},
                lower = -Inf,upper = Inf)$value
pi0 <- computePi(se,sp,prev)$pi0
pi1 <- computePi(se,sp,prev)$pi1

nseq <- c(100,300,500)
nu.seq <- 3*10^seq(from=-1.5,to=2.5,length.out=10)


RUN <- 1000
perf.allcase <- NULL;tune.allcase <- NULL
for(j1 in 1:length(rhoseq)){
  rho <- rhoseq[j1]
  for(j2 in 1:length(nseq)){
    n <- m <- nseq[j2]
    set.mat <- data.frame(rho=rep(rho,4),n=rep(n,4),
                          method=c("AIC","CV5","NP","naive"))
    cluster <- makeCluster(min(77,detectCores()-2))
    registerDoParallel(cluster)
    res.l <- list()
    res.l <- foreach(run = 1:RUN,
                     .packages = c("mgcv","spatstat.univar","rootSolve","mvtnorm"))%dopar%{
                       pack.uni.asmp1(run,nu.seq=nu.seq,n=n,m=m,
                                      pi0=pi0,pi1=pi1,K=50,knots="q",s=0.2,s.range=c(0.1,0.3),
                                      rho=rho,varX1=varX1,varX2=varX2,vareps=vareps,prev=prev,coef=coef)
                     }
    stopCluster(cl = cluster)
    sumout <- sum.pack.uni.assump1(res.l)
    print(sumout$conv)
    save(sumout,file=paste0("n",n,"_RUN",RUN,"_rho",rho,"_norm-sensi.Rdata"))
    perf.allcase <- rbind(perf.allcase,cbind(set.mat,sumout$perf))
    tune.allcase <- rbind(tune.allcase,cbind(set.mat[1:2,],rbind(sumout$id.aic,sumout$id.cv5)))
  }
}
proc.time() - ptm

perf.allcase
#tune.allcase
write_xlsx(perf.allcase,paste0("rho",rho,"-perf-normal-sensi.xlsx")) 
