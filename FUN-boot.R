###functions for bootstrapp intervals


pack.uni.boot <- function(run,nu.seq,n,m,pi0,pi1,
                          K=50,knots="q",datatype="norm",
                          s=0.2,s.range=c(0.1,0.3),B=500,btune=F){
  if(datatype=="norm"){
    ###settings
    mu <- c(0,1)#mean of R=0 & 1; the second element should be larger!
    sigma <- c(1,1)#sd of R=0 & 1
    eval.points <- seq(from=-1,to=2,length.out=20)
    ROC.true <- 1-pnorm(qnorm(1-s,mu[1],sigma[1]),mu[2],sigma[2])
    AUC.true <- integrate(function(c){(1-pnorm(c,mu[2],sigma[2]))*dnorm(c,mu[1],sigma[1])},
                          lower=-50,upper=50)$value
    true.cutoff <- uniroot.all(function(t){log(dnorm(t,mu[2],sigma[2])/dnorm(t,mu[1],sigma[1]))},c(-5,5))
    YI.true <- pnorm(true.cutoff,mu[1],sigma[1])-pnorm(true.cutoff,mu[2],sigma[2])
    ROC.fun0 <- function(s){1-pnorm(qnorm(1-s,mu[1],1),mu[2],1)}
    pAUC.true <- (integrate(ROC.fun0,s.range[1],s.range[2])$value)/(s.range[2]-s.range[1])
    true.all <- c(ROC.true,AUC.true,true.cutoff,YI.true,pAUC.true)
    ###data
    set.seed(run*77+7)
    mydt <- gnormmix(n,m,pi0,pi1,mu=mu,sigma=sigma)
  }else if(datatype=="gamma"){
    shape <- c(2,3)
    rate <- c(1,1)
    eval.points <- seq(from=0.5,to=5.5,length.out=20)
    ROC.true <- 1-pgamma(qgamma(1-s,shape=shape[1],rate=rate[1]),shape=shape[2],rate=rate[2])
    AUC.true <- integrate(function(c){(1-pgamma(c,shape=shape[2],rate=rate[2]))*dgamma(c,shape=shape[1],rate=rate[1])},
                          lower=0,upper=50)$value
    eval.true <- log(dgamma(eval.points,shape=shape[2],rate=rate[2])/dgamma(eval.points,shape=shape[1],rate=rate[1]))
    true.cutoff <- uniroot.all(function(t){log(dgamma(t,shape=shape[2],rate=rate[2])/dgamma(t,shape=shape[1],rate=rate[1]))},c(0,5))
    YI.true <- pgamma(true.cutoff,shape=shape[1],rate=rate[1])-pgamma(true.cutoff,shape=shape[2],rate=rate[2])
    ROC.fun0 <- function(s){1-pgamma(qgamma(1-s,shape=shape[1],rate=rate[1]),shape=shape[2],rate=rate[2])}
    pAUC.true <- (integrate(ROC.fun0,s.range[1],s.range[2])$value)/(s.range[2]-s.range[1])
    true.all <- c(ROC.true,AUC.true,true.cutoff,YI.true,pAUC.true)
    ###data
    set.seed(run*77+7)
    mydt <- ggammamix(n,m,pi0,pi1,shape=shape,rate=rate)
  } 
  
  ###estimation
  ini.G <- NULL#no initial
  ic.out <- IC.DRest.EM(mydt,pi0,pi1,
                        nu=nu.seq,k=K,ord=4,ord.pen=2,thres = 1e-5,
                        knots=knots,ini.G=ini.G)
  id.aic <- ic.out$tune.selection$aic[1]#is.na(id.aic)?
  conv.aic <- ic.out$tune.res$conv[id.aic]
  if(isTRUE(conv.aic)){
    sumout <- sum.DRest0.EM(mydt,ic.out$all.list[[id.aic]],s=s,s.range=s.range)
    stats.aic <- sumout$stats.df
  }else{
    stats.aic <- data.frame(ROC=NA,AUC=NA,cutoff=NA,Youden=NA,pAUC=NA)
  }
  paras <- list();paras$true <- true.all;paras$aic <- stats.aic
  
  #5 cv
  cvout <- cv.DRest.EM(mydt,pi0,pi1,nfolds=5,seed=123,
                       nu=nu.seq,k=K,ord=4,ord.pen=2,
                       maxit=500,thres=1e-5,knots=knots,ini.G=ini.G)
  id.cv5 <- cvout$tune.id
  if(length(id.cv5)>0){
    conv.cv5 <- cvout$res$conv[id.cv5]
  }else{
    conv.cv5 <- F
    id.cv5 <- NA
  }
  if(isTRUE(conv.cv5)){
    stats.cv5 <- try(sum.DRest0.EM(mydt,ic.out$all.list[[id.cv5]],
                                   s=s,s.range=s.range)$stats.df,silent = T)
    if(class(stats.cv5)=="try-error"){
      stats.cv5 <- data.frame(ROC=NA,AUC=NA,cutoff=NA,Youden=NA,pAUC=NA)
    }
  }else{
    stats.cv5 <- data.frame(ROC=NA,AUC=NA,cutoff=NA,Youden=NA,pAUC=NA)
  }
  paras$cv5 <- stats.cv5
  #naive method (using cv5 tuning)
  if(isTRUE(conv.cv5)){
    test.naive <- DRest0.EM(mydt,pi0=1,pi1=1,
                            nu=nu.seq[id.cv5],k=K,ord=4,ord.pen=2,
                            maxit=500,thres=1e-5,knots=knots)
    if(isTRUE(test.naive$conv)){
      stats.naive <- sum.DRest0.EM(mydt,test.naive,s=s,s.range=s.range)$stats.df
    }else{
      stats.naive <- data.frame(ROC=NA,AUC=NA,cutoff=NA,Youden=NA,pAUC=NA)
    }
    paras$naive <- stats.naive
  }else{
    paras$naive <- data.frame(ROC=NA,AUC=NA,cutoff=NA,Youden=NA,pAUC=NA)
  }
  #np method
  mydt2 <- mydt
  if(datatype=="norm"){
    mydt2$biomarker <- exp(mydt2$biomarker)/(1+exp(mydt2$biomarker))
    npout <- np.roc(mydt2,pi0,pi1,t=s,s.range=s.range)
    npout$cutoff <- log(npout$cutoff/(1-npout$cutoff))
  }else if(datatype=="gamma"){
    mydt2$biomarker <- (mydt2$biomarker)/(1+(mydt2$biomarker))
    npout <- np.roc(mydt2,pi0,pi1,t=s,s.range=s.range)
    npout$cutoff <- npout$cutoff/(1-npout$cutoff)
  }
  paras$np <- npout 
  
  finout <- list(paras=paras,conv=data.frame(aic=conv.aic,cv5=conv.cv5),
                 tune.id=data.frame(aic=id.aic,cv5=id.cv5))
  
  ###bootstrap
  if(isFALSE(finout$conv$cv5)){
    finout$conv.boot <- 0
  }else{
    if(isFALSE(btune)){
      tune.work <- nu.seq[id.cv5]
    }else{tune.work <- nu.seq}
    
    allest <- allest.np <- allest.naive <- NULL
    conv.prop.boot <- vector(length = B)
    for (jj in 1:B) {
      #bootstrap sample
      set.seed(jj)
      R0.selec.id <- sample(c(1:n),n,T)
      R1.selec.id <- sample(c((n+1):(n+m)),m,T)
      mydt.selec <- mydt[c(R0.selec.id,R1.selec.id),]
      
      #estimation: proposed
      ic.out <- IC.DRest.EM(mydt.selec,pi0,pi1,
                            nu=tune.work,k=K,ord=4,ord.pen=2,thres = 1e-5,
                            knots=knots,ini.G=NULL)#no initial#knots are not the same as full data
      id.aic <- ic.out$tune.selection$aic[1]
      if(is.na(id.aic)){
        conv.aic <- F
      }else{
        conv.aic <- ic.out$tune.res$conv[id.aic]
        conv.prop.boot[jj] <- conv.aic
      }
      if(isTRUE(conv.aic)){
        sumout <- sum.DRest0.EM(mydt.selec,ic.out$all.list[[id.aic]],s=s,s.range=s.range)
        stats.aic <- sumout$stats.df
      }else{
        stats.aic <- data.frame(ROC=NA,AUC=NA,cutoff=NA,Youden=NA,pAUC=NA)
      }
      allest <- rbind(allest,stats.aic)
      
      #estimation: naive
      test.naive <- try(DRest0.EM(mydt.selec,pi0=1,pi1=1,
                                  nu=nu.seq[id.cv5],k=K,ord=4,ord.pen=2,
                                  maxit=500,thres=1e-5,knots=knots),silent = T)
      if(class(test.naive)!="try-error"){
        if(isTRUE(test.naive$conv)){
          stats.naive <- sum.DRest0.EM(mydt.selec,test.naive,s=s,s.range=s.range)$stats.df
        }
      }else{
        stats.naive <- data.frame(ROC=NA,AUC=NA,cutoff=NA,Youden=NA,pAUC=NA)
      }
      allest.naive <- rbind(allest.naive,stats.naive)
      
      #estimation: np
      if(datatype=="norm"){
        mydt.selec$biomarker <- exp(mydt.selec$biomarker)/(1+exp(mydt.selec$biomarker))
        npout <- np.roc(mydt.selec,pi0,pi1,t=s,s.range=s.range)
        npout$cutoff <- log(npout$cutoff/(1-npout$cutoff))
      }else if(datatype=="gamma"){
        mydt.selec$biomarker <- (mydt.selec$biomarker)/(1+(mydt.selec$biomarker))
        npout <- np.roc(mydt.selec,pi0,pi1,t=s,s.range=s.range)
        npout$cutoff <- npout$cutoff/(1-npout$cutoff)
      }
      allest.np <- rbind(allest.np,npout)
    }
    
    #summary
    finout$conv.boot <- mean(conv.prop.boot)
    bootsum <- rbind(apply(allest,2,quantile,na.rm=T,probs=c(0.025,0.975,0.05,0.95)),
                     apply(allest,2,mean,na.rm=T),apply(allest,2,sd,na.rm=T))
    bootsum.naive <- rbind(apply(allest.naive,2,quantile,na.rm=T,probs=c(0.025,0.975,0.05,0.95)),
                     apply(allest.naive,2,mean,na.rm=T),apply(allest.naive,2,sd,na.rm=T))
    bootsum.np <- rbind(apply(allest.np,2,quantile,na.rm=T,probs=c(0.025,0.975,0.05,0.95)),
                     apply(allest.np,2,mean,na.rm=T),apply(allest.np,2,sd,na.rm=T))
    rownames(bootsum)[c(5,6)] <- rownames(bootsum.naive)[c(5,6)] <- rownames(bootsum.np)[c(5,6)] <- c("mean","sd")
    finout$bootsum <- bootsum
    finout$bootsum.naive <- bootsum.naive
    finout$bootsum.np <- bootsum.np
    #CI 95%
    finout$cover <- cover.analysis(true.all,finout$paras$cv5,bootsum,debias=F)
    finout$cover.naive <- cover.analysis(true.all,finout$paras$naive,bootsum.naive,debias=F)
    finout$cover.np <- cover.analysis(true.all,finout$paras$np,bootsum.np,debias=F)
  }
  
  return(finout)
}

cover.analysis <- function(true.all,est.all,bootsum,debias=F){
  if(isFALSE(debias)){
    tominus <- 0
  }else{
    tominus <- bootsum['mean',]-est.all
  }
  ll <- 2*est.all-bootsum['97.5%',]-tominus
  uu <- 2*est.all-bootsum['2.5%',]-tominus
  out <- ifelse(ll<=true.all & uu>true.all,T,F)
}

sum.pack.uni.boot <- function(res.l){
  RUN <- length(res.l)
  true.mat <- do.call(rbind,lapply(res.l, function(x){x$paras$true}))
  all.true <- apply(true.mat, 2, median)#
  out <- list(para.true=all.true,
              res.df.cv5=do.call(rbind,lapply(res.l, function(x){x$paras$cv5})),
              res.df.np=do.call(rbind,lapply(res.l, function(x){x$paras$np})),
              res.df.naive=do.call(rbind,lapply(res.l, function(x){x$paras$naive})))
  id.cv5 <- sapply(res.l,function(x){x$tune.id$cv5})
  out$id.cv5 <- table(factor(id.cv5,levels = 1:length(nu.seq)))/RUN
  conv_cv5 <- sapply(res.l, function(x){x$conv$cv5})
  out$conv.cv5 <- mean(conv_cv5)
  out$perf <- rbind(perf.15(out$res.df.cv5,all.true),
                    perf.15(out$res.df.np,all.true),
                    perf.15(out$res.df.naive,all.true))
  rownames(out$perf) <- c("cv5","np","naive")
  
  convprop <- sapply(res.l, function(x){x$conv.boot})
  out$conv.prop.Boot <- data.frame(mean=mean(convprop),median=median(convprop))
  
  cover <- rbind(apply(do.call(rbind,lapply(res.l, function(x){x$cover})),2,mean),
                 apply(do.call(rbind,lapply(res.l, function(x){x$cover.np})),2,mean),
                 apply(do.call(rbind,lapply(res.l, function(x){x$cover.naive})),2,mean))
  rownames(cover) <- c("proposed","np","naive")
  out$cover <- cover
  return(out)
}


pack.bi.boot <- function(run,nu.seq,n,m,pi0,pi1,rho,
                         K=50,knots="q",datatype="norm",
                         s=0.2,s.range=c(0.1,0.3),B=500,btune=F){
  if(datatype=="norm"){
    ###settings
    require(mvtnorm)
    mu1 <- c(2,1);mu0 <- c(0,0)
    #1: diff of ROC
    ROC.true.1 <- 1-pnorm(qnorm(1-s,mu0[1],1),mu1[1],1)
    ROC.true.2 <- 1-pnorm(qnorm(1-s,mu0[2],1),mu1[2],1)
    DROC.true <- ROC.true.1-ROC.true.2
    #2: diff of AUC
    AUC.true.1 <- integrate(function(c){(1-pnorm(c,mu1[1],1))*dnorm(c,mu0[1],1)},lower=-50,upper=50)$value
    AUC.true.2 <- integrate(function(c){(1-pnorm(c,mu1[2],1))*dnorm(c,mu0[2],1)},lower=-50,upper=50)$value
    DAUC.true <- AUC.true.1-AUC.true.2
    #3: diff of Youden indices
    true.cutoff.1 <- uniroot.all(function(t){log(dnorm(t,mu1[1],1)/dnorm(t,mu0[1],1))},c(-5,5))
    true.cutoff.2 <- uniroot.all(function(t){log(dnorm(t,mu1[2],1)/dnorm(t,mu0[2],1))},c(-5,5))
    YI.true.1 <- pnorm(true.cutoff.1,mu0[1],1)-pnorm(true.cutoff.1,mu1[1],1)
    YI.true.2 <- pnorm(true.cutoff.2,mu0[2],1)-pnorm(true.cutoff.2,mu1[2],1)
    DYI.true <- YI.true.1-YI.true.2
    #4: diff of pAUC
    ROC.fun0 <- function(s){1-pnorm(qnorm(1-s,mu0[1],1),mu1[1],1)}
    pAUC.true.1 <- (integrate(ROC.fun0,s.range[1],s.range[2])$value)/(s.range[2]-s.range[1])
    ROC.fun0 <- function(s){1-pnorm(qnorm(1-s,mu0[2],1),mu1[2],1)}
    pAUC.true.2 <- (integrate(ROC.fun0,s.range[1],s.range[2])$value)/(s.range[2]-s.range[1])
    DpAUC.true <- pAUC.true.1-pAUC.true.2
    ###data
    set.seed(run*77+7)
    mydt <- gnormmix2(n,m,pi0,pi1,
                      mu1=mu1,mu0=mu0,
                      var1=matrix(c(1,rho,rho,1),2,2),
                      var0=matrix(c(1,rho,rho,1),2,2))
  }else if(datatype=="gamma"){
    ###settings
    #1: diff of ROC
    ROC.true.1 <- 1-pgamma(qgamma(1-s,2,1),4,1)
    ROC.true.2 <- 1-pgamma(qgamma(1-s,2,1),3,1)
    DROC.true <- ROC.true.1-ROC.true.2
    #2: diff of AUC
    AUC.true.1 <- integrate(function(c){(1-pgamma(c,4,1))*dgamma(c,2,1)},lower=0,upper=100)$value
    AUC.true.2 <- integrate(function(c){(1-pgamma(c,3,1))*dgamma(c,2,1)},lower=0,upper=100)$value
    DAUC.true <- AUC.true.1-AUC.true.2
    #3: diff of Youden indices
    true.cutoff.1 <- uniroot.all(function(t){log(dgamma(t,4,1)/dgamma(t,2,1))},c(0,50))
    true.cutoff.2 <- uniroot.all(function(t){log(dgamma(t,3,1)/dgamma(t,2,1))},c(0,50))
    YI.true.1 <- pgamma(true.cutoff.1,2,1)-pgamma(true.cutoff.1,4,1)
    YI.true.2 <- pgamma(true.cutoff.2,2,1)-pgamma(true.cutoff.2,3,1)
    DYI.true <- YI.true.1-YI.true.2
    #4: diff of pAUC
    ROC.fun0 <- function(s){1-pgamma(qgamma(1-s,2,1),4,1)}
    pAUC.true.1 <- (integrate(ROC.fun0,s.range[1],s.range[2])$value)/(s.range[2]-s.range[1])
    ROC.fun0 <- function(s){1-pgamma(qgamma(1-s,2,1),3,1)}
    pAUC.true.2 <- (integrate(ROC.fun0,s.range[1],s.range[2])$value)/(s.range[2]-s.range[1])
    DpAUC.true <- pAUC.true.1-pAUC.true.2
    ###data
    set.seed(run*77+7)
    mydt <- ggammamix2(n,m,pi0,pi1,shape1=c(4,3),shape0=c(2,2),
                       rho=rho)
  }
  true.all <- c(DROC.true,DAUC.true,DYI.true,DpAUC.true)
  mydt.1 <- mydt[,1:2];mydt.1$biomarker <- mydt$biomarker1
  mydt.2 <- mydt[,1:2];mydt.2$biomarker <- mydt$biomarker2
  paras <- list();paras$true <- true.all
  

  ###estimation (part1)
  ini.G <- NULL
  cvout <- cv.DRest.EM(mydt.1,pi0,pi1,nfolds=5,seed=123,
                       nu=nu.seq,k=K,ord=4,ord.pen=2,
                       maxit=500,thres=1e-5,knots=knots,ini.G=ini.G)
  id.cv5.1 <- cvout$tune.id#
  if(length(id.cv5.1)>0){
    conv.cv5.1 <- cvout$res$conv[id.cv5.1]
  }else{
    conv.cv5.1 <- F
    id.cv5.1 <- NA
  }
  if(isTRUE(conv.cv5.1)){
    test <- try(DRest0.EM(mydt.1,pi0,pi1,
                          nu=nu.seq[id.cv5.1],k=K,ord=4,ord.pen=2,
                          maxit=500,thres=1e-5,knots=knots,ini.G=ini.G),silent = T)
    if(class(test)!="try-error"){
      stats.cv5.1 <- sum.DRest0.EM(mydt.1,test,s=s,s.range=s.range)$stats.df
    }else{
      stats.cv5.1 <- data.frame(ROC=NA,AUC=NA,cutoff=NA,Youden=NA,pAUC=NA)
    }
  }else{
    stats.cv5.1 <- data.frame(ROC=NA,AUC=NA,cutoff=NA,Youden=NA,pAUC=NA)
  }
  ###estimation (part2)
  cvout <- cv.DRest.EM(mydt.2,pi0,pi1,nfolds=5,seed=123,
                       nu=nu.seq,k=K,ord=4,ord.pen=2,
                       maxit=500,thres=1e-5,knots=knots,ini.G=ini.G)
  id.cv5.2 <- cvout$tune.id#
  if(length(id.cv5.2)>0){
    conv.cv5.2 <- cvout$res$conv[id.cv5.2]
  }else{
    conv.cv5.2 <- F
    id.cv5.2 <- NA
  }
  if(isTRUE(conv.cv5.2)){
    test <- try(DRest0.EM(mydt.2,pi0,pi1,
                          nu=nu.seq[id.cv5.2],k=K,ord=4,ord.pen=2,
                          maxit=500,thres=1e-5,knots=knots,ini.G=ini.G),silent = T)
    if(class(test)!="try-error"){
      stats.cv5.2 <- sum.DRest0.EM(mydt.2,test,s=s,s.range=s.range)$stats.df
    }else{
      stats.cv5.2 <- data.frame(ROC=NA,AUC=NA,cutoff=NA,Youden=NA,pAUC=NA)
    }
  }else{
    stats.cv5.2 <- data.frame(ROC=NA,AUC=NA,cutoff=NA,Youden=NA,pAUC=NA)
  }
  paras$cv5 <- (stats.cv5.1-stats.cv5.2)[-3]
  
  #naive method (using cv5 tuning)
  if(conv.cv5.1*conv.cv5.2==1){
    test.naive <- DRest0.EM(mydt.1,pi0=1,pi1=1,
                            nu=nu.seq[id.cv5.1],k=K,ord=4,ord.pen=2,
                            maxit=500,thres=1e-5,knots=knots)#或许也得用try
    if(isTRUE(test.naive$conv)){
      stats.naive.1 <- sum.DRest0.EM(mydt.1,test.naive,s=s,s.range=s.range)$stats.df
    }else{
      stats.naive.1 <- data.frame(ROC=NA,AUC=NA,cutoff=NA,Youden=NA,pAUC=NA)
    }
    test.naive <- DRest0.EM(mydt.2,pi0=1,pi1=1,
                            nu=nu.seq[id.cv5.2],k=K,ord=4,ord.pen=2,
                            maxit=500,thres=1e-5,knots=knots)
    if(isTRUE(test.naive$conv)){
      stats.naive.2 <- sum.DRest0.EM(mydt.2,test.naive,s=s,s.range=s.range)$stats.df
    }else{
      stats.naive.2 <- data.frame(ROC=NA,AUC=NA,cutoff=NA,Youden=NA,pAUC=NA)
    }
    paras$naive <- (stats.naive.1-stats.naive.2)[-3]
  }else{
    paras$naive <- data.frame(ROC=NA,AUC=NA,Youden=NA,pAUC=NA)
  }
  
  #np method
  mydt2.1 <- mydt.1; mydt2.2 <- mydt.2
  if(datatype=="norm"){
    mydt2.1$biomarker <- exp(mydt2.1$biomarker)/(1+exp(mydt2.1$biomarker))
    mydt2.2$biomarker <- exp(mydt2.2$biomarker)/(1+exp(mydt2.2$biomarker))
  }else if(datatype=="gamma"){
    mydt2.1$biomarker <- (mydt2.1$biomarker)/(1+(mydt2.1$biomarker))
    mydt2.2$biomarker <- (mydt2.2$biomarker)/(1+(mydt2.2$biomarker))
  }
  npout1 <- np.roc(mydt2.1,pi0,pi1,t=s,s.range=s.range)
  npout2 <- np.roc(mydt2.2,pi0,pi1,t=s,s.range=s.range)
  paras$np <- (npout1-npout2)[-3]
  finout <- list(paras=paras,conv=data.frame(cv5.1=conv.cv5.1,cv5.2=conv.cv5.2),
                 tune.id=data.frame(cv5.1=id.cv5.1,cv5.2=id.cv5.2))
  
  
  ###bootstrap
  convtem <- (finout$conv$cv5.1)*(finout$conv$cv5.2)
  if(convtem==0){
    finout$conv.boot <- 0
  }else{
    if(isFALSE(btune)){
      tune.work1 <- nu.seq[id.cv5.1];tune.work2 <- nu.seq[id.cv5.2]
    }else{tune.work1 <- tune.work2 <- nu.seq}
    
    allest <- allest.np <- allest.naive <- NULL
    conv.prop.boot <- vector(length = B)
    for (jj in 1:B){
      #bootstrap sample
      set.seed(jj)
      R0.selec.id <- sample(c(1:n),n,T)
      R1.selec.id <- sample(c((n+1):(n+m)),m,T)
      mydt.selec1 <- mydt.1[c(R0.selec.id,R1.selec.id),]; mydt.selec2 <- mydt.2[c(R0.selec.id,R1.selec.id),]
      
      #estimation: proposed
      ic.out1 <- IC.DRest.EM(mydt.selec1,pi0,pi1,
                            nu=tune.work1,k=K,ord=4,ord.pen=2,thres = 1e-5,
                            knots=knots,ini.G=NULL)#no initial#knots are not the same as full data
      id.aic1 <- ic.out1$tune.selection$aic[1]
      ic.out2 <- IC.DRest.EM(mydt.selec2,pi0,pi1,
                            nu=tune.work2,k=K,ord=4,ord.pen=2,thres = 1e-5,
                            knots=knots,ini.G=NULL)#no initial#knots are not the same as full data
      id.aic2 <- ic.out2$tune.selection$aic[1]
      if(is.na(id.aic1*id.aic2)){
        conv.aic <- F
      }else{
        conv.aic <- (ic.out1$tune.res$conv[id.aic1])*(ic.out2$tune.res$conv[id.aic2])
        conv.prop.boot[jj] <- conv.aic
      }
      if(conv.aic==1){
        sumout1 <- sum.DRest0.EM(mydt.selec1,ic.out1$all.list[[id.aic1]],s=s,s.range=s.range)
        sumout2 <- sum.DRest0.EM(mydt.selec2,ic.out2$all.list[[id.aic2]],s=s,s.range=s.range)
        stats.aic <- (sumout1$stats.df-sumout2$stats.df)[-3]
      }else{
        stats.aic <- data.frame(ROC=NA,AUC=NA,Youden=NA,pAUC=NA)
      }
      allest <- rbind(allest,stats.aic)
      
      #estimation: naive
      test.naive1 <- try(DRest0.EM(mydt.selec1,pi0=1,pi1=1,
                                  nu=nu.seq[id.cv5.1],k=K,ord=4,ord.pen=2,
                                  maxit=500,thres=1e-5,knots=knots),silent = T)
      if(class(test.naive1)!="try-error"){
        if(isTRUE(test.naive1$conv)){
          stats.naive1 <- sum.DRest0.EM(mydt.selec1,test.naive1,s=s,s.range=s.range)$stats.df
        }
      }else{
        stats.naive1 <- data.frame(ROC=NA,AUC=NA,cutoff=NA,Youden=NA,pAUC=NA)
      }
      test.naive2 <- try(DRest0.EM(mydt.selec2,pi0=1,pi1=1,
                                   nu=nu.seq[id.cv5.2],k=K,ord=4,ord.pen=2,
                                   maxit=500,thres=1e-5,knots=knots),silent = T)
      if(class(test.naive2)!="try-error"){
        if(isTRUE(test.naive2$conv)){
          stats.naive2 <- sum.DRest0.EM(mydt.selec2,test.naive2,s=s,s.range=s.range)$stats.df
        }
      }else{
        stats.naive2 <- data.frame(ROC=NA,AUC=NA,cutoff=NA,Youden=NA,pAUC=NA)
      }
      stats.naive <- (stats.naive1-stats.naive2)[-3]
      allest.naive <- rbind(allest.naive,stats.naive)
      
      #estimation: np
      if(datatype=="norm"){
        mydt.selec1$biomarker <- exp(mydt.selec1$biomarker)/(1+exp(mydt.selec1$biomarker))
        mydt.selec2$biomarker <- exp(mydt.selec2$biomarker)/(1+exp(mydt.selec2$biomarker))
      }else if(datatype=="gamma"){
        mydt.selec1$biomarker <- (mydt.selec1$biomarker)/(1+(mydt.selec1$biomarker))
        mydt.selec2$biomarker <- (mydt.selec2$biomarker)/(1+(mydt.selec2$biomarker))
      }
      npout1 <- np.roc(mydt.selec1,pi0,pi1,t=s,s.range=s.range)
      npout2 <- np.roc(mydt.selec2,pi0,pi1,t=s,s.range=s.range)
      allest.np <- rbind(allest.np,(npout1-npout2)[-3])
    }
    #summary
    finout$conv.boot <- mean(conv.prop.boot)
    bootsum <- rbind(apply(allest,2,quantile,na.rm=T,probs=c(0.025,0.975,0.05,0.95)),
                     apply(allest,2,mean,na.rm=T),apply(allest,2,sd,na.rm=T))
    bootsum.naive <- rbind(apply(allest.naive,2,quantile,na.rm=T,probs=c(0.025,0.975,0.05,0.95)),
                           apply(allest.naive,2,mean,na.rm=T),apply(allest.naive,2,sd,na.rm=T))
    bootsum.np <- rbind(apply(allest.np,2,quantile,na.rm=T,probs=c(0.025,0.975,0.05,0.95)),
                        apply(allest.np,2,mean,na.rm=T),apply(allest.np,2,sd,na.rm=T))
    rownames(bootsum)[c(5,6)] <- rownames(bootsum.naive)[c(5,6)] <- rownames(bootsum.np)[c(5,6)] <- c("mean","sd")
    finout$bootsum <- bootsum
    finout$bootsum.naive <- bootsum.naive
    finout$bootsum.np <- bootsum.np
    #CI 95%
    finout$cover <- cover.analysis(true.all,finout$paras$cv5,bootsum,debias=F)
    finout$cover.naive <- cover.analysis(true.all,finout$paras$naive,bootsum.naive,debias=F)
    finout$cover.np <- cover.analysis(true.all,finout$paras$np,bootsum.np,debias=F)
  }
  return(finout)
}

sum.pack.bi.boot <- function(res.l){
  RUN <- length(res.l)
  true.mat <- do.call(rbind,lapply(res.l, function(x){x$paras$true}))
  all.true <- apply(true.mat, 2, median)#
  out <- list(para.true=all.true,
              res.df.cv5=do.call(rbind,lapply(res.l, function(x){x$paras$cv5})),
              res.df.np=do.call(rbind,lapply(res.l, function(x){x$paras$np})),
              res.df.naive=do.call(rbind,lapply(res.l, function(x){x$paras$naive})))
  
  id.cv5.1 <- sapply(res.l,function(x){x$tune.id$cv5.1})
  id.cv5.2 <- sapply(res.l,function(x){x$tune.id$cv5.2})
  out$id.cv5.1 <- table(factor(id.cv5.1,levels = 1:length(nu.seq)))/RUN
  out$id.cv5.2 <- table(factor(id.cv5.2,levels = 1:length(nu.seq)))/RUN
  
  conv.T1 <- sapply(res.l, function(x){x$conv$cv5.1})
  conv.T2 <- sapply(res.l, function(x){x$conv$cv5.2})
  out$conv.cv5 <- data.frame(cv5.1=mean(conv.T1),cv5.2=mean(conv.T2))
  
  out$perf <- rbind(perf.12(out$res.df.cv5,all.true),
                    perf.12(out$res.df.np,all.true),
                    perf.12(out$res.df.naive,all.true))
  rownames(out$perf) <- c("cv5","np","naive")

  convprop <- sapply(res.l, function(x){x$conv.boot})
  out$conv.prop.Boot <- data.frame(mean=mean(convprop),median=median(convprop))
  
  cover <- rbind(apply(do.call(rbind,lapply(res.l, function(x){x$cover})),2,mean),
                 apply(do.call(rbind,lapply(res.l, function(x){x$cover.np})),2,mean),
                 apply(do.call(rbind,lapply(res.l, function(x){x$cover.naive})),2,mean))
  rownames(cover) <- c("proposed","np","naive")
  out$cover <- cover
  return(out)
}

