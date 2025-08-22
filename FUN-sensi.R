###functions for sensitivity analysis

gnormmix.assump1 <- function(n0,n1,rho=0,varX1=0.5,varX2=0.5,vareps=0.5,prev=0.4,
                             coef=list(bX1=1,bGT=1,bX2=1,bGR=4.78,aR=-2.39)){
  countR0 <- countR1 <- 0
  dtR0 <- dtR1 <- NULL
  while (countR0<n0 | countR1<n1) {
    Gi <- rbinom(1,1,prev)#true label
    Xi <- mvtnorm::rmvnorm(1,mean=c(0,0),sigma=matrix(c(varX1,rho,rho,varX2),2,2))#X
    Ti <- (coef$bX1)*Xi[1]+(coef$bGT)*Gi+rnorm(1,mean=0,sd=sqrt(vareps))#T
    prb <- 1/(1+exp((-coef$bX2)*Xi[2]-(coef$bGR)*Gi-coef$aR))
    Ri <- rbinom(1,1,prb)#reference label
    if(Ri==0){
      dtR0 <- rbind(dtR0,data.frame(R=Ri,biomarker=Ti,G=Gi))
      countR0 <- countR0+1
    }else{
      dtR1 <- rbind(dtR1,data.frame(R=Ri,biomarker=Ti,G=Gi))
      countR1 <- countR1+1
    }
  }
  return(rbind(dtR0[c(1:n0),],dtR1[c(1:n1),]))
}



pack.uni.asmp1 <- function(run,nu.seq,n,m,pi0,pi1,
                           K=50,knots="q",s=0.2,s.range=c(0.1,0.3),
                           rho=0,varX1=0.5,varX2=0.5,vareps=0.5,prev=0.4,
                           coef=list(bX1=1,bGT=1,bX2=1,bGR=4.78,aR=-2.39)){
  ###settings
  mu <- c(0,coef$bGT)#mean of R=0 & 1; the second element should be larger!
  sigma <- c(varX1+vareps,varX1+vareps)#sd of R=0 & 1
  ROC.true <- 1-pnorm(qnorm(1-s,mu[1],sigma[1]),mu[2],sigma[2])
  AUC.true <- integrate(function(c){(1-pnorm(c,mu[2],sigma[2]))*dnorm(c,mu[1],sigma[1])},
                        lower=-50,upper=50)$value
  true.cutoff <- rootSolve::uniroot.all(function(t){log(dnorm(t,mu[2],sigma[2])/dnorm(t,mu[1],sigma[1]))},c(-5,5))
  YI.true <- pnorm(true.cutoff,mu[1],sigma[1])-pnorm(true.cutoff,mu[2],sigma[2])
  ROC.fun0 <- function(s){1-pnorm(qnorm(1-s,mu[1],1),mu[2],1)}
  pAUC.true <- (integrate(ROC.fun0,s.range[1],s.range[2])$value)/(s.range[2]-s.range[1])
  true.all <- c(ROC.true,AUC.true,true.cutoff,YI.true,pAUC.true)
  ###data
  set.seed(run*77+7)
  mydt <- gnormmix.assump1(n,m,rho=rho,varX1=varX1,varX2=varX2,vareps=vareps,prev=prev,coef=coef)
  
  ###estimation
  ini.G <- NULL
  
  #aic
  ic.out <- IC.DRest.EM(mydt,pi0,pi1,
                        nu=nu.seq,k=K,ord=4,ord.pen=2,thres = 1e-5,
                        knots=knots,ini.G=ini.G)
  id.aic <- ic.out$tune.selection$aic[1]
  if(is.na(id.aic)){
    conv.aic <- F
  }else{
    conv.aic <- ic.out$tune.res$conv[id.aic]
  }
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
  mydt$biomarker <- exp(mydt$biomarker)/(1+exp(mydt$biomarker))
  npout <- np.roc(mydt,pi0,pi1,t=s,s.range=s.range)
  npout$cutoff <- log(npout$cutoff/(1-npout$cutoff))
  paras$np <- npout
  
  return(list(paras=paras,conv=data.frame(aic=conv.aic,cv5=conv.cv5),
              tune.id=data.frame(aic=id.aic,cv5=id.cv5)))
}



sum.pack.uni.assump1 <- function(res.l){
  RUN <- length(res.l)
  true.mat <- do.call(rbind,lapply(res.l, function(x){x$paras$true}))
  all.true <- apply(true.mat, 2, median)
  out <- list(para.true=all.true,
              res.df.aic=do.call(rbind,lapply(res.l, function(x){x$paras$aic})),
              res.df.cv5=do.call(rbind,lapply(res.l, function(x){x$paras$cv5})),
              res.df.np=do.call(rbind,lapply(res.l, function(x){x$paras$np})),
              res.df.naive=do.call(rbind,lapply(res.l, function(x){x$paras$naive})))
  id.aic <- sapply(res.l,function(x){x$tune.id$aic})
  id.cv5 <- sapply(res.l,function(x){x$tune.id$cv5})
  out$id.aic <- table(factor(id.aic,levels = 1:length(nu.seq)))/RUN
  out$id.cv5 <- table(factor(id.cv5,levels = 1:length(nu.seq)))/RUN
  
  conv_aic <- sapply(res.l, function(x){x$conv$aic})
  conv_cv5 <- sapply(res.l, function(x){x$conv$cv5})
  out$conv <- data.frame(aic=mean(conv_aic),cv5=mean(conv_cv5))
  fail.run1 <- which(conv_aic==F);fail.run2 <- which(conv_cv5==F)
  out$fail.run <- list(aic=fail.run1,cv5=fail.run2)
  
  #performance:
  out$perf <- rbind(perf.15(out$res.df.aic,all.true),
                    perf.15(out$res.df.cv5,all.true),
                    perf.15(out$res.df.np,all.true),
                    perf.15(out$res.df.naive,all.true))
  rownames(out$perf) <- c("aic","cv5","np","naive")
  return(out)
}


pack.uni.asmp2 <- function(run,nu.seq,n,m,pi0,pi1,
                           K=50,knots="q",datatype="norm",
                           s=0.2,s.range=c(0.1,0.3),
                           pi0.true=0.9,pi1.true=0.9){
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
    mydt <- gnormmix(n,m,pi0.true,pi1.true,mu=mu,sigma=sigma)
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
    mydt <- ggammamix(n,m,pi0.true,pi1.true,shape=shape,rate=rate)
  } 
  
  ###estimation
  ini.G <- NULL
  #aic
  ic.out <- IC.DRest.EM(mydt,pi0,pi1,
                        nu=nu.seq,k=K,ord=4,ord.pen=2,thres = 1e-5,
                        knots=knots,ini.G=ini.G)
  id.aic <- ic.out$tune.selection$aic[1]#
  if(is.na(id.aic)){
    conv.aic <- F
  }else{
    conv.aic <- ic.out$tune.res$conv[id.aic]
  }
  
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
  id.cv5 <- cvout$tune.id#
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
  if(datatype=="norm"){
    mydt$biomarker <- exp(mydt$biomarker)/(1+exp(mydt$biomarker))
    npout <- np.roc(mydt,pi0,pi1,t=s,s.range=s.range)
    npout$cutoff <- log(npout$cutoff/(1-npout$cutoff))
  }else if(datatype=="gamma"){
    mydt$biomarker <- (mydt$biomarker)/(1+(mydt$biomarker))
    npout <- np.roc(mydt,pi0,pi1,t=s,s.range=s.range)
    npout$cutoff <- npout$cutoff/(1-npout$cutoff)
  }
  paras$np <- npout
  
  return(list(paras=paras,conv=data.frame(aic=conv.aic,cv5=conv.cv5),
              tune.id=data.frame(aic=id.aic,cv5=id.cv5)))
}