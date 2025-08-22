###functions for paralell computing 

########################univariate case
pack.uni <- function(run,nu.seq,n,m,pi0,pi1,
                     K=50,knots="q",datatype="norm",
                     s=0.2,s.range=c(0.1,0.3),
                     InputInitG=F,InitGcon=list(degree=1,log1=F,log2=F),
                     alt.auc=F){
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
  ini.G <- NULL
  if(isTRUE(InputInitG)){
    init.tes <- get.init.G(mydt,pi0,pi1,
                           degree=InitGcon$degree,log1=InitGcon$log1,log2=InitGcon$log2,
                           maxit=500,thres=1e-5,ini.G=NULL)
    if(isTRUE(init.tes$conv)){
      ini.G <- init.tes$init.G
    }else{ini.G <- NULL}
  }
  #aic
  ic.out <- IC.DRest.EM(mydt,pi0,pi1,
                        nu=nu.seq,k=K,ord=4,ord.pen=2,thres = 1e-5,
                        knots=knots,ini.G=ini.G)
  id.aic <- ic.out$tune.selection$aic[1]#should check is.na(id.aic)
  conv.aic <- ic.out$tune.res$conv[id.aic]
  if(isTRUE(conv.aic)){
    sumout <- sum.DRest0.EM(mydt,ic.out$all.list[[id.aic]],s=s,s.range=s.range,alt.auc=alt.auc)
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
                                   s=s,s.range=s.range,alt.auc=alt.auc)$stats.df,silent = T)
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
      stats.naive <- sum.DRest0.EM(mydt,test.naive,s=s,s.range=s.range,alt.auc=alt.auc)$stats.df
    }else{
      stats.naive <- data.frame(ROC=NA,AUC=NA,cutoff=NA,Youden=NA,pAUC=NA)
    }
    paras$naive <- stats.naive
  }else{
    paras$naive <- data.frame(ROC=NA,AUC=NA,cutoff=NA,Youden=NA,pAUC=NA)
  }
  
  #np method
  np.d <- np.cutoff.di(mydt,pi0,pi1)
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
  
  #direct youden estimations
  paras.d <- paras;paras.d$aic <- NULL#aic is not included
  paras.d$np$cutoff <- np.d$cutoff;paras.d$np$Youden <- np.d$Youden
  paras.d$cv5$cutoff <- paras.d$naive$cutoff
  paras.d$cv5$Youden <- (paras.d$naive$Youden)/(pi0+pi1-1)
  
  return(list(paras=paras,conv=data.frame(aic=conv.aic,cv5=conv.cv5),
              tune.id=data.frame(aic=id.aic,cv5=id.cv5),
              paras.d=paras.d))
}
sum.pack.uni <- function(res.l){
  RUN <- length(res.l)
  true.mat <- do.call(rbind,lapply(res.l, function(x){x$paras$true}))
  all.true <- apply(true.mat, 2, median)#
  out <- list(para.true=all.true,
              res.df.aic=do.call(rbind,lapply(res.l, function(x){x$paras$aic})),
              res.df.cv5=do.call(rbind,lapply(res.l, function(x){x$paras$cv5})),
              res.df.np=do.call(rbind,lapply(res.l, function(x){x$paras$np})),
              res.df.naive=do.call(rbind,lapply(res.l, function(x){x$paras$naive})),
              res.df.cv5.d=do.call(rbind,lapply(res.l, function(x){x$paras.d$cv5})),
              res.df.np.d=do.call(rbind,lapply(res.l, function(x){x$paras.d$np})))
  id.aic <- sapply(res.l,function(x){x$tune.id$aic})
  id.cv5 <- sapply(res.l,function(x){x$tune.id$cv5})
  out$id.aic <- table(factor(id.aic,levels = 1:length(nu.seq)))/RUN
  out$id.cv5 <- table(factor(id.cv5,levels = 1:length(nu.seq)))/RUN
  
  conv_aic <- sapply(res.l, function(x){x$conv$aic})
  conv_cv5 <- sapply(res.l, function(x){x$conv$cv5})
  out$conv <- data.frame(aic=mean(conv_aic),
                         cv5=mean(conv_cv5))
  fail.run1 <- which(conv_aic==F);fail.run2 <- which(conv_cv5==F)
  out$fail.run <- list(aic=fail.run1,cv5=fail.run2)
  
  #performance:
  out$perf <- rbind(perf.15(out$res.df.aic,all.true),
                    perf.15(out$res.df.cv5,all.true),
                    perf.15(out$res.df.np,all.true),
                    perf.15(out$res.df.naive,all.true))
  rownames(out$perf) <- c("aic","cv5","np","naive")
  
  #comparison to direct method:
  out$comp.d <- rbind(perf.15(out$res.df.cv5,all.true),
                      perf.15(out$res.df.cv5.d,all.true),
                      perf.15(out$res.df.np,all.true),
                      perf.15(out$res.df.np.d,all.true))
  rownames(out$comp.d) <- c("cv5","cv5.d","np","np.d")
  out$comp.d <- out$comp.d[,c(7:12)]
  return(out)
}
perf.15 <- function(res.df,all.true){
  mat.out <- matrix(nrow=1,ncol=15)
  colnames(mat.out) <- c("ROC.b","ROC.sd","ROC.mse",
                         "AUC.b","AUC.sd","AUC.mse",
                         "cutoff.b","cutoff.sd","cutoff.mse",
                         "YI.b","YI.sd","YI.mse",
                         "pAUC.b","pAUC.sd","pAUC.mse")
  mat.out[1,c(1,4,7,10,13)] <- (apply(res.df, 2, mean,na.rm=T)-all.true)
  mat.out[1,c(2,5,8,11,14)] <- apply(res.df, 2, sd,na.rm=T)
  mat.out[1,c(3,6,9,12,15)] <- mat.out[1,c(2,5,8,11,14)]^2+mat.out[1,c(1,4,7,10,13)]^2
  #mat.out[1,c(1,4,7,10,13)] <- mat.out[1,c(1,4,7,10,13)]/all.true
  return(mat.out)
}




########################bivariate case
pack.bi <- function(run,nu.seq,n,m,pi0,pi1,rho,
                    K=50,knots="q",datatype="norm",
                    s=0.2,s.range=c(0.1,0.3),
                    InputInitG=F,InitGcon=list(degree=1,log1=F,log2=F)){
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
  true.all.1 <- c(ROC.true.1,AUC.true.1,YI.true.1,pAUC.true.1)
  true.all.2 <- c(ROC.true.2,AUC.true.2,YI.true.2,pAUC.true.2)
  mydt.1 <- mydt[,1:2];mydt.1$biomarker <- mydt$biomarker1
  mydt.2 <- mydt[,1:2];mydt.2$biomarker <- mydt$biomarker2
  
  paras <- list();paras$true <- true.all
  paras1 <- list();paras1$true <- true.all.1
  paras2 <- list();paras2$true <- true.all.2
  ###estimation (part1)
  ini.G <- NULL
  if(isTRUE(InputInitG)){
    init.tes <- get.init.G(mydt.1,pi0,pi1,
                           degree=InitGcon$degree,log1=InitGcon$log1,log2=InitGcon$log2,
                           maxit=500,thres=1e-5,ini.G=NULL)
    if(isTRUE(init.tes$conv)){
      ini.G <- init.tes$init.G
    }else{ini.G <- NULL}
  }
  cvout <- cv.DRest.EM(mydt.1,pi0,pi1,nfolds=5,seed=123,
                       nu=nu.seq,k=K,ord=4,ord.pen=2,
                       maxit=500,thres=1e-5,knots=knots,ini.G=ini.G)
  id.cv5.1 <- cvout$tune.id
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
  ini.G <- NULL
  if(isTRUE(InputInitG)){
    init.tes <- get.init.G(mydt.2,pi0,pi1,
                           degree=InitGcon$degree,log1=InitGcon$log1,log2=InitGcon$log2,
                           maxit=500,thres=1e-5,ini.G=NULL)
    if(isTRUE(init.tes$conv)){
      ini.G <- init.tes$init.G
    }else{ini.G <- NULL}
  }
  cvout <- cv.DRest.EM(mydt.2,pi0,pi1,nfolds=5,seed=123,
                       nu=nu.seq,k=K,ord=4,ord.pen=2,
                       maxit=500,thres=1e-5,knots=knots,ini.G=ini.G)
  id.cv5.2 <- cvout$tune.id
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
  paras1$cv5 <- stats.cv5.1[-3];paras2$cv5 <- stats.cv5.2[-3]
  
  
  #naive method (using cv5 tuning)
  if(conv.cv5.1*conv.cv5.2==1){
    test.naive <- DRest0.EM(mydt.1,pi0=1,pi1=1,
                            nu=nu.seq[id.cv5.1],k=K,ord=4,ord.pen=2,
                            maxit=500,thres=1e-5,knots=knots)#need to use 'try'
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
    paras1$naive <- stats.naive.1[-3];paras2$naive <- stats.naive.2[-3]
  }else{
    paras$naive <- data.frame(ROC=NA,AUC=NA,Youden=NA,pAUC=NA)
    paras2$naive <- paras1$naive <- data.frame(ROC=NA,AUC=NA,Youden=NA,pAUC=NA)
  }
  
  #np method
  #np.d <- np.cutoff.di(mydt,pi0,pi1)
  if(datatype=="norm"){
    mydt.1$biomarker <- exp(mydt.1$biomarker)/(1+exp(mydt.1$biomarker))
    mydt.2$biomarker <- exp(mydt.2$biomarker)/(1+exp(mydt.2$biomarker))
    #npout$cutoff <- log(npout$cutoff/(1-npout$cutoff))
  }else if(datatype=="gamma"){
    mydt.1$biomarker <- (mydt.1$biomarker)/(1+(mydt.1$biomarker))
    mydt.2$biomarker <- (mydt.2$biomarker)/(1+(mydt.2$biomarker))
    #npout$cutoff <- npout$cutoff/(1-npout$cutoff)
  }
  npout1 <- np.roc(mydt.1,pi0,pi1,t=s,s.range=s.range)
  npout2 <- np.roc(mydt.2,pi0,pi1,t=s,s.range=s.range)
  paras$np <- (npout1-npout2)[-3]
  paras1$np <- npout1[-3];paras2$np <- npout2[-3]
  
  return(list(paras=paras,paras1=paras1,paras2=paras2,
              conv=data.frame(cv5.1=conv.cv5.1,cv5.2=conv.cv5.2),
              tune.id=data.frame(cv5.1=id.cv5.1,cv5.2=id.cv5.2)))
}
sum.pack.bi <- function(res.l){
  RUN <- length(res.l)
  true.mat <- do.call(rbind,lapply(res.l, function(x){x$paras$true}))
  all.true <- apply(true.mat, 2, median)#
  all.true.part1 <- apply(do.call(rbind,lapply(res.l, function(x){x$paras1$true})),
                          2, median)
  all.true.part2 <- apply(do.call(rbind,lapply(res.l, function(x){x$paras2$true})),
                          2, median)
  
  out <- list(para.true=all.true,
              para.true.1=all.true.part1,para.true.2=all.true.part2,
              res.df.cv5=do.call(rbind,lapply(res.l, function(x){x$paras$cv5})),
              res.df.np=do.call(rbind,lapply(res.l, function(x){x$paras$np})),
              res.df.naive=do.call(rbind,lapply(res.l, function(x){x$paras$naive})),
              res.df.cv5.1=do.call(rbind,lapply(res.l, function(x){x$paras1$cv5})),
              res.df.np.1=do.call(rbind,lapply(res.l, function(x){x$paras1$np})),
              res.df.naive.1=do.call(rbind,lapply(res.l, function(x){x$paras1$naive})),
              res.df.cv5.2=do.call(rbind,lapply(res.l, function(x){x$paras2$cv5})),
              res.df.np.2=do.call(rbind,lapply(res.l, function(x){x$paras2$np})),
              res.df.naive.2=do.call(rbind,lapply(res.l, function(x){x$paras2$naive}))
              )
  id.cv5.1 <- sapply(res.l,function(x){x$tune.id$cv5.1})
  id.cv5.2 <- sapply(res.l,function(x){x$tune.id$cv5.2})
  out$id.cv5.1 <- table(factor(id.cv5.1,levels = 1:length(nu.seq)))/RUN
  out$id.cv5.2 <- table(factor(id.cv5.2,levels = 1:length(nu.seq)))/RUN
  
  conv.T1 <- sapply(res.l, function(x){x$conv$cv5.1})
  conv.T2 <- sapply(res.l, function(x){x$conv$cv5.2})
  out$conv <- data.frame(cv5.1=mean(conv.T1),cv5.2=mean(conv.T2))
  fail.run1 <- which(conv.T1==F);fail.run2 <- which(conv.T2==F)
  out$fail.run <- list(T1=fail.run1,T2=fail.run2)
  
  #performance:
  out$perf <- rbind(perf.12(out$res.df.cv5,all.true),
                    perf.12(out$res.df.np,all.true),
                    perf.12(out$res.df.naive,all.true))
  out$perf1 <- rbind(perf.12(out$res.df.cv5.1,all.true.part1),
                          perf.12(out$res.df.np.1,all.true.part1),
                          perf.12(out$res.df.naive.1,all.true.part1))
  out$perf2 <- rbind(perf.12(out$res.df.cv5.2,all.true.part2),
                          perf.12(out$res.df.np.2,all.true.part2),
                          perf.12(out$res.df.naive.2,all.true.part2))
  rownames(out$perf) <- c("cv5","np","naive")
  rownames(out$perf2) <- rownames(out$perf1) <- c("cv5","np","naive")
  return(out)
}
perf.12 <- function(res.df,all.true){
  mat.out <- matrix(nrow=1,ncol=12)
  colnames(mat.out) <- c("ROC.b","ROC.sd","ROC.mse",
                         "AUC.b","AUC.sd","AUC.mse",
                         "YI.b","YI.sd","YI.mse",
                         "pAUC.b","pAUC.sd","pAUC.mse")
  mat.out[1,c(1,4,7,10)] <- (apply(res.df, 2, mean,na.rm=T)-all.true)
  mat.out[1,c(2,5,8,11)] <- apply(res.df, 2, sd,na.rm=T)
  mat.out[1,c(3,6,9,12)] <- mat.out[1,c(2,5,8,11)]^2+mat.out[1,c(1,4,7,10)]^2
  #mat.out[1,c(1,4,7,10)] <- mat.out[1,c(1,4,7,10)]/all.true
  return(mat.out)
}





########################test effect of K (number of B-splines): univariate settings, EM with 5-fold cv
pack.Kstudy <- function(run,nu.seq,n,m,pi0,pi1,
                        Kseq,knots="q",datatype="norm",
                        s=0.2,s.range=c(0.1,0.3)){
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
  paras <- list();paras$true <- true.all
  
  stats.cv5.l <- list();id <- numeric(length = length(Kseq))
  for(jj in 1:length(Kseq)){
    cvout <- cv.DRest.EM(mydt,pi0,pi1,nfolds=5,seed=123,
                         nu=nu.seq,k=Kseq[jj],ord=4,ord.pen=2,
                         maxit=500,thres=1e-5,knots=knots)
    id.cv5 <- cvout$tune.id
    if(length(id.cv5)>0){
      conv.cv5 <- cvout$res$conv[id.cv5]
    }else{
      conv.cv5 <- F
      id.cv5 <- NA
    }
    if(isTRUE(conv.cv5)){
      test <- DRest0.EM(mydt,pi0,pi1,nu=nu.seq[id.cv5],k=Kseq[jj],
                        ord=4,ord.pen=2,maxit=500,thres=1e-5,knots=knots)
      stats.cv5.l[[jj]] <- sum.DRest0.EM(mydt,test,s=s,s.range=s.range)$stats.df
    }else{
      stats.cv5.l[[jj]] <- data.frame(ROC=NA,AUC=NA,cutoff=NA,Youden=NA,pAUC=NA)
    }
    id[jj] <- id.cv5
  }
  paras$stats <- do.call(rbind,stats.cv5.l)
  return(list(paras=paras,id=id,Kseq=Kseq))
}
sum.pack.Kstudy <- function(res.l,Kseq){
  RUN <- length(res.l)
  true.mat <- do.call(rbind,lapply(res.l, function(x){x$paras$true}))
  all.true <- apply(true.mat, 2, median)#
  stats.l <- list();id.l <- list()
  for(j in 1:length(Kseq)){
    stats.l[[j]] <- do.call(rbind,lapply(res.l,function(x){x$paras$stats[j,]}))
    id.j <- sapply(res.l,function(x){x$id[j]})
    id.l[[j]] <- table(factor(id.j,levels = 1:length(nu.seq)))/RUN
  }
  perf <- lapply(stats.l, perf.15,all.true=all.true)
  return(list(id=do.call(rbind,id.l),
              perf=do.call(rbind,perf)))
}




