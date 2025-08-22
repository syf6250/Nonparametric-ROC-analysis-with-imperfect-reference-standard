###methods from Sun J, Tang C, Xie W, Zhou XH. Nonparametric receiver operating characteristic curve analysis with an imperfect gold standard. Biometrics. 2024;80(3):ujae063.

np.roc <- function(dat,pi0=1,pi1=1,t=0.2,s.range=c(0.1,0.3)){
  ###CAUTION: Y should be scaled into [0,1] (see EP below)
  R <- dat$R
  n0<-sum(R==0)
  n1<-sum(R==1)
  Y <- dat$biomarker
  ##AUC
  AUC_<-0
  for(i in 1:(n0+n1)){
    for(j in 1:(n0+n1)){
      if(R[i]==1&R[j]==0){
        AUC_<-AUC_+(Y[i]>Y[j])+1/2*(Y[i]==Y[j])
      }
    }
  }
  AUC_<-AUC_/(n0*n1)
  AUC_est<-(AUC_+(pi1+pi0)/2-1)/(pi1+pi0-1)
  
  ##ROC
  EP<-sort(unique(Y))
  F1<-rep(0,length(EP)+1)
  F0<-rep(0,length(EP)+1)
  for(i in 1:length(EP)){
    G1_t<-sum(Y<=EP[i]&R==1)/sum(R==1)
    G0_t<-sum(Y<=EP[i]&R==0)/sum(R==0)
    F1[i+1]<-pi0/(pi1+pi0-1)*G1_t-(1-pi1)/(pi1+pi0-1)*G0_t
    F0[i+1]<-pi1/(pi1+pi0-1)*G0_t-(1-pi0)/(pi1+pi0-1)*G1_t 
  }
  
  EP_<-c(0,EP,1)
  delta_EP<-rep(0,length(EP)+1)
  for(i in 1:(length(EP)+1)){
    delta_EP[i]<-EP_[i+1]-EP_[i]
  }
  
  F_points<-c(F1,F0)
  F_points<-F_points*(F_points>=0)*(F_points<=1)
  F_points<-sort(unique(F_points))
  
  Q1<-rep(0,length(F_points))
  Q0<-rep(0,length(F_points))
  for(i in 1:length(F_points)){
    Q1[i]<-sum((F1<=F_points[i])*delta_EP)
    Q0[i]<-sum((F0<=F_points[i])*delta_EP)
  }
  
  q0<-sum((F0<=(1-t))*delta_EP)
  q1<-1-F_points[max(c(which(Q1<=q0),1))]
  
  ###newly added functions (not included in the original paper from Sun et al.):
  t.eval <- sort(unique(c(Q0,Q1)))
  delta.F <- numeric(length(t.eval))
  for(i in 1:length(t.eval)){
    delta.F[i] <- F_points[max(c(which(Q0<=t.eval[i]),1))]-F_points[max(c(which(Q1<=t.eval[i]),1))]
  }
  # #pAUC:
  # #upper and lower limits
  # lo <- sum((F0<=(1-s.range[2]))*delta_EP)
  # up <- sum((F0<=(1-s.range[1]))*delta_EP)
  # 
  # pAUC_<-0
  # for(i in 1:(n0+n1)){
  #   for(j in 1:(n0+n1)){
  #     if(R[i]==1&R[j]==0&Y[j]>=lo&Y[j]<=up){
  #       pAUC_<-pAUC_+(Y[i]>Y[j])+1/2*(Y[i]==Y[j])
  #     }
  #   }
  # }
  # pAUC_<-pAUC_/(n0*n1)/(s.range[2]-s.range[1])
  # 
  # ROC.s1<-1-F_points[max(c(which(Q1<=up),1))]
  # ROC.s2<-1-F_points[max(c(which(Q1<=lo),1))]
  # modify1 <- (1-pi0)*pi1/2*(ROC.s2^2-ROC.s1^2)/(s.range[2]-s.range[1])
  # modify2 <- (1-pi0)*(1-pi1)*(s.range[2]*ROC.s2-s.range[1]*ROC.s1)/(s.range[2]-s.range[1])
  # modify0 <- (1-pi1)*pi0/2*(s.range[2]+s.range[1])
  # pAUC_est <- (pAUC_-modify0-modify1-modify2)/(pi1+pi0-1)
  
  #pAUC: direct integration:
  ROC.fun <- function(s){
    tem<-sum((F0<=(1-s))*delta_EP)
    return(1-F_points[max(c(which(Q1<=tem),1))])
  }
  pAUC_est2 <- try((integrate(Vectorize(ROC.fun),
                              s.range[1],s.range[2])$value)/(s.range[2]-s.range[1]),
                   silent = T)
  if(class(pAUC_est2)=="try-error"){
    x.eval <- seq(from=s.range[1],to=s.range[2],length.out=1000)
    pAUC_est2 <- trapz(x.eval,Vectorize(ROC.fun)(x.eval))/(s.range[2]-s.range[1])
  }
  return(data.frame(ROC=q1,AUC=AUC_est,
                    cutoff=t.eval[which.max(delta.F)],
                    Youden=max(delta.F),pAUC2=pAUC_est2#,pAUC=pAUC_est
                    ))
}

###direct method to estimate the cutoff and Youden:
np.cutoff.di <- function(dat,pi0=1,pi1=1){
  id0 <- which(dat$R==0);id1 <- which(dat$R==1)
  require(spatstat.univar)
  F0 <- ewcdf(dat$biomarker[id0])
  F1 <- ewcdf(dat$biomarker[id1])
  
  eval <- sort(unique(dat$biomarker))
  delta.F <- F0(eval)-F1(eval)
  return(data.frame(cutoff=eval[which.max(delta.F)],
                    Youden=max(delta.F)/(pi0+pi1-1)))
}



