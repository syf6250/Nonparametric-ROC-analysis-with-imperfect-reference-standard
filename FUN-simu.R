###simulated data generation


###generate pi0 pi1 from sensitivity, specificity and prevalence of R 
computePi <- function(Sen,Spe,prev){
  ppv<-prev*Sen/(prev*Sen+(1-prev)*(1-Spe))
  npv<-(1-prev)*Spe/((1-prev)*Spe+prev*(1-Sen))
  return(data.frame(pi0=npv,pi1=ppv))
}


#########single biomarker
###generate normal mixture data (one dimensional biomarker)
gnormmix <- function(n0,n1,pi0=1,pi1=1,mu=c(0,1),sigma=c(1,1),trans=F){
  ##################
  ###INPUTS:
  #n0,n1: No. of subjects for R=0 and R=1 groups
  #pi0=P(G=0|R=0); pi1=P(G=1|R=1)
  #mu=(mean of G=0,mean of G=1)
  #sigma=(sd of G=0,sd of G=1)
  #trans: if True, x is modified by exp(x)/(1+exp(x))!
  ###OUTPUTS: a dataframe containing T (biomarker values), R & G.
  ##################
  df <- data.frame(R=c(rep(0,n0),rep(1,n1)))
  ###generation of R=0 group
  z0 <- sample(2,n0, replace = TRUE, prob = c(pi0,1-pi0))#z=1:G=0;z=2:G=1
  T0 <- rnorm(n0,mean=mu[z0],sd=sigma[z0])
  ###generation of R=1 group
  z1 <- sample(2,n1, replace = TRUE, prob = c(1-pi1,pi1))#z=1:G=0;z=2:G=1
  T1 <- rnorm(n1,mean=mu[z1],sd=sigma[z1])
  
  df$biomarker <- c(T0,T1)
  if(isTRUE(trans)){
    df$biomarker <- exp(df$biomarker)/(1+exp(df$biomarker))
  }
  df$G <- c(z0-1,z1-1)
  return(df)
}

#determine interval/points to evaluation h:
gen.eval.pts <- function(data.type='norm',
                         n0,n1,pi0=1,pi1=1,mu=c(0,1),sigma=c(1,1),
                         shape=c(2,3),rate=c(1,1),
                         Repeat=100,p=c(0.05,0.95),no.pts=30){
  Tdata <- list()
  for(j in 1:Repeat){
    if(data.type=='norm'){
      Tdata[[j]] <- gnormmix(n0,n1,pi0,pi1, mu=mu,sigma=sigma)$biomarker
    }else if(data.type=='gamma'){
      Tdata[[j]] <- ggammamix(n0,n1,pi0,pi1, shape=shape,rate=rate)$biomarker
    }
    #else?
  }
  Tdata <- unlist(Tdata)
  eval.interval <- quantile(Tdata,probs=p)
  eval.pts <- seq(from=eval.interval[1],to=eval.interval[2],length.out=no.pts)
  return(list(eval.interval=eval.interval,eval.pts=eval.pts))
}

###generate gamma mixture data (one dimensional biomarker)
ggammamix <- function(n0,n1,pi0=1,pi1=1,shape=c(2,3),rate=c(1,1),trans=F){
  ##################
  ###INPUTS:
  #n0,n1: No. of subjects for R=0 and R=1 groups
  #pi0=P(G=0|R=0); pi1=P(G=1|R=1)
  #shape=(shape of G=0,mean of G=1)
  #rate=(rate of G=0,sd of G=1)
  #trans: if True, x is modified by exp(x)/(1+exp(x))!
  ###OUTPUTS: a dataframe containing T (biomarker values), R & G.
  ##################
  df <- data.frame(R=c(rep(0,n0),rep(1,n1)))
  ###generation of R=0 group
  z0 <- sample(2,n0, replace = TRUE, prob = c(pi0,1-pi0))#z=1:G=0;z=2:G=1
  T0 <- rgamma(n0,shape=shape[z0],rate=rate[z0])
  ###generation of R=1 group
  z1 <- sample(2,n1, replace = TRUE, prob = c(1-pi1,pi1))#z=1:G=0;z=2:G=1
  T1 <- rgamma(n1,shape=shape[z1],rate=rate[z1])
  
  df$biomarker <- c(T0,T1)
  if(isTRUE(trans)){
    df$biomarker <- exp(df$biomarker)/(1+exp(df$biomarker))
  }
  df$G <- c(z0-1,z1-1)
  return(df)
}




#########multiple biomarkers
gnormmix2 <- function(n0,n1,pi0=1,pi1=1,
                      mu1=c(1.5,1),mu0=c(0,0),
                      var1=matrix(c(1,0.2,0.2,1),2,2),
                      var0=matrix(c(1,0.2,0.2,1),2,2)){
  ##################
  ###INPUTS:
  #n0,n1: No. of subjects for R=0 and R=1 groups
  #pi0=P(G=0|R=0); pi1=P(G=1|R=1)
  #mu0:mean of G=0; mu1:mean of G=1
  #var0:variance of G=0; var1:variance of G=1
  ###OUTPUTS: a dataframe containing T (biomarker values), R & G.
  ##################
  df <- data.frame(R=c(rep(0,n0),rep(1,n1)))
  ###generation of R=0 group
  z0 <- sample(2,n0, replace = TRUE, prob = c(pi0,1-pi0))#z=1:G=0;z=2:G=1
  ###generation of R=1 group
  z1 <- sample(2,n1, replace = TRUE, prob = c(1-pi1,pi1))#z=1:G=0;z=2:G=1
  df$G <- c(z0-1,z1-1)
  G.ind0 <- which(df$G==0);G.ind1 <- which(df$G==1)
  
  ###biomarkers
  df$biomarker2 <- df$biomarker1 <- 0
  df[G.ind0,3:4] <- mvtnorm::rmvnorm(length(G.ind0),mean=mu0,sigma=var0)
  df[G.ind1,3:4] <- mvtnorm::rmvnorm(length(G.ind1),mean=mu1,sigma=var1)
  return(df)
}

ggammamix2 <- function(n0,n1,pi0=1,pi1=1,
                       shape1=c(4,3),shape0=c(2,2),
                       rho=0.5){
  ##################
  ###INPUTS:
  #n0,n1: No. of subjects for R=0 and R=1 groups
  #pi0=P(G=0|R=0); pi1=P(G=1|R=1)
  #shape1: shape of G=1 (for two biomarkers)
  #rho: cov(T_1,T_2), given rate==1 (should be smaller than shapes)
  ###OUTPUTS: a dataframe containing T (biomarker values), R & G.
  ##################
  df <- data.frame(R=c(rep(0,n0),rep(1,n1)))
  ###generation of R=0 group
  z0 <- sample(2,n0, replace = TRUE, prob = c(pi0,1-pi0))#z=1:G=0;z=2:G=1
  ###generation of R=1 group
  z1 <- sample(2,n1, replace = TRUE, prob = c(1-pi1,pi1))#z=1:G=0;z=2:G=1
  df$G <- c(z0-1,z1-1)
  G.ind0 <- which(df$G==0);G.ind1 <- which(df$G==1)
  
  ###biomarkers
  df$biomarker2 <- df$biomarker1 <- rgamma(n0+n1,shape=rho,rate=1)
  df$biomarker1[G.ind0] <- df$biomarker1[G.ind0]+rgamma(length(G.ind0),shape=shape0[1]-rho,rate=1)
  df$biomarker1[G.ind1] <- df$biomarker1[G.ind1]+rgamma(length(G.ind1),shape=shape1[1]-rho,rate=1)
  df$biomarker2[G.ind0] <- df$biomarker2[G.ind0]+rgamma(length(G.ind0),shape=shape0[2]-rho,rate=1)
  df$biomarker2[G.ind1] <- df$biomarker2[G.ind1]+rgamma(length(G.ind1),shape=shape1[2]-rho,rate=1)
  return(df)
}


