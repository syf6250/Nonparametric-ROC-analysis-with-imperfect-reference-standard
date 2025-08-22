library(mgcv)
library(ggplot2)
library(spatstat.univar)
library(rootSolve)
source('FUN-simu.R')
source('FUN-np.R')
source('FUN.R')
rawdata <- read.table("data.txt")
colnames(rawdata) <- c("R","biomarker")
mydt <- rawdata[rawdata$biomarker>0,]
mydt$biomarker <- log(mydt$biomarker)

#id2rm <- which(mydt$biomarker==min(mydt$biomarker))
#mydt <- mydt[-id2rm,]


### n:R=0; m:R=1
n <- sum(mydt$R==0);R0id <- which(mydt$R==0)
m <- sum(mydt$R==1);R1id <- which(mydt$R==1)
### pi0=P(G=0|R=0); pi1=P(G=1|R=1)
pi0 <- 1;pi1 <- 0.677



###tuning selection:
#AIC
nu.seq <- 10^seq(from=-3,to=2,length.out=11)
ic.out <- IC.DRest.EM(mydt,pi0,pi1,
                      nu=nu.seq,k=50,ord=4,ord.pen=2,maxit=1000,
                      thres = 1e-5,knots="q")#check edf to see if the range of tunings is appropriate
ic.out$tune.res


#CV
cvout <- cv.DRest.EM(mydt,pi0,pi1,nfolds=5,seed=123,
                     nu=nu.seq,k=50,ord=4,ord.pen=2,
                     maxit=1000,thres=1e-5,knots="q")
#cvout$tune
cvout$tune.id#the number of the picked tuning (see also application-cv5.R for cross validation results)
#cvout$res

#stats 
s <- 0.2;s.range=c(0.1,0.3)
stats <- sum.DRest0.EM(mydt,ic.out$all.list[[cvout$tune.id]])
stats <- stats$stats.df
stats


###fully nonparametric method:
mydt2 <- mydt
a <- min(mydt2$biomarker);b <- max(mydt2$biomarker)
mydt2$biomarker <- (mydt2$biomarker-a)/(b-a)
npout <- np.roc(mydt2,pi0,pi1,s,s.range)
npout$cutoff <- npout$cutoff*(b-a)+a
npout

###naive method:
test0 <- DRest0.EM(mydt,pi0=1,pi1=1,
                   nu=nu.seq[cvout$tune.id],k=50,ord=4,ord.pen=2,
                   maxit=1000,thres=1e-5,
                   knots="q",ini.G=NULL)
stats0 <- sum.DRest0.EM(mydt,test0)
stats0 <- stats0$stats.df
stats0


###plots: ROC
pts <- seq(0,1,0.01)
ROC <- est.ROC(pts,mydt,ic.out$all.list[[cvout$tune.id]])
s <- 0.2;s.range=c(0.1,0.3)
mydt2 <- mydt
a <- min(mydt2$biomarker);b <- max(mydt2$biomarker)
mydt2$biomarker <- (mydt2$biomarker-a)/(b-a)
ROCnp <- numeric(length(pts))
for(i in 1:length(ROCnp)){
  npout <- np.roc(mydt2,pi0,pi1,pts[i],s.range)
  ROCnp[i] <- npout$ROC
}
ROCnaive <- est.ROC(pts,mydt,test0)
df2draw <- data.frame(s=pts,EM=ROC,NP=ROCnp,naive=ROCnaive)
library(reshape2)
df2draw.melt <- melt(df2draw,id.vars="s")
colnames(df2draw.melt)[2] <- "method"
p.roc <- ggplot(data=df2draw.melt,mapping=aes(x=s,y=value,group=method))+
  geom_line(aes(linetype=method,color=method))+
  scale_y_continuous(name="Sensitivity")+scale_x_continuous(name="1 - specificity")+
  theme_bw()+theme(legend.position='none',
                   axis.title.x = element_text(size=15),axis.title.y = element_text(size=15))
#ggsave("roc.pdf",p.roc,device = "pdf")





###draw distributions for R=0,1: 
#R=0
mydt_R0 <- subset(mydt,R==0)
emp0 <- ewcdf(mydt_R0$biomarker)
EM0 <- ewcdf(mydt$biomarker,weights = ic.out$all.list[[cvout$tune.id]]$weight.Dis$G0)
#R=1
mydt_R1 <- subset(mydt,R==1)
emp1 <- ewcdf(mydt_R1$biomarker)
EM1 <- ewcdf(mydt$biomarker,
             weights = ic.out$all.list[[cvout$tune.id]]$weight.Dis$G1)
EM.mix <- function(t){pi1*EM1(t)+(1-pi1)*EM0(t)}
id2del <- which(mydt_R1$biomarker==0)
mydt_R1 <- mydt_R1[-id2del,]
Tall <- c(mydt_R0$biomarker,mydt_R1$biomarker)

cdfall <- ggplot()+scale_x_continuous(name="log positive parasite",limits = c(min(Tall),max(Tall)))+
  theme_bw()+scale_y_continuous(name="Distribution function value")+
  geom_function(fun=EM0,aes(colour = "proposed0",linetype = "proposed0"))+
  geom_function(fun=emp0,aes(colour = "empirical0",linetype = "empirical0"))+
  geom_function(fun=EM.mix,aes(colour = "proposed1",linetype = "proposed1"))+
  geom_function(fun=emp1,aes(colour = "empirical1",linetype = "empirical1"))+
  scale_colour_manual(name="method",values =c('blue','#00BA38','#F8766D','black'), 
                      labels = c('empirical0','empirical1','proposed0','proposed1'))+
  scale_linetype_manual(name="method",values = c(3,4,1,2), 
                        labels = c('empirical0','empirical1','proposed0','proposed1'))+
  theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),
        legend.position = 'none')

roccdf <- gridExtra::grid.arrange(p.roc,cdfall,ncol=2)
#ggsave("roccdf.pdf",roccdf,device = "pdf")






