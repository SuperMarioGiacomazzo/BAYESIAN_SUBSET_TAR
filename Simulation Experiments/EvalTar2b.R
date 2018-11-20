#Libraries
library(tidyverse)

#Set Working Directory for Loading Information
setwd("D:/Mario Documents/Research/JAS/BAYESIAN_SUBSET_TAR/Simulation Experiments")

#Loop to Gather Results
FULL.OUT=list()
for(k in 1:100){
  load(paste("SimTar2_",k,".Rdata",sep=""))
  FULL.OUT[[k]]=OUT3
}

rm("OUT2","OUT3")

#Set Working Director for Saving Output
setwd("D:/Mario Documents/Graduate School/Ph.D Dissertation/TAR Simulation")

#Functions
CUMSUMCOEF.func=function(k){
  OUT=FULL.OUT[[k]]$mean.coef
  COEF.MAT=matrix(OUT,nrow=51,ncol=4,byrow=T)
  COEF.MAT2=apply(COEF.MAT,2,cumsum)
  THRESH=FULL.OUT[[k]]$seq.thresh
  FINAL=as.tibble(cbind(k,THRESH,COEF.MAT2[-1,]))
  FINAL$Regime=1:50
  names(FINAL)=c("Sim","Threshold","Intercept","Lag 1","Lag 2","Lag 3","Regime")
  FINAL2=FINAL %>% gather(3:6,key="Lag",value="Coefficient")
  FINAL2$Lag=factor(FINAL2$Lag,levels=c("Intercept","Lag 1","Lag 2","Lag 3"))
  return(FINAL2)
}

#Loop to Get the Change Values of 
FULL.COEF.CHANGE=NULL
for(k in 1:100){
  FULL.COEF.CHANGE=rbind(FULL.COEF.CHANGE,CUMSUMCOEF.func(k))
}

#Create Vector of True Values for Figure
FULL.COEF.CHANGE$ACTUAL = 0
for(k in 1:20000){
  FULL.COEF.CHANGE$ACTUAL[k] =if(FULL.COEF.CHANGE$Lag[k]=="Intercept"){
    0
  }else if(FULL.COEF.CHANGE$Lag[k]=="Lag 1" & FULL.COEF.CHANGE$Threshold[k] < -0.4){
      0.5
  }else if(FULL.COEF.CHANGE$Lag[k]=="Lag 1" & FULL.COEF.CHANGE$Threshold[k] > -0.4){
      -0.5
  }else if(FULL.COEF.CHANGE$Lag[k]=="Lag 2"){
    0
  }else if(FULL.COEF.CHANGE$Lag[k]=="Lag 3"){
    0
  }
}

#Plot Showing the Coefficients At All Thresholds Considered
p1=ggplot(data=FULL.COEF.CHANGE) + 
  geom_point(aes(x=Threshold,y=Coefficient),alpha=0.3,color="gray") + theme_minimal() +
  geom_line(aes(x=Threshold,y=ACTUAL),color="black")+
  facet_wrap(Lag~.) + theme(legend.position="Null")
ggsave("SimTar2bCoeff.png",plot=p1,width=6,height=4)

#Count Regimes
REGIME.COUNT=rep(NA,100)
for(k in 1:100){
  OUT=FULL.OUT[[k]]$mean.coef
  COEF.MAT=matrix(OUT,nrow=51,ncol=4,byrow=T)
  COEF.VEC=apply(COEF.MAT,1,function(x){sum(x==0)})
  REGIME.COUNT[k]=sum(COEF.VEC!=4)
}

PROP.TRUE=sum(REGIME.COUNT==2)/100
PROP.OVER=sum(REGIME.COUNT>2)/100
PROP.UNDER=sum(REGIME.COUNT<2)/100

#Estimates of Significant Coefficients of Correctly Identified Models
SIG.PARAM.OUT=matrix(NA,100,4)
Accurate.Thresh=which(REGIME.COUNT==2)
for(k in Accurate.Thresh){
  OUT1A=FULL.OUT[[k]]$mean.coef
  OUT1B=matrix(OUT1A,nrow=51,ncol=4,byrow=T)
  idnonzero=which(apply(OUT1B,1,sum)!=0)
  OUT1C=OUT1B[idnonzero,]
  OUT1D=OUT1C[1,2]
  OUT2A=OUT1C[2,2]+OUT1D
  OUT3A=sqrt(FULL.OUT[[k]]$mean.s2)
  OUT4A=FULL.OUT[[k]]$seq.thresh[idnonzero[2]-1]
  SIG.PARAM.OUT[k,]=c(OUT1D,OUT2A,OUT3A,OUT4A)
}

SIG.PARAM.OUT2=SIG.PARAM.OUT[which(REGIME.COUNT==2),]

bias.func=function(x){
  x-c(0.5,-0.5,2,-0.4)
}
sqerror.func=function(x){
  (x-c(0.5,-0.5,2,-0.4))^2
}

BIAS.OUT=t(apply(SIG.PARAM.OUT2,1,bias.func))
MSE.OUT=t(apply(SIG.PARAM.OUT2,1,sqerror.func))

#Get Table For Paper
OUT.TABLE=cbind(colMeans(SIG.PARAM.OUT2),
                apply(SIG.PARAM.OUT2,2,sd),
                colMeans(BIAS.OUT),
                sqrt(colMeans(MSE.OUT)))
