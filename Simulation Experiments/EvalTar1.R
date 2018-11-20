#Libraries
library(tidyverse)

#Set Working Directory for Loading Information
setwd("D:/Mario Documents/Research/JAS/BAYESIAN_SUBSET_TAR/Simulation Experiments")

#Loop to Gather Results
FULL.OUT=list()
for(k in 1:100){
  load(paste("SimTar1_",k,".Rdata",sep=""))
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

#Loop to Get the Change Values of Coefficients
FULL.COEF.CHANGE=NULL
for(k in 1:100){
  FULL.COEF.CHANGE=rbind(FULL.COEF.CHANGE,CUMSUMCOEF.func(k))
}
FULL.COEF.CHANGE$ACTUAL = 0
for(k in 1:20000){
  FULL.COEF.CHANGE$ACTUAL[k] =if(FULL.COEF.CHANGE$Lag[k]=="Intercept"){
                                0
                              }else if(FULL.COEF.CHANGE$Lag[k]=="Lag 1"){
                                0.5
                              }else if(FULL.COEF.CHANGE$Lag[k]=="Lag 2"){
                                0
                              }else if(FULL.COEF.CHANGE$Lag[k]=="Lag 3"){
                                0
                              }
}
  
p1=ggplot(data=FULL.COEF.CHANGE) + 
  geom_point(aes(x=Threshold,y=Coefficient),alpha=0.3,color="gray") + theme_minimal() +
  geom_line(aes(x=Threshold,y=ACTUAL),color="black")+
  facet_wrap(Lag~.) + theme(legend.position="Null")
ggsave("SimTar1Coeff.png",plot=p1,width=6,height=4)
  

#Count Regimes
REGIME.COUNT=rep(NA,100)
for(k in 1:100){
  OUT=FULL.COEF
}














