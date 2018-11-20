#Libraries Used
library(doParallel)
library(forecast)
library(foreach)
library(bayesreg)
library(coda)
library(MASS)
library(ggplot2)
library(fanplot)
library(stargazer)
library(xtable)

options(scipen=999)

#Working Directory
setwd("D:/Mario Documents/Research/JAS/BAYESIAN_SUBSET_TAR/Traffic Analysis/Horizon 5")

# Function To Output Sequence of Thresholds
seq.thresh.func=function(day,series,P,h,nquant){
  #Obtain Train Data for Detector 103 including Before (L108) and After (L101)
  L106.TR=April.3.Day[[day]]$L106_occupancy[1:1440]/100
  L101.TR=April.3.Day[[day]]$L101_occupancy[1:1440]/100
  L108.TR=April.3.Day[[day]]$L108_occupancy[1:1440]/100
  L104.TR=April.3.Day[[day]]$L104_occupancy[1:1440]/100
  L102.TR=April.3.Day[[day]]$L102_occupancy[1:1440]/100
  L107.TR=April.3.Day[[day]]$L107_occupancy[1:1440]/100
  L103.TR=April.3.Day[[day]]$L103_occupancy[1:1440]/100
  
  #Adjust Data for 0 and 1 to be 0.0001 and 0.9999
  L106.TR[L106.TR==1]=0.9999
  L106.TR[L106.TR==0]=0.0001
  L101.TR[L101.TR==1]=0.9999
  L101.TR[L101.TR==0]=0.0001
  L108.TR[L108.TR==1]=0.9999
  L108.TR[L108.TR==0]=0.0001
  L104.TR[L104.TR==1]=0.9999
  L104.TR[L104.TR==0]=0.0001
  L102.TR[L102.TR==1]=0.9999
  L102.TR[L102.TR==0]=0.0001
  L107.TR[L107.TR==1]=0.9999
  L107.TR[L107.TR==0]=0.0001
  L103.TR[L103.TR==1]=0.9999
  L103.TR[L103.TR==0]=0.0001
  
  #Logit Transformed Data
  yA=logit.func(L101.TR)
  yB=logit.func(L106.TR)
  yC=logit.func(L108.TR)
  yD=logit.func(L102.TR)
  yE=logit.func(L104.TR)
  yF=logit.func(L107.TR)
  yG=logit.func(L103.TR)
  
  #List of Time Series
  y=list(yA,yB,yC,yD,yE,yF,yG)[[series]]
  
  #Creation of Matrix
  N=length(y)
  X=matrix(NA,N,P)
  for (j in 1:P){
    X[,j]=lag.func(y,k=j+h-1)
  }
  
  #Remove Missing Data
  X=X[-(1:(P+h-1)),]
  y=y[-(1:(P+h-1))]
  N=length(y)
  
  #Reorder X and y matrix According to Delay Parameter
  X.order=X[order(X[,1]),]
  y.order=y[order(X[,1])]
  t.order=(1:N)[order(X[,1])]
  
  #Determine Possible Thresholds To Consider
  seq.quant=seq(0.15,0.85,length=nquant)
  seq.thresh=revlogit.func(as.numeric(quantile(X[,1],seq.quant)))

  #Create Variable Groups Based on Quantiles
  Var.Groups=rep(0:nquant,each=(P+1))

  #Out List
  out=list(seq.quant=seq.quant,seq.thresh=seq.thresh,Var.Groups=Var.Groups)
  
  return(out)
}

#Load Data From Results
load("TARTRAFFIC5C.Rdata")

FINAL.MOD.NREGIME=matrix(NA,length(DAY.ID),length(LOOP.ID))

for(day in DAY.ID){
  for(series in LOOP.ID){
    OUT3=TAR.RESULTS3[[day]][[series]]
    THRESH=seq.thresh.func(day=day,series=series,P=7,h=5,nquant=50)
    REG.IN=rep(NA,51)
    for(k in unique(THRESH$Var.Groups)){
      ID=which(THRESH$Var.Groups==k)
      COEF=OUT3$mean.coef[ID]
      REG.IN[k+1]=sum(COEF==0)!=8
    }
    FINAL.MOD.NREGIME[day,which(LOOP.ID==series)]=sum(REG.IN)
  }
}  

rownames(FINAL.MOD.NREGIME)=DAYS
colnames(FINAL.MOD.NREGIME)=LOOPS
write.csv(FINAL.MOD.NREGIME,file="FINAL.MOD.NREGIME5.csv")
xtable(FINAL.MOD.NREGIME)

DELTA=NULL
count=0
for(day in DAY.ID){
  for(series in LOOP.ID){
    count=count+1
    OUT3=TAR.RESULTS3[[day]][[series]]
    THRESH=seq.thresh.func(day=day,series=series,P=7,h=5,nquant=50)
    REG.IN=rep(NA,51)
    for(k in unique(THRESH$Var.Groups)){
      ID=which(THRESH$Var.Groups==k)
      COEF=OUT3$mean.coef[ID]
      REG.IN[k+1]=sum(COEF==0)!=8
    }
    REG.IN2=REG.IN[-1]
    INIT.DELTA=c(0,rep(NA,(8-1)))
    if(sum(REG.IN2)>0) INIT.DELTA[2:(1+sum(REG.IN2))]=THRESH$seq.thresh[REG.IN2]
    INIT.DELTA[is.na(INIT.DELTA)]=max(INIT.DELTA,na.rm=T)
    jump=1:8
    DELTA=rbind(DELTA,cbind(count,day,series,jump,INIT.DELTA))
  }
}

DELTA=as.data.frame(DELTA)
names(DELTA)=c("n","day","series","jump","delta")
for(day in DAY.ID){
  for(series in LOOP.ID){
    DELTA$day[which(DELTA$day==day)]=DAYS[day]
    DELTA$series[which(DELTA$series==series)]=LOOPS[which(LOOP.ID==series)]
  }
}
DELTA$day=factor(DELTA$day,levels=DAYS)

png("Threshold5.png",width=800,height=1000)
ggplot(DELTA) +
  geom_step(aes(x=jump,y=delta),color="maroon") +
  facet_grid(series~day,switch="y") +
  xlab("Number of Regimes") +
  ylab("Occupancy Thresholds") +
  scale_y_continuous(position="right")  +
  scale_x_continuous(breaks=c(1,3,5,7,9))+
  theme_light()+
  theme(strip.background =element_rect(fill=rgb(0.5,0.5,0.5,1/4)))+
  theme(strip.text = element_text(colour = 'maroon',size=16)) +
  theme(axis.text.x=element_text(size=14)) +
  theme(axis.text.y=element_text(size=14)) +
  theme(axis.title.x=element_text(size=18)) +
  theme(axis.title.y=element_text(size=18))
dev.off()

#Specific Models
EX3=TAR.RESULTS3[[1]][[1]]$mean.coef
seq.thresh.func(day=1,series=1,P=7,h=5,nquant=50)$seq.thresh[c(0,33,48)]
EX3[seq.thresh.func(day=1,series=1,P=7,h=5,nquant=50)$Var.Groups %in% c(0,33,48)]
EX3[1:8]
EX3[]

logit.func(0.48)
logit.func(0.73)
