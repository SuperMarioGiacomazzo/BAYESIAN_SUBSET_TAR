#Libraries Used
library(doParallel)
library(forecast)
library(foreach)
library(runjags)
library(bayesreg)
library(coda)
library("plot3D")
library(fanplot)
library(xtable)

options(scipen=999)

setwd("D:/Mario Documents/Research/JAS/BAYESIAN_SUBSET_TAR/Traffic Analysis/Seasonal")

load("SEASTRAFFIC.Rdata")

#Coding of Locations Based on Graphic
LOOPS=c("A","B","C","D","E","F","G")
LOOP.ID=c(1,7,2,3,4,5,6)

#Codeing Days
DAYS=day.names
DAY.ID=1:5

#Function to Obtain Actual Series for Forecasting Period
pre.actual.data.func<-function(day,series){
  L106.TR = April.3.Day[[day]]$L106_occupancy[(1:1440)]/100
  L101.TR = April.3.Day[[day]]$L101_occupancy[(1:1440)]/100
  L108.TR = April.3.Day[[day]]$L108_occupancy[(1:1440)]/100
  L104.TR = April.3.Day[[day]]$L104_occupancy[(1:1440)]/100
  L102.TR = April.3.Day[[day]]$L102_occupancy[(1:1440)]/100
  L107.TR = April.3.Day[[day]]$L107_occupancy[(1:1440)]/100
  L103.TR = April.3.Day[[day]]$L103_occupancy[(1:1440)]/100
  yA = L101.TR
  yB = L106.TR
  yC = L108.TR
  yD = L102.TR
  yE = L104.TR
  yF = L107.TR
  yG = L103.TR
  y=list(yA,yB,yC,yD,yE,yF,yG)[[series]]
  return(y)
}
actual.data.func<-function(day,series){
  L106.TR = April.3.Day[[day]]$L106_occupancy[-(1:1440)]/100
  L101.TR = April.3.Day[[day]]$L101_occupancy[-(1:1440)]/100
  L108.TR = April.3.Day[[day]]$L108_occupancy[-(1:1440)]/100
  L104.TR = April.3.Day[[day]]$L104_occupancy[-(1:1440)]/100
  L102.TR = April.3.Day[[day]]$L102_occupancy[-(1:1440)]/100
  L107.TR = April.3.Day[[day]]$L107_occupancy[-(1:1440)]/100
  L103.TR = April.3.Day[[day]]$L103_occupancy[-(1:1440)]/100
  yA = L101.TR
  yB = L106.TR
  yC = L108.TR
  yD = L102.TR
  yE = L104.TR
  yF = L107.TR
  yG = L103.TR
  y=list(yA,yB,yC,yD,yE,yF,yG)[[series]]
  return(y)
}

#Function to Obtain Point Estimates and Density Forecasts
results.func<-function(mu,coef,s2,X){
  S=length(mu)
  N=dim(X)[1]
  
  dist=data.frame(matrix(NA,12,N))
  est=rep(NA,N)
  for(n in 1:N){
    MU=mu+as.numeric(X[n,]%*%coef)
    PREDICT=rnorm(S,mean=MU,sd=sqrt(s2))
    REV.PREDICT=revlogit.func(PREDICT)
    
    Q=quantile(REV.PREDICT,c(.05,.1,.15,.20,.25,.40,.60,.75,.80,.85,.90,.95))
    REV.MU=mean(REV.PREDICT)
    
    dist[,n]=Q
    est[n]=REV.MU
  }
  dist2=t(dist)
  est2=est
  return(list(dist2,est2))
}

POINT.EST=NULL
MASE.SEAS1=matrix(NA,5,7)
MASE.SEAS3=matrix(NA,5,7)
MASE.SEAS5=matrix(NA,5,7)

for(day in DAY.ID){
  for(series in LOOP.ID){
    
    #Day and Loop Name for Particular Iteration
    DAY.NAME=DAYS[day]
    LOOP.NAME=LOOPS[series]
    
    #Actual Traffic Occupancy for Forecast Period
    PREACTUAL=pre.actual.data.func(day=day,series=series)
    ACTUAL=actual.data.func(day=day,series=series)
    
    #MAE for Training Period
    mae.1=mean( abs(PREACTUAL-lag.func(PREACTUAL,k=1)),na.rm=T )
    mae.3=mean( abs(PREACTUAL-lag.func(PREACTUAL,k=3)),na.rm=T )
    mae.5=mean( abs(PREACTUAL-lag.func(PREACTUAL,k=5)),na.rm=T )
    
    mu=SEASMOD.RESULTS[[day]][[series]]$s.mu.p
    coef=SEASMOD.RESULTS[[day]][[series]]$s.coef.p
    s2=SEASMOD.RESULTS[[day]][[series]]$s.s2.p
    
    #Data Matrix for Day and Loop For Forecast Period
    X=as.matrix(seas.data.func2(day=day,nfreq=150,series=series))[,-1]
    
    #Output Prediction Results on Original Scale
    SUMMARY3=results.func(mu=mu,coef=coef,s2=s2,X=X)
    
    
    POINT.EST=rbind(POINT.EST,cbind(day,series,ACTUAL,SUMMARY3[[1]],SUMMARY3[[2]]))
    MASE.SEAS1[which(DAY.ID==day),which(LOOP.ID==series)]=mean( abs((ACTUAL-SUMMARY3[[2]]))/mae.1)  
    MASE.SEAS3[which(DAY.ID==day),which(LOOP.ID==series)]=mean( abs((ACTUAL-SUMMARY3[[2]]))/mae.3)
    MASE.SEAS5[which(DAY.ID==day),which(LOOP.ID==series)]=mean( abs((ACTUAL-SUMMARY3[[2]]))/mae.5)
    
  }
}


save.image("SEASTRAFFIC4.Rdata")



POINT.EST=data.frame(POINT.EST)
names(POINT.EST)=c("DAY","SERIES","ACTUAL",c("5","10","15","20","25","40","60","75","80","85","90","95"),"MEAN")



series=1
png(filename="SEASESTPlots.png",width=1200,height=1200)
par(mfrow=c(7,1),mar=c(1.9,3.9,2.3,.5))
for(series in LOOP.ID){
  plot(subset(POINT.EST,SERIES==series)[,"ACTUAL"],col="gray50",ylim=c(0,1),
       xlab="",xaxt="n",cex.main=2,cex.axis=2,ylab="",yaxt="n",pch=16,cex=1,
       main=paste("Location",LOOPS[which(LOOP.ID==series)]))
  points(subset(POINT.EST,SERIES==series)[,"MEAN"],type="l",col="maroon",lwd=2)
  abline(v=((length(subset(POINT.EST,SERIES==series)[,"ACTUAL"])/5)*(1:4)),lwd=2,lty=2)
  axis(side=1, labels=DAYS,cex.axis=2,
       at=(1:5)*(length(subset(POINT.EST,SERIES==series)[,"ACTUAL"])/5)-480/2)
  axis(2,at=c(0,0.5,1),label=c(0,0.5,1),cex.axis=2)
}
dev.off()

png(filename="SEASDENSPlots.png",width=1200,height=1200)
par(mfrow=c(7,1),mar=c(1.9,3.9,2.3,.5))
for(series in LOOP.ID){
  POINT.EST2=(POINT.EST[POINT.EST$SERIES==series,])$MEAN
  DENSITY2=t((POINT.EST[POINT.EST$SERIES==series,])[,c("5","10","15","20","25","40","60","75","80","85","90","95")])
  ACTUAL2=(POINT.EST[POINT.EST$SERIES==series,])$ACTUAL
  
  plot(POINT.EST2,type="n",main=paste("Location",LOOPS[which(LOOP.ID==series)]),
       xlab="",ylab="",xaxt="n",yaxt="n",cex.main=2,cex.axis=2)
  fan(DENSITY2,data.type="values",probs=c(.05,.1,.15,.20,.25,.40,.60,.75,.80,.85,.90,.95),
      alpha=0.2,llab=F,rlab=F)
  points(ACTUAL2,type="p",pch=16,col="gray50",cex=1)
  abline(v=((length(subset(POINT.EST,SERIES==series)[,"ACTUAL"])/5)*(1:4)),lwd=2,lty=2)
  axis(side=1, labels=DAYS,cex.axis=2,
       at=(1:5)*(length(subset(POINT.EST,SERIES==series)[,"ACTUAL"])/5)-469/2)
  axis(2,at=c(0,0.5,1),label=c(0,0.5,1),cex.axis=2)
}
dev.off()

xtable(MASE.SEAS1)
write.csv(round(MASE.SEAS1,2),file="MASE.SEAS1.csv")
xtable(MASE.SEAS3)
write.csv(round(MASE.SEAS3,2),file="MASE.SEAS3.csv")
xtable(MASE.SEAS5)
write.csv(round(MASE.SEAS5,2),file="MASE.SEAS5.csv")


