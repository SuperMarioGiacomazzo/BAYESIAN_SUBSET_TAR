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
setwd("D:/Mario Documents/Research/JAS/BAYESIAN_SUBSET_TAR/Traffic Analysis/Horizon 3")


#Load Data
load("TARTRAFFIC3.Rdata")

#Coding of Locations Based on Graphic
LOOPS=c("A","B","C","D","E","F","G")
LOOP.ID=c(1,7,2,3,4,5,6)

#Codeing Days
DAYS=day.names
DAY.ID=1:5

#Function to Obtain Actual Series for Forecasting Period
pre.actual.data.func<-function(day,series,P,h){
  L106.TR = April.3.Day[[day]]$L106_occupancy[(1:1440)]/100
  L101.TR = April.3.Day[[day]]$L101_occupancy[(1:1440)]/100
  L108.TR = April.3.Day[[day]]$L108_occupancy[(1:1440)]/100
  L104.TR = April.3.Day[[day]]$L104_occupancy[(1:1440)]/100
  L102.TR = April.3.Day[[day]]$L102_occupancy[(1:1440)]/100
  L107.TR = April.3.Day[[day]]$L107_occupancy[(1:1440)]/100
  L103.TR = April.3.Day[[day]]$L103_occupancy[(1:1440)]/100
  yA = L101.TR[-(1:(P + h - 1))]
  yB = L106.TR[-(1:(P + h - 1))]
  yC = L108.TR[-(1:(P + h - 1))]
  yD = L102.TR[-(1:(P + h - 1))]
  yE = L104.TR[-(1:(P + h - 1))]
  yF = L107.TR[-(1:(P + h - 1))]
  yG = L103.TR[-(1:(P + h - 1))]
  y=list(yA,yB,yC,yD,yE,yF,yG)[[series]]
  return(y)
}
actual.data.func<-function(day,series,P,h){
  L106.TR = April.3.Day[[day]]$L106_occupancy[-(1:1440)]/100
  L101.TR = April.3.Day[[day]]$L101_occupancy[-(1:1440)]/100
  L108.TR = April.3.Day[[day]]$L108_occupancy[-(1:1440)]/100
  L104.TR = April.3.Day[[day]]$L104_occupancy[-(1:1440)]/100
  L102.TR = April.3.Day[[day]]$L102_occupancy[-(1:1440)]/100
  L107.TR = April.3.Day[[day]]$L107_occupancy[-(1:1440)]/100
  L103.TR = April.3.Day[[day]]$L103_occupancy[-(1:1440)]/100
  yA = L101.TR[-(1:(P + h - 1))]
  yB = L106.TR[-(1:(P + h - 1))]
  yC = L108.TR[-(1:(P + h - 1))]
  yD = L102.TR[-(1:(P + h - 1))]
  yE = L104.TR[-(1:(P + h - 1))]
  yF = L107.TR[-(1:(P + h - 1))]
  yG = L103.TR[-(1:(P + h - 1))]
  y=list(yA,yB,yC,yD,yE,yF,yG)[[series]]
  return(y)
}

#Function to Obtain Point Estimates and Density Forecasts
results.func<-function(data,X,order){
  POST.COEF=data$fullpost.coef
  POST.S2=data$fullpost.s2
  S=length(POST.S2)
  N=dim(X)[1]
  
  dist=data.frame(matrix(NA,12,N))
  est=rep(NA,N)
  for(n in 1:N){
    MU=as.numeric(X[n,]%*%POST.COEF)
    PREDICT=rnorm(S,mean=MU,sd=sqrt(POST.S2))
    REV.PREDICT=revlogit.func(PREDICT)
    
    Q=quantile(REV.PREDICT,c(.05,.1,.15,.20,.25,.40,.60,.75,.80,.85,.90,.95))
    REV.MU=mean(REV.PREDICT)
    
    dist[,n]=Q
    est[n]=REV.MU
  }
  dist2=t(dist[,ORDER])
  est2=est[ORDER]
  return(list(dist2,est2))
}

#Possible Set of Regimes(Group Parameters Based on Regime)
nquant=50                           #Max Number of Regimes Considered
P=7                                 #AR Order in Each Regime
Var.Groups=c(rep(1:(nquant+1),each=(P+1))) #Variable Groups

#Forecasting Results
#POINT.EST1.3=NULL
#POINT.EST2.3=NULL
POINT.EST3.3=NULL
#RMSFE.TAR1.3=matrix(NA,5,7)
#RMSFE.TAR2.3=matrix(NA,5,7)
#RMSFE.TAR3.3=matrix(NA,5,7)
#MASE.TAR1.3=matrix(NA,5,7)
#MASE.TAR2.3=matrix(NA,5,7)
MASE.TAR3.3=matrix(NA,5,7)

for(day in DAY.ID){
  for(series in LOOP.ID){
  
  #Day and Loop Name for Particular Iteration
  DAY.NAME=DAYS[day]
  LOOP.NAME=LOOPS[series]
  
  #Actual Traffic Occupancy for Forecast Period
  PREACTUAL=pre.actual.data.func(day=day,series=series,P=7,h=3)
  ACTUAL=actual.data.func(day=day,series=series,P=7,h=3)
  
  #MAE for Training Period
  mae.3=mean( abs(PREACTUAL-lag.func(PREACTUAL,k=3)),na.rm=T )
  
  #Data Matrix for Day and Loop For Forecast Period
  X=as.matrix(cbind(1,tar.data.func2(day=day,series=series,P=7,h=3,nquant=50)$data[,-1]))
  REVORDER=tar.data.func2(day=day,series=series,P=7,h=3,nquant=50)$t.order
  ORDER=order(REVORDER)
  
  #Get Results from 3-stage Method
  #OUT1=TAR.RESULTS[[day]][[series]]
  #OUT2=TAR.RESULTS2[[day]][[series]]
  #OUT2$fullpost.coef=t(OUT2$fullpost.coef)
  OUT3=TAR.RESULTS3[[day]][[series]]
  OUT3$fullpost.coef=t(OUT3$fullpost.coef)

  #Get Point Estimates and Forecast Distributions
  #SUMMARY1=results.func(data=OUT1,X=X,order=ORDER)
  #SUMMARY2=results.func(data=OUT2,X=X,order=ORDER)
  SUMMARY3=results.func(data=OUT3,X=X,order=ORDER)
  
  #Create Data Frames for All 3-Stages
  #POINT.EST1.3=rbind(POINT.EST1.3,cbind(day,series,ACTUAL,SUMMARY1[[1]],SUMMARY1[[2]]))
  #RMSFE.TAR1.3[which(DAY.ID==day),which(LOOP.ID==series)]=sqrt(mean((ACTUAL-SUMMARY1[[2]])^2))
  #MASE.TAR1.3[which(DAY.ID==day),which(LOOP.ID==series)]=mean( abs((ACTUAL-SUMMARY1[[2]]))/mae.3)
  
  #POINT.EST2.3=rbind(POINT.EST2.3,cbind(day,series,ACTUAL,SUMMARY2[[1]],SUMMARY2[[2]]))
  #RMSFE.TAR2.3[which(DAY.ID==day),which(LOOP.ID==series)]=sqrt(mean((ACTUAL-SUMMARY2[[2]])^2))
  #MASE.TAR2.3[which(DAY.ID==day),which(LOOP.ID==series)]=mean( abs((ACTUAL-SUMMARY2[[2]]))/mae.3)
  
  POINT.EST3.3=rbind(POINT.EST3.3,cbind(day,series,ACTUAL,SUMMARY3[[1]],SUMMARY3[[2]]))
  #RMSFE.TAR3.3[which(DAY.ID==day),which(LOOP.ID==series)]=sqrt(mean((ACTUAL-SUMMARY3[[2]])^2))  
  MASE.TAR3.3[which(DAY.ID==day),which(LOOP.ID==series)]=mean( abs((ACTUAL-SUMMARY3[[2]]))/mae.3)
    
  }
}

save.image("TARTRAFFIC3C.Rdata")
load("TARTRAFFIC3C.Rdata")

#POINT.EST1.3=data.frame(POINT.EST1.3)
#names(POINT.EST1.3)=c("DAY","SERIES","ACTUAL",c("5","10","15","20","25","40","60","75","80","85","90","95"),"MEAN")
#POINT.EST2.3=data.frame(POINT.EST2.3)
#names(POINT.EST2.3)=c("DAY","SERIES","ACTUAL",c("5","10","15","20","25","40","60","75","80","85","90","95"),"MEAN")
POINT.EST3.3=data.frame(POINT.EST3.3)
names(POINT.EST3.3)=c("DAY","SERIES","ACTUAL",c("5","10","15","20","25","40","60","75","80","85","90","95"),"MEAN")

#plot(subset(POINT.EST1.3,SERIES==1)[,"ACTUAL"],col="gray50",xlab="",ylab="Occupancy")
#points(subset(POINT.EST1.3,SERIES==1)[,"MEAN"],type="l",col="black")
#abline(v=((length(subset(POINT.EST1.3,SERIES==1)[,"ACTUAL"])/5)*(1:4)),lwd=2,lty=2)

#plot(subset(POINT.EST2.3,SERIES==1)[,"ACTUAL"],col="gray50",xlab="",ylab="Occupancy")
#points(subset(POINT.EST2.3,SERIES==1)[,"MEAN"],type="l",col="black")
#abline(v=((length(subset(POINT.EST2.3,SERIES==1)[,"ACTUAL"])/5)*(1:4)),lwd=2,lty=2)




png(filename="EST3Plots.png",width=1200,height=1200)
par(mfrow=c(7,1),mar=c(1.9,3.9,2.3,.5))
for(series in LOOP.ID){
  plot(subset(POINT.EST3.3,SERIES==series)[,"ACTUAL"],col="gray50",ylim=c(0,1),
       xlab="",xaxt="n",cex.main=2,cex.axis=2,ylab="",yaxt="n",pch=16,cex=1,
       main=paste("Location",LOOPS[which(LOOP.ID==series)]))
  points(subset(POINT.EST3.3,SERIES==series)[,"MEAN"],type="l",col="maroon",lwd=2)
  abline(v=((length(subset(POINT.EST3.3,SERIES==series)[,"ACTUAL"])/5)*(1:4)),lwd=2,lty=2)
  axis(side=1, labels=DAYS,cex.axis=2,
       at=(1:5)*(length(subset(POINT.EST3.3,SERIES==series)[,"ACTUAL"])/5)-471/2)
  axis(2,at=c(0,0.5,1),label=c(0,0.5,1),cex.axis=2)
}
dev.off()

png(filename="DENS3Plots.png",width=1200,height=1200)
par(mfrow=c(7,1),mar=c(1.9,3.9,2.3,.5))
for(series in LOOP.ID){
  POINT.EST=(POINT.EST3.3[POINT.EST3.3$SERIES==series,])$MEAN
  DENSITY=t((POINT.EST3.3[POINT.EST3.3$SERIES==series,])[,c("5","10","15","20","25","40","60","75","80","85","90","95")])
  ACTUAL=(POINT.EST3.3[POINT.EST3.3$SERIES==series,])$ACTUAL
  
  plot(POINT.EST,type="n",main=paste("Location",LOOPS[which(LOOP.ID==series)]),
       xlab="",ylab="",xaxt="n",yaxt="n",cex.main=2,cex.axis=2)
  fan(DENSITY,data.type="values",probs=c(.05,.1,.15,.20,.25,.40,.60,.75,.80,.85,.90,.95),
      alpha=0.2,llab=F,rlab=F)
  points(ACTUAL,type="p",pch=16,col="gray50",cex=1)
  abline(v=((length(subset(POINT.EST3.3,SERIES==series)[,"ACTUAL"])/5)*(1:4)),lwd=2,lty=2)
  axis(side=1, labels=DAYS,cex.axis=2,
       at=(1:5)*(length(subset(POINT.EST3.3,SERIES==series)[,"ACTUAL"])/5)-471/2)
  axis(2,at=c(0,0.5,1),label=c(0,0.5,1),cex.axis=2)
}
dev.off()


xtable(MASE.TAR3.3)
write.csv(round(MASE.TAR3.3,2),file="MASETAR33.csv")













