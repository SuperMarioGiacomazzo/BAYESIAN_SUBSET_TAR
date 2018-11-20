options(scipen=999)

#Open Data Directory
setwd("D:/Mario Documents/Research/JAS/BAYESIAN_SUBSET_TAR/Traffic Analysis/Source Code")

#Gather Data from Source Code
source("APRIL_SOURCE.R")

#Directory for Saving
setwd("D:/Mario Documents/Research/JAS/BAYESIAN_SUBSET_TAR/Traffic Analysis/Traffic Plots")

#Coding of Locations Based on Graphic
LOOPS=c("A","B","C","D","E","F","G")
LOOP.ID=c(1,7,2,3,4,5,6)

#Codeing Days
DAYS=day.names
DAY.ID=1:5

#Logit Transformation Function
logit.func<-function(x) return(log(x/(1-x)))

#Function to Obtain Actual Series for Forecasting Period
data.func<-function(day,series){
  L106.TR = April.3.Day[[day]]$L106_occupancy/100
  L101.TR = April.3.Day[[day]]$L101_occupancy/100
  L108.TR = April.3.Day[[day]]$L108_occupancy/100
  L104.TR = April.3.Day[[day]]$L104_occupancy/100
  L102.TR = April.3.Day[[day]]$L102_occupancy/100
  L107.TR = April.3.Day[[day]]$L107_occupancy/100
  L103.TR = April.3.Day[[day]]$L103_occupancy/100
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

data2.func<-function(day,series){
  L106.TR = April.3.Day[[day]]$L106_occupancy/100
  L101.TR = April.3.Day[[day]]$L101_occupancy/100
  L108.TR = April.3.Day[[day]]$L108_occupancy/100
  L104.TR = April.3.Day[[day]]$L104_occupancy/100
  L102.TR = April.3.Day[[day]]$L102_occupancy/100
  L107.TR = April.3.Day[[day]]$L107_occupancy/100
  L103.TR = April.3.Day[[day]]$L103_occupancy/100
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
  y=list(yA,yB,yC,yD,yE,yF,yG)[[series]]
  return(y)
}

#Plot of All Data
png(filename="RawPlots.png",width=1000,height=1200)
par(mfrow=c(7,1),mar=c(2.1,3.9,2.1,.5))
for(series in LOOP.ID){
  data.1=data.func(day=1,series=series)
  data.2=data.func(day=2,series=series)
  data.3=data.func(day=3,series=series)
  data.4=data.func(day=4,series=series)
  data.5=data.func(day=5,series=series)
  length(full.data)/480
  full.data=c(data.1[1:480],data.2[1:480],data.3[1:480],data.4[1:480],data.5[1:480],
              data.1[481:960],data.2[481:960],data.3[481:960],data.4[481:960],data.5[481:960],
              data.1[961:1440],data.2[961:1440],data.3[961:1440],data.4[961:1440],data.5[961:1440],
              data.1[-(1:1440)],data.2[-(1:1440)],data.3[-(1:1440)],data.4[-(1:1440)],data.5[-(1:1440)])
  
  plot(full.data,xlab="",ylab="",ylim=c(0,1),xaxt="n",yaxt="n",type="l",col="maroon",
       main=paste("Raw Occupancy for Location",LOOPS[which(LOOP.ID==series)]),cex.main=2)
  axis(1,at=240+480*(0:19),label=rep(DAYS,4),cex.axis=2)
  axis(2,at=c(0,0.25,0.5,0.75,1),label=c(0,0.25,0.5,0.75,1),cex.axis=2)
  rect(7201,0,9600,1,col = rgb(0.5,0.5,0.5,1/4))
}
dev.off()

#Plot of All Transformed Data
png(filename="TransPlots.png",width=1000,height=1200)
par(mfrow=c(7,1),mar=c(2.1,3.9,2.1,.5))
for(series in LOOP.ID){
  data.1=data2.func(day=1,series=series)
  data.2=data2.func(day=2,series=series)
  data.3=data2.func(day=3,series=series)
  data.4=data2.func(day=4,series=series)
  data.5=data2.func(day=5,series=series)
  length(full.data)/480
  full.data=c(data.1[1:480],data.2[1:480],data.3[1:480],data.4[1:480],data.5[1:480],
              data.1[481:960],data.2[481:960],data.3[481:960],data.4[481:960],data.5[481:960],
              data.1[961:1440],data.2[961:1440],data.3[961:1440],data.4[961:1440],data.5[961:1440],
              data.1[-(1:1440)],data.2[-(1:1440)],data.3[-(1:1440)],data.4[-(1:1440)],data.5[-(1:1440)])
  
  plot(full.data,xlab="",ylab="",xaxt="n",type="l",col="maroon",cex.axis=2,
       main=paste("Transformed Occupancy for Location",LOOPS[which(LOOP.ID==series)]),cex.main=2)
  axis(1,at=240+480*(0:19),label=rep(DAYS,4),cex.axis=2)
  rect(7201,min(full.data),9600,max(full.data),col = rgb(0.5,0.5,0.5,1/4))
}
dev.off()








