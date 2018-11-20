#Libraries Used
library(doParallel)
library(forecast)
library(foreach)
library(runjags)
library(bayesreg)
library(coda)
library("plot3D")

options(scipen=999)

setwd("D:/Mario Documents/Research/JAS/BAYESIAN_SUBSET_TAR/Traffic Analysis/Seasonal"))

load("SEASTRAFFIC.Rdata")

#Loop Names
loop.names=c("L106","L101","L108","L104","L102","L107","L103")

#Plots of Seasonal Time Series
png(filename="seasplots.png",width=2000,height=1800)
par(mfrow=c(5,7))
for(day in 1:5){
  for(series in 1:7){
    subdat=SEASMOD.RESULTS[[day]][[series]]$seas.data
    plot(subdat$Actual[1:960],ylim=c(min(subdat$Actual),max(subdat$Actual)),col="maroon",
         type="n",xlab="Time",ylab="Transformed Occpancy")
    points(1:480,subdat$Actual[1:480],col="blue",type="l")
    points(1:480,subdat$Actual[481:960],col="maroon",type="l")
    points(1:480,subdat$Actual[961:1440],col="gold",type="l")
    points(481:960,subdat$Actual[-(1:1440)],col="gray50",type="l")
    points(subdat$Seas[1:960],type="l",lwd=4,col="black")
    abline(v=480,lwd=4)
    }
}
dev.off()

png(filename="seasplots2.png",width=2000,height=1800)
par(mfrow=c(5,7),mar=c(0.2,4,1.5,0.2),mgp=c(0.5,0,0))
for(day in 1:5){
  for(series in 1:7){
    subdat=SEASMOD.RESULTS[[day]][[series]]$seas.data
    plot(subdat$Actual[1:1440],ylim=c(min(subdat$Actual),max(subdat$Actual)),col="maroon",
         type="n",main=loop.names[series],ylab=day.names[day],xlab="",
         cex.lab=2,cex.axis=2,xaxt="n",yaxt="n",cex.main=2)
    points(1:1440,subdat$Actual[1:1440],col="maroon")
    points(subdat$Seas[1:1440],type="l",lwd=4,col="black")
  }
}
dev.off()

png(filename="seasplots3.png",width=2000,height=1800)
par(mfrow=c(5,7),mar=c(0.2,4,1.5,0.2),mgp=c(0.5,0,0))
for(day in 1:5){
  for(series in 1:7){
    subdat=SEASMOD.RESULTS[[day]][[series]]$seas.data
    plot(subdat$Actual[1:480],ylim=c(min(subdat$Actual),max(subdat$Actual)),col="maroon",
         type="n",main=loop.names[series],ylab=day.names[day],xlab="",
         cex.lab=2,cex.axis=2,xaxt="n",yaxt="n",cex.main=2)
    points(1:480,subdat$Actual[1:480],col="maroon",pch=16,cex=1)
    points(481:960,subdat$Actual[1:480],col="maroon",pch=16,cex=1)
    points(961:1440,subdat$Actual[1:480],col="maroon",pch=16,cex=1)
    points(subdat$Seas[1:480],type="l",lwd=4,col="black")
  }
}
dev.off()


#Gather Deviations from Seasonal Profiles
png(filename="devplots1.png",width=2000,height=1800)
par(mfrow=c(5,7),mar=c(0.2,4,1.5,0.2),mgp=c(0.5,0,0))
for(day in 1:5){
  for(series in 1:7){
    subdat=SEASMOD.RESULTS[[day]][[series]]$seas.data$Seas.Dev
    plot(subdat,type="l",xaxt="n",yaxt="n",
         main=loop.names[series],ylab=day.names[day],xlab="",
         cex.lab=2,cex.main=2)
    points(subdat$Seas[1:480],type="l",lwd=4,col="black")
  }
}
dev.off()


occ.names2=occ.names[c(5,1,7,4,2,6,3)]
vol.names2=vol.names[c(5,1,7,4,2,6,3)]

png(filename="3Ddevplots.png",width=2000,height=1800)
par(mfrow=c(5,7),mar=c(0.5,0.5,0.5,0.1),mgp=c(3,1,0),cex.lab=3)
for(day in 1:5){
  for(series in 1:7){
    y=SEASMOD.RESULTS[[day]][[series]]$seas.data$Seas.Dev2[1:1440]
    x1=lag.func(April.3.Day[[day]][1:1440,vol.names2[series]],k=1)
    x2=lag.func(April.3.Day[[day]][1:1440,occ.names2[series]],k=1)
    scatter3D(x=x2,y=x1,z=y,phi=25,theta=-45,pch=16,xlab="V",ylab="O",zlab="D",
              col = ramp.col(c("maroon", "gold")),colkey = F)
  }
}
dev.off()

png(filename="3Ddev2plots.png",width=2000,height=1800)
par(mfrow=c(5,7),mar=c(0.5,0.5,0.5,0.1),mgp=c(3,1,0),cex.lab=3)
for(day in 1:5){
  for(series in 1:7){
    y=SEASMOD.RESULTS[[day]][[series]]$seas.data$Seas.Dev2[1:1440]
    x1=lag.func(April.3.Day[[day]][1:1440,vol.names2[series]],k=1)
    x2=lag.func(April.3.Day[[day]][1:1440,occ.names2[series]],k=1)
    scatter3D(x=x2,y=x1,z=y,phi=25,theta=-45,pch=16,xlab="V",ylab="O",zlab="D2",
              col = ramp.col(c("maroon", "gold")),colkey = F)
  }
}
dev.off()

png(filename="DevPlots.png",width=2000,height=1800)
par(mfrow=c(5,7),mar=c(0.2,4,1.5,0.2),mgp=c(0.5,0,0))
for(day in 1:5){
  for(series in 1:7){
    y=SEASMOD.RESULTS[[day]][[series]]$seas.data$Seas.Dev[1:1440]
    plot(y,type="l",xaxt="n",yaxt="n",
         main=loop.names[series],ylab=day.names[day],xlab="",
         cex.lab=2,cex.main=2)
  }
}
dev.off()

png(filename="Dev2Plots.png",width=2000,height=1800)
par(mfrow=c(5,7),mar=c(0.2,4,1.5,0.2),mgp=c(0.5,0,0))
for(day in 1:5){
  for(series in 1:7){
    y=SEASMOD.RESULTS[[day]][[series]]$seas.data$Seas.Dev2[1:1440]
    plot(y,type="l",xaxt="n",yaxt="n",
         main=loop.names[series],ylab=day.names[day],xlab="",
         cex.lab=2,cex.main=2)
  }
}
dev.off()

png(filename="DevVolPlots.png",width=2000,height=1800)
par(mfrow=c(5,7),mar=c(0.2,4,1.5,0.2),mgp=c(0.5,0,0))
for(day in 1:5){
  for(series in 1:7){
    y=SEASMOD.RESULTS[[day]][[series]]$seas.data$Seas.Dev[1:1440]
    x=lag.func(April.3.Day[[day]][1:1440,vol.names2[series]],k=1)
    plot(x,y,xaxt="n",yaxt="n",
         main=loop.names[series],ylab=day.names[day],xlab="",
         cex.lab=2,cex.main=2)
  }
}
dev.off()

png(filename="Dev2VolPlots.png",width=2000,height=1800)
par(mfrow=c(5,7),mar=c(0.2,4,1.5,0.2),mgp=c(0.5,0,0))
for(day in 1:5){
  for(series in 1:7){
    y=SEASMOD.RESULTS[[day]][[series]]$seas.data$Seas.Dev2[1:1440]
    x=lag.func(April.3.Day[[day]][1:1440,vol.names2[series]],k=1)
    plot(x,y,xaxt="n",yaxt="n",
         main=loop.names[series],ylab=day.names[day],xlab="",
         cex.lab=2,cex.main=2)
  }
}
dev.off()

png(filename="DevPACFPlots.png",width=2000,height=1800)
par(mfrow=c(5,7),mar=c(0.2,4,1.5,0.2),mgp=c(0.5,0,0))
for(day in 1:5){
  for(series in 1:7){
    y=SEASMOD.RESULTS[[day]][[series]]$seas.data$Seas.Dev[1:1440]/
      lag.func(SEASMOD.RESULTS[[day]][[series]]$seas.data$Seas.Dev2[1:1440],k=1)
    plot(y,xaxt="n",yaxt="n",type="l",
         main=loop.names[series],ylab=day.names[day],xlab="",
         cex.lab=2,cex.main=2)
    abline(h=0,col="maroon")
    pacf(na.omit(y))
  }
}
dev.off()



