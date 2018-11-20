#Libraries Used
library(doParallel)
library(forecast)
library(foreach)
library(runjags)
library(bayesreg)
library(coda)

options(scipen=999)


#Open Data Directory
setwd("D:/Mario Documents/Research/JAS/BAYESIAN_SUBSET_TAR/Traffic Analysis/Source Code")

#Gather Data from Source Code
source("APRIL_SOURCE.R")

#Directory for Saving
setwd("D:/Mario Documents/Research/JAS/BAYESIAN_SUBSET_TAR/Traffic Analysis/Horizon 5")

##############################
# Lag Function
##############################
lag.func<-function(x,k=1){
  t=length(x)
  y=c(rep(NA,t))
  for(i in (k+1):t){
    y[i]=x[i-k]
  }
  return(y)
}
seasdiff.func<-function(x,d=1){
  t=length(x)
  y=diff(x,differences=1)
  y=c(rep(NA,d),y)
  return(y)
}

##############################
# Logit Functions
##############################
logit.func<-function(x) return(log(x/(1-x)))
revlogit.func<-function(x) return(exp(x)/(1+exp(x)))

###################################
# Data Functions
###################################

#For Training Data
tar.data.func=function(day,series,P,h,nquant){
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
  seq.thresh=as.numeric(quantile(X[,1],seq.quant))
  
  #Add Column of 1s for the Mean
  X.order=cbind(1,X.order)
  
  #Create Full Model Matrix
  Large.X=X.order
  for(Q in seq.thresh){
    X.edit=X.order
    j=max(which(X.order[,2]<Q))
    X.edit[(1:j),]=0
    Large.X=cbind(Large.X,X.edit)
  }
  
  #Create Variable Groups Based on Quantiles
  Var.Groups=rep(0:nquant,each=P)
  
  #Remove the First Column of 1s To Ensure Overall Mean is Not Shrunk
  Large.X=Large.X[,-1]
  
  #Output Data Frame
  data=data.frame(y=y.order,X=Large.X)
  
  #Out List
  out=list(data=data,Var.Groups=Var.Groups,t.order=t.order)
  
  return(out)
}
#For Testing Data
tar.data.func2=function(day,series,P,h,nquant){
  #Obtain Train Data for Detector 103 including Before (L108) and After (L101)
  L106.TR=April.3.Day[[day]]$L106_occupancy[-(1:1440)]/100
  L101.TR=April.3.Day[[day]]$L101_occupancy[-(1:1440)]/100
  L108.TR=April.3.Day[[day]]$L108_occupancy[-(1:1440)]/100
  L104.TR=April.3.Day[[day]]$L104_occupancy[-(1:1440)]/100
  L102.TR=April.3.Day[[day]]$L102_occupancy[-(1:1440)]/100
  L107.TR=April.3.Day[[day]]$L107_occupancy[-(1:1440)]/100
  L103.TR=April.3.Day[[day]]$L103_occupancy[-(1:1440)]/100
  
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
  seq.thresh=as.numeric(quantile(X[,1],seq.quant))
  
  #Add Column of 1s for the Mean
  X.order=cbind(1,X.order)
  
  #Create Full Model Matrix
  Large.X=X.order
  for(Q in seq.thresh){
    X.edit=X.order
    j=max(which(X.order[,2]<Q))
    X.edit[(1:j),]=0
    Large.X=cbind(Large.X,X.edit)
  }
  
  #Create Variable Groups Based on Quantiles
  Var.Groups=rep(0:nquant,each=P)
  
  #Remove the First Column of 1s To Ensure Overall Mean is Not Shrunk
  Large.X=Large.X[,-1]
  
  #Output Data Frame
  data=data.frame(y=y.order,X=Large.X)
  
  #Out List
  out=list(data=data,Var.Groups=Var.Groups,t.order=t.order)
  
  return(out)
}

####################################
# Create Empty Lists to Save Results
####################################
TAR.RESULTS=list()

for(day in 1:5){

  TAR.RESULTS[[day]]=foreach(series=1:7)%do%{
    
    #Obtain Full Matrix of Data
    out=tar.data.func(day=day,series=series,P=7,h=5,nquant=50)
    
    #Initial Seasonal Model Using JAGS
    tarmod=bayesreg(y~.,data=out$data,prior="hs+",nsamples=5000,burnin=5000,thin=10)
    
    #Obtain Full Posterior of Coefficients
    fullpost.coef=rbind(tarmod$beta0,tarmod$beta)
    fullpost.s2=as.numeric(tarmod$sigma2)
    
    #Number of Params
    nparam=dim(fullpost.coef)[1]
    
    eff.size=rep(NA,nparam)
    for(k in 1:nparam){
      eff.size[k]=effectiveSize(fullpost.coef[k,])
    }
    eff.size=c(eff.size,effectiveSize(fullpost.s2))
    
    min.eff.size=min(eff.size)
    
    tar.coef=rowMeans(fullpost.coef)
    tar.s2=mean(fullpost.s2)
    
    train.data=tar.data.func(day=day,series=series,P=7,h=5,nquant=50)
    train.predict=as.numeric(cbind(1,as.matrix(train.data$data[,-1]))%*%tar.coef)
    
    order.rev1=train.data$t.order
    order.forward1=order(order.rev1)
    
    png(filename=paste("Fit5D",day,"S",series,".png",sep=""),width=800,height=500)
    plot(train.data$data$y[order.forward1],col="gray50")
    points(train.predict[order.forward1],type="l")
    dev.off()
    
    test.data=tar.data.func2(day=day,series=series,P=7,h=5,nquant=50)
    test.predict=as.numeric(cbind(1,as.matrix(test.data$data[,-1]))%*%tar.coef)

    order.rev2=test.data$t.order
    order.forward2=order(order.rev2)
    
    png(filename=paste("Predict5D",day,"S",series,".png",sep=""),width=800,height=500)
    plot(test.data$data$y[order.forward2],col="gray50")
    points(test.predict[order.forward2],type="l")
    dev.off()
    
    tar.profile=c(train.predict[order.forward1],test.predict[order.forward2])
    actual.data=c(train.data$data$y[order.forward1],test.data$data$y[order.forward2])
    tar.dev=actual.data-tar.profile
    tar.dev2=tar.dev^2
    
    tar.data=data.frame(cbind(actual.data,tar.profile,tar.dev,tar.dev2))
    names(tar.data)=c("Actual","Predict","Dev","Dev2")
    
    #Optimal Seasonal Model
    OUT=list(tar.data=tar.data, min.eff.size=min.eff.size,
             fullpost.coef=fullpost.coef,fullpost.s2=fullpost.s2)
    OUT
  } 
  save.image("TARTRAFFIC5.Rdata")
}




