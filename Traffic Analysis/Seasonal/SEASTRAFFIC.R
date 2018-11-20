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
setwd("D:/Mario Documents/Research/JAS/BAYESIAN_SUBSET_TAR/Traffic Analysis/Seasonal")

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
seas.data.func=function(day,nfreq,series){
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
  
  #Time Variable
  time=(1:length(yA))
  
  #Creation of Harmonic Matrices
  N=length(yA)
  X=matrix(NA,N,nfreq*2)
  for (j in 1:nfreq){
    X[,(2*j-1):(2*j)]=cbind(cos(2*pi*time*j/480),sin(2*pi*time*j/480))
  }
  P=dim(X)[2]
  y=list(yA,yB,yC,yD,yE,yF,yG)[[series]]
  data=data.frame(y=y,X=X)
  return(data)
}
seas.data.func2=function(day,nfreq,series){
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
  
  #Time Variable
  time=(1:length(yA))
  
  #Creation of Harmonic Matrices
  N=length(yA)
  X=matrix(NA,N,nfreq*2)
  for (j in 1:nfreq){
    X[,(2*j-1):(2*j)]=cbind(cos(2*pi*time*j/480),sin(2*pi*time*j/480))
  }
  P=dim(X)[2]
  y=list(yA,yB,yC,yD,yE,yF,yG)[[series]]
  data=data.frame(y=y,X=X)
  return(data)
}

####################################
# Create Empty Lists to Save Results
####################################
SEASMOD.RESULTS=list()

for(day in 1:5){

  SEASMOD.RESULTS[[day]]=foreach(series=1:7)%do%{
    
    #Initial Seasonal Model Using JAGS
    seasmod=bayesreg(y~.,data=seas.data.func(day=day,nfreq=150,series=series),
                     prior="hs+",nsamples=2000,burnin=5000,thin=10)
    
    eff.size=rep(NA,300)
    for(k in 1:300){
      eff.size[k]=effectiveSize(c(seasmod$beta[k,]))
    }
    min.eff.size=min(eff.size)
    
    seas.mean=seasmod$muBeta0
    seas.coef=seasmod$muBeta
    
    train.data=seas.data.func(day=day,nfreq=150,series=series)
    train.predict=as.numeric(c(seas.mean)+as.matrix(train.data[,-1])%*%seas.coef)
    
    test.data=seas.data.func2(day=day,nfreq=150,series=series)
    test.predict=as.numeric(c(seas.mean)+as.matrix(test.data[,-1])%*%seas.coef)
    
    seas.profile=c(train.predict,test.predict)
    actual.data=c(train.data$y,test.data$y)
    seas.dev=actual.data-seas.profile
    seas.dev2=seas.dev^2
    
    seas.data=data.frame(cbind(actual.data,seas.profile,seas.dev,seas.dev2))
    names(seas.data)=c("Actual","Seas","Seas.Dev","Seas.Dev2")
    
    #Optimal Seasonal Model
    out=list(seas.data=seas.data,s.mu.p=as.numeric(seasmod$beta0),s.coef.p=seasmod$beta,s.s2.p=as.numeric(seasmod$sigma2))
    out
  } 
  save.image("SEASTRAFFIC.Rdata")
}






