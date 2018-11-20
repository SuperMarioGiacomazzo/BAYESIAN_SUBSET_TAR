#Libraries Used
library(doParallel)
library(forecast)
library(foreach)
library(bayesreg)
library(coda)

options(scipen=999)

#Working Directory
setwd("D:/Mario Documents/Research/JAS/BAYESIAN_SUBSET_TAR/Traffic Analysis/Horizon 5")

#Load Data
load("TARTRAFFIC5.Rdata")

#Required Functions for Regime Projection Model Selection
regime.proj.func<-function(post.coef,post.var,X,indproj){
  S=length(post.var)
  N=dim(X)[1]
  P=dim(X)[2]
  
  COEF.PROJ=matrix(0,S,P)
  X.proj=X[,indproj]
  pred.proj=X%*%post.coef
  
  coef.proj=solve(t(X.proj)%*%X.proj)%*%t(X.proj)%*%pred.proj
  var.proj=post.var+colMeans((pred.proj-X.proj%*%coef.proj)^2)
  KL.PROJ=0.5*log(var.proj/post.var)
  COEF.PROJ[,indproj]=t(coef.proj)
  KL.MEAN=mean(KL.PROJ)
  COEF.MEAN=colMeans(COEF.PROJ)
  VAR.PROJ=var.proj
  VAR.MEAN=mean(var.proj)
  return(list(KL.MEAN=KL.MEAN,COEF.MEAN=COEF.MEAN,COEF.PROJ=COEF.PROJ,VAR.MEAN=VAR.MEAN,VAR.PROJ=VAR.PROJ))
}

#Possible Set of Regimes(Group Parameters Based on Regime)
nquant=50                           #Max Number of Regimes Considered
P=7                                 #AR Order in Each Regime
Var.Groups=c(rep(1:(nquant+1),each=(P+1))) #Variable Groups

#Projection Model Selection to Select Regimes

TAR.RESULTS2=list()
TAR.RESULTS3=list()
for(day in 1:5){
  TAR.RESULTS2[[day]]=list()
  TAR.RESULTS3[[day]]=list()
  for(series in 1:7){
    
    #Retrieve Results From Full HS+ Estimated Model
    OUT=TAR.RESULTS[[day]][[series]]
    FULL.POST.COEF=OUT$fullpost.coef
    FULL.POST.S2=OUT$fullpost.s2
    
    DATA.TRAIN=tar.data.func(day=day,series=series,P=P,h=5,nquant=nquant)
    X=as.matrix(cbind(1,DATA.TRAIN$data[,-1]))
    
    
    #Regime Selection
    n.regimes=length(unique(Var.Groups))
    KL.REGIME=rep(NA,n.regimes)
    chosen=1
    notchosen=setdiff(1:n.regimes,chosen)
  
    FIRST=regime.proj.func(post.coef=FULL.POST.COEF,post.var=FULL.POST.S2,X=X,indproj=which(Var.Groups %in% chosen))
    KL.REGIME[1]=FIRST$KL.MEAN
    modnum=2
    check=0
    
    cl<-makeCluster(5)
    registerDoParallel(cl)
    while(check<0.95 & modnum<=n.regimes){
      nleft<-length(notchosen)
      val=foreach(j=1:nleft,.combine=c,.export="regime.proj.func")%dopar%{
        ind<-sort(c(chosen,notchosen[j]))
        NEXT<-tryCatch({regime.proj.func(post.coef=FULL.POST.COEF,post.var=FULL.POST.S2,X=X,
                        indproj=which(Var.Groups %in% ind))$KL.MEAN},
                        error=function(e){return(NA)})
        NEXT
      }
      if(sum(is.na(val))!=length(val)){
        minval<-which.min(val)
        chosen<-c(chosen,notchosen[minval])
        notchosen<-setdiff(1:n.regimes,chosen)
        KL.REGIME[modnum]<-val[minval]
        check=1-KL.REGIME[modnum]/KL.REGIME[1]
        modnum=modnum+1
      }else{
        check=9999
        modnum=9999
      }
    }
    stopCluster(cl)
    
    
    #Plots After Regime Selection
    tar.coef=regime.proj.func(post.coef=FULL.POST.COEF,post.var=FULL.POST.S2,X=X,
                              indproj=which(Var.Groups %in% chosen))$COEF.MEAN
    tar.s2=regime.proj.func(post.coef=FULL.POST.COEF,post.var=FULL.POST.S2,X=X,
                            indproj=which(Var.Groups %in% chosen))$VAR.MEAN
    
    train.data=tar.data.func(day=day,series=series,P=7,h=5,nquant=50)
    train.predict=as.numeric(cbind(1,as.matrix(train.data$data[,-1]))%*%tar.coef)
    
    order.rev1=train.data$t.order
    order.forward1=order(order.rev1)
    
    png(filename=paste("SELECT1Fit5D",day,"S",series,".png",sep=""),width=800,height=500)
    plot(train.data$data$y[order.forward1],col="gray50")
    points(train.predict[order.forward1],type="l")
    dev.off()
    
    test.data=tar.data.func2(day=day,series=series,P=7,h=5,nquant=50)
    test.predict=as.numeric(cbind(1,as.matrix(test.data$data[,-1]))%*%tar.coef)
    
    order.rev2=test.data$t.order
    order.forward2=order(order.rev2)
    
    png(filename=paste("SELECT1Predict5D",day,"S",series,".png",sep=""),width=800,height=500)
    plot(test.data$data$y[order.forward2],col="gray50")
    points(test.predict[order.forward2],type="l")
    dev.off()
    
    #Save Results From Regime Selection
    tar.profile=c(train.predict[order.forward1],test.predict[order.forward2])
    actual.data=c(train.data$data$y[order.forward1],test.data$data$y[order.forward2])
    tar.dev=actual.data-tar.profile
    tar.dev2=tar.dev^2
    
    tar.data=data.frame(cbind(actual.data,tar.profile,tar.dev,tar.dev2))
    names(tar.data)=c("Actual","Predict","Dev","Dev2")

    TAR.RESULTS2[[day]][[series]]=list(tar.data=tar.data,
                                  regimes=chosen,
                                  KL=KL.REGIME[1:length(chosen)],
                                  chosen=which(Var.Groups %in% chosen),
                                  fullpost.coef=regime.proj.func(post.coef=FULL.POST.COEF,post.var=FULL.POST.S2,X=X,
                                                                 indproj=which(Var.Groups %in% chosen))$COEF.PROJ,
                                  fullpost.s2=regime.proj.func(post.coef=FULL.POST.COEF,post.var=FULL.POST.S2,X=X,
                                                               indproj=which(Var.Groups %in% chosen))$VAR.PROJ,
                                  mean.coef=tar.coef,mean.s2=tar.s2)

    
    #Variable Selection
    possible.vars=which(Var.Groups %in% chosen)
    n.vars=length(possible.vars)
    KL.REGIME2=rep(NA,n.vars)
    chosen2=1
    notchosen2=setdiff(1:n.vars,chosen2)
    
    FIRST2=regime.proj.func(post.coef=FULL.POST.COEF,post.var=FULL.POST.S2,X=X,indproj=possible.vars[chosen2])
    KL.REGIME2[1]=FIRST2$KL.MEAN
    modnum2=2
    check2=0
    
    cl<-makeCluster(5)
    registerDoParallel(cl)
    while(check2<0.95 & modnum2<=n.vars){
      nleft2<-length(notchosen2)
      val2=foreach(j=1:nleft2,.combine=c,.export="regime.proj.func")%dopar%{
        ind<-sort(c(chosen2,notchosen2[j]))
        NEXT<-tryCatch({regime.proj.func(post.coef=FULL.POST.COEF,post.var=FULL.POST.S2,X=X,
                                         indproj=possible.vars[ind])$KL.MEAN},
                       error=function(e){return(NA)})
        NEXT
      }
      minval2<-which.min(val2)
      chosen2<-c(chosen2,notchosen2[minval2])
      notchosen2<-setdiff(1:n.vars,chosen2)
      KL.REGIME2[modnum2]<-val2[minval2]
      check2=1-KL.REGIME2[modnum2]/KL.REGIME2[1]
      modnum2=modnum2+1
    }
    stopCluster(cl)
    
    chosen2=possible.vars[chosen2]
    
    #Plots After Variable Selection
    tar2.coef=regime.proj.func(post.coef=FULL.POST.COEF,post.var=FULL.POST.S2,X=X,
                              indproj=chosen2)$COEF.MEAN
    tar2.s2=regime.proj.func(post.coef=FULL.POST.COEF,post.var=FULL.POST.S2,X=X,
                            indproj=chosen2)$VAR.MEAN
    
    train.data=tar.data.func(day=day,series=series,P=7,h=5,nquant=50)
    train.predict2=as.numeric(cbind(1,as.matrix(train.data$data[,-1]))%*%tar2.coef)
    
    order.rev1=train.data$t.order
    order.forward1=order(order.rev1)
    
    png(filename=paste("SELECT2Fit5D",day,"S",series,".png",sep=""),width=800,height=500)
    plot(train.data$data$y[order.forward1],col="gray50")
    points(train.predict2[order.forward1],type="l")
    dev.off()
    
    test.data=tar.data.func2(day=day,series=series,P=7,h=5,nquant=50)
    test.predict2=as.numeric(cbind(1,as.matrix(test.data$data[,-1]))%*%tar2.coef)
    
    order.rev2=test.data$t.order
    order.forward2=order(order.rev2)
    
    png(filename=paste("SELECT2Predict5D",day,"S",series,".png",sep=""),width=800,height=500)
    plot(test.data$data$y[order.forward2],col="gray50")
    points(test.predict2[order.forward2],type="l")
    dev.off()
    
    
    #Save Results From Variable Selection
    tar.profile2=c(train.predict2[order.forward1],test.predict2[order.forward2])
    actual.data=c(train.data$data$y[order.forward1],test.data$data$y[order.forward2])
    tar.dev2=actual.data-tar.profile2
    tar.dev22=tar.dev2^2
    
    tar.data2=data.frame(cbind(actual.data,tar.profile2,tar.dev2,tar.dev22))
    names(tar.data2)=c("Actual","Predict","Dev","Dev2")
    
    TAR.RESULTS3[[day]][[series]]=list(tar.data=tar.data2,
                                       KL=KL.REGIME2[1:length(chosen2)],
                                       chosen=chosen2,
                                       fullpost.coef=regime.proj.func(post.coef=FULL.POST.COEF,post.var=FULL.POST.S2,X=X,
                                                                      indproj=chosen2)$COEF.PROJ,
                                       fullpost.s2=regime.proj.func(post.coef=FULL.POST.COEF,post.var=FULL.POST.S2,X=X,
                                                                    indproj=chosen2)$VAR.PROJ,
                                       mean.coef=tar2.coef,mean.s2=tar2.s2)
  
    save.image("TARTRAFFIC5B.Rdata")
  }
}


TAR.RESULTS=NULL
fullpost.coef=NULL
save.image("TARTRAFFIC5B.Rdata")










