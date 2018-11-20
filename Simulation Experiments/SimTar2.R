#Args from Slurm File
argument <- commandArgs(trailingOnly = TRUE)

#Load Important Libraries
library(horseshoe)
library(coda)
library(foreach)
library(doParallel)

options(scipen=999)

#Functions for Projection Model Selection
lag.func<-function(x,k=1){
  t=length(x)
  y=c(rep(NA,t))
  for(i in (k+1):t){
    y[i]=x[i-k]
  }
  return(y)
}
tar.data.func=function(y,P,h,nquant){
  
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
  Large.X=Large.X
  
  #Output Data Frame
  data=data.frame(y=y.order,X=Large.X)
  
  #Out List
  out=list(data=data,Var.Groups=Var.Groups,t.order=t.order,seq.thresh=seq.thresh)
  
  return(out)
}
regime.proj.func<-function(post.coef,post.var,X,indproj){
  S=length(post.var)
  N=dim(X)[1]
  P=dim(X)[2]
  
  COEF.PROJ=matrix(0,S,P)
  X.proj=X[,indproj]
  pred.proj=X%*%post.coef
  
  #coef.proj=solve(t(X.proj)%*%X.proj)%*%t(X.proj)%*%pred.proj
  coef.proj=chol2inv(chol(t(X.proj)%*%X.proj))%*%t(X.proj)%*%pred.proj
  var.proj=post.var+colMeans((pred.proj-X.proj%*%coef.proj)^2)
  KL.PROJ=0.5*log(var.proj/post.var)
  COEF.PROJ[,indproj]=t(coef.proj)
  KL.MEAN=mean(KL.PROJ)
  COEF.MEAN=colMeans(COEF.PROJ)
  VAR.PROJ=var.proj
  VAR.MEAN=mean(var.proj)
  return(list(KL.MEAN=KL.MEAN,COEF.MEAN=COEF.MEAN,COEF.PROJ=COEF.PROJ,VAR.MEAN=VAR.MEAN,VAR.PROJ=VAR.PROJ))
}

#Simulation 1
a1<-0.5
a2<--0.5
delta<--0.4
sigma=2

N=1000
burn=1000

set.seed(argument[1])
y=rnorm((N+burn),0,1)
e=rnorm((N+burn),0,sigma)
for(i in 2:(N+burn)){
  if(y[i-1]<delta){
    y[i]=a1*y[i-1]+e[i]
  }else{
    y[i]=a2*y[i-1]+e[i]
  }
}

y=y[-(1:burn)]

#Create High Dimensional Model Matrix
DATA=tar.data.func(y=y,P=3,h=1,nquant=50)

#Initial Seasonal Model Using JAGS
tarmod=horseshoe(y=as.vector(DATA$data$y),X=as.matrix(DATA$data[,-1]),method.tau="halfCauchy",
                 method.sigma="Jeffreys",burn=5000,nmc=10000,thin=20)
fullpost.coef=tarmod$BetaSamples
fullpost.s2=as.numeric(tarmod$Sigma2Samples)

#Number of Params
nparam=dim(fullpost.coef)[1]

eff.size=rep(NA,nparam)
for(k in 1:nparam){
  eff.size[k]=effectiveSize(fullpost.coef[k,])
}
eff.size=c(eff.size,effectiveSize(fullpost.s2))

min.eff.size=min(eff.size)

#Initial Bayesian Shrinkage Regression
OUT=list(tar.data=DATA, min.eff.size=min.eff.size,
         fullpost.coef=fullpost.coef,fullpost.s2=fullpost.s2)

#Possible Set of Regimes(Group Parameters Based on Regime)
nquant=50                                      #Max Number of Regimes Considered
P=3                                            #AR Order in Each Regime
Var.Groups=c(rep(1:(nquant+1),each=(P+1))) #Variable Groups

#Retrieve Results From Full HS+ Estimated Model
FULL.POST.COEF=OUT$fullpost.coef
FULL.POST.S2=OUT$fullpost.s2

#Get Model Matrix
X=as.matrix(DATA$data[,-1])

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
  val=foreach(j=1:nleft,.combine=c,.export="regime.proj.func")%do%{
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

tar.coef=regime.proj.func(post.coef=FULL.POST.COEF,post.var=FULL.POST.S2,X=X,
                          indproj=which(Var.Groups %in% chosen))$COEF.MEAN
tar.s2=regime.proj.func(post.coef=FULL.POST.COEF,post.var=FULL.POST.S2,X=X,
                        indproj=which(Var.Groups %in% chosen))$VAR.MEAN

OUT2=list(tar.data=DATA,
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

tar2.coef=regime.proj.func(post.coef=FULL.POST.COEF,post.var=FULL.POST.S2,X=X,
                           indproj=chosen2)$COEF.MEAN
tar2.s2=regime.proj.func(post.coef=FULL.POST.COEF,post.var=FULL.POST.S2,X=X,
                         indproj=chosen2)$VAR.MEAN

train.data=DATA
train.predict=as.numeric(as.matrix(train.data$data[,-1])%*%tar2.coef)
order.rev1=train.data$t.order
order.forward1=order(order.rev1)

#Save Results From Variable Selection
tar.profile=train.predict[order.forward1]
actual.data=train.data$data$y[order.forward1]
sq.dev=(actual.data-tar.profile)^2

tar.data2=data.frame(cbind(actual.data,tar.profile,sq.dev))
names(tar.data2)=c("Actual","Predict","Dev")

OUT3=list(tar.data=tar.data2,
          seq.thresh=DATA$seq.thresh,
          KL=KL.REGIME2[1:length(chosen2)],
          chosen=chosen2,
          fullpost.coef=regime.proj.func(post.coef=FULL.POST.COEF,post.var=FULL.POST.S2,X=X,
                                         indproj=chosen2)$COEF.PROJ,
          fullpost.s2=regime.proj.func(post.coef=FULL.POST.COEF,post.var=FULL.POST.S2,X=X,
                                       indproj=chosen2)$VAR.PROJ,
          mean.coef=tar2.coef,mean.s2=tar2.s2)

rm(list=setdiff(ls(), c("argument","OUT2","OUT3")))

save.image(paste("SimTar2_",argument[1],".Rdata",sep=""))