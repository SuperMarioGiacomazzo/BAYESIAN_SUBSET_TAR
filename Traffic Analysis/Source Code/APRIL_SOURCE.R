options(scipen=999)

#Source BayesReg Package for Functions
#Package is Currently Obsolete
source("bayesregsource.R")

##################
#Get Data
##################

#Read in All Data for 3 Min
April.3=read.csv(file="APRIL_MIN_3.csv",header=T)

#Split Data Into Days
April.3.Day=list()
for(k in 1:5){
  days=(2+k)+c(0,1,2,3)*7
  April.3.Day[[k]]=subset(April.3,DAY %in% days)
}

#Read in All Data for 6 Min
April.6=read.csv(file="APRIL_MIN_6.csv",header=T)

#Split Data Into Days
April.6.Day=list()
for(k in 1:5){
  days=(2+k)+c(0,1,2,3)*7
  April.6.Day[[k]]=subset(April.6,DAY %in% days)
}

#Read in All Data for 15 Min
April.15=read.csv(file="APRIL_MIN_15.csv",header=T)

#Split Data Into Days
April.15.Day=list()
for(k in 1:5){
  days=(2+k)+c(0,1,2,3)*7
  April.15.Day[[k]]=subset(April.15,DAY %in% days)
}

############################################
#Create Vectors for Different Traffic Checks
############################################
vol.names=names(April.15)[1:7]
occ.names=names(April.15)[1:7+7]
day.names=c("M","T","W","Th","F")

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

##############################
# Logit Functions
##############################
logit.func<-function(x) return(log(x/(1-x)))
revlogit.func<-function(x) return(exp(x)/(1+exp(x)))

