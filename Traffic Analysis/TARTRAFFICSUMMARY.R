library(ggplot2)

setwd("D:/Mario Documents/Graduate School/Ph.D Dissertation/Traffic Analysis/TAR")

FINAL.MOD.NREGIME1=as.matrix(read.csv("FINAL.MOD.NREGIME1.csv")[,-1])
FINAL.MOD.NREGIME3=as.matrix(read.csv("FINAL.MOD.NREGIME3.csv")[,-1])
FINAL.MOD.NREGIME5=as.matrix(read.csv("FINAL.MOD.NREGIME5.csv")[,-1])

ALL.VEC=c(as.vector(FINAL.MOD.NREGIME1),as.vector(FINAL.MOD.NREGIME3),as.vector(FINAL.MOD.NREGIME5))
ALL.VEC=data.frame(x=as.factor(ALL.VEC))

png("TOTALREGIME.png",width=600,height=400)
ggplot(ALL.VEC)+
  geom_bar(aes(x=x,y=..count../sum(..count..)),fill="maroon") +
  theme_light() + xlab("Number of Regimes") + ylab("Relative Frequency")+
  theme(axis.text.x=element_text(size=14)) +
  theme(axis.text.y=element_text(size=14)) +
  theme(axis.title.x=element_text(size=18)) +
  theme(axis.title.y=element_text(size=18))
dev.off()

