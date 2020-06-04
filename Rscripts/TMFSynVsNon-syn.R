#Transition mutation Frequency Synonymous v. Non-Synonymous
#Uses TransFreqGLM to graph comparison between synonymous or nonsynonymous transition frequencies
library(ape)
library(seqinr)
library(pegas)
library(sfsmisc)
library(ggplot2)
library(scales)
library(plotrix)
library(RColorBrewer)
x1<-50
x2<--50
sdcard<- 10
limitx1<- 2*x1
Limitx2<- 2*x2


TMFData<-read.csv("./Output/TransFreqGLM.csv")
TMFData$Columns[which(TMFData$nonsyn == 0 & TMFData$stop == 0 & TMFData$Res == 0)]<- x1 #assign 1 if syn assign 2 if nonsyn
TMFData$Columns[which(TMFData$nonsyn == 1 & TMFData$stop == 0 & TMFData$Res == 0)]<- x2


plot(TMFData$Columns, TMFData$TransMutFreq, main="SynVSNon-Syn", ylab="Frequencies", type = "p", xlim = c(limitx2, limitx1))
r=TMFData$TransMutFreq[TMFData$Columns == x1]
r1<-rnorm(n = length(r), mean = x1, sd = sdcard)
s=TMFData$TransMutFreq[TMFData$Columns == x2]
s1<-rnorm(n = length(s), mean = x2, sd = sdcard)
points(r1,r,col= alpha("#97CC04",1))
points(s1,s,col= alpha("#7E7ED1",1))


dev.off()
