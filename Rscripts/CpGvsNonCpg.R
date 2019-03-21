library(ape)
library(seqinr)
library(pegas)
library(sfsmisc)
library(ggplot2)
library(scales)
library(plotrix)
library(RColorBrewer)
#Makes graphs for T->C and A->G

OverviewDFZanini<-read.table("Output/OverviewSelCoeffZanini.csv",sep=",",header=TRUE,row.names=1)

cpg.y<-subset(OverviewDFZanini, makesCpG==1)
cpg.n<-subset(OverviewDFZanini, makesCpG==0)

AC<-subset(cpg.y, WTnt=='a') #Does make cpg
ANC<-subset(cpg.n, WTnt=='a') #Doesn't make cpg
GC<-subset(cpg.y, WTnt=='g') #Does make cpg
GNC<-subset(cpg.n, WTnt=='g') #Doesn't make cpg
TC<-subset(cpg.y, WTnt=='t') #Does make cpg
TNC<-subset(cpg.n, WTnt=='t') #Doesn't make cpg
CC<-subset(cpg.y, WTnt=='c') #Does make cpg
CNC<-subset(cpg.n, WTnt=='c') #Doesn't make cpg


#Function to help create errorbars
sem<-function(x){
  return(sd(x,na.rm = FALSE)/sqrt(length(x)))
}

#making the data frames with all information about a, t, c, g 
AllA = rbind(AC,ANC)
AllT = rbind(TC,TNC)
AllC = rbind(CC,CNC)
AllG = rbind(GC,GNC)

AllA$meanofmeans_value<-0
AllT$meanofmeans_value<-0
AllC$meanofmeans_value<-0
AllG$meanofmeans_value<-0

# for loops to caculate mean ans errorbars and 1, 2, 3, 4 for position
for (i in 1:length(AllA$makesCpG)) {
  if (AllA$makesCpG[i] == 1 && AllA$TypeOfSite[i] == "syn") {
    AllA$graphit[i] <- 2
  }
  if (AllA$makesCpG[i] == 0 && AllA$TypeOfSite[i] == "syn") {
    AllA$graphit[i] <- 1
  }
}

AllA[(which(AllA$makesCpG == 1 & AllA$TypeOfSite == "syn") ),13] <- median(AllA$colMeansTsZanini[(which(AllA$makesCpG == 1 & AllA$TypeOfSite == "syn") )])
AllA[(which(AllA$makesCpG == 0 & AllA$TypeOfSite == "syn") ),13] <- median(AllA$colMeansTsZanini[(which(AllA$makesCpG == 0 & AllA$TypeOfSite == "syn") )])

for (i in 1:length(AllT$makesCpG)) {
  if (AllT$makesCpG[i] == 1 && AllT$TypeOfSite[i] == "syn") {
    AllT$graphit[i] <- 2
  }
  if (AllT$makesCpG[i] == 0 && AllT$TypeOfSite[i] == "syn") {
    AllT$graphit[i] <- 1
  }
}
AllT[(which(AllT$makesCpG == 1 & AllT$TypeOfSite == "syn") ),13] <- median(AllT$colMeansTsZanini[(which(AllT$makesCpG == 1 & AllT$TypeOfSite == "syn") )])
AllT[(which(AllT$makesCpG == 0 & AllT$TypeOfSite == "syn") ),13] <- median(AllT$colMeansTsZanini[(which(AllT$makesCpG == 0 & AllT$TypeOfSite == "syn") )])

for (i in 1:length(AllC$makesCpG)) {
  if (AllC$makesCpG[i] == 1 && AllC$TypeOfSite[i] == "syn") {
    AllC$graphit[i] <- 2
  }
  if (AllC$makesCpG[i] == 0 && AllC$TypeOfSite[i] == "syn") {
    AllC$graphit[i] <- 1
  }
}
AllC[(which(AllC$makesCpG == 1 & AllC$TypeOfSite == "syn") ),13] <- median(AllC$colMeansTsZanini[(which(AllC$makesCpG == 1 & AllC$TypeOfSite == "syn") )])
AllC[(which(AllC$makesCpG == 0 & AllC$TypeOfSite == "syn") ),13] <- median(AllC$colMeansTsZanini[(which(AllC$makesCpG == 0 & AllC$TypeOfSite == "syn") )])

for (i in 1:length(AllG$makesCpG)) {
  if (AllG$makesCpG[i] == 1 && AllG$TypeOfSite[i] == "syn") {
    AllG$graphit[i] <- 2
  }
  if (AllG$makesCpG[i] == 0 && AllG$TypeOfSite[i] == "syn") {
    AllG$graphit[i] <- 1
  }
}
AllG[(which(AllG$makesCpG == 1 & AllG$TypeOfSite == "syn") ),13] <- median(AllG$colMeansTsZanini[(which(AllG$makesCpG == 1 & AllG$TypeOfSite == "syn") )])
AllG[(which(AllG$makesCpG == 0 & AllG$TypeOfSite == "syn") ),13] <- median(AllG$colMeansTsZanini[(which(AllG$makesCpG == 0 & AllG$TypeOfSite == "syn") )])





#####################################################################

truenamepng = paste("Output/CpG_vs_NonCpG.png",sep="")
png(truenamepng, width = 6.75, height = 6.75, units = "in", res= 300)
par(mfrow=c(2,2), mar=c(4.1, 4.1, 1.9, 0.8),oma=c(0.1,0.1,1.5,0.1)) 
palette(alpha(c("#99FF99","#9999FF","#FF9900","#FF3300"),0.3))

plot(AllA$graphit, AllA$meanofmeans_value,log='y',col=factor(AllA$graphit),pch=19,cex = 3, main="A->G",xlab = " ", ylab = "Mutation Frequency", yaxt="n", xaxt="n", ylim=c(0.0001, 0.5),xlim=c(.5,2.5))
#plot all points used from median function onto 1 & 2
#make table for p-values
p=AllA$colMeansTsZanini[AllA$makesCpG == 0 & AllA$TypeOfSite == "syn"]
p1<-rep(1,length(p))
p1<-rnorm(n = length(p), mean = 1, sd = 0.05)
q=AllA$colMeansTsZanini[AllA$makesCpG == 1 & AllA$TypeOfSite == "syn"] 
q1<-rnorm(n = length(q), mean = 2, sd = 0.05)
points(p1,p,col="black")
points(q1,q,col="black")

eaxis(2,at=c(10^-3,10^-2,10^-1,10^0),cex.axis=1.1)
axis(1, at= c(1:4),labels = c("No CpG \n Non-Coding", " CpG \n Non-Coding", "No CpG \n NonSyn", "CpG \n NonSyn"), mgp=c(3, 1.5, 0))
axis.break(2,0.0002*(1+0.02),breakcol="black",style="slash")
mtext('0', side=2, line=1.5, at=0.0001, las=1.1)
mtext("5'UTR", outer=TRUE, adj=0.55, cex=1.7, line=0.01)



plot(AllT$graphit, AllT$meanofmeans_value,log='y',col=factor(AllT$graphit),pch=19,cex = 3,main="T->C",xlab = " ", ylab = "Mutation Frequency", yaxt="n", xaxt = "n", ylim=c(0.0001, .5), xlim=c(.5,2.4))
r=AllT$colMeansTsZanini[AllT$makesCpG == 0 & AllT$TypeOfSite == "syn"]
r1<-rnorm(n = length(r), mean = 1, sd = 0.05)
s=AllT$colMeansTsZanini[AllT$makesCpG == 1 & AllT$TypeOfSite == "syn"] 
s1<-rnorm(n = length(s), mean = 2, sd = 0.05)
points(r1,r,col="black")
points(s1,s,col="black")

#plot all points used from median function onto 1 & 2
#AllT$colMeansTsZanini[(which(AllT$makesCpG == 1 & AllT$TypeOfSite == "syn") )] on 1
#AllT$colMeansTsZanini[(which(AllT$makesCpG == 0 & AllT$TypeOfSite == "syn") )] on 2

eaxis(2,at=c(10^-3,10^-2,10^-1,10^0),cex.axis=1.1)
axis(1, at= c(1:4),labels = c("No CpG \n Non-Coding", " CpG \n Non-Coding", "No CpG \n NonSyn", "CpG \n NonSyn"),mgp=c(3, 1.5, 0))
axis.break(2,0.0002*(1+0.02),breakcol="black",style="slash")
mtext('0', side=2, line=1.5, at=0.0001, las=1.1)




#

plot(AllC$graphit, AllC$meanofmeans_value,log='y',col=factor(AllC$graphit),pch=19,cex = 3,main="G->C",xlab = " ", ylab = "Mutation Frequency", yaxt="n", xaxt = "n", ylim=c(0.0001, .5), xlim=c(.5,2.4))
q=AllC$colMeansTsZanini[AllC$makesCpG == 0 & AllC$TypeOfSite == "syn"]
q1<-rnorm(n = length(q), mean = 1, sd = 0.05)
w=AllC$colMeansTsZanini[AllC$makesCpG == 1 & AllC$TypeOfSite == "syn"] 
w1<-rnorm(n = length(w), mean = 2, sd = 0.05)
points(r1,r,col="black")
points(s1,s,col="black")

#plot all points used from median function onto 1 & 2
#AllT$colMeansTsZanini[(which(AllT$makesCpG == 1 & AllT$TypeOfSite == "syn") )] on 1
#AllT$colMeansTsZanini[(which(AllT$makesCpG == 0 & AllT$TypeOfSite == "syn") )] on 2

eaxis(2,at=c(10^-3,10^-2,10^-1,10^0),cex.axis=1.1)
axis(1, at= c(1:4),labels = c("No CpG \n Non-Coding", " CpG \n Non-Coding", "No CpG \n NonSyn", "CpG \n NonSyn"),mgp=c(3, 1.5, 0))
axis.break(2,0.0002*(1+0.02),breakcol="black",style="slash")
mtext('0', side=2, line=1.5, at=0.0001, las=1.1)




plot(AllG$graphit, AllG$meanofmeans_value,log='y',col=factor(AllG$graphit),pch=19,cex = 3,main="C->G",xlab = " ", ylab = "Mutation Frequency", yaxt="n", xaxt = "n", ylim=c(0.0001, .5), xlim=c(.5,2.4))
f=AllG$colMeansTsZanini[AllG$makesCpG == 0 & AllG$TypeOfSite == "syn"]
f1<-rnorm(n = length(f), mean = 1, sd = 0.05)
g=AllG$colMeansTsZanini[AllG$makesCpG == 1 & AllG$TypeOfSite == "syn"] 
g1<-rnorm(n = length(g), mean = 2, sd = 0.05)
points(r1,r,col="black")
points(s1,s,col="black")

#plot all points used from median function onto 1 & 2
#AllT$colMeansTsZanini[(which(AllT$makesCpG == 1 & AllT$TypeOfSite == "syn") )] on 1
#AllT$colMeansTsZanini[(which(AllT$makesCpG == 0 & AllT$TypeOfSite == "syn") )] on 2

eaxis(2,at=c(10^-3,10^-2,10^-1,10^0),cex.axis=1.1)
axis(1, at= c(1:4),labels = c("No CpG \n Non-Coding", " CpG \n Non-Coding", "No CpG \n NonSyn", "CpG \n NonSyn"),mgp=c(3, 1.5, 0))
axis.break(2,0.0002*(1+0.02),breakcol="black",style="slash")
mtext('0', side=2, line=1.5, at=0.0001, las=1.1)



#

palette(alpha(c("#99FF99","#FF9900"),0.3))
dev.off()




#This is the wilcox file that needs to become its own file (just lazy)
truenamepng1 = paste("Output/P-ValueTable.png",sep="")
png(truenamepng1, width = 6.75, height = 6.75, units = "in", res= 300)




wilcox.test(q,p,alternative='less')

wilcox.test(r,s,alternative='less')

dev.off()



