library(ape)
library(seqinr)
library(pegas)
library(sfsmisc)
library(ggplot2)
library(scales)
library(plotrix)
library(RColorBrewer)
#Makes graphs for T->C and A->G
#need to make this into a function to make a .png for each gene

OverviewDFZanini<-read.table("Output/OverviewSelCoeffZanini.csv",sep=",",header=TRUE,row.names=1)
dim(OverviewDFZanini)
#548 is when data is not zero
#9591 is when zero's appear again

OverviewDFZanini1<- OverviewDFZanini[548:634,1:10] #5LTR
OverviewDFZanini2<- OverviewDFZanini[790:2292,1:10] #Gag
OverviewDFZanini3<- OverviewDFZanini[2085:5096,1:10] #pol
OverviewDFZanini4<- OverviewDFZanini[5041:6062,1:10] #vif, vpr, vpu
OverviewDFZanini5<- OverviewDFZanini[6225:8795,1:10] #env
OverviewDFZanini6<- OverviewDFZanini[8797:9417,1:10] #nef
OverviewDFZanini7<- OverviewDFZanini[9086:9591,1:10] #3LTR

for (j in 1:7){
print(j)
cpg.y<-subset(get(paste("OverviewDFZanini",j,sep = "")), makesCpG==1)
#print(length(row.names(cpg.y)))
cpg.n<-subset(get(paste("OverviewDFZanini",j,sep = "")), makesCpG==0)
#print(length(row.names(cpg.n)))


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
AllC = rbind(CNC)
AllG = rbind(GNC)

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

AllA[(which(AllA$makesCpG == 1 & AllA$TypeOfSite == "syn") ),13] <- mean(AllA$colMeansTsZanini[(which(AllA$makesCpG == 1 & AllA$TypeOfSite == "syn") )])
AllA[(which(AllA$makesCpG == 0 & AllA$TypeOfSite == "syn") ),13] <- mean(AllA$colMeansTsZanini[(which(AllA$makesCpG == 0 & AllA$TypeOfSite == "syn") )])

print(mean(AllA$colMeansTsZanini[(which(AllA$makesCpG == 1 & AllA$TypeOfSite == "syn") )]))
print(mean(AllA$colMeansTsZanini[(which(AllA$makesCpG == 0 & AllA$TypeOfSite == "syn") )]))

for (i in 1:length(AllT$makesCpG)) {
  if (AllT$makesCpG[i] == 1 && AllT$TypeOfSite[i] == "syn") {
    AllT$graphit[i] <- 2
  }
  if (AllT$makesCpG[i] == 0 && AllT$TypeOfSite[i] == "syn") {
    AllT$graphit[i] <- 1
  }
}
AllT[(which(AllT$makesCpG == 1 & AllT$TypeOfSite == "syn") ),13] <- mean(AllT$colMeansTsZanini[(which(AllT$makesCpG == 1 & AllT$TypeOfSite == "syn") )])
AllT[(which(AllT$makesCpG == 0 & AllT$TypeOfSite == "syn") ),13] <- mean(AllT$colMeansTsZanini[(which(AllT$makesCpG == 0 & AllT$TypeOfSite == "syn") )])

print(mean(AllT$colMeansTsZanini[(which(AllT$makesCpG == 1 & AllT$TypeOfSite == "syn") )]))
print(mean(AllT$colMeansTsZanini[(which(AllT$makesCpG == 0 & AllT$TypeOfSite == "syn") )]))

for (i in 1:length(AllC$makesCpG)) {
  if (AllC$makesCpG[i] == 1 && AllC$TypeOfSite[i] == "syn") {
    AllC$graphit[i] <- 2
  }
  if (AllC$makesCpG[i] == 0 && AllC$TypeOfSite[i] == "syn") {
    AllC$graphit[i] <- 1
  }
}
AllC[(which(AllC$makesCpG == 1 & AllC$TypeOfSite == "syn") ),13] <- mean(AllC$colMeansTsZanini[(which(AllC$makesCpG == 1 & AllC$TypeOfSite == "syn") )])
AllC[(which(AllC$makesCpG == 0 & AllC$TypeOfSite == "syn") ),13] <- mean(AllC$colMeansTsZanini[(which(AllC$makesCpG == 0 & AllC$TypeOfSite == "syn") )])

for (i in 1:length(AllG$makesCpG)) {
  if (AllG$makesCpG[i] == 1 && AllG$TypeOfSite[i] == "syn") {
    AllG$graphit[i] <- 2
  }
  if (AllG$makesCpG[i] == 0 && AllG$TypeOfSite[i] == "syn") {
    AllG$graphit[i] <- 1
  }
}
AllG[(which(AllG$makesCpG == 1 & AllG$TypeOfSite == "syn") ),13] <- mean(AllG$colMeansTsZanini[(which(AllG$makesCpG == 1 & AllG$TypeOfSite == "syn") )])
AllG[(which(AllG$makesCpG == 0 & AllG$TypeOfSite == "syn") ),13] <- mean(AllG$colMeansTsZanini[(which(AllG$makesCpG == 0 & AllG$TypeOfSite == "syn") )])





#####################################################################

truenamepng = paste("Output/CpG_vs_NonCpG",j,".png",sep="")
png(truenamepng, width = 6.75, height = 6.75, units = "in", res= 300)
par(mfrow=c(2,2), mar=c(4.1, 4.1, 1.9, 0.8),oma=c(0.1,0.1,1.5,0.1)) 
palette(alpha(c("#2E4C2E","#38385D","#FF9900","#FF3300"),0.3))

plot(AllA$graphit, AllA$meanofmeans_value,log='y',col=factor(AllA$graphit),pch=19,cex = 3, main="A->G",xlab = " ", ylab = "Mutation Frequency", yaxt="n", xaxt="n", ylim=c(0.0001, 0.5),xlim=c(.5,2.5))
#plot all points used from mean function onto 1 & 2
#make table for p-values
p=AllA$colMeansTsZanini[AllA$makesCpG == 0 & AllA$TypeOfSite == "syn"]
p1<-rep(1,length(p))
p1<-rnorm(n = length(p), mean = 1, sd = 0.25)
q=AllA$colMeansTsZanini[AllA$makesCpG == 1 & AllA$TypeOfSite == "syn"] 
q1<-rnorm(n = length(q), mean = 2, sd = 0.25)
points(p1,p,col= alpha("#97CC04",1))
points(q1,q,col= alpha("#7E7ED1",1))

eaxis(2,at=c(10^-3,10^-2,10^-1,10^0),cex.axis=1.1)
axis(1, at= c(1:4),labels = c("No CpG ", " CpG ", "No CpG \n NonSyn", "CpG \n NonSyn"), mgp=c(3, 1.5, 0))
axis.break(2,0.0002*(1+0.02),breakcol="black",style="slash")
mtext('0', side=2, line=1.5, at=0.0001, las=1.1)
mtext("5'UTR", outer=TRUE, adj=0.55, cex=1.7, line=0.01)



plot(AllT$graphit, AllT$meanofmeans_value,log='y',col=factor(AllT$graphit),pch=19,cex = 3,main="T->C",xlab = " ", ylab = "Mutation Frequency", yaxt="n", xaxt = "n", ylim=c(0.0001, .5), xlim=c(.5,2.4))
r=AllT$colMeansTsZanini[AllT$makesCpG == 0 & AllT$TypeOfSite == "syn"]
r1<-rnorm(n = length(r), mean = 1, sd = 0.25)
s=AllT$colMeansTsZanini[AllT$makesCpG == 1 & AllT$TypeOfSite == "syn"] 
s1<-rnorm(n = length(s), mean = 2, sd = 0.25)
points(r1,r,col= alpha("#97CC04",1))
points(s1,s,col= alpha("#7E7ED1",1))

#plot all points used from mean function onto 1 & 2
#AllT$colMeansTsZanini[(which(AllT$makesCpG == 1 & AllT$TypeOfSite == "syn") )] on 1
#AllT$colMeansTsZanini[(which(AllT$makesCpG == 0 & AllT$TypeOfSite == "syn") )] on 2

eaxis(2,at=c(10^-3,10^-2,10^-1,10^0),cex.axis=1.1)
axis(1, at= c(1:4),labels = c("No CpG ", " CpG ", "No CpG \n NonSyn", "CpG \n NonSyn"),mgp=c(3, 1.5, 0))
axis.break(2,0.0002*(1+0.02),breakcol="black",style="slash")
mtext('0', side=2, line=1.5, at=0.0001, las=1.1)




#

plot(AllG$graphit, AllG$meanofmeans_value,log='y',col=factor(AllG$graphit),pch=19,cex = 3,main="G->A",xlab = " ", ylab = "Mutation Frequency", yaxt="n", xaxt = "n", ylim=c(0.0001, .5), xlim=c(.5,2.4))
f=AllG$colMeansTsZanini[AllG$makesCpG == 0 & AllG$TypeOfSite == "syn"]
f1<-rnorm(n = length(f), mean = 1, sd = 0.25)
g=AllG$colMeansTsZanini[AllG$makesCpG == 1 & AllG$TypeOfSite == "syn"] 
g1<-rnorm(n = length(g), mean = 2, sd = 0.25)
points(f1,f,col= alpha("#97CC04",1))
points(g1,g,col= alpha("#7E7ED1",1))

#plot all points used from mean function onto 1 & 2
#AllT$colMeansTsZanini[(which(AllT$makesCpG == 1 & AllT$TypeOfSite == "syn") )] on 1
#AllT$colMeansTsZanini[(which(AllT$makesCpG == 0 & AllT$TypeOfSite == "syn") )] on 2

eaxis(2,at=c(10^-3,10^-2,10^-1,10^0),cex.axis=1.1)
axis(1, at= c(1:4),labels = c("No CpG ", " CpG ", "No CpG \n NonSyn", "CpG \n NonSyn"),mgp=c(3, 1.5, 0))
axis.break(2,0.0002*(1+0.02),breakcol="black",style="slash")
mtext('0', side=2, line=1.5, at=0.0001, las=1.1)



plot(AllC$graphit, AllC$meanofmeans_value,log='y',col=factor(AllC$graphit),pch=19,cex = 3,main="C->T",xlab = " ", ylab = "Mutation Frequency", yaxt="n", xaxt = "n", ylim=c(0.0001, .5), xlim=c(.5,2.4))
q=AllC$colMeansTsZanini[AllC$makesCpG == 0 & AllC$TypeOfSite == "syn"]
q1<-rnorm(n = length(q), mean = 1, sd = 0.25)
w=AllC$colMeansTsZanini[AllC$makesCpG == 1 & AllC$TypeOfSite == "syn"] 
w1<-rnorm(n = length(w), mean = 2, sd = 0.25)
points(q1,q,col= alpha("#97CC04",1))
points(w1,w,col= alpha("#7E7ED1",1))

#plot all points used from mean function onto 1 & 2
#AllT$colMeansTsZanini[(which(AllT$makesCpG == 1 & AllT$TypeOfSite == "syn") )] on 1
#AllT$colMeansTsZanini[(which(AllT$makesCpG == 0 & AllT$TypeOfSite == "syn") )] on 2

eaxis(2,at=c(10^-3,10^-2,10^-1,10^0),cex.axis=1.1)
axis(1, at= c(1:4),labels = c("No CpG ", " CpG ", "No CpG \n NonSyn", "CpG \n NonSyn"),mgp=c(3, 1.5, 0))
axis.break(2,0.0002*(1+0.02),breakcol="black",style="slash")
mtext('0', side=2, line=1.5, at=0.0001, las=1.1)

#

palette(alpha(c("#99FF99","#FF9900"),0.3))
dev.off()
}


