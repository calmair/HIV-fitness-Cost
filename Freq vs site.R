library(ape)
library(seqinr)
library(pegas)
library(sfsmisc)
library(ggplot2)
library(scales)
library(plotrix)
library(RColorBrewer)


OverviewDFZanini<-read.table("Output/OverviewSelCoeffZanini.csv",sep=",",header=TRUE,row.names=1)
for(n in 1:nrow(OverviewDFZanini)){
OverviewDFZanini[n,1]<-OverviewDFZanini[n,1]+550
}

truenamepng = paste("Output/FreqvsSite.png",sep="")
png(truenamepng, width = 6.75, height = 6.75, units = "in", res= 300)
#edit to give each nucliec acid a different color
OverviewDFZanini$Colour="#8b9dff"
OverviewDFZanini$Colour[OverviewDFZanini$WTnt=="t"]="#ffdf8f"
OverviewDFZanini$Colour[OverviewDFZanini$WTnt=="c"]="#ff83a2"
OverviewDFZanini$Colour[OverviewDFZanini$WTnt=="g"]="#74ff9b"





# Set new column values to appropriate colours

# Plot all points at once, using newly generated colours
plot(OverviewDFZanini$num,OverviewDFZanini$colMeansTsZanini, col=OverviewDFZanini$Colour,log="y",xlab="Sites",ylab="Estimated Frequencies",main="5UTR Region",pch=19)


dev.off()
