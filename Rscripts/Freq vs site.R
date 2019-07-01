library(ape)
library(seqinr)
library(pegas)
library(sfsmisc)
library(ggplot2)
library(scales)
library(plotrix)
library(RColorBrewer)


OverviewDFZanini<-read.table("Output/OverviewSelCoeffZanini.csv",sep=",",header=TRUE,row.names=1)

truenamepng = paste("Output/FreqvsSite.png",sep="")
png(truenamepng, width = 30, height = 7.5, units = "in",pointsize = 5, res= 1000)
#edit to give each nucliec acid a different color
OverviewDFZanini$Colour="#8b9dff"
OverviewDFZanini$Colour[OverviewDFZanini$WTnt=="t"]="#ffdf8f"
OverviewDFZanini$Colour[OverviewDFZanini$WTnt=="c"]="#ff83a2"
OverviewDFZanini$Colour[OverviewDFZanini$WTnt=="g"]="#74ff9b"





# Set new column values to appropriate colours

# Plot all points at once, using newly generated colours
# add section that illustrates genome underneath xlab
plot(OverviewDFZanini$num,OverviewDFZanini$colMeansTsZanini, col=OverviewDFZanini$Colour,log="y",xlab="Sites",ylab="Estimated Frequencies",main="Genome-Wide",pch=19, xaxt  = "n")
axis(1,at = c(0,500,1000,1500,2000,2500,3000,3500,4000,4500,5000,5500,6000,6500,7000,7500,8000,8500,9000,9500,9719), labels = c(0,500,1000,1500,2000,2500,3000,3500,4000,4500,5000,5500,6000,6500,7000,7500,8000,8500,9000,9500,9719))

dev.off()
