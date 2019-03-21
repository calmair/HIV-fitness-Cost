  library(graphics)
  library(plyr)
  library(dplyr)
  
  data<-read.table("Output/OverviewSelCoeffZanini.csv",sep=",",header=TRUE,row.names=1)
  #This function is going to read the data from the csv files 
  
  
  pVals = c()
  shrtval = 0
  options(scipen=999)
  #prevents pvalues from becoming scientific notation. 
  
  array1 = data$colMeansTsZanini[data$WTnt =="a" & data$TypeOfSite == 'syn' & data$makesCpG == 1]
  array2 = data$colMeansTsZanini[data$WTnt =="a" & data$TypeOfSite == 'syn' & data$makesCpG == 0]
  CpGa = data$colMeansTsZanini[data$WTnt =="a" & data$makesCpG == 1]
  nonCpGa = data$colMeansTsZanini[data$WTnt =="a" &  data$makesCpG == 0]
  
  print("For a: Comparing makes CpG with noCpG (syn). Wilcox test less: red/blue")
  print(wilcox.test(array1, array2, alternative='less'))
  pVals = c(pVals,format(wilcox.test(array1, array2, alternative='less')$p.value, nsmall = 6))
  print(pVals)
  
  
  
  
  array5 = data$colMeansTsZanini[data$WTnt =="t" & data$TypeOfSite == 'syn' & data$makesCpG == 1]
  print(array5)
  array6 = data$colMeansTsZanini[data$WTnt =="t" & data$TypeOfSite == 'syn' & data$makesCpG == 0]
  print(array6)
  CpGt = data$colMeansTsZanini[data$WTnt =="t" & data$makesCpG == 1]
  nonCpGt = data$colMeansTsZanini[data$WTnt =="t" &  data$makesCpG == 0]
  
  print("For t: Comparing makes CpG with noCpG (syn). Wilcox test less: red/blue")
  print(wilcox.test(array5, array6, alternative='less'))
  pVals = c(pVals,format(wilcox.test(array5, array6, alternative='less')$p.value, nsmall = 6))
  
  Pvalues= c(pVals)
  #save Pvalues into list
print(Pvalues)


  options(scipen = 999)
  #setwd("output/redeploy/")
  truenamepdf= paste("Output/","P-values",".pdf",sep="")
  truenamepng= paste("P-values","tables", ".png", sep="")
  #print(truenamepdf)
  #prevents pvalues from becoming scientific notation
  options(warn=-1)
  #suppress warnings
  
  #setwd("~/Desktop/Something_Cool-CpG_Sites-/Tables")
  #table construct
  pdf(truenamepdf, width = 7, height= 5)
  #png(truenamepng, width = 6.75, height = 6.75, units = "in", res= 300)
  col1 <- c("A->G", "T->C")
  col2 <- c("Syn: CpG v NonCpG")
  ycoor <- c(4*100/5+.7 , 3*100/5 + 5.1, 3*100/5 - 10, 2*100/5 -6.4, 1*100/5-1.3, 100/5-14- 2.9, 100)
  ycoorb <- c(4*100/5+.7 , 3*100/5 + 5.1, 3*100/5 - 10.6, 2*100/5 -6.4, 1*100/5-1.3, 100/5-14- 2.9, 100)
  df = data.frame(col1, col2, Pvalues)
  
  par(xpd=F)
  plot(1, 2, xlim=c(0,100),ylim=c(0,100), col=0, xaxt="n", yaxt="n", xlab="", ylab="")
  title(main = "CpGvsNoCpG", family = "Times", adj = 0.5, cex.main= 2)
  abline(v = 100/5)
  abline(v = 2*100/3)
  abline(h = 100-100/7 + 3)
  abline(h = 42)

  
  
  text(x=100/7- 6, y= 5*100/5-3, "Mutation Type")
  text(x=3*100/7, y = 5*100/5-3, "Comparison")
  text(x= 6*100/7, y=5*100/5-3, "P-Value")
  rect(xleft = -4, xright = 100/5, ybottom =42, ytop =100-100/7+3 , col = "white")
  text(x= 100/12, y= 3*100/5+5, "A->G", cex = 1.7, family = "Times")
  rect(xleft = -4, xright = 100/5, ybottom =-4, ytop =42 , col = "white")
  text(x= 100/12, y = 1*100/5 - 1, "T->C", cex = 1.7, family ='Times')
  
  text(x= 3*100/7, y =  3*100/5+5, labels= col2[1])#cpgvsnoncpg

  text(x= 3*100/7, y =  1*100/5, labels= col2[1])#cpgvsnoncpg

    text(x= 6*100/7, y =3*100/5+5, labels = Pvalues[1],cex=1.2)
    text(x= 6*100/7, y =1*100/5, labels = Pvalues[2],cex = 1.2)
    
    
  print("end")
  #dev.copy(pdf, truenamepng)
  dev.off()

