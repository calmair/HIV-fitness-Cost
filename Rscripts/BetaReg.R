#Prep Data for GLM analysis
#need CountData & OverviewDF

#USE BTR instead of GLM
library(betareg)
library(tidyverse)
require(xtable)

TMFD<-read.csv("./Output/TransFreqGLM.csv")
print("done4")

TMFD<-na.omit(TMFD) #Removes rows that contain NA's found from Shape Column
#table(is.na(TMFD$TransMutFreq))
###Edit Estimate code block
#See notebook for corrections


############################################################################################################
#tempsample<-TMFD[TMFD$Res == 0 & TMFD$stop == 0 & TMFD$Section1 ==1 & TMFD$TransMutFreq >0,] #used to test if all rows are correctly called for

#Section 1 Analysis 5utr
# fullmodel.int <- betareg(TransMutFreq ~ shape + t + c + g + CpG + CpG*t  + nonsyn + t*nonsyn + c*nonsyn + g*nonsyn + nonsyn*CpG + t:nonsyn:CpG + bigAAchange, data = TMFD[TMFD$Res == 0 & TMFD$stop == 0 & TMFD$Section1 ==1 & TMFD$TransMutFreq >0,])
# sumOfModel <- summary(fullmodel.int)$coeff
# print(sumOfModel)
# OverallOfModel<-xtable(sumOfModel$mean, digits = 3)
# OverallOfModel$Effect<-NA
# for (g in 1:length(row.names(OverallOfModel)) ){
#   print(g)
# if (g==1){
#   OverallOfModel$Effect[1]<- round(exp(OverallOfModel[1,g]), digits = 3)
# }
# else{
# OverallOfModel$Effect[g]<- round((((exp(OverallOfModel[1,1] + OverallOfModel$Estimate[g]) - exp(OverallOfModel[1,1])) /exp(OverallOfModel[1,1]))*100)) #add estimate % column
# }
# }
# print(xtable(OverallOfModel, digits = 3),type="html",file="Output/OverallOfGLM5UTR.html")

#Section 2 Analysis GAG
fullmodel.int <- betareg(TransMutFreq ~ shape + t + c + g + CpG + CpG*t  + t*nonsyn + c*nonsyn + g*nonsyn + nonsyn*CpG + t:nonsyn:CpG + bigAAchange, data = TMFD[TMFD$Res == 0 & TMFD$stop == 0 & TMFD$Section2 ==1 & TMFD$TransMutFreq >0,]) # Not returning the correct positions from gene section 2
sumOfModel <- summary(fullmodel.int)$coeff
print(sumOfModel)
OverallOfModel<-xtable(sumOfModel$mean, digits = 3)
OverallOfModel$Effect<-NA
for (g in 1:length(row.names(OverallOfModel)) ){
  print(g)
  if (g==1){
    OverallOfModel$Effect[1]<- round(exp(OverallOfModel[1,g]), digits = 3)
  }
  else{
    OverallOfModel$Effect[g]<- round((((exp(OverallOfModel[1,1] + OverallOfModel$Estimate[g]) - exp(OverallOfModel[1,1])) /exp(OverallOfModel[1,1]))*100)) #add estimate % column
  }
}
print(xtable(OverallOfModel, digits = 3),type="html",file="Output/OverallOfGLMGAG.html")

#Section 3 POL
fullmodel.int <- betareg(TransMutFreq ~ shape + t + c + g + CpG + CpG*t  + t*nonsyn + c*nonsyn + g*nonsyn + nonsyn*CpG + t:nonsyn:CpG + bigAAchange,   data = TMFD[TMFD$Res == 0 & TMFD$stop == 0 & TMFD$Section3 ==1 & TMFD$TransMutFreq >0,])
#fullmodel.int <- betareg(TransMutFreq ~ nonsyn,   data = TMFD[TMFD$Res == 0 & TMFD$stop == 0 & TMFD$Section3 ==1 & TMFD$TransMutFreq >0,])
sumOfModel <- summary(fullmodel.int)$coeff
print(sumOfModel)
OverallOfModel<-xtable(sumOfModel$mean, digits = 3)
OverallOfModel$Effect<-NA
for (g in 1:length(row.names(OverallOfModel)) ){
  print(g)
  if (g==1){
    OverallOfModel$Effect[1]<- round(exp(OverallOfModel[1,g]), digits = 3)
  }
  else{
    OverallOfModel$Effect[g]<- round((((exp(OverallOfModel[1,1] + OverallOfModel$Estimate[g]) - exp(OverallOfModel[1,1])) /exp(OverallOfModel[1,1]))*100)) #add estimate % column
  }
}
print(xtable(OverallOfModel, digits = 3),type="html",file="Output/OverallOfGLMPol.html")

#Section 4 Analysis VIF
fullmodel.int <- betareg(TransMutFreq ~ shape + t + c + g + CpG + CpG*t  + t*nonsyn + c*nonsyn + g*nonsyn + nonsyn*CpG + t:nonsyn:CpG + bigAAchange,   data = TMFD[TMFD$Res == 0 & TMFD$stop == 0 & TMFD$Section4 ==1 & TMFD$TransMutFreq >0,])
sumOfModel <- summary(fullmodel.int)$coeff
print(sumOfModel)
OverallOfModel<-xtable(sumOfModel$mean, digits = 3)
OverallOfModel$Effect<-NA
for (g in 1:length(row.names(OverallOfModel)) ){
  print(g)
  if (g==1){
    OverallOfModel$Effect[1]<- round(exp(OverallOfModel[1,g]), digits = 3)
  }
  else{
    OverallOfModel$Effect[g]<- round((((exp(OverallOfModel[1,1] + OverallOfModel$Estimate[g]) - exp(OverallOfModel[1,1])) /exp(OverallOfModel[1,1]))*100)) #add estimate % column
  }
}
print(xtable(OverallOfModel, digits = 3),type="html",file="Output/OverallOfGLMVIF.html")


#Section 5 Analysis VPR
fullmodel.int <- betareg(TransMutFreq ~ shape + t + c + g + CpG + CpG*t  + t*nonsyn + c*nonsyn + g*nonsyn + nonsyn*CpG + t:nonsyn:CpG + bigAAchange,   data = TMFD[TMFD$Res == 0 & TMFD$stop == 0 & TMFD$Section5 ==1 & TMFD$TransMutFreq >0,])
sumOfModel <- summary(fullmodel.int)$coeff
print(sumOfModel)
OverallOfModel<-xtable(sumOfModel$mean, digits = 3)
OverallOfModel$Effect<-NA
for (g in 1:length(row.names(OverallOfModel)) ){
  print(g)
  if (g==1){
    OverallOfModel$Effect[1]<- round(exp(OverallOfModel[1,g]), digits = 3)
  }
  else{
    OverallOfModel$Effect[g]<- round((((exp(OverallOfModel[1,1] + OverallOfModel$Estimate[g]) - exp(OverallOfModel[1,1])) /exp(OverallOfModel[1,1]))*100)) #add estimate % column
  }
}
print(xtable(OverallOfModel, digits = 3),type="html",file="Output/OverallOfGLMVPR.html")

#Section 6 Analysis ENV
fullmodel.int <- betareg(TransMutFreq ~ shape + t + c + g + CpG + CpG*t  + t*nonsyn + c*nonsyn + g*nonsyn + nonsyn*CpG + t:nonsyn:CpG + bigAAchange,   data = TMFD[TMFD$Res == 0 & TMFD$stop == 0 & TMFD$Section6 ==1 & TMFD$TransMutFreq >0,])
sumOfModel <- summary(fullmodel.int)$coeff
print(sumOfModel)
OverallOfModel<-xtable(sumOfModel$mean, digits = 3)
OverallOfModel$Effect<-NA
for (g in 1:length(row.names(OverallOfModel)) ){
  print(g)
  if (g==1){
    OverallOfModel$Effect[1]<- round(exp(OverallOfModel[1,g]), digits = 3)
  }
  else{
    OverallOfModel$Effect[g]<- round((((exp(OverallOfModel[1,1] + OverallOfModel$Estimate[g]) - exp(OverallOfModel[1,1])) /exp(OverallOfModel[1,1]))*100)) #add estimate % column
  }
}
print(xtable(OverallOfModel, digits = 3),type="html",file="Output/OverallOfGLMENV.html")

#Section 7 Analysis NEF
fullmodel.int <- betareg(TransMutFreq ~ shape + t + c + g + CpG + CpG*t  + t*nonsyn + c*nonsyn + g*nonsyn + nonsyn*CpG + t:nonsyn:CpG + bigAAchange,   data = TMFD[TMFD$Res == 0 & TMFD$stop == 0 & TMFD$Section7 ==1 & TMFD$TransMutFreq >0,])
sumOfModel <- summary(fullmodel.int)$coeff
print(sumOfModel)
OverallOfModel<-xtable(sumOfModel$mean, digits = 3)
OverallOfModel$Effect<-NA
for (g in 1:length(row.names(OverallOfModel)) ){
  print(g)
  if (g==1){
    OverallOfModel$Effect[1]<- round(exp(OverallOfModel[1,g]), digits = 3)
  }
  else{
    OverallOfModel$Effect[g]<- round((((exp(OverallOfModel[1,1] + OverallOfModel$Estimate[g]) - exp(OverallOfModel[1,1])) /exp(OverallOfModel[1,1]))*100)) #add estimate % column
  }
}
print(xtable(OverallOfModel, digits = 3),type="html",file="Output/OverallOfGLMNEF.html")

#Overall Analysis
#Need More RAM to Run LMAO
names(TMFD)[16:22]<-c("FiveUTR","Gag", "Pol", "Vif", "Vpr", "Env", "Nef")
fullmodel.int <- betareg(TransMutFreq ~ shape + t + c + g + CpG + CpG*t  + t*nonsyn + c*nonsyn + g*nonsyn + nonsyn*CpG + t:nonsyn:CpG + bigAAchange + FiveUTR + Gag + Pol + Vif + Vpr + Env + Nef,   data = TMFD[TMFD$Res == 0 & TMFD$stop == 0 & TMFD$TransMutFreq >0,])
sumOfModel <- summary(fullmodel.int)$coeff
print(sumOfModel)
OverallOfModel<-xtable(sumOfModel$mean, digits = 3)
OverallOfModel$Effect<-NA
for (g in 1:length(row.names(OverallOfModel)) ){
  print(g)
  if (g==1){
    OverallOfModel$Effect[1]<- round(exp(OverallOfModel[1,g]), digits = 3)
  }
  else{
    OverallOfModel$Effect[g]<- round((((exp(OverallOfModel[1,1] + OverallOfModel$Estimate[g]) - exp(OverallOfModel[1,1])) /exp(OverallOfModel[1,1]))*100)) #add estimate % column
  }
}
print(xtable(OverallOfModel, digits = 3),type="html",file="Output/OverallOfGLM.html")
#Use AIC function and fiddle with column exclusion for Overall Analysis
