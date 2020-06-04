#Resistances are from Wensing et al. 2017
#BOLDED Only
#Mult-NRTI's excluded
#Creates dataframes for Resistances to Entry inhibitors, Integrase inhibitors, Nucleoside Reverse Transcriptase inhibitors,
# Non-Nucleoside Reverse Transcriptase inhibitors, and Protease Inhibitors
Entrymuts<- data.frame("pos"=numeric(7),"wt"=numeric(7),"mut"=numeric(7))
Intmuts<- data.frame("pos"=numeric(7),"wt"=numeric(7),"mut"=numeric(7))
PImuts<-data.frame("pos"=numeric(14),"wt"=numeric(14),"mut"=numeric(14))
NRTImuts<-data.frame("pos"=numeric(15),"wt"=numeric(15),"mut"=numeric(15))
NNRTImuts<-data.frame("pos"=numeric(14),"wt"=numeric(14),"mut"=numeric(14))

Intmuts$pos<-c(66,92,121,143,147,148,155)
Intmuts$mut<-c("IAK","QG","Y","RHC","G","HKR","H")
Entrymuts$pos<-c(36,37,38,39,40,42,43)
Entrymuts$mut<-c("DS","V","AME","R","H","T","D")
PImuts$pos<-c(30,32,46,47,48,50,58,74,76,82,83,84,88,90)
PImuts$mut<-c("N","I","IL","VA","V","VL","E","P","V","AFLTS","D","V","S","M")
NRTImuts$pos<-c(41,62,65,67,70,74,75,77,115,116,151,184,210,215,219)
NRTImuts$mut<-c("L","V", "REN","N","RE","V","I","L","F","Y","M","VI","W","YF","QE")
NNRTImuts$pos<-c(100,101,103,106,108,138,179,181,188,190,221,225,227,230)
NNRTImuts$mut<-c("I","P","NS","MA","I","AGKQR","L","CIV","CLH","SA","Y","H","C","IL")

AllMuts<-rbind(Entrymuts,Intmuts,PImuts,NRTImuts,NNRTImuts)
AllMuts<-AllMuts[order(AllMuts$pos),]
