#Find Theys and Victoria data
#Format them to correct data frames
#Calculate CpG ratio for my dataset(8 genes), Theys(Pol gene), Victoria(8 genes)

Victoria<-read.csv(file="Data/VictoriaHIV/HIV1_FLT_2017_pol_DNA.csv")
names(Victoria)[3]<-"MeanFreq"
names(Victoria)[4]<-"WTnt"
#I'm currently using the CSV file provided on Victoria's github
#This data was extracted from http://www.hiv.lanl.gov/

#Theys data
Theys1<-read.csv(file="Data/TheysHIV/OverviewSelCoeff_BachelerFilter.csv")
Theys2<-read.csv(file="Data/TheysHIV/OverviewSelCoeffLehman.csv")
names(Theys2)[3]<-"MeanFreq"
#Insert CpG site column
Theys2$makesCpG <- 0
for(i in 1:(nrow(Theys2)-1)){
  trip <- Theys2$WTnt[c(i, i+1)]
  if(trip[1] == "c" & trip[2] == "a"){
    Theys2$makesCpG[i+1] <- 1
  }else if(trip[1] == "t" & trip[2] == "g"){
    Theys2$makesCpG[i] <- 1
  }
  else{
    next()
  }
}
Theys3<-read.csv(file="Data/TheysHIV/OverviewSelCoeffZanini.csv")
names(Theys3)[3]<-"MeanFreq"
#Insert CpG site column
Theys3$makesCpG <- 0
for(i in 1:(nrow(Theys3)-1)){
  trip <- Theys3$WTnt[c(i, i+1)]
  if(trip[1] == "c" & trip[2] == "a"){
    Theys3$makesCpG[i+1] <- 1
  }else if(trip[1] == "t" & trip[2] == "g"){
    Theys3$makesCpG[i] <- 1
  }
  else{
    next()
  }
}


#My data
Zanini<-read.table("Output/OverviewSelCoeffZanini.csv",sep=",",header=TRUE,row.names=1)
names(Zanini)[1]<-"MeanFreq"
Zanini1<- Zanini[1:634,1:10] #5UTR
Zanini2<- Zanini[790:2292,1:10] #Gag
Zanini3<- Zanini[2358:5040,1:10] #pol
Zanini4<- Zanini[5041:5557,1:10] #vif
Zanini5<- Zanini[5559:5795,1:10] #vpr
Zanini6<- Zanini[6225:8795,1:10] #env
Zanini7<- Zanini[8797:9168,1:10] #nef

dataNames<-c("Victoria","Theys1","Theys2","Theys3","Zanini1","Zanini2","Zanini3","Zanini4","Zanini5","Zanini6","Zanini7")
OfficName<-c("Victoria Pol","TheysPol-1","TheysPol-2","TheysPol-3","5UTR","Gag","Pol","Vif","Vpr","Env","Nef")
Ratio<-data.frame(Num=1:length(dataNames),dataset= dataNames,AG=0,TC=0)
for ( i in 1:length(dataNames)){
  print(dataNames[i])
#Mean frequencies of A->G syn that creates CpG
AGC<- get(paste(dataNames[i],sep = ""))$MeanFreq[which(get(paste(dataNames[i],sep = ""))$WTnt == "a" & get(paste(dataNames[i],sep = ""))$TypeOfSite == "syn" & get(paste(dataNames[i],sep = ""))$makesCpG == 1)]
AGC[which(AGC == 0)]<-NA
#A->G syn that does not create CpG
AGNC<-get(paste(dataNames[i],sep = ""))$MeanFreq[which(get(paste(dataNames[i],sep = ""))$WTnt == "a" & get(paste(dataNames[i],sep = ""))$TypeOfSite == "syn" & get(paste(dataNames[i],sep = ""))$makesCpG == 0)]
AGNC[which(AGNC == 0)]<-NA
#A->G Ratio
Ratio$AG[i]<-mean(AGNC,na.rm = TRUE)/mean(AGC,na.rm = TRUE)
#T->C syn that creates CpG
TCC<-get(paste(dataNames[i],sep = ""))$MeanFreq[which(get(paste(dataNames[i],sep = ""))$WTnt == "t" & get(paste(dataNames[i],sep = ""))$TypeOfSite == "syn" & get(paste(dataNames[i],sep = ""))$makesCpG == 1)]
TCC[which(TCC == 0)]<-NA
#T->C syn that does not create CpG
TCNC<-get(paste(dataNames[i],sep = ""))$MeanFreq[which(get(paste(dataNames[i],sep = ""))$WTnt == "t" & get(paste(dataNames[i],sep = ""))$TypeOfSite == "syn" & get(paste(dataNames[i],sep = ""))$makesCpG == 0)]
TCNC[which(TCNC == 0)]<-NA
#T->C Ratio
Ratio$TC[i]<-mean(TCNC,na.rm = TRUE)/mean(TCC,na.rm = TRUE)
}
png("Output/CpGRatio.png", width = 6.75, height = 6.75, units = "in", res= 300)
plot(x = 1,
     type = "n", xaxt="n",
     xlim = c(0.9, 11.1), xlab = "" ,
     ylim = c(-5, 10), ylab = "Cost of CpG Mutations",
     pch = 16)
ticks = c(0, 1, 5, 10)
axis(side = 2, at = ticks, labels = TRUE)
axis(1,at=Ratio$Num,labels = OfficName,las=2, cex.axis=0.85)
points(Ratio$Num+0.3, Ratio$AG, pch=16 ,col="red")
points(Ratio$Num-0.3, Ratio$TC, pch=16 ,col="blue")
abline(h=1,lty=1,lwd=.5)

abline(v=seq(0.5,11,1),lty=2,col="grey")
dev.off()


