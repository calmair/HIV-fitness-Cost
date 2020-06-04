library(ape)
library(seqinr)
library(pegas)
library(sfsmisc)
library(ggplot2)
library(scales)
library(plotrix)
library(RColorBrewer)
#Script to analyse the frequency data and associate with features for Zanini data. 
#Read the csv files 

#Overview is assigning nonsyn, when it's actually synonymous
#Example Number 9162
source("Rscripts/ResistanceMutations.r")

consensusfasta<-read.dna("./Data/Modified_Sequences_Aligned2_Consensus/modifiedconsensusBsequence", format = "fasta",as.character=TRUE)
consensusB<-data.frame(Num = 1:length(consensusfasta),WTnt = consensusfasta[1:9719]) #changed consensus B to include whole genome

transition<-function(nuc){
  if (nuc=="a") {return("g")
  }
  if (nuc=="g") {return("a")
  }
  if (nuc=="c") {return("t")
  }
  if (nuc =="t") {return("c")
  }
}

#Insert typeofsitefunction
typeofsitefunction<-function(WTcodon, mutantcodon){
  
  WTAA<-seqinr::translate(WTcodon)
  MUTAA<-seqinr::translate(mutantcodon)
  if (WTAA == MUTAA) return ("syn")
  else if (MUTAA == "*") return ("stop")
  else return ("nonsyn")
}


start1<- 1
end1<-636

start2<-790
end2<-2292

start3<-2358
end3<-5040 #was 5096 

start4<-5041 
end4<-5557 #* was 5619

start5<-5559
end5<-5795

start6<-6225
end6<-8795

start7<-8797
end7<-9168

OverviewDFZanini1<- consensusB[start1:end1,] #5UTR
OverviewDFZanini2<- consensusB[start2:end2,] #Gag
OverviewDFZanini3<- consensusB[start3:end3,] #pol ############Problem that makes all nonsyn for this region and on
OverviewDFZanini4<- consensusB[start4:end4,] #vif
OverviewDFZanini5<- consensusB[start5:end5,] #vpr
OverviewDFZanini6<- consensusB[start6:end6,] #env
OverviewDFZanini7<- consensusB[start7:end7,] #nef

TypeOfSite<-c()

for (j in 1:7){
for (codon in get(paste("start",j,sep = "")):get(paste("end",j,sep = ""))){#for each Open reading frame for each gene
  inlistpos<-which(get(paste("OverviewDFZanini",j,sep = "")) == codon)
  if (inlistpos%%3 == 1){
  WTcodon <- as.character(get(paste("OverviewDFZanini",j,sep = ""))[(inlistpos):(inlistpos+2),2])
  mutantcodon <- c(transition(WTcodon[1]), WTcodon[2:3])
  TypeOfSite[codon]<-typeofsitefunction(WTcodon,mutantcodon)
}
if (inlistpos%%3 == 2){
  WTcodon <- as.character(get(paste("OverviewDFZanini",j,sep = ""))[(inlistpos-1):(inlistpos+1),2])
  mutantcodon <- c(WTcodon[1],transition(WTcodon[2]), WTcodon[3])
  TypeOfSite[codon]<-typeofsitefunction(WTcodon,mutantcodon)
}
if (inlistpos%%3 == 0){
  WTcodon <- as.character(get(paste("OverviewDFZanini",j,sep = ""))[(inlistpos-2):(inlistpos),2])
  mutantcodon <- c(WTcodon[1:2], transition(WTcodon[3]))
  TypeOfSite[codon]<-typeofsitefunction(WTcodon,mutantcodon)
}

  #Adds Res tag if which function greater than 0
 if (length(which(codon == AllMuts$pos))>0){
   TypeOfSite[codon]<-"Res"
 }
}

  if (j!=7){
    if (get(paste("end",j,sep = "")) != (get(paste("start",j+1,sep = ""))-1)){
  for(B in (get(paste("end",j,sep = ""))+1):(get(paste("start",j+1,sep = ""))-1)){
    TypeOfSite[B]<-"NA"
  }
    }
  }
  if (j==7){
    for(Z in ((get(paste("end",j,sep = "")))+1):9719){
      TypeOfSite[Z]<-"NA"
    }
  }
  }
print(TypeOfSite)



EstimatedS <- function(mu, meanfreq){
  if (meanfreq == 0) return (NA)
  else return (min(c(mu/meanfreq,1)))
}

#Read the stored frequencies rather than calculating frequencies again
read.table("Output/freqPatTs_Zanini.csv",sep=",",header=TRUE,row.names=1)->freqPatTsZanini

freqPatTsZanini<-freqPatTsZanini[,-c(1,2)]
#has 69 observations
#CHANGE ALL 0'S IN freqPatTsZanini into NA just so mean values are not skewed

#freqPatTsZanini[freqPatTsZanini == 0] <- NA



colMeansTsZanini<-apply(freqPatTsZanini, 2 , mean, na.rm=TRUE) #FINDS THE MEAN for all 69 rows



#what does the 2 mean?

## Create overview dataframe and plot site frequency spectra 
#Only synonymous, non-synomous and stop codons are considered
#- for each mutation, determine whether it is synonymous, non-synonymous or creates a stop
#- add information on resistance  positions

#PSP Nov 11 2015 I renamed x OverviewDFZanini and newdata OverviewDFZaniniOrderedByFreq
numsitesZanini<-length(colMeansTsZanini)
print(length(colMeansTsZanini))
OverviewDFZanini<-data.frame(num=1:numsitesZanini,colMeansTsZanini)


for (i in 1:length(TypeOfSite)){
  OverviewDFZanini$TypeOfSite[i]<-TypeOfSite[i]
}




OverviewDFZanini$WTnt<-as.character(consensusB[1:numsitesZanini,2])

#Mut rates and sel coefficients from Abrams paper
read.csv("Data/HIVMutRates/HIVMutRates.csv")->mutrates

OverviewDFZanini$TSmutrate<-0
OverviewDFZanini$TSmutrate[OverviewDFZanini$WTnt=="a"]<-mutrates$Probability[mutrates$Nucleotide.substitution=="AG"]
OverviewDFZanini$TSmutrate[OverviewDFZanini$WTnt=="c"]<-mutrates$Probability[mutrates$Nucleotide.substitution=="CU"]
OverviewDFZanini$TSmutrate[OverviewDFZanini$WTnt=="g"]<-mutrates$Probability[mutrates$Nucleotide.substitution=="GA"]
OverviewDFZanini$TSmutrate[OverviewDFZanini$WTnt=="t"]<-mutrates$Probability[mutrates$Nucleotide.substitution=="UC"]

OverviewDFZanini$TSmutZan<-0
OverviewDFZanini$TSmutZan[OverviewDFZanini$WTnt=="a"]<-mutrates$ZaniniProb[mutrates$Nucleotide.substitution=="AG"]
OverviewDFZanini$TSmutZan[OverviewDFZanini$WTnt=="c"]<-mutrates$ZaniniProb[mutrates$Nucleotide.substitution=="CU"]
OverviewDFZanini$TSmutZan[OverviewDFZanini$WTnt=="g"]<-mutrates$ZaniniProb[mutrates$Nucleotide.substitution=="GA"]
OverviewDFZanini$TSmutZan[OverviewDFZanini$WTnt=="t"]<-mutrates$ZaniniProb[mutrates$Nucleotide.substitution=="UC"]

for (i in 1:length(OverviewDFZanini$num)){
  OverviewDFZanini$EstSelCoeff[i] = EstimatedS(OverviewDFZanini$TSmutrate[i],OverviewDFZanini$colMeansTsZanini[i])
  OverviewDFZanini$EstSelCoeffZan[i] = EstimatedS(OverviewDFZanini$TSmutZan[i],OverviewDFZanini$colMeansTsZanini[i])}

#OverviewDFZanini$EstSelCoeff= OverviewDFZanini$TSmutrate/OverviewDFZanini$colMeansTsZanini
#OverviewDFZanini$EstSelCoeff[OverviewDFZanini$EstSelCoeff>1]<-1


#WT AAs 

OverviewDFZanini$WTAA<-""
OverviewDFZanini$MUTAA<-""

for (z in 1:7){
  #Need to fix this so it correctly analyzes codons in each section!

for (i in get(paste("start",z,sep = "")):get(paste("end",z,sep = ""))){ #must be in correct ORF
  inlistpos<-which(get(paste("OverviewDFZanini",z,sep = "")) == i)
  if (inlistpos%%3==1) OverviewDFZanini$WTAA[i] = seqinr::translate(OverviewDFZanini$WTnt[c(i,i+1,i+2)])
  if (inlistpos%%3==2) OverviewDFZanini$WTAA[i] = seqinr::translate(OverviewDFZanini$WTnt[c(i-1,i,i+1)])
  if (inlistpos%%3==0) OverviewDFZanini$WTAA[i] = seqinr::translate(OverviewDFZanini$WTnt[c(i-2,i-1,i)])
}

#MUT AAs
  for (i in get(paste("start",z,sep = "")):get(paste("end",z,sep = ""))){ #must be in correct ORF
  inlistpos<-which(get(paste("OverviewDFZanini",z,sep = "")) == i)
  if (inlistpos%%3==1) OverviewDFZanini$MUTAA[i] = seqinr::translate(c(transition(OverviewDFZanini$WTnt[i]),OverviewDFZanini$WTnt[c(i+1,i+2)]))
  if (inlistpos%%3==2) OverviewDFZanini$MUTAA[i] = seqinr::translate(c(OverviewDFZanini$WTnt[c(i-1)],transition(OverviewDFZanini$WTnt[i]),OverviewDFZanini$WTnt[c(i+1)]))
  if (inlistpos%%3==0) OverviewDFZanini$MUTAA[i] = seqinr::translate(c(OverviewDFZanini$WTnt[c(i-2,i-1)],transition(OverviewDFZanini$WTnt[i])))
}
}

#Add whether CpG sites
#Need to edit to add C&G
#Current ones determined are A&T
OverviewDFZanini$makesCpG <- 0
for(i in 1:nrow(OverviewDFZanini)){
  trip <- OverviewDFZanini$WTnt[c(i, i+1)]
  print(i)
  if(trip[1] == "c" & trip[2] == "a" ){
    OverviewDFZanini$makesCpG[i+1] <- 1
  }else if(trip[1] == "t" & trip[2] == "g"){
    OverviewDFZanini$makesCpG[i] <- 1
  }
  else{
    next()
  }
}



write.csv(OverviewDFZanini,"Output/OverviewSelCoeffZanini.csv",row.names = FALSE)
