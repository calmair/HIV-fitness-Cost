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

consensusfasta<-read.dna("./Data/HIV1_CON_2004_POL_DNA.fasta", format = "fasta",as.character=TRUE)
consensusB<-consensusfasta$CONSENSUS_B[551:873]

transition<-function(nuc){
  if (nuc=="a") {return("g")
  }
  if (nuc=="g") {return("a")
  }
  if (nuc=="c") {return("t")
  }
  if (nuc=="t") {return("c")
  }
}


TypeOfSite<-c()
for (codon in 1:(length(consensusB)/3)){#for each codon in the sequence
  positions <- c(codon*3-2,codon*3-1, codon*3)
  WTcodon <- consensusB[positions]
  mutant1codon <- c(transition(WTcodon[1]), WTcodon[2:3])
  mutant2codon <- c(WTcodon[1],transition(WTcodon[2]), WTcodon[3])
  mutant3codon <- c(WTcodon[1:2], transition(WTcodon[3]))
  TypeOfSite<-c(TypeOfSite,typeofsitefunction(WTcodon,mutant1codon))
  TypeOfSite<-c(TypeOfSite,typeofsitefunction(WTcodon,mutant2codon))
  TypeOfSite<-c(TypeOfSite,typeofsitefunction(WTcodon,mutant3codon))
}
print(TypeOfSite)

EstimatedS <- function(mu, meanfreq){
  if (meanfreq == 0) return (1)
  else return (min(c(mu/meanfreq,1)))
}

#Read the stored frequencies rather than calculating frequencies again
read.table("Output/freqPatTs_Zanini.csv",sep=",",header=TRUE,row.names=1)->freqPatTsZanini
#has 69 observations
colMeansTsZanini<-apply(freqPatTsZanini, 2 , mean, na.rm=TRUE) #FINDS THE MEAN for all 69 rows
## Create overview dataframe and plot site frequency spectra 
#Only synonymous, non-synomous and stop codons are considered
#- for each mutation, determine whether it is synonymous, non-synonymous or creates a stop
#- add information on resistance  positions

#PSP Nov 11 2015 I renamed x OverviewDFZanini and newdata OverviewDFZaniniOrderedByFreq
numsitesZanini<-length(colMeansTsZanini)
OverviewDFZanini<-data.frame(num=1:numsitesZanini,colMeansTsZanini)

OverviewDFZanini$TypeOfSite<-"syn"
#why does it give me NA?
OverviewDFZanini$WTnt<-consensusB[1:numsitesZanini]

#Mut rates and sel coefficients from Abrams paper
#ARE THESE MUTATION RATES UNIVERSAL?
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
#DONT KNOW HOW THIS WORKS
OverviewDFZanini$WTAA<-""
for (i in 1:numsitesZanini){
  if (i%%3==1) OverviewDFZanini$WTAA[i] = seqinr::translate(OverviewDFZanini$WTnt[c(i,i+1,i+2)])
  if (i%%3==2) OverviewDFZanini$WTAA[i] = seqinr::translate(OverviewDFZanini$WTnt[c(i-1,i,i+1)])
  if (i%%3==0) OverviewDFZanini$WTAA[i] = seqinr::translate(OverviewDFZanini$WTnt[c(i-2,i-1,i)])
}

OverviewDFZanini$MUTAA<-""
#MUT AAs
#DONT KNOW HOW THIS WORKS
for (i in 1:numsitesZanini){
  if (i%%3==1) OverviewDFZanini$MUTAA[i] = seqinr::translate(c(transition(OverviewDFZanini$WTnt[i]),OverviewDFZanini$WTnt[c(i+1,i+2)]))
  if (i%%3==2) OverviewDFZanini$MUTAA[i] = seqinr::translate(c(OverviewDFZanini$WTnt[c(i-1)],transition(OverviewDFZanini$WTnt[i]),OverviewDFZanini$WTnt[c(i+1)]))
  if (i%%3==0) OverviewDFZanini$MUTAA[i] = seqinr::translate(c(OverviewDFZanini$WTnt[c(i-2,i-1)],transition(OverviewDFZanini$WTnt[i])))
}

#Add whether CpG sites
OverviewDFZanini$makesCpG <- 0
for(i in 1:nrow(OverviewDFZanini)){
  trip <- OverviewDFZanini$WTnt[c(i-1, i, i + 1)]
  if(trip[1] == "c" & trip[2] == "a" ){
    OverviewDFZanini$makesCpG[i] <- 1
  }
  if(trip[2] == "t" & trip[3] == "g"){
    OverviewDFZanini$makesCpG[i] <- 1
  }
}

OverviewDFZanini$NonCpG <- 0
for(i in 1:nrow(OverviewDFZanini)){
  trip <- OverviewDFZanini$WTnt[c(i-1, i, i + 1)]
  if(trip[1] != "c" & trip[2] == "a" ){
    OverviewDFZanini$NonCpG[i] <- 1
  }
  if(trip[2] == "t" & trip[3] != "g"){
    OverviewDFZanini$NonCpG[i] <- 1
  }
}


write.csv(OverviewDFZanini,"Output/OverviewSelCoeffZanini.csv")
