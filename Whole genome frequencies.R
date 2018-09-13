#Working on this Sep 2015, preparing code for Marion and Kristof

#Working on this again in December 2014. Want to test how well we can estimate the mean frequency of a mutation. 
#I will determine mean freq for each site and then do bootstrapping to get 95% conf intervals
#This may be useful as prelim data for NSF proposal. 

#I plan to plot the frequencies of mutations WITHIN all patients. 
#So for each site, I determine the B consensus base, this will be WT
#Next, for each patient, I determine whether the seqs on day 1 were WT. 
#If that is the case then I will look at the freq of the non-WT bases on all days after day 1. 
#Accross all patients, I will have around one hundred frequencies for each base and each site. Now I look at the freq dist for each site. 
#If all 4-fold sites are truly neutral, then they should all show the same frequency distribution. 
#The distribution will not look like a neutral one because of the conditioning on starting off entirely WT. But it should still be OK to compare between sites, because I do the conditioning for each site.

#DEPENDS ON source("/Users/pleuni/Documents/Research/HIV/SoftSweepsInHIV/Bacheler2000/RResistanceMutations.r")
#DEPENDS ON "/Users/pleuni/Documents/Research/HIV/SoftSweepsInHIV/HowToGenbank/HIV1_CON_2004_POL_DNA.fasta"
#DEPENDS ON OR CREATES "freqPatSite.csv" (No longer, now this is made in createFrequencies-Bacheler.R)

#load relevant libraries and read consensusfasta file

library(ape)
library(seqinr)
library(pegas)
library(sfsmisc)
library(ggplot2)
library(scales)
library(plotrix)
library(RColorBrewer)

#read the file with the resistance mutations
source("Rscripts/RResistanceMutations.r")
#read the fasta file 
consensusfasta<-read.dna("./Data/HIV1_CON_2004_POL_DNA.fasta", format = "fasta",as.character=TRUE)	

#where is the start of POL? 
polstart=regexpr("cctca",paste(consensusfasta[which(row.names(consensusfasta)=="CONSENSUS_B"),],collapse=""))[1]
consensusB<-consensusfasta$CONSENSUS_B[551:9717] #changed consensus B to have whole genome
consensusA<-consensusfasta[which(row.names(consensusfasta)=="CONSENSUS_A1"), polstart:(polstart+983)]
consensusC<-consensusfasta[which(row.names(consensusfasta)=="CONSENSUS_C"), polstart:(polstart+983)]
consofcons<-consensusfasta[which(row.names(consensusfasta)=="CON_OF_CONS"), polstart:(polstart+983)]
consensus01AE<-consensusfasta[which(row.names(consensusfasta)=="CONSENSUS_01_AE"), polstart:(polstart+983)]
list.files(path="./Data/BachelerFiles/FASTAfiles/")->listfastafiles


#* Transition function*
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

typeofsitefunction<-function(WTcodon, mutantcodon){
  WTAA<-seqinr::translate(WTcodon)
  MUTAA<-seqinr::translate(mutantcodon)
  if (WTAA == MUTAA) return ("syn")
  else if (MUTAA == "*") return ("stop")
  else return ("nonsyn")
}

TypeOfSite<-c()
for (codon in 1:13) TypeOfSite<-c(TypeOfSite,c("overlap","overlap","overlap"))
for (codon in 14:(984/3)){#for each codon in the sequence
  positions <- c(codon*3-2,codon*3-1, codon*3)
  WTcodon <- consensusB[positions]
  mutant1codon <- c(transition(WTcodon[1]), WTcodon[2:3])
  mutant2codon <- c(WTcodon[1],transition(WTcodon[2]), WTcodon[3])
  mutant3codon <- c(WTcodon[1:2], transition(WTcodon[3]))
  TypeOfSite<-c(TypeOfSite,typeofsitefunction(WTcodon,mutant1codon))
  TypeOfSite<-c(TypeOfSite,typeofsitefunction(WTcodon,mutant2codon))
  TypeOfSite<-c(TypeOfSite,typeofsitefunction(WTcodon,mutant3codon))
}
#make sure that resistance sites in RT have a diff type of site
TypeOfSite[sort(c((RTImuts$pos*3)-2,(RTImuts$pos*3)-1,(RTImuts$pos*3)))+297]<-"res"
TypeOfSite[sort(c((Indinavirmuts$pos*3)-2,(Indinavirmuts$pos*3)-1,(Indinavirmuts$pos*3)))]<-"res"

EstimatedS <- function(mu, meanfreq){
  if (meanfreq == 0) return (1)
  else return (min(c(mu/meanfreq,1)))
}

#Amino acid changes
pos <- "R|H|K"
neg <- "D|E"
unc <- "S|T|N|Q"
spe <- "C|U|G|P"
hyd <- "A|I|L|F|M|W|Y|V"
amCat <- function(AA){
  if(regexpr(pos, AA) > 0){ return(0) }
  if(regexpr(neg, AA) > 0){ return(1) }
  if(regexpr(unc, AA) > 0){ return(2) }
  if(regexpr(spe, AA) > 0){ return(3) }
  if(regexpr(hyd, AA) > 0){ return(4) }
  return(5)
}

#Colors
cols <- c("#66CCEE","#228833","#CCBB44","#EE6677","#AA3377","#4477AA","#BBBBBB")
# from https://personal.sron.nl/~pault/

NonSynSites<-which(TypeOfSite=="nonsyn")
SynSites<-which(TypeOfSite=="syn")
#--------------------------------------------------------------------------------------


#Will read Zanini files, determine where POL is. 

#This is just for one patient. Later add others.
ZaniniFiles<-list.files("Data/ZaniniNeherData/",recursive = TRUE,pattern="tsv")


#SeqData<-read.csv(paste("..//Data/ZaniniNeherData/",ZaniniFiles[1],sep=""),sep="\t",skip=1)
#SeqData$MajNt<-""
#SeqData$MajNt<-apply(X[,1:4],1,function(x) c("a","c","g","t")[which.max(x)])

#i=0
#while(i , 10000){
#    print(i)
#    polstart=regexpr("cctca",paste(X$MajNt[i:9060],collapse="")) [[1]]
#    print(polstart)
#    print(seqinr::translate(X$MajNt[(i+polstart-1):(i+polstart+100)]))
#    i=i+polstart+1
#}

#start of POL is 1719 :-) in patient 1
#seqinr::translate(X$MajNt[(1719):(1821)])
#seqinr::translate(X$MajNt[(1719+297):(1821+297)])

#make dataframe with frequencies for all non-muts for all patients for all sites filtered with the WT threshold.
freqPatTs_Zanini<-data.frame(row.names=ZaniniFiles)

#formats data for next use
Seq<-vector()
fomattedseq<-vector()
for (i in 1:length(ZaniniFiles)){
  Seqpoldata<-read.csv(paste("Data/ZaniniNeherData/",ZaniniFiles[i],sep=""),sep="\t",skip=1)
  Seqpoldata$MajNt<-""
  Seqpoldata$MajNt<-apply(Seqpoldata[,1:4],1,function(x) c("a","c","g","t")[which.max(x)])
  Seq[i]<-toString(Seqpoldata$MajNt, width = NULL)
  fomattedseq[i] = gsub(", ","", Seq[i])
  #write.fasta(Seq[i], names = paste("SEQUENCE", i) , file.out = paste("pat_", i), open = "w", nbchar = 60, as.string = FALSE)
}

#first for loop that is has faulty inner for loop
for (i in 1:length(ZaniniFiles)){
  print(i)
  print(ZaniniFiles[i])
  SeqData<-read.csv(paste("Data/ZaniniNeherData/",ZaniniFiles[i],sep=""),sep="\t",skip=1)

  pat = substr(ZaniniFiles[i],1,regexpr("/",ZaniniFiles[i])[[1]]-1)
  print(pat)
  #whole genome analysis starts at location 4
  SeqData<-SeqData[4:nchar(fomattedseq[i]),] #this line sets which data points will be analyzed
  SeqData$MajNt<-""
  SeqData$MajNt<-apply(SeqData[,1:4],1,function(x) c("a","c","g","t")[which.max(x)])
  #check that the right position is read  in the right reading frame
  print(seqinr::translate(SeqData$MajNt[1:30]))
  print(seqinr::translate(SeqData$MajNt[298:(298+29)]))  
  #What is transition mut?
  SeqData$consensusB<-consensusB[1:length(SeqData[,1])]
  
  for (j in 1:length(SeqData$consensusB)) {
    SeqData$transition[j]<-transition(SeqData$consensusB[j])
  }
  
  
  #determine Ts freq of every site. 
  SeqData$freq<-0
  for (k in 1:length(SeqData$consensusB)){#for each site in the sequence
    MutNum<- SeqData [k,which(c("a","c","g","t")==SeqData$transition[k])]
    WTNum <- SeqData [k,which(c("a","c","g","t")==SeqData$consensusB[k])]
    #check wether the neighboring sequences are the same / WE CANT DO THIS FOR THE ZANINI DATA
    SeqData$freq[k]<-MutNum/(MutNum+WTNum) 
    if (MutNum>=WTNum)SeqData$freq[k]<-NA #filter majority manority out 
    freqPatTs_Zanini[i,k]<-SeqData$freq[k]
  }
}

write.csv(freqPatTs_Zanini,file="Output/freqPatTs_Zanini.csv")    

pdf("TRY.pdf")
par(mfrow=c(3,3))
for (i in 1:100) hist(freqPatTs_Zanini[,i],main=i,xlim=c(0,1),breaks=seq(0,1,by=0.01),col=2)
dev.off()