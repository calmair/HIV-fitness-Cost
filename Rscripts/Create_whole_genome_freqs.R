library(ape)
library(seqinr)
library(pegas)
library(sfsmisc)
library(ggplot2)
library(scales)
library(plotrix)
library(RColorBrewer)

#NEED TO FURTHER ANALYZE ON WHY THERE IS AN ABUNDANCE OF 0.5 FREQS THIS MAY BE LEADING TO THE FINAL PRODUCT NOT BEING SIGNIFICANT
#0.5 was created when all row were all zero's

consensusfasta<-read.dna("./Data/Modified_Sequences_Aligned2_Consensus/modifiedconsensusBsequence", format = "fasta",as.character=TRUE)
consensusB<-consensusfasta[1:9719] #changed consensus B to include whole genome
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

freqPatTs_Zanini<-data.frame(row.names=ZaniniFiles)
#import modified_act_files
ZaniniFiles<-list.files("Data/Modified_Frequencies_Aligned2_Consensus/",recursive = TRUE,pattern="tsv")
print(ZaniniFiles)


freqPatTs_Zanini<- data.frame()
for (i in 1:length(ZaniniFiles)){
  print(i)
  print(ZaniniFiles[i])
  SeqData<-read.csv(paste("Data/Modified_Frequencies_Aligned2_Consensus/",ZaniniFiles[i],sep=""),sep=",", stringsAsFactors = TRUE, header = TRUE)
  #SOLVED
  
  colnames(SeqData)[which(names(SeqData) == "X..A")] <- "A"
  
  #What is transition mut?
  SeqData$consensusB<-consensusB
  
  #Why use transition then compare transition to Consensus? 
  #Why not compare Consensus to the Zanini Sequence?
  
  for (j in 1:length(SeqData$consensusB)) SeqData$transition[j]<-transition(SeqData$consensusB[j])
  
  #Why use transition then compare transition to Consensus? 
  #Why not compare Consensus to the Zanini Sequence?
  
  
  #determine Ts freq of every site. 
  SeqData$freq<-0
  #check out freq at 6658
  #for i == 13:16
  for (k in 1:length(SeqData$consensusB)){#for each site in the sequence
    
    MutNum<- SeqData [k,which(c("a","c","g","t")==SeqData$transition[k])]
    WTNum <- SeqData [k,which(c("a","c","g","t")==SeqData$consensusB[k])]
    
    if(WTNum =="-"){SeqData$freq[k]<-0}
    if(WTNum !="-"){
      SeqData$freq[k]<-as.numeric(MutNum)/(as.numeric(MutNum)+as.numeric(WTNum)) 
      #if (as.numeric(MutNum)>=as.numeric(WTNum))SeqData$freq[k]<-NA #filter majority manority out 
    }
    if(WTNum == 0){SeqData$freq[k]<-0}
    freqPatTs_Zanini[i,k]<-SeqData$freq[k]
  }
  
  }

write.csv(freqPatTs_Zanini,file="Output/freqPatTs_Zanini.csv")

pdf("Output/ZaniniHistograms.pdf")
par(mfrow=c(3,3)) #rows of 3 in each pdf page
for (i in 1:length(freqPatTs_Zanini)) hist(freqPatTs_Zanini[,i],main=paste("Site",i,"WT:" ,SeqData$MajNt[i]),xlim=c(0,.01),ylim = c(0,10),breaks=seq(0,1,by=0.0001),xlab = "Frequency",ylab="Number of Patients",col=2)
dev.off()
