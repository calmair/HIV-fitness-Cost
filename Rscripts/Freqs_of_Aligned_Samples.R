library(ape)
library(seqinr)
library(pegas)
library(sfsmisc)
library(ggplot2)
library(scales)
library(plotrix)
library(RColorBrewer)
library(stringr)

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

#import modified_act_files
ZaniniFiles<-list.files("Data/Modified_Frequencies_Aligned2_Consensus/",recursive = TRUE,pattern="tsv")
freqPatTs_Zanini<-data.frame(Patient = str_extract(ZaniniFiles, "p[[:digit:]]+"), Day= str_extract(ZaniniFiles, "[[:digit:]]+_days"))

for (i in 1:length(ZaniniFiles)){ #Access correct files
  print(i)
  print(ZaniniFiles[i])
  SeqData<-read.csv(paste("Data/Modified_Frequencies_Aligned2_Consensus/",ZaniniFiles[i],sep=""),sep=",", stringsAsFactors = TRUE, header = TRUE)

  
  colnames(SeqData)[which(names(SeqData) == "X..A")] <- "A"
  

  SeqData$consensusB<-consensusB
  
  for (j in 1:length(SeqData$consensusB)) SeqData$transition[j]<-transition(SeqData$consensusB[j])
  
  
  #determine Ts freq of every site. 
  SeqData$freq<-NA
  #ESTIMATES FREQUENCIES
  for (k in 1:length(SeqData$consensusB)){#for each site in the sequence

    #if(!is.na(as.numeric(as.character(SeqData[k,1]))+as.numeric(as.character(SeqData[k,2]))+as.numeric(as.character(SeqData[k,3]))+as.numeric(as.character(SeqData[k,4]))>1000)){ that's over 1000 reads
      MutNum<- SeqData [k,which(c("a","c","g","t")==SeqData$transition[k])]
      WTNum <- SeqData [k,which(c("a","c","g","t")==SeqData$consensusB[k])]
      #if(WTNum =="-"){SeqData$freq[k]<-NA}
      if(WTNum !="-"){
        SeqData$freq[k]<-as.numeric(as.character(MutNum))/(as.numeric(as.character(MutNum))+as.numeric(as.character(WTNum))) 
        #if (as.numeric(MutNum)>=as.numeric(WTNum))SeqData$freq[k]<-NA #filter majority manority out 
      }
      if(WTNum == 0){SeqData$freq[k]<-0} #removes 1's

    freqPatTs_Zanini[i,k+2]<-SeqData$freq[k]
  }
  #which(freqPatTs_Zanini %in% 1) #checks for 1's
  #which(freqPatTs_Zanini %in% "NA")
  
  }

write.csv(freqPatTs_Zanini,file="Output/freqPatTs_Zanini.csv")

pdf("Output/ZaniniHistograms.pdf")
par(mfrow=c(3,3)) #rows of 3 in each pdf page
for (i in 3:length(freqPatTs_Zanini)){
  if (length(which(is.na(freqPatTs_Zanini[,i])))==0){#exclude NA's
  hist(freqPatTs_Zanini[,i],main=paste("Site",i-2,"WT:" ,consensusB[i]),xlim=c(0,.01),ylim = c(0,10),breaks=seq(0,1,by=0.0001),xlab = "Frequency",ylab="Number of Patients",col=2)
  }
  }
  dev.off()
