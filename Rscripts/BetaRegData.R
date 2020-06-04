library(ape)
library(dplyr)

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
#Will give final dataframe needed for GLM.R
AlleleConsensus<- read.csv("./Output/freqPatTs_Zanini.csv",stringsAsFactors = FALSE)
ZaniniFiles<-list.files("Data/Modified_Frequencies_Aligned2_Consensus/",recursive = TRUE,pattern="tsv")
consensusfasta<-read.dna("./Data/Modified_Sequences_Aligned2_Consensus/modifiedconsensusBsequence", format = "fasta",as.character=TRUE)
consensusB<-consensusfasta[1:9719]
TotalSeqData<-data.frame("Reference" = 1:(9719*50), "Patient" = '',"Day" = "","Position" = "",a=0,c=0,g=0,t=0,bigAAchange=0,nonsyn=0,stop=0,Res=0,shape=0,paringprob=0,CpG=0,TransMutFreq=0, Section1=0, Section2 = 0, Section3 = 0, Section4=0, Section5 = 0, Section6 = 0, Section7=0,stringsAsFactors = FALSE)
TotalSeqData<-TotalSeqData[,-1]
OverviewZaniniDF<-read.csv("./Output/OverviewSelCoeffZanini.csv") #why use OverviewZaniniDF
GenomeShape<-read.csv("./Output/shape-parameters.csv") #Use ShapeParameters.txt
start1 = 1
end1= 50
#function to judge amino acid categories
pos <- "R|H|K"
neg <- "D|E"
unc <- "S|T|N|Q"
spe <- "C|U|G|P"
hyd <- "A|I|L|F|M|W|Y|V"

amCat <- function(AA){
  if(regexpr(pos, AA) > 0){ return(0) }
  if(regexpr(neg, AA) > 0){ return(1) }
  if(regexpr(unc, AA) > 0){ return(2) }
  if(regexpr(hyd, AA) > 0){ return(3) }
  if(regexpr("C", AA) > 0){ return(4) }
  if(regexpr("U", AA) > 0){ return(6) }
  if(regexpr("G", AA) > 0){ return(7) }
  if(regexpr("P", AA) > 0){ return(8) }
  return(5)
}
Gene1<-1:636 #5UTR
Gene2<- 790:2292 #Gag
Gene3<-2358:5040 #Pol
Gene4<-5041:5557 #Vif
Gene5<-5559:5795 #Vpr
Gene6<-6225:8795 #Env
Gene7<-8797:9168 #Nef


#Creates a row called bigChange that will be used in next for loop
bigChange <- rep(NA, nrow(OverviewZaniniDF))
for(i in 1:nrow(OverviewZaniniDF)){
  
  WT <- amCat(OverviewZaniniDF[i,'WTAA'])
  MUT <- amCat(OverviewZaniniDF[i,'MUTAA'])
  if(WT == MUT){ bigChange[i] <- 0
  }
  else{
    bigChange[i] <- 1
  }
}

for (i in 1:length(consensusB)){
  TotalSeqData$Patient[start1:end1]<-AlleleConsensus[1:50,2] #Patient
  TotalSeqData$Day[start1:end1]<-AlleleConsensus[1:50,3] #day
  TotalSeqData$Position[start1:end1]<-i #Position
  
  #Adjust this to show Actuall numbers of WT
  #TotalSeqData[start1:end1,which(colnames(TotalSeqData)==consensusB[i])]<-1 #WT
  
  
  TotalSeqData$bigAAchange[start1:end1]<-bigChange[i] #bigAAchange
  TotalSeqData[start1:end1,which(colnames(TotalSeqData)==OverviewZaniniDF$TypeOfSite[i])]<-1 #Type of Site
  TotalSeqData$shape[start1:end1]<- GenomeShape$Shape[i] #shape
  TotalSeqData$paringprob[start1:end1]<-GenomeShape$Prob[i] #probability
  TotalSeqData$CpG[start1:end1]<- OverviewZaniniDF$makesCpG[i] #cpg
  TotalSeqData[start1:end1, which(colnames(TotalSeqData) == consensusB[i])]<-1
  TotalSeqData$TransMutFreq[start1:end1]<-AlleleConsensus[1:50,i+3]
  for (g in 1:7){
    if (length(which (i == get(paste("Gene",g,sep = "")))) > 0){ #need to modify
      TotalSeqData[start1:end1,paste("Section",g,sep = "")]<-1#assign 1 if position is in Gene 1,2,3 etc.
      break
    }
  }
  
  
  start1=start1+50
  end1=end1+50
}




# for (z in 1:length(ZaniniFiles)){ #Access correct files
#   print(z)
#   print(ZaniniFiles[z])
#   SeqData<-read.csv(paste("Data/Modified_Frequencies_Aligned2_Consensus/",ZaniniFiles[z],sep=""),sep=",", stringsAsFactors = FALSE, header = TRUE)
#   colnames(SeqData)[which(names(SeqData) == "X..A")] <- "a"
#   colnames(SeqData)[which(names(SeqData) == "C")] <- "c"
#   colnames(SeqData)[which(names(SeqData) == "G")] <- "g"
#   colnames(SeqData)[which(names(SeqData) == "T")] <- "t"
#   
#   multipos<-which(TotalSeqData$Patient == str_extract(ZaniniFiles[z], "p[[:digit:]]+") & TotalSeqData$Day == str_extract(ZaniniFiles[z], "[[:digit:]]+_days")) #find correct position for each file in TotalSeq variable
#   start2 = 1
#   for (h in multipos){
#     if (SeqData$MajNt[start2] == "-"){
#       TotalSeqData[h,4:7]<- "NA"
#     }
#     else{
#       TotalSeqData[h,4:7]<- SeqData[start2,1:4] 
#     }
#     start2= start2+1
# }
# }




write.csv(TotalSeqData,"./Output/TransFreqGLM.csv",row.names = FALSE)