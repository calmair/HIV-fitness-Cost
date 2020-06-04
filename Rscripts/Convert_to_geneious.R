library(seqinr)
filename<-vector()
Seq<-vector()
ZaniniFiles<-list.files("Data/ZaniniNeherData/",recursive = TRUE,pattern="tsv")

for (i in 1:length(ZaniniFiles)){
SeqData<-read.csv(paste("Data/ZaniniNeherData/",ZaniniFiles[i],sep=""),sep="\t",skip=1)
#there is a problem with reading.csv which doesnt read up to a certain extent of the data
SeqData$MajNt<-""
colnames(SeqData)[which(names(SeqData) == "X..A")] <- "A"
#max(SeqData[1,1:4])
#SeqData$MajNt[2]
for (j in 1:nrow(SeqData)){
  #print(j)
  #print(max(SeqData[1,1:4]))
if (max(SeqData[j,1:4])!=0){
  SeqData$MajNt[j]<- noquote(names(which.max(SeqData[j,1:4])))
} else{
  SeqData[j,7]<-"N"
}
}
print(SeqData$MajNt)
Seq<-toString(SeqData$MajNt, width = NULL)
Seq<-gsub(", ","", Seq)
print(ZaniniFiles[i])
ZaniniFiles[i]<-gsub(".tsv","",ZaniniFiles[i])
filename[i]<-gsub("act_","",ZaniniFiles[i])
filename[i]<-gsub("/","_",filename[i])
filename[i]
write.fasta(Seq, names = paste("SEQUENCE", ZaniniFiles[i]) , file.out = paste("/Users/Calmair/HIV/Data/Zanini_Consensus_Sequences_Fastas",filename[i],sep="/"), open = "w", nbchar = 60, as.string = FALSE)
}
