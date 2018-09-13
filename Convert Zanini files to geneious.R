library(seqinr)
ZaniniFiles<-list.files("Data/ZaniniNeherData/",recursive = TRUE,pattern="tsv")
for (i in 1:length(ZaniniFiles)){
SeqData<-read.csv(paste("Data/ZaniniNeherData/",ZaniniFiles[i],sep=""),sep="\t",skip=1)
SeqData$MajNt<-""
SeqData$MajNt<-apply(SeqData[,1:4],1,function(x) c("a","c","g","t")[which.max(x)])
Seq[i]<-toString(SeqData$MajNt, width = NULL)
write.fasta(Seq[i], names = paste("SEQUENCE", i) , file.out = paste("pat_", i), open = "w", nbchar = 60, as.string = FALSE)

}
