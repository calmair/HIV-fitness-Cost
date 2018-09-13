library(seqinr)
Seq<-vector()
ZaniniFiles<-list.files("Data/ZaniniNeherData/",recursive = TRUE,pattern="tsv")
for (i in 1:length(ZaniniFiles)){
  SeqData<-read.csv(paste("Data/ZaniniNeherData/",ZaniniFiles[i],sep=""),sep="\t",skip=1)
  SeqData$MajNt<-""
  SeqData$MajNt<-apply(SeqData[,1:4],1,function(x) c("a","c","g","t")[which.max(x)])
  Seq[i]<-toString(SeqData$MajNt, width = NULL)
  #write.fasta(Seq[i], names = paste("SEQUENCE", i) , file.out = paste("pat_", i), open = "w", nbchar = 60, as.string = FALSE)
}

#Notes: Need to change variable Seq into readable string data for matchPattern function
polstart<-vector()
for (i in 1:length(ZaniniFiles)){
  fomattedseq = gsub(", ","", Seq[i])    #Changes Seq[i] to remove comma
  total <- regexpr(pattern = "cctcagatcactcttt", fomattedseq)[[1]] #assigns the start location of the sequence in to total
  totalx <- regexpr(pattern = "cctcaaatcactcttt", fomattedseq)[[1]]
  if(total >= 0) { #these if statements make sure that each sequence returns a location
    polstart[i] <- total
  }else{
    polstart[i] <- totalx
  }
}
polstart

?regexpr
