library(ape)
library(seqinr)
library(pegas)
library(sfsmisc)
library(ggplot2)
library(scales)
library(plotrix)
library(RColorBrewer)
library(DataCombine)

#creates modified sequences and frequency tables

#add more extensive comments
AlignedCount<-data.frame()
ZaniniFiles<-list.files("Data/ZaniniNeherData/",recursive = TRUE,pattern="tsv")
Aligned_Geneious<-read.fasta(file = "Data/AlignedGeneiousData.fasta", seqtype = "DNA", as.string = FALSE, forceDNAtolower = TRUE, set.attributes = TRUE, legacy.mode = TRUE, seqonly = FALSE, strip.desc = FALSE, bfa = FALSE, sizeof.longlong = .Machine$sizeof.longlong, endian = .Platform$endian, apply.mask = TRUE)
Realnames<-vector()
modiseq<-vector()
New <- matrix("-", ncol = 7, nrow = 1)
for(a in 1:length(Aligned_Geneious)){
Realnames[a]<-attr(Aligned_Geneious[[a]], "Annot")
}
Seqnames<-names(Aligned_Geneious)
Realnames_2<- gsub(">SEQUENCE_[0-9] ","",Realnames)
Realnames_2<- gsub(">SEQUENCE ","",Realnames_2)
Realnames_2<- gsub(">SEQUENCE_[0-9][0-9] ","",Realnames_2)

#Make a list to remove insertions in Consensus_B later in the program
addressdeletion<-c()
u<-1
for (h in 1:length(Aligned_Geneious$`HIV_Reference_(Consensus_B)`)){#10655
  if (Aligned_Geneious[[1]][h]=="-"){
    addressdeletion[u]<- h
    u<-u+1
  }
}
#length(Aligned_Geneious[[1]])


for (i in 1:length(ZaniniFiles)){ #rows 72
  UnalignedCount<- read.csv(paste("Data/ZaniniNeherData/",ZaniniFiles[i],sep=""),sep="\t",skip=1)#change the one to i in ZaniniFiles[]
  UnalignedCount$MajNt<-""
  ZaniniModifiedName<- gsub(".tsv","",ZaniniFiles[i])

  #Find correct row from Aligned_Geneious to match the correct Zaninifile due to unorganized Aligned_Geneious
  for (k in 1:length(Aligned_Geneious)){#70

    if (ZaniniModifiedName == Realnames_2[k]){
      print(Realnames_2[k])
      #Makes Consensus sequence if Match found
      
      #Kaho writing code atm to make sure it inserts N/A for two sequences that have 0,0,0,0
      
      UnalignedCount$MajNt<-apply(UnalignedCount[,1:4],1,function(x){ if (max(x)==0) "N"
        else c("a","c","g","t")[which.max(x)]})
      
          
      
      
      #print("Done1")
      
      #for (f in 1:length(Aligned_Geneious)){#70
      
        for(b in 1:length(Aligned_Geneious$`HIV_Reference_(Consensus_B)`)){#10655
          if (Aligned_Geneious[[k]][b]=="-"){
              UnalignedCount<-InsertRow(UnalignedCount, NewRow = New, RowNum = b)
          #Inserts a blank row into the UnalignedCount this ensures that it aligns
          #Now there should be at least 10655 rows in UnalignedCount
          
              
              }
        }
      
      #}
      #print("Done2")
     

      #Uses a vector of index values to remove rows
      #Needs work
      length(addressdeletion)
      #for(s in 1:length(addressdeletion)){
      AlignedCount<-UnalignedCount[-addressdeletion,]
      Aligned_Geneious[[k]]<-Aligned_Geneious[[k]][-addressdeletion]
      #length of AlignedCount should be 9719
      
      AlignedCount<-na.omit(AlignedCount)
      print(length(row.names(AlignedCount))) # checks if all rows are the same length
      
      # also make one HIV-B reference sequence for next geneious step
      
      #write fasta file of sequences to use for geneious
      #need to rename correctly for next rscript (freqs)
      write.fasta(Aligned_Geneious[k], names = gsub("/","_",ZaniniFiles[i])  ,file.out = paste("/Users/Calmair/HIV/Data/Modified_Sequences_Aligned2_Consensus/modifiedsequence",gsub("/","_",ZaniniFiles[i] ),sep = "_"), open = "w")

      write.csv(AlignedCount, file = paste("/Users/Calmair/HIV/Data/Modified_Frequencies_Aligned2_Consensus/modifiedfrequency",gsub("/","_",ZaniniFiles[i] ),sep = "_"),row.names = FALSE,quote = FALSE)
    
      }
  }
}
Aligned_Geneious[[1]]<-Aligned_Geneious[[1]][-addressdeletion]
write.fasta(Aligned_Geneious[1], names ="consensus_B"  ,file.out = paste("/Users/Calmair/HIV/Data/Modified_Sequences_Aligned2_Consensus/modifiedconsensusBsequence"), open = "w")


