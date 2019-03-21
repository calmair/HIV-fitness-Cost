#New version that only delete's "-" in column 1

library(ape)
library(seqinr)
library(pegas)
library(sfsmisc)
library(ggplot2)
library(scales)
library(plotrix)
library(RColorBrewer)

Geneiousfasta<-read.dna("./Data/HIV_geneious_alignment.fasta", format = "fasta",as.character=TRUE)
m<-list()
k<-1
for (i in 1:length(Geneiousfasta[1,])){
  
  if (Geneiousfasta[1,i] == "-"){
    #Identifies which rows in column 1 contain "-"
    m[k]<-i
    k<-k+1
  }
  #makes sure there is no error due to overexceeding parameters
  if (i==length(Geneiousfasta[1,]))
  {
    break
  }
}
m<- unlist(m)
#Uses a vector of index values to remove rows
Geneiousfasta<-Geneiousfasta[,-m]
#Check to see if any "-" left
match("-",Geneiousfasta[1,])
#Creates new fasta file to be used in 5UTR_Freqs.R
write.csv(Geneiousfasta, file = "edited_geneious_alignement.csv")