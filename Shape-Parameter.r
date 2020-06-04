#Shape-Parameter.csv

GenomeShape<-read.delim("./Data/shape-parameters.txt", header = TRUE, sep = "\t", dec = ".")
GenomeShape<-GenomeShape[,-c(4,5,6,7,8,9)]
write.csv(GenomeShape,"Output/shape-parameters.csv",row.names = FALSE)