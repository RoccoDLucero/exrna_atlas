data.files <- list.files("~/miRNA_Table", pattern="*.txt", ignore.case=FALSE)
ncolsData <- 0
nrowsData <- 0
miRNAsData <- NULL
sampleNames <- NULL
for (i in 1:length(data.files)) {
  currFile <- read.table(data.files[i], header=TRUE, sep="\t")
  miRNAsData <- c(miRNAsData, as.character(currFile[,"X"]))
  sampleNames <- c(sampleNames, colnames(currFile)[2:length(colnames(currFile))])
  ncolsData <- ncolsData + ncol(currFile)
  nrowsData <- nrowsData + nrow(currFile)
  print(sprintf("Unique miRNAs in integrated dataset = %d", length(unique(miRNAsData))))
  print(sprintf("Unique sample names in integrated dataset = %d", length(unique(sampleNames))))
  print(sprintf("Total rows = %d", nrowsData))
}


integratedData <- matrix(0, nrow=length(unique(miRNAsData)), ncol=length(unique(sampleNames)))
rownames(integratedData) <- unique(miRNAsData)
colnames(integratedData) <- unique(sampleNames)
for (i in 1:length(data.files)) {
  currFile <- read.table(data.files[i], header=TRUE, sep="\t")
  rownames(currFile) <- currFile[,"X"]
  currFile <- currFile[,-1]
  integratedData[rownames(currFile),colnames(currFile)] <- as.matrix(currFile)
}

write.table(integratedData, file="Atlas_Human_miRNA_Read_Counts.txt", sep="\t", quote=FALSE, row.names=TRUE)
