## Current approach to generate Venn Diagrams
## Samples: Healthy Control, Meets ERCC QC
## RNA types: miRNA, piRNA, tRNA, gencode
## Biofluids: CSF, Plasma, Saliva, Serum
## Threshold: Mean and Median >= 0, 25, 100

library(ggplot2)
library(gplots)
library(VennDiagram)
library(miscTools)

setwd("~/Documents/BCM/Lab/exRNA/miRNA_Table")
source('~/Documents/BCM/Lab/exRNA/miRNA_Table/functionsForMean.R')

miRNA = read.table("~/Documents/BCM/Lab/exRNA/miRNA_Table/Merged_Meta_Atlas_miRNA_AllReadsPerMillion.txt", header=TRUE, stringsAsFactors = TRUE, check.names = FALSE, comment.char = "", quote = "\"", sep = "\t", row.names = 1)

#############################################
## Generate necessary counts for Venn Diagram
## MEAN, Healthy Control, Biofluids
#############################################

levels.venn = c(25)
for (x in levels.venn){
  print(x)
  RNA.healthy.CSF.mean = meanFunction(miRNA, "Healthy Control", "Cerebrospinal fluid", x)
  RNA.healthy.plasma.mean = meanFunction(miRNA, "Healthy Control", "Plasma", x)
  RNA.healthy.saliva.mean = meanFunction(miRNA, "Healthy Control", "Saliva", x)
  RNA.healthy.serum.mean = meanFunction(miRNA, "Healthy Control", "Serum", x)
  
  maxlength.venn = max(length(RNA.healthy.CSF.mean), length(RNA.healthy.plasma.mean), length(RNA.healthy.saliva.mean), length(RNA.healthy.serum.mean))
  length(RNA.healthy.CSF.mean) = maxlength.venn
  length(RNA.healthy.plasma.mean) = maxlength.venn
  length(RNA.healthy.saliva.mean) = maxlength.venn
  length(RNA.healthy.serum.mean) = maxlength.venn
  
  cbind.col.venn = as.data.frame(cbind(RNA.healthy.CSF.mean, RNA.healthy.plasma.mean, RNA.healthy.saliva.mean, RNA.healthy.serum.mean))
  csv.name = paste(x,"_Threshold_miRNA.csv", sep = '')
  write.csv(cbind.col.venn, file = csv.name, quote = FALSE, na = "", row.names = FALSE)
  
  RNA.list = read.csv(csv.name, stringsAsFactors = FALSE, header = TRUE)
  RNA.list.venn = lapply(RNA.list, removeEMPTYstrings)
  names(RNA.list.venn) = c("CSF", "Plasma", "Saliva", "Serum")
  venn.title = paste(x,"_Threshold_miRNA", sep = '')
  
  VENN.LIST = RNA.list.venn
  venn.plot = venn.diagram(VENN.LIST, NULL, cex = 2, fill=c("darkmagenta", "darkblue", "yellow", "orange"), cat.fontface = 4, category.names = names(VENN.LIST), main = venn.title)
  pdf(paste(venn.title,"_Mean_Healthy_Biofluids.pdf", sep = ''))
  grid.draw(venn.plot)
  dev.off()
  
  venn.string = venn(VENN.LIST, show.plot=FALSE)
  str(venn.string)
  inters <- attr(venn.string,"intersections")
}
