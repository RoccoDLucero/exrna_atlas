###############################################################################
#Load libraries and functions
#
source('./exrna_atlas_functions.R',echo = F)
################################################################################
####################################################################################################
#Load the Atlas Raw miRNA counts and associated sample metaData for ~750 samples of various
# conditions, rectify data types, and merge into a single data frame

Atlas.miRNA.nrmlzd <- read.delim2("./input/Atlas_Human_miRNA_Normalized_RPM_Counts.txt")
Atlas.miRNA.meta <- read.delim2("./input/Atlas_Human_Sample_Metadata2.txt",sep ="\t",row.names = 1)
Atlas.miRNA.raw <- read.delim2("./input/Atlas_Human_miRNA_Raw_Read_Counts.txt",sep = "\t",row.names = 1)
rn <- rownames(Atlas.miRNA.raw)
Atlas.miRNA.raw <- sapply(Atlas.miRNA.raw,as.numeric)
rownames(Atlas.miRNA.raw) <- rn
Atlas.miRNA.raw <- as.data.frame(Atlas.miRNA.raw)
Atlas.miRNA <- merge.data.frame(Atlas.miRNA.meta, t(Atlas.miRNA.raw), by.x = 0, by.y = 0 )
rm(Atlas.miRNA.meta,Atlas.miRNA.raw)
################################################################################
colnames(Atlas.miRNA[1:7])
levels(Atlas.miRNA$Anatomical.Location)
levels(Atlas.miRNA$Biofluid.Name)
sapply(Atlas.miRNA[1:10],class)

tmp <- prcomp(Atlas.miRNA[,8:ncol(Atlas.miRNA)],scale. = T)
#print(tmp)
tmp1 <- predict(tmp)
tmp1 <- as.data.frame(tmp1)
levels(Atlas.miRNA$Condition)
library(ggplot2)
qplot(PC1, PC2, data=tmp1, colour= Atlas.miRNA$Condition, shape = Atlas.miRNA$Biofluid.Name )

qplot(PC1, PC2, data=tmp1, colour= Atlas.miRNA$Biofluid.Name, shape = Atlas.miRNA$exRNA.Source )
qplot(PC1, PC2, data=tmp1, colour= Atlas.miRNA$Anatomical.Location)
qplot(PC1, PC2, data=tmp1, colour= Atlas.miRNA$exRNA.Source)
qplot(PC1, PC2, data=tmp1, colour= Atlas.miRNA$RNA.Isolation.Kit)

qplot(PC1, PC3, data=tmp1, colour= Atlas.miRNA$Condition)
qplot(PC1, PC3, data=tmp1, colour= Atlas.miRNA$Biofluid.Name)
qplot(PC1, PC3, data=tmp1, colour= Atlas.miRNA$Anatomical.Location)
qplot(PC1, PC3, data=tmp1, colour= Atlas.miRNA$exRNA.Source)
qplot(PC1, PC3, data=tmp1, colour= Atlas.miRNA$RNA.Isolation.Kit)



