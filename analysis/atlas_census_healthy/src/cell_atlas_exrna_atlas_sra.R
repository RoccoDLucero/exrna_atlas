################################################################################
#Load libraries and functions
#
source('./exrna_atlas_functions.R',echo = F)
################################################################################
#Load Atlas miRNA readcounts
exRna.Atlas.miR <- read.delim2('./input/William/Atlas_Human_miRNA_Raw_Read_Counts.txt',
                               sep = "\t", row.names = 1)
exRna.Atlas.miR <- to.numeric.df(tFrame(exRna.Atlas.miR))

Atlas.miR.meta <- read.delim2("./input/William/Atlas_Human_Sample_Metadata2.txt",
                              sep ="\t",row.names = 1)
##########
##########
Tissue.Atlas <- read.delim2("./input/Merged_results_joel/exceRpt_miRNA_ReadsPerMillion.txt",
                          sep = "\t",row.names = 1,stringsAsFactors = F)
Tissue.Atlas <- to.numeric.df(tFrame(Tissue.Atlas))

Tissue.group.defs <- read.delim2("./input/Merged_results_joel/exceRpt_sampleGroupDefinitions.txt",
                             sep = "\t",row.names = 1,stringsAsFactors = T)
tss.rlvt.smp <- grep("CSF|brain|Plamsa|Serum|prostate|pancreas|colon", Tissue.group.defs$sampleGroup,value = F)
Tissue.group.defs1 <-Tissue.group.defs[tss.rlvt.smp,]

##########
##########
BrainSpan_SRA <- read.delim2("./input/Brainspan+SRA_merged_results/exceRpt_miRNA_ReadsPerMillion.txt",
                          sep = "\t",row.names = 1,stringsAsFactors = F)
BrainSpan_SRA <- to.numeric.df(tFrame(BrainSpan_SRA))
SRA.meta <- read.delim2("./input/Brainspan+SRA_merged_results/run_annotations_full--fixed.txt",
                        sep ="\t",blank.lines.skip = T)
#Until we get the full meta data uses onlt the SRA samples (for which we have metadata)
SRA.miR <- BrainSpan_SRA[SRA.meta$ID,]

shared.mirna <- intersect( colnames(Tissue.Atlas),colnames(SRA.miR) )
shared.mirna <- intersect( colnames(exRna.Atlas.miR),shared.mirna )

SRA.miR <- SRA.miR[,shared.mirna]
Tissue.Atlas <- Tissue.Atlas[,shared.mirna]
exRna.Atlas.miR <- exRna.Atlas.miR[,shared.mirna]
sapply(X = list(SRA.miR,Tissue.Atlas,exRna.Atlas.miR),FUN = dim)

if (   all(colnames(SRA.miR) == colnames(exRna.Atlas.miR))
    && all(colnames(SRA.miR) == colnames(Tissue.Atlas)) ){
    miRNA.RPM <- rbind.data.frame(Tissue.Atlas,SRA.miR,exRna.Atlas.miR)
}

#my.excluded.subset <- (grep("Tissue_",rownames(Cell.Atlas)))
#my.cell.atlas <- Cell.Atlas[-my.excluded.subset,]
#group.defs <- group.defs[-my.excluded.subset,]
#Cell.Atlas.QN <- Cell.Atlas.QN[-my.excluded.subset,]
#Cell.Atlas <- Cell.Atlas[(grep("Tissue_",rownames(Cell.Atlas))),]
#rownames(Cell.Atlas.QN)
#colnames(Cell.Atlas)

################################################################################
##Look at overview of PCA results
#res.pca <- prcomp(my.cell.atlas[,-2141])
if(F){
res.pca <- prcomp(Cell.Atlas.QN)
eig <- (res.pca$sdev)^2
component.variance <- eig/sum(eig)*100
cumvar <- cumsum(component.variance)

pdf('./output/Cell_Atlas/CellAtlas_PCA_Scree.pdf')
PC4scree <- 12
barplot(component.variance[1:PC4scree], names.arg=1:length(component.variance[1:PC4scree]), 
        main = "Cell Atlas summary statistics:\nVariance of PCs 1-12",
        sub = paste("(",round(cumvar[PC4scree],1),"% variance explained )"),
        xlab = "Principal Components",
        ylab = "Percentage of variances",
        col ="skyblue")
dev.off()

#circ.var <- get_pca_var(res.pca = res.pca)
########################
##This plots the the feature vectors on the unit circle with their labels
#Do this when the metadata are rich and informative
#pdf('./output/Cell_Atlas/CellAtlas_FACMAP_1_2.pdf')
#fviz_pca_var(res.pca,labelsize = 1,repel = T,axes = c(1,2) )
#dev.off()
# Helper function : 
# Correlation between variables and principal components
var_cor_func <- function(var.loadings, comp.sdev){
    var.loadings*comp.sdev
}
# Variable correlation/coordinates
loadings <- res.pca$rotation
sdev <- res.pca$sdev
var.coord <- var.cor <- t(apply(loadings, 1, var_cor_func, sdev))
head(var.coord[, 1:4])

# Plot the correlation circle
if(F){
a <- seq(0, 2*pi, length = 100)
plot( cos(a), sin(a), type = 'l', col="gray",
      xlab = "PC2",  ylab = "PC3")
abline(h = 0, v = 0, lty = 2)
# Add active variables
arrows(0, 0, var.coord[, 2][1:10], var.coord[, 3][1:10], 
       length = 0.1, angle = 15, code = 2)
# Add labels
text(var.coord, labels=rownames(var.coord), cex = 1, adj=1)
}
###############################################

###########################


#print(tmp)
tmp1 <- predict(res.pca)
tmp1 <- as.data.frame(tmp1)
library(ggplot2)
pdf('./output/Cell_Atlas/CellAtlas_PCA.pdf',onefile = T)
qplot(PC1, PC2, data=tmp1, colour=group.defs) #, shape = Cell.Atlas$Group )
qplot(PC2, PC3, data=tmp1, colour=group.defs) #, shape = KJ1.meta$Condition )
qplot(PC3, PC4, data=tmp1, colour=group.defs) #, shape = KJ1.meta$Condition )
dev.off()
}
####################