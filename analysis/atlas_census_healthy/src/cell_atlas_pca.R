################################################################################
#Load libraries and functions
#
source('./exrna_atlas_functions.R',echo = F)
################################################################################
#Load the CEll Atlas miRNA Read counts
Cell.Atlas <- read.delim2("./input/Merged_results_joel/exceRpt_miRNA_ReadCounts.txt",
                          sep = "\t",row.names = 1,stringsAsFactors = F)

Cell.Atlas <- read.delim2("./input/Merged_results_joel/exceRpt_miRNA_ReadsPerMillion.txt",
                          sep = "\t",row.names = 1,stringsAsFactors = F)

Cell.Atlas <- to.numeric.df(tFrame(Cell.Atlas))

group.defs <- read.delim2("./input/Merged_results_joel/exceRpt_sampleGroupDefinitions.txt",
                             sep = "\t",row.names = 1,stringsAsFactors = T)

CA.Qnorm <- normalize.quantiles( as.matrix(Cell.Atlas) )
Cell.Atlas.QN <- as.data.frame(CA.Qnorm)
dimnames(Cell.Atlas.QN) <- dimnames(Cell.Atlas)

if( all(rownames(Cell.Atlas) ==  rownames(group.defs)) ){
    Cell.Atlas$Group <- as.factor(group.defs$sampleGroup)
}

my.excluded.subset <- (grep("Tissue_",rownames(Cell.Atlas)))
my.cell.atlas <- Cell.Atlas[my.excluded.subset,]
group.defs <- group.defs[my.excluded.subset,]
#Cell.Atlas <- Cell.Atlas[(grep("Tissue_",rownames(Cell.Atlas))),]
rownames(Cell.Atlas)
colnames(Cell.Atlas)
################################################################################
##Look at overview of PCA results
res.pca <- prcomp(my.cell.atlas[,-2141])
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

####################
