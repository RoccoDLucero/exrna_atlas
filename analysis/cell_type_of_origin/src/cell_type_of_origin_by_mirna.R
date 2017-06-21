##cell_type_of_origin_by_mirna.R
##by: Rocco Lucero
##started: June 19, 2017
##last updated: June 20, 2017
##
##Objective: Identify tissue of origin, or cell type of origin signatures
##    in extracellular RNA component of various biofluids
##Input: ExceRpt processed EnTEX small RNA-seq data; exRNA Atlas Data
##Output:

source("./src/cell_type_of_origin_functions.R")
library(tsne)
library(gplots)
library(ggplot2)
library(huge)
source("https://bioconductor.org/biocLite.R")
biocLite("sva")
biocLite("preprocessCore")
get.qn <- preprocessCore::normalize.quantiles
devtools::install_github("BRL-BCM/EDec")
library(EDec)
set.seed(20170620)

##Function definitions
filterBySD <- function(x, min_sd){
    s <- sd(x)
    if(s < min_sd){return(TRUE)}else{return(FALSE)}

}
###################################################################################################
##Load the ENTEX data then tidy the global environment
load(file = "./input/getex/exceRpt_ENCODE_smallRNAseq_merged_results/exceRpt_smallRNAQuants_ReadsPerMillion.RData")

rm_dat <- grep(pattern = "exogenous", value = T, x = ls())

do.call(rm, as.list(rm_dat))

###################################################################################################
##create metadata dataframe from metadata embedded in ENTEX sample names (any table will suffice)
meta_all <- make.meta.from.entex.names(exprs.miRNA.rpm)


##For this analysis we only care to use data from adult donors (no cell line of fetal samples)
adult_samples <-  which((meta_all$group == "adult") & (meta_all$smp_src != 'GM12878'))

meta_adult <- meta_all[adult_samples,]

group <-  sapply(X = meta_adult$group, FUN = function(x){paste(x, collapse = " & ")})

sex <-  sapply(X = meta_adult$sex, FUN = function(x){paste(x, collapse = " & ")})

age <-  sapply(X = meta_adult$age, FUN = function(x){paste(x, collapse = " & ")})

tissue <- sapply(X = meta_adult$smp_src, FUN = function(x){paste(x, collapse = " & ")})

tissue_subset <- sample(x = tissue, size = 8)
samp_subset <- which(meta_adult$smp_src %in% tissue_subset)

###################################################################################################
##Preprocess the data: normalization and surrogate variable analysis
##to remove unmodeled sources of variance
mirna <- t(exprs.miRNA.rpm)


keep_mirna <- !(apply(X = mirna, MARGIN = 2, FUN = filterBySD, min_sd = 1))

mirna_adult <- mirna[adult_samples, keep_mirna]
rownames(mirna_adult) <- meta_adult$smp_src
mirna_adult <- mirna_adult[samp_subset, ]


mirna_adult_qn <- t(get.qn(t(mirna_adult)))

mirna_adult_npn <- huge.npn(mirna_adult)


#mirna_adult_s <- sweep(mirna_adult, 2, colSums(mirna_adult), FUN="/")
#mirna_adult_s <- scale(mirna_adult_s, center=FALSE, scale=colSums(mirna_adult_s))
#Run SVA and Combat...
#
# p <- prcomp(mirna_adult_npn)
# p <- p$x
# qplot(x = p[,1], y = p[,2], color = sex[samp_subset])
# qplot(x = p[,1], y = p[,2], color = age[samp_subset])
# qplot(x = p[,3], y = p[,4], color = tissue[samp_subset])

heatmap.2(x = mirna_adult_npn, trace = 'none')

# Create a vector of colors representing the class of each reference
tiss_colors <- as.factor(tissue[samp_subset])
levels(tiss_colors) <- RColorBrewer::brewer.pal(length(unique(tiss_colors)),"Accent")
tiss_colors <- as.character(tiss_colors)

# Create a color gradient to be used in a heatmap of correlations
color_gradient <- colorRampPalette(c("white","steelblue"))

# Compute correlation matrix and draw a heatmap
cors_tiss_mirna <- cor(x = t(mirna_adult))
heatmap.2(cors_tiss_mirna,
                  trace="none",
                  col=color_gradient(10),
                  breaks=seq(0,1,0.1),
                  margins=c(15,15),
                  ColSideColors = tiss_colors,
                  RowSideColors = tiss_colors)


