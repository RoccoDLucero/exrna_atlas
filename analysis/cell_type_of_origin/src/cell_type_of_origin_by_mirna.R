##cell_type_of_origin_by_mirna.R
##by: Rocco Lucero
##started: June 19, 2017
##last updated: June 20, 2017
##
##Objective: Identify tissue of origin, or cell type of origin signatures
##    in extracellular RNA component of various biofluids
##Input: ExceRpt processed EnTEX small RNA-seq data; exRNA Atlas Data
##Output:

source("../cell_type_of_origin/src/cell_type_of_origin_functions.R")
library(tsne)
library(gplots)
library(ggplot2)
library(huge)
source("https://bioconductor.org/biocLite.R")
biocLite( c("sva", "preprocessCore"), suppressUpdates = T)
get.qn <- preprocessCore::normalize.quantiles
devtools::install_github("BRL-BCM/EDec")
library(EDec)
devtools::install_github("BRL-BCM/EDecExampleData")
library(EDecExampleData)
set.seed(20170620)

##Function definitions
filterBySD <- function(x, min_sd){
    s <- sd(x)
    if(s < min_sd){return(TRUE)}else{return(FALSE)}

}

txfmQN  <- function(df){preprocessCore::normalize.quantiles(df)} #columns will share a distribution
txfmNPN <- function(df){huge::huge.npn(x = df)}
txfmLinear <- function(x){(x-min(x))/(max(x)-min(x))}
txfmLogistic1 <- function(x, a = 1/max(x)){ 1 / (1 + exp(1)^(-(a * x )) ) }
txfmLogistic2 <- function(x, a = 1/max(x)){ 1.5*( 1 - exp(1)^(-(a * x)) ) }
txfmLogistic3 <- function(x, a = 1/max(x)){ 1 / (1 + exp(1)^(-(a * (x-mean(x)) ) )) }
txfmLogit <- function(x){log( x / (1 - x) )}

####################################################################################################
####################################################################################################
##Preprocess the data: normalization and surrogate variable analysis
##to remove unmodeled sources of variance

###################
#######ENTEX######

##Load ENTEX rpm data##
load(file = "./input/getex/exceRpt_ENCODE_smallRNAseq_merged_results/exceRpt_smallRNAQuants_ReadsPerMillion.RData")
entex_mirna <- t(exprs.miRNA.rpm)

##tidy the global environment##
rm_dat <- grep(pattern = "exogenous", value = T, x = ls())
do.call(rm, as.list(rm_dat))

##create metadata##
entex_meta_all <- make.meta.from.entex.names(exprs.miRNA.rpm)

##exclude cell line and fetal samples)##
adult_samples <-  which((entex_meta_all$group == "adult") & (entex_meta_all$smp_src != 'GM12878'))

entex_meta_adult <- entex_meta_all[adult_samples,]

group <-  sapply(X = entex_meta_adult$group, FUN = function(x){paste(x, collapse = " & ")})

sex <-  sapply(X = entex_meta_adult$sex, FUN = function(x){paste(x, collapse = " & ")})

age <-  sapply(X = entex_meta_adult$age, FUN = function(x){paste(x, collapse = " & ")})

tissue <- sapply(X = entex_meta_adult$smp_src, FUN = function(x){paste(x, collapse = " & ")})

tissue_subset <- sample(x = tissue, size = 7)

samp_subset <- which(entex_meta_adult$smp_src %in% tissue_subset)

table(unlist(entex_meta_adult[which(entex_meta_adult$smp_src %in% tissue_subset),][,'smp_src']))

##subset entex reads##
keep_mirna <- !(apply(X = entex_mirna, MARGIN = 2, FUN = filterBySD, min_sd = 1))
table(keep_mirna)

mirna_adult <- entex_mirna[adult_samples, keep_mirna]
rownames(mirna_adult) <- entex_meta_adult$smp_src
mirna_adult <- mirna_adult[samp_subset, ]

mirna_adult_qn <- t(get.qn(t(mirna_adult)))

mirna_adult_npn <- huge.npn(mirna_adult)

mirna_adult_npn_lgstc <- apply(X = mirna_adult_npn, MARGIN = 2, FUN = txfmLogistic1)
mirna_adult_npn_lgstc <- mirna_adult_npn_lgstc[,complete.cases(t(mirna_adult_npn_lgstc))]


###################
#######ATLAS######

##Load ExRNA Atlas rpm data##
exrna_atlas_1 <- readRDS("../get_atlas_data_in_R/interim/exrna_atlas_readcounts_non_gencode.RDS")
atlas_mirna <- exrna_atlas_1$miRNA
rm(exrna_atlas_1)

exrna_atlas_2 <- readRDS("./input/exrna_atlas_split_from_gencode_readcounts.Rdat")
atlas_linc <- exrna_atlas_2$lincRNA
rm(exrna_atlas_2)

##Load ExRNA Atlas metadata##
exrna_atlas_meta <- readRDS(file = '../cell_type_of_origin/input/exrna_atlas_metadata_and_map_to_samples.RDS')
exrna_atlas_meta

##subset atlas reads##
keep_mirna <- !(apply(X = atlas_mirna, MARGIN = 2, FUN = filterBySD, min_sd = 1))
table(keep_mirna)

atlas_mirna <- as.matrix(atlas_mirna[ , keep_mirna],
                         row.names = rownames(atlas_mirna),
                         col.names = colnames(atlas_mirna))

atlas_mirna_qn <- get.qn(atlas_mirna)
rownames(atlas_mirna_qn) <- rownames(atlas_mirna)
colnames(atlas_mirna_qn) <- colnames(atlas_mirna)

atlas_mirna_npn <- huge.npn(atlas_mirna_qn)

atlas_mirna_npn_lgstc <- apply(X = atlas_mirna_npn, MARGIN = 2, FUN = txfmLogistic1)
atlas_mirna_npn_lgstc <- atlas_mirna_npn_lgstc[,complete.cases(t(atlas_mirna_npn_lgstc))]



####################################################################################################
# Create a vector of colors representing the class of each reference
tiss_colors <- as.factor(tissue[samp_subset])
levels(tiss_colors) <- RColorBrewer::brewer.pal(length(unique(tiss_colors)),"Accent")
tiss_colors <- as.character(tiss_colors)

# Create a color gradient to be used in a heatmap of correlations
color_gradient <- colorRampPalette(c("white","steelblue"))
cbreaks <- 40
col_step <- 2/cbreaks

# Compute correlation matrix and draw a heatmap
cors_tiss_mirna <- cor(x = t(mirna_adult))
heatmap.2(cors_tiss_mirna,
                  trace="none",
                  #col=color_gradient(10),
                  #breaks=seq(0,1,0.1),
                  col= redgreen(cbreaks),
                  breaks=seq(-1,1,col_step),
                  margins=c(15,15),
                  ColSideColors = tiss_colors,
                  RowSideColors = tiss_colors)
#This demonstrates that the tissues have consistent and distinct miRNA signatures
# prior to any normalization and probe selection


cors_tiss_mirna <- cor(x = t(mirna_adult_npn))
heatmap.2(cors_tiss_mirna,
          trace="none",
          #col=color_gradient(10),
          #breaks=seq(0,1,0.1),
          col= redgreen(cbreaks),
          breaks=seq(-1,1,col_step),
          margins=c(15,15),
          ColSideColors = tiss_colors,
          RowSideColors = tiss_colors)
#This demonstrates that the normalized data shows more distinct tissue signatures
# prior to probe selection

cors_tiss_mirna <- cor(x = t(mirna_adult_npn_lgstc))
heatmap.2(cors_tiss_mirna,
          trace="none",
          #col=color_gradient(10),
          #breaks=seq(0,1,0.1),
          col= redgreen(cbreaks ),
          breaks=seq(-1,1,col_step),
          margins=c(15,15),
          ColSideColors = tiss_colors,
          RowSideColors = tiss_colors)



# Selecting marker loci based on comparisons of each class of reference against
# all other samples
mpv <- 1e-3
markers_ovr <- run_edec_stage_0(
                            reference_meth = t(mirna_adult_npn_lgstc),
                            reference_classes = tissue[samp_subset],
                            max_p_value = mpv,
                            num_markers = 100,
                            version = "one.vs.rest")

# Selecting marker loci based on comparisons of between each pair of
# reference classes
markers_ep <- run_edec_stage_0(
                            reference_meth = t(mirna_adult_npn_lgstc),
                            reference_classes = tissue[samp_subset],
                            max_p_value = mpv,
                            num_markers = 100,
                            version = "each.pair")



cors_tiss_markers_ovr <- cor(t(mirna_adult_npn_lgstc)[markers_ovr,])
heatmap.2(cors_tiss_markers_ovr,
                  trace="none",
                  #col=color_gradient(10),
                  #breaks=seq(0,1,0.1),
                  col= redgreen(cbreaks ),
                  breaks=seq(-1,1,col_step),
                  margins=c(10,10),
                  ColSideColors = tiss_colors,
                  RowSideColors = tiss_colors)

cors_tiss_markers_ep <- cor(t(mirna_adult_npn_lgstc)[markers_ep,])
heatmap.2(cors_tiss_markers_ep,
          trace="none",
          #col=color_gradient(10),
          #breaks=seq(0,1,0.1),
          col= redgreen(cbreaks ),
          breaks=seq(-1,1,col_step),
          margins=c(10,10),
          ColSideColors = tiss_colors,
          RowSideColors = tiss_colors)


cors_tiss_markers_all <- cor(t(mirna_adult_npn_lgstc)[c(markers_ep,markers_ovr),])
heatmap.2(cors_tiss_markers_all,
          trace="none",
          #col=color_gradient(10),
          #breaks=seq(0,1,0.1),
          col= redgreen(cbreaks ),
          breaks=seq(-1,1,col_step),
          margins=c(10,10),
          ColSideColors = tiss_colors,
          RowSideColors = tiss_colors)

####################################################################################################
## EDEC STAGE 1 ENTEX DATA
####################################################################################################
my_probes <- c(markers_ep, markers_ovr)
my_probes <-(intersect(my_probes, colnames(atlas_mirna_npn_lgstc)))

stage1_result_entex_5 = EDec::run_edec_stage_1(meth_bulk_samples = t(mirna_adult_npn_lgstc),
                           informative_loci = my_probes,
                           num_cell_types = 5)

dim(stage1_result_entex_5$methylation)


# Compute correlation between estimated methylation profiles,
# and reference methylation profiles

cors_deconv_refs_5ct = cor(t(mirna_adult_npn_lgstc)[my_probes,],
                           stage1_result_entex_5$methylation[my_probes,])

# Check what references had the highest correlation with each
# of the estimated methylation profiles
best_cors = rbind(apply(cors_deconv_refs_5ct,2,which.max),
                  apply(cors_deconv_refs_5ct,2,max))

best_cor_labels = matrix("",nrow=nrow(cors_deconv_refs_5ct),
                         ncol=ncol(cors_deconv_refs_5ct))
for (i in 1:5){
    best_cor_labels[best_cors[1,i],i] = as.character(round(best_cors[2,i],2))
}

# Plot correlation matrix
gplots::heatmap.2(cors_deconv_refs_5ct,
                  trace="none",
                  #col=color_gradient(10),
                  #breaks=seq(0,1,0.1),
                  col= redgreen(cbreaks ),
                  breaks=seq(-1,1,col_step),
                  margins=c(4,12),
                  RowSideColors = tiss_colors,
                  cellnote = best_cor_labels,
                  notecol="black")
#This shows that the selected probes can identify tissues of origin

colnames(stage1_result_entex_5$methylation) <- c("spleen", "gastro_sphincter", "stomach",
                                             "Peyer's", "esophageal mucosa")

colnames(stage1_result_entex_5$proportions) <- colnames(stage1_result_entex_5$methylation)


# Plot a heat map of proportions of constituent cell types in
# each sample with a side bar indicating whether the sample is
# tumor or normal
gplots::heatmap.2( stage1_result_entex_5$proportions,
                  trace="none",
                  #col=color_gradient(10),
                  #breaks=seq(0,1,0.1),
                  col= redgreen(cbreaks ),
                  breaks=seq(-1,1,col_step),
                  margins=c(10,4),
                  labRow = FALSE ,
                  RowSideColors = tiss_colors)


####################################################################################################
## EDEC STAGE 1 : ON EXRNA ATLAS DATA
####################################################################################################
my_probes <- c(markers_ep, markers_ovr)
my_probes <-(intersect(my_probes, colnames(atlas_mirna_npn_lgstc)))

stage1_result_atlas_5 = EDec::run_edec_stage_1(meth_bulk_samples = t(atlas_mirna_npn_lgstc),
                                           informative_loci = my_probes,
                                           num_cell_types = 5)

dim(stage1_result_atlas_5$methylation)


# Compute correlation between estimated methylation profiles,
# and reference methylation profiles
cors_deconv_refs_5ct = cor(t(mirna_adult_npn_lgstc)[my_probes,],
                           stage1_result_atlas_5$methylation[my_probes,])

# Check what references had the highest correlation with each
# of the estimated methylation profiles
best_cors = rbind(apply(cors_deconv_refs_5ct,2,which.max),
                  apply(cors_deconv_refs_5ct,2,max))

best_cor_labels = matrix("",nrow=nrow(cors_deconv_refs_5ct),
                         ncol=ncol(cors_deconv_refs_5ct))
for (i in 1:5){
    best_cor_labels[best_cors[1,i],i] = as.character(round(best_cors[2,i],2))
}

# Plot correlation matrix
gplots::heatmap.2(cors_deconv_refs_5ct,
                  trace="none",
                  #col=color_gradient(10),
                  #breaks=seq(0,1,0.1),
                  col= redgreen(cbreaks),
                  breaks=seq(-1,1,col_step),
                  margins=c(4,12),
                  #RowSideColors = tiss_colors,
                  cellnote = best_cor_labels,
                  notecol="white")


# Create a vector of colors indicating whether a sample represents a
# tumor or a normal control

colnames(stage1_result_entex_5$methylation) <- c("spleen", "gastro_sphincter", "stomach",
                                                 "Peyer's", "esophageal mucosa")[c(4,2,1,3,5)]

colnames(stage1_result_entex_5$proportions) <- colnames(stage1_result_entex_5$methylation)

#biofluid_colors =


# Plot a heat map of proportions of constituent cell types in
# each sample with a side bar indicating whether the sample is
# tumor or normal
gplots::heatmap.2(stage1_result_atlas_5$proportions,
                  trace="none",
                  col=color_gradient(10),
                  breaks=seq(0,1,0.1),
                  margins=c(10,4),
                  labRow = FALSE,
                  RowSideColors = biofluid_colors)

