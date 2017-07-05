##cell_type_of_origin_by_mirna.R
##by: Rocco Lucero
##started: June 19, 2017
##last updated: July 3, 2017
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


min_sd <- 1
####################################################################################################
####################################################################################################
####################################################################################################
###################
#######ENTEX######

##Load ENTEX rpm data##
entex_data_file <- paste("./input/getex/exceRpt_ENCODE_smallRNAseq_merged_results/",
                         "exceRpt_smallRNAQuants_ReadsPerMillion.RData", sep ='')

load(file = entex_data_file)

##Select mirnatable
entex_mirna <- t(exprs.miRNA.rpm)

##create metadata##
entex_meta_all <- make.meta.from.entex.names(exprs.miRNA.rpm)

##exclude cell line and fetal samples)##
adult_samples <- with(entex_meta_all, (group == 'adult' & smp_src != 'GM12878'))

entex_meta_adult <- entex_meta_all[adult_samples,]
entex_meta_adult  <- droplevels.data.frame(x = entex_meta_adult)

tissue_subset <- grep('spleen|liver|lung|stomach', entex_meta_adult$smp_src, value = T )
samp_subset <- which(entex_meta_adult$smp_src %in% tissue_subset)
table(unlist(entex_meta_adult[which(entex_meta_adult$smp_src %in% tissue_subset),][,'smp_src']))

##subset entex reads##
keep_mirna <- !(apply(X = entex_mirna, MARGIN = 2, FUN = filterBySD, min_sd = min_sd))
table(keep_mirna)

mirna_adult <- entex_mirna[adult_samples, keep_mirna]
rownames(mirna_adult) <- rownames(entex_meta_adult)
mirna_adult <- mirna_adult[samp_subset, ]
entex_meta_adult <- entex_meta_adult[samp_subset, ]

tissue_refs_nrmlzd <- get.normalized.readcounts.obj(mirna_adult)

####################################################################################################
###################
#######ATLAS######

##Load ExRNA Atlas rpm data and metadata##
exrna_atlas_1 <- readRDS("../get_atlas_data_in_R/interim/exrna_atlas_readcounts_non_gencode.RDS")
exrna_atlas_meta <- readRDS(file = '../cell_type_of_origin/input/comprehensive_atlas_metadata_draft_1.RDS')

#select miRNA table
atlas_mirna <- exrna_atlas_1$miRNA

##subset samples on metadata properties##
blood_fractions <- grep("[Pp]lasma|Serum",exrna_atlas_meta$Biofluid.Name )

blood_samples <- exrna_atlas_meta[blood_fractions,]$Sample.Name

atlas_mirna <- atlas_mirna[blood_samples,]

##subset atlas reads and normalize##
keep_mirna <- !(apply(X = atlas_mirna, MARGIN = 2, FUN = filterBySD, min_sd = min_sd))
table(keep_mirna)

atlas_mirna <- atlas_mirna[ , keep_mirna]

atlas_mirna_nrmlzd <- get.normalized.readcounts.obj(atlas_mirna)


##tidy the global environment##
rm_dat <- grep(pattern = "exogenous", value = T, x = ls())
do.call(rm, as.list(rm_dat))
rm(exrna_atlas_1)

####################################################################################################

####################################################################################################

####################################################################################################
tiss_refs <- tissue_refs_nrmlzd$qn.npn.lgstc
tiss_meta <- droplevels.data.frame(entex_meta_adult)
rownames(tiss_refs) <- tiss_meta$smp_src

atlas_blood_mirna <- atlas_mirna_nrmlzd$qn.npn.lgstc
atlas_blood_meta <- droplevels.data.frame(exrna_atlas_meta[blood_fractions,])

refs_meta #add later
target_meta #add later

edec_input_data <- subset.on.shared.probes(ref_df = refs, exp_df = dat)



################################################################################
# Create a vector of colors representing the tissue type of each reference
factor_vars_list <- list( refs_meta$ )

factor_colors <-  lapply(factor_vars_list, make.colors.vector)
# Create a color gradient to be used in a heatmap of correlations

cbreaks <- 40
col_step <- 2/cbreaks
heat_breaks = seq(-1,1,col_step)
heat_colors = list(red_green = redgreen(cbreaks),
                   col_grad = colorRampPalette(c("white","steelblue")))

heat_margins = c(15,15)


################################################################################





# Selecting marker loci based on comparisons of each class of reference against
# all other samples
informative_probes <- perform.probe.selection( )

cors_ref <- cor(x = ref)
cors_tiss_markers_ovr <- cor(t(tiss_refs)[markers_ovr,])
cors_tiss_markers_ep <- cor(t(tiss_refs)[markers_ep,])

my_probes <- intersect(colnames(atlas_blood_mirna), union(markers_ep, markers_ovr))
cors_tiss_markers_all <- cor(t(tiss_refs)[my_probes,])

sapply(X = probe_sets_list, FUN = get_heatmaps)
heatmap.2(cors_ref,
          trace="none",
          col = heat_colors,
          breaks = heat_breaks,
          margins = heat_margins,
          ColSideColors = col_colors,
          RowSideColors = row_colors)

####################################################################################################
## EDEC STAGE 1 ENTEX DATA
####################################################################################################
n_ct <- 4

stage1_res_tiss_refs = EDec::run_edec_stage_1(meth_bulk_samples = t(tiss_refs),
                           informative_loci = my_probes,
                           num_cell_types = n_ct)
check_cors_tiss <-T
if(check_cors_tiss){
    # Compute correlation between estimated methylation profiles,
    # and reference methylation profiles
    cors_deconv_refs = cor(t(tiss_refs)[my_probes,],
                               stage1_res_tiss_refs$methylation[my_probes,])

    # Check what references had the highest correlation with each
    # of the estimated methylation profiles
    best_cors = rbind(apply(cors_deconv_refs,2,which.max),
                      apply(cors_deconv_refs,2,max))

    best_cor_labels = matrix("",nrow=nrow(cors_deconv_refs),
                             ncol=ncol(cors_deconv_refs))
    for (i in seq(n_ct)){
        best_cor_labels[best_cors[1,i],i] = as.character(round(best_cors[2,i],2))
    }

    # Plot correlation matrix
    gplots::heatmap.2(cors_deconv_refs,
                      trace="none",
                      #col=color_gradient(10),
                      #breaks=seq(0,1,0.1),
                      col= redgreen(cbreaks ),
                      breaks=seq(-1,1,col_step),
                      margins=c(4,12),
                      RowSideColors = tiss_colors,
                      cellnote = best_cor_labels,
                      notecol="black")
}

colnames(stage1_res_tiss_refs$proportions) <- colnames(stage1_res_tiss_refs$methylation)


gplots::heatmap.2( stage1_res_tiss_refs$proportions,
                  trace="none",
                  col=color_gradient(10),
                  breaks=seq(0,1,0.1),
                  margins=c(10,4),
                  labRow = FALSE ,
                  RowSideColors = tiss_colors)
####################################################################################################




####################################################################################################
## EDEC STAGE 1 : ON EXRNA ATLAS DATA no reference profiles
####################################################################################################
atlas_rc_subset <- atlas_blood_mirna
spike_in <- tiss_refs[1:4, ]

atlas_blood_pls_refs_mirna <- rbind(atlas_blood_mirna, spike_in)

atlas_rc_subset <- atlas_blood_pls_refs_mirna


stage1_res_atlas_pls_refs = EDec::run_edec_stage_1(meth_bulk_samples = t(atlas_rc_subset),
                                                  informative_loci = my_probes,
                                                  num_cell_types = n_ct)

colnames(stage1_res_atlas_no_refs$methylation) <- LETTERS[seq(n_ct)]

atlas_deconv_cors <- T
prof_a <- t(tiss_refs)[my_probes,]
prof_b <- stage1_res_atlas_pls_refs$methylation[my_probes, ]

if(atlas_deconv_cors){
    cors_deconv = cor(prof_a, prof_b)

    best_cors = rbind(apply(cors_deconv, 2, which.max), apply(cors_deconv, 2, max))

    best_cor_labels = matrix("", nrow=nrow(cors_deconv), ncol=ncol(cors_deconv))

    for (i in seq(n_ct)){
        best_cor_labels[best_cors[1, i], i] = as.character(round(best_cors[2, i], 2))
    }

    # Plot correlation matrix
    gplots::heatmap.2(cors_deconv,
                      trace="none",
                      #col=color_gradient(10),
                      #breaks=seq(0,1,0.1),
                      col= redgreen(cbreaks),
                      breaks=seq(-1,1,col_step),
                      margins=c(4,12),
                      RowSideColors = tiss_colors,
                      cellnote = best_cor_labels,
                      notecol="black")


}

if(T){
    # Plot a heat map of proportions of constituent cell types in
    # each sample with a side bar indicating whether the sample is
    # tumor or normal


    gplots::heatmap.2(stage1_res_atlas_pls_refs$proportions,
                      trace="none",
                      col=color_gradient(10),
                      breaks=seq(0,1,0.1),
                      margins=c(10,4),
                      labRow = FALSE,
                      RowSideColors = c(donor_type_colors, tiss_colors[1:4]))
}
