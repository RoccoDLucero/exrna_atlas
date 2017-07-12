##cell_type_of_origin_by_generic.R
##by: Rocco Lucero
##started: June 19, 2017
##last updated: July 6, 2017
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


min_sd <- 0.07
normalization_types <- list("qn","qn.npn","qn.npn.lgstc","original")
names(normalization_types) <- normalization_types
normalization <- normalization_types$qn.npn.lgstc

if(T){
####################################################################################################
####################################################################################################
####################################################################################################
###################
#######ENTEX######

##Load ENTEX rpm data##
entex_data_file <- paste("./input/getex/exceRpt_ENCODE_smallRNAseq_merged_results/",
                         "exceRpt_smallRNAQuants_ReadsPerMillion.RData", sep ='')

load(file = entex_data_file)

##create metadata##
entex_meta_all <- make.meta.from.entex.names(exprs.miRNA.rpm)

################################################################################
##Filter on Metadata Properties of Reference Profiles
################################################################################
adult_samples <- with(entex_meta_all, (group == 'adult' & smp_src != 'GM12878'))
entex_meta <- entex_meta_all[adult_samples,]

#entex_meta$age <- factor_to_numeric(entex_meta$age)
#entex_meta <- subset.data.frame(x = entex_meta, subset = entex_meta$age > 37)

tissues_for_analysis <- 'spleen|lung|stomach|adrenal'
tissue_subset <- grep(tissues_for_analysis, entex_meta$smp_src, value = T )

samp_subset <- which(entex_meta$smp_src %in% tissue_subset)
entex_meta <- droplevels.data.frame(entex_meta[samp_subset, ])



table(unlist(entex_meta[which(entex_meta$smp_src %in% tissue_subset),][,'smp_src']))


##subset entex reads##
entex_rna_types <- list(exprs.miRNA.rpm, exprs.piRNA.rpm, exprs.tRNA.rpm)
entex_rna <- t(Reduce(f = rbind, entex_rna_types))

entex_rna <- entex_rna[adult_samples, ]

entex_rna <- entex_rna[samp_subset, ]

rownames(entex_rna) <- rownames(entex_meta)

tissue_refs_nrmlzd <- get.normalized.readcounts.obj(df = entex_rna, na_rm = T)
entex_rna <- tissue_refs_nrmlzd[[normalization]]


quantile(apply(X = entex_rna, MARGIN = 2, FUN = sd))
keep_rna <- !(apply(X = entex_rna, MARGIN = 2, FUN = filterBySD, min_sd = min_sd))
table(keep_rna)

entex_rna <- entex_rna[, keep_rna]
####################################################################################################
###################
#######ATLAS######

##Load ExRNA Atlas rpm data and metadata##
#atlas_exrna_tables_genc  <- readRDS("../get_atlas_data_in_R/interim/exrna_atlas_readcounts_gencode_split.RDS")
#atlas_exrna_tables_genc  <- readRDS("../get_atlas_data_in_R/interim/exrna_atlas_readcounts_gencode_mixed.RDS")
atlas_exrna_tables   <- readRDS("../get_atlas_data_in_R/interim/exrna_atlas_readcounts_non_gencode.RDS")
atlas_exrna <- Reduce(f = cbind, atlas_exrna_tables)

atlas_meta <- readRDS(file = '../cell_type_of_origin/input/comprehensive_atlas_metadata_draft_1.RDS')

##############################
blood_samples <- grep("plasma|serum", atlas_meta$Biofluid.Name, ignore.case = T)
healthy_samples <- grep("control", atlas_meta$Donor.Type, ignore.case = T)
kjens_study <- grep("KJENS", atlas_meta$Study, ignore.case = T)
sex <- grep("male", atlas_meta$Sex, ignore.case = T)

subset_options <- list(blood_samples, healthy_samples, kjens_study, sex)[c(1,2)]

atlas_subset <- Reduce(f = intersect, x = subset_options) # <===============
############################
atlas_meta   <- atlas_meta[atlas_subset,]

atlas_exrna <- atlas_exrna[as.character(atlas_meta$Sample.Name),]
all(rownames(atlas_exrna) == as.character(atlas_meta$Sample.Name))


##subset atlas reads and normalize##
atlas_nrmlzd <- get.normalized.readcounts.obj(df = atlas_exrna, na_rm = T)
atlas_exrna <- atlas_nrmlzd[[normalization]]

quantile(apply(X = entex_rna, MARGIN = 2, FUN = sd))
keep_exrna <- !(apply(X = atlas_exrna, MARGIN = 2, FUN = filterBySD, min_sd = min_sd))
table(keep_exrna)

atlas_exrna <- atlas_exrna[, keep_exrna]


##tidy the global environment##
rm_dat <- grep(pattern = "exogenous", value = T, x = ls())
do.call(rm, as.list(rm_dat))


####################################################################################################
refs <- entex_rna
refs_meta <- entex_meta
rownames(refs) <- make.names(refs_meta$smp_src,unique = T)

target <- atlas_exrna
target_meta <- atlas_meta

###############################################################################
probe_sel_max_p <- 1e-2
edec_input_data <- lapply(subset.on.shared.probes(ref_df = refs, exp_df = target),t)
lapply(edec_input_data, dim)

informative_probes <- perform.probe.selection(ref_probes =  edec_input_data$ref,
                                              ref_classes = refs_meta$smp_src,
                                              n_markers = 40,
                                              max_p_val = probe_sel_max_p)




#selected_probes <- select.best.probes(informative_probes, probes_per_comparison = 15)
tmp1 <- lapply(c(1:4), function(x){head(sort(informative_probes$one.vs.rest$t_test$p.values[,x]),20)})
tmp2 <- unique(names(unlist(tmp1)))
informative_probes1 <- tmp2

tmp11 <- lapply(c(1:4), function(x){head(sort(informative_probes$each.pair$t_test$p.values[,x]),20)})
tmp22 <- unique(names(unlist(tmp11)))
informative_probes2 <- tmp22
informative_probes <- union(informative_probes1, informative_probes2)
expression_correlations_refs <- list(test = cor(edec_input_data$ref[informative_probes,]))

#expression_correlations_refs <- lapply(X = informative_probes,
#                                       FUN = function(probes){cor(edec_input_data$ref[probes,],use = 'complete.obs')})


###############################################################################
if(T){
factor_vars_list <- list(refs_classes = refs_meta$smp_src, biofluid = target_meta$Biofluid.Name)
factor_colors <-  lapply(factor_vars_list, make.colors.vector)
col_colors <- factor_colors$refs_classes
row_colors <- factor_colors$refs_classes

#Function to set heatmap parameters {cbreaks heatmin heatmax colorstyle heatmargins, col_colors, row_colors }
cbreaks <- 20
heat_min <- 0
heat_max <- 1
heat_range <- (heat_max - heat_min)
col_step <- heat_range / cbreaks
heat_breaks <- seq(heat_min, heat_max, col_step)
n_colors    <- length(heat_breaks) - 1
heat_colors <- list(red_green = redgreen(n_colors), col_grad = colorRampPalette(c("white","steelblue"))(n_colors))
heat_colors <- heat_colors[[2]]
heat_margins <- c(15,15)

}

sapply(X = names(expression_correlations_refs),
       FUN = function(probeset){
           cors <- expression_correlations_refs
           main = probeset
           print(probeset)
           heatmap.2(x = cors[[probeset]], trace="none", col = heat_colors,
                      breaks = heat_breaks, margins = heat_margins, ColSideColors = col_colors, RowSideColors = row_colors,
                      main = main, key.xlab = "Expression correlation" )
})

 }
####################################################################################################
n_ct <- 5
my_probes <- informative_probes #$max_set
####################################################################################################
## EDEC STAGE 1 ENTEX DATA
####################################################################################################
stage1_refs <- run_edec_stage_1(meth_bulk_samples = edec_input_data$ref,
                           informative_loci = my_probes,
                           num_cell_types = n_ct)

check_stage1_results(stage1_result = stage1_refs, probeset = my_probes,
                     sample_meta_colors = factor_colors$refs_classes,
                     reference_profiles = edec_input_data$ref,
                     ref_colors = factor_colors$refs_classes)

####################################################################################################
## EDEC STAGE 1 : ON EXRNA ATLAS DATA with reference profiles
####################################################################################################
refs_spike_in_idx <-  c(1:5,8) #seq_along(colnames(edec_input_data$ref ))
refs_spike_in <- edec_input_data$ref[,refs_spike_in_idx]
exp_w_refs <- cbind(edec_input_data$exp, refs_spike_in)


stage1_exp_w_refs = EDec::run_edec_stage_1(meth_bulk_samples = exp_w_refs,
                                                  informative_loci = my_probes,
                                                  num_cell_types = n_ct)

colnames(stage1_exp_w_refs$methylation) <- LETTERS[seq(n_ct)]

check_stage1_results(stage1_result = stage1_exp_w_refs, probeset = my_probes, show_props = T,
                     sample_meta_colors = c(factor_colors$biofluid,
                                            factor_colors$refs_classes[refs_spike_in_idx]),
                     reference_profiles = edec_input_data$ref,
                     ref_colors =  factor_colors$refs_classes)


####################################################################################################
## EDEC STAGE 1 : ON EXRNA ATLAS DATA without reference profiles
####################################################################################################
if(F){
exp_wo_refs <- edec_input_data$exp

stage1_exp_wo_refs = EDec::run_edec_stage_1(meth_bulk_samples = exp_wo_refs,
                                            informative_loci = my_probes,
                                            num_cell_types = n_ct)

colnames(stage1_exp_wo_refs$methylation) <- LETTERS[seq(n_ct)]

check_stage1_results(stage1_result = stage1_exp_wo_refs, probeset = my_probes,
                     sample_meta_colors = factor_colors$biofluid,
                     reference_profiles = edec_input_data$ref,
                     ref_colors = col_colors)
}
