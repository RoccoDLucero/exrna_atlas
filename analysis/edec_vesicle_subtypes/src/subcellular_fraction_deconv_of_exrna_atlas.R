###############################
##  subcellular_fraction_deconv_of_exrna_atlas.R
##  by:Rocco Lucero
##  started: 8.16.2017
##    
###############################

################################################################################
##  LOAD PACKAGES & SOURCE FILES
################################################################################
rm(list=ls())

pkgs <- list('gplots', 'class', 'RColorBrewer', 'stringr', 'ggplot2', 'reshape2',
             'rms', 'devtools', 'clue', 'EDec', 'tidyverse', 'plot3D', 'tsne', 'huge')

lapply(pkgs, library, character.only = T)

source("https://bioconductor.org/biocLite.R")
biocLite( c("sva", "preprocessCore"), suppressUpdates = T)
get.qn <- preprocessCore::normalize.quantiles

source("../edec_vesicle_subtypes/src/edec_vesicle_subtypes_functions.R")


################################################################################
##  GLOBAL PARAMETERS; I/O LOCATIONS;
################################################################################
set.seed(20170816)

## INPUT LOCATIONS ##
ref_path  <- "../edec_vesicle_subtypes/input/rd_cnt/lasser_atlas_gingeras/ref/"
targ_path <- "../edec_vesicle_subtypes/input/rd_cnt/lasser_atlas_gingeras/target/"
prbs_path <- "../edec_vesicle_subtypes/input/rd_cnt/lasser_atlas_gingeras/prbs/"

## OUTPUT LOCATIONS ##
figs_dir <- "../edec_vesicle_subtypes/interim/figs/"
rdat_dir <- '../edec_vesicle_subtypes/interim/rdat/'

################################################################################
##  LOAD DATA AND METADATA
################################################################################
## FIRST GET FILENAMES AND PATHS FOR EXCERPT RESULT .RDATA ###
## FROM EACH NON-ATLAS STUDY:
## --HD/LD DATA FROM LASSER ET AL (later...)
## --GINGERAS CELLS AND SUBCELLULAR FRACTIONS FROM ENCODE
## ------SRA_STUDY_ACCESSION: SRP003754 
## **SEE get_gingeras_subcell_refs_rna_seq_raw_dat.R
##
## **ATLAS EXCERPT DATA IS PROCESSED SEPARATELY AND ELSEWHERE
## **SEE analysis "get_atlas_data_in_R" 
################################################################################



######################################
## LOAD REFERENCE DATA AND METADATA ##
gingeras_rpm <- condense.post.proc.results(results_rdata_path = ref_path,
                                         per_million = T,
                                         include_exogenous = F)

in_file <- "gingeras_subcell_meta.RDS"
gingeras_meta <- readRDS(paste(ref_path, in_file,sep = ''))

#######################################
##LOAD TARGET DATA AND METADATA FOR ALL STUDIES: ##
in_file <- "atlas_rpm_all_nc.RDS"
#atlas_rpm_all_nc <- read.delim(file = paste(targ_path, in_file, sep = ''), sep = "\t")
atlas_rpm_all_nc <- readRDS(paste(targ_path, in_file))

in_file <- "comprehensive_atlas_metadata_draft_1.RDS"
atlas_meta <-readRDS(file = paste(targ_path, in_file, sep = ''))

atlas_all_data_list_obj_path <- "./input/exrna_atlas/exrna_atlas_readcounts_non_gencode.RDS"
non_gc_names <- get.rna.species.by.type(atlas_all_data_list_obj_path)

atlas_all_data_list_obj_path <- "./input/exrna_atlas/exrna_atlas_readcounts_gencode_split.RDS"
gc_split_names <- get.rna.species.by.type(atlas_all_data_list_obj_path)

targ_probes_by_type <- c(non_gc_names,
                         list( snRNA = gc_split_names$snRNA,
                         lincRNA = gc_split_names$lincRNA,
                         snoRNA = gc_split_names$snoRNA))
##################################
## LOAD PRE_SELECTED PROBE SETS ##
## --OSCAR'S HD/LD PROBES
in_file <- "originalProbes_ncRNA.txt"
murillo_probes_hd_ld <- read.delim(paste(prbs_path, in_file, sep = ''), sep = "\t")

################################################################################
## COMPUTE AND OR ADD NEW METADATA PROPERTIES ##

## EXTRA METADATA PROPERTIES FOR ENCODE REFERENCES ##
#cell_stats <- list(cancer = c(hepg2, hela, k562, a549, sk-n-sh),
#                   normal= c(nhek, imr90, gm12878))

################################################################################
## PREPROCESS READCOUNT DATA OR PROBESETS ##
##--FIX SAMPLE AND RNA NAMES
##--RNA SUBSET IS DETERMINED BY RNA OSCAR SELECTED FOR 'atlas_all_ncrna_rpm.txt'

###############
## FIX NAMES ##
gingeras_dimnames <- process.dimnames(gingeras_rpm)
dimnames(gingeras_rpm) <- gingeras_dimnames

atlas_dimnames <- process.dimnames(atlas_rpm_all_nc)
dimnames(atlas_rpm_all_nc) <- atlas_dimnames

rownames(murillo_probes_hd_ld) <- murillo_probes_hd_ld$Val
murillo_probes_hd_ld <- process.dimnames(murillo_probes_hd_ld)[[1]]

targ_probes_by_type <- process.names.list(targ_probes_by_type)
################################################################################
## SUBSET READCOUNT DATA OR PROBESETS
################################################################################

###############################################
## FIRST CREATE VECTORS USED TO SUBSET DATA AND METADATA:
## --MATCH DATA TO METADATA BASED ON RUN ACCESSION (SRR ID)
## **--FOR NOW... NO ATLAS DATA SUBSETTING

## TO MATCH METADATA TO SAMPLE DATA
included_ref_samps <- gingeras_meta$run_accession %in% colnames(gingeras_rpm)
sum(included_ref_samps)

######################################################
## TO SELECT ONLY PROBES DETECTED IN APPROPRIATE SETS:

## PROBES IDENTIFIED IN BOTH THE REFERENCE AND TARGET SET ##
univ_detected_rpm <- subset.on.shared.probes(gingeras_rpm, atlas_rpm_all_nc)
univ_probes <- rownames(univ_detected_rpm$ref)

## PROBES IDENTIFIED IN EITHER THE REFERENCE AND/OR TARGET SET ##
max_probes <- union(gingeras_dimnames[[1]], atlas_dimnames[[1]])

## PROBES FOUND IN THE REFERENCE SET AND OSCARS PRE-SELECTED SET
murillo_probes <- intersect(murillo_probes_hd_ld, univ_probes)

random_probes_500 <- sample(x = univ_probes, size = 500)



probesets <- list(univ = univ_probes,
                  max = max_probes,
                  murillo = murillo_probes,
                  random_500 = random_probes_500)



#################################################
## NOW CREATE THE MAIN DATA AND METADATA OBJECTS:

## REFERENCE METADATA SUBSET
gingeras_meta <- gingeras_meta[included_ref_samps,]

## SELECT SMALL RNA ONLY ##
#gingeras_rpm <- subset.small.rna(gingeras_rpm,rna_in_colnames = F)
#atlas_rpm_all_nc <- subset.small.rna(atlas_rpm_all_nc, rna_in_colnames = F)

################################################################################
################################################################################
################################################################################
## CREATE REFERENCE SAMPLE LABELS FOR EDEC
## MAKE THIS EXPANDABLE TO MULTPLE METADATA VARIABLES IE A FUNCTION ##
if(T){
if(T){
meta_field <-gingeras_meta$sample_name
lbls_frac <- c("cytosol","nucleus")#,"cell")
use_lbls <- lbls_frac

ref_lbl <- sapply(X = use_lbls,
              FUN =  function(frac){
                  frac <- gsub(' ','',frac)
                  pos <- grep(frac, meta_field)
                  
                  clss <- rep(frac, length(pos))
                  
                  
                  nm <- gingeras_meta$run_accession[pos]
                  
                  
                  lbl_map <- cbind(nm,clss)
                  
                  return(lbl_map)},
              simplify = T)

ref_lbl <- Reduce(f = rbind, ref_lbl)
}

if(T){
meta_field <- gingeras_meta$sample_attribute
lbls_cell_ty <- unique(gsub("^source_name:|\\|\\|.*$","", meta_field))
use_lbls <- lbls_cell_ty

alt_lbl <- sapply(X = use_lbls,
                  FUN =  function(frac){
                      frac <- gsub(' ','',frac)
                      pos <- grep(frac, meta_field)
                      
                      clss <- rep(frac, length(pos))
                      
                      
                      nm <- gingeras_meta$run_accession[pos]
                      
                      
                      lbl_map <- cbind(nm,clss)
                      
                      return(lbl_map)},
                  simplify = T)

alt_lbl <- Reduce(f = rbind, alt_lbl)
}

if(F){
meta_field <- gingeras_meta$instrument_model
lbls_instr_ty <- unique(meta_field)
use_lbls <- lbls_instr_ty

alt_lbl <- sapply(X = use_lbls,
                  FUN =  function(frac){
                      
                      pos <- grep(paste("^",frac,"$",sep = ''), meta_field)
                      
                      clss <- rep(frac, length(pos))
                      
                      
                      nm <- gingeras_meta$run_accession[pos]
                      
                      
                      lbl_map <- cbind(nm,clss)
                      
                      return(lbl_map)},
                  simplify = T)

alt_lbl <- Reduce(f = rbind, alt_lbl)
}
}
#####################################################
## ENSURE THAT OTHER METADATA COME FROM
## SAMPLES SELECTED ON THE PRIMARY METATDATA PROPERTY
meta_labels <- merge.data.frame(ref_lbl, alt_lbl, by = 'nm' )
meta_labels <- sapply(X = meta_labels, as.character)

################################################################################
################################################################################
################################################################################

################################################################################
## END OF PRELIMINARIES
################################################################################

################################################################################
## BEGIN ANALYSIS
################################################################################

ref_samples <- meta_labels[,'nm']
ref_classes <- meta_labels[, 2]
refs <- gingeras_rpm[, ref_samples]

## ENSURE THAT REFERENCE METADATA MATCHES REFERENCE DATA 
meta_reorder <- match(x = ref_samples, table = gingeras_meta$run_accession)
refs_meta <- gingeras_meta[meta_reorder, ]

target <- atlas_rpm_all_nc
target_meta <- atlas_meta

target_samples <- colnames(target)

####################################
## MERGE AND NORMALIZE THE DATA
###################################

if(F){
## NORMALIZE DATA WITHOUT MERGING REF AND TARG
## SELECT NORMALIZATION METHOD
normalization_types <- list("qn","qn.npn","qn.npn.lgstc","original")
names(normalization_types) <- normalization_types
normalization <- normalization_types$qn.npn.lgstc

refs_nrmlzd_objs <- get.normalized.readcounts.obj(df = refs, na_rm = F)
normalized_refs <- refs_nrmlzd_objs[[normalization]]

targ_nrmlzd_objs <- get.normalized.readcounts.obj(df = target, na_rm = F)
normalized_targ <- targ_nrmlzd_objs[[normalization]]

}

if(F){
## NORMALIZE DATA FOR EACH BATCH SEPARATELY,


## SELECT NORMALIZATION METHOD
normalization_types <- list("qn","qn.npn","qn.npn.lgstc","original")
names(normalization_types) <- normalization_types
normalization <- normalization_types$qn.npn.lgstc
normalization_batch <- normalization_types$qn


## NORMALIZE EACH 'BATCH/STUDY' IN THE REFERENCE SET
batch1_samps <- grep("SRR527", refs_meta$run_accession)
batch2_samps <- grep("SRR527", refs_meta$run_accession, invert = T)
refs_batch1 <- refs[, batch1_samps]
refs_batch2 <- refs[, batch2_samps]

refs_nrmlzd_objs1 <- get.normalized.readcounts.obj(df = refs_batch1, na_rm = F)
normalized_refs1 <- refs_nrmlzd_objs1[[normalization_batch]]

refs_nrmlzd_objs2 <- get.normalized.readcounts.obj(df = refs_batch2, na_rm = F)
normalized_refs2 <- refs_nrmlzd_objs2[[normalization_batch]]

refs_batch_merged <- combine.rpm.tables(normalized_refs1, normalized_refs2)[, ref_samples]

refs_nrmlzd_objs <- get.normalized.readcounts.obj(df = refs_batch_merged, na_rm = F)
normalized_refs <- refs_nrmlzd_objs[[normalization]]

## NORMALIZE EACH 'BATCH/STUDY' IN THE EXRNA ATLAS
#targ_nrmlzd_objs <- get.normalized.readcounts.obj(df = target, na_rm = F)
#normalized_targ <- targ_nrmlzd_objs[[normalization]]

    if(F){
    #merge refs and target then normalize agian
    }
}

if(T){
## COMBINE DATA TABLES BEFORE NORMALIZATION
save_file <- 'all_nrmlzd_obj_cmbn_before_nrmlz.RDS'
nmz_obj <- paste(rdat_dir,save_file,sep = '')

## SELECT NORMALIZATION METHOD
normalization_types <- list("qn","qn.npn", "qn.lgstc","qn.npn.lgstc","original")
names(normalization_types) <- normalization_types
#normalization <- normalization_types$qn.npn.lgstc
normalization <- normalization_types$qn.lgstc

if(file.exists(nmz_obj)){
    
    all_nrmlzd_objs <- readRDS(nmz_obj)

}else{    
    merged_ref_target <- combine.rpm.tables(refs, target)

    all_nrmlzd_objs <- get.normalized.readcounts.obj(df = merged_ref_target, na_rm = T)
    
    saveRDS(all_nrmlzd_objs,file = paste(rdat_dir,save_file,sep = ''))
    
}

all_nrmlzd <- all_nrmlzd_objs[[normalization]]

#some NA values are present: "Neptune_Urine_33.11" "Neptune_Urine_33.24"
#setdiff(colnames(all_nrmlzd_objs$original), colnames(all_nrmlzd_objs$qn.npn.lgstc))
target_samples <- intersect(colnames(all_nrmlzd), target_samples)

normalized_refs <- all_nrmlzd[, ref_samples]
normalized_targ <- all_nrmlzd[, target_samples]

}

###############################################################################
##  PROBE SELECTION
###############################################################################

## SET PROBE SELECTION PARAMETERS ##
input_probes <- probesets$univ

probe_sel_max_p <- 1e-4
probe_sel_refs <- normalized_refs[input_probes, ]
n_classes <- length(unique(ref_classes))
n_pairs   <- ncol(combn(x = unique(ref_classes), m = 2, simplify = T))
probes_ovr <- 100
probes_ep <- 100
n_mrkrs <- sum(probes_ovr, probes_ep)

####


tmp <- perform.probe.select.by.probeset(ref_df = normalized_refs,
                                        ref_classes = ref_classes,
                                        probe_sets_for_selection = probesets,
                                        max_p_val = probe_sel_max_p,
                                        n_markers = n_mrkrs)

probesets_by_rna_type <- lapply(X = targ_probes_by_type,
                                FUN = function(ps){intersect(ps, probesets$univ)}) 

names(probesets_by_rna_type) <- names(targ_probes_by_type)
tmp2 <- perform.probe.select.by.probeset(ref_df = normalized_refs,
                                        ref_classes = ref_classes,
                                        probe_sets_for_selection = probesets_by_rna_type,
                                        max_p_val = probe_sel_max_p,
                                        n_markers = n_mrkrs)
####
    

informative_probes_obj <- perform.probe.selection(ref_probes =  probe_sel_refs,
                                              ref_classes = ref_classes,
                                              n_markers = n_mrkrs,
                                              max_p_val = probe_sel_max_p)

#selected_probes <- select.best.probes(informative_probes, probes_per_comparison = 15)

tmp1 <- lapply(c(1:n_classes),
               function(x){head(sort(informative_probes_obj$one.vs.rest$t_test$p.values[,x]), probes_ovr)})

tmp2 <- unique(names(unlist(tmp1)))

informative_probes1 <- tmp2

tmp11 <- lapply(c(1:n_pairs),
                function(x){head(sort(informative_probes_obj$each.pair$t_test$p.values[,x]), probes_ep)})

tmp22 <- unique(names(unlist(tmp11)))

informative_probes2 <- tmp22
chosen_probes <- unlist(union(informative_probes1, informative_probes2))

################################################################################

################################################################################
## PRODUCE A HEATMEAP TO SHOW CLUSTERING OF REFERECNCE SAMPLES
## BY CORRELATION OVER CHOSEN PROBES
################################################################################
## SET PLOTTING PARAMETERS
if(T){
    factor_vars_list <- list(subcell_fraction = ref_classes,
                             cell_type = meta_labels[, 3])
    
    factor_colors <- lapply(X = factor_vars_list, FUN = mk.col)  
    
    col_colors <- factor_colors$subcell_fraction
    row_colors <- factor_colors$cell_type
    
    #Function to set heatmap parameters {cbreaks heatmin heatmax colorstyle heatmargins, col_colors, row_colors }
    cbreaks <- 20
    heat_min <- 0
    heat_max <- 1
    heat_range <- (heat_max - heat_min)
    col_step <- heat_range / cbreaks
    heat_breaks <- seq(heat_min, heat_max, col_step)
    n_colors    <- length(heat_breaks) - 1
    heat_colors <- list(red_green = redgreen(n_colors), col_grad = colorRampPalette(c("white","steelblue"))(n_colors))
    heat_colors <- heat_colors[[1]]
    heat_margins <- c(15,15)
    
}

optim_ref_set <- probe_sel_refs[chosen_probes,]
#optim_ref_set <- probe_sel_refs[murillo_probes,]
#optim_ref_set <- probe_sel_refs[sample(x = 1:dim(probe_sel_refs)[[1]], size = 100),]
#rand_probes <- rownames(optim_ref_set)

colnames(optim_ref_set) <- paste(meta_labels[,1], meta_labels[,2],sep = '_')
colnames(optim_ref_set) <- paste(meta_labels[,1], meta_labels[,3],sep = '_')

expression_correlations_refs <- list(test = cor(optim_ref_set))

pdf(paste(figs_dir,"hmap_refs_only_univ_probes.pdf", sep = ''))
sapply(X = names(expression_correlations_refs),
       FUN = function(probeset){
           cors <- expression_correlations_refs
           main = probeset
           print(probeset)
           heatmap.2(x = cors[[probeset]], trace="none", col = heat_colors,
                     breaks = heat_breaks, margins = heat_margins,
                     ColSideColors = col_colors, RowSideColors = row_colors,
                     main = main, key.xlab = "Expression correlation" )
           
           
       })
dev.off()

## MAKE METHOD 'PRINT.PARAMS' THAT CAN SOMEHOW ENCAPSULATE THE SETTINGS THAT PRODUCE A GIVEN PLOT
## USE THE DO.CALL METHOD ON AN OBJECT THAT HOLDS NAMES FOR ALL THE SALIENT PARAMETERS FOR THE PLOT
## AS INPUT TO EACH
##  TAKES:  OUTFILE NAMES LIST;
##          PROBE SELECTION PARAMETERS [MAXPVAL, PROBESET, REF_CLASSES];
##          PLOT TITLE; ...

## WHEN PRODUCING A HEATMAP WITH SIDE
## COLORS MAKE A LEGEND MACHING COLORS TO LABELS
#color.legend(-.5,.2,-.1,.2,legend = c("a", "M"),
#rect.col = unique(factor_colors$subcell_fraction))

####################################################################################################

####################################################################################################
## EDEC STAGE 1 TEST ON REFERENCE DATA
####################################################################################################

## SET EDEC PARAMETERS ##
n_ct <- 5

optim_ref_set <- probe_sel_refs[chosen_probes,]
my_probes <- chosen_probes

stage1_refs <- run_edec_stage_1(meth_bulk_samples = optim_ref_set,
                                informative_loci = my_probes,
                                num_cell_types = n_ct)


check_stage1_results(stage1_result = stage1_refs, probeset = my_probes,
                     sample_meta_colors = factor_colors$subcell_fraction,
                     reference_profiles = optim_ref_set,
                     ref_colors = factor_colors$subcell_fraction)

####################################################################################################
## EDEC STAGE 1 : ON EXRNA ATLAS DATA with reference profiles
####################################################################################################
## OVER ALL PROBES ##
if(F){
refs_spike_in <- sample(x = colnames(normalized_refs),size = 20,replace = F)
refs_spike_in <- normalized_refs[,refs_spike_in]

#some NA values are present: "Neptune_Urine_33.11" "Neptune_Urine_33.24"
good_samps <- grep("Neptune_Urine_33.11|Neptune_Urine_33.24",
                   colnames(normalized_targ),invert = T, value = T)

normalized_targ <- normalized_targ[,good_samps]

smp <- sample(x = colnames(normalized_targ),size = round(x = (.8*(dim(normalized_targ)[[2]]))) )

normalized_targ <- normalized_targ[, smp]

exp_w_refs <- cbind(normalized_targ, refs_spike_in)


stage1_exp_w_refs = EDec::run_edec_stage_1(meth_bulk_samples = exp_w_refs,
                                           informative_loci = my_probes,
                                           num_cell_types = n_ct)

colnames(stage1_exp_w_refs$methylation) <- LETTERS[seq(n_ct)]

##IMPROVE THE CHECK STAGE 1 RESULTS FUNCTION TO RETURN QUANTITATIVE INFORMATION
## IN A LIST OBJECT
check_stage1_results(stage1_result = stage1_exp_w_refs,
                     probeset = my_probes,
                     show_props = T,
                     #These should be colors for target metadata
                     sample_meta_colors = rainbow(ncol(exp_w_refs)), 
                     reference_profiles = normalized_refs,
                     ref_colors =  factor_colors$subcell_fraction)
}
## OVER CHOSEN PROBES ONLY ##
my_probes <- murillo_probes
my_probes <- chosen_probes
exp_w_refs_chopro <- cbind(normalized_targ, refs_spike_in)[my_probes, ]

stage1_exp_w_refs_chopro = EDec::run_edec_stage_1(meth_bulk_samples = exp_w_refs_chopro,
                                           informative_loci = my_probes,
                                           num_cell_types = n_ct)

colnames(stage1_exp_w_refs$methylation) <- LETTERS[seq(n_ct)]

check_stage1_results(stage1_result = stage1_exp_w_refs_chopro,
                     probeset = my_probes,
                     show_props = T,
                     sample_meta_colors = rainbow(ncol(exp_w_refs_chopro)), 
                     reference_profiles = normalized_refs[my_probes,],
                     ref_colors =  factor_colors$subcell_fraction)

View(stage1_exp_w_refs_chopro$methylation)

intersect(chosen_probes, murillo_probes)


###########################
## Normalization - allncRNA
###########################
## MERGE REFERENCE AND TARGET SETS
## BEFORE NORMALIZATION???
input_allncRNA = (MERGE atlas AND lasser) 
input_allncRNA_.QN = quantile_normalisation(input.allncRNA)
## Plot 1
raw.allncRNA = 
    QN.allncRNA =
    
    
    pdf(paste(figs_dir,"Raw_vs_QN_Atlas_allncRNA.pdf", sep = ''))

plot(raw.allncRNA, QN.allncRNA)
dev.off()
## Transformations with alpha
input.allncRNA.trans.QN.logistic = t(apply(t(input.allncRNA.trans.QN),2, logistic))
input.allncRNA.trans.QN.logistic.nonZero = input.allncRNA.trans.QN.logistic[rowSums(input.allncRNA.trans.QN.logistic[, -1]) > 0, ]
input.allncRNA.trans.QN.logistic.nonZero.max = (1/max(input.allncRNA.trans.QN.logistic.nonZero)) * input.allncRNA.trans.QN.logistic.nonZero
## Plot 2
Logist.allncRNA = unlist(input.allncRNA.trans.QN.logistic.nonZero.max)
pdf("~/Documents/BCM/Lab/LouiseData/EDec_Analysis_Laurent_09/Stage_0/Raw_vs_QNandLogistic_Atlas_allncRNA.pdf")
plot(raw.allncRNA, Logist.allncRNA)
dev.off()

##############
## Subset Data
##############



#################################
## EDec - DWONG1, Saliva, Healthy
#################################
DWONG1 = c(1:length(DWONG1.Saliva.ncRNA.data))
Set1.DWONG1 = sample(x = DWONG1, size = length(DWONG1)/2, replace = FALSE)
Set2.DWONG1 = DWONG1[!DWONG1 %in% Set1.DWONG1]

DWONG1.Saliva.ncRNA.data.Set1 = DWONG1.Saliva.ncRNA.data[,Set1.DWONG1]
DWONG1.Saliva.ncRNA.data.Set2 = DWONG1.Saliva.ncRNA.data[,Set2.DWONG1]

Samples.DWONG1 = cor(DWONG1.Saliva.ncRNA.data.Set1, DWONG1.Saliva.ncRNA.data.Set2)
png("~/Documents/BCM/Lab/LouiseData/EDec_Analysis_Laurent_09/Stage_1/Correlation.chosenProbes.DWONG1.Saliva.Samples.png",1000,1000)
heatmap.2(as.matrix(Samples.DWONG1),
          trace = "none",
          col=my_palette,
          breaks = col_breaks,
          cexCol = 4,
          Colv = TRUE,
          margins = c(10,10))
dev.off()

EDec.DWONG1.ncRNA.Set1 = EDecStage1(methMixtureSamples = DWONG1.Saliva.ncRNA.data.Set1, cellTypeSpecificLoci = chosenProbes.allncRNA, nCts = 3, maxIts=2000, rssDiffStop=1e-10)
EDec.DWONG1.ncRNA.Set1.methy = EDec.DWONG1.ncRNA.Set1$methylation
EDec.DWONG1.ncRNA.Set1.methy.chosenProbes = EDec.DWONG1.ncRNA.Set1.methy[chosenProbes.allncRNA,]
EDec.DWONG1.ncRNA.Set1.prop = as.data.frame(t(EDec.DWONG1.ncRNA.Set1$proportions))

CorMatrix.DWONG1.ncRNA.Set1 = cor(EDec.DWONG1.ncRNA.Set1.methy.chosenProbes, Lasser.ncRNA.data.chosenProbes)
png("~/Documents/BCM/Lab/LouiseData/EDec_Analysis_Laurent_09/Stage_1/Profiles.DWONG1.Saliva.chosenProbes.ncRNA.3.Set1.png",1000,1000)
heatmap.2(as.matrix(t(CorMatrix.DWONG1.ncRNA.Set1)),
          trace = "none",
          col=my_palette,
          breaks = col_breaks,
          cexCol = 4,
          Colv = TRUE,
          margins = c(10,10))
dev.off()

png("~/Documents/BCM/Lab/LouiseData/EDec_Analysis_Laurent_09/Stage_1/Proportions.DWONG1.Saliva.chosenProbes.ncRNA.3.Set1.png",1000,1000)
heatmap.2(as.matrix(EDec.DWONG1.ncRNA.Set1.prop),
          trace="none",
          col=my_palette,
          breaks = col_breaks,
          labCol = FALSE,
          cexRow = 4,
          margins=c(10,10))
dev.off()

EDec.DWONG1.ncRNA.Set2 = EDecStage1(methMixtureSamples = DWONG1.Saliva.ncRNA.data.Set2, cellTypeSpecificLoci = chosenProbes.allncRNA, nCts = 3, maxIts=2000, rssDiffStop=1e-10)
EDec.DWONG1.ncRNA.Set2.methy = EDec.DWONG1.ncRNA.Set2$methylation
EDec.DWONG1.ncRNA.Set2.methy.chosenProbes = EDec.DWONG1.ncRNA.Set2.methy[chosenProbes.allncRNA,]
EDec.DWONG1.ncRNA.Set2.prop = as.data.frame(t(EDec.DWONG1.ncRNA.Set2$proportions))

CorMatrix.DWONG1.ncRNA.Set2 = cor(EDec.DWONG1.ncRNA.Set2.methy.chosenProbes, Lasser.ncRNA.data.chosenProbes)
png("~/Documents/BCM/Lab/LouiseData/EDec_Analysis_Laurent_09/Stage_1/Profiles.DWONG1.Saliva.chosenProbes.ncRNA.3.Set2.png",1000,1000)
heatmap.2(as.matrix(t(CorMatrix.DWONG1.ncRNA.Set2)),
          trace = "none",
          col=my_palette,
          breaks = col_breaks,
          cexCol = 4,
          Colv = TRUE,
          margins = c(10,10))
dev.off()

png("~/Documents/BCM/Lab/LouiseData/EDec_Analysis_Laurent_09/Stage_1/Proportions.DWONG1.Saliva.chosenProbes.ncRNA.3.Set2.png",1000,1000)
heatmap.2(as.matrix(EDec.DWONG1.ncRNA.Set2.prop),
          trace="none",
          col=my_palette,
          breaks = col_breaks,
          labCol = FALSE,
          cexRow = 4,
          margins=c(10,10))
dev.off()

EDec.DWONG1.ncRNA = EDecStage1(methMixtureSamples = DWONG1.Saliva.ncRNA.data, cellTypeSpecificLoci = chosenProbes.allncRNA, nCts = 3, maxIts=2000, rssDiffStop=1e-10)
EDec.DWONG1.ncRNA.methy = EDec.DWONG1.ncRNA$methylation
EDec.DWONG1.ncRNA.methy.chosenProbes = EDec.DWONG1.ncRNA.methy[chosenProbes.allncRNA,]
EDec.DWONG1.ncRNA.prop = as.data.frame(t(EDec.DWONG1.ncRNA$proportions))

CorMatrix.DWONG1.ncRNA = cor(EDec.DWONG1.ncRNA.methy.chosenProbes, Lasser.ncRNA.data.chosenProbes)
png("~/Documents/BCM/Lab/LouiseData/EDec_Analysis_Laurent_09/Stage_1/Profiles.DWONG1.Saliva.chosenProbes.ncRNA.3.All.png",1000,1000)
heatmap.2(as.matrix(t(CorMatrix.DWONG1.ncRNA)),
          trace = "none",
          col=my_palette,
          breaks = col_breaks,
          cexCol = 4,
          Colv = TRUE,
          margins = c(10,10))
dev.off()

png("~/Documents/BCM/Lab/LouiseData/EDec_Analysis_Laurent_09/Stage_1/Proportions.DWONG1.Saliva.chosenProbes.ncRNA.3.All.png",1000,1000)
heatmap.2(as.matrix(EDec.DWONG1.ncRNA.prop),
          trace="none",
          col=my_palette,
          breaks = col_breaks,
          labCol = FALSE,
          cexRow = 4,
          margins=c(10,10))
dev.off()


#####################################
## Name all predicted samples (by hand??)
#####################################

colnames(EDec.KJENS1AP.Serum.ncRNA.methy) = c("LD","Other","HD")
colnames(EDec.KJENS1AP.Serum.ncRNA.methy.chosenProbes) = c("LD","Other","HD")
colnames(EDec.KJENS1AP.Serum.ncRNA.Set1.methy.chosenProbes) = c("HD","LD","Other")
colnames(EDec.KJENS1AP.Serum.ncRNA.Set2.methy.chosenProbes) = c("HD","Other","LD")

colnames(EDec.KJENS1RID.Saliva.ncRNA.methy) = c("HD","Other","LD")
colnames(EDec.KJENS1RID.Saliva.ncRNA.methy.chosenProbes) = c("HD","Other","LD")
colnames(EDec.KJENS1RID.Saliva.ncRNA.Set1.methy.chosenProbes) = c("Other","LD","HD")
colnames(EDec.KJENS1RID.Saliva.ncRNA.Set2.methy.chosenProbes) = c("Other","HD","LD")

colnames(EDec.KJENS1RID.Urine.ncRNA.methy) = c("LD","Other","HD")
colnames(EDec.KJENS1RID.Urine.ncRNA.methy.chosenProbes) = c("LD","Other","HD")
colnames(EDec.KJENS1RID.Urine.ncRNA.Set1.methy.chosenProbes) = c("LD","Other","HD")
colnames(EDec.KJENS1RID.Urine.ncRNA.Set2.methy.chosenProbes) = c("HD","LD","Other")

all.species = c(colnames(EDec.KJENS1AP.Serum.ncRNA.methy), colnames(EDec.KJENS1RID.Saliva.ncRNA.methy), colnames(EDec.KJENS1RID.Urine.ncRNA.methy))

######################
## Correlations - Sets
######################

Sets.KJENS1RID.Urine = round(cor(EDec.KJENS1RID.Urine.ncRNA.Set1.methy.chosenProbes, EDec.KJENS1RID.Urine.ncRNA.Set2.methy.chosenProbes, use="complete"), digits=2)
png("~/Documents/BCM/Lab/LouiseData/EDec_Analysis_Laurent_09/Stage_1/1.Correlation.chosenProbes.KJENS1RID.Urine.Sets.png",1000,1000)
qplot(x=Var1, y=Var2, data=melt(Sets.KJENS1RID.Urine), fill=value, geom="tile") + theme(axis.text = element_text(size = 24)) + geom_text(aes(Var2, Var1, label = value), color = "black", size = 8) + scale_fill_gradient2(low="white", high="darkgreen")
dev.off()

########################
## Correlation - Studies
########################
## all Probes
Sets.KJENS1AP.Serum.KJENS1RID.Saliva = round(cor(EDec.KJENS1AP.Serum.ncRNA.methy, EDec.KJENS1RID.Saliva.ncRNA.methy, use="complete"), digits=2)
png("~/Documents/BCM/Lab/LouiseData/EDec_Analysis_Laurent_09/Stage_1/2.a.Correlation.KJENS1AP.Serum.KJENS1RID.Saliva.png",1000,1000)
qplot(x=Var1, y=Var2, data=melt(Sets.KJENS1AP.Serum.KJENS1RID.Saliva), fill=value, geom="tile") + theme(axis.text = element_text(size = 24)) + geom_text(aes(Var2, Var1, label = value), color = "black", size = 8) + scale_fill_gradient2(low="white", high="darkgreen")
dev.off()

## chosen Probes
Sets.KJENS1AP.Serum.KJENS1RID.Saliva = round(cor(EDec.KJENS1AP.Serum.ncRNA.methy.chosenProbes, EDec.KJENS1RID.Saliva.ncRNA.methy.chosenProbes, use="complete"), digits=2)
png("~/Documents/BCM/Lab/LouiseData/EDec_Analysis_Laurent_09/Stage_1/2.b.Correlation.chosenProbes.KJENS1AP.Serum.KJENS1RID.Saliva.png",1000,1000)
qplot(x=Var1, y=Var2, data=melt(Sets.KJENS1AP.Serum.KJENS1RID.Saliva), fill=value, geom="tile") + theme(axis.text = element_text(size = 24)) + geom_text(aes(Var2, Var1, label = value), color = "black", size = 8) + scale_fill_gradient2(low="white", high="darkgreen")
dev.off()

####################################
## Clustering/PCA for all references
####################################
colnames(EDec.KJENS1AP.Serum.ncRNA.methy) = c("LD_KJENS1AP_Serum","Other_KJENS1AP_Serum","HD_KJENS1AP_Serum")
colnames(EDec.KJENS1AP.Serum.ncRNA.methy.chosenProbes) = c("LD_KJENS1AP_Serum","Other_KJENS1AP_Serum","HD_KJENS1AP_Serum")

colnames(EDec.KJENS1RID.Saliva.ncRNA.methy) = c("HD_KJENS1RID_Saliva","Other_KJENS1RID_Saliva","LD_KJENS1RID_Saliva")
colnames(EDec.KJENS1RID.Saliva.ncRNA.methy.chosenProbes) = c("HD_KJENS1RID_Saliva","Other_KJENS1RID_Saliva","LD_KJENS1RID_Saliva")

colnames(EDec.KJENS1RID.Urine.ncRNA.methy) = c("LD_KJENS1RID_Urine","Other_KJENS1RID_Urine","HD_KJENS1RID_Urine")
colnames(EDec.KJENS1RID.Urine.ncRNA.methy.chosenProbes) = c("LD_KJENS1RID_Urine","Other_KJENS1RID_Urine","HD_KJENS1RID_Urine")

pca.input.1 = cbind(EDec.KJENS1AP.Serum.ncRNA.methy, EDec.KJENS1RID.Saliva.ncRNA.methy, EDec.KJENS1RID.Urine.ncRNA.methy)

all.pca.1 = prcomp(t(pca.input.1), center = TRUE)
all.plot.1 = ggbiplot(all.pca.1, obs.scale = 1, var.scale = 1, groups = all.species, ellipse = TRUE, circle = TRUE, var.axes = FALSE) 
all.plot.1 = all.plot.1 + scale_color_discrete(name = '')
all.plot.1 = all.plot.1 + geom_point(shape = 1, size = 5)
all.plot.1 = all.plot.1 + theme(legend.direction = 'horizontal', legend.position = 'top', legend.text = element_text(size = 20), axis.text = element_text(colour = "black", size = 20))

png("~/Documents/BCM/Lab/LouiseData/EDec_Analysis_Laurent_09/Stage_1/3.PCA.CarrierTypes.KJENS1AP.Serum.KJENS1RID.Saliva.KJENS1RID.Urine.png",1000,1000)
plot(all.plot.1)
dev.off()
remove(all.pca.1)

pca.input.2 = cbind(EDec.KJENS1AP.Serum.ncRNA.methy.chosenProbes, EDec.KJENS1RID.Saliva.ncRNA.methy.chosenProbes, EDec.KJENS1RID.Urine.ncRNA.methy.chosenProbes)

all.pca.2 = prcomp(t(pca.input.2), center = TRUE)
all.plot.2 = ggbiplot(all.pca.2, obs.scale = 1, var.scale = 1, groups = all.species, ellipse = TRUE, circle = TRUE, var.axes = FALSE) 
all.plot.2 = all.plot.2 + scale_color_discrete(name = '')
all.plot.2 = all.plot.2 + geom_point(shape = 1, size = 5)
all.plot.2 = all.plot.2 + theme(legend.direction = 'horizontal', legend.position = 'top', legend.text = element_text(size = 20), axis.text = element_text(colour = "black", size = 20))

png("~/Documents/BCM/Lab/LouiseData/EDec_Analysis_Laurent_09/Stage_1/3.PCA.CarrierTypes.KJENS1AP.Serum.KJENS1RID.Saliva.KJENS1RID.Urine.chosenProbes.png",1000,1000)
plot(all.plot.2)
dev.off()
remove(all.pca.2)
