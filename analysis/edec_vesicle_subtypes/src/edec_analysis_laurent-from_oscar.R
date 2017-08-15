###############################
##  Based on: EDec_OscarMurillo
###############################

################################################################################
##  LOAD PACKAGES & SOURCE FILES
################################################################################
rm(list=ls())

pkgs <- list('gplots', 'class', 'RColorBrewer', 'stringr', 'ggplot2', 'reshape2',
          'rms', 'devtools', 'clue', 'EDec', 'tidyverse')
lapply(pkgs, library, character.only = T)
session_info()

#library('ggbiplot')

#source("~/Documents/BCM/Lab/LouiseData/sources/PAM50ExerciseFunctions.R")
#source("~/Documents/BCM/Lab/LouiseData/sources/EDec.R")
source("../edec_vesicle_subtypes/src/K_cell_number.R")

################################################################################
##  FUNCTION DEFINITIONS
################################################################################
source("../edec_vesicle_subtypes/src/edec_vesicle_subtypes_functions.R")

################################################################################
##  LOAD DATA
################################################################################
################################################
##  READ IN METADATA
################################################
## EXRNA ATLAS METADATA 


################################################
##  READ IN EXCERPT READS PER MILLION OUTPUT DATA
################################################

######  FIRST GET FILENAMES AND PATHS FOR EXCERPT RESULT .RDATA FROM EACH STUDY  ######
## HD/LD DATA FROM LASSER ET AL:
## LAURENT DATA FROM KIT X CONDITIONED MEDIA EXPERIMENTS:
## GINGERAS CELLS AND SUBCELLULAR FRACTIONS
#path <- "../edec_vesicle_subtypes/input/hd_ld_atlas_gingeras_laurent"
path <- "../edec_vesicle_subtypes/input/lasser_atlas_laurent"

all_combined_rpm <-t(condense.post.proc.results(results_rdata_path = path,
                                               per_million = T,
                                               include_exogenous = F))

## PROCESS DIMNAMES (MAKE THIS A FUNCTION)
rownames(all_combined_rpm) <- gsub("sample_|_fastq","",rownames(all_combined_rpm))

gsub(pattern = "\\|.*$",replacement = "",
     x = colnames(all_combined_rpm)[grep("hsa_piR",colnames(all_combined_rpm))])


small_rna_keys <- paste("Gly|Ala|Leu|Met|Phe|Trp|Lys|Gln|Glu|Ser",
                        "Pro|Val|Ile|Cys|Tyr|His|Arg|Asn|Asp|Thr|Mt_tRNA",
                        "hsa-let|hsa-mi[Rr]|hsa_piR|snoRNA|snRNA|^Y_RNA",
                        #"hsa_circ|circ",
                        "3prime_overlapping|vaultRNA",
                        "miRNA|sRNA|scRNA|VTRNA",
                        sep = "|")

##############

dim(all_combined_rpm[,grep(small_rna_keys,colnames(all_combined_rpm))])

## EXRNA ATLAS DATA FOR ALL STUDIES:
path <- "./input/exrna_atlas/Atlas_allncRNA_RPM.txt"
atlas_rpm_all_nc <- t(read.delim(file = path, sep = "\t"))
dim(atlas_rpm_all_nc)
dim(atlas_rpm_all_nc[,grep(small_rna_keys,colnames(atlas_rpm_all_nc))])

if(F){ ## OSCAR'S APPROACH ##
    
Lasser_allncRNA_RPM <- read.delim("P:/brl/proj/exrna_atlas/analysis/edec_vesicle_subtypes/input/Oscar/input/input/Lasser_HD_LD/Lasser_allncRNA_RPM.txt", 
                                      header = T, sep = "\t", row.names = 1)
    
Atlas_allncRNA_RPM <- read.delim("P:/brl/proj/exrna_atlas/analysis/edec_vesicle_subtypes/input/Oscar/input/input/Atlas_Data_20170804/Atlas_allncRNA_RPM.txt", 
                                     header = T, sep = "\t", row.names = 1)
    
lapply(FUN = dim, X = list(Lasser_allncRNA_RPM, Atlas_allncRNA_RPM))    
    

common.allncRNA = intersect(row.names(Lasser.allncRNA.RPM), row.names(Atlas.allncRNA.RPM.common))

## GET THE RNA SPECIES COMMON TO ALL THREE DATASETS:
## MERGE ON THESE...


}

################################################################################
##  SET GLOBAL PARAMETERS
################################################################################

my_palette = colorRampPalette(c("white","lightgreen","darkgreen"))(n = 29)
col_breaks = c(seq(0.0000,0.3330,length=10), # for white
               seq(0.3331,0.6660,length=10), # for lightgreen
               seq(0.6661,1.0000,length=10)) # for darkgreen

my_palette = colorRampPalette(c("white","lightgreen","darkgreen"))(n = 29)
col_breaks_2 = c(seq(0.50,0.6699,length=10), # for white
                 seq(0.67,0.8399,length=10), # for lightgreen
                 seq(0.84,1.0000,length=10)) # for darkgreen



####################################
## SELECT RNA THAT CAN SEPARATE LD/HD
## SIGNATURES BASED ON EXPRESSION
####################################
## T-test
Lasser.subtypes = c("HD","HD","LD","LD")
colors.1 = rep("white",length(Lasser.subtypes))
colors.1[Lasser.subtypes=="HD"] = "gold"
colors.1[Lasser.subtypes=="LD"] = "violet"

tTestResultOneVsRest.ncRNA = perform_t_tests_all_classes_one_vs_rest(dataMatrix = Lasser.allncRNA.RPM, classVector = Lasser.subtypes)

chosenProbes = c()
for(i in 1:2){
  sigProbes = row.names(tTestResultOneVsRest.ncRNA$P.Values[tTestResultOneVsRest.ncRNA$P.Values[,i] < 0.005,])
  sigDiff = tTestResultOneVsRest.ncRNA$Difference.Between.Means[which(row.names(tTestResultOneVsRest.ncRNA$P.Values) %in% sigProbes), i]
  highDiffProbes = names(tail(sort(sigDiff),100))
  lowDiffProbes = names(head(sort(sigDiff),100))
  chosenProbes = c(chosenProbes, highDiffProbes, lowDiffProbes)
}
chosenProbes.allncRNA = unique(c(chosenProbes))
Lasser.allncRNA.RPM.chosenProbes = Lasser.allncRNA.RPM[chosenProbes.allncRNA,]

png("~/Documents/BCM/Lab/LouiseData/EDec_Analysis_Laurent_09/Stage_0/Clustering.chosenProbes.ncRNA.png",1000,1000)
heatmap.2(as.matrix(Lasser.allncRNA.RPM.chosenProbes),
          trace="none",
          col=my_palette,
          labCol = FALSE,
          margins=c(10,10))
dev.off()

ProbesFile = file("~/Documents/BCM/Lab/LouiseData/EDec_Analysis_Laurent_09/Stage_0/chosenProbes_ncRNA.txt")
writeLines(c(chosenProbes.allncRNA), ProbesFile)
close(ProbesFile)

###########################
## Normalization - allncRNA
###########################
input.allncRNA = merge.Lasser.Atlas.allncRNA
input.allncRNA.trans = t(input.allncRNA)
input.allncRNA.trans.QN = as.data.frame(t(quantile_normalisation(input.allncRNA.trans)))
## Plot 1
raw.allncRNA = unlist(input.allncRNA)
QN.allncRNA = unlist(input.allncRNA.trans.QN)
pdf("~/Documents/BCM/Lab/LouiseData/EDec_Analysis_Laurent_09/Stage_0/Raw_vs_QN_Atlas_allncRNA.pdf")
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
meta.Lasser = colnames(Lasser.allncRNA.RPM)
Lasser.ncRNA.data = input.allncRNA.trans.QN.logistic.nonZero.max[,meta.Lasser]
Lasser.ncRNA.data.chosenProbes = Lasser.ncRNA.data[chosenProbes.allncRNA,]

condition.1 = "Healthy Control"
biofluid.1 = "Serum"
biofluid.2 = "Plasma"
biofluid.3 = "Cerebrospinal fluid"
biofluid.4 = "Saliva"
biofluid.5 = "Urine"

study.1 = "DWONG1-gastric-cancer-and-controls-2016-10-25" #Saliva
study.2 = "JFREE1-Framingham-UH2-40Samples-2016-10-17" #Plasma
study.3 = "LLAUR1-placental-dysfunction-biomarkers-2016-10-13" #Serum
study.4 = "TPATE1-human-plasma-healthyVsCancer-2016-10-17" #Plasma
study.5 = "KJENS1-Alzheimers_Parkinsons-2016-10-17" # Serum, CSF
study.6 = "KJENS1-RIDProject-2016-06-27" #Plasma, Saliva, Urine

DWONG1.Saliva.temp = as.data.frame(t(Atlas.meta[which(Atlas.meta$biofluid_name == biofluid.4 & Atlas.meta$condition == condition.1 & Atlas.meta$Study == study.1),]))
DWONG1.Saliva = c(colnames(DWONG1.Saliva.temp))
DWONG1.Saliva.ncRNA.data = as.data.frame(input.allncRNA.trans.QN.logistic.nonZero.max[,DWONG1.Saliva])

JFREE1.Plasma.temp = as.data.frame(t(Atlas.meta[which(Atlas.meta$biofluid_name == biofluid.2 & Atlas.meta$condition == condition.1 & Atlas.meta$Study == study.2),]))
JFREE1.Plasma = c(colnames(JFREE1.Plasma.temp))
JFREE1.Plasma.ncRNA.data = as.data.frame(input.allncRNA.trans.QN.logistic.nonZero.max[,JFREE1.Plasma])

LLAUR1.Serum.temp = as.data.frame(t(Atlas.meta[which(Atlas.meta$biofluid_name == biofluid.1 & Atlas.meta$condition == condition.1 & Atlas.meta$Study == study.3),]))
LLAUR1.Serum = c(colnames(LLAUR1.Serum.temp))
LLAUR1.Serum.ncRNA.data = as.data.frame(input.allncRNA.trans.QN.logistic.nonZero.max[,LLAUR1.Serum])

TPATE1.Plasma.temp = as.data.frame(t(Atlas.meta[which(Atlas.meta$biofluid_name == biofluid.2 & Atlas.meta$condition == condition.1 & Atlas.meta$Study == study.4),]))
TPATE1.Plasma = c(colnames(TPATE1.Plasma.temp))
TPATE1.Plasma.ncRNA.data = as.data.frame(input.allncRNA.trans.QN.logistic.nonZero.max[,TPATE1.Plasma])

KJENS1AP.Serum.temp = as.data.frame(t(Atlas.meta[which(Atlas.meta$biofluid_name == biofluid.1 & Atlas.meta$condition == condition.1 & Atlas.meta$Study == study.5),]))
KJENS1AP.Serum = c(colnames(KJENS1AP.Serum.temp))
KJENS1AP.Serum.ncRNA.data = as.data.frame(input.allncRNA.trans.QN.logistic.nonZero.max[,KJENS1AP.Serum])

KJENS1AP.CSF.temp = as.data.frame(t(Atlas.meta[which(Atlas.meta$biofluid_name == biofluid.3 & Atlas.meta$condition == condition.1 & Atlas.meta$Study == study.5),]))
KJENS1AP.CSF = c(colnames(KJENS1AP.CSF.temp))
KJENS1AP.CSF.ncRNA.data = as.data.frame(input.allncRNA.trans.QN.logistic.nonZero.max[,KJENS1AP.CSF])

KJENS1RID.Plasma.temp = as.data.frame(t(Atlas.meta[which(Atlas.meta$biofluid_name == biofluid.2 & Atlas.meta$condition == condition.1 & Atlas.meta$Study == study.6),]))
KJENS1RID.Plasma = c(colnames(KJENS1RID.Plasma.temp))
KJENS1RID.Plasma.ncRNA.data = as.data.frame(input.allncRNA.trans.QN.logistic.nonZero.max[,KJENS1RID.Plasma])

KJENS1RID.Saliva.temp = as.data.frame(t(Atlas.meta[which(Atlas.meta$biofluid_name == biofluid.4 & Atlas.meta$condition == condition.1 & Atlas.meta$Study == study.6),]))
KJENS1RID.Saliva = c(colnames(KJENS1RID.Saliva.temp))
KJENS1RID.Saliva.ncRNA.data = as.data.frame(input.allncRNA.trans.QN.logistic.nonZero.max[,KJENS1RID.Saliva])

KJENS1RID.Urine.temp = as.data.frame(t(Atlas.meta[which(Atlas.meta$biofluid_name == biofluid.5 & Atlas.meta$condition == condition.1 & Atlas.meta$Study == study.6),]))
KJENS1RID.Urine = c(colnames(KJENS1RID.Urine.temp))
KJENS1RID.Urine.ncRNA.data = as.data.frame(input.allncRNA.trans.QN.logistic.nonZero.max[,KJENS1RID.Urine])

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

#################################
## EDec - JFREE1, Plasma, Healthy
#################################
JFREE1 = c(1:length(JFREE1.Plasma.ncRNA.data))
Set1.JFREE1 = sample(x = JFREE1, size = length(JFREE1)/2, replace = FALSE)
Set2.JFREE1 = JFREE1[!JFREE1 %in% Set1.JFREE1]

JFREE1.Plasma.ncRNA.data.Set1 = JFREE1.Plasma.ncRNA.data[,Set1.JFREE1]
JFREE1.Plasma.ncRNA.data.Set2 = JFREE1.Plasma.ncRNA.data[,Set2.JFREE1]

Samples.JFREE1 = cor(JFREE1.Plasma.ncRNA.data.Set1, JFREE1.Plasma.ncRNA.data.Set2)
png("~/Documents/BCM/Lab/LouiseData/EDec_Analysis_Laurent_09/Stage_1/Correlation.chosenProbes.JFREE1.Plasma.Samples.png",1000,1000)
heatmap.2(as.matrix(Samples.JFREE1),
          trace = "none",
          col=my_palette,
          breaks = col_breaks,
          cexCol = 4,
          Colv = TRUE,
          margins = c(10,10))
dev.off()

EDec.JFREE1.ncRNA.Set1 = EDecStage1(methMixtureSamples = JFREE1.Plasma.ncRNA.data.Set1, cellTypeSpecificLoci = chosenProbes.allncRNA, nCts = 3, maxIts=2000, rssDiffStop=1e-10)
EDec.JFREE1.ncRNA.Set1.methy = EDec.JFREE1.ncRNA.Set1$methylation
EDec.JFREE1.ncRNA.Set1.methy.chosenProbes = EDec.JFREE1.ncRNA.Set1.methy[chosenProbes.allncRNA,]
EDec.JFREE1.ncRNA.Set1.prop = as.data.frame(t(EDec.JFREE1.ncRNA.Set1$proportions))

CorMatrix.JFREE1.ncRNA.Set1 = cor(EDec.JFREE1.ncRNA.Set1.methy.chosenProbes, Lasser.ncRNA.data.chosenProbes)
png("~/Documents/BCM/Lab/LouiseData/EDec_Analysis_Laurent_09/Stage_1/Profiles.JFREE1.Plasma.chosenProbes.ncRNA.3.Set1.png",1000,1000)
heatmap.2(as.matrix(t(CorMatrix.JFREE1.ncRNA.Set1)),
          trace = "none",
          col=my_palette,
          breaks = col_breaks,
          cexCol = 4,
          Colv = TRUE,
          margins = c(10,10))
dev.off()

png("~/Documents/BCM/Lab/LouiseData/EDec_Analysis_Laurent_09/Stage_1/Proportions.JFREE1.Plasma.chosenProbes.ncRNA.3.Set1.png",1000,1000)
heatmap.2(as.matrix(EDec.JFREE1.ncRNA.Set1.prop),
          trace="none",
          col=my_palette,
          breaks = col_breaks,
          labCol = FALSE,
          cexRow = 4,
          margins=c(10,10))
dev.off()

EDec.JFREE1.ncRNA.Set2 = EDecStage1(methMixtureSamples = JFREE1.Plasma.ncRNA.data.Set2, cellTypeSpecificLoci = chosenProbes.allncRNA, nCts = 3, maxIts=2000, rssDiffStop=1e-10)
EDec.JFREE1.ncRNA.Set2.methy = EDec.JFREE1.ncRNA.Set2$methylation
EDec.JFREE1.ncRNA.Set2.methy.chosenProbes = EDec.JFREE1.ncRNA.Set2.methy[chosenProbes.allncRNA,]
EDec.JFREE1.ncRNA.Set2.prop = as.data.frame(t(EDec.JFREE1.ncRNA.Set2$proportions))

CorMatrix.JFREE1.ncRNA.Set2 = cor(EDec.JFREE1.ncRNA.Set2.methy.chosenProbes, Lasser.ncRNA.data.chosenProbes)
png("~/Documents/BCM/Lab/LouiseData/EDec_Analysis_Laurent_09/Stage_1/Profiles.JFREE1.Plasma.chosenProbes.ncRNA.3.Set2.png",1000,1000)
heatmap.2(as.matrix(t(CorMatrix.JFREE1.ncRNA.Set2)),
          trace = "none",
          col=my_palette,
          breaks = col_breaks,
          cexCol = 4,
          Colv = TRUE,
          margins = c(10,10))
dev.off()

png("~/Documents/BCM/Lab/LouiseData/EDec_Analysis_Laurent_09/Stage_1/Proportions.JFREE1.Plasma.chosenProbes.ncRNA.3.Set2.png",1000,1000)
heatmap.2(as.matrix(EDec.JFREE1.ncRNA.Set2.prop),
          trace="none",
          col=my_palette,
          breaks = col_breaks,
          labCol = FALSE,
          cexRow = 4,
          margins=c(10,10))
dev.off()

EDec.JFREE1.ncRNA = EDecStage1(methMixtureSamples = JFREE1.Plasma.ncRNA.data, cellTypeSpecificLoci = chosenProbes.allncRNA, nCts = 3, maxIts=2000, rssDiffStop=1e-10)
EDec.JFREE1.ncRNA.methy = EDec.JFREE1.ncRNA$methylation
EDec.JFREE1.ncRNA.methy.chosenProbes = EDec.JFREE1.ncRNA.methy[chosenProbes.allncRNA,]
EDec.JFREE1.ncRNA.prop = as.data.frame(t(EDec.JFREE1.ncRNA$proportions))

CorMatrix.JFREE1.ncRNA = cor(EDec.JFREE1.ncRNA.methy.chosenProbes, Lasser.ncRNA.data.chosenProbes)
png("~/Documents/BCM/Lab/LouiseData/EDec_Analysis_Laurent_09/Stage_1/Profiles.JFREE1.Plasma.chosenProbes.ncRNA.3.All.png",1000,1000)
heatmap.2(as.matrix(t(CorMatrix.JFREE1.ncRNA)),
          trace = "none",
          col=my_palette,
          breaks = col_breaks,
          cexCol = 4,
          Colv = TRUE,
          margins = c(10,10))
dev.off()

png("~/Documents/BCM/Lab/LouiseData/EDec_Analysis_Laurent_09/Stage_1/Proportions.JFREE1.Plasma.chosenProbes.ncRNA.3.All.png",1000,1000)
heatmap.2(as.matrix(EDec.JFREE1.ncRNA.prop),
          trace="none",
          col=my_palette,
          breaks = col_breaks,
          labCol = FALSE,
          cexRow = 4,
          margins=c(10,10))
dev.off()

#################################
## EDec - LLAUR1, Serum, Healthy
#################################
LLAUR1 = c(1:length(LLAUR1.Serum.ncRNA.data))
Set1.LLAUR1 = sample(x = LLAUR1, size = length(LLAUR1)/2, replace = FALSE)
Set2.LLAUR1 = LLAUR1[!LLAUR1 %in% Set1.LLAUR1]

LLAUR1.Serum.ncRNA.data.Set1 = LLAUR1.Serum.ncRNA.data[,Set1.LLAUR1]
LLAUR1.Serum.ncRNA.data.Set2 = LLAUR1.Serum.ncRNA.data[,Set2.LLAUR1]

Samples.LLAUR1 = cor(LLAUR1.Serum.ncRNA.data.Set1, LLAUR1.Serum.ncRNA.data.Set2)
png("~/Documents/BCM/Lab/LouiseData/EDec_Analysis_Laurent_09/Stage_1/Correlation.chosenProbes.LLAUR1.Serum.Samples.png",1000,1000)
heatmap.2(as.matrix(Samples.LLAUR1),
          trace = "none",
          col=my_palette,
          breaks = col_breaks,
          cexCol = 4,
          Colv = TRUE,
          margins = c(10,10))
dev.off()

EDec.LLAUR1.ncRNA.Set1 = EDecStage1(methMixtureSamples = LLAUR1.Serum.ncRNA.data.Set1, cellTypeSpecificLoci = chosenProbes.allncRNA, nCts = 3, maxIts=2000, rssDiffStop=1e-10)
EDec.LLAUR1.ncRNA.Set1.methy = EDec.LLAUR1.ncRNA.Set1$methylation
EDec.LLAUR1.ncRNA.Set1.methy.chosenProbes = EDec.LLAUR1.ncRNA.Set1.methy[chosenProbes.allncRNA,]
EDec.LLAUR1.ncRNA.Set1.prop = as.data.frame(t(EDec.LLAUR1.ncRNA.Set1$proportions))

CorMatrix.LLAUR1.ncRNA.Set1 = cor(EDec.LLAUR1.ncRNA.Set1.methy.chosenProbes, Lasser.ncRNA.data.chosenProbes)
png("~/Documents/BCM/Lab/LouiseData/EDec_Analysis_Laurent_09/Stage_1/Profiles.LLAUR1.Serum.chosenProbes.ncRNA.3.Set1.png",1000,1000)
heatmap.2(as.matrix(t(CorMatrix.LLAUR1.ncRNA.Set1)),
          trace = "none",
          col=my_palette,
          breaks = col_breaks,
          cexCol = 4,
          Colv = TRUE,
          margins = c(10,10))
dev.off()

png("~/Documents/BCM/Lab/LouiseData/EDec_Analysis_Laurent_09/Stage_1/Proportions.LLAUR1.Serum.chosenProbes.ncRNA.3.Set1.png",1000,1000)
heatmap.2(as.matrix(EDec.LLAUR1.ncRNA.Set1.prop),
          trace="none",
          col=my_palette,
          breaks = col_breaks,
          labCol = FALSE,
          cexRow = 4,
          margins=c(10,10))
dev.off()

EDec.LLAUR1.ncRNA.Set2 = EDecStage1(methMixtureSamples = LLAUR1.Serum.ncRNA.data.Set2, cellTypeSpecificLoci = chosenProbes.allncRNA, nCts = 3, maxIts=2000, rssDiffStop=1e-10)
EDec.LLAUR1.ncRNA.Set2.methy = EDec.LLAUR1.ncRNA.Set2$methylation
EDec.LLAUR1.ncRNA.Set2.methy.chosenProbes = EDec.LLAUR1.ncRNA.Set2.methy[chosenProbes.allncRNA,]
EDec.LLAUR1.ncRNA.Set2.prop = as.data.frame(t(EDec.LLAUR1.ncRNA.Set2$proportions))

CorMatrix.LLAUR1.ncRNA.Set2 = cor(EDec.LLAUR1.ncRNA.Set2.methy.chosenProbes, Lasser.ncRNA.data.chosenProbes)
png("~/Documents/BCM/Lab/LouiseData/EDec_Analysis_Laurent_09/Stage_1/Profiles.LLAUR1.Serum.chosenProbes.ncRNA.3.Set2.png",1000,1000)
heatmap.2(as.matrix(t(CorMatrix.LLAUR1.ncRNA.Set2)),
          trace = "none",
          col=my_palette,
          breaks = col_breaks,
          cexCol = 4,
          Colv = TRUE,
          margins = c(10,10))
dev.off()

png("~/Documents/BCM/Lab/LouiseData/EDec_Analysis_Laurent_09/Stage_1/Proportions.LLAUR1.Serum.chosenProbes.ncRNA.3.Set2.png",1000,1000)
heatmap.2(as.matrix(EDec.LLAUR1.ncRNA.Set2.prop),
          trace="none",
          col=my_palette,
          breaks = col_breaks,
          labCol = FALSE,
          cexRow = 4,
          margins=c(10,10))
dev.off()

EDec.LLAUR1.ncRNA = EDecStage1(methMixtureSamples = LLAUR1.Serum.ncRNA.data, cellTypeSpecificLoci = chosenProbes.allncRNA, nCts = 3, maxIts=2000, rssDiffStop=1e-10)
EDec.LLAUR1.ncRNA.methy = EDec.LLAUR1.ncRNA$methylation
EDec.LLAUR1.ncRNA.methy.chosenProbes = EDec.LLAUR1.ncRNA.methy[chosenProbes.allncRNA,]
EDec.LLAUR1.ncRNA.prop = as.data.frame(t(EDec.LLAUR1.ncRNA$proportions))

CorMatrix.LLAUR1.ncRNA = cor(EDec.LLAUR1.ncRNA.methy.chosenProbes, Lasser.ncRNA.data.chosenProbes)
png("~/Documents/BCM/Lab/LouiseData/EDec_Analysis_Laurent_09/Stage_1/Profiles.LLAUR1.Serum.chosenProbes.ncRNA.3.All.png",1000,1000)
heatmap.2(as.matrix(t(CorMatrix.LLAUR1.ncRNA)),
          trace = "none",
          col=my_palette,
          breaks = col_breaks,
          cexCol = 4,
          Colv = TRUE,
          margins = c(10,10))
dev.off()

png("~/Documents/BCM/Lab/LouiseData/EDec_Analysis_Laurent_09/Stage_1/Proportions.LLAUR1.Serum.chosenProbes.ncRNA.3.All.png",1000,1000)
heatmap.2(as.matrix(EDec.LLAUR1.ncRNA.prop),
          trace="none",
          col=my_palette,
          breaks = col_breaks,
          labCol = FALSE,
          cexRow = 4,
          margins=c(10,10))
dev.off()

#################################
## EDec - TPATE1, Plasma, Healthy
#################################
TPATE1 = c(1:length(TPATE1.Plasma.ncRNA.data))
Set1.TPATE1 = sample(x = TPATE1, size = length(TPATE1)/2, replace = FALSE)
Set2.TPATE1 = TPATE1[!TPATE1 %in% Set1.TPATE1]

TPATE1.Plasma.ncRNA.data.Set1 = TPATE1.Plasma.ncRNA.data[,Set1.TPATE1]
TPATE1.Plasma.ncRNA.data.Set2 = TPATE1.Plasma.ncRNA.data[,Set2.TPATE1]

Samples.TPATE1 = cor(TPATE1.Plasma.ncRNA.data.Set1, TPATE1.Plasma.ncRNA.data.Set2)
png("~/Documents/BCM/Lab/LouiseData/EDec_Analysis_Laurent_09/Stage_1/Correlation.chosenProbes.TPATE1.Plasma.Samples.png",1000,1000)
heatmap.2(as.matrix(Samples.TPATE1),
          trace = "none",
          col=my_palette,
          breaks = col_breaks,
          cexCol = 4,
          Colv = TRUE,
          margins = c(10,10))
dev.off()

EDec.TPATE1.ncRNA.Set1 = EDecStage1(methMixtureSamples = TPATE1.Plasma.ncRNA.data.Set1, cellTypeSpecificLoci = chosenProbes.allncRNA, nCts = 3, maxIts=2000, rssDiffStop=1e-10)
EDec.TPATE1.ncRNA.Set1.methy = EDec.TPATE1.ncRNA.Set1$methylation
EDec.TPATE1.ncRNA.Set1.methy.chosenProbes = EDec.TPATE1.ncRNA.Set1.methy[chosenProbes.allncRNA,]
EDec.TPATE1.ncRNA.Set1.prop = as.data.frame(t(EDec.TPATE1.ncRNA.Set1$proportions))

CorMatrix.TPATE1.ncRNA.Set1 = cor(EDec.TPATE1.ncRNA.Set1.methy.chosenProbes, Lasser.ncRNA.data.chosenProbes)
png("~/Documents/BCM/Lab/LouiseData/EDec_Analysis_Laurent_09/Stage_1/Profiles.TPATE1.Plasma.chosenProbes.ncRNA.3.Set1.png",1000,1000)
heatmap.2(as.matrix(t(CorMatrix.TPATE1.ncRNA.Set1)),
          trace = "none",
          col=my_palette,
          breaks = col_breaks,
          cexCol = 4,
          Colv = TRUE,
          margins = c(10,10))
dev.off()

png("~/Documents/BCM/Lab/LouiseData/EDec_Analysis_Laurent_09/Stage_1/Proportions.TPATE1.Plasma.chosenProbes.ncRNA.3.Set1.png",1000,1000)
heatmap.2(as.matrix(EDec.TPATE1.ncRNA.Set1.prop),
          trace="none",
          col=my_palette,
          breaks = col_breaks,
          labCol = FALSE,
          cexRow = 4,
          margins=c(10,10))
dev.off()

EDec.TPATE1.ncRNA.Set2 = EDecStage1(methMixtureSamples = TPATE1.Plasma.ncRNA.data.Set2, cellTypeSpecificLoci = chosenProbes.allncRNA, nCts = 3, maxIts=2000, rssDiffStop=1e-10)
EDec.TPATE1.ncRNA.Set2.methy = EDec.TPATE1.ncRNA.Set2$methylation
EDec.TPATE1.ncRNA.Set2.methy.chosenProbes = EDec.TPATE1.ncRNA.Set2.methy[chosenProbes.allncRNA,]
EDec.TPATE1.ncRNA.Set2.prop = as.data.frame(t(EDec.TPATE1.ncRNA.Set2$proportions))

CorMatrix.TPATE1.ncRNA.Set2 = cor(EDec.TPATE1.ncRNA.Set2.methy.chosenProbes, Lasser.ncRNA.data.chosenProbes)
png("~/Documents/BCM/Lab/LouiseData/EDec_Analysis_Laurent_09/Stage_1/Profiles.TPATE1.Plasma.chosenProbes.ncRNA.3.Set2.png",1000,1000)
heatmap.2(as.matrix(t(CorMatrix.TPATE1.ncRNA.Set2)),
          trace = "none",
          col=my_palette,
          breaks = col_breaks,
          cexCol = 4,
          Colv = TRUE,
          margins = c(10,10))
dev.off()

png("~/Documents/BCM/Lab/LouiseData/EDec_Analysis_Laurent_09/Stage_1/Proportions.TPATE1.Plasma.chosenProbes.ncRNA.3.Set2.png",1000,1000)
heatmap.2(as.matrix(EDec.TPATE1.ncRNA.Set2.prop),
          trace="none",
          col=my_palette,
          breaks = col_breaks,
          labCol = FALSE,
          cexRow = 4,
          margins=c(10,10))
dev.off()

EDec.TPATE1.ncRNA = EDecStage1(methMixtureSamples = TPATE1.Plasma.ncRNA.data, cellTypeSpecificLoci = chosenProbes.allncRNA, nCts = 3, maxIts=2000, rssDiffStop=1e-10)
EDec.TPATE1.ncRNA.methy = EDec.TPATE1.ncRNA$methylation
EDec.TPATE1.ncRNA.methy.chosenProbes = EDec.TPATE1.ncRNA.methy[chosenProbes.allncRNA,]
EDec.TPATE1.ncRNA.prop = as.data.frame(t(EDec.TPATE1.ncRNA$proportions))

CorMatrix.TPATE1.ncRNA = cor(EDec.TPATE1.ncRNA.methy.chosenProbes, Lasser.ncRNA.data.chosenProbes)
png("~/Documents/BCM/Lab/LouiseData/EDec_Analysis_Laurent_09/Stage_1/Profiles.TPATE1.Plasma.chosenProbes.ncRNA.3.All.png",1000,1000)
heatmap.2(as.matrix(t(CorMatrix.TPATE1.ncRNA)),
          trace = "none",
          col=my_palette,
          breaks = col_breaks,
          cexCol = 4,
          Colv = TRUE,
          margins = c(10,10))
dev.off()

png("~/Documents/BCM/Lab/LouiseData/EDec_Analysis_Laurent_09/Stage_1/Proportions.TPATE1.Plasma.chosenProbes.ncRNA.3.All.png",1000,1000)
heatmap.2(as.matrix(EDec.TPATE1.ncRNA.prop),
          trace="none",
          col=my_palette,
          breaks = col_breaks,
          labCol = FALSE,
          cexRow = 4,
          margins=c(10,10))
dev.off()

#################################
## EDec - KJENS1AP, Serum, Healthy
#################################
KJENS1AP.Serum = c(1:length(KJENS1AP.Serum.ncRNA.data))
Set1.KJENS1AP.Serum = sample(x = KJENS1AP.Serum, size = length(KJENS1AP.Serum)/2, replace = FALSE)
Set2.KJENS1AP.Serum = KJENS1AP.Serum[!KJENS1AP.Serum %in% Set1.KJENS1AP.Serum]

KJENS1AP.Serum.ncRNA.data.Set1 = KJENS1AP.Serum.ncRNA.data[,Set1.KJENS1AP.Serum]
KJENS1AP.Serum.ncRNA.data.Set2 = KJENS1AP.Serum.ncRNA.data[,Set2.KJENS1AP.Serum]

Samples.KJENS1AP.Serum = cor(KJENS1AP.Serum.ncRNA.data.Set1, KJENS1AP.Serum.ncRNA.data.Set2)
png("~/Documents/BCM/Lab/LouiseData/EDec_Analysis_Laurent_09/Stage_1/Correlation.chosenProbes.KJENS1AP.Serum.Samples.png",1000,1000)
heatmap.2(as.matrix(Samples.KJENS1AP.Serum),
          trace = "none",
          col=my_palette,
          breaks = col_breaks,
          cexCol = 4,
          Colv = TRUE,
          margins = c(10,10))
dev.off()

EDec.KJENS1AP.Serum.ncRNA.Set1 = EDecStage1(methMixtureSamples = KJENS1AP.Serum.ncRNA.data.Set1, cellTypeSpecificLoci = chosenProbes.allncRNA, nCts = 3, maxIts=2000, rssDiffStop=1e-10)
EDec.KJENS1AP.Serum.ncRNA.Set1.methy = EDec.KJENS1AP.Serum.ncRNA.Set1$methylation
EDec.KJENS1AP.Serum.ncRNA.Set1.methy.chosenProbes = EDec.KJENS1AP.Serum.ncRNA.Set1.methy[chosenProbes.allncRNA,]
EDec.KJENS1AP.Serum.ncRNA.Set1.prop = as.data.frame(t(EDec.KJENS1AP.Serum.ncRNA.Set1$proportions))

CorMatrix.KJENS1AP.Serum.ncRNA.Set1 = cor(EDec.KJENS1AP.Serum.ncRNA.Set1.methy.chosenProbes, Lasser.ncRNA.data.chosenProbes)
png("~/Documents/BCM/Lab/LouiseData/EDec_Analysis_Laurent_09/Stage_1/Profiles.KJENS1AP.Serum.chosenProbes.ncRNA.3.Set1.png",1000,1000)
heatmap.2(as.matrix(t(CorMatrix.KJENS1AP.Serum.ncRNA.Set1)),
          trace = "none",
          col=my_palette,
          breaks = col_breaks,
          cexCol = 4,
          Colv = TRUE,
          margins = c(10,10))
dev.off()

png("~/Documents/BCM/Lab/LouiseData/EDec_Analysis_Laurent_09/Stage_1/Proportions.KJENS1AP.Serum.chosenProbes.ncRNA.3.Set1.png",1000,1000)
heatmap.2(as.matrix(EDec.KJENS1AP.Serum.ncRNA.Set1.prop),
          trace="none",
          col=my_palette,
          breaks = col_breaks,
          labCol = FALSE,
          cexRow = 4,
          margins=c(10,10))
dev.off()

EDec.KJENS1AP.Serum.ncRNA.Set2 = EDecStage1(methMixtureSamples = KJENS1AP.Serum.ncRNA.data.Set2, cellTypeSpecificLoci = chosenProbes.allncRNA, nCts = 3, maxIts=2000, rssDiffStop=1e-10)
EDec.KJENS1AP.Serum.ncRNA.Set2.methy = EDec.KJENS1AP.Serum.ncRNA.Set2$methylation
EDec.KJENS1AP.Serum.ncRNA.Set2.methy.chosenProbes = EDec.KJENS1AP.Serum.ncRNA.Set2.methy[chosenProbes.allncRNA,]
EDec.KJENS1AP.Serum.ncRNA.Set2.prop = as.data.frame(t(EDec.KJENS1AP.Serum.ncRNA.Set2$proportions))

CorMatrix.KJENS1AP.Serum.ncRNA.Set2 = cor(EDec.KJENS1AP.Serum.ncRNA.Set2.methy.chosenProbes, Lasser.ncRNA.data.chosenProbes)
png("~/Documents/BCM/Lab/LouiseData/EDec_Analysis_Laurent_09/Stage_1/Profiles.KJENS1AP.Serum.chosenProbes.ncRNA.3.Set2.png",1000,1000)
heatmap.2(as.matrix(t(CorMatrix.KJENS1AP.Serum.ncRNA.Set2)),
          trace = "none",
          col=my_palette,
          breaks = col_breaks,
          cexCol = 4,
          Colv = TRUE,
          margins = c(10,10))
dev.off()

png("~/Documents/BCM/Lab/LouiseData/EDec_Analysis_Laurent_09/Stage_1/Proportions.KJENS1AP.Serum.chosenProbes.ncRNA.3.Set2.png",1000,1000)
heatmap.2(as.matrix(EDec.KJENS1AP.Serum.ncRNA.Set2.prop),
          trace="none",
          col=my_palette,
          breaks = col_breaks,
          labCol = FALSE,
          cexRow = 4,
          margins=c(10,10))
dev.off()

EDec.KJENS1AP.Serum.ncRNA = EDecStage1(methMixtureSamples = KJENS1AP.Serum.ncRNA.data, cellTypeSpecificLoci = chosenProbes.allncRNA, nCts = 3, maxIts=2000, rssDiffStop=1e-10)
EDec.KJENS1AP.Serum.ncRNA.methy = EDec.KJENS1AP.Serum.ncRNA$methylation
EDec.KJENS1AP.Serum.ncRNA.methy.chosenProbes = EDec.KJENS1AP.Serum.ncRNA.methy[chosenProbes.allncRNA,]
EDec.KJENS1AP.Serum.ncRNA.prop = as.data.frame(t(EDec.KJENS1AP.Serum.ncRNA$proportions))

CorMatrix.KJENS1AP.Serum.ncRNA = cor(EDec.KJENS1AP.Serum.ncRNA.methy.chosenProbes, Lasser.ncRNA.data.chosenProbes)
png("~/Documents/BCM/Lab/LouiseData/EDec_Analysis_Laurent_09/Stage_1/Profiles.KJENS1AP.Serum.chosenProbes.ncRNA.3.png",1000,1000)
heatmap.2(as.matrix(t(CorMatrix.KJENS1AP.Serum.ncRNA)),
          trace = "none",
          col=my_palette,
          breaks = col_breaks,
          cexCol = 4,
          Colv = TRUE,
          margins = c(10,10))
dev.off()

png("~/Documents/BCM/Lab/LouiseData/EDec_Analysis_Laurent_09/Stage_1/Proportions.KJENS1AP.Serum.chosenProbes.ncRNA.Kit.3.png",1000,1000)
heatmap.2(as.matrix(EDec.KJENS1AP.Serum.ncRNA.prop),
          trace="none",
          col=my_palette,
          breaks = col_breaks,
          labCol = FALSE,
          cexRow = 4,
          margins=c(10,10))
dev.off()

#################################
## EDec - KJENS1AP, CSF, Healthy
#################################
KJENS1AP.CSF = c(1:length(KJENS1AP.CSF.ncRNA.data))
Set1.KJENS1AP.CSF = sample(x = KJENS1AP.CSF, size = length(KJENS1AP.CSF)/2, replace = FALSE)
Set2.KJENS1AP.CSF = KJENS1AP.CSF[!KJENS1AP.CSF %in% Set1.KJENS1AP.CSF]

KJENS1AP.CSF.ncRNA.data.Set1 = KJENS1AP.CSF.ncRNA.data[,Set1.KJENS1AP.CSF]
KJENS1AP.CSF.ncRNA.data.Set2 = KJENS1AP.CSF.ncRNA.data[,Set2.KJENS1AP.CSF]

Samples.KJENS1AP.CSF = cor(KJENS1AP.CSF.ncRNA.data.Set1, KJENS1AP.CSF.ncRNA.data.Set2)
png("~/Documents/BCM/Lab/LouiseData/EDec_Analysis_Laurent_09/Stage_1/Correlation.chosenProbes.KJENS1AP.CSF.Samples.png",1000,1000)
heatmap.2(as.matrix(Samples.KJENS1AP.CSF),
          trace = "none",
          col=my_palette,
          breaks = col_breaks,
          cexCol = 4,
          Colv = TRUE,
          margins = c(10,10))
dev.off()

EDec.KJENS1AP.CSF.ncRNA.Set1 = EDecStage1(methMixtureSamples = KJENS1AP.CSF.ncRNA.data.Set1, cellTypeSpecificLoci = chosenProbes.allncRNA, nCts = 3, maxIts=2000, rssDiffStop=1e-10)
EDec.KJENS1AP.CSF.ncRNA.Set1.methy = EDec.KJENS1AP.CSF.ncRNA.Set1$methylation
EDec.KJENS1AP.CSF.ncRNA.Set1.methy.chosenProbes = EDec.KJENS1AP.CSF.ncRNA.Set1.methy[chosenProbes.allncRNA,]
EDec.KJENS1AP.CSF.ncRNA.Set1.prop = as.data.frame(t(EDec.KJENS1AP.CSF.ncRNA.Set1$proportions))

CorMatrix.KJENS1AP.CSF.ncRNA.Set1 = cor(EDec.KJENS1AP.CSF.ncRNA.Set1.methy.chosenProbes, Lasser.ncRNA.data.chosenProbes)
png("~/Documents/BCM/Lab/LouiseData/EDec_Analysis_Laurent_09/Stage_1/Profiles.KJENS1AP.CSF.chosenProbes.ncRNA.3.Set1.png",1000,1000)
heatmap.2(as.matrix(t(CorMatrix.KJENS1AP.CSF.ncRNA.Set1)),
          trace = "none",
          col=my_palette,
          breaks = col_breaks,
          cexCol = 4,
          Colv = TRUE,
          margins = c(10,10))
dev.off()

png("~/Documents/BCM/Lab/LouiseData/EDec_Analysis_Laurent_09/Stage_1/Proportions.KJENS1AP.CSF.chosenProbes.ncRNA.3.Set1.png",1000,1000)
heatmap.2(as.matrix(EDec.KJENS1AP.CSF.ncRNA.Set1.prop),
          trace="none",
          col=my_palette,
          breaks = col_breaks,
          labCol = FALSE,
          cexRow = 4,
          margins=c(10,10))
dev.off()

EDec.KJENS1AP.CSF.ncRNA.Set2 = EDecStage1(methMixtureSamples = KJENS1AP.CSF.ncRNA.data.Set2, cellTypeSpecificLoci = chosenProbes.allncRNA, nCts = 3, maxIts=2000, rssDiffStop=1e-10)
EDec.KJENS1AP.CSF.ncRNA.Set2.methy = EDec.KJENS1AP.CSF.ncRNA.Set2$methylation
EDec.KJENS1AP.CSF.ncRNA.Set2.methy.chosenProbes = EDec.KJENS1AP.CSF.ncRNA.Set2.methy[chosenProbes.allncRNA,]
EDec.KJENS1AP.CSF.ncRNA.Set2.prop = as.data.frame(t(EDec.KJENS1AP.CSF.ncRNA.Set2$proportions))

CorMatrix.KJENS1AP.CSF.ncRNA.Set2 = cor(EDec.KJENS1AP.CSF.ncRNA.Set2.methy.chosenProbes, Lasser.ncRNA.data.chosenProbes)
png("~/Documents/BCM/Lab/LouiseData/EDec_Analysis_Laurent_09/Stage_1/Profiles.KJENS1AP.CSF.chosenProbes.ncRNA.3.Set2.png",1000,1000)
heatmap.2(as.matrix(t(CorMatrix.KJENS1AP.CSF.ncRNA.Set2)),
          trace = "none",
          col=my_palette,
          breaks = col_breaks,
          cexCol = 4,
          Colv = TRUE,
          margins = c(10,10))
dev.off()

png("~/Documents/BCM/Lab/LouiseData/EDec_Analysis_Laurent_09/Stage_1/Proportions.KJENS1AP.CSF.chosenProbes.ncRNA.3.Set2.png",1000,1000)
heatmap.2(as.matrix(EDec.KJENS1AP.CSF.ncRNA.Set2.prop),
          trace="none",
          col=my_palette,
          breaks = col_breaks,
          labCol = FALSE,
          cexRow = 4,
          margins=c(10,10))
dev.off()

EDec.KJENS1AP.CSF.ncRNA = EDecStage1(methMixtureSamples = KJENS1AP.CSF.ncRNA.data, cellTypeSpecificLoci = chosenProbes.allncRNA, nCts = 3, maxIts=2000, rssDiffStop=1e-10)
EDec.KJENS1AP.CSF.ncRNA.methy = EDec.KJENS1AP.CSF.ncRNA$methylation
EDec.KJENS1AP.CSF.ncRNA.methy.chosenProbes = EDec.KJENS1AP.CSF.ncRNA.methy[chosenProbes.allncRNA,]
EDec.KJENS1AP.CSF.ncRNA.prop = as.data.frame(t(EDec.KJENS1AP.CSF.ncRNA$proportions))

CorMatrix.KJENS1AP.CSF.ncRNA = cor(EDec.KJENS1AP.CSF.ncRNA.methy.chosenProbes, Lasser.ncRNA.data.chosenProbes)
png("~/Documents/BCM/Lab/LouiseData/EDec_Analysis_Laurent_09/Stage_1/Profiles.KJENS1AP.CSF.chosenProbes.ncRNA.3.png",1000,1000)
heatmap.2(as.matrix(t(CorMatrix.KJENS1AP.CSF.ncRNA)),
          trace = "none",
          col=my_palette,
          breaks = col_breaks,
          cexCol = 4,
          Colv = TRUE,
          margins = c(10,10))
dev.off()

png("~/Documents/BCM/Lab/LouiseData/EDec_Analysis_Laurent_09/Stage_1/Proportions.KJENS1AP.CSF.chosenProbes.ncRNA.3.png",1000,1000)
heatmap.2(as.matrix(EDec.KJENS1AP.CSF.ncRNA.prop),
          trace="none",
          col=my_palette,
          breaks = col_breaks,
          labCol = FALSE,
          cexRow = 4,
          margins=c(10,10))
dev.off()

#################################
## EDec - KJENS1RID, Plasma, Healthy
#################################
KJENS1RID.Plasma = c(1:length(KJENS1RID.Plasma.ncRNA.data))
Set1.KJENS1RID.Plasma = sample(x = KJENS1RID.Plasma, size = length(KJENS1RID.Plasma)/2, replace = FALSE)
Set2.KJENS1RID.Plasma = KJENS1RID.Plasma[!KJENS1RID.Plasma %in% Set1.KJENS1RID.Plasma]

KJENS1RID.Plasma.ncRNA.data.Set1 = KJENS1RID.Plasma.ncRNA.data[,Set1.KJENS1RID.Plasma]
KJENS1RID.Plasma.ncRNA.data.Set2 = KJENS1RID.Plasma.ncRNA.data[,Set2.KJENS1RID.Plasma]

Samples.KJENS1RID.Plasma = cor(KJENS1RID.Plasma.ncRNA.data.Set1, KJENS1RID.Plasma.ncRNA.data.Set2)
png("~/Documents/BCM/Lab/LouiseData/EDec_Analysis_Laurent_09/Stage_1/Correlation.chosenProbes.KJENS1RID.Plasma.Samples.png",1000,1000)
heatmap.2(as.matrix(Samples.KJENS1RID.Plasma),
          trace = "none",
          col=my_palette,
          breaks = col_breaks,
          cexCol = 4,
          Colv = TRUE,
          margins = c(10,10))
dev.off()

EDec.KJENS1RID.Plasma.ncRNA.Set1 = EDecStage1(methMixtureSamples = KJENS1RID.Plasma.ncRNA.data.Set1, cellTypeSpecificLoci = chosenProbes.allncRNA, nCts = 3, maxIts=2000, rssDiffStop=1e-10)
EDec.KJENS1RID.Plasma.ncRNA.Set1.methy = EDec.KJENS1RID.Plasma.ncRNA.Set1$methylation
EDec.KJENS1RID.Plasma.ncRNA.Set1.methy.chosenProbes = EDec.KJENS1RID.Plasma.ncRNA.Set1.methy[chosenProbes.allncRNA,]
EDec.KJENS1RID.Plasma.ncRNA.Set1.prop = as.data.frame(t(EDec.KJENS1RID.Plasma.ncRNA.Set1$proportions))

CorMatrix.KJENS1RID.Plasma.ncRNA.Set1 = cor(EDec.KJENS1RID.Plasma.ncRNA.Set1.methy.chosenProbes, Lasser.ncRNA.data.chosenProbes)
png("~/Documents/BCM/Lab/LouiseData/EDec_Analysis_Laurent_09/Stage_1/Profiles.KJENS1RID.Plasma.chosenProbes.ncRNA.3.Set1.png",1000,1000)
heatmap.2(as.matrix(t(CorMatrix.KJENS1RID.Plasma.ncRNA.Set1)),
          trace = "none",
          col=my_palette,
          breaks = col_breaks,
          cexCol = 4,
          Colv = TRUE,
          margins = c(10,10))
dev.off()

png("~/Documents/BCM/Lab/LouiseData/EDec_Analysis_Laurent_09/Stage_1/Proportions.KJENS1RID.Plasma.chosenProbes.ncRNA.3.Set1.png",1000,1000)
heatmap.2(as.matrix(EDec.KJENS1RID.Plasma.ncRNA.Set1.prop),
          trace="none",
          col=my_palette,
          breaks = col_breaks,
          labCol = FALSE,
          cexRow = 4,
          margins=c(10,10))
dev.off()

EDec.KJENS1RID.Plasma.ncRNA.Set2 = EDecStage1(methMixtureSamples = KJENS1RID.Plasma.ncRNA.data.Set2, cellTypeSpecificLoci = chosenProbes.allncRNA, nCts = 3, maxIts=2000, rssDiffStop=1e-10)
EDec.KJENS1RID.Plasma.ncRNA.Set2.methy = EDec.KJENS1RID.Plasma.ncRNA.Set2$methylation
EDec.KJENS1RID.Plasma.ncRNA.Set2.methy.chosenProbes = EDec.KJENS1RID.Plasma.ncRNA.Set2.methy[chosenProbes.allncRNA,]
EDec.KJENS1RID.Plasma.ncRNA.Set2.prop = as.data.frame(t(EDec.KJENS1RID.Plasma.ncRNA.Set2$proportions))

CorMatrix.KJENS1RID.Plasma.ncRNA.Set2 = cor(EDec.KJENS1RID.Plasma.ncRNA.Set2.methy.chosenProbes, Lasser.ncRNA.data.chosenProbes)
png("~/Documents/BCM/Lab/LouiseData/EDec_Analysis_Laurent_09/Stage_1/Profiles.KJENS1RID.Plasma.chosenProbes.ncRNA.3.Set2.png",1000,1000)
heatmap.2(as.matrix(t(CorMatrix.KJENS1RID.Plasma.ncRNA.Set2)),
          trace = "none",
          col=my_palette,
          breaks = col_breaks,
          cexCol = 4,
          Colv = TRUE,
          margins = c(10,10))
dev.off()

png("~/Documents/BCM/Lab/LouiseData/EDec_Analysis_Laurent_09/Stage_1/Proportions.KJENS1RID.Plasma.chosenProbes.ncRNA.3.Set2.png",1000,1000)
heatmap.2(as.matrix(EDec.KJENS1RID.Plasma.ncRNA.Set2.prop),
          trace="none",
          col=my_palette,
          breaks = col_breaks,
          labCol = FALSE,
          cexRow = 4,
          margins=c(10,10))
dev.off()

EDec.KJENS1RID.Plasma.ncRNA = EDecStage1(methMixtureSamples = KJENS1RID.Plasma.ncRNA.data, cellTypeSpecificLoci = chosenProbes.allncRNA, nCts = 3, maxIts=2000, rssDiffStop=1e-10)
EDec.KJENS1RID.Plasma.ncRNA.methy = EDec.KJENS1RID.Plasma.ncRNA$methylation
EDec.KJENS1RID.Plasma.ncRNA.methy.chosenProbes = EDec.KJENS1RID.Plasma.ncRNA.methy[chosenProbes.allncRNA,]
EDec.KJENS1RID.Plasma.ncRNA.prop = as.data.frame(t(EDec.KJENS1RID.Plasma.ncRNA$proportions))

CorMatrix.KJENS1RID.Plasma.ncRNA = cor(EDec.KJENS1RID.Plasma.ncRNA.methy.chosenProbes, Lasser.ncRNA.data.chosenProbes)
png("~/Documents/BCM/Lab/LouiseData/EDec_Analysis_Laurent_09/Stage_1/Profiles.KJENS1RID.Plasma.chosenProbes.ncRNA.3.png",1000,1000)
heatmap.2(as.matrix(t(CorMatrix.KJENS1RID.Plasma.ncRNA)),
          trace = "none",
          col=my_palette,
          breaks = col_breaks,
          cexCol = 4,
          Colv = TRUE,
          margins = c(10,10))
dev.off()

png("~/Documents/BCM/Lab/LouiseData/EDec_Analysis_Laurent_09/Stage_1/Proportions.KJENS1RID.Plasma.chosenProbes.ncRNA.3.png",1000,1000)
heatmap.2(as.matrix(EDec.KJENS1RID.Plasma.ncRNA.prop),
          trace="none",
          col=my_palette,
          breaks = col_breaks,
          labCol = FALSE,
          cexRow = 4,
          margins=c(10,10))
dev.off()

#################################
## EDec - KJENS1RID, Saliva, Healthy
#################################
KJENS1RID.Saliva = c(1:length(KJENS1RID.Saliva.ncRNA.data))
Set1.KJENS1RID.Saliva = sample(x = KJENS1RID.Saliva, size = length(KJENS1RID.Saliva)/2, replace = FALSE)
Set2.KJENS1RID.Saliva = KJENS1RID.Saliva[!KJENS1RID.Saliva %in% Set1.KJENS1RID.Saliva]

KJENS1RID.Saliva.ncRNA.data.Set1 = KJENS1RID.Saliva.ncRNA.data[,Set1.KJENS1RID.Saliva]
KJENS1RID.Saliva.ncRNA.data.Set2 = KJENS1RID.Saliva.ncRNA.data[,Set2.KJENS1RID.Saliva]

Samples.KJENS1RID.Saliva = cor(KJENS1RID.Saliva.ncRNA.data.Set1, KJENS1RID.Saliva.ncRNA.data.Set2)
png("~/Documents/BCM/Lab/LouiseData/EDec_Analysis_Laurent_09/Stage_1/Correlation.chosenProbes.KJENS1RID.Saliva.Samples.png",1000,1000)
heatmap.2(as.matrix(Samples.KJENS1RID.Saliva),
          trace = "none",
          col=my_palette,
          breaks = col_breaks,
          cexCol = 4,
          Colv = TRUE,
          margins = c(10,10))
dev.off()

EDec.KJENS1RID.Saliva.ncRNA.Set1 = EDecStage1(methMixtureSamples = KJENS1RID.Saliva.ncRNA.data.Set1, cellTypeSpecificLoci = chosenProbes.allncRNA, nCts = 3, maxIts=2000, rssDiffStop=1e-10)
EDec.KJENS1RID.Saliva.ncRNA.Set1.methy = EDec.KJENS1RID.Saliva.ncRNA.Set1$methylation
EDec.KJENS1RID.Saliva.ncRNA.Set1.methy.chosenProbes = EDec.KJENS1RID.Saliva.ncRNA.Set1.methy[chosenProbes.allncRNA,]
EDec.KJENS1RID.Saliva.ncRNA.Set1.prop = as.data.frame(t(EDec.KJENS1RID.Saliva.ncRNA.Set1$proportions))

CorMatrix.KJENS1RID.Saliva.ncRNA.Set1 = cor(EDec.KJENS1RID.Saliva.ncRNA.Set1.methy.chosenProbes, Lasser.ncRNA.data.chosenProbes)
png("~/Documents/BCM/Lab/LouiseData/EDec_Analysis_Laurent_09/Stage_1/Profiles.KJENS1RID.Saliva.chosenProbes.ncRNA.3.Set1.png",1000,1000)
heatmap.2(as.matrix(t(CorMatrix.KJENS1RID.Saliva.ncRNA.Set1)),
          trace = "none",
          col=my_palette,
          breaks = col_breaks,
          cexCol = 4,
          Colv = TRUE,
          margins = c(10,10))
dev.off()

png("~/Documents/BCM/Lab/LouiseData/EDec_Analysis_Laurent_09/Stage_1/Proportions.KJENS1RID.Saliva.chosenProbes.ncRNA.3.Set1.png",1000,1000)
heatmap.2(as.matrix(EDec.KJENS1RID.Saliva.ncRNA.Set1.prop),
          trace="none",
          col=my_palette,
          breaks = col_breaks,
          labCol = FALSE,
          cexRow = 4,
          margins=c(10,10))
dev.off()

EDec.KJENS1RID.Saliva.ncRNA.Set2 = EDecStage1(methMixtureSamples = KJENS1RID.Saliva.ncRNA.data.Set2, cellTypeSpecificLoci = chosenProbes.allncRNA, nCts = 3, maxIts=2000, rssDiffStop=1e-10)
EDec.KJENS1RID.Saliva.ncRNA.Set2.methy = EDec.KJENS1RID.Saliva.ncRNA.Set2$methylation
EDec.KJENS1RID.Saliva.ncRNA.Set2.methy.chosenProbes = EDec.KJENS1RID.Saliva.ncRNA.Set2.methy[chosenProbes.allncRNA,]
EDec.KJENS1RID.Saliva.ncRNA.Set2.prop = as.data.frame(t(EDec.KJENS1RID.Saliva.ncRNA.Set2$proportions))

CorMatrix.KJENS1RID.Saliva.ncRNA.Set2 = cor(EDec.KJENS1RID.Saliva.ncRNA.Set2.methy.chosenProbes, Lasser.ncRNA.data.chosenProbes)
png("~/Documents/BCM/Lab/LouiseData/EDec_Analysis_Laurent_09/Stage_1/Profiles.KJENS1RID.Saliva.chosenProbes.ncRNA.3.Set2.png",1000,1000)
heatmap.2(as.matrix(t(CorMatrix.KJENS1RID.Saliva.ncRNA.Set2)),
          trace = "none",
          col=my_palette,
          breaks = col_breaks,
          cexCol = 4,
          Colv = TRUE,
          margins = c(10,10))
dev.off()

png("~/Documents/BCM/Lab/LouiseData/EDec_Analysis_Laurent_09/Stage_1/Proportions.KJENS1RID.Saliva.chosenProbes.ncRNA.3.Set2.png",1000,1000)
heatmap.2(as.matrix(EDec.KJENS1RID.Saliva.ncRNA.Set2.prop),
          trace="none",
          col=my_palette,
          breaks = col_breaks,
          labCol = FALSE,
          cexRow = 4,
          margins=c(10,10))
dev.off()

EDec.KJENS1RID.Saliva.ncRNA = EDecStage1(methMixtureSamples = KJENS1RID.Saliva.ncRNA.data, cellTypeSpecificLoci = chosenProbes.allncRNA, nCts = 3, maxIts=2000, rssDiffStop=1e-10)
EDec.KJENS1RID.Saliva.ncRNA.methy = EDec.KJENS1RID.Saliva.ncRNA$methylation
EDec.KJENS1RID.Saliva.ncRNA.methy.chosenProbes = EDec.KJENS1RID.Saliva.ncRNA.methy[chosenProbes.allncRNA,]
EDec.KJENS1RID.Saliva.ncRNA.prop = as.data.frame(t(EDec.KJENS1RID.Saliva.ncRNA$proportions))

CorMatrix.KJENS1RID.Saliva.ncRNA = cor(EDec.KJENS1RID.Saliva.ncRNA.methy.chosenProbes, Lasser.ncRNA.data.chosenProbes)
png("~/Documents/BCM/Lab/LouiseData/EDec_Analysis_Laurent_09/Stage_1/Profiles.KJENS1RID.Saliva.chosenProbes.ncRNA.3.png",1000,1000)
heatmap.2(as.matrix(t(CorMatrix.KJENS1RID.Saliva.ncRNA)),
          trace = "none",
          col=my_palette,
          breaks = col_breaks,
          cexCol = 4,
          Colv = TRUE,
          margins = c(10,10))
dev.off()

png("~/Documents/BCM/Lab/LouiseData/EDec_Analysis_Laurent_09/Stage_1/Proportions.KJENS1RID.Saliva.chosenProbes.ncRNA.3.png",1000,1000)
heatmap.2(as.matrix(EDec.KJENS1RID.Saliva.ncRNA.prop),
          trace="none",
          col=my_palette,
          breaks = col_breaks,
          labCol = FALSE,
          cexRow = 4,
          margins=c(10,10))
dev.off()

#################################
## EDec - KJENS1RID, Urine, Healthy
#################################
KJENS1RID.Urine = c(1:length(KJENS1RID.Urine.ncRNA.data))
Set1.KJENS1RID.Urine = sample(x = KJENS1RID.Urine, size = length(KJENS1RID.Urine)/2, replace = FALSE)
Set2.KJENS1RID.Urine = KJENS1RID.Urine[!KJENS1RID.Urine %in% Set1.KJENS1RID.Urine]

KJENS1RID.Urine.ncRNA.data.Set1 = KJENS1RID.Urine.ncRNA.data[,Set1.KJENS1RID.Urine]
KJENS1RID.Urine.ncRNA.data.Set2 = KJENS1RID.Urine.ncRNA.data[,Set2.KJENS1RID.Urine]

Samples.KJENS1RID.Urine = cor(KJENS1RID.Urine.ncRNA.data.Set1, KJENS1RID.Urine.ncRNA.data.Set2)
png("~/Documents/BCM/Lab/LouiseData/EDec_Analysis_Laurent_09/Stage_1/Correlation.chosenProbes.KJENS1RID.Urine.Samples.png",1000,1000)
heatmap.2(as.matrix(Samples.KJENS1RID.Urine),
          trace = "none",
          col=my_palette,
          breaks = col_breaks,
          cexCol = 4,
          Colv = TRUE,
          margins = c(10,10))
dev.off()

EDec.KJENS1RID.Urine.ncRNA.Set1 = EDecStage1(methMixtureSamples = KJENS1RID.Urine.ncRNA.data.Set1, cellTypeSpecificLoci = chosenProbes.allncRNA, nCts = 3, maxIts=2000, rssDiffStop=1e-10)
EDec.KJENS1RID.Urine.ncRNA.Set1.methy = EDec.KJENS1RID.Urine.ncRNA.Set1$methylation
EDec.KJENS1RID.Urine.ncRNA.Set1.methy.chosenProbes = EDec.KJENS1RID.Urine.ncRNA.Set1.methy[chosenProbes.allncRNA,]
EDec.KJENS1RID.Urine.ncRNA.Set1.prop = as.data.frame(t(EDec.KJENS1RID.Urine.ncRNA.Set1$proportions))

CorMatrix.KJENS1RID.Urine.ncRNA.Set1 = cor(EDec.KJENS1RID.Urine.ncRNA.Set1.methy.chosenProbes, Lasser.ncRNA.data.chosenProbes)
png("~/Documents/BCM/Lab/LouiseData/EDec_Analysis_Laurent_09/Stage_1/Profiles.KJENS1RID.Urine.chosenProbes.ncRNA.3.Set1.png",1000,1000)
heatmap.2(as.matrix(t(CorMatrix.KJENS1RID.Urine.ncRNA.Set1)),
          trace = "none",
          col=my_palette,
          breaks = col_breaks,
          cexCol = 4,
          Colv = TRUE,
          margins = c(10,10))
dev.off()

png("~/Documents/BCM/Lab/LouiseData/EDec_Analysis_Laurent_09/Stage_1/Proportions.KJENS1RID.Urine.chosenProbes.ncRNA.3.Set1.png",1000,1000)
heatmap.2(as.matrix(EDec.KJENS1RID.Urine.ncRNA.Set1.prop),
          trace="none",
          col=my_palette,
          breaks = col_breaks,
          labCol = FALSE,
          cexRow = 4,
          margins=c(10,10))
dev.off()

EDec.KJENS1RID.Urine.ncRNA.Set2 = EDecStage1(methMixtureSamples = KJENS1RID.Urine.ncRNA.data.Set2, cellTypeSpecificLoci = chosenProbes.allncRNA, nCts = 3, maxIts=2000, rssDiffStop=1e-10)
EDec.KJENS1RID.Urine.ncRNA.Set2.methy = EDec.KJENS1RID.Urine.ncRNA.Set2$methylation
EDec.KJENS1RID.Urine.ncRNA.Set2.methy.chosenProbes = EDec.KJENS1RID.Urine.ncRNA.Set2.methy[chosenProbes.allncRNA,]
EDec.KJENS1RID.Urine.ncRNA.Set2.prop = as.data.frame(t(EDec.KJENS1RID.Urine.ncRNA.Set2$proportions))

CorMatrix.KJENS1RID.Urine.ncRNA.Set2 = cor(EDec.KJENS1RID.Urine.ncRNA.Set2.methy.chosenProbes, Lasser.ncRNA.data.chosenProbes)
png("~/Documents/BCM/Lab/LouiseData/EDec_Analysis_Laurent_09/Stage_1/Profiles.KJENS1RID.Urine.chosenProbes.ncRNA.3.Set2.png",1000,1000)
heatmap.2(as.matrix(t(CorMatrix.KJENS1RID.Urine.ncRNA.Set2)),
          trace = "none",
          col=my_palette,
          breaks = col_breaks,
          cexCol = 4,
          Colv = TRUE,
          margins = c(10,10))
dev.off()

png("~/Documents/BCM/Lab/LouiseData/EDec_Analysis_Laurent_09/Stage_1/Proportions.KJENS1RID.Urine.chosenProbes.ncRNA.3.Set2.png",1000,1000)
heatmap.2(as.matrix(EDec.KJENS1RID.Urine.ncRNA.Set2.prop),
          trace="none",
          col=my_palette,
          breaks = col_breaks,
          labCol = FALSE,
          cexRow = 4,
          margins=c(10,10))
dev.off()

EDec.KJENS1RID.Urine.ncRNA = EDecStage1(methMixtureSamples = KJENS1RID.Urine.ncRNA.data, cellTypeSpecificLoci = chosenProbes.allncRNA, nCts = 3, maxIts=2000, rssDiffStop=1e-10)
EDec.KJENS1RID.Urine.ncRNA.methy = EDec.KJENS1RID.Urine.ncRNA$methylation
EDec.KJENS1RID.Urine.ncRNA.methy.chosenProbes = EDec.KJENS1RID.Urine.ncRNA.methy[chosenProbes.allncRNA,]
EDec.KJENS1RID.Urine.ncRNA.prop = as.data.frame(t(EDec.KJENS1RID.Urine.ncRNA$proportions))

CorMatrix.KJENS1RID.Urine.ncRNA = cor(EDec.KJENS1RID.Urine.ncRNA.methy.chosenProbes, Lasser.ncRNA.data.chosenProbes)
png("~/Documents/BCM/Lab/LouiseData/EDec_Analysis_Laurent_09/Stage_1/Profiles.KJENS1RID.Urine.chosenProbes.ncRNA.3.png",1000,1000)
heatmap.2(as.matrix(t(CorMatrix.KJENS1RID.Urine.ncRNA)),
          trace = "none",
          col=my_palette,
          breaks = col_breaks,
          cexCol = 4,
          Colv = TRUE,
          margins = c(10,10))
dev.off()

png("~/Documents/BCM/Lab/LouiseData/EDec_Analysis_Laurent_09/Stage_1/Proportions.KJENS1RID.Urine.chosenProbes.ncRNA.3.png",1000,1000)
heatmap.2(as.matrix(EDec.KJENS1RID.Urine.ncRNA.prop),
          trace="none",
          col=my_palette,
          breaks = col_breaks,
          labCol = FALSE,
          cexRow = 4,
          margins=c(10,10))
dev.off()

#####################################
## Name all predicted samples by hand
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
Sets.KJENS1AP.Serum = round(cor(EDec.KJENS1AP.Serum.ncRNA.Set1.methy.chosenProbes, EDec.KJENS1AP.Serum.ncRNA.Set2.methy.chosenProbes, use="complete"), digits=2)
png("~/Documents/BCM/Lab/LouiseData/EDec_Analysis_Laurent_09/Stage_1/1.Correlation.chosenProbes.KJENS1AP.Serum.Sets.png",1000,1000)
qplot(x=Var1, y=Var2, data=melt(Sets.KJENS1AP.Serum), fill=value, geom="tile") + theme(axis.text = element_text(size = 24)) + geom_text(aes(Var2, Var1, label = value), color = "black", size = 8) + scale_fill_gradient2(low="white", high="darkgreen")
dev.off()

Sets.KJENS1RID.Saliva = round(cor(EDec.KJENS1RID.Saliva.ncRNA.Set1.methy.chosenProbes, EDec.KJENS1RID.Saliva.ncRNA.Set2.methy.chosenProbes, use="complete"), digits=2)
png("~/Documents/BCM/Lab/LouiseData/EDec_Analysis_Laurent_09/Stage_1/1.Correlation.chosenProbes.KJENS1RID.Saliva.Sets.png",1000,1000)
qplot(x=Var1, y=Var2, data=melt(Sets.KJENS1RID.Saliva), fill=value, geom="tile") + theme(axis.text = element_text(size = 24)) + geom_text(aes(Var2, Var1, label = value), color = "black", size = 8) + scale_fill_gradient2(low="white", high="darkgreen")
dev.off()

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

Sets.KJENS1AP.Serum.KJENS1RID.Urine = round(cor(EDec.KJENS1AP.Serum.ncRNA.methy, EDec.KJENS1RID.Urine.ncRNA.methy, use="complete"), digits=2)
png("~/Documents/BCM/Lab/LouiseData/EDec_Analysis_Laurent_09/Stage_1/2.a.Correlation.KJENS1AP.Serum.KJENS1RID.Urine.png",1000,1000)
qplot(x=Var1, y=Var2, data=melt(Sets.KJENS1AP.Serum.KJENS1RID.Urine), fill=value, geom="tile") + theme(axis.text = element_text(size = 24)) + geom_text(aes(Var2, Var1, label = value), color = "black", size = 8) + scale_fill_gradient2(low="white", high="darkgreen")
dev.off()

Sets.KJENS1RID.Saliva.KJENS1RID.Urine = round(cor(EDec.KJENS1RID.Saliva.ncRNA.methy, EDec.KJENS1RID.Urine.ncRNA.methy, use="complete"), digits=2)
png("~/Documents/BCM/Lab/LouiseData/EDec_Analysis_Laurent_09/Stage_1/2.a.Correlation.KJENS1RID.Saliva.KJENS1RID.Urine.png",1000,1000)
qplot(x=Var1, y=Var2, data=melt(Sets.KJENS1RID.Saliva.KJENS1RID.Urine), fill=value, geom="tile") + theme(axis.text = element_text(size = 24)) + geom_text(aes(Var2, Var1, label = value), color = "black", size = 8) + scale_fill_gradient2(low="white", high="darkgreen")
dev.off()

## chosen Probes
Sets.KJENS1AP.Serum.KJENS1RID.Saliva = round(cor(EDec.KJENS1AP.Serum.ncRNA.methy.chosenProbes, EDec.KJENS1RID.Saliva.ncRNA.methy.chosenProbes, use="complete"), digits=2)
png("~/Documents/BCM/Lab/LouiseData/EDec_Analysis_Laurent_09/Stage_1/2.b.Correlation.chosenProbes.KJENS1AP.Serum.KJENS1RID.Saliva.png",1000,1000)
qplot(x=Var1, y=Var2, data=melt(Sets.KJENS1AP.Serum.KJENS1RID.Saliva), fill=value, geom="tile") + theme(axis.text = element_text(size = 24)) + geom_text(aes(Var2, Var1, label = value), color = "black", size = 8) + scale_fill_gradient2(low="white", high="darkgreen")
dev.off()

Sets.KJENS1AP.Serum.KJENS1RID.Urine = round(cor(EDec.KJENS1AP.Serum.ncRNA.methy.chosenProbes, EDec.KJENS1RID.Urine.ncRNA.methy.chosenProbes, use="complete"), digits=2)
png("~/Documents/BCM/Lab/LouiseData/EDec_Analysis_Laurent_09/Stage_1/2.b.Correlation.chosenProbes.KJENS1AP.Serum.KJENS1RID.Urine.png",1000,1000)
qplot(x=Var1, y=Var2, data=melt(Sets.KJENS1AP.Serum.KJENS1RID.Urine), fill=value, geom="tile") + theme(axis.text = element_text(size = 24)) + geom_text(aes(Var2, Var1, label = value), color = "black", size = 8) + scale_fill_gradient2(low="white", high="darkgreen")
dev.off()

Sets.KJENS1RID.Saliva.KJENS1RID.Urine = round(cor(EDec.KJENS1RID.Saliva.ncRNA.methy.chosenProbes, EDec.KJENS1RID.Urine.ncRNA.methy.chosenProbes, use="complete"), digits=2)
png("~/Documents/BCM/Lab/LouiseData/EDec_Analysis_Laurent_09/Stage_1/2.b.Correlation.chosenProbes.KJENS1RID.Saliva.KJENS1RID.Urine.png",1000,1000)
qplot(x=Var1, y=Var2, data=melt(Sets.KJENS1RID.Saliva.KJENS1RID.Urine), fill=value, geom="tile") + theme(axis.text = element_text(size = 24)) + geom_text(aes(Var2, Var1, label = value), color = "black", size = 8) + scale_fill_gradient2(low="white", high="darkgreen")
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
