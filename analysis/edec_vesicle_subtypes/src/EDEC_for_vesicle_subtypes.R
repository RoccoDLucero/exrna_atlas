source("./exrna_atlas_get_data_functions.R")




#This list is generated manually by inspection of ExceRpt Biotype counts outputs 
gencode_types <- list("Y_RNA", "snoRNA", "snRNA", "lincRNA", "miRNA",
                      "protein_coding","processed_transcript" )

#Identify exRNA detected in a given propotion of samples at a given rpm value
# so that we can reduce the size of the data
filterByPctDtect <- function(x, min_detect_proportion, min_rpm){
    detect_prop <- length(x[x >= min_rpm]) / length(x)
    
    if(detect_prop >= min_detect_proportion){return(TRUE)
    }else{ return(FALSE)}
    
}
#We want our chosen miRNA to show variability across conditions so we may be able to 
# remove miRNA that do not vary enough
filterByVariance <- function(x, min_variance){
    v <- var(x)
    if(v < min_variance){return(TRUE)}else{return(FALSE)}
    
}

filterBySD <- function(x, min_sd){
    s <- sd(x)
    if(s < min_sd){return(TRUE)}else{return(FALSE)}
    
}

filterByIQR <- function(x, min_iqr){
    i <- IQR(x)
    if(i < min_iqr){return(TRUE)}else{return(FALSE)}
    
}


removeRareRNA <- function(df, min_detect_proportion = .05, min_rpm = 1){
    r <- apply(df, 2, filterByPctDtect, min_detect_proportion, min_rpm)
    
    return(df[,r])
    
}

#Functions for various data transformations that will be neeed to ensure compatibility
# With EDec
txfmQN  <- function(df){preprocessCore::normalize.quantiles(df)} #columns will share a distribution
txfmNPN <- function(df){huge::huge.npn(x = df)}
txfmLin <- function(x){(x-min(x))/(max(x)-min(x))}
txfmLogistic1 <- function(x, a = 1/max(x)){ 1 / (1 + exp(1)^(-(a * x )) ) }
txfmLogistic2 <- function(x, a = 1/max(x)){ 1.5*( 1 - exp(1)^(-(a * x)) ) }
txfmLogistic3 <- function(x, a = 1/max(x)){ 1 / (1 + exp(1)^(-(a * (x-mean(x)) ) )) } 
txfmLogit <- function(x){log( x / (1 - x) )}

#######################################################
##Functions for selecting informative probe molecules##
perform_t_tests_all_vars <- function(dataGroup1,dataGroup2){
    nGroup1 <- nrow(dataGroup1)
    nGroup2 <- nrow(dataGroup2)
    dataAll <- rbind(dataGroup1,dataGroup2)
    tTestWithErrorHandling = function(x){
        testResult <- try(t.test(x[1:nGroup1],x[(nGroup1+1):(nGroup1+nGroup2)]),silent=TRUE);
        if(is.character(testResult)){
            warning(testResult)
            c(NA,NA,NA)
        }else{
            c(testResult$p.value,testResult$estimate)
        }
    }
    results <- matrix(unlist(apply(dataAll,2,tTestWithErrorHandling)),ncol=3,byrow=TRUE)
    
    results <- as.data.frame(results)
    colnames(results) = c("P.value","Mean.group.1","Mean.group.2")
    rownames(results) = colnames(dataGroup1)
    return(results)
}

perform_t_tests_all_classes_one_vs_rest = function(dataMatrix,classVector){
    if(nrow(dataMatrix)!=length(classVector)){
        stop("Number of samples in data matrix must be equal to the length of the class vector")
    }
    possibleClasses <- unique(classVector)
    nClasses <- length(possibleClasses)
    
    allPValues <- matrix(NA,nrow = ncol(dataMatrix),ncol=nClasses)
    allDiffMeans <- matrix(NA,nrow=ncol(dataMatrix),ncol=nClasses)
    
    colnames(allPValues) <- sapply(possibleClasses, paste, "vs.others", sep = "." )
    rownames(allPValues) <- colnames(dataMatrix)
    colnames(allDiffMeans) <- sapply(possibleClasses, paste, "vs.others", sep = "." )
    rownames(allDiffMeans) <- colnames(dataMatrix)
    
    for(i in 1:nClasses){
        class <- possibleClasses[i]
        resultTest <- perform_t_tests_all_vars(dataMatrix[classVector==class,],
                                              dataMatrix[classVector!=class,])
        allPValues[,i] <- resultTest[,1]
        allDiffMeans[,i] <- resultTest[,2] - resultTest[,3]
    }
    
    result = lapply(X = list(allPValues, allDiffMeans), as.data.frame)
    names(result) <- c("P.Values","Difference.Between.Means")
    return(result)
}

perform_t_tests_all_classes_each_pair = function(dataMatrix, classVector){
    if(nrow(dataMatrix) != length(classVector)){
        stop("Number of samples in data matrix must be equal to the length of the class vector")
    }
    possibleClasses <- unique(classVector)
    nClasses <- length(possibleClasses)
    
    allPValues <- NULL
    allDiffMeans <- NULL
    names <- NULL 
    for(i in 1:(nClasses-1)){
        for(j in (i+1):nClasses){
            class1 <- possibleClasses[i]
            class2 <- possibleClasses[j]
            names <- c(names,paste(class1,class2,sep="."))
            result <- perform_t_tests_all_vars(dataMatrix[classVector==class1,],
                                              dataMatrix[classVector==class2,])
            allPValues <- cbind(allPValues,result[,1])
            allDiffMeans <- cbind(allDiffMeans, abs((result[,2] - result[,3])))
        }
    }
    colnames(allPValues) <- names
    rownames(allPValues) <- colnames(dataMatrix)
    colnames(allDiffMeans) <- names
    rownames(allDiffMeans) <- colnames(dataMatrix)
    result <- lapply(X = list(allPValues, allDiffMeans), as.data.frame)
    names(result) <- c("P.Values", "Difference.Between.Means")
    return(result)
}



################################################################################
################################################################################
#Load the ExceRpt Post Processed results for Vesicle Subtype data from Lasser
#and Kit/prep method data set from Laurent RNA_ISO_LARGE (From Meenu) 

##Generic suffix for Excerpt Post processed result RDATA##
pth_tail_rpm        <- "_exceRpt_smallRNAQuants_ReadsPerMillion.RData"
pth_tail_raw_counts <- "_exceRpt_smallRNAQuants_ReadCounts.Rdata"

##LASSER HD LD##
pth_root <-file.path(".", "input", "Lasser")
pth_root <-file.path(pth_root, list.files(pth_root,recursive = F)[1])

inp_1 <-file.path(pth_root, paste("Lasser3_2017-5-7", pth_tail_rpm, sep = ""))
inp_2 <-file.path(pth_root, paste("Lasser1_Large_2017-5-9", pth_tail_rpm, sep = ""))

lasser_rpm_1 <- my.load.excerpt.rdata(local_rdat_file = inp_1)
lasser_rpm_2 <- my.load.excerpt.rdata(local_rdat_file = inp_2)

#Taken from the Lasser Paper:
#{Lab, CellLine, Biofluid, KIT, IsolationRep, SequencingRep, Vol_Biofluid}
lasser_meta <- (c("Lotvall", "HMC1_MastCell", "CellSup", "MiRCury", "No_Data", "1", "No_Data"))


lasser_meta_df <- as.data.frame(matrix(data = rep(lasser_meta,4), nrow = 4, byrow = T),
                                stringsAsFactors = F)
lasser_meta_df[,5] <- as.character(c(1,2,1,2))
rownames(lasser_meta_df) <- rownames(lasser_all$miRNA.rpm)

##LAURENT RNA_ISO_LARGE##
pth_root <-file.path(".", "input", "Laurent")

inp_3 <-file.path(pth_root, paste("Laurent_2017-5-4-11_RNAIsoLarge", pth_tail_rpm, sep = ""))

laurent_rpm_1 <- my.load.excerpt.rdata(local_rdat_file = inp_3)

laurent_meta <- read.delim("./input/Laurent/20170508_RNAisolation_RMsummary_annot_Meenu.txt",
                           header = T, sep = "\t", row.names = 1, stringsAsFactors = F)
laurent_meta <- laurent_meta[,1:7]
laurent_meta_sup <- laurent_meta[grep("Sup",laurent_meta$Biofluid),] 

setdiff(rownames(laurent_meta_sup), rownames(laurent_rpm_1$miRNA.rpm))
fix_rn <- function(df){rownames(df) <- gsub("ExoQuicN", "ExoQuick", rownames(df)); df}
laurent_rpm_1 <- lapply(laurent_rpm_1, fix_rn)

##MERGE Lasser and Laurent Data Sets
lasser_all <- mapply(my.combine.frames, lasser_rpm_1, lasser_rpm_2)
lasser_laurent_all <- mapply(my.combine.frames, lasser_all, laurent_rpm_1)

##MERGE Lasser and Laurent MetaData
colnames(lasser_meta_df) <- colnames(laurent_meta_sup)
lasser_laurent_meta <- rbind(lasser_meta_df, laurent_meta_sup)
lasser_laurent_meta <- as.data.frame(unclass(lasser_laurent_meta),
                                     stringsAsFactors = T, row.names = rownames(lasser_laurent_meta) )
##Remove unused ExceRpt output
lasser_laurent_all <- lasser_laurent_all[-grep("exogenous|endogenous",names(lasser_laurent_all))]

#split Gencode mapped output into multiple data frames by RNA Type
lasser_laurent_all <- c(lasser_laurent_all[-grep("gencode",names(lasser_laurent_all))],
                        sep_df_by_colname_pttn(lasser_laurent_all$gencode.rpm,
                                               gencode_types))

rm(inp_1, inp_2, inp_3, lasser_rpm_1, lasser_rpm_2)
rm(pth_root, pth_tail_raw_counts, pth_tail_rpm)

#####################################################################################
#^^^^ THE ExceRpt ouput for the Lasser and Laurent Studies have been combined
# and made more granular. They are now ready for transformation to a form that
# should be compatible with EDec.
################################################################################################

#This is to quickly confirm that Excerpt maps similar proportions of reads from HD LD
# as the "tuxedo" mapper used in the Lasser paper.
for(n in 1:length(lasser_laurent_all)){
    rna_class <- lasser_laurent_all[[n]][1:4,]
    slices <- apply(rna_class,1,sum) 
    pie(x = slices, main = names(lasser_laurent_all)[n])
    rm(rna_class, slices)
}

########################################################################################
#Filter and transform the raw exRNA data to 0-1 range so we meet
#assumptions of constrained matrix optimization employed by EDec

mir_dat <- lasser_laurent_all$miRNA.rpm

#Subset to include only cell supernatant samples and
#Select only those records with both readmap data and metadata
mir_dat <- mir_dat[intersect(rownames(mir_dat), rownames(lasser_laurent_meta)),]

##Whole matrix transformations##
mir_mtx <- as.matrix(mir_dat, drop = F)
mir_mtx_qn_var <- txfmQN(mir_mtx) 
mir_mtx_qn_obs <- t(txfmQN(t(mir_mtx))) #samples will have same distribution. This may eliminate 'crosstalk' between mapped reads 
mir_mtx_qn_obs.var <- txfmQN(mir_mtx_qn_obs) #samples and vars sequentially normalized within their own dimension
mir_mtx_qn_o.v_npn <- txfmNPN(mir_mtx_qn_obs.var)
#Within variables transformations
mir_mtx_final <- apply(mir_mtx_qn_obs.var , 2, txfmLogistic3)
#mir_mtx_final <- apply(mir_mtx_qn_o.v_npn , 2, txfmLogistic3)
#mir_mtx_final <- apply(mir_mtx_qn_obs , 2, txfmLogistic3)
#mir_mtx_final <- apply(mir_mtx_qn_var , 2, txfmLogistic3)
dimnames(mir_mtx_final) <- dimnames(mir_mtx)



boxplot(mir_mtx_final[1:30,],use.cols = F,
       col = c(rep("lightcoral",2),rep("lightpink",2),rep("lightcyan",26)))

boxplot(mir_mtx_final[1:20,1:60],use.cols = T)

##################
## Selection of informative exRNA
##################
#===============================================================================

#g <- apply(mir_mtx_final, 2, IQR)
#qtl <- quantile(g,probs = seq(0,1,.05))
#h <- apply(mir_mtx_final, 2, filterByIQR, min_iqr = qtl["80%"])
#var_filtr_mir_mtx <- mir_mtx_final[,!h]
#dim(var_filtr_mir_mtx)
#boxplot(var_filtr_mir_mtx[,sample(1:ncol(var_filtr_mir_mtx),size = 50)], use.cols = T,ylim = c(0,1))

hd_dif <- abs(mapply(`-`,mir_mtx_final[1,], mir_mtx_final[2,]))
ld_dif <- abs(mapply(`-`,mir_mtx_final[3,], mir_mtx_final[4,]))
sum_difs <- mapply(sum, hd_dif, ld_dif)
hd_avg <- apply(mir_mtx_final[1:2,], 2, mean)
ld_avg <- apply(mir_mtx_final[3:4,], 2, mean)
hd_ld_diff <- abs(ld_avg - hd_avg)
head(hd_ld_diff)
q <- quantile(hd_ld_diff)
length(hd_ld_diff[which(hd_ld_diff > q["75%"])])
tst <- rbind(hd_dif,ld_dif, sum_difs, hd_ld_diff)
tst <- tst[,order(-tst[4,],tst[3,],decreasing = F)]
tst <- rbind(tst,(tst[4,] / abs(range(mir_mtx_final)[2] - range(mir_mtx_final)[1])))
tst1 <- tst[,tst[5,] >= .3]
head(tst1[,1:12])
dim(mir_mtx_final[,colnames(tst1)])

pair_by_kit <- perform_t_tests_all_classes_each_pair(dataMatrix = mir_mtx_final,
                                    classVector = lasser_laurent_meta$Kit) 
one_v_all_by_kit <- perform_t_tests_all_classes_one_vs_rest(dataMatrix = mir_mtx_final,
                                    classVector = lasser_laurent_meta$Kit)

good_probes <- unique(as.vector((apply(pair_by_kit$P.Values, 2, function(x){head(order(x),20)} ))))

tmp <- mir_mtx_final[,good_probes]
head(tmp[,1:3])

#===============================================================================
library(RColorBrewer)
my_colors_fct <- with(laurent_meta_sup, data.frame(M_F.CellLine = levels(M_F.CellLine),
                                                   color = rainbow(nlevels(M_F.CellLine))))
my_colors_f <- merge(laurent_meta_sup, my_colors_fct)
my_colors <- as.character(my_colors_f$color)
gplots::heatmap.2(tmp, trace = "none") #RowSideColors = my_colors )
gplots::heatmap.2(tmp, trace = "none",Rowv = F,Colv = T) #RowSideColors = my_colors )

levels(laurent_meta_sup$Kit)
tsne_mir_mat <- tsne(mir_mtx_final[,colnames(tst1)],k = 8,perplexity = 8)
par(bg = "lightgrey")
plot(tsne_mir_mat[,1],tsne_mir_mat[,2], pch = 20, cex = 1,
     col = c(rep("magenta",2),rep("yellow",2),laurent_meta_sup$Lab))
cor(mir_dat[1,colnames(tst1)],mir_dat[2,colnames(tst1)])
cor(mir_mtx_final[3,],mir_mtx_final[4,])
cor(mir_mtx_final[1,],mir_mtx_final[2,])
cor(mir_mtx_final[1,],mir_mtx_final[2,])


#Load EDec Souce files provided by Oscar
dr <- "./input/Oscar/Glioblastoma/sources/"
fl <- as.list(list.files(dr))
fl <- sapply(dr,paste,fl, sep = "")
lapply(fl,source,echo = F)
#source("~/Documents/BCM/Lab/Glioblastoma/sources/diffExpEDec.R")
