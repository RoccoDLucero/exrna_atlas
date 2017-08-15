

####################################################################################################
source("https://bioconductor.org/biocLite.R")
biocLite('GEOquery', suppressUpdates = T)
library('GEOquery')
####################################################################################################
####################################################################################################
#####GET NANOSTRING DATA FROM BLOOD CELLS GSE57649##
blood_nano_dat <- getGEO(GEO = "GSE57649")

blood_mirna_expression <- blood_nano_dat$GSE57649_series_matrix.txt.gz@assayData$exprs

blood_mirna_meta <- blood_nano_dat$GSE57649_series_matrix.txt.gz@phenoData@data

nano_informative_probes <- perform.probe.selection(ref_probes =  blood_mirna_expression,
                                              ref_classes = blood_mirna_meta$source_name_ch1,
                                              n_markers = 60,
                                              max_p_val = probe_sel_max_p)


#selected_probes <- select.best.probes(informative_probes, probes_per_comparison = 15)
tmp1 <- lapply(c(1:4), function(x){head(sort(nano_informative_probes$one.vs.rest$t_test$p.values[,x]),30)})
tmp2 <- unique(names(unlist(tmp1)))
informative_probes1 <- tmp2

tmp11 <- lapply(c(1:4), function(x){head(sort(nano_informative_probes$each.pair$t_test$p.values[,x]),30)})
tmp22 <- unique(names(unlist(tmp11)))
informative_probes2 <- tmp22

nano_informative_probes <- unlist(union(informative_probes1, informative_probes2))

#This is a stop-gap that removes probes not listed in the rna-seq-data for the tissue references
rem_prbs <- setdiff(nano_informative_probes, rownames(edec_input_data$ref) )
nano_prbs <- setdiff(nano_informative_probes, rem_prbs)

inf_prbs <- union(nano_prbs, informative_probes)
expression_correlations_refs <- list(test = cor(edec_input_data$ref[inf_prbs,]))


