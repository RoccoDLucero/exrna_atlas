
## PRE-PROCESS ATLAS METADATA TABLE TO ENSURE QUALITY
mta <- readRDS("../get_atlas_data_in_R/interim/comprehensive_atlas_metadata_draft_1.RDS")
mta <- mta[-(which(is.na(as.character(mta$Sample.Name)))),]
samps_mta <- sort(as.character(mta$Sample.Name))
#WHY is there a sample.name listed as NA in the metadata table at row 1397???

## PRE-PROCESS EACH ATLAS RPM TABLE TO ENSURE QUALITY AND RECTIFY IDENTIFIERS
b <- readRDS("../get_atlas_data_in_R/interim/exrna_atlas_readcounts_non_gencode.RDS")

##MIRNA##
rpm_tab <- b$miRNA
samps_b <- sort(gsub("\\.","-",rownames(rpm_tab)))
rownames(rpm_tab) <- samps_b
mir <- rpm_tab

#TRNA##
rpm_tab <- b$tRNA
samps_b <- sort(gsub("\\.","-",rownames(rpm_tab)))
rownames(rpm_tab) <- samps_b
tr <- rpm_tab

##PIRNA##
rpm_tab <- b$piRNA
samps_b <- sort(gsub("\\.","-",rownames(rpm_tab)))
rownames(rpm_tab) <- samps_b
rna_b <- gsub(pattern = "\\|.*$",replacement = "", colnames(rpm_tab))
colnames(rpm_tab) <- rna_b
pir <- rpm_tab

if(all(rownames(pir) == rownames(mir)) & all(rownames(pir) == rownames(tr))){samps_good <- T}

## ENSURE THAT ONLY SAMPLES WITH DATA AND METADATA ARE PASSED ON
if(samps_good){ all_samps <- intersect(rownames(mir), samps_mta)}

mir <- mir[all_samps,]
pir <- pir[all_samps,]
tr <- tr[all_samps,]

if(all(rownames(pir) == sort(as.character(mta$Sample.Name)))){

    mir <- cbind(as.character(mta$Biosample), rownames(mir), mir)
    colnames(mir)[1:2] <- c('Biosample','Study.Sample.Name')

    pir <- cbind(as.character(mta$Biosample), rownames(pir), pir)
    colnames(pir)[1:2] <- c('Biosample','Study.Sample.Name')

    tr <- cbind(as.character(mta$Biosample), rownames(tr), tr)
    colnames(tr)[1:2] <- c('Biosample','Study.Sample.Name')

}


write.table(a[,1:33],"../get_atlas_data_in_R/interim/for_william/atlas_meta_test.txt",
            quote = F, sep = "\t", eol = "\r", row.names = F, col.names = T)

write.table(mir, "../get_atlas_data_in_R/interim/for_william/atlas_mirna_rpm_test.txt",
            quote = F, sep = "\t", eol = "\r", row.names = F, col.names = T)

write.table(pir, "../get_atlas_data_in_R/interim/for_william/atlas_pirna_rpm_test.txt",
            quote = F, sep = "\t", eol = "\r", row.names = F, col.names = T)

write.table(tr, "../get_atlas_data_in_R/interim/for_william/atlas_trna_rpm_test.txt",
            quote = F, sep = "\t", eol = "\r", row.names = F, col.names = T)






