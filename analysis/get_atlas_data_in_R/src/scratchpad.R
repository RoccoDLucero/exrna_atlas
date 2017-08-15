meta_d  <- readRDS("../get_atlas_data_in_R/interim/comprehensive_atlas_metadata_draft_1.RDS")

kj_meta <- meta_d[grep("kjens",meta_d$Biosample,ignore.case = T),]

p <- levels(droplevels(kj_meta$Study))
alz_prk <- head(kj_meta[kj_meta$Study == p[1],])
sah <- head(kj_meta[kj_meta$Study == p[2],])
ivh <- head(kj_meta[kj_meta$Study == p[3],])
rid <- head(kj_meta[kj_meta$Study == p[4],])
