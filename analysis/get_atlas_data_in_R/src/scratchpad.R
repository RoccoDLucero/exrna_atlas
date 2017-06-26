

gs <- readRDS('../get_atlas_data_in_R/interim/exrna_atlas_readcounts_gencode_split.RDS')
sapply(X = gs, FUN = dim)
head(gs$snRNA)[,1:4]

gs <- readRDS('../get_atlas_data_in_R/interim/exrna_atlas_readcounts_gencode_mixed.RDS')
sapply(X = gs, FUN = dim)
head(gs$gencode)[,1:4]

gs <- readRDS('../get_atlas_data_in_R/interim/exrna_atlas_readcounts_non_gencode.RDS')
sapply(X = gs, FUN = dim)
head(gs$miRNA)[,1:4]

mt <- readRDS('../get_atlas_data_in_R/interim/exrna_atlas_metadata_and_map_to_samples.RDS')
sapply(X = mt, FUN = length)


rm(gs)
gc()
