



bsid_mp <- get.combined.bsid.to.smp.name.maps(studies_url,studies_dirs = dirs_studies)

mt[["bsid_sample_map"]] <- bsid_mp

saveRDS(object = mt, file = '../get_atlas_data_in_R/interim/exrna_atlas_metadata_and_map_to_samples.2.RDS')

head(bsid_mp)

####################################################################################################


meta <- readRDS(file = '../get_atlas_data_in_R/interim/exrna_atlas_metadata_and_map_to_samples.RDS')
cmbnd_meta <- lapply(X = meta_types,
                   FUN = function(dt){
                       d <- get.doctype.tables(atlas_metadata = bsid ,studies = dirs_studies, doc_type = dt)
                   }
              )
names(cmbnd_meta) <- meta_types



