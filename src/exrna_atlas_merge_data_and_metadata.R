###################################################################################################
##title:                 exrna_atlas_merge_data_and_metadata.R
##started_date:          unknown
##started_by:            Rocco Lucero
##last_updated_date:     June 8, 2017
##last_updated_by:       Rocco Lucero
###################################################################################################
##description:   Pull together the exrna_atlas readcount data and metatdata into
##               a single comprehensive R object
##
##inputs:                <What goes in>
##
##outputs:               <What are the outputs>
##
##dependencies:          <Name any dependencies not publicly available>
###################################################################################################


source(file = "./exrna_atlas_get_read_counts_functions.R")
source(file = "./exrna_atlas_get_meta_data_functions.R")


studies.url <- paste('ftp://ftp.genboree.org/exRNA-atlas/grp/',
                          'Extracellular RNA Atlas/db/exRNA Repository - hg19/file/',
                          'exRNA-atlas/exceRptPipeline_v4.6.2/', sep = '')

dirs.studies <- my.get.print.studies(studies_url = studies.url)

rna.types1 <- c('miRNA', 'tRNA', 'piRNA')
rna.types2 <- c('gencode')
rna.types3 <- c('circularRNA')

meta.types <- c('BS', 'DO', 'ST', 'SU')

get_rc <- F
if(get_rc){
exrna_atlas_readcounts1 <- my.get.atlas.readcounts(rna_types = rna.types1,
                                                  studies_url = studies.url,
                                                  dirs_studies = dirs.studies)

saveRDS(object = exrna_atlas_readcounts1, file = "./exrna_atlas_readcounts_non_gencode.Rdat")
rm(exrna_atlas_readcounts1)

exrna_atlas_readcounts2 <- my.get.atlas.readcounts(rna_types = rna.types2,
                                                  studies_url = studies.url,
                                                  dirs_studies = dirs.studies)

saveRDS(object = exrna_atlas_readcounts2,
        file = "./exrna_atlas_gencode_readcounts.Rdat")

exrna_atlas_readcounts2_sep <- sep.df.by.colname.pttn(exrna_atlas_readcounts2$gencode)

rm(exrna_atlas_readcounts2)

saveRDS(object = exrna_atlas_readcounts2_sep,
        file = "./exrna_atlas_split_from_gencode_readcounts.Rdat")

rm(exrna_atlas_readcounts2_sep)
}

get_meta <- T
if(get_meta){

mm <- my.get.study.metadata(studies_url = studies.url,
                      study_dirname = dirs.studies[1])


#exrna_atlas_metadata <- my.merge.metadata(studies_url = studies.url,
#                                          studies_dirs = dirs.studies[1:2])

saveRDS(object = exrna_atlas_metadata, file = "./exrna_atlas_metadata.Rdat")
rm(exrna_atlas_metadata)

exrna_atlas_bsid_to_sample_names <- lapply(X = dirs.studies,
                                           FUN = function(s){})
}


