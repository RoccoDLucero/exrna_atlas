###################################################################################################
##title:                 exrna_atlas_merge_data_and_metadata.R
##associated_project:    exrna_atlas
##associated_analysis:   get_atlas_data_in_R
##started_date:          April 2017
##started_by:            Rocco Lucero
##last_updated_date:     June 21, 2017
##last_updated_by:       Rocco Lucero
###################################################################################################
##description:   Pull together the exrna_atlas readcount data and metatdata into
##               a single comprehensive R object
##
##inputs:                <What goes in>
##
##outputs:               <What are the outputs>
##
##dependencies:          <exrna_atlas_get_data_functions.R; exrna_atlas_get_metadata_functions.R
##                        get_atlas_data_in_R_generic_functions.R>
###################################################################################################


source(file = "./src/exrna_atlas_get_data_functions.R")
source(file = "./src/exrna_atlas_get_metadata_functions.R")


studies_url <- paste('ftp://ftp.genboree.org/exRNA-atlas/grp/',
                          'Extracellular RNA Atlas/db/exRNA Repository - hg19/file/',
                          'exRNA-atlas/exceRptPipeline_v4.6.2/', sep = '')

dirs_studies <- my.get.print.studies(studies_url = studies_url, print_only = F)

rna_types1 <- c('miRNA', 'tRNA', 'piRNA')
rna_types2 <- c('gencode')
rna_types3 <- c('circularRNA')

meta_types <- c('BS', 'DO', 'ST', 'SU')

get_rc <- F
if(get_rc){

    obj <- do.call(what = my.get.atlas.readcounts,
                   args = list(rna_types = rna_types1,
                               studies_url = studies_url,
                               dirs_studies = dirs_studies[1]))

exrna_atlas_readcounts1 <- my.get.atlas.readcounts(rna_types = rna_types1,
                                                  studies_url = studies_url,
                                                  dirs_studies = dirs_studies)

saveRDS(object = exrna_atlas_readcounts1, file = "./exrna_atlas_readcounts_non_gencode.Rdat")
rm(exrna_atlas_readcounts1)

exrna_atlas_readcounts2 <- my.get.atlas.readcounts(rna_types = rna_types2,
                                                  studies_url = studies_url,
                                                  dirs_studies = dirs_studies)

saveRDS(object = exrna_atlas_readcounts2,
        file = "./output/exrna_atlas_gencode_readcounts.RDS")

exrna_atlas_readcounts2_sep <- sep.df.by.colname.pttn(exrna_atlas_readcounts2$gencode)

rm(exrna_atlas_readcounts2)

saveRDS(object = exrna_atlas_readcounts2_sep,
        file = "./output/exrna_atlas_split_from_gencode_readcounts.RDS")

rm(exrna_atlas_readcounts2_sep)
}

get_meta <- T
if(get_meta){

    exrna_atlas_meta <- get.atlas.metadata(studies_url = studies_url,
                                           studies_dirs = dirs_studies,
                                           meta_types = meta_types)


}

s <- get.atlas.metadata(studies_url = studies_url, studies_dirs = dirs_studies[1:2], meta_types = meta_types)
