###################################################################################################
##title:                 exrna_atlas_merge_data_and_metadata.R
##associated_project:    exrna_atlas
##associated_analysis:   get_atlas_data_in_R
##started_date:          April 2017
##started_by:            Rocco Lucero
##last_updated_date:     June 24, 2017
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
source(file = "../get_atlas_data_in_R/src/exrna_atlas_get_data_functions.R")
source(file = "../get_atlas_data_in_R/src/exrna_atlas_get_metadata_functions.R")


#####SET PRELIMINARY PARAMETERS#####
save_path <- '../get_atlas_data_in_R/interim/'

studies_url <- paste('ftp://ftp.genboree.org/exRNA-atlas/grp/',
                          'Extracellular RNA Atlas/db/exRNA Repository - hg19/file/',
                          'exRNA-atlas/exceRptPipeline_v4.6.2/', sep = '')

dirs_studies <- my.get.print.studies(studies_url = studies_url, print_only = F)

meta_types <- c('BS', 'DO', 'ST', 'SU')

gencode_biotypes <- get.atlas.biotype.names(studies_url, dirs_studies)

################################################################################
get_meta <- F
if(get_meta){

    get_meta_args <- list(studies_url = studies_url,
                          studies_dirs = dirs_studies,
                          meta_types = meta_types,
                          get_map = F )

    rds_save_output(fun = get.atlas.metadata, args = get_meta_args,
                    save_path = save_path,
                    save_as_name = 'exrna_atlas_metadata_and_map_to_samples.RDS')
}

get_rc <- F
if(get_rc){

    #Download non-gencode readcount data and store in a single R object
    get_data_args <- list(rna_types = NULL, studies_url = studies_url, dirs_studies = dirs_studies)

    get_data_args$rna_types <- c('miRNA', 'tRNA', 'piRNA')

    rds_save_output(fun = my.get.atlas.readcounts, args = get_data_args, save_path = save_path,
                    save_as_name = 'exrna_atlas_readcounts_non_gencode.RDS')



    #Download gencode mapped readcount data and store in a single R object
    get_data_args$rna_types <- c('gencode')

    rds_save_output(fun = my.get.atlas.readcounts, args = get_data_args, save_path = save_path,
                    save_as_name = 'exrna_atlas_readcounts_gencode.RDS')

    gencode_mixed <- readRDS('./interim/exrna_atlas_readcounts_gencode.RDS')

    gencode_mixed <- gencode_mixed$gencode

    split_gencode_args <-  list(df = gencode_mixed, patterns = gencode_biotypes)

    rds_save_output(fun = sep.df.by.colname.pttn, args = split_gencode_args,
                    save_path = save_path,
                    save_as_name = 'exrna_atlas_readcounts_gencode_split.RDS')

    gencode_split <- readRDS('./interim/exrna_atlas_readcounts_gencode_split.RDS')

    rm(gencode_mixed, split_gencode_args, get_data_args)

}




