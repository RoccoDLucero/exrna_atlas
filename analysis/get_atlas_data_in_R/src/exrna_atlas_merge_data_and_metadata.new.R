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
source(file = "../get_atlas_data_in_R/src/get_atlas_data_in_r_src_functions.R")


#####SET PRELIMINARY PARAMETERS#####
tmp_fldr <- '../get_atlas_data_in_R/interim/exrna_atlas_rc_tmp/'
save_path <- '../get_atlas_data_in_R/interim/'

studies_url <- paste('ftp://ftp.genboree.org/exRNA-atlas/grp/',
                          'Extracellular RNA Atlas/db/exRNA Repository - hg19/file/',
                          'exRNA-atlas/exceRptPipeline_v4.6.2/', sep = '')

dirs_studies <- my.get.print.studies(studies_url = studies_url, print_only = F)

################################################################################
################################################################################
meta_types <- c('BS', 'DO', 'EX', 'RU', 'ST', 'SU', 'AN')

meta_rds_file <- 'exrna_atlas_metadata_no_map_to_samples.RDS'
get_meta <- F
if(get_meta){



    get_meta_args <- list(studies_url = studies_url,
                          studies_dirs = dirs_studies,
                          meta_types = meta_types,
                          get_map = T )

    rds_save_output(fun = get.atlas.metadata, args = get_meta_args,
                    save_path = save_path,
                    save_as_name = meta_rds_file)

}

process_meta <- F
if(process_meta){

    meta_path <- paste(save_path, meta_rds_file, sep = "")

    meta <- readRDS(file = meta_path)

    #meta_map <- meta$bsid_sample_map

    #meta_map <- reformat.names(meta_map,2)

    cmbnd_meta <- lapply(X = meta_types,
                     FUN = function(dt){


                         d <- get.doctype.tables(atlas_metadata = meta,
                                                 studies = dirs_studies,
                                                 doc_type = dt, universal_fields = T)
                         dim <- 1

                         d <- drop.field(d, dim, 'DocURL')
                         d <- drop.field(d, dim, '-+ Type')
                         d <- drop.field(d, dim, 'Status')
                         d <- drop.field(d, dim, 'MD5.*')
                         d <- drop.field(d, dim, 'File.Name')

                         d <- reformat.names(d, dim)



                         return(t(d))

                     }
    )

    names(cmbnd_meta) <- meta_types


    ds <- as.data.frame(dirs_studies)
    colnames(ds) <- "Study.Name.FTP.Dir"
    cmbnd_meta$ST <- cbind(cmbnd_meta$ST, ds)


    colnames(cmbnd_meta$BS) <- change.colnames(cmbnd_meta$BS,
                                               "Donor.ID", "Donor")

    colnames(cmbnd_meta$BS) <- change.colnames(cmbnd_meta$BS,
                                               "Related.Experiment$", "Experiment")

    colnames(cmbnd_meta$RU) <- change.colnames(cmbnd_meta$RU,
                                               "Related.Study", "Study")

    colnames(cmbnd_meta$RU) <- change.colnames(cmbnd_meta$RU,
                                               "Biosample.ID", "Biosample")

    colnames(cmbnd_meta$ST) <- change.colnames(cmbnd_meta$ST,
                                               "Related.Submission$", "Submission")

    colnames(meta_map) <- change.colnames(meta_map,
                                          "BS.ID", "Biosample")

    sapply(X = cmbnd_meta, colnames)

    atlas_core_metadata <- merge.data.frame(x = cmbnd_meta$BS,
                                            y = cmbnd_meta$DO,
                                            by = 'Donor',
                                            all.x = T)

    atlas_core_metadata <- merge.data.frame(x = atlas_core_metadata,
                                            y = cmbnd_meta$EX,
                                            by = 'Experiment',
                                            all.x = T)

    atlas_core_metadata <- merge.data.frame(x = atlas_core_metadata,
                                            y = cmbnd_meta$RU,
                                            by = "Biosample",
                                            all.x = T)

    atlas_core_metadata <- merge.data.frame(x = atlas_core_metadata,
                                            y = cmbnd_meta$ST,
                                            by = 'Study',
                                            all.x = T)

    atlas_core_metadata <- merge.data.frame(x = atlas_core_metadata,
                                            y = cmbnd_meta$SU,
                                            by = 'Submission',
                                            all.x = T)

    atlas_core_metadata <- merge.data.frame(x = atlas_core_metadata,
                                            y = meta_map[,c(1,2,4)],
                                            by = 'Biosample',
                                            all.x = T)

    colnames(atlas_core_metadata) <- gsub(pattern = '^X[.]+',
                                          replacement = '',
                                          colnames(atlas_core_metadata))


    #atlas_core_metadata <-  atlas_core_metadata[ ,c(1, 23, 4, 24, 2, 19, 15, 21, 22, 11, 13, 14, 17, 18, 25) ]


    saveRDS(atlas_core_metadata,
            "../get_atlas_data_in_R/interim/comprehensive_atlas_metadata_draft_1.RDS")
}



################################################################################

################################################################################
get_rc <- F
get_non_gencode <- F
get_gencode <- F
split_gencode <- F
get_circ <- F
if(get_rc){

    if(get_non_gencode){

        rna_types <- c('miRNA', 'piRNA', 'tRNA')
        get_data_args <- list(studies_url = studies_url, studies_dirs = dirs_studies,
                              rna_types = rna_types,
                              tmp_fldr = tmp_fldr )

        rds_save_output(fun = aggregate.exrna.atlas.readcounts,
                        args = get_data_args, save_path = save_path,
                        save_as_name = 'exrna_atlas_readcounts_non_gencode.RDS')

        rm(rna_types, get_data_args)
    }


    if(get_gencode){

        rna_types <- c('gencode')
        get_data_args <- list(studies_url = studies_url, studies_dirs = dirs_studies,
                              rna_types = rna_types,
                              tmp_fldr = tmp_fldr )

        rds_save_output(fun = aggregate.exrna.atlas.readcounts,
                        args = get_data_args, save_path = save_path,
                        save_as_name = 'exrna_atlas_readcounts_gencode_mixed.RDS')

        rm(rna_types, get_data_args)
    }

    if(split_gencode){

        gencode_mixed <- readRDS('./interim/exrna_atlas_readcounts_gencode_mixed.RDS')
        gencode_mixed <- gencode_mixed$gencode

        gencode_biotypes <- get.atlas.biotype.names(studies_url, dirs_studies)

        split_gencode_args <-  list(df = gencode_mixed, patterns = gencode_biotypes)

        rds_save_output(fun = sep.df.by.colname.pttn, args = split_gencode_args,
                        save_path = save_path,
                        save_as_name = 'exrna_atlas_readcounts_gencode_split.RDS')

        rm(gencode_mixed, split_gencode_args, gencode_biotypes)
    }

    if(get_circ){

        rna_types <- c('circRNA')
        get_data_args <- list(studies_url = studies_url, studies_dirs = dirs_studies,
                              rna_types = rna_types,
                              tmp_fldr = tmp_fldr )

        rds_save_output(fun = aggregate.exrna.atlas.readcounts,
                        args = get_data_args, save_path = save_path,
                        save_as_name = 'exrna_atlas_readcounts_circRNA_mixed.RDS')

        rm(rna_types, get_data_args)
    }

}







