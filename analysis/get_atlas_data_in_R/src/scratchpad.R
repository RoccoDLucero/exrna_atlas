

#rna_types <- c('gencode')


tmp_fldr <- '../get_atlas_data_in_R/interim/exrna_atlas_rc_tmp/'
save_path <- '../get_atlas_data_in_R/interim/'



rna_types <- c('miRNA', 'piRNA')
get_data_args <- list(studies_url = studies_url, studies_dirs = dirs_studies,
                      rna_types = rna_types,
                      tmp_fldr = tmp_fldr )

rds_save_output(fun = aggregate.exrna.atlas.readcounts,
                args = get_data_args, save_path = save_path,
                save_as_name = 'exrna_atlas_readcounts_non_gencode.RDS')



rna_types <- c('gencode')
get_data_args <- list(studies_url = studies_url, studies_dirs = dirs_studies,
                      rna_types = rna_types,
                      tmp_fldr = tmp_fldr )

rds_save_output(fun = aggregate.exrna.atlas.readcounts,
                args = get_data_args, save_path = save_path,
                save_as_name = 'exrna_atlas_readcounts_non_gencode.RDS')
