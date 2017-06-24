################################################################################
#Figure out how to pass arguments without triggering recursive reference warnings



write.study.readcounts.obj.by.rna <-function(save_fldr, study_dirname,
                                             rna_type, ...){


    study_ppr <- my.get.study.post.processed.results(studies_url = studies_url,
                                                     study_dirname = study_dirname,
                                                     postprocess_version = postprocess_version,
                                                     rna_type = rna_type,
                                                     check_names = check_names,
                                                     per_million = per_million)

    save_name <- paste(study_dirname, "_", rna_type, '_readcounts.RDS')
    save_name <- gsub(" ", "", save_name)

    save_path <- paste(save_fldr, save_name, sep = '')

    saveRDS(study_ppr, save_path)

    return(NULL)
}


combine.on.disk.study.readcounts.by.rna <- function(rna_type, path = save_fldr){

    filter.empty.df <- function(df){ is.null(df)}

    #path <- '../get_atlas_data_in_R/interim/exrna_atlas_rc_tmp/'

    rc_tables <- list.files(path = path)

    tmp_list <- vector("list", length(rc_tables))

    names(tmp_list) <- rc_tables

    while(length(tmp_list) > 1){

        if(names(tmp_list)[1] == rc_tables[1]){

            df_1 <- readRDS(paste(path,names(tmp_list)[1],sep = ''))

        }else{

            df_1 <- tmp_list[[1]]

        }

        df_2 <- readRDS(paste(path,names(tmp_list)[2],sep = ''))

        if(filter.empty.df(df_1)){ tmp_list <- tmp_list[-1]; next }

        if(filter.empty.df(df_2)){ tmp_list <- tmp_list[-2]; next }

        new_df <- my.combine.atlas.rpm.tables(df_1, df_2)

        tmp_list[[1]] <- new_df

        names(tmp_list)[1] <- "new_df"

        tmp_list <- tmp_list[-2]

    }

    df_cmbnd <- tmp_list[[1]]

    df_cmbnd <- t(df_cmbnd)

    return(df_cmbnd)

}


get.and.combine.readcounts.by.rna <- function(studies_dirs = studies_dirs,
                                                 rna_type = rna_type,
                                                 save_fldr = save_fldr){
    dir.create(save_fldr, showWarnings = F)

    sapply(X = studies_dirs,
           FUN = function(st){
               write.study.readcounts.obj.by.rna(study_dirname = st, rna_type = rna_type)})

    rt_combnd <- combine.on.disk.study.readcounts.by.rna(rna_type = rna_type)

    unlink(x = save_fldr, recursive = T)

    return(rt_combnd)
}


aggregate.exrna.atlas.readcounts <- function(studies_dirs, rna_types, save_fldr,
                                             studies_url, per_million = T, check_names = T,
                                             postprocess_version = "v4.6.3" ){

    dir.create(save_fldr, showWarnings = F)

    rc_lst <- lapply(X = rna_types,
                     FUN = function(rt){ get.and.combine.readcounts.by.rna(rt)},
                     studies_dirs = studies_dirs)

    names(rc_lst) <- rna_types

    rc_lst <- lapply(X = rc_lst, as.data.frame )

}

#rna_types <- c('gencode')
rna_types <- c('miRNA', 'tRNA', 'piRNA')
save_fldr <- '../get_atlas_data_in_R/interim/exrna_atlas_rc_tmp/'


write.study.readcounts.obj.by.rna(save_fldr = save_fldr, study_dirname = dirs_studies[1],
                                  rna_type = rna_types[1])

get.and.combine.readcounts.by.rna(studies_dirs = dirs_studies[1:2], rna_type = rna_type, save_fldr = save_fldr)










