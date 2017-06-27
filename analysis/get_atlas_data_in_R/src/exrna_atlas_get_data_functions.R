###################################################################################################
##title:                 exrna_atlas_get_data_functions.R
##associated_project:    exrna_atlas
##associated_analysis:   get_atlas_data_in_R
##started_date:          April 2017
##started_by:            Rocco Lucero
##last_updated_date:     June 22, 2017
##last_updated_by:       Rocco Lucero
###################################################################################################
##description:   functions for retrieval and processing exrna atlas read count data to return
##               a single comprehensive R object
##
##inputs:                <NA>
##
##outputs:               <NA>
##
##dependencies:          <exrna_atlas_generic_functions.R>
###################################################################################################
##Table of Contents:
##.....GET AND PROCESS READ COUNTS.....
##--my.get.study.post.processed.results
##--my.combine.study.counts.by.rna
##--my.get.atlas.readcounts (MAIN)
################################################################################
################################################################################

required_packages <- c("RCurl", "stringr")

lapply(X = required_packages, FUN = library, character.only = T)

source("../get_atlas_data_in_R/src/get_atlas_data_in_R_generic_functions.R")

################################################################################
#Functions to interface with the ExRNA Atlas Data on the BRL FTP site
################################################################################

################################################################################
#Get RNA-seq readcount data from a given study in the ExRNA-Atlas stored on BRL FTP-site.

my.get.study.post.processed.results <- function(studies_url, study_dirname,
                                                postprocess_version = "v4.6.3",
                                                rna_type = 'miRNA', check_names = F,
                                                per_million = T){

    ############################################################################
    #Input: (4 character)
    #       The URL to Genboree FTP StudiesExceRpt Jobs;
    #       name of study; name of an RNA type mapped in the ExceRpt Run;
    #       version number for Excerpt post processing scripts
    #
    #       e.g. \"ftp://ftp.genboree.org/exRNA-atlas/grp/
    #                   Extracellular RNA Atlas/db/exRNA Repository - hg19/
    #                   file/exRNA-atlas/exceRptPipeline_v4.6.2/\";
    #             \"KJENS1-RIDProject-2016-06-27\" ; \"gencode\"; \"v4.6.3\"
    #
    #Output: (single R object)
    #       A dataframe of RNA readcounts for all samples in the given study
    #
    #Called by: my.combine.study.reads
    ############################################################################

    if(per_million){

        normalization <- "ReadsPerMillion.txt"

    }else{

        normalization <- 'ReadCounts.txt'

    }

    study_url <- paste(studies_url, study_dirname,'/', sep = '')

    dirs_data <- retry.connection(target_url = study_url, dirs_only = T, text_connect = F)
    dirs_data <- unlist(strsplit(x = dirs_data, split = '\r\n'))

    results_ver <- paste('postProcessedResults_', postprocess_version, sep = "")

    results_dir <- grep(results_ver, dirs_data, value = T)

    ppr_url <- paste(study_url, results_dir,'/', sep = '')

    dirs_data <- retry.connection(target_url = ppr_url, dirs_only = T, text_connect = F)
    dirs_data <- unlist(strsplit(dirs_data,'\r\n'))

    ptt <- paste("exceRpt", rna_type, normalization, sep = "_")

    data_file <- grep(ptt, dirs_data, value = T)

    targ_file_url <- paste(ppr_url, data_file, sep = '')

    cat("\r\nRetrieving: ", data_file)

    cat("\r\nFrom: ", targ_file_url)

    reads_dat <- retry.connection(targ_file_url, text_connect = T)

    if( any(dim(reads_dat) == 0) ){

        cat(study_dirname, " had no data in ", data_file )

        return(NULL)

    }

    rownames(reads_dat) <- reads_dat[,1]

    reads_dat <- remove.columns(df = reads_dat, remove_cols = c(1))

    Sys.sleep(.1)

    return( reads_dat )

}

################################################################################
################################################################################
##THESE FUNCTIONS GET READCOUND DATA BUT ONLY USE MEMORY TO PROCESS FILES
################################################################################
################################################################################
my.combine.study.counts.by.rna <- function( dirs_studies, rna_type, studies_url,
                                            per_million = T, check_names = T){

    get.ppr <- function(st){
        my.get.study.post.processed.results(studies_url = studies_url,
                                            rna_type = rna_type,
                                            study_dirname = st)
    }

    filter.empty.df <- function(df){ is.null(df)}

    tmp_list <- vector("list", length(dirs_studies))
    names(tmp_list) <- dirs_studies

    while(length(tmp_list) > 1){

        if(names(tmp_list)[1] == dirs_studies[1]){

            df_1 <- get.ppr(st = names(tmp_list)[1])

        }else{

            df_1 <- tmp_list[[1]]

        }

        df_2 <- get.ppr(st = names(tmp_list)[2])

        if(filter.empty.df(df_1)){ tmp_list <- tmp_list[-1]; next }

        if(filter.empty.df(df_2)){ tmp_list <- tmp_list[-2]; next }

        new_df <- my.combine.atlas.rpm.tables(df_1, df_2)

        tmp_list[[1]] <- new_df

        names(tmp_list)[1] <- "new_df"

        tmp_list <- tmp_list[-2]

        #print(names(tmp_list))

    }

    df_cmbnd <- tmp_list[[1]]
    #df_cmbnd <- remove.columns(df_cmbnd, remove_cols = 1)
    df_cmbnd <- t(df_cmbnd)

    #rn <- rownames(df_cmbnd)

    #df_cmbnd <- as.data.frame(df_cmbnd)
    #df_cmbnd <- apply(X = df_cmbnd, MARGIN = 2, FUN = as.numeric)
    #rownames(df_cmbnd) <- rn

    return(df_cmbnd)
    #return( list(df_cmbnd, df_1, df_2 ))

}

################################################################################
#THIS IS THE MAIN FUNCTION USED TO GET READCOUND DATA
################################################################################
my.get.atlas.readcounts <- function(rna_types, studies_dirs, studies_url,
                                    per_million = T, check_names = T){

    cat("\n", rep("=", 60), "\n", sep = "")
    cat("\n", "Acquiring read count for:", "\n")
    print(dirs_studies)
    cat("\n", "with selected RNA types:", "\n")
    cat(unlist(rna_types))
    cat("\n", rep("=", 60), "\n", sep = "")

    rc <- lapply(X = rna_types,
                 FUN = my.combine.study.counts.by.rna, dirs_studies = dirs_studies,
                 studies_url = studies_url, per_million = per_million,
                 check_names = check_names)

    names(rc) <- rna_types

    rc <- lapply(X = rc, as.data.frame )

    rc["Access Date"] <- paste( "This record was created on:", Sys.time(), sep = " ")

    cat("Timestamp Added.")

    gc()

    return(rc)

}

################################################################################
################################################################################
################################################################################
##THE FUNCTIONS BELOW ACCOMPLISH THE SAME AS THOSE ABOVE, BUT THEY USE THE DISK
##AS A TEMPORARY STORAGE OF UNPROCESSED FILES
################################################################################
write.study.readcounts.obj.by.rna <-function(save_fldr, study_dirname,
                                             rna_type, ...){


    study_ppr <- my.get.study.post.processed.results(studies_url, study_dirname,
                                                     rna_type,
                                                     postprocess_version = "v4.6.3",
                                                     check_names = F,
                                                     per_million = T)

    save_name <- paste(study_dirname, "_", rna_type, '_readcounts.RDS')
    save_name <- gsub(" ", "", save_name)

    save_path <- paste(save_fldr, save_name, sep = '')

    saveRDS(study_ppr, save_path)

    return(NULL)
}


combine.on.disk.study.readcounts.by.rna <- function(rna_type,
                                                    save_fldr, ...){
    path <- save_fldr

    filter.empty.df <- function(df){ is.null(df)}

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


get.and.combine.readcounts.by.rna <- function(studies_dirs,
                                              rna_type,
                                              save_fldr,
                                              unlink = T, ...){
    dir.create(save_fldr, showWarnings = F)

    sapply(X = studies_dirs,
           FUN = function(st){
               write.study.readcounts.obj.by.rna(study_dirname = st,
                                                 rna_type = rna_type,
                                                 save_fldr = save_fldr)
           })

    rt_combnd <- combine.on.disk.study.readcounts.by.rna(rna_type, save_fldr, unlink)

    unlink_fldr <- gsub(pattern = "/$", replacement = '', x = save_fldr)

    if(unlink){unlink(x = unlink_fldr, recursive = T)

        return(rt_combnd)

    }else{return(NULL)}


}

################################################################################
#THIS IS THE ALTERNATIVE MAIN FUNCTION USED TO GET READCOUND DATA
################################################################################
aggregate.exrna.atlas.readcounts <- function( studies_url,
                                              studies_dirs,
                                              rna_types,
                                              tmp_fldr, ...){
    save_fldr <- tmp_fldr

    dir.create(save_fldr, showWarnings = F)

    rc_lst <- lapply(X = rna_types,
                     FUN = function(rt){
                         get.and.combine.readcounts.by.rna(rna_type = rt,
                                                           studies_dirs,
                                                           save_fldr)})

    names(rc_lst) <- rna_types

    rc_lst <- lapply(X = rc_lst, as.data.frame )

    print(sapply(X = rc_lst, FUN = dim))

    return(rc_lst)

}




