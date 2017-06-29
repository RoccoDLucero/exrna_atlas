###################################################################################################
##title:                 exrna_atlas_get_metadata_functions.R
##associated_project:    exrna_atlas
##associated_analysis:   get_atlas_data_in_R
##started_date:          April 2017
##started_by:            Rocco Lucero
##last_updated_date:     June 22, 2017
##last_updated_by:       Rocco Lucero
###################################################################################################
##description:   functions for retrieving and processing exrna atlas metadata to return
##               a single comprehensive R object
##
##inputs:                <NA>
##
##outputs:               <NA>
##
##dependencies:          <get_atlas_data_in_R_generic_functions.R>
###################################################################################################
##Table of Contents:
##.....GET AND PROCESS METADATA.....
##--my.get.study.bioGPS.metadata
##--my.map.BSIDtoSampleName
##       --get.BS.smp.ID
##--get.combined.bsid.to.smp.name.maps(MAIN)
##--my.get.study.metadata
##       --get.BS.smp.ID
##       --extract.metadata
##               --extract.meta
##--my.merge.metadata (MAIN)

source("../get_atlas_data_in_R/src/get_atlas_data_in_R_generic_functions.R")

################################################################################
#FUNCTIONS FOR RETRIEVING AND PROCESSING ATLAS METADATA
################################################################################
################################################################################

get.study.biotypes <- function(studies_url, study_dirname){

    study_url <- paste(studies_url, study_dirname, '/', sep = '')

    dirs_data <- retry.connection(target_url = study_url, dirs_only = T, text_connect = F)
    dirs_data <- unlist(strsplit(dirs_data,'\r\n'))

    ppr_dir <- grep('postProcessedResults_', dirs_data, value = T)
    ppr_url  <- paste(study_url, ppr_dir,'/', sep = '')

    dirs_data <- retry.connection(target_url =  ppr_url, dirs_only = T, text_connect = F)
    dirs_data <- unlist(strsplit(dirs_data,'\r\n'))

    data_file <- grep('biotypeCounts.t[xs][tv]', dirs_data, value = T)

    targ_file_url <- paste(ppr_url, data_file, sep = '')

    cat("\r\nRetrieving: ", data_file, "\r\n")

    reads_dat <- retry.connection(targ_file_url, text_connect = T)

    if( any(dim(reads_dat) == 0) ){

        cat(study_dirname, " had no data in ", data_file )

        return(NULL)

    }

    return(reads_dat[,1])
}

get.atlas.biotype.names <- function(studies_url, studies_dirs){
    bt <-   sapply(X = studies_dirs,
            FUN = get.study.biotypes, studies_url = studies_url)

    bt <- unique(unlist(bt))

    return(bt)
}
################################################################################

my.get.study.biogps.metadata <- function(studies.url, study.dirname){
    ############################################################################
    #Input: (2 character)   URL to all public Atlas studies;
    #                       name of specific study
    #
    #Output: (single R object) A dataframe containing the BIOGPS metadata
    ############################################################################

    my.url <- paste(studies.url, study.dirname,'/', sep = '')
    dirs.data <- retry.connection(target_url = my.url, dirs_only = T, text_connect = F)
    dirs.data <- unlist(strsplit(dirs.data,'\r\n'))

    BioGPS.dir <- grep('BioGPS', dirs.data, value = T)
    my.url <- paste(my.url, BioGPS.dir,'/', sep = '')
    dirs.data <- retry.connection(target_url = my.url, dirs_only = T, text_connect = F)
    dirs.data <- unlist(strsplit(dirs.data,'\r\n'))

    metadata.file <- grep('metadata.txt', dirs.data, value = T)
    my.url <- paste(my.url,metadata.file, sep = '')

    meta.dat <- retry.connection(target_url = my.url, text_connect = T)

    meta.dat[,1] <- rep(study.dirname,nrow(meta.dat))
    colnames(meta.dat)[1] <- 'Study'

    gc()
    Sys.sleep(1)
    return(meta.dat)
}

make.studies.biogps.df <- function(studies_url, studies_dirs){
    gps_lst <- lapply(X = studies_dirs,
                      FUN =  my.get.study.biogps.metadata,
                      studies.url = studies_url)

    df <- Reduce(f = rbind, x = gps_lst)

    return(df)
}
################################################################################

my.map.BSIDtoSampleName <- function(studies.url, study.dirname){
    ############################################################################
    #Input: (2 character)   URL to all public Atlas studies;
    #                       name of specific study,
    #
    #Output: (single R object) A dataframes which  for the given study,
    #           contains metadata for: Biosample, Study, Donor, Experiment, Run.
    ############################################################################

    #Load the ExRNA-Atlas Result files into memory and extract the
    #BiosampleID sample name pairings
    #Return a dataframe of pairings

    ############################################################################
    #SUB-FUNCTIONS
    ############################################################################
    get.BS.smp.ID <- function(study.dirname, rf.string){
        #Find and retrieve the sample name within an RF file for the study:
        smp.qry.string <- paste(study.dirname,'.*', '/CORE_RESULTS/',sep = '')

        sm <- str_extract(string = rf.string, pattern = smp.qry.string)

        sm <- gsub(pattern = paste(study.dirname,'/',sep = ''), replacement = '', x = sm)

        sm <- gsub(pattern = '/CORE_RESULTS/','',x = sm)

        BS.qry.string <- paste('Biosample ID\tEXR','.*','-BS',sep = '')

        bs <- str_extract(string = rf.string, pattern = BS.qry.string)

        bs <- gsub(pattern = 'Biosample ID\t', replacement = '', x = bs)

        gc()

        return(c(bs,sm))
    }
    ############################################################################
    #Begin my.map.BSIDtoSampleName
    ############################################################################

    my.url <- paste(studies.url, study.dirname, '/', sep = '')
    dirs.data <- retry.connection(target_url = my.url, dirs_only = T, text_connect = F)
    dirs.data <- unlist(strsplit(dirs.data, '\r\n'))

    meta.dir <- grep('metadataFiles', dirs.data, value = T)
    my.url <- paste(my.url, meta.dir, '/', sep = '')
    met.files <- retry.connection(target_url = my.url, dirs_only = T, text_connect = F)
    met.files <- unlist(strsplit(met.files,'\r\n'))

    #Now for each results file (-RF.metadata)
    rf.files <- grep('-RF.metadata.tsv', met.files, value = T)

    mappings <- NULL
    mpp <- NULL

    for(mf in rf.files){

        my.url1 <- paste(my.url, mf, sep = '')

        meta.dat <- NULL

        attempt <- 0

        while( is.null(meta.dat) && attempt <= 10 ) {
            attempt <- attempt + 1
            if(attempt > 1){
                print(paste('retrying ', mf, ': attempt #', attempt, sep = ''))
            }
            try(meta.dat <- getURL(url = my.url1))
            Sys.sleep(1)
        }

        Sys.sleep(1)
        if(is.null(meta.dat) && attempt >= 10){
            mpp <- c("Failed Download", "Failed Dowload")
        }else{
            mpp <- get.BS.smp.ID(study.dirname,meta.dat)
        }
        mpp <- c(mf,mpp,study.dirname)
        mappings <- c(mappings,mpp)
        my.url1 <- NULL


    }

    mappings <- as.data.frame(matrix(mappings,ncol = 4,byrow = T))
    colnames(mappings) <-c("RF File","BS ID", "Sample Name", "Study")
    return( t(mappings))
}

get.combined.bsid.to.smp.name.maps <- function(studies_url, studies_dirs){

    get.maps <- my.map.BSIDtoSampleName

    maps <- lapply(X = studies_dirs,
                   FUN = function(st){print(st); get.maps(studies_url,st)})
    names(maps) <- studies_dirs

    maps <- Reduce(f = cbind, maps)

    maps <- as.data.frame(t(maps))

    return(maps[,c(2,3,4,1)])

}

my.combine.atlas.meta.tables <- function(df1, df2){
    ############################################################################
    #Input: (2 Data Frames)
    #
    #
    #Output: (single R object) A dataframe with values for all rows and columns
    #           of both inputs. NA values are filled in where data is absent
    ############################################################################

    no_data_fill_val <- "NA"

    df1 <- as.matrix( df1 )
    df2 <- as.matrix( df2 )

    r_names_1 <- df1[,1]
    r_names_2 <- df2[,1]

    properties_not_in_df_2 <- setdiff(r_names_1, r_names_2)
    properties_not_in_df_1 <- setdiff(r_names_2, r_names_1)

    n_novel_properties <- length(properties_not_in_df_1)
    n_old_properties   <- length(properties_not_in_df_2)

    add_to_df1 <- (matrix(data = no_data_fill_val, ncol = ncol(df1),
                          nrow = n_novel_properties ))

    add_to_df2 <- (matrix(data = no_data_fill_val, ncol = ncol(df2),
                          nrow = n_old_properties ))


    colnames(add_to_df1) <- colnames(df1)
    colnames(add_to_df2) <- colnames(df2)


    rn_combined_1 <- c(r_names_1, properties_not_in_df_1)
    rn_combined_2 <- c(r_names_2, properties_not_in_df_2)

    df1_expanded <- rbind(df1, add_to_df1)
    rownames(df1_expanded) <- rn_combined_1

    df2_expanded  <- rbind(df2, add_to_df2)
    rownames(df2_expanded) <- rn_combined_2

    df2_expanded <- df2_expanded[rn_combined_1,]

    df_merged <- cbind(df1_expanded, df2_expanded[,-1])
    #df_merged <- cbind(rownames(df_merged), df_merged)
    return( df_merged )

}
################################################################################
my.get.study.metadata <- function(studies_url, study_dirname, meta_types = NULL){
    ############################################################################
    #Input: (3 character)   URL to all public Atlas studies;
    #                       name of specific study;
    #                       vector of meta data types to fetch from the atlas
    #
    #Output: (single R object) A list of dataframes which, for the given study,
    #           contains metadata for: Biosample, Study, Donor, Experiment, Run.
    ############################################################################


    ############################################################################
    #SUB-FUNCTIONS
    ############################################################################
    extract.metadata <- function(meta_url, file_list, merge = T,
                                 file_type_unique_string, sleep = 0.1){

        #Process a single file
        extract.meta <- function(meta_url, metadata_file){
            text_only_meta_types <- c("RF.meta")

            if(file_type_unique_string %in% text_only_meta_types){

                print("Unsupported metadata file.")
                print("Extract metadata by auxiliary process.")
                return(NULL)

            }

            meta_url <- paste(meta_url, metadata_file, sep = '')

            meta_dat <- retry.connection(target_url = meta_url, text_connect = T)

            #cat(class(meta_dat))
            #print(head(rownames(meta_dat)))
            #print(head(colnames(meta_dat)))

            Sys.sleep(sleep)

            return(as.matrix(meta_dat))
        }

        #Get and place metadata into a matrix
        m_files <- grep(file_type_unique_string, file_list, value = T)

        m_dat <- lapply(FUN = extract.meta, X = m_files, meta_url = meta_url)

        gc()

        if(!merge){
            print("Returning unmerged study metadata.")
            return(m_dat)
        }

        ####
        m_dat_mtx <- Reduce(f = function(dat1,dat2){my.combine.atlas.meta.tables(dat1,dat2)},
                            x = m_dat)
        ####
        print(paste("Returning merged study metadata for ", file_type_unique_string, "data",
                    sep = "" ))

        return(m_dat_mtx)
    }
    ############################################################################
    #(END) SUB-FUNCTIONS (END)
    ############################################################################

    #Get the list of metadata files for the study
    study_url <- paste(studies_url, study_dirname, '/', sep = '')

    dirs_data <- getURL(url = study_url, dirlistonly = T)
    dirs_data <- unlist(strsplit(dirs_data,'\r\n'))

    meta_dir <- grep('metadataFiles', dirs_data, value = T)
    meta_url  <- paste(study_url, meta_dir,'/', sep = '')

    dirs_data <- getURL(url = meta_url, dirlistonly = T)
    dirs_data <- unlist(strsplit(dirs_data,'\r\n'))
    meta_files <- grep('metadata.t[xs][tv]', dirs_data, value = T)


    if(is.null(meta_types)){
        cat("\nUsing default metadata types:")
        cat("BS, DO, ST, SU\n")
        meta_types <-  as.list(unlist(strsplit("BS,DO,ST,SU", ",")))
    }

    names(meta_types) <- meta_types
    meta_types <- lapply(FUN = paste, X = meta_types, ".meta", sep = "")

    #Get the metadata
    combined_meta_all <- lapply(X = meta_types,
                                FUN = function(mt){
                                    print(paste(study_dirname, mt))
                                    extract.metadata(meta_url = meta_url,
                                                     file_list = meta_files, merge = T,
                                                     file_type_unique_string = mt)
                                })


    gc()
    Sys.sleep(1)
    return(combined_meta_all)
}

################################################################################
my.get.study.metadata.old <- function(studies_url, study_dirname, meta_types = NULL){
    ############################################################################
    #Input: (3 character)   URL to all public Atlas studies;
    #                       name of specific study;
    #                       vector of meta data types to fetch from the atlas
    #
    #Output: (single R object) A list of dataframes which, for the given study,
    #           contains metadata for: Biosample, Study, Donor, Experiment, Run.
    ############################################################################


    ############################################################################
    #SUB-FUNCTIONS
    ############################################################################
    extract.metadata <- function(meta_url, file_list, merge = T,
                                 file_type_unique_string, sleep = 0.1){

        #Process a single file
        extract.meta <- function(meta_url, metadata_file){
            text_only_meta_types <- c("EX.meta","RF.meta")

            if(file_type_unique_string %in% text_only_meta_types){

                print("Unsupported metadata file.")
                print("Extract metadata by auxiliary process.")
                return(NULL)

            }

            meta_url <- paste(meta_url, metadata_file, sep = '')

            meta_dat <- retry.connection(target_url = meta_url, text_connect = T)

            Sys.sleep(sleep)

            #return(meta_dat)
            return(as.matrix(meta_dat))
        }

        #Get and place metadata into a data frame
        m_files <- grep(file_type_unique_string, file_list, value = T)

        m_dat <- lapply(FUN = extract.meta, X = m_files, meta_url = meta_url)

        gc()

        if(!merge){
            print("Returning unmerged study metadata.")
            return(m_dat)
        }


        #m_dat_df <- Reduce(function(x,y){merge.data.frame(x, y,by = 1, all = T)}, m_dat)
        m_dat_mtx <- Reduce(f = function(x,y){cbind(x,y[,-1])},x = m_dat)

        print(paste("Returning merged study metadata for ", file_type_unique_string, "data",
                    sep = "" ))

        #return(m_dat_df)
        return(m_dat_mtx)
    }
    ############################################################################
    #(END) SUB-FUNCTIONS (END)
    ############################################################################

    #Get the list of metadata files for the study
    study_url <- paste(studies_url, study_dirname, '/', sep = '')

    dirs_data <- getURL(url = study_url, dirlistonly = T)
    dirs_data <- unlist(strsplit(dirs_data,'\r\n'))

    meta_dir <- grep('metadataFiles', dirs_data, value = T)
    meta_url  <- paste(study_url, meta_dir,'/', sep = '')

    dirs_data <- getURL(url = meta_url, dirlistonly = T)
    dirs_data <- unlist(strsplit(dirs_data,'\r\n'))
    meta_files <- grep('metadata.t[xs][tv]', dirs_data, value = T)


    if(is.null(meta_types)){
        cat("\nUsing default metadata types:")
        cat("BS, DO, ST, SU\n")
        meta_types <-  as.list(unlist(strsplit("BS,DO,ST,SU", ",")))
    }

    names(meta_types) <- meta_types
    meta_types <- lapply(FUN = paste, X = meta_types, ".meta", sep = "")

    #Get the metadata
    combined_meta_all <- lapply(X = meta_types,
                                FUN = function(mt){
                                    print(paste(study_dirname, mt))
                                    extract.metadata(meta_url = meta_url,
                                                     file_list = meta_files, merge = T,
                                                     file_type_unique_string = mt)
                                })


    gc()
    Sys.sleep(1)
    return(combined_meta_all)
}

################################################################################
write.study.meta.obj <-function(st){

    st_all_meta <- my.get.study.metadata(studies_url = studies_url,
                                         study_dirname = st,
                                         meta_types)

    save_fldr <- './interim/exrna_atlas_meta/'

    save_name <- paste(st,'_all_meta.RDS')

    save_name <- gsub(" ","",save_name)

    save_path <- paste(save_fldr,save_name,sep = '')

    saveRDS(st_all_meta, save_path)
}

combine.studies.meta.obj <- function(studies_dirs){

    objs_path <- '../get_atlas_data_in_R/interim/exrna_atlas_meta/'

    dir.create(objs_path,showWarnings = F)

    meta_obj  <- lapply(X = studies_dirs,
                    FUN = function(st){
                        file_pth <- paste(objs_path,st,'_all_meta.RDS', sep = '')
                        readRDS(file = file_pth)
                })
    names(meta_obj) <- studies_dirs

    return(meta_obj)

}

create.aggregate.study.metadata.obj <- function(studies_url, studies_dirs, meta_types){

    sapply(X = studies_dirs, FUN =  write.study.meta.obj)

    combine.studies.meta.obj(studies_dirs)


}

################################################################################
get.atlas.metadata <- function(studies_url, studies_dirs, meta_types = FALSE, get_map = F, update){
    ############################################################################
    #Input: (3 character)   URL to all public Atlas studies;
    #                       name of desired studies;
    #                       name of desired meta data doc types (BS, DO, SU, ST)
    #
    #Output: (single R object) A list of dataframes which, for the doc type,
    #           contains metadata for: all sudies in studies_dirs
    ############################################################################

    cat("\n", rep("=",60), "\n", sep = "")
    cat("\n", "Acquiring metadata for:", "\n")
    print(studies_dirs)
    cat("\n", rep("=",60), "\n", sep = "")

    #meta_list <<- vector(mode = 'list', length = length(studies_dirs))
    #names(meta_list) <- studies_dirs

    meta_list <-  create.aggregate.study.metadata.obj(studies_url, studies_dirs, meta_types)

    cat("\n", rep("=",60), "\n", sep = "")
    cat("\n", "Generating Biosample_ID to sample_name map for:", "\n")
    print(studies_dirs)
    cat("\n", rep("=",60), "\n", sep = "")

    if(get_map){
        meta_list[["bsid_sample_map"]] <- get.combined.bsid.to.smp.name.maps(studies_url,
                                                                             studies_dirs)

    }

    meta_list[["Access_Date"]] <- paste( "This record was created on:",
                                                      Sys.time(),sep = " " )
    cat("Timestamp Added.")

    gc()

    return(meta_list)
}






