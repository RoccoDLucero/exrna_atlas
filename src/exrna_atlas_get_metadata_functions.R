################################################################################
#FUNCTIONS FOR RETRIEVING AND PROCESSING ATLAS METADATA
################################################################################

my.get.study.bioGPS.metadata <- function(studies.url, study.dirname){
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
    #meta.dat <- getURL(url = my.url)
    #meta.con <- textConnection(meta.dat,"r")
    #meta.dat <- read.table(meta.con, sep = '\t', header = T, fill = F,
    #                       stringsAsFactors = F,quote = "")
    meta.dat <- retry.connection(target_url = my.url, text_connect = T)

    meta.dat[,1] <- rep(study.dirname,nrow(meta.dat))
    colnames(meta.dat)[1] <- 'Study'

    #close.connection(meta.con)

    gc()
    Sys.sleep(1)
    return(meta.dat)
}


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
    #mf <- rf.files[2]
    mappings <- NULL
    mpp <- NULL

    #retry.connection(target_url = my.url, dirs_only = F, text_connect = F)
    for(mf in rf.files){
        print(mf)

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
            text_only_meta_types <- c("EX.meta","RF.meta")

            if(file_type_unique_string %in% text_only_meta_types){

                print("Unsupported metadata file.")
                print("Extract metadata by auxiliary process.")
                return(NULL)

            }

            meta_url <- paste(meta_url, metadata_file, sep = '')

            meta_dat <- retry.connection(target_url = meta_url, text_connect = T)

            Sys.sleep(sleep)

            return(meta_dat)
        }

        #Get and place metadata into a data frame
        m_files <- grep(file_type_unique_string, file_list, value = T)

        m_dat <- lapply(FUN = extract.meta, X = m_files, meta_url = meta_url)

        gc()

        if(!merge){
            print("Returning unmerged study metadata.")
            return(m_dat)
        }

        m_dat_df <- Reduce(function(x,y){merge.data.frame(x, y,by = 1, all = T)}, m_dat)
        #m_dat_df <- Reduce(merge, m_dat)

        print(paste("Returning merged study metadata for ", file_type_unique_string, "data",
                    sep = "" ))

        return(m_dat_df)
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
my.merge.metadata <- function(studies_url, studies_dirs, meta_types = FALSE){
    ############################################################################
    #Input: (3 character)   URL to all public Atlas studies;
    #                       name of desired studies;
    #                       name of desired meta data doc types (BS, DO, SU, ST)
    #
    #Output: (single R object) A list of dataframes which, for the doc type,
    #           contains metadata for: all sudies in studies_dirs
    ############################################################################


    ############################################################################
    #SUB-FUNCTIONS
    ############################################################################

    cat("\n", rep("=",60), "\n", sep = "")
    cat("\n", "Acquiring metadata for:", "\n")
    print(studies_dirs)
    cat("\n", rep("=",60), "\n", sep = "")

    all_meta_by_study <- lapply(studies_dirs,
                                function(st){
                                    my.get.study.metadata(studies_url = studies_url,
                                                          study_dirname = st,
                                                          meta_types)
                                })

    all_meta_by_doc <-  unlist(all_meta_by_study, recursive = F, use.names = T)

    md_doc_groups <- grp.by.list.names( all_meta_by_doc )

    all_meta_by_doc_group <- apply(md_doc_groups, 2,
                                   function(grp){all_meta_by_doc[grp]})


    merged_meta_by_doc_group <- lapply(all_meta_by_doc_group,
                                       function(doc_lst){
                                           Reduce(function(x, y){
                                               merge.data.frame(x, y, by = 1, all = T)},
                                               doc_lst)
                                       })

    merged_meta_by_doc_group <- lapply( merged_meta_by_doc_group,
                                        function(df){
                                            df <- t(df)

                                            df[1,] <- gsub("^X\\.+", "", make.names(df[1,]))

                                            return(df)
                                        })

    merged_meta_by_doc_group["Access Date"] <- paste( "This record was created on:",
                                                      Sys.time(),sep = " " )
    cat("Timestamp Added.")

    gc()

    return( merged_meta_by_doc_group )
}
