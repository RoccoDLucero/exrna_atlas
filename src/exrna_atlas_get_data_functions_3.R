###################################################################################################
##title:                 exrna_atlas_get_data_functions
##associated_project:    exrna_atlas
##associated_analysis:   XXXXX
##started_date:          October 20, 2016
##started_by:            Rocco Lucero
##last_updated_date:     June 8, 2017
##last_updated_by:       Rocco Lucero
###################################################################################################
##description:           <What task does this script purport to accomplish??>
##
##inputs:                <What goes in>
##
##outputs:               <What are the outputs>
##
##dependencies:          <Name any dependencies not publicly available>
###################################################################################################
#These functions should be used to retrieve public ExRNA Atlas and process it
#Into R objects that can then be used to efficently analyse all/any of the studies
#the ExRNA Atlas.
#
#Table of Contents:
#.....GENERIC FUNCTIONS
#--my.print.studies
#--retry.connection
#--grp.by.list.names
#--remove.columns
#--my.combine.atlas.rpm.tables
#--sep.df.by.colname.pttn
#
#.....GET AND PROCESS READ COUNTS
#--my.get.study.post.processed.results
#--my.combine.study.reads
#--my.get.atlas.readcounts (MAIN)
#
#.....GET AND PROCESS METADATA
#--my.get.study.bioGPS.metadata
#--my.map.BSIDtoSampleName
#--my.get.study.metadata
#       --get.BS.smp.ID
#       --extract.metadata
#               --extract.meta
#--my.merge.metadata (MAIN)
#
#.....PROCESS RESULTS FROM WORKBENCH EXCERPT JOBS
#--my.load.excerpt.rdata
#
#These functions work with scripts including:
#--generate.dfs.from.ftp.R
#--metadataprocessing.R
#(You're welcome, future Rocco...)
################################################################################
################################################################################
####################
required_packages <- c("RCurl", "stringr")

lapply(required_packages, library, character.only = T)


set.seed(20161017)
####################




################################################################################
#Functions to interface with the ExRNA Atlas Data on the BRL FTP site
################################################################################
my.get.print.studies <- function(studies_url){

        dirs_data <- getURL(url = studies_url, dirlistonly = T)
        dirs_data <- unlist(strsplit(dirs_data,'\r\n'))

        cat( length(dirs_data), " studies found : ")
        cat(rep("=",60), "\r\n", sep = "")
        print( dirs_data )
        cat(rep("=",60), "\r\n", sep = "")

        return(dirs_data)
}

################################################################################
retry.connection <- function(target_url, text_connect = T, dirs_only = F){

        dat <- NULL
        attempt <- 0

        while( is.null(dat) && attempt <= 10 ) {

            attempt <- attempt + 1

            if(attempt > 1){

                cat("\r\nRetrying connection:","\r\n","\t",
                    "attempt #", attempt, "\r\n", sep = '')

                cat("\r\n", target_url)
            }

            try(dat <- getURL(url = target_url, dirlistonly = dirs_only))

            if(text_connect){

                try(dat_con <- textConnection(dat, "r"))

                try(dat <- read.table(dat_con, sep = '\t', header = T, fill = F,
                                      stringsAsFactors = F, quote = "", comment.char = ""))

                try(close.connection(dat_con))
            }

            Sys.sleep(.1)
        }

        return(dat)
}

################################################################################
grp.by.list.names <- function(lst){
        labels <- unique(names(lst))

        grouping <- sapply(X   = labels, lst = lst,
                           FUN = function(lab,lst){which(names(lst) == lab)})

        return(grouping)
}

################################################################################
remove.columns <- function(df, remove_cols){ df[, setdiff(seq(ncol(df)), remove_cols)]}

################################################################################
my.combine.atlas.rpm.tables <- function(df1, df2){
        #Inputs: (2 Data Frames)

        #Outputs: (1 Data Frame)

        #Suggestions: Needs Error Handling

        no_data_fill_val <- "0"

        df1 <- as.matrix( df1 )
        df2 <- as.matrix( df2 )

        r_names_1 <- df1[,1]
        r_names_2 <- df2[,1]

        rna_not_in_df_2 <- setdiff(r_names_1, r_names_2)
        rna_not_in_df_1 <- setdiff(r_names_2, r_names_1)
        n_novel_ids <- length(rna_not_in_df_1)
        n_old_ids   <- length(rna_not_in_df_2)

        df1 <- ( remove.columns(df1, 1))
        df2 <- ( remove.columns(df2, 1))

        n_old_samples <- ncol(df1)
        n_new_samples <- ncol(df2)

        cat("\nAdding ", n_novel_ids, "RNA species to ", length(r_names_1), "indexed RNA." )
        cat("\nAdding ", n_new_samples, "new samples to ", n_old_samples, "existing records.\n")


        add_to_df1 <- (matrix(data = no_data_fill_val, ncol = ncol(df1),
                              nrow = n_novel_ids ))

        add_to_df2 <- (matrix(data = no_data_fill_val, ncol = ncol(df2),
                              nrow = n_old_ids ))


        colnames(add_to_df1) <- colnames(df1)
        colnames(add_to_df2) <- colnames(df2)


        rn_combined_1 <- c(r_names_1, rna_not_in_df_1)
        rn_combined_2 <- c(r_names_2, rna_not_in_df_2)

        df1_expanded <- rbind(df1, add_to_df1)
        rownames(df1_expanded) <- rn_combined_1

        df2_expanded  <- rbind(df2, add_to_df2)
        rownames(df2_expanded) <- rn_combined_2

        df2_expanded <- df2_expanded[rn_combined_1,]

        df_merged <- cbind(df1_expanded, df2_expanded)
        df_merged <- cbind(rownames(df_merged), df_merged)

        return( df_merged )

}

################################################################################
sep.df.by.colname.pttn <- function(df, patterns = NULL){

        gencode_types <- list("Y_RNA", "snoRNA", "snRNA", "lincRNA", "miRNA",
                              "protein_coding", "processed_transcript" )

        if(is.null(patterns)){patterns <- gencode_types}

            make_new_df <-  function(rna_type){
                                df_cols <- grep(rna_type, colnames(df))
                                return(df[,df_cols])

                            }

        new_dfs <- lapply(X = patterns, FUN = make_new_df)

        names(new_dfs) <- patterns

        return(new_dfs)

}

################################################################################
exclude.if.duplicated <- function(x){
        if(is.vector(x)){

            return(x[!(duplicated(x) | duplicated(x, fromLast = TRUE))])

        }

        if(is.data.frame(x)){

            rn <- rownames(x)

            return( x[!(duplicated(rn) | duplicated(rn, fromLast = TRUE)), ] )

        }

}

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

        Sys.sleep(.1)

        return( reads_dat )

}

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
        df_cmbnd <- remove.columns(df_cmbnd, remove_cols = 1)
        df_cmbnd <- t(df_cmbnd)
        rn <- rownames(df_cmbnd)

        df_cmbnd <- as.data.frame(df_cmbnd)
        df_cmbnd <- apply(X = df_cmbnd, MARGIN = 2, FUN = as.numeric)
        rownames(df_cmbnd) <- rn

        gc()

        return( df_cmbnd )

}

################################################################################
#THIS IS THE MAIN FUNCTION USED TO GET READCOUND DATA
################################################################################
my.get.atlas.readcounts <- function(rna_types, dirs_studies, studies_url,
                                    per_million = T, check_names = T){

        cat("\n", rep("=", 60), "\n", sep = "")
        cat("\n", "Acquiring read count for:", "\n")
        cat(dirs_studies)
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


