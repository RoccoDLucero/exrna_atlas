###################################################################################################
##title:                 get_atlas_data_in_R_generic_functions.R
##associated_project:    exrna_atlas
##associated_analysis:   get_atlas_data_in_R
##started_date:          April 2017
##started_by:            Rocco Lucero
##last_updated_date:     June 22, 2017
##last_updated_by:       Rocco Lucero
###################################################################################################
##description:   generic function to assist in retrieval and processing of exrna atlas contents
##
##inputs:                <NA>
##
##outputs:               <NA>
##
##dependencies:          <NA>
###################################################################################################
##Table of Contents:
##--my.get.print.studies
##--retry.connection
##--grp.by.list.names
##--remove.columns
##--my.combine.atlas.rpm.tables
##--sep.df.by.colname.pttn
##--exclude.if.duplicated
##--rds_save_output
################################################################################
#Functions to interface with the ExRNA Atlas Data on the BRL FTP site
################################################################################
memory.limit(size = 4095)


################################################################################
my.get.print.studies <- function(studies_url, print_only = T){

    dirs_data <- getURL(url = studies_url, dirlistonly = T)
    dirs_data <- unlist(strsplit(dirs_data,'\r\n'))

    cat(rep("=",60), "\r\n", sep = "")
    cat(rep("=",60), "\r\n", sep = "")
    cat( length(dirs_data), " studies found : \r\n")
    print(dirs_data)
    cat(rep("=",60), "\r\n", sep = "")
    cat(rep("=",60), "\r\n", sep = "")

    if(!print_only){return(dirs_data)}
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

    #Suggestions: Needs Error Handling, make this generic

    no_data_fill_val <- 0

    df1 <- as.matrix( df1 )
    df2 <- as.matrix( df2 )

    r_names_1 <- rownames(df1)
    r_names_2 <- rownames(df2)

    rna_not_in_df_2 <- setdiff(r_names_1, r_names_2)
    rna_not_in_df_1 <- setdiff(r_names_2, r_names_1)

    n_novel_ids <- length(rna_not_in_df_1)
    n_old_ids   <- length(rna_not_in_df_2)

    #df1 <- ( remove.columns(df1, 1))
    #df2 <- ( remove.columns(df2, 1))

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
rds_save_output <- function(fun, args, save_path, save_as_name){

    output_name <- paste(save_path, save_as_name, sep = "")

    obj <- do.call(what = fun, args = args)

    saveRDS(object = obj, file = output_name)

    return(NULL)

}
