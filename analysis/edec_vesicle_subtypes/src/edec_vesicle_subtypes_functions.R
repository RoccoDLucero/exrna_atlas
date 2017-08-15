#edec_vesicle_subtypes_functions.R



#combine.atlas.rpm.tables
#condense.post.proc.results
#make_excerpt_results_list

################################################################################
combine.atlas.rpm.tables <- function(df1, df2){
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
    
   
    df2_expanded <- df2_expanded[rn_combined_1,,drop = F]

    #df1 <<- dimnames(df1_expanded)
    #df2 <<- dimnames(df2_expanded)
        
    df_merged <- cbind(df1_expanded, df2_expanded)
    
    return( df_merged )
    
}

################################################################################
condense.post.proc.results <- function(results_rdata_path, 
                               per_million = T, include_exogenous = F,
                               exclude_tables, exclude_gencode_types){
    
    ## REQUIRES FUNCTIONS: combine.atlas.rpm.tables, make_excerpt_results_list
    
    ## FIRST PUT ALL THE RESULTS FILES TO BE PROCESSED IN A SINGLE FOLDER
    ## PROVIDE THE PATH TO THIS FOLDER AS THE FIRST ARGUMENT
    ## OUTPUTS A SINGLE TABLE OF READCOUNTS FOR ALL RNA AND ALL STUDIES 
    
    path <- results_rdata_path
    
    if(per_million){ptt <- 'PerMillion.RData' }else{ppt <- 'ReadCounts.RData'}
    
    files <- list.files(path, pattern = ptt, full.names = T)
    
    names <- list.files(path, pattern = ptt, full.names = F)
    
    res_list <- mapply(FUN = make_excerpt_results_list, excerpt_rdata_files = files,
                       names_for_output = names, exogenous = include_exogenous)
    names(res_list) <- names
    
    merged_runs_list <- lapply(X = res_list, function(excerpt_run){Reduce(f = rbind, excerpt_run) })
    
    all_merged <- Reduce(f = combine.atlas.rpm.tables, merged_runs_list)
}

################################################################################
logistic = function(x){
    a = 1/max(x)
    1 - exp(1)^(-a*x)
}

################################################################################
make_excerpt_results_list <- function(excerpt_rdata_files, names_for_output, exogenous = F){
    
    load_to_list <- function(rdata, exo = exogenous){
        
        load(rdata)
        
        if(!exo){rm(list = ls(pattern = 'exogenous'))}
        
        vars <- ls(pattern = 'exprs.*rpm')
        
        print(vars)
        
        mget(vars)
        
    }
    
    fls <- excerpt_rdata_files
    
    res <- lapply(X = fls, FUN = load_to_list, exo = exogenous)
    
    names(res) <- names_for_output
    
    return(res)
    
}

################################################################################
quantile_normalisation <- function(df){
    df_rank <- apply(df,2,rank,ties.method="min")
    df_sorted <- data.frame(apply(df, 2, sort))
    df_mean <- apply(df_sorted, 1, mean)
    
    index_to_mean <- function(my_index, my_mean){
        return(my_mean[my_index])
    }
    
    df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
    rownames(df_final) <- rownames(df)
    return(df_final)
}


