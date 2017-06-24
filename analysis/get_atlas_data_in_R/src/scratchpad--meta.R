

my.combine.atlas.meta.tables <- function(df1, df2){
    #Inputs: (2 Data Frames)

    #Outputs: (1 Data Frame)

    #Suggestions: Needs Error Handling, make this generic

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


get_meta_args <- list(studies_url = studies_url, studies_dirs = dirs_studies[1],
                      meta_types = meta_types, get_map = F)

rds_save_output(fun = get.atlas.metadata, args = get_meta_args,
                save_path = save_path,
                save_as_name = 'exrna_atlas_metadata_and_map_to_samples.RDS')

#b <- readRDS(file = '../get_atlas_data_in_R/interim/exrna_atlas_meta/KJENS1-aSAH_Project-2016-10-13_all_meta.RDS')
b <- readRDS(file = '../get_atlas_data_in_R/interim/exrna_atlas_metadata_and_map_to_samples.RDS')
sapply(X = b, FUN = length)
sapply(X = b, FUN = function(c){sapply(X = c, FUN = dim)})
dirs_studies


