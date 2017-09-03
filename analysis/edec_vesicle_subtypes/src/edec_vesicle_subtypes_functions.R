#edec_vesicle_subtypes_functions.R



#combine.rpm.tables
#condense.post.proc.results
#make_excerpt_results_list

################################################################################

## 3' ADAPTER FOR GINGERAS/ENCODE SUBCELLULAR COMPARTMENT DATA
## AAAAAAAAAAAAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGA

combine.rpm.tables <- function(df1, df2){
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
    
    ## REQUIRES FUNCTIONS: combine.rpm.tables, make_excerpt_results_list
    
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
    
    all_merged <- Reduce(f = combine.rpm.tables, merged_runs_list)
    
    print("Removing Duplicate records for: ")
    print(colnames(all_merged)[duplicated(colnames(all_merged))])
    
    return( all_merged[, unique(colnames(all_merged))])
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

subset.small.rna <- function(df, rna_in_colnames = T){
    small_rna_keys <- paste("Gly|Ala|Leu|Met|Phe|Trp|Lys|Gln|Glu|Ser",
                            "Pro|Val|Ile|Cys|Tyr|His|Arg|Asn|Asp|Thr|Mt_tRNA",
                            "hsa-let|hsa-mi[Rr]|hsa_piR|snoRNA|snRNA|^Y_RNA",
                            #"hsa_circ|circ",
                            "3prime_overlapping|vaultRNA",
                            "miRNA|sRNA|scRNA|VTRNA",
                            sep = "|")
    if(rna_in_colnames){
        dim <- 2
        return( df[, grep(small_rna_keys, dimnames(df)[[dim]]) ])
    
    }else{
        dim <- 1
        return( df[ grep(small_rna_keys, dimnames(df)[[dim]]) ,])
    
    }
    
}

get.study.specific.probes <- function(results_rdata_path, 
                                      per_million = T, include_exogenous = F,
                                      exclude_tables, exclude_gencode_types){
    
    ## REQUIRES FUNCTIONS: combine.tables, make_excerpt_results_list
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
    
    probes_lists <- lapply(X =  merged_runs_list,
                           
                           FUN = function(x){dn <- dimnames(x)
                           names(dn) <- c("row","col")
                           return(dn)
                           })
    
    names(probes_lists) <- names
    
    return(probes_lists)
}

get.sra.expt.info <- function(study_identifiers, query_col = 'study_accession'){
    
    study_name <- paste("'", study_identifiers[1], "'", sep = "")
    
    exp_rs_ID <- dbGetQuery(conn = sra_con,
                            statement = paste( "select * from experiment",
                                               "where", query_col, "=", study_name,
                                               sep=" "))
    
    #GET ADDTIONAL INFROMATION ABOUT THE SAMPLES
    smp_rs_ID <- dbGetQuery(conn = sra_con,
                            statement = paste("select * from sample inner join experiment on",
                                              "sample.sample_accession=experiment.sample_accession",
                                              "where", query_col, "=", study_name, sep = " "))
    
    res <- list(sra_experiment = exp_rs_ID,
                sra_expt_samples = smp_rs_ID,
                geo_series_mtx = NULL)
    
    res$geo_series_mtx <- getGEO(GEO = study_identifiers[2],destdir = "./input/geo_dat", GSEMatrix = F)
    
    return(res)
    
}

move_files_to_subfolder <- function(files_vec, from_path){
    
    sub <- deparse(substitute(files_vec))
    
    sub <- paste(from_path, sub, sep = '')
    
    if(!dir.exists(sub)){dir.create(sub)}
    
    files_avail <- list.files(from_path)
    print(files_avail)
    
    found_files <- intersect(files_avail, files_vec)
    print(found_files)
    
    not_found <- setdiff(files_vec, files_avail)
    
    to_files <- paste(sub, '/', found_files, sep = '')
    
    from_files <- paste(from_path, found_files, sep = '')
    
    file.rename(from_files, to_files)
    
}

multigrep <- function(patterns, target, ...){
    
    xx <- sapply(X = patterns, FUN = grep, x = target)
    Reduce(f = intersect, x = xx)
}

mk.col <-function(factor_vec, color_fun = rainbow){
    
    color.function <- color_fun
    
    colors <- as.factor(factor_vec)
    
    levels(colors) <- color.function(length(levels(colors)))
    
    return(as.character(colors))
    
}


process.dimnames <- function(excerpt_counts_table){
    
    dn <- dimnames(excerpt_counts_table)
    
    ## PROCESS COLNAMES (SAMPLES)
    dn[[2]] <- gsub("sample_|_fastq","", dn[[2]])
    
    ## PROCESS ROWNAMES (RNA)
    fx <- grep("hsa_piR", dn[[1]])
    
    dn[[1]][fx] <- gsub(pattern = "\\|.*$", replacement = "", x = dn[[1]][fx])
    
    return(dn)
    
}

process.names.list <- function(rna_names_list){
    
    lst <- rna_names_list
    
    ## PROCESS ROWNAMES (RNA)
    names_out <- lapply(X = lst,
                        FUN = function(x){
            
                            fx <- grep("hsa_piR", x)
                
                            x[fx] <- gsub(pattern = "\\|.*$",
                                          replacement = "",
                                          x = x[fx])
                            
                            return(x)
                            }
                        
                        )
    
    return(names_out)
    
}


################################################################################
#Identify exRNA detected in a given propotion of samples at a given rpm value
# so that we can reduce the size of the data
filterByPctDtect <- function(x, min_detect_proportion, min_rpm){
    detect_prop <- length(x[x >= min_rpm]) / length(x)
    
    if(detect_prop >= min_detect_proportion){return(TRUE)
    }else{ return(FALSE)}
    
}
#We want our chosen miRNA to show variability across conditions so we may be able to
# remove miRNA that do not vary enough
filterByVariance <- function(x, min_variance){
    v <- var(x)
    if(v < min_variance){return(TRUE)}else{return(FALSE)}
    
}

filterBySD <- function(x, min_sd){
    s <- sd(x)
    if(s < min_sd){return(TRUE)}else{return(FALSE)}
    
}

filterByIQR <- function(x, min_iqr){
    i <- IQR(x)
    if(i < min_iqr){return(TRUE)}else{return(FALSE)}
    
}


removeRareRNA <- function(df, min_detect_proportion = .05, min_rpm = 1){
    r <- apply(df, 2, filterByPctDtect, min_detect_proportion, min_rpm)
    
    return(df[,r])
    
}

#Functions for various data transformations that will be neeed to ensure compatibility
# With EDec
txfmQN  <- function(df){preprocessCore::normalize.quantiles(df)} #columns will share a distribution
txfmNPN <- function(df){huge::huge.npn(x = df)}
txfmLinear <- function(x){(x-min(x))/(max(x)-min(x))}
txfmLogistic1 <- function(x, a = 1/max(x)){ 1 / (1 + exp(1)^(-(a * x )) ) }
txfmLogistic2 <- function(x, a = 1/max(x)){ 1.5*( 1 - exp(1)^(-(a * x)) ) }
txfmLogistic3 <- function(x, a = 1/max(x)){ 1 / (1 + exp(1)^(-(a * (x-mean(x)) ) )) }
txfmLogit <- function(x){log( x / (1 - x) )}

factor_to_numeric <- function(fac){as.numeric(as.character(fac))}


################################################################################
get.normalized.readcounts.obj <- function(df, na_rm = T, vars_dim = 1){
    
    df <- as.matrix(df, row.names = rownames(df), col.names = colnames(df))
    
    df_qn <- get.qn(df)
    rownames(df_qn) <- rownames(df)
    colnames(df_qn) <- colnames(df)
    
    df_npn <- huge.npn(df_qn)
    
    df_qn_lgstc <- apply(X = df_qn, MARGIN = vars_dim, FUN = txfmLogistic1)
    
    df_npn_lgstc <- apply(X = df_npn, MARGIN = vars_dim, FUN = txfmLogistic1)
    
    if(na_rm){
        df_qn_lgstc <- df_qn_lgstc[ , complete.cases(df_qn_lgstc)]
        
        df_npn_lgstc <- df_npn_lgstc[ , complete.cases(df_npn_lgstc)]
        
    }
    
    normalized_lst <- list(qn = df_qn,
                           qn.npn = df_npn,
                           qn.lgstc = t(df_qn_lgstc),
                           qn.npn.lgstc =  t(df_npn_lgstc),
                           original = df)
    
    return(normalized_lst)
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


subset.on.shared.probes <- function(ref_df, exp_df, samples_in_cols = T){
    
    if(samples_in_cols){probe_dim <- 1}else{probe_dim <- 2}
    
    avail_probes <- intersect(dimnames(ref_df)[[probe_dim]],
                              dimnames(exp_df)[[probe_dim]])
    
    cat(length(avail_probes), "probes are shared." )
    
    if(probe_dim == 2 ){
        ref_df <- ref_df[, avail_probes]
        exp_df <- exp_df[, avail_probes]
    }else{
        ref_df <- ref_df[avail_probes,]
        exp_df <- exp_df[avail_probes,]
        
    }
    
    return(list(ref = ref_df, exp = exp_df))
}

make.colors.vector <- function(factor_vec, color_fun = rainbow){
    
    color.function <- color_fun
    
    colors <- as.factor(factor_vec)
    
    levels(colors) <- color.function(length(colors))
    
    return(as.character(colors))
    
}


my.run.edec.stage.0 <- function (reference_meth, reference_classes,
                                 max_p_value, num_markers, version = "one.vs.rest"){
    num_comps = NULL
    num_hyper_markers_per_comp = NULL
    t_test_results = NULL
    if (version == "one.vs.rest") {
        num_comps <- length(unique(reference_classes))
        t_test_results <- EDec:::perform_t_tests_all_classes_one_vs_rest(reference_meth,
                                                                         reference_classes)
        num_hyper_markers_per_comp <- ceiling((num_markers/num_comps)/2)
    }
    else if (version == "each.pair") {
        t_test_results <- EDec:::perform_t_tests_all_classes_each_pair(reference_meth,
                                                                       reference_classes)
        num_comps <- ncol(t_test_results$p.values)
        num_hyper_markers_per_comp <- ceiling((num_markers/num_comps)/2)
    }
    else {
        stop("Error: type has to be either 'one.vs.rest' or 'each.pair'")
    }
    
    markers <- NULL
    for (i in 1:num_comps) {
        significant_loci <- names(which(t_test_results$p.values[,i] <= max_p_value))
        
        significant_loci <- setdiff(significant_loci, markers)
        
        hyper_loci <- names(utils::tail(sort(t_test_results$diff.means[significant_loci, i]),
                                        num_hyper_markers_per_comp))
        
        hypo_loci <- names(utils::head(sort(t_test_results$diff.means[significant_loci, i]),
                                       num_hyper_markers_per_comp))
        
        markers <- c(markers, hyper_loci, hypo_loci)
    }
    
    num_extra_markers <- length(markers) - num_markers
    
    if (num_extra_markers > 0) {
        min_p_values <- apply(t_test_results$p.values[markers,], 1, min)
        
        markers_to_drop <- utils::tail(sort(min_p_values, index.return = TRUE)$ix,
                                       num_extra_markers)
        
        markers <- markers[-markers_to_drop]
    }
    return(list(probes = markers, t_test = t_test_results))
}


perform.probe.selection <- function(ref_probes, ref_classes, max_p_val,
                                    n_markers, method = list("one.vs.rest", "each.pair" )){
    
    info_probes <- lapply(X = method,
                          FUN = function(mth){
                              print(mth)
                              my.run.edec.stage.0(version = mth, reference_meth = ref_probes,
                                                  reference_classes = ref_classes,
                                                  max_p_value = max_p_val, num_markers = n_markers)}
                          
    )
    
    names(info_probes) <- unlist(method)
    
    if(length(info_probes) > 1){info_probes[['max_set']] <- Reduce(f = union, x = info_probes)}
    
    info_probes[['no_selection']] <- rownames(ref_probes)
    
    return(info_probes)
    
}


check_stage1_results <- function(stage1_result, probeset = 'all', reference_profiles,
                                 sample_meta_colors = F,
                                 ref_colors = F, show_corrs = T, show_props = T){
    
    n_ct <- ncol(stage1_result$methylation)
    use_all_probes <-  (length(probeset) == 1  && probeset == "all")
    if(use_all_probes){probeset <- seq(nrow(stage1_result$methylation))}
    
    cors_deconv_refs = cor(reference_profiles[probeset, ], stage1_result$methylation[probeset,] )
    
    # Check what references had the highest correlation with each
    # of the estimated methylation profiles
    best_cors = rbind(apply(cors_deconv_refs,2,which.max),
                      apply(cors_deconv_refs,2,max))
    
    best_cor_labels = matrix("",nrow=nrow(cors_deconv_refs),
                             ncol=ncol(cors_deconv_refs))
    for (i in seq(n_ct)){
        best_cor_labels[best_cors[1,i],i] = as.character(round(best_cors[2,i],2))
    }
    
    if(show_corrs){
        # Plot correlation matrix
        #main = "Inferred vs. reference\nexpression profile correlations"
        main = 'title'
        gplots::heatmap.2(cors_deconv_refs,
                          trace="none", col = heat_colors,
                          breaks = heat_breaks, margins = heat_margins,
                          RowSideColors = ref_colors,
                          main = main, key.xlab = "Expression correlation",
                          cellnote = best_cor_labels,
                          notecol="black")
    }
    
    if(show_props){
        colnames(stage1_result$proportions) <- colnames(stage1_result$methylation)
        
        col_grad = colorRampPalette(c("white","steelblue"))
        
        #main = "Per sample Composition\nof inferred profiles"
        gplots::heatmap.2( stage1_result$proportions,
                           trace="none",
                           col = col_grad(10),
                           breaks=seq(0,1,0.1),
                           margins=c(10,4),
                           labRow = FALSE ,
                           RowSideColors = sample_meta_colors)
    }
}

############################################################################################################

get.rna.species.by.type <- function(atlas_all_data_list_obj_path,
                                       rna_names_dim = 2){
    
    r_dim <-rna_names_dim
    
    atls <- readRDS(atlas_all_data_list_obj_path)
    
    r_names <- lapply(FUN = dimnames, atls)
    
    r_names <- lapply(X = r_names, FUN = function(x){return(x[[r_dim]])})
    
    names(r_names) <- names(atls)
    
    return(r_names)
    
}

################################################################################
################################################################################
################################################################################
#functions from K_cell_number.R

pMatrix.min <- function(A, B) {
    # finds the permutation P of A such that ||PA - B|| is minimum
    # in Frobenius norm
    # Uses the linear-sum assignment problem (LSAP) solver
    # in the "clue" package
    # Returns P%*%A and the permutation vector `pvec' such that
    # A[pvec, ] is the permutation of A closest to B
    n <- nrow(A)
    D <- matrix(NA, n, n)
    for (i in 1:n) {
        for (j in 1:n) {
            #           D[j, i] <- sqrt(sum((B[j, ] - A[i, ])^2))
            D[j, i] <- (sum((B[j, ] - A[i, ])^2))  # this is better
        } }
    vec <- c(solve_LSAP(D))
    list(A=A[vec,], pvec=vec)
}


meth.cor.min <- function(E1, E2,i) {
    d = diag(1, i, i)
    c = cor(E1,E2)
    f= pMatrix.min(c,d)
    min(diag(f$A)) }


prop.cor.min <- function(E1, E2,i) {
    d = diag(1, i, i)
    com = intersect(row.names(E1),row.names(E2))
    c = cor(E1[com,],E2[com,])
    f= pMatrix.min(c,d)
    min(diag(f$A))}

#find cell number

find.cell.number <- function(lower,upper,reps,probes,tum_df){
    measure = data.frame(0,0,0)
    for(i in lower:upper){
        prop = list()
        meth = list()
        for(k in 1:reps){
            num <- sample(1:ncol(tum_df), round((0.8*ncol(tum_df)),0), replace=F)
            tum = tum_df[,num]
            tum = as.matrix(tum)
            x = EDecStage1(methMixtureSamples = tum,cellTypeSpecificLoci = probes,nCts = i)
            prop[[k]] <- as.data.frame(x$proportions)
            meth[[k]] <- as.data.frame(x$methylation)}
        pro = c()
        met = c()
        comb = combn( 1:length(prop), 2)
        for(col in 1:ncol(comb)){
            met[[length(met)+1]] <- meth.cor.min(meth[[comb[1,col]]],meth[[comb[2,col]]],i)
            pro[[length(pro)+1]] <- prop.cor.min(prop[[comb[1,col]]],prop[[comb[2,col]]],i)
        }
        
        results = c(i,round(min(met),3),round(min(pro),3))
        measure = rbind(measure,results)
    }
    measure = measure[2:nrow(measure),]
    colnames(measure) = c("Cell_Number","Meth_Cor","Prop_Cor")
    return(measure)
}

find.cell.number.Return <- function(lower,upper,reps,probes,tum_df){
    measure = data.frame(0,0,0)
    for(i in lower:upper){
        prop = list()
        meth = list()
        for(k in 1:reps){
            num <- sample(1:ncol(tum_df), round((0.8*ncol(tum_df)),0), replace=F)
            tum = tum_df[,num]
            tum = as.matrix(tum)
            x = EDecStage1(methMixtureSamples = tum,cellTypeSpecificLoci = probes,nCts = i)
            prop[[k]] <- as.data.frame(x$proportions)
            meth[[k]] <- as.data.frame(x$methylation)}
        pro = c()
        met = c()
        comb = combn( 1:length(prop), 2)
        for(col in 1:ncol(comb)){
            met[[length(met)+1]] <- meth.cor.min(meth[[comb[1,col]]],meth[[comb[2,col]]],i)
            pro[[length(pro)+1]] <- prop.cor.min(prop[[comb[1,col]]],prop[[comb[2,col]]],i)
        }
        
        results = c(i,round(min(met),3),round(min(pro),3))
        measure = rbind(measure,results)
    }
    measure = measure[2:nrow(measure),]
    colnames(measure) = c("Cell_Number","Meth_Cor","Prop_Cor")
    assign("Prop_80",prop,.GlobalEnv)
    assign("Meth_80",meth,.GlobalEnv)
    return(measure)
}
############################################################################################################

get.df.rna.species.by.type <- function(atlas_all_data_list_obj_path,
                                          rna_names_dim = 2){
    
    r_dim <-rna_names_dim
    
    atls <- readRDS(atlas_all_data_list_obj_path)
    
    r_names <- lapply(FUN = dimnames, atls)
    
    r_names <- lapply(X = r_names, FUN = function(x){return(x[[r_dim]])})
    
    names(r_names) <- names(atls)
    
    return(r_names)
    
}

