#cell_type_of_origin_functions.R



library(stringr)


make.meta.from.entex.names <- function(data){

    flatten <- function(col_list){

        sapply(X = col_list, FUN = function(x){paste(x, collapse = " & ")})

    }
    ############################################################################
    sn <- colnames(data)

    group <- str_extract(string = sn, pattern = "adult|fetal")

    smp_id <- str_extract(string = sn, pattern = "^([^_]+)")

    smp_src <- gsub(pattern = "^.*(year|week|unknown)_", replacement = "", x = sn  )
    smp_src <- gsub(pattern = "_(adult|fetal)", replacement = "", x = smp_src)

    time_unit <- str_extract(string = sn, pattern = "year|week|unknown")
    age <- gsub(pattern = "^.*male_", replacement = "", x = sn)
    age <- str_extract(string = age, pattern = "^.*(year|week|unknown)")
    age <- str_extract_all(string = age, pattern = "[0-9]+")


    sex <- str_extract_all(string = sn, pattern = "(male|female)")
    #sex <- sapply(X = sex, FUN = function(x){paste(x, collapse = " & ")})

    meta_vars <- c('group', 'sex', 'age', 'time_unit', 'smp_src')
    #colnames(meta_dat)  <- meta_vars
    #rownames(meta_dat)  <-  smp_id

    meta_dat            <- vector(mode = "list", length = length(meta_vars))
    names(meta_dat)     <- meta_vars

    meta_dat$group      <- group
    meta_dat$sex        <- sex
    meta_dat$age        <- age
    meta_dat$time_unit  <- time_unit
    meta_dat$smp_src    <- smp_src

    meta_dat <- lapply(meta_dat,flatten)

    meta_df <- as.data.frame(Reduce(cbind, meta_dat), stringsAsFactors = T)
    colnames(meta_df) <- meta_vars
    rownames(meta_df) <- smp_id



return(meta_df)

}




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



get.normalized.readcounts.obj <- function(df){

    df <- as.matrix(df, row.names = rownames(df), col.names = colnames(df))

    df_qn <- get.qn(df)
    rownames(df_qn) <- rownames(df)
    colnames(df_qn) <- colnames(df)

    df_npn <- huge.npn(df_qn)

    df_npn_lgstc <- apply(X = df_npn, MARGIN = 2, FUN = txfmLogistic1)
    df_npn_lgstc <- df_npn_lgstc[ , complete.cases(t(df_npn_lgstc))]

    normalized_lst <- list(qn = df_qn,
                           qn.npn = df_npn,
                           qn.npn.lgstc =  df_npn_lgstc,
                           original = df)

    return(normalized_lst)
}




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


subset.on.shared.probes <- function(ref_df, exp_df){
    avail_probes <- intersect(colnames(ref_df), colnames(exp_df))
    ref_df <- ref_df[, avail_probes]
    exp_df <- exp_df[, avail_probes]

    return(list(ref =  ref_df, exp = exp_df))
}

make.colors.vector <- function(factor_vec, color_fun = rainbow){

    color.function <- color_fun

    colors <- as.factor(factor_vec)

    levels(colors) <- color.function(length(colors))

    return(as.character(colors))

}

perform.probe.selection <- function(ref_probes, ref_classes, max_p_val = 1e-3,
                                    n_markers = 100, method = c("one.vs.rest", "each.pair" )){

    info_probes <- lapply(X = method, run_edec_stage_0)

    names(info_probes) <- method

    if(length(info_probes) > 1){info_probes[['max_set']] <- Reduce(f = union, x = info_probes)}

    return(info_probes)

}

get.correlation.heatmaps <- ( )





######################################################################################################
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



