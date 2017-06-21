#cell_type_of_origin_functions.R



library(stringr)


make.meta.from.entex.names <- function(data){

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

    meta_df <- as.data.frame(Reduce(cbind,meta_dat), stringsAsFactors = T)
    colnames(meta_df) <- meta_vars
    rownames(meta_df) <- smp_id

return(meta_df)

}

