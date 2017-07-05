
###################################################################################################
##title:                 exrna_atlas_metadata_formatting_functions.R
##associated_project:    exrna_atlas
##associated_analysis:   get_atlas_data_in_R
##started_date:          June 20 2017
##started_by:            Rocco Lucero
##last_updated_date:     June 26, 2017
##last_updated_by:       Rocco Lucero
###################################################################################################
##description:   functions for formatting and processing exrna atlas metadata to return
##               a single comprehensive R object
##
##inputs:                <NA>
##
##outputs:               <NA>
##
##dependencies:          <NA>
################################################################################
##Table of Contents:
##
##--extract.dt.doc
##--extract.field.names
##--extract.univ.data
##--get.doctype.tables
################################################################################
library(dplyr)
library(reshape2)
source('../get_atlas_data_in_R/src/exrna_atlas_get_metadata_functions.R')
################################################################################
extract.dt.doc <- function(st){

    st_idx <- grep(pattern = st, x = names(studies_dat))

    st_docs_lst <- studies_dat[[st_idx]]

    dt_doc_idx <- grep(pattern = doc_type, names(st_docs_lst))

    dt_doc <- st_docs_lst[[dt_doc_idx]]

    return(dt_doc)
}

################################################################################
extract.field.names <- function(dat_table){return(dat_table[, 1])}

get.univ.fields <- function(doc_type_tables_list){

    dtt <- doc_type_tables_list

    fn_lst <- lapply(X = dtt, FUN = extract.field.names)

    univ <- Reduce(f = intersect, x = fn_lst)

    return(univ)
}

################################################################################
extract.univ.data <- function(dat_table, univ_fields){

    rownames(dat_table) <- dat_table[, 1]

    dat_table <- dat_table[univ_fields, -1]

    return(dat_table)

}

################################################################################
get.ru.table <- function(atlas_metadata, studies){

    extract.dt.doc <- function(st){

        st_idx <- grep(pattern = st, x = names(studies_dat))

        st_docs_lst <- studies_dat[[st_idx]]

        dt_doc_idx <- grep(pattern = doc_type, names(st_docs_lst))

        dt_doc <- st_docs_lst[[dt_doc_idx]]

        return(dt_doc)
    }

    ############################################################################
    get.univ.fields.ru <- function(doc_type_tables_list){

        dtt <- doc_type_tables_list

        fn_lst <- lapply(X = dtt, FUN = colnames)

        univ <- Reduce(f = intersect, x = fn_lst)

        return(univ)
    }

    ############################################################################
    extract.univ.data.ru <- function(dat_table, univ_fields){

        return(dat_table[,univ_fields])

    }

    ############################################################################

    doc_type <- "RU"

    studies_dat <- atlas_metadata[studies]

    tables <- lapply(X = names(studies_dat), extract.dt.doc)

    tables <- lapply(tables, process.ru.tables)

    univ <- get.univ.fields.ru(tables)

    tables <- lapply(tables, extract.univ.data.ru, univ_fields = univ)

    names(tables) <- names(studies_dat)

    ru_cmbd <- Reduce(x = tables, f = rbind)

    return(ru_cmbd)

}
################################################################################
get.doctype.tables <- function(atlas_metadata, studies, doc_type,
                               universal_fields = T, proc_run_tables = T){

    extract.dt.doc <- function(st){

        st_idx <- grep(pattern = st, x = names(studies_dat))

        st_docs_lst <- studies_dat[[st_idx]]

        dt_doc_idx <- grep(pattern = doc_type, names(st_docs_lst))

        dt_doc <- st_docs_lst[[dt_doc_idx]]

        return(dt_doc)
    }

    ############################################################################
    extract.field.names <- function(dat_table){return(dat_table[, 1])}

    get.univ.fields <- function(doc_type_tables_list){

        dtt <- doc_type_tables_list

        fn_lst <- lapply(X = dtt, FUN = extract.field.names)

        univ <- Reduce(f = intersect, x = fn_lst)

        return(univ)
    }

    ############################################################################
    extract.univ.data <- function(dat_table, univ_fields){

        rownames(dat_table) <- dat_table[, 1]

        dat_table <- dat_table[univ_fields, -1]

        return(dat_table)

    }

    ############################################################################
    if(doc_type == 'RU' & proc_run_tables){

        ru_tbl <- get.ru.table(atlas_metadata, studies)

        return(t(ru_tbl))

    }

    studies_dat <- atlas_metadata[studies]
    print(names(studies_dat))

    #get doc_type tables for all studies
    tables <- lapply(X = names(studies_dat), extract.dt.doc)
    print(sum(sapply(tables,dim)[2,]))

    names(tables) <- names(studies_dat)

    if(!universal_fields){

        tables <- lapply(tables,t)

        return(tables)

    }

    univ <- get.univ.fields(tables)

    univ_dat <- lapply(X = tables, extract.univ.data, univ_fields = univ)

    univ_cmbd <- Reduce(f = cbind, x = univ_dat)
    print(dim(univ_cmbd))

    return(univ_cmbd)
}




################################################################################
drop.field <- function(table, dim, field_name){

    nms <- dimnames(table)[[dim]]

    fix <- grep(pattern = field_name, x = nms, invert = T)

    fix <- nms[fix]

    return(table[fix,])

}

################################################################################
reformat.names <- function(obj, dim){

    obj <- as.data.frame(obj, stringsAsFactors = F)

    dimnames(obj)[[dim]] <- make.names(names = dimnames(obj)[[dim]], unique = T)

    dimnames(obj)[[dim]] <- gsub(pattern = '^X[.]+',
                                 replacement = '',
                                 dimnames(obj)[[dim]])

    return(obj)

}

################################################################################
change.colnames <- function(table, from_name, to_name){

    cn <- colnames(table)

    cn[grep(pattern = from_name, x = cn)] <- to_name

    return(cn)

}

################################################################################
fac.to.char <- function(fac_vec){ as.character(unlist(fac_vec))}

################################################################################
process.ru.tables <- function(run_table){

    rt_vals <-  run_table[,1]

    bs_cols <-  list(  "--*- Biosample ID",
                       "--*-- File URL",
                       "--*-- File Name",
                       "--*-- MD5 Checksum",
                       "--*-- Type",
                       "--*-- DocURL",
                       "- Status")

    nsamp_row   <- grep("Raw Data Files", rt_vals, fixed = T)

    nsamp       <- as.integer(run_table[nsamp_row, 2])

    bs_rows     <- unlist(sapply(X = bs_cols, grep, x = rt_vals, fixed = T))


    run_dat <- run_table[-c(bs_rows), ]

    rownames(run_dat) <- run_dat[,1]

    run_dat <- as.data.frame(t(run_dat), stringsAsFactors = F)[2,]


    bs_tab <- run_table[bs_rows,][1:nsamp,]

    bs_tab <- cbind(bs_tab, run_dat[1,1])[,2:3]

    colnames(bs_tab) <- c("Biosample", "Run")


    processed_rt <- merge(x = bs_tab, y = run_dat, by = 'Run', all.x = T)

    return(processed_rt)
}
