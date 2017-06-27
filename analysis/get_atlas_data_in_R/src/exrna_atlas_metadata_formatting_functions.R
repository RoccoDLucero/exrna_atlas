
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
get.doctype.tables <- function(atlas_metadata, studies, doc_type, universal_fields = T){

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
    studies_dat <- atlas_metadata[studies]

    tables <- lapply(X = names(studies_dat), extract.dt.doc)

    names(tables) <- names(studies_dat)

    if(!universal_fields){return(tables)}

    univ <- get.univ.fields(tables)

    univ_dat <- lapply(X = tables, extract.univ.data, univ_fields = univ)

    univ_cmbd <- Reduce(f = cbind, x = univ_dat)

    return(univ_cmbd)
}
