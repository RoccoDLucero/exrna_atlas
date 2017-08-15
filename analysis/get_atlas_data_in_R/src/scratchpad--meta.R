

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

    if(doc_type == 'AN' & proc_run_tables){

        an_tbl <- get.ru.table(atlas_metadata, studies)

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



