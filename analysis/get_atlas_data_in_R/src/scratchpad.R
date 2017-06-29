



ex_ru_meta <- readRDS('../get_atlas_data_in_R/interim/exrna_atlas_metadata_EX_RU.RDS')

ex_ru_meta




####################################################################################################

process.ru.tables <- function(run_table){

        rt_vals <-  run_table[,1]

        bs_cols <-  list(  "--*- Biosample ID",
                            "--*-- File URL",
                            "--*-- File Name",
                            "--*-- MD5 Checksum",
                            "--*-- Type",
                            "--*-- DocURL")

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

s <- process.ru.tables(meta$`AKRIC1-GBM_exosome-2016-10-17`$RU)
t <- (meta$`AKRIC1-GBM_exosome-2016-10-17`$RU)



d <- get.ru.table(atlas_metadata = meta,
                        studies = dirs_studies,
                        doc_type = 'RU')


sapply(d, function(z){head(z)[,-(grep('Experimental Design', colnames(z)))]})
