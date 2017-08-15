#compare_evs_to_their_cells_of_origin.R
#Data come from study GSE99430

####################################################################################################
source("https://bioconductor.org/biocLite.R")
biocLite('SRAdb', suppressUpdates = T)
library('SRAdb')
biocLite('GEOquery', suppressUpdates = T)
library('GEOquery')
####################################################################################################

sqlfile <- '../cell_type_of_origin/input/SRAmetadb.sqlite'

if(!file.exists('../cell_type_of_origin/input/SRAmetadb.sqlite')){

    sqlfile <<- getSRAdbFile(destdir = '../cell_type_of_origin/input/')
}

sra_con    <- dbConnect(drv = SQLite(), sqlfile)

#Get all table names
sra_tables <- dbListTables(sra_con)

#For tables of interest get all field names
study_fields <- dbListFields(sra_con,"study")
samp_fields  <- dbListFields(sra_con,"sample")
sra_fields  <- dbListFields(sra_con,"sra")
exp_fields  <- dbListFields(sra_con,"experiment")

colDesc <- colDescriptions(sra_con = sra_con)

################################################################################
#GET RAW SMALL-RNA-SEQ READ DATA FROM SRA FOR GEO STUDY:
study_name <- "'GSE99430'"

exp_rs_ID <- dbGetQuery(conn = sra_con,
                        statement =  paste( "select * from experiment",
                                            "where study_name = ", study_name,
                                            sep=" "))

#GET ADDTIONAL INFROMATION ABOUT THE SAMPLES
smp_rs_ID <- dbGetQuery(conn = sra_con,
                        statement =  paste( "select * from sample inner join experiment",
                                            "on sample.sample_accession=experiment.sample_accession",
                                            "where study_name=",study_name, sep = " "))

accsns <- smp_rs_ID$sample_accession

study_sra_info <- getSRAinfo ( accsns, sra_con, sraType = "sra" )

getSRAfile(in_acc =  accsns, sra_con =  sra_con,
           destDir = "../cell_type_of_origin/input/sm_RNA_Cells_and_EVs/",
           fileType = 'sra')
###########################################################################################################
#GET METADATA TABLE FROM GEO STUDY
geo_dat <- getGEO(GEO = gsub(pattern = "'", replacement = "", study_name))
View(geo_dat$GSE99430_series_matrix.txt.gz)

