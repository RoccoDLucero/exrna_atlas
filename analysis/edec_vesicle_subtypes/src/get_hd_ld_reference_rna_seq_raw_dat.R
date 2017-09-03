#get_cell_type_reference_rna_seq_raw_dat.R


################################################################################
##  LASSER HD/LD  ## 
################################################################################
# 
#
#
################################################################################

################################################################################
source("https://bioconductor.org/biocLite.R")
biocLite('SRAdb', suppressUpdates = T)
library('SRAdb')
biocLite('GEOquery', suppressUpdates = T)
library('GEOquery')
biocLite('GEOmetadb', suppressUpdates = T)
library('GEOmetadb')
library('openssl')

################################################################################

####################################################################################################
######  Connect to local SRA db copy and explore data tables                                  ######
####################################################################################################
######
sqlfile <- '../cell_type_of_origin/input/SRAmetadb.sqlite'

if(!file.exists(sqlfile)){

    sqlfile <<- getSRAdbFile(destdir = '../cell_type_of_origin/input/')
}

sra_con    <- dbConnect(drv = SQLite(), sqlfile)
######

######
sqlfile <- '../cell_type_of_origin/input/GEOmetadb.sqlite'

if(!file.exists(sqlfile)){
    getSQLiteFile(destdir = "../cell_type_of_origin/input")}

geo_meta_con <- dbConnect(SQLite(),sqlfile)
######

if(F){
#Get all table names
sra_tables <- dbListTables(sra_con)

#For tables of interest get all field names
study_fields <- dbListFields(sra_con,"study")
samp_fields  <- dbListFields(sra_con,"sample")
sra_fields  <- dbListFields(sra_con,"sra")
exp_fields  <- dbListFields(sra_con,"experiment")
run_fields  <- dbListFields(sra_con,"run")

############################
#Get all table names
geo_tables <- dbListTables(geo_meta_con)

#For tables of interest get all field names
geoConv_fields  <- dbListFields(geo_meta_con, "geoConvert")
gse_gsm_fields  <- dbListFields(geo_meta_con, "gse_gsm")
smatrix_fields  <- dbListFields(geo_meta_con, "sMatrix")
gsm_fields  <- dbListFields(geo_meta_con, "gsm")


View(dbGetQuery(geo_meta_con, "select * from gsm limit 200"))

}

query <- "select * from sra where study_accession = 'SRP090496'"

lasser_all_rna_seq_sra <- (dbGetQuery(sra_con, query ))
View(lasser_all_rna_seq_sra)

##################

####################################################################################################
######       DOWNLOAD THE SRA DATA FILES
####################################################################################################

if(F){
#dest_dir <- "../cell_type_of_origin/input/sm_RNA_Cells_and_EVs_refs/"

#getSRAfile(in_acc = GSE74759_cell_samps_sra$run_accession,
#           sra_con = sra_con, fileType = 'sra', destDir = dest_dir)

}

####################################################################################################
######      PREPARE ARCHIVES FOR EXCERPT RUNS
####################################################################################################

### CREATE ZIP ARCHIVES FOR UPLOAD TO WORKBENCH ###
#zip(zipfile = paste(dest_dir,"blood_smrna_refs_1.zip",sep = ''),
#    files = paste(dest_dir,"*.sra",sep = ''))
#the zip function needs to be pointed somehow at the 7zip executable
####################################################################################################
######  Disconnect from local SRA db copy                                                     ######
####################################################################################################
dbDisconnect(conn = sra_con)
dbDisconnect(conn = geo_meta_con)

