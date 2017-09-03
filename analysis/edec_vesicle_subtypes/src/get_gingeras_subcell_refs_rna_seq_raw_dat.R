## get_gingeras_subcell_refs_rna_seq_raw_dat.R
## by: Rocco Lucero
## date created: August 14, 2017
################################################################################
##  Sub-cellular fractions small rna-seq data  ##
################################################################################
## SmallRNA-seq:
##   Tom Gingeras Encode GSE24565
##
####################################################################################################
source("https://bioconductor.org/biocLite.R")
biocLite('SRAdb', suppressUpdates = T)
library('SRAdb')
biocLite('GEOquery', suppressUpdates = T)
library('GEOquery')
biocLite('GEOmetadb', suppressUpdates = T)
library('GEOmetadb')
library('openssl')

####################################################################################################

####################################################################################################
######  Connect to local SRA db copy and explore data tables                                  ######
####################################################################################################

sqlfile <- '../edec_vesicle_subtypes/input/SRAmetadb.sqlite'

if(!file.exists(sqlfile)){

    sqlfile <<- getSRAdbFile(destdir = '../edec_vesicle_subtypes/input/')
}

sra_con    <- dbConnect(drv = SQLite(), sqlfile)
######

######
sqlfile <- '../edec_vesicle_subtypes/input/GEOmetadb.sqlite'

if(!file.exists(sqlfile)){
    getSQLiteFile(destdir = "../edec_vesicle_subtypes/input")}

geo_meta_con <- dbConnect(SQLite(), sqlfile)
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
    
    View(dbGetQuery(sra_con, "select * from sra limit 12"))

    #Get all table names
    geo_tables <- dbListTables(geo_meta_con)

    #For tables of interest get all field names
    geoConv_fields  <- dbListFields(geo_meta_con, "geoConvert")
    gse_gsm_fields  <- dbListFields(geo_meta_con, "gse_gsm")
    smatrix_fields  <- dbListFields(geo_meta_con, "sMatrix")
    View(dbGetQuery(geo_meta_con, "select * from gsm limit 12"))

}


####################################################################################################
###### GET INFORMATION ABOUT SMALL-RNA-SEQ READ DATA FROM SRA FOR GEO STUDIES
####################################################################################################
#study_accession_ids <- list(gingeras_encode = list(geo ="GSE24565", sra = "SRP003754"))

################################################################################
##### MANUALLY EXPLORE THE SRA AND GEO TABLES TO IDENTIFY DESIRED SAMPLES  #####
################################################################################
##########       SUBSET DATA WITHIN A GIVEN STUDY                ###############

##################
#Explore the available samples with:
#gsm_subset <- (dbGetQuery(geo_meta_con, "select * from gsm where series_id = 'GSE24565' "))
#View(gsm_subset)

query <- paste(sep = " ",
               "SELECT * ",
               "FROM",
               "sra join run on sra.run_accession=run.run_accession",
               "WHERE",
               "sra.study_accession = 'SRP003754'",
               "AND (run.experiment_name like '%cell%'",
               " OR  run.experiment_name like '%cytosol%'",
               " OR  run.experiment_name like '%nucleus%')",
               "AND sample_attribute like '%datatype: RnaSeq%'",
               "AND (sample_attribute like '%RNA < 200%' ",
               " OR  sample_attribute like '%shorter than 200 nt%')" )

GSE24565_all_samps_sra <- dbGetQuery(conn = sra_con, statement = query)

nuc <- grep("localization: nuc|localization description: Large membrane bound", GSE24565_all_samps_sra$sample_attribute)
cel.1 <- grep("localization description: Whole cell|CD34", GSE24565_all_samps_sra$sample_attribute)
cel.2 <- grep("BJ_cell", GSE24565_all_samps_sra$sample_name)
cel <- c(cel.1, cel.2)
cyt <- grep("localization: cyto|description: The fluid between", GSE24565_all_samps_sra$sample_attribute)

plya <- "rnaextract description: Rna shorter than 200 nt that has not been seperated based on Poly Adenalyation"
tap <- grep("TAP-Only", GSE24565_all_samps_sra$sample_attribute)
cip_tap <- grep("CIP-TAP", GSE24565_all_samps_sra$sample_attribute)
no_txt <- grep("No special", GSE24565_all_samps_sra$sample_attribute)

att <- GSE24565_all_samps_sra$sample_attribute
all_treatments <- unique(gsub("^.*protocol: |protocol description.*$", "", att))

## CHECK THAT ALL DESIRED SAMPLES ARE INCLUDED ##
unaccounted_samps <- GSE24565_all_samps_sra[-(c(nuc,cel,cyt)),]
#View(unaccounted_samps)
any(Reduce(intersect, list(nuc,cel,cyt)))

nuc_tap     <- GSE24565_all_samps_sra[intersect(nuc,tap),]
nuc_cip_tap <- GSE24565_all_samps_sra[intersect(nuc,cip_tap),]
nuc_no_txt  <- GSE24565_all_samps_sra[intersect(nuc,no_txt),]

cyt_tap     <- GSE24565_all_samps_sra[intersect(cyt,tap),]
cyt_cip_tap <- GSE24565_all_samps_sra[intersect(cyt,cip_tap),]
cyt_no_txt  <- GSE24565_all_samps_sra[intersect(cyt,no_txt),]

cel_tap     <- GSE24565_all_samps_sra[intersect(cel,tap),]
cel_cip_tap <- GSE24565_all_samps_sra[intersect(cel,cip_tap),]
cel_no_txt  <- GSE24565_all_samps_sra[intersect(cel,no_txt),]

sra_meta <- rbind(nuc_tap, nuc_cip_tap, nuc_no_txt,
                  cyt_tap, cyt_cip_tap, cyt_no_txt,
                  cel_tap, cel_cip_tap, cel_no_txt)

if(F){
meta_out <-"../edec_vesicle_subtypes/input/gingeras_subcell_meta.RDS"
saveRDS(object = sra_meta, file = meta_out)
}
View(sra_meta)

##################

####################################################################################################
######  DOWNLOAD THE SRA DATA FILES & PREPARE ARCHIVES
####################################################################################################
dest_dir <- "../edec_vesicle_subtypes/input/gingeras_gse24565/"
 
if(F){

    ## TAP TREATED SAMPLES
    range <- c(1:15)
    getSRAfile(in_acc = nuc_tap$run_accession[range],  #1-15 downloaded
               sra_con = sra_con, fileType = 'fastq', destDir = dest_dir)

    getSRAfile(in_acc = cyt_tap$run_accession[range],  #1-15 downloaded
               sra_con = sra_con, fileType = 'fastq', destDir = dest_dir)
    #grep("SRR527629|SRR527630|SRR527638",cyt_tap$run_accession)

    getSRAfile(in_acc = cel_tap$run_accession[range],  #1-15 downloaded
               sra_con = sra_con, fileType = 'fastq', destDir = dest_dir)

}

if(F){
    
    ## GET ALL SAMPLES FOR REMAINING SRR527* RUNS:
    get.remainingSRR537 <- function(run_accs, start_from = 1){
        #run_accs <- nuc_tap$run_accession
        rng <- start_from:length(run_accs)
        
        rst <- grep('SRR527', run_accs[rng], value = T)
        
        run_accs <- intersect(run_accs,rst)
        
        return(run_accs)
        
    }
    
    rst <- get.remainingSRR537(nuc_tap$run_accession,16)
    getSRAfile(in_acc = rst ,sra_con = sra_con, fileType = 'fastq', destDir = dest_dir)
    
    rst <- get.remainingSRR537(cyt_tap$run_accession,16)
    getSRAfile(in_acc = rst ,sra_con = sra_con, fileType = 'fastq', destDir = dest_dir)
    
    rst <- get.remainingSRR537(nuc_no_txt$run_accession)
    getSRAfile(in_acc = rst ,sra_con = sra_con, fileType = 'fastq', destDir = dest_dir)
    
    rst <- get.remainingSRR537(cyt_no_txt$run_accession)
    getSRAfile(in_acc = rst ,sra_con = sra_con, fileType = 'fastq', destDir = dest_dir)
    
    rst <- get.remainingSRR537(nuc_cip_tap$run_accession)
    getSRAfile(in_acc = rst ,sra_con = sra_con, fileType = 'fastq', destDir = dest_dir)
    
    rst <- get.remainingSRR537(cyt_cip_tap$run_accession)
    getSRAfile(in_acc = rst ,sra_con = sra_con, fileType = 'fastq', destDir = dest_dir)
  
    
}



####################################################################################################
######      PREPARE ARCHIVES FOR EXCERPT RUNS
####################################################################################################
##sAMPLES 1-15 FOR EACH
## ORGANIZE DOWNLOADED FILES BY FRACTION TYPE
sra_files <- list.files(path = dest_dir, pattern = "*.fastq.gz")

cyt_tap_files <- paste(cyt_tap$run_accession,'.fastq.gz',sep = '')
cyt_tap_files_1_15 <- intersect(cyt_tap_files, sra_files)

cel_tap_files <- paste(cel_tap$run_accession,'.fastq.gz',sep = '')
cel_tap_files_1_15 <- intersect(cel_tap_files, sra_files)

nuc_tap_files <- paste(nuc_tap$run_accession,'.fastq.gz',sep = '')
nuc_tap_files_1_15 <- intersect(nuc_tap_files, sra_files)

## REORGANIZE THE FILES INTO SUBFOLDERS BY TREATMENT AND CELLULAR FRACTION
move_files_to_subfolder(files_vec = cel_tap_files_1_15, from_path = dest_dir)

move_files_to_subfolder(files_vec = cyt_tap_files_1_15, from_path = dest_dir)

move_files_to_subfolder(files_vec = nuc_tap_files_1_15, from_path = dest_dir)


################################################################################
######CREATE ARCHIVES FOR UPLOAD TO WORKBENCH######
## CURRENTLY THIS IS DONE MANUALLY IN WINDOWS
#zip(zipfile = paste(dest_dir,"gingeras_refs_1.zip",sep = ''),
#    files = paste(dest_dir,"*.sra",sep = ''))

####################################################################################################
######  Disconnect from local SRA db copy                                                     ######
####################################################################################################
dbDisconnect(conn = sra_con)
dbDisconnect(conn = geo_meta_con)

