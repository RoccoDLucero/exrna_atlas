#get_cell_type_reference_rna_seq_raw_dat.R

################################################################################
##  BLOOD  ##
################################################################################
## SmallRNA-seq:
#   Tom Gingeras Encode GSE24565
#   B-cells 2 samples
#   CD-14 positive monocyte 2 samples
#   HSaVEC: Human Vein Endothelial Cell
#
## miRNA-seq:
#   White_blood_cells: GSE16368 #GSM669582-GSM669585
#   Erythrocytes: GSE63703
#   EBV transformed peripheral blood B lymphocytes: GSE74759; SX1421941-SRX1421974
#
#  *B-CELL GSE15229,SRX015696
#
#   SRR5223045,SRR5223046
################################################################################

################################################################################
##  LIVER  ## (see http://www.wikilectures.eu/index.php/Cells_of_Liver)
################################################################################
#### Hepatocyte
#
#### Kuppfer Cell
#
#### Hepatic Stellate Cell
#
#### Liver endothelial Cell
#
################################################################################

################################################################################
##  SPLEEN  ##
################################################################################
#### B-LYMPHOCYTE
#
#### T-LYMPHOCYTE
#
####
#
################################################################################

################################################################################
##  LUNG  ##
################################################################################
#### TYPE I ALVEOLAR EPITHELIUM
#
#### TYPE II ALVEOLAR EPITHELIUM
#
####
#
################################################################################


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

multigrep <- function(patterns, target, ...){

    xx <- sapply(X = patterns, FUN = grep, x = target)
    Reduce(f = intersect, x = xx)
}


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
study_geo_ids <- c("GSE74759","GSE24565","GSE16368","GSE63703")
study_sra_accessions <- c("SRP065890","SRP003754", "SRP046944","SRP050333")
study_identifiers_mtx <- rbind(study_sra_accessions, study_geo_ids)

samples_info_list <- apply(X = study_identifiers_mtx, MARGIN = 2,FUN = get.sra.expt.info)
names(samples_info_list) <- study_identifiers_mtx[1,]


################################################################################
##### MANUALLY EXPLORE THE SRA AND GEO TABLES TO IDENTIFY DESIRED SAMPLES  #####
################################################################################
##########       SUBSET DATA WITHIN A GIVEN STUDY                ###############

##################
#Explore the available samples for "GSE74759; SRP065890" (1 of 4):
sra_samps <- samples_info_list$SRP065890$sra_expt_samples

geo_samps <- dbGetQuery(geo_meta_con, "SELECT * FROM gsm WHERE series_id = 'GSE74759'")


#Select the desired samples (Take all HiSeq 2000 samples)
#Subset on Cellular samples
query <- paste(sep = " ",
               "SELECT *",
               "FROM sra ",
               "WHERE",
               "study_accession = 'SRP065890'",
               "AND sample_attribute LIKE '%compartment: cell%'",
               "AND instrument_model = 'Illumina HiSeq 2000'")

GSE74759_cell_samps_sra <- dbGetQuery(sra_con, query)
GSE74759_cell_samps_sra <- GSE74759_cell_samps_sra[1:6,]

#Subset on Exosomal samples
query <- paste(sep = " ",
               "SELECT *",
               "FROM sra ",
               "WHERE",
               "study_accession = 'SRP065890'",
               "AND sample_attribute LIKE '%compartment: exosome%'",
               "AND instrument_model = 'Illumina HiSeq 2000'")

GSE74759_exosome_samps_sra <- dbGetQuery(sra_con, query)
GSE74759_exosome_samps_sra <- GSE74759_exosome_samps_sra[1:6,]

##################

##################
#Explore the available samples for "GSE24565; SRP003754" (2 of 4):
sra_samps <- samples_info_list$SRP003754$sra_expt_samples

geo_samps <- dbGetQuery(geo_meta_con, "SELECT * FROM gsm WHERE series_id = 'GSE24565'")


#Select the desired samples
query <- paste(sep = " ",
               "SELECT * ",
               "FROM",
               "sra",
               "WHERE",
               "study_accession = 'SRP003754'",
               "AND sample_attribute like '%RNA < 200%'",
               "AND (sample_attribute like '%source_name: CD%'",
               "OR sample_attribute like '%cell_type: CD%'",
               "OR UPPER(sample_attribute) like '%HSAVEC%' )")

GSE24565_samps_sra <- dbGetQuery(conn = sra_con, statement = query)

##################

##################
#Explore the available samples for "GSE16368; SRP046944" (3 of 4):
sra_samps <- samples_info_list$SRP046944$sra_expt_samples

geo_samps <- dbGetQuery(geo_meta_con, "SELECT * FROM gsm WHERE series_id = 'GSE16368'")


#Select the desired samples
query <- paste(sep = " ",
               "SELECT * ",
               "FROM",
               "sra",
               "WHERE",
               "study_accession = 'SRP046944'",
               "AND library_strategy like '%RNA%'",
               "AND sample_attribute like '%blood%'",
               "AND library_layout like '%SINGLE%'")


GSE16368_samps_sra <- dbGetQuery(sra_con, query)

##################

##################
#Explore the available samples for "GSE63703; SRP050333" (4 of 4):
sra_samps <- samples_info_list$SRP050333$sra_expt_samples

geo_samps <- dbGetQuery(geo_meta_con, "SELECT * FROM gsm WHERE series_id = 'GSE63703'")


#Select the desired samples
query <- paste(sep = " ",
               "SELECT * ",
               "FROM",
               "sra",
               "WHERE",
               "study_accession = 'SRP050333'",
               "AND library_strategy = 'ncRNA-Seq'")


GSE63703_samps_sra <- dbGetQuery(sra_con, query)

##################

####################################################################################################
######       DOWNLOAD THE SRA DATA FILES
####################################################################################################

if(F){
dest_dir <- "../cell_type_of_origin/input/sm_RNA_Cells_and_EVs_refs/"

#getSRAfile(in_acc = GSE74759_cell_samps_sra$run_accession,
#           sra_con = sra_con, fileType = 'sra', destDir = dest_dir)

#getSRAfile(in_acc = GSE74759_exosome_samps_sra$run_accession,
           #sra_con = sra_con, fileType = 'sra', destDir = dest_dir)


#getSRAfile(in_acc = GSE24565_samps_sra$run  ,
#           sra_con = sra_con, fileType = 'sra', destDir = dest_dir)

##NEED to request these through DBGAP...
#getSRAfile(in_acc = GSE16368_samps_sra$run_accession,
#           sra_con = sra_con, fileType = 'sra', destDir = dest_dir)

#getSRAfile(in_acc = GSE63703_samps_sra$run_accession ,
#           sra_con = sra_con, fileType = 'sra', destDir = dest_dir)
}

####################################################################################################
######      PREPARE ARCHIVES FOR EXCERPT RUNS
####################################################################################################
######FIRST SEPARATE SAMPLES BY PLATFORM/ ADAPTER SEQUENCE######


dnld_runs <- rbind(GSE74759_cell_samps_sra,
              GSE74759_exosome_samps_sra,
              GSE24565_samps_sra,
              GSE16368_samps_sra,
              GSE63703_samps_sra)

sra_files <- list.files(path = dest_dir, pattern = "*.sra")

yet_to_dnld <- setdiff(dnld_runs$run_accession, gsub('.sra', '', sra_files) )
exclude_from_zip <- setdiff( gsub('.sra', '', sra_files), dnld_runs$run_accession )




if(F){
run_accs <- c(GSE74759_cell_samps_sra$run_accession,
              GSE74759_exosome_samps_sra$run_accession,
              GSE24565_samps_sra$run_accession,
              GSE16368_samps_sra$run_accession,
              GSE63703_samps_sra$run_accession)

run_accs <- paste(run_accs,collapse = "','")

query <-  paste(sep = "",
                "SELECT * ",
                "FROM ",
                "sra ",
                "WHERE ",
                "sample_accession IN ('", run_accs ,"')")


GSE24565_samps_dld <- dbGetQuery(conn = sra_con, statement = query)
}

######CREATE ZIP ARCHIVES FOR UPLOAD TO WORKBENCH######
zip(zipfile = paste(dest_dir,"blood_smrna_refs_1.zip",sep = ''),
    files = paste(dest_dir,"*.sra",sep = ''))

####################################################################################################
######  Disconnect from local SRA db copy                                                     ######
####################################################################################################
dbDisconnect(conn = sra_con)
dbDisconnect(conn = geo_meta_con)

