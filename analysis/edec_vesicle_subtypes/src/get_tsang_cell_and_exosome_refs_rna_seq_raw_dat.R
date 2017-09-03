## get_tsang_cell_and_exosome_refs_rna_seq_raw_dat.R
## by: Rocco Lucero
## date created: August 21, 2017
####################################################################################################
source("https://bioconductor.org/biocLite.R")
biocLite('SRAdb', suppressUpdates = T)
library('SRAdb')
biocLite('GEOquery', suppressUpdates = T)
library('GEOquery')
biocLite('GEOmetadb', suppressUpdates = T)
library('GEOmetadb')
library('openssl')



dest_dir <- "../edec_vesicle_subtypes/input/tsang/"
####################################################################################################
######  Connect to local SRA db copy and explore data tables                                  ######
####################################################################################################

sqlfile <- "../edec_vesicle_subtypes/input/SRAmetadb.sqlite"
sra_con <- dbConnect(drv = SQLite(), sqlfile)


##################
#Explore the available samples with:
query <- paste(sep = " ",
               "SELECT * ",
               "FROM",
               "sra join run on sra.run_accession=run.run_accession",
               "WHERE",
               "sra.study_accession = 'SRP065890'"
               )

all_samps_sra <- dbGetQuery(conn = sra_con, statement = query)
View(all_samps_sra)

meta_out <- paste(dest_dir,"tsang_meta.RDS", sep = '')
saveRDS(object = all_samps_sra, file = meta_out)
#View(all_samps_sra)

####################################################################################################
######  DOWNLOAD THE SRA DATA FILES & PREPARE ARCHIVES
####################################################################################################

mi_seq_samps <- with(data = all_samps_sra, expr = instrument_model == 'Illumina MiSeq')
mi_seq_samps <- all_samps_sra[mi_seq_samps,]

mi_seq_dir <- paste(dest_dir, "mi_seq", sep = '')
getSRAfile(in_acc = mi_seq_samps$run_accession, sra_con = sra_con, fileType = 'fastq', destDir = mi_seq_dir)

hi_seq_samps <- with(data = all_samps_sra, expr = instrument_model == 'Illumina HiSeq 2000') 
hi_seq_samps <- all_samps_sra[hi_seq_samps,]

hi_seq_dir <- paste(dest_dir, "hi_seq", sep = '')
getSRAfile(in_acc = hi_seq_samps$run_accession, sra_con = sra_con, fileType = 'fastq', destDir = hi_seq_dir)

####################################################################################################
######  Disconnect from local SRA db copy                                                     ######
####################################################################################################
dbDisconnect(conn = sra_con)
dbDisconnect(conn = geo_meta_con)

