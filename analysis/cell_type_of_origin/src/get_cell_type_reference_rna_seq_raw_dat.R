#get_cell_type_reference_rna_seq.R

##dim(smp_rs_homo[grep('PRJNA388448',smp_rs_homo$xref_link),]) EV and Cell DATA...KEEP##

#See:
#http://www.bioconductor.org/packages/release/bioc/vignettes/SRAdb/inst/doc/SRAdb.pdf
####################################################################################################
source("https://bioconductor.org/biocLite.R")
biocLite('SRAdb', suppressUpdates = T)
library('SRAdb')

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

#dbGetQuery(sra_con,'PRAGMA TABLE_INFO(study)')

colDesc <- colDescriptions(sra_con = sra_con)
colnames(colDesc)
colDesc[,1:4]

#Select from the sra_table Study three records
std_rs <- dbGetQuery(conn = sra_con, statement = "select * from study limit 3")
colnames(std_rs)
std_rs[,5]

rs <- dbGetQuery(conn = sra_con,
                 statement =  paste( "select study_accession",
                                     "study_title from study where",
                                     "study_description like 'Transcriptome%'", sep=" "))
rs


getTableCounts <- function(tableName,conn){
    sql <- sprintf("select count(*) from %s",tableName)
    return(dbGetQuery(conn,sql)[1,1])
}

c <- do.call(what = rbind, args = sapply(sra_tables[c(2,4,5,11,12)],
                                     getTableCounts, sra_con, simplify=FALSE))

rs <- dbGetQuery(conn = sra_con,
                 statement = paste( "SELECT study_type AS StudyType, ",
                        "count( * ) AS Number FROM `study` ",
                        "GROUP BY study_type order by Number DESC", sep=""))
rs

##########
smp_rs <- dbGetQuery(conn = sra_con,
                     statement = "select * from sample limit 15")
run_rs <- dbGetQuery(conn = sra_con,
                     statement = "select * from run limit 5")

expt_rs <- dbGetQuery(conn = sra_con,
                     statement = "select * from experiment limit 5")
colnames(smp_rs)

smp_rs_homo <- dbGetQuery(conn = sra_con,
                 statement =  paste( "select *",
                                     "from sample where",
                                     "scientific_name like 'Homo%'", sep=" "))


exp_rs_lib_strat <- dbGetQuery(conn = sra_con,
                          statement =  paste( "select *",
                                              "from experiment where",
                                              "library_strategy like 'RNA-seq%'", sep=" "))

dim(smp_rs_homo)
head(smp_rs_homo)
dim(exp_rs_lib_strat)
head(exp_rs_lib_strat)
head(exp_rs_lib_strat$library_strategy)
homo_rseq <- intersect(smp_rs_homo$sample_accession, exp_rs_lib_strat$sample_accession)

SAMN07177409


rs1 <- getSRA (search_terms ='liver* NEAR/2 cell*', out_types=c('run','study', 'sample', 'experiment'), sra_con=sra_con)

liver_samps_homo <- intersect(smp_rs_homo$sample_alias, rs1$sample_alias)
liver_samps_r_seq <- intersect(exp_rs_lib_strat$sample_name, rs1$sample_name)
homo_r_sreq_samps <- intersect(liver_samps_homo, liver_samps_r_seq)
tt <- rs1[which(rs1$sample_name %in% homo_r_sreq_samps),]
rr <- tt[which(tt$study_type == "Transcriptome Analysis"),]

rs2 <- getSRA (search_terms ='kidney* NEAR/2 cell*', out_types=c('run','study', 'sample', 'experiment'), sra_con=sra_con)
kdny_samps_homo <- intersect(smp_rs_homo$sample_alias, rs2$sample_alias)
kdny_samps_r_seq <- intersect(exp_rs_lib_strat$sample_name, rs2$sample_name)
homo_r_sreq_samps <- intersect(kdny_samps_homo, kdny_samps_r_seq)
tt <- rs1[which(rs2$sample_name %in% homo_r_sreq_samps),]
rr <- tt[which(tt$study_type == "Transcriptome Analysis"),]


kidney_cl_and_ev <- (smp_rs_homo[grep('PRJNA388448',smp_rs_homo$xref_link),])

accsns <- sapply( (seq(25,42)), function(x){paste("SRS22391", x, sep = '')})

rs = listSRAfile( accsns, sra_con, fileType = 'sra' )
rs
rs = getSRAinfo ( accsns, sra_con, sraType = "sra" )
rs

#getSRAfile( accsns, sra_con, fileType = 'sra' )

