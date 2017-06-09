#Download ExRNA Atlas Metadata from the Genboree KB before proceeding#
################################################################################
my.KB.to.DF <- function(file, metadat.type){
    md <- read.table(file = file, sep = "\t", quote = '')
    #create a list of entities with all possible field names as properties
    cnames <- unique(md[,1])
    ptt <-paste("^",metadat.type,"$",sep = '')
    md.records <- md[grep(ptt,md[,1]) , ][,2]
    metadat.mtx <- matrix(data = "", nrow = length(md.records),
                          ncol = length(cnames),dimnames = list(md.records,cnames))
    
    for(rec in md.records[-length(md.records)]){
        stp <- md.records[which(md.records == rec)+1]
        a <- grep(rec,md[,2])
        b <- grep(stp,md[,2])
        for(i in a:(b-1)){
            fld <- md[i,]
            metadat.mtx[rec, as.character(fld[,1]) ] <- as.character(fld[,2])
            
            
        }
        
       
    }
    return(metadat.mtx)
}
################################################################################


Study.mtd <- my.KB.to.DF(file = "./input/Atlas_Metadata/Studies.docs.tsv",metadat.type = "Study")
Donor.mtd <- my.KB.to.DF(file = "./input/Atlas_Metadata/Donors.docs.tsv",metadat.type = "Donor")
Biosample.mtd <- my.KB.to.DF(file = "./input/Atlas_Metadata/Biosamples.docs.tsv",metadat.type = "Biosample")

Study.mtd <- as.data.frame(Study.mtd)
Biosample.mtd <- as.data.frame(Biosample.mtd)
Biosample.mtd$`Biosample.Donor ID.Age at Sampling` <- "use Donor.Age" 

Donor.mtd <- as.data.frame(Donor.mtd)
Donor.mtd$Donor.Age <- as.numeric(gsub(" y","",Donor.mtd$Donor.Age))

#Try to merge all of the metadata into a single flat table
#Combine Donor and Biosample metadata
ad <- merge.data.frame(Biosample.mtd,Donor.mtd,by.x = "Biosample.Donor ID" , by.y = "Donor" )

#Try to add Study metadata
####DO LATER####

#Combine RPM data with extended metadata
atlas.rpm.lst <- readRDS("./atlas.all.counts.and.biogps.meta")
z <- atlas.rpm.lst$miRNA
rm(atlas.rpm.lst)

ad.z <- merge(z,ad,by.x = "biosample_metadata_id" ,by.y = "Biosample" )


