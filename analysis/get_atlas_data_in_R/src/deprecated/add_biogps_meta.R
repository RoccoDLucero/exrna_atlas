#Rocco Lucero February 8 2017#
#This script should produce a single data frame containing data & metadata for all samples
#from any subset of studies in the Atlas

source('./exrna_atlas_get_data_functions.R')

#Go to Genboree FTP and look for the human studies in ExRNA-Atlas Version 4
#and display the options
ftp.site <- 'ftp://ftp.genboree.org'
subdir1 <- '/exRNA-atlas/grp/Extracellular RNA Atlas/db/exRNA Repository - hg19/file/exRNA-atlas/exceRptPipeline_v4.6.2'
base.url <- paste(ftp.site,subdir1,'/',sep = '')
dirs.studies <- getURL(url = base.url, dirlistonly = T)
dirs.studies <- unlist(strsplit(dirs.studies,'\r\n'))
study.dirname <- dirs.studies[8]
#######################################
#Get BIOGPS METADATA
#######################################
if(file.exists('./input/Atlas.all.BioGPS.metadata')){
    Atlas.all.BioGPS.metadata <- readRDS('./input/Atlas.all.BioGPS.metadata')
}else{ if(!exists('Atlas.all.BioGPS.metadata')){
            Atlas.all.BioGPS.metadata <- NULL
            for(cur.study in dirs.studies ){
                df <- my.get.study.bioGPS.metadata(base.url,cur.study)
                Atlas.all.BioGPS.metadata <- rbind(Atlas.all.BioGPS.metadata,df)
            }
            saveRDS(object = Atlas.all.BioGPS.metadata, file = './input/Atlas.all.BioGPS.metadata')
    }

}

#Atlas.all.BioGPS.metadata holds the BioGPS metadata for all included studies
dirs.studies <- setdiff(dirs.studies, unique(Atlas.all.BioGPS.metadata$Study))
for(cur.study in dirs.studies ){
    df <- my.get.study.bioGPS.metadata(base.url,cur.study)
    Atlas.all.BioGPS.metadata <- rbind(Atlas.all.BioGPS.metadata,df)
    df
}





################################################################################
#######################################
#Get and merge Read count data 
#######################################
#First put the data into a list in memory
RNA.types <- list("miRNA", "piRNA", "circularRNA", "gencode" , "tRNA")
names(RNA.types) <- RNA.types
Atlas.all.reads.per.million <- vector("list",length(RNA.types))
names(Atlas.all.reads.per.million) <- RNA.types
if(file.exists('./input/Atlas.all.reads.per.million')){
    Atlas.all.reads.per.million <- readRDS('./input/Atlas.all.reads.per.million')
}else{
    if(!exists('Atlas.all.reads.per.million')|any(unlist(lapply(Atlas.all.reads.per.million,is.null)))){
        
        
        for(ty in RNA.types){
            if(is.null(Atlas.all.reads.per.million[[ty]])){
                print(paste('start',ty))
                
                Atlas.all.reads.per.million[[ty]] <- my.combine.study.reads(dirs.studies,
                                                                            RNA.type = ty, 
                                                                            base.url,
                                                                            per.million = T,
                                                                            check.names = T)
                print(paste(ty, "completed."))
            }else{print(paste(ty,"is already populated"))
            }
        }
    }
    saveRDS(object = Atlas.all.reads.per.million, file = './input/Atlas.all.reads.per.million')
}
################################################################################
#######################################
#Get and merge Biosample
# to Sample Name mappings 
#######################################
Atlas.all.BStoSample <- vector("list",length(dirs.studies))
names(Atlas.all.BStoSample) <- dirs.studies
if(file.exists('./input/Atlas.all.BStoSample')){
    Atlas.all.BStoSample <- readRDS('./input/Atlas.all.BStoSample')
}else{
    
    for(st in dirs.studies){
        if(is.null(Atlas.all.BStoSample[[st]])){
            print(paste("starting",st))
            Atlas.all.BStoSample[[st]] <- my.map.BSIDtoSampleName(base.url,st)
        }else{print(paste(st,"is already populated"))}
    }
    saveRDS(Atlas.all.BStoSample, './Atlas.all.BStoSample')
    
    
    Atlas.all.BStoSample.merged <- as.data.frame(NULL)
    for(st.mp in Atlas.all.BStoSample){
        Atlas.all.BStoSample.merged <- rbind(Atlas.all.BStoSample.merged,st.mp)
    }
    saveRDS(Atlas.all.BStoSample.merged,'./Atlas.all.BStoSample.merged')
    write.table(Atlas.all.BStoSample.merged,'./output/Atlas.all.BStoSample.merged.txt',
                quote = F,sep = '\t',row.names = F,col.names = T)
}
    mp.renamed <- Atlas.all.BStoSample.merged
    mp.renamed$`Sample Name` <- gsub("['.']|[-]","_",mp.renamed$`Sample Name`)
    meta.map.to.samples <- merge(x = Atlas.all.BioGPS.metadata,y = mp.renamed,
               by.x = "biosample_metadata_id", by.y = "BS ID"  )
    
    ren.cols <- c(grep("sample_name",colnames(meta.map.to.samples)),
                  grep("Sample Name",colnames(meta.map.to.samples)))
    colnames(meta.map.to.samples)[ren.cols] <- c("BioGPS_Sample_NAME", "RF_Sample_NAME")
    rownames(meta.map.to.samples) <- make.unique(meta.map.to.samples$RF_Sample_NAME)
    meta.RF.names <- gsub("^sample_|_fastq$|_fq$","", meta.map.to.samples$RF_Sample_NAME, ignore.case = T)
    meta.map.to.samples$RF_Sample_NAME <- meta.RF.names
    meta.map.to.samples <- data.frame(lapply(meta.map.to.samples,factor))
    
################################################################################    
#######################################
#Merge Biosample to Sample Name mappings
#with readcount data
#######################################
atlas.all.counts.and.biogps.meta <- vector("list",length(Atlas.all.reads.per.million))
names(atlas.all.counts.and.biogps.meta) <- names(Atlas.all.reads.per.million)   
                                                   
   
for(r.type in names(Atlas.all.reads.per.million)){
    df <- NULL
    print(r.type)
    print(dim(Atlas.all.reads.per.million[[r.type]]))

    df <- Atlas.all.reads.per.million[[r.type]]
    df[,grep("^X$",colnames(df))] <- NULL
    
    df <- as.data.frame(t(df)) #Now samples are in rows, variables in columns
    
    rownames(df) <- gsub("['.']|[-]","_",rownames(df))
    colnames(df) <- gsub("^circRNA#","",colnames(df))
    df$`Sample Name` <- gsub("^sample_|_fastq$|_fq$","", rownames(df), ignore.case = T)
    df$`Sample Name` <- factor(df$`Sample Name`)
    
    df <- (merge(x =  df, y = meta.map.to.samples,
              by.x = 'Sample Name', by.y = 'RF_Sample_NAME'))
    
    atlas.all.counts.and.biogps.meta[[r.type]] <- df
    print(paste(r.type,"read count data and BioGPS metadata have been merged."))
}

lapply(atlas.all.counts.and.biogps.meta, dim)

saveRDS(object = atlas.all.counts.and.biogps.meta,file = "./atlas.all.counts.and.biogps.meta")

