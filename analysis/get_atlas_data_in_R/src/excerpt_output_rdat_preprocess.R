#excerpt_output_rdat_preprocess.R
#started: unknown
#last_updated: June 8, 2017

################################################################################
################################################################################
#FUNCTIONS FOR WORKING WITH EXCERPT DATA GENERATED ON THE GENBOREE WORKBENCH
################################################################################
atlas_load_excerpt_readcounts <- function(rna.types.lst = NULL, per.million = T,
                                          add.new.studies = T, rdcnt.obj.pth){
    
    #This script needs to be more generic so that it can handle readcount
    #objects of non-RDS type
    
    add.atlas.readcounts <- function(readcounts.lst, rna.types){
        readcounts.lst <- vector("list", length(rna.types))
        
        for(ty in rna.types){
            if(is.null(readcounts.lst[[ty]])){
                
                print(paste('start',ty))
                
                x <- my.combine.study.reads(dirs.studies, rna.type = ty,
                                            studies.url, per.million = T, check.names = T)
                
                readcounts.lst[[ty]] <- x
                
                print(paste(ty, "completed."))
                
            }else{print(paste(ty,"is already populated"))}
            
        }
        
    }
    
    rna.types <- list("miRNA", "piRNA", "circularRNA", "gencode" , "tRNA")
    
    if(!is.null(rna.types.lst)){rna.types <- rna.types.lst
    print("Custom RNA types selected.")
    
    }
    
    names(rna.types) <- rna.types
    print("Selected RNA types:")
    print(names(rna.types))
    
    
    if(file.exists(rdcnt.obj.pth)){
        
        atlas.readcounts <- readRDS(rdcnt.obj.pth)
        print("Loading pre compiled Atlas read counts data")
        
        if( any(sapply(atlas.readcounts, is.null)) ){
            rna.to.update <- setdiff(names(rna.types.lst),names(atlas.readcounts)) 
            atlas.readcounts <- add.atlas.readcounts(atlas.readcounts, rna.types)
            
        }
        
        
    }
    
    if(!exists('atlas.readcounts') | any(sapply(atlas.readcounts, is.null))){
        
        atlas.readcounts <- add.atlas.readcounts(atlas.readcounts, rna.types)
        
    }
    
    saveRDS(atlas.readcounts, rdcnt.obj.pth)
}

################################################################################
my.load.excerpt.rdata <- function(local_rdat_file){
    #Input: (Single File) An RDATA Structure Produced by Genboree ExceRpt Jobs.
    #       Downloaded from the Genboree Workbench
    #       e.g. Lasser1_Large_2017-5-9_exceRpt_smallRNAQuants_ReadCounts.RData
    #
    #Output: (single R object) A list of dataframes of RNA readcounts,
    #       each data frame corresponding to one of the databases used 
    #       for mapping (Typically a single RNA type)
    #       e.g. miRNA, piRNA, gencode(contains several RNA types), exogenous...    
    ############################################################################
    if(!file.exists(local_rdat_file)) {
        
        print("File not found.")
        
        return()
    }
    
    loaded <- load(local_rdat_file)
    
    e <- environment()
    d <- lapply(X = loaded,FUN = get, envir = e)
    
    names(d) <- sapply(loaded,gsub,replacement = "", pattern = "exprs.")
    d <- lapply(d,t)
    d <- lapply(d,as.data.frame)
    
    return(d)
}





