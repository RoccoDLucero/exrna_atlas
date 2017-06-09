
###################################################################################################
#Prepare an R data structrure that can facilitate the analysis
#These sub-lists can br nested in different ways to make data more
# amenable to analyses of different variables
clear.templates <- function(){
    kits.lst <- list( exoRneasy = NULL, miRvana = NULL )
    
    bio.conditions.lst <- list(AD = NULL, Control = NULL, GBM = NULL,
                                LGG = NULL, PD = NULL, SAH = NULL)
    
    RNA.type.lst <- list(circRNA_sense = NULL, circRNA_antisense = NULL,
                         gencode_sense = NULL, gencode_antisense = NULL,
                         gencode_sense_geneLevel = NULL, gencode_antisense_geneLevel = NULL,
                         miRNAmature_sense = NULL,
                         miRNAprecursor_sense = NULL, miRNAprecursor_antisense = NULL,
                         tRNA_sense = NULL, tRNA_antisense = NULL,
                         piRNA_sense = NULL, piRNA_antisense = NULL)
    
    assign('kits.lst',kits.lst, envir = .GlobalEnv)
    assign('bio.conditions.lst',bio.conditions.lst, envir = .GlobalEnv)
    assign('RNA.type.lst',RNA.type.lst, envir = .GlobalEnv)
}
#KIT > CONDITION > RNA.TYPE
clear.templates()
#Place RNA type under biocondition, then condition under kit
for(cond in 1:length(bio.conditions.lst)){bio.conditions.lst[[cond]] <- RNA.type.lst} 
for(kit in 1:length(kits.lst)){kits.lst[[kit]] <- bio.conditions.lst}
UH2.CSF.K.C.R.rdct <-kits.lst
clear.templates()
#str(UH2.CSF.K.C.R.rdct)

#RNA.TYPE > KIT > CONDITION
#Place condition under kit, then kit under RNA type
clear.templates()
for(kit in 1:length(kits.lst)){kits.lst[[kit]] <- bio.conditions.lst}
for(type in 1:length(RNA.type.lst)){RNA.type.lst[[type]] <- kits.lst}
UH2.CSF.R.K.C.rdct <- RNA.type.lst
#str(UH2.CSF.R.K.C.rdct)
clear.templates()
###################################################################################################
#Populate lists with DFs
#targ.lst <- UH2.CSF.R.K.C.rdct 
populate.lst <- function(targ.lst, df.path = df.pathRoot){
    success = 0
    fails = 0
    for(cur.kit in names(kits.lst)){
        fl.1 <- list.files(path = df.path, pattern = cur.kit,ignore.case = T)
        for(cur.cond in names(bio.conditions.lst)){
            dat.fldr <- fl.1[grep(cur.cond,fl.1,ignore.case = T)]
            for(cur.rna.ty in names(RNA.type.lst)){
                #cur.rna.ty.1 <- gsub(pattern = '\\.',replacement = '_',cur.rna.ty)
                ptt <- paste(cur.rna.ty,'.',sep = '')
                pth <- paste(df.path,dat.fldr,sep = '')
                
                f.cur <- list.files(path = pth,pattern = ptt,ignore.case = T)
                
                
                print(paste(cur.kit,cur.cond,cur.rna.ty))#
            
                if(length(f.cur) == 1){
                    df.cur <- read.delim(paste(pth,f.cur,sep = '/'),row.names = 1)
                
                    heir <- c(cur.kit,cur.cond,cur.rna.ty)
                    h1 <- heir[which(heir %in% names(targ.lst))]
                    h2 <- heir[which(heir %in% names(targ.lst[[1]]))]
                    h3 <- heir[which(heir %in% names(targ.lst[[1]][[1]]))]
                    
                    targ.lst[[h1]][[h2]][[h3]] <- df.cur
                    print("Added to list")
                    success = success +1
                }
                else{   fails = fails +1
                        print("File not found")
                    }
                    
            }
        }
    }
    print(paste('Succ:',success))
    print(paste('Fail:',fails))
    return(targ.lst)
}

###################################################################################################
#Clean up the folder names output by exceRpt...
#For some reason this fails on the data directly downloaded from genboree, but works if I move
#the core results by hand to a folder with a shorter path
#I will try to circumvent this later
clean.foldernames <- function(df.path = df.pathRoot){
    orig.wd <- getwd()
    setwd(df.path)
    
    my.names.orig <- list.files(path = './',pattern = '.*fastq')
    my.names.repl <- gsub(pattern = '_[ATGC]{6}.*fastq',replacement = '', x = my.names.orig)
    if(length(my.names.orig) != 0){
        for(i in 1:length(my.names.orig)){
            file.rename(from = my.names.orig[i],to = my.names.repl[i])
            
        }
    }    
    setwd(orig.wd)  
}    
  
###################################################################################################
#Use this function to add rowIDs based on RNAreferenceIDs so that the
#RNA species can be compared according to a unique identifier
addRowIdsToDFList <- function(my.inputs){
    if(!is.list(my.inputs)){
        print("Required data is a list of data frames(1)")
        return(1)}
    outputs <- list()
    for(i in 1:length(my.inputs)){
        if(!is.data.frame(my.inputs[[i]])){
            print("Required data is a list of data frames(2)")
            return(2)}
        else{       
            refID.list <- strsplit(as.character(my.inputs[[i]]$ReferenceID),':')
            #instead of rebuilding this as a matrix, only select the unique identifiers
            #This will be handled wit a switch statement
            #e.g take the MIMAT###(column3) for mature microRNA,
            #ENST for gencode
            #Gene symbol for gencode gene level
            #column1 for circular RNA
            #MI### for microRNA precursors
            #??? for piRNA
            #eg. for TRNA't RNA-Val-AAC-1-2''
            #id.idx <- length(refID.list[[1]])
            #id.mtx <- matrix(data = unlist(refID.list), ncol = id.idx,byrow = T)
            #rownames(input) <- id.mtx[,1]
            outputs[[i]] <- input
        }
    }
    return(outputs)
}
###################################################################################################