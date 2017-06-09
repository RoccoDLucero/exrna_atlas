####################################################################################################
#Here we will use Lilli's modules and the the ExRNA Atlas miRNA expression raw counts to test
# the hypothesis that miRNA within the covariance modules for PD or AD
# will be correlated more highly than those outside of modules
# in a different contition e.g. cancer.
################################################################################
#Experiment 1: Graph built from healthy control count data
#Experiment 2: Graph built from alzheimer's patient data
#Experiment 3: Graph built from parkinson's patient data
#Experiment 4: Graph built from alzheimer's and healthy control patient count data
#Experiment 5: Graph built from parkinson's and healthy control patient count data
#Experiment 6: Graph built from alzheimer's and parkinson's patient count data
#Experiment 7: Graph from experiment 1, with node names permuted
#Experiment 8: Graph from experiment 2, with node names permuted
#Experiment 9: Graph from experiment 3, with node names permuted
#Experiment 10: Graph from experiment 4, with node names permuted
#Experiment 11: Graph from experiment 5, with node names permuted
#Experiment 12: Graph from experiment 6, with node names permuted
################################################################################
#Load libraries
library(lsr)
library(naturalsort)

####################################################################################################
####################################################################################################
#Load the Atlas Raw miRNA counts and associated sample metaData for ~750 samples of various
# conditions, rectify data types, and merge into a single data frame
Atlas.miRNA.meta <- read.delim2("./input/William/Atlas_Human_Sample_Metadata2.txt",sep ="\t",row.names = 1)
Atlas.miRNA.raw <- read.delim2("./input/William/Atlas_Human_miRNA_Raw_Read_Counts.txt",sep = "\t",row.names = 1)
rn <- rownames(Atlas.miRNA.raw)
Atlas.miRNA.raw <- sapply(Atlas.miRNA.raw,as.numeric)
rownames(Atlas.miRNA.raw) <- rn
Atlas.miRNA.raw <- as.data.frame(Atlas.miRNA.raw)
Atlas.miRNA <- merge.data.frame(Atlas.miRNA.meta, t(Atlas.miRNA.raw), by.x = 0, by.y = 0 )
rm(Atlas.miRNA.meta,Atlas.miRNA.raw)
apply(X = Atlas.miRNA,MARGIN = 2,FUN = class)

################################################################################

#Load Lilli's R environment from the graphical analysis:
exp.files <- list.files(path = './input', pattern = '^exrna' )

for( my.file in exp.files){  
    #my.file <- exp.files[2]
    load(file = paste('./input/', my.file,sep = '' ))
    #The vector called 'optimize.results' contains Fisher information scores for 
    # modules identified in the graph for the given comparison
    #A given miRNA may participate in multiple modules
    
    #Extract nodes from "significant" modules 
    mods.pval <- info2p(optimize.results,p.cutoff = .01)
    module.members <- get.mod.membs(mods.pval)
    all.nodes <- naturalsort(unique(unlist(module.members)))
    
    valid.miRs <- grep("hsa-mir",colnames(Atlas.miRNA),ignore.case = T)
    valid.miRs <- colnames(Atlas.miRNA[valid.miRs])
    valid.miRs <- naturalsort(valid.miRs)
    
    length(intersect(all.nodes,valid.miRs))/length(all.nodes)
    
    non.module.miRs <- setdiff(valid.miRs,all.nodes)
####################################################################################################
####################################################################################################
#For those miRNA that appear in at least one module, for each such miRNA
# identify all of the miRNA that participate in a module with it, and
# all of thos that do not appear in any modules with it.
mirs.in.ex <- per.node.int.v.ext.nodes(all.nodes,module.members)

###################################################################################################
#For those miRNA that appear in at least one significant module...     
        
    all.im.prs <- list()
    all.ex.prs <- list()
    for(m in miRTst){
        if(length(m$intramod) > 1){
            all.im.prs <- c( all.im.prs, combn(m$intramod,2,simplify = F) )
        }
        if(length(m$extramod) > 1){
        all.ex.prs <- c( all.ex.prs, combn(m$extramod,2,simplify = F) )
        }
    
    }
    
    im <- as.data.frame(t(sapply(all.im.prs,naturalsort)))
    im <- im[!duplicated(im),]    
    im <- im[naturalorder(im[,1]),]
    
    ex <- as.data.frame(t(sapply(all.ex.prs,naturalsort)))
    ex <- ex[!duplicated(ex),]
    ex <- ex[naturalorder(ex[,1]),]
    
######################
#FOR THE MIRNA THAT ARE NOT OBSERVED IN ANY MODULE
    set.seed(20161017)
    #all.bg.prs <- list()
    all.bg.prs <- combn(non.module.miRs,2,simplify = F) 
    bg <- sample(all.bg.prs, 10000)
    bg <- as.data.frame(t(sapply(bg,naturalsort)))
    
    bg <- bg[!duplicated(bg),]
    bg <- bg[naturalorder(bg[,1]),]

######################    
    
    prs <-ex
    cors.ex <- c()
    for(i in 1:nrow(prs)){
        m1 <- as.character(prs[i,1])
        m2 <- as.character(prs[i,2])
        cors.ex <- c(cors.ex,cor(Atlas.miRNA[,m1],Atlas.miRNA[,m2]))
    }
    
    prs <-im
    cors.im <- c()
    for(i in 1:nrow(prs)){
        m1 <- as.character(prs[i,1])
        m2 <- as.character(prs[i,2])
        cors.im <- c(cors.im,cor(Atlas.miRNA[,m1],Atlas.miRNA[,m2]))
    }
    
    
    prs <-bg
    cors.bg <- c()
    for(i in 1:nrow(prs)){
        m1 <- as.character(prs[i,1])
        m2 <- as.character(prs[i,2])
        cors.bg <- c(cors.bg,cor(Atlas.miRNA[,m1],Atlas.miRNA[,m2]))
    }
    txp <-.3
    #hist(cors.ex,col = rgb(1,0,0,txp),breaks = 100,xlim = c(-1,1),ylim = c(0,10),freq = F)
    #hist(cors.im,col = rgb(0,1,0,txp),add = T, breaks = 100, freq = F)
    #hist(cors.bg,col = rgb(0,0,1,txp),add = T, breaks = 100, freq = F)
    
    pdf(file = paste('./output/test.',gsub(".RData","",my.file),'.pdf',sep = '' ) )
    plot(ecdf(x = cors.ex),col = rgb(1,0,0,txp))
    plot(ecdf(x = cors.im),col = rgb(0,1,0,txp),add = T)
    plot(ecdf(x = cors.bg),col = rgb(0,0,1,txp), add =T)
    dev.off()
}
