####################################################################################################
#Here we will use Lilli's modules and the the Cell ExRNA Atlas miRNA expression raw counts to test
# the hypothesis that miRNA within the covariance modules for PD or AD
# will be differentially correlated (but more highly than those outside of modules)
# in a different cell types
#The idea here is that if a given cell type contributes to a given disease profile/phenotype
##that cell type will "resonate" more with with a network learned from patients with that disease
##Alternatively differnent modules in general may be associated with different cell types.
##
####################################################################################################
####################################################################################################
#Load libraries and functions
#
source('./exrna_atlas_functions.R',echo = F)
####################################################################################################
##Import miRBASE 21 and select only human miRNA (should be 1881 miRNA)
miRBASE21 <- read.delim("./input/miRNA.txt", stringsAsFactors=FALSE)
miRBASE21.hsa <- (miRBASE21[grep('hsa-',miRBASE21$ID),grep('ID',colnames(miRBASE21))])
miR.NAMES.hsa <- union(union(miRBASE21.hsa$ID,miRBASE21.hsa$Mature1_ID),miRBASE21.hsa$Mature2_ID)
miR.NAMES.hsa <- (sort(miR.NAMES.hsa))

#Load ExRNA Atlas RAW readcounts from William (11.28.2016 1465 samples)
atlas <- read.table('./input/William/Atlas_Human_miRNA_Read_Counts_11_28_16.txt')
atlas <- as.matrix(atlas)
colnames(atlas) <- gsub('sample_','',colnames(atlas),ignore.case = T)
colnames(atlas) <- gsub('sample_','',colnames(atlas),ignore.case = T)
colnames(atlas) <- gsub('_fastq','',colnames(atlas),ignore.case = T)
colnames(atlas) <- gsub('_sequence','',colnames(atlas),ignore.case = T)
colnames(atlas) <- gsub('_clean_fq','',colnames(atlas),ignore.case = T)

QN <- T
if(QN){
    atlas.qn <- normalize.quantiles(atlas)
    dimnames(atlas.qn) <- dimnames(atlas)
    #To Reassure myself that the data is properly normalized:
    #pdf(file = './output/qnbox.pdf',width = 20,height = 10,onefile = T)
    #boxplot.matrix(atlas.qn[,1:20], use.cols = T)
    #boxplot.matrix(atlas[,1:20], use.cols = T)
    #dev.off()
}


my.query.dat <- atlas.qn
####################################################################################################
#Load Lilli's R environment from the graphical analysis:
#Here we use the network modules built with k = 10 depth for nov_update

mods.dat.dir <- './input/Lilli/nov_update/'
exp.files <- list.files(path = mods.dat.dir, pattern = '^[ek][xr][ra][ns]' )
exp.files <- naturalsort(exp.files)
rn <- gsub('Prostrate','Prostate', (gsub('.RData','',(make.names(exp.files)))))


out.put.lst <- vector('list', length(exp.files) )
mir.mod.mtx <- as.data.frame(matrix(nrow = length(rownames(atlas)),ncol = length(exp.files) ))
rownames(mir.mod.mtx) <- rownames(atlas) # update colnames as each expt is processed

for(my.file in exp.files){
    load(file = paste(mods.dat.dir, my.file, sep = '' ))
    expt.name <- gsub('.RData','', make.names(my.file) )
    expt.name <-gsub('Prostrate','Prostate',expt.name)
    #The vector called 'optimize.results' contains Fisher information scores for 
    # modules identified in the graph for the given comparison
    # A given miRNA may participate in multiple modules
    ###########################################
    #Extract nodes from "significant" modules
    p.cut <- .05
    
    ##For nodes in significant modules
    sig.mods <- info2p(optimize.results, p.cutoff=p.cut, '<=')
    module.members <- get.mod.membs(sig.mods)
    module.members <- lapply(module.members,naturalsort)
    all.sig.nodes <- naturalsort(unique(unlist(module.members)))
    
    mod.vec <- vector(length=length(rownames(atlas)))
    for(nd in rownames(mir.mod.mtx)){

        if(nd %in% all.sig.nodes){mod.vec[which(rownames(mir.mod.mtx)==nd)] <- 1 }
        else{mod.vec[which(rownames(mir.mod.mtx)==nd)] <- 0 }
    }
    
    mir.mod.mtx[,which(exp.files == my.file)] <- mod.vec
    colnames(mir.mod.mtx)[which(exp.files == my.file)] <- expt.name
    
    
    out.put.lst[[expt.name]] <- rmv.sub.modules(module.members)
    out.put.lst <- out.put.lst[!sapply(out.put.lst,is.null)]
    
}


#For the miRNA that appear in modules do we observe clustering by experiments
colnames( mir.mod.mtx )
mir.mod.mtx <- mir.mod.mtx[unique(unlist(out.put.lst)),]
mir.mod.mtx <- mir.mod.mtx[complete.cases(mir.mod.mtx),]
rownames(mir.mod.mtx)[1:20]

#Based on the non-subset modules identified in each experiment create the 
#list of unique modules identified in all experiments
unq.modules <- unique(unlist(out.put.lst,recursive = F))
unq.modules <- unique(sapply(unq.modules,paste,collapse = '.'))

modules.table <- as.data.frame(matrix(nrow = length(unq.modules),ncol = length(out.put.lst)))
rownames(modules.table) <- unq.modules # update colnames as each expt is processed
colnames(modules.table) <- names(out.put.lst)
idx <-0
for(ex in out.put.lst){
    idx <- idx+1
    my.t <- sapply(ex,paste,collapse = '.')
    my.sum <- c()
    for(m in unq.modules){
        if( m %in% my.t){my.sum <-c(my.sum,1)}else{my.sum <-c(my.sum,0)}
    }
    if(sum(my.sum) == length(my.t)){
        modules.table[,idx] <- my.sum
    }
}

data <- as.matrix(mir.mod.mtx)
data <- as.matrix(modules.table)
data <- as.matrix(modules.table[,grep("AD|PD",colnames(modules.table))])
data <- as.matrix(modules.table[,grep("ras",colnames(modules.table))])
heatmap.2(data, trace="none", scale="row",cexRow = .1,cexCol = .8)


module.sizes <- 1+sapply(rownames(modules.table),str_count,'\\.')
module.recurrence <- rowSums(modules.table)
df <- as.data.frame(cbind(module.sizes,module.recurrence))

pdf(file = paste(plots.dir,'size.v.recurrence.pdf'),width = 11,height = 8.5)
ggplot(df,aes(x=module.sizes,y=module.recurrence)) + stat_binhex() + ggtitle('Module recurrence vs. module size\n in K-JENS and TPATEL datasets (module size 1-10)')
dev.off()

x <- rownames(modules.table[which(module.recurrence > 1),])
z <- modules.table[x,]
module.sizes.z <- 1+sapply(rownames(z),str_count,'\\.')
module.recurrence.z <- rowSums(z)

if(F){
    my.file <- './output/ExRNA_Atlas/recurrentMods.txt'
    for(i in 1:nrow(z)){
        cat('============',file = my.file, sep = '\n',append = T)
        cat(rownames(z)[i],file = my.file, sep = '\n',append = T)
        cat(colnames(z[i,z[i,]>0]),file = my.file, sep = ', ',append = T)
        cat('\n',file = my.file,append = T)
    }
}

hist(module.recurrence,breaks = seq(0,12))
hist(module.sizes)
png(file="./output/mygraphic.png",width=4000,height=5500)
heatmap.2(atlas.qn[1:30,1:20], trace = 'no')
dev.off()

#Test the correlation of module miRNA in KJENS data as positive control for
# Alzheimer's and Parkinson's modules:
jensen.subset <- atlas.qn[,grep('AD|PD|CONTROL',colnames(atlas.qn))]



colnames(atlas.qn)

