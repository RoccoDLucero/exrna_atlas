####################################################################################################
#Here we will use Lilli's modules and the the Cell ExRNA Atlas miRNA expression raw counts to test
# the hypothesis that miRNA within the covariance modules for PD or AD
# will be differentially correlated (but more highly than those outside of modules)
# in a different cell types
#The idea here is that if a given cell type contributes to a given disease profile/phenotype
##that cell type will "resonate" more with with a network learned from patients with that disease
##Alternatively differnent modules in general may be associated with different cell types.
##
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
#Load libraries and functions
#
source('./exrna_atlas_functions.R',echo = F)
####################################################################################################
#Load the Cell Atlas Raw miRNA counts and rectify data types.
#
Cell.Atlas <- read.delim2("./input/Merged_results_joel/exceRpt_miRNA_ReadCounts.txt",
                          sep = "\t",row.names = 1,stringsAsFactors = F)

Cell.Atlas <- to.numeric.df(tFrame(Cell.Atlas))
head(Cell.Atlas[,1:4])
####################################################################################################
#Load Lilli's R environment from the graphical analysis:
#
exp.files <- list.files(path = './input/Lilli/', pattern = '^[ek][xr][ra][ns]' )
#smpl <- sample(1:length(exp.files),size = 4)
for(my.file in exp.files){
    #my.file <- exp.files[1]
    load(file = paste('./input/Lilli/', my.file, sep = '' ))
    #The vector called 'optimize.results' contains Fisher information scores for 
    # modules identified in the graph for the given comparison
    # A given miRNA may participate in multiple modules
    print('get.sigmods')
    ###########################################
    #Extract nodes from "significant" modules
    p.cut <- .05
    ##For nodes in significant modules
    sig.mods <- info2p(optimize.results, p.cutoff=p.cut, '<=')
    module.members <- get.mod.membs(sig.mods)
    all.sig.nodes <- naturalsort(unique(unlist(module.members)))
    
    ##For non-significant nodes in the graph
    nonsig.mods <- info2p(optimize.results, p.cutoff=p.cut, '>')
    non.sig.nodes <- get.mod.membs(sig.mods)
    non.sig.nodes<- naturalsort(unique(unlist(module.members)))
    
    
    ###########################################
    #Identify all miRNA which do not appear
    # as nodes in the "significant" modules
    #but for which we have mapped reads
    non.module.miRs <- get.non.module.vars(Cell.Atlas, all.sig.nodes,
                                           non.sig.nodes,"hsa-mir")
    if( (length(all.sig.nodes) == 0) ){
        print(paste( my.file,"Has Failed"))
        next
    } 
    
    print('in.vs.ex')
    #Identify, for all nodes that participate in a 'significant' module
    # those nodes that are in a module with that node
    # those nodes that are not in a module with that node
    sig.mod.in.ex.pr.lst <- per.node.int.v.ext.nodes(all.sig.nodes,module.members)
    
    module.node.pairings <- get.reduced.pairings(sig.mod.in.ex.pr.lst)
    
    bg.smpl <- 1000
    background.node.pairings <- get.bg.pairings(non.module.miRs, bg.smpl)
    
    cors.im <- get.nodepair.cors(Cell.Atlas, module.node.pairings$intramodule.prs)
    cors.ex <- get.nodepair.cors(Cell.Atlas, module.node.pairings$extramodule.prs)
    cors.bg <- get.nodepair.cors(Cell.Atlas, background.node.pairings)
    if( length(cors.im) == 0){
        paste( my.file,"Has Failed")
        next
    }
    if( length(cors.bg) == 0){
        paste( my.file,"Has Failed")
        next
    }
    
    print('plotting')
    #Plotting Details
    cur.expt <-  gsub('.RData','', my.file)
    txp <-.5
    hst.brks <- 50
    ks.res <- ks.test(cors.im, cors.bg,alternative = "less")
    options(scipen = 5)
    options(digits=4)
    #pdf(file = paste('./output/Cell_Atlas/.',cur.expt,'.pdf',sep = '' ) )
    pdf(file = paste('./output/Cell_Atlas/',cur.expt,'.pdf',sep = '' ) )
    
    hist(cors.im,col = rgb(1,0,0,txp),breaks = hst.brks,
         xlim = c(-0.5,1),ylim = c(0,15),freq = F,
         main = paste('Correlations of module vs. non-module node pairs\nfor',cur.expt),
         xlab = 'between nodes correlation', ylab = 'Probability density as %' )
    
    hist(cors.bg,col = rgb(0,0,1,txp),add = T, breaks = hst.brks, freq = F)
    
    legend(x = 'topleft',legend = c(paste('intramodule n =',length(cors.im)),
                                    paste('non-module n =',bg.smpl)) ,
           fill = c(rgb(1,0,0,txp),rgb(0,0,1,txp)) )
    mtext(paste('KS-test D:',round(ks.res$statistic,2)), line = -1,side = 3,adj = 1)
    mtext(paste('p-val:',signif(ks.res$p.value,3)), line = -2,side = 3,adj = 1)
    #mtext(paste('nod:',signif(ks.res$p.value,3)), line = -3,side = 3,adj = 1)

    dev.off()
    
}

    
################################################################################
if(F){    
#For each miRNA in current module, identify all miRNA that are not in a 
# module with ANY of those module members"
valid.miRs <- colnames(Atlas.miRNA)
tested.pairs <- list()
my.cors <- c()
for(cur.mod in module.members){
    non.mod <- setdiff(all.nodes, cur.mod)
    #my.cors <- c()
    for(cur.miR.intr in cur.mod){
        
        for(cur.miR.extr in non.mod){
            if(cur.miR.extr %in% tested.pairs[[cur.miR.intr]]){next}
            if(!(cur.miR.extr %in% valid.miRs && cur.miR.extr %in% valid.miRs)){next}
            cr <- cor(Atlas.miRNA[,cur.miR.intr],Atlas.miRNA[,cur.miR.extr])
            my.cors <- c( my.cors, cr)
            tested.pairs[[deparse(cur.miR.intr)]] <- c(tested.pairs[[deparse(cur.miR.intr)]],cur.miR.extr)
        }
    }
    #print(c(tst.mod,length(module.members)))
    #if(tst.mod > length(module.members)){break}
}
length(my.cors)
hist(my.cors, col = "skyblue", xlim = c(-1,1))
print(length(module.members))

for(cur.mod in module.members[1:100]){
    tst.miRNA <- c()
    for(tst.mod in module.members[1:100]){
        if(cur.mod == tst.mod){next}
        if( length(intersect(cur.mod,tst.mod)) == 0 ){
            tst.miRNA <- c(tst.miRNA,tst.mod)
        }
    }
    #print(cur.mod)
    print(length(tst.miRNA))
    #print(tst.miRNA)
    
}

length(module.members[grep(pattern = 'hsa-miR-374a-3p',x = module.members)])

#To select inter-module
####################################################################################################
#For the miRNA identified in a given module, check the correlations within and between
#members of the module.
#
#Subset the data to include only the samples for the condtions in which we are interested.
table( (Atlas.miRNA$Condition) )
tst.conditions <- c('Carcinoma','Alzheimer', 'Parkinson', 'Healthy Control')

conds.lst <- list()
for(my.condition in tst.conditions){ 
    my.Atlas.subset <- Atlas.miRNA[grep( my.condition ,Atlas.miRNA$Condition,fixed = T ),]
    Cor.lst <- list()
    #Cor.lst[['Query.Condition']] <- my.condition
    Cor.lst[['Pairwise.Correlations']] <- list()
    for( mod in module.members ){
        
        if( length(mod) > 1 ){
            intra.module.pairs <- combn(mod,2)
            
            for( pair in 1:ncol(intra.module.pairs) ){
                cur.pair <- intra.module.pairs[,pair]
                cur.pair <- sapply(cur.pair,strsplit,"\\|",simplify = T)
                cur.pair <- c(cur.pair[[1]][1],cur.pair[[2]][1])
                if( !all(cur.pair %in% colnames(Atlas.miRNA)) ){ next }
                pair.cor <- cor( my.Atlas.subset[,cur.pair[1]], my.Atlas.subset[,cur.pair[2]],use = 'pairwise' )
                pairname <- paste(cur.pair[1],cur.pair[2],sep = '&')
                if( !is.na(pair.cor) ){ Cor.lst[['Pairwise.Correlations']][[pairname]] <- pair.cor }
            }   
        }
        
    }
    
    conds.lst[[my.condition]] <- Cor.lst
}

rm(cor.df)
for( cond in names(conds.lst) ){
    cor.tmp <- as.numeric(unlist(conds.lst[[cond]]$Pairwise.Correlations))
    cat( cond, ": ",length(cor.tmp)," - ",length(names(conds.lst[[cond]]$Pairwise.Correlations)) )
    names(cor.tmp) <- names(conds.lst[[cond]]$Pairwise.Correlations)
    cor.tmp <- as.data.frame(cor.tmp)
    colnames(cor.tmp) <- cond
    if( !exists('cor.df')){
        cor.df<- cor.tmp
        next
    }
    if( ncol(cor.df)==1 ){
        cor.df <- merge.data.frame(cor.df,cor.tmp,by.x = 0, by.y = 0 )
        next
    }
    if( ncol(cor.df)>1 ){
        cor.df <- merge.data.frame(cor.df,as.data.frame(cor.tmp),by.x = 1, by.y = 0 )
    }
    
}
rm(cor.tmp)
head(cor.df)
}
################################################################################
