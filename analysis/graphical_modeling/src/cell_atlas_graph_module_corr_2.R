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
#Load the Cell Atlas Raw miRNA counts and rectify data types.
#
Cell.Atlas <- read.delim2("./input/Merged_results_joel/exceRpt_miRNA_ReadCounts.txt",
                          sep = "\t",row.names = 1,stringsAsFactors = F)

Cell.Atlas <- to.numeric.df(tFrame(Cell.Atlas))


QN <- TRUE
if(QN){
    CA.Qnorm <- normalize.quantiles( as.matrix(Cell.Atlas) )
    Cell.Atlas.QN <- as.data.frame(CA.Qnorm)
    dimnames(Cell.Atlas.QN) <- dimnames(Cell.Atlas)
    Cell.Atlas <- Cell.Atlas.QN
}

head(Cell.Atlas[,1:4])
####################################################################################################
#Load Lilli's R environment from the graphical analysis:
matchedAnalysis <- F
exp.files <- list.files(path = './input/Lilli/', pattern = '^[ek][xr][ra][ns]' )
exp.files <- naturalsort(exp.files)
#ecxlude permuted network experiments
exp.files <- exp.files[-grep('[789]|[1][012]|serum',exp.files,value = F)]

sum.stat.cols <- c('p.cut','n.sigMods', 'largest.sigMod', 'n.sigNodes')
rn <- gsub('Prostrate','Prostate', (gsub('.RData','',(make.names(exp.files)))))
my.summary.stats <- data.frame(matrix(nrow = length(exp.files),ncol = length(sum.stat.cols)),
                               row.names = rn )
colnames(my.summary.stats) <- sum.stat.cols

out.put.lst <- list()
#smpl <- sample(1:length(exp.files),size = 6)
smpl <- 1:length( exp.files )
for(my.file in exp.files[smpl]){
    load(file = paste('./input/Lilli/', my.file, sep = '' ))
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
        print(paste( my.file,"Has Failed; no significant nodes"))
        next
    } 

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
    
    #CALL GET.INFORMATIVE.MODULES HERE  <<===
    ##################################
    makeplot <-function(mods,mod.membs,expt.name){
        if(all( rownames(mods) == names(mod.membs)) ){
            if(length(mod.membs)>1){
                #print('plottable')
                a <- sapply(mod.membs,length)
                b <- info2p(mods,reverse = T)
                p.plot <- qplot(x = a, y = b, asp = 9/16,
                                main = paste("Module Information by Size for\n experiment", expt.name),
                                xlab = "Module Size", ylab = "Module Fisher Information")
                
            }
        }
        
    }
    p.plot <- makeplot(mods = sig.mods,mod.membs = module.members,expt.name = expt.name)
    pdf(paste('./output/Cell_Atlas/p.vs.modSize.',expt.name,'.pdf'))
    print(p.plot)
    dev.off()
    
    #print("**")
    
    modsInfo <- get.mods.and.info(sig.mods,module.members)
    
    out.put.lst[[expt.name]] <- list(intra.cors = cors.im,
                                   extra.cors = cors.ex,
                                   bg.cors = cors.bg,
                                   bg.smpl = bg.smpl,
                                   mods.and.infos = modsInfo,
                                   node.pairs = module.node.pairings
                                   )
    
    my.summary.stats[expt.name,'p.cut'] <- p.cut
    my.summary.stats[expt.name,'n.sigMods'] <- length(sig.mods)
    my.summary.stats[expt.name,'largest.sigMod'] <- max(sapply(module.members,length))
    my.summary.stats[expt.name,'n.sigNodes'] <- length(all.sig.nodes)
    #my.summary.stats[expt.name,'Nodes/Mods'] <- length(all.sig.nodes)/length(sig.mods)
}

short.stats <- apply(my.summary.stats,MARGIN = 2,quantile, na.rm = T)

source('./graphNames.R')
graphSamples <- graphSamples[-grep("serum",graphSamples$`Derived from`),]
rownames(graphSamples) <- graphSamples[,1]

graphSamples <- graphSamples[rownames(my.summary.stats),]
all(rownames(graphSamples) == rownames(my.summary.stats))
cmbd.stats <- cbind(graphSamples,my.summary.stats)
dim(cmbd.stats)
write.table(cmbd.stats,'./output/Cell_Atlas/10.26.ModuleSummaryStats.tab',
            quote = F,row.names = F,col.names = T,sep = '\t', na = 'NA')
write.table(short.stats,'./output/Cell_Atlas/OverallpropertiesStats.csv',
            quote = F,row.names = T,col.names = T,sep = ',')
saveRDS(object = out.put.lst,file = './CellAtlasMod_OutLst10.26.rdat')
#out.put.lst <- readRDS('./CellAtlasModQN_OutLst.rdat')
if(matchedAnalysis){ExpPrs <- matched.expts(names(out.put.lst)) }
#source('./graphNames.R')
####################################################################################################
if(F){
dta.st = 0
source('./graphNames.R')
for(expt in names(ExpPrs)){
    dta.st = dta.st + 1
    print(expt)
    tru.ex <- ExpPrs[[expt]]$true
    rnd.ex <- ExpPrs[[expt]]$randomized
    cors.im <- NULL
    cors.rnd <- NULL
    cors.bg <- NULL
    if(!is.null(tru.ex) ){
        cors.im <- out.put.lst[[tru.ex]]$intra.cors
    }else{ print(paste("true module NULL error:",expt));next }
    
    if(!is.null(rnd.ex) ){
        cors.rnd <- out.put.lst[[rnd.ex]]$intra.cors
    }else{ print(paste("random module NULL error:",expt));next } 
    cors.bg <- c(out.put.lst[[tru.ex]]$bg.cors)#, out.put.lst[[rnd.ex]]$bg.cors)
 
    
    print('plotting')
    txp <-.3
    hst.brks <- 50
    ks.tru <- NULL
    ks.rnd <- NULL
    if(!is.null(cors.im)){ ks.tru <- ks.test(cors.im, cors.bg,alternative = "less")}
    if(!is.null(cors.rnd)){ ks.rnd <- ks.test(cors.rnd, cors.bg,alternative = "less")}
    options(scipen = 5)
    options(digits=4)
    
expt.lbl <- graphSamples[grep]
    out.file.name <- paste(graphSamples[dta.st,2],'--',graphSamples[dta.st,3],sep = '')
    pdf(file = paste('./output/Cell_Atlas/QN',out.file.name,'.pdf',sep = '' ) )
    
    hist(cors.im,col = rgb(1,0,0,txp),breaks = hst.brks,
         xlim = c(-0.5,1),ylim = c(0,15),freq = F,
         main = paste('Cell-Atlas miRNA intra-module expression correlations\n',
                      'Graph: ', graphSamples[dta.st,2],'\n',
                      'Applied to: ', graphSamples[dta.st,3],' patients', sep = ''),
         xlab = 'Pairwise correlation', ylab = 'Probability density as %' )
    
    hist(cors.rnd,col = rgb(0,1,0,txp),add = T, breaks = hst.brks, freq = F)
    hist(cors.bg,col = rgb(0,0,1,txp),add = T, breaks = hst.brks, freq = F)
    
    legend(x = 'topleft',legend = c(paste('true-module n =',length(cors.im)),
                                    paste('permuted-module n =',length(cors.rnd)),
                                    paste('non-module n =',bg.smpl)) ,
           fill = c(rgb(1,0,0,txp),rgb(0,1,0,txp), rgb(0,0,1,txp)) )
    mtext(paste('True KS-test D:',round(ks.tru$statistic,2)), line = -1,side = 3,adj = 1)
    mtext(paste('p-val:',signif(ks.tru$p.value,3)), line = -2,side = 3,adj = 1)
    mtext(paste('Random KS-test D:',round(ks.rnd$statistic,2)), line = -3,side = 3,adj = 1)
    mtext(paste('p-val:',signif(ks.rnd$p.value,3)), line = -4,side = 3,adj = 1)
    #mtext(paste('nod:',signif(ks.res$p.value,3)), line = -3,side = 3,adj = 1)
    
    dev.off()
}

}
####################################################################################################
#get the plots for the true-intra vs. true-extramodule
if(T){
dta.st = 0
out.file.name <- "CellATLASTSTJUMK"
pdf(file = paste('./output/Cell_Atlas/in_ex',out.file.name,'.pdf',sep = '' ))

expt.names <- gsub('.RData','', make.names(exp.files[smpl]) )
expt.names <-gsub('Prostrate','Prostate',expt.names)    


for(expt in expt.names){
    dta.st = dta.st + 1
    cors.im <- NULL
    cors.ex <- NULL
    cors.bg <- NULL
    if(!is.null(out.put.lst[[expt]]) ){
        cors.im <- out.put.lst[[expt]]$intra.cors
        cors.ex <- out.put.lst[[expt]]$extra.cors
    }else{ print(paste("true module NULL error:",expt));next }
    cors.bg <- c(out.put.lst[[expt]]$bg.cors)
    txp <-.3
    hst.brks <- 50
    ks.im <- NULL
    ks.ex <- NULL
    if(!is.null(cors.im)){ ks.tru <- ks.test(cors.im, cors.ex,alternative = "less")}
    #if(!is.null(cors.im)){ wlcx.tru <- wilcox.test(cors.im, cors.ex,alternative = "greater",
    #                                               paired = F,exact = T)}
    #if(!is.null(cors.ex)){ ks.ex <- ks.test(cors.ex, cors.bg,alternative = "less")}
    options(scipen = 5)
    options(digits=4)
    
    hist(cors.im,col = rgb(1,0,0,txp),breaks = hst.brks,
         xlim = c(-0.5,1),ylim = c(0,15),freq = F,
         main = paste('Cell-Atlas pairwise miRNA expression correlations\n',
                      'Graph: ', graphSamples[dta.st,2],'\n',
                      'Applied to: ', graphSamples[dta.st,3],' patients', sep = ''),
         xlab = 'Pairwise correlation', ylab = 'Probability density as %' )
    
    hist(cors.ex,col = rgb(0,0,1,txp),add = T, breaks = hst.brks, freq = F)
    legend(x = 'topleft',legend = c(paste('within-modules n =',length(cors.im)),
                                    paste('between-modules n =',length(cors.ex)) ),
           fill = c(rgb(1,0,0,txp), rgb(0,0,1,txp)) )
    
    mtext(paste('True KS-test D:',round(ks.tru$statistic,2)), line = -1,side = 3,adj = 1)
    mtext(paste('p-val:',signif(ks.tru$p.value,3)), line = -2,side = 3,adj = 1)
    #mtext(paste('Permutation test:',round(tic,2)), line = -3,side = 3,adj = 1)
    #mtext(paste('p-val:',signif(wlcx.tru$p.value,3)), line = -4,side = 3,adj = 1)
    
    
}
dev.off()
}


