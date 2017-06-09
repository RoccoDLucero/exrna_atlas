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
####################################################################################################
#Load the Atlas Raw miRNA counts and associated sample metaData for ~750 samples of various
# conditions, rectify data types, and merge into a single data frame
Atlas.miRNA.meta <- read.delim2("./input/Atlas_Human_Sample_Metadata2.txt",sep ="\t",row.names = 1)
Atlas.miRNA.raw <- read.delim2("./input/Atlas_Human_miRNA_Raw_Read_Counts.txt",sep = "\t",row.names = 1)
rn <- rownames(Atlas.miRNA.raw)
Atlas.miRNA.raw <- sapply(Atlas.miRNA.raw,as.numeric)
rownames(Atlas.miRNA.raw) <- rn
Atlas.miRNA.raw <- as.data.frame(Atlas.miRNA.raw)
Atlas.miRNA <- merge.data.frame(Atlas.miRNA.meta, t(Atlas.miRNA.raw), by.x = 0, by.y = 0 )
rm(Atlas.miRNA.meta,Atlas.miRNA.raw)
####################################################################################################
#Load Lilli's R environment from the graphical analysis:
exp.files <- list.files(path = './input', pattern = '^exrna' )

for( my.file in exp.files){
    
    load(file = paste('./input/', my.file,sep = '' ))
    #The matrix called 'optimize.results' contains Fisher information scores for 
    # modules identified in the graph for the given comparison
    
    a <- optimize.results[order(optimize.results),,drop = F]
    a <-2^-a 
    a<- a[a<.01,,drop = F]
    
    module.members <- sapply(rownames(a),strsplit, "/")
    head(sort(table(unlist(module.members)),decreasing = T))
    
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
    
    
    #Later add the biofluid meta data and condition name ('condition:biofluid')
    # For each of the conditions represented, and add a color key...
    #Also print the total number of intranode correlations used
    pdf(file = paste('./output/',gsub('.RData',"",my.file),'intramod.by.cond.fxd.pdf'),
        width = 4,height = 4)
    bw <- 0.05
    txp <- .2
    my.cols <- c(rgb(1,0,0,txp),rgb(0,1,0,txp),rgb(.05,0,1,txp),rgb(.6,.5,.2,txp))
    h1 <- hist(cor.df$Carcinoma, breaks = seq(-1,1,bw), col = my.cols[1],
               freq = T , ylim = c(0,100),
               main = paste('Graph:',gsub('.RData',"",my.file)),
               xlab ='Atlas miRNA: intranode pairwise correlations')
    h2 <- hist(cor.df$Alzheimer, breaks = seq(-1,1,bw), col = my.cols[2], freq = T,add = T)
    h3 <- hist(cor.df$Parkinson, breaks = seq(-1,1,bw), col = my.cols[3], freq = T, add = T)
    h4 <- hist(cor.df$`Healthy Control`, breaks = seq(-1,1,bw), col = my.cols[4], freq = T, add = T)
    legend("topleft", legend = c('Carc','Alzh', 'Park', 'Healthy'),lwd = 5,
           col =  my.cols, cex = .6)
    
    dev.off()
    scatters <- F
    if(scatters){
    plot(mod.cnd.corr.df$prk.corr,mod.cnd.corr.df$carc.corr, xlim = c(0,1),ylim = c(0,1))
    fit <- lsfit(mod.cnd.corr.df$prk.corr,mod.cnd.corr.df$carc.corr)
    abline(fit)
    
    plot(mod.cnd.corr.df$prk.corr,mod.cnd.corr.df$alz.corr, xlim = c(0,1),ylim = c(0,1))
    fit <- lsfit(mod.cnd.corr.df$prk.corr,mod.cnd.corr.df$alz.corr)
    abline(fit)
    
    plot(mod.cnd.corr.df$alz.corr,mod.cnd.corr.df$carc.corr, xlim = c(0,1),ylim = c(0,1))
    fit <- lsfit(mod.cnd.corr.df$prk.corr,mod.cnd.corr.df$alz.corr)
    abline(fit)
    
    
    qqplot(mod.cnd.corr.df$prk.corr,mod.cnd.corr.df$carc.corr, xlim = c(0,1),ylim = c(0,1))
    qqplot(mod.cnd.corr.df$alz.corr,mod.cnd.corr.df$carc.corr, xlim = c(0,1),ylim = c(0,1))
    qqplot(mod.cnd.corr.df$prk.corr,mod.cnd.corr.df$alz.corr, xlim = c(0,1),ylim = c(0,1))
    }
}