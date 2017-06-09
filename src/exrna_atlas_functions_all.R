#All functions for the exrna_atlas R project

################################################################################
#Load required libraries:
library(perm)
library(lsr)
library(stringr)
library(RColorBrewer)
library(colorRamps)
library(naturalsort)
library(ggplot2)
library(gplots)
library(devtools)
library(graphics)
library(infotheo)
#devtools::install_github("kassambara/factoextra")
library(factoextra)
source("http://www.bioconductor.org/biocLite.R")
biocLite(c("Biobase","preprocessCore"),suppressUpdates = T)
library(Biobase)
library(preprocessCore)
library(huge)
library(igraph)
library(VennDiagram)
library(RCurl)
library(tsne)
set.seed(20161017)
################################################################################
#Define custom functinos for Atlas analysis workflows
#
#####Processing graph diffusion results from Lilli's current outputs------######
#Converts information scores from 'optimize.results' matrix 
# to matirx of p-values at or below some threshold
info2p <- function(gph.info.Mtx, p.cutoff = .05,cut.style = '<=', reverse = F){
    if(!is.numeric(gph.info.Mtx)){return(NULL)}
    if(reverse){a <- -log2(gph.info.Mtx); return(a)}
    else{a <- gph.info.Mtx[order(gph.info.Mtx),,drop = F]
        a <-2^-a
        switch(cut.style,
               '<'  = { a<- a[ a <  p.cutoff, ,drop = F]},
               '<=' = { a<- a[ a <= p.cutoff, ,drop = F]},
               '>'  = { a<- a[ a >  p.cutoff, ,drop = F]},
               '>=' = { a<- a[ a >= p.cutoff, ,drop = F]}
        )
        return(a)
    }
}

my.get.mode.info <- function(my.vec){
    a <- names(rev(sort(table(my.vec)))[1])
    b <- (rev(sort(table(my.vec)))[1])
    return(paste(a,' (',b,'x)',sep = ''))
}

#Splits modules as rownames to list of nodes in each module: 
get.mod.membs <- function(gph.Mtx,delimiter = "/"){
    mod.membs <- sapply(rownames(gph.Mtx), strsplit, delimiter)
    return(mod.membs)
}

#converts imported data to data frame with numeric columns
to.numeric.df <- function(df){
    for( cl in 1:ncol(df) ){
        df[,cl] <- as.numeric(as.character(df[,cl]))
    }
    chk <- all((apply(df,2,is.numeric)))
    if(chk){return(df)}else{print('problem with data conversion')}
}    


get.non.module.vars <- function(query.df,sig,non.sig,var.ptt = "hsa-mir"){
    #The valid set of nodes are those which are present in the query
    vld <- grep(var.ptt,colnames(query.df),ignore.case = T)
    vld <- colnames(query.df[,vld])
    vld <- naturalsort(vld)
    
    excl.sig    <- (length(sig)-length(intersect(sig,vld)))
    excl.nonsig <- (length(non.sig)-length(intersect(non.sig,vld)))
    
    #print( paste(excl.sig, "significant nodes not in query" ))
    #print( paste(excl.nonsig, "non-significant nodes not in query" ))

    non.module.vars <- setdiff( vld,sig )
    return(non.module.vars)
}

#This function should, for each node, give the list of all nodes
# 1.that participate in ANY significant module with the given node.
# 2.that never appear in a significant module with the given node.
per.node.int.v.ext.nodes <- function(allNodes,modMembs){
    #For every node
    int.ext.lst <- list()
    for(cur.nd in allNodes){
        intramod <- c()
        extramod <- c()
        for(tst.mod in modMembs){
            if( cur.nd %in% tst.mod ){
                intramod <- c(intramod,tst.mod)
            }
            if( !(cur.nd %in% tst.mod) ){
                extramod <- c(extramod,tst.mod)
            }
        }
        intramod <- unique(intramod)
        extramod <- unique(extramod)
        extramod <- setdiff(extramod,intramod)
        int.ext.lst[[make.names(cur.nd)]] <- list(node = cur.nd,
                                                  intramod = naturalsort(intramod),
                                                  extramod = naturalsort(extramod))
    }
    return(int.ext.lst)
}


get.reduced.pairings <- function(in.ex.pairs.lst){                         
        all.im.prs <- list()
        all.ex.prs <- list()
        for(m in in.ex.pairs.lst){
            if(length(m$intramod) > 1){
                #For each significant node enumerate all pairs of
                # nodes that appear together in a 'significant' module
                all.im.prs <- c( all.im.prs, combn(m$intramod,2,simplify = F) )
            }
            if(length(m$extramod) > 1){
                #For each significant node enumerate all pairs of
                # nodes that never appear together in a 'significant' module
                # though both nodes appear in at least one 'significant' module
                all.ex.prs <- c( all.ex.prs, combn(m$extramod,2,simplify = F) )
            }
            
        }
        #Get rid of any node pairs in the extramodule group that are
        # present in the intramodule pairings
        all.im.prs <- all.im.prs[!duplicated(all.im.prs)]
        all.ex.prs <- all.ex.prs[!duplicated(all.ex.prs)]
        all.ex.prs <- all.ex.prs[!(all.ex.prs %in% all.im.prs) ]
        #Collapse the module pairs to a table that records a given
        # pairing only once, and sort the table
        im <- as.data.frame(t(sapply(all.im.prs,naturalsort)))
        im <- im[!duplicated(im),]    
        im <- im[naturalorder(im[,1]),]
        
        ex <- as.data.frame(t(sapply(all.ex.prs,naturalsort)))
        ex <- ex[!duplicated(ex),]
        ex <- ex[naturalorder(ex[,1]),]
        
        im.ex <- list(intramodule.prs = im, extramodule.prs = ex)
        return(im.ex)
}    

get.bg.pairings <- function(non.mod.vars,bg.sample = 1000, bg.seed = 20161017){
    set.seed(bg.seed)
    bg.prs <- combn(non.mod.vars,2,simplify = F) 
    bg <- sample(bg.prs, bg.sample)
    
    bg <- as.data.frame(t(sapply(bg,naturalsort)))
    bg <- bg[!duplicated(bg),]
    bg <- bg[naturalorder(bg[,1]),]
    return(bg)
}

get.nodepair.cors <- function(query.df, node.pairings){
    prs <-node.pairings
    cors <- c()
    for(i in 1:nrow(prs)){
        if(length(prs[1,])==2){
            m1 <- as.character(prs[i,1])
            m2 <- as.character(prs[i,2])
            if((m1 %in% colnames(query.df)) && (m2 %in% colnames(query.df)) ){
                cors <- c(cors, cor(query.df[,m1], query.df[,m2])) 
            }else{cors <- c(cors,NA)} 
        }else{ cors <- c(cors,NA)}
    }
    return(cors)
}

###################################################################################################
#This function should look at the names in the out.put.lst
# and match the true and random network experiments
# currently this is # and #+6 based on Lilli's convention for non- and permuted- networks
matched.expts <- function(expt.names){
    exp.categories <- c('exrna','kras')
    conditions <- c('AD$','AD-serum','PD$','PD-serum','Colo', 'Pros','Panc')
    nms <- naturalsort(names(out.put.lst))
    lst <- list()
    t <- 0
    for(cat in exp.categories){
        nn <- grep(cat,nms,value = T)
        for(cnd in conditions){
            ex <- grep(cnd,nn,value = T)
            if(length(ex)>0){
                for(p in 1:6){
                    a <- grep(paste('[as]',p,'-',sep = ''),ex,fixed = F,value = T)
                    b <- grep(p+6,ex,fixed = T,value = T)
                    if(length(a) == 0 ){a <- NULL}
                    if(length(b) == 0 ){b <- NULL} 
                    t <- t+1
                    res <-  list(true = a, randomized = b)
                    lst[[paste(cat,cnd,p,sep = '')]] <- res
                }
            }    
        }
    }
    for(i in names(lst)){
        if(all(sapply(lst[[i]],is.null))){lst[[i]] <-NULL }
    }
    return(lst)
}

harmonic.mean <- function(x){
    if(is.numeric(x)){
    hm <-  1/mean(1/x)
    }
    else{"check input vector to harmonic mean"
        hm <-NULL
    }
    return(hm)
}


get.mods.and.info <- function(sig.mods,mod.membs){
    #Returns the top modules and plots all input modules vs. fisher informatoin
    #Plots the p-value vs size for each module identified in a given experiment
    #Displays cutoffs for size and pvalue
    #This should help us to identify modules with the most potential
    # for biomarkers and will help in prioritizing per module analysis
    #Weights for information and module size can be used to fine tune the model
    # these may later be determined by a linaear model
    #This funcion is a hackjob at this point..should be streamlined later...
    if(all( rownames(sig.mods) == names(mod.membs)) ){
        a <- sapply(mod.membs,length)
        b <- info2p(sig.mods,reverse = T)
        tm <- as.data.frame(cbind(rownames(b),b),row.names = F)
        tm[,1] <- as.character(tm[,1])
        tm[,2] <- as.numeric(as.character(tm[,2]))
        tm <- tm[order(tm[,2],decreasing = T),]
        #tm <- tm[1:n.top,]
        x  <- sapply(tm[,1], strsplit, split = "/")
        y  <- tm[,2]
        tm2 <-list()
        for(s in 1:length(x)){
            nm <- gsub('hsa','',make.names(names(x)[s]))
            tm2[[ nm ]] <- list(members = x[[s]], information = y[s])
            
        }
    }
    return(tm2)
}

################################################################################
#This function removes modules which are subsets of other modules. 
rmv.sub.modules <- function(mds.lst){
    mds <- module.members#mds.lst
    md.red <- vector("list",length(mds))
    i <- 0
    for(mod.1 in mds){
        i <- i+1
        sub <- F
        for(mod.2 in mds[-i]){
            if(all( mod.1 %in% mod.2)){
                sub <- T
                if(sub){break}
            }
        }
        if(!sub){md.red[[i]] <- naturalsort(mod.1)}
    }
    md.red <- md.red[!sapply(md.red,is.null)]
    return( unique(md.red) )
}


#Remove isolated nodes from a graph
delete.isolated <- function(graph, mode = 'all') {
    isolates <- which(degree(graph, mode = mode) == 0) - 1
    delete.vertices(graph, isolates)
}

